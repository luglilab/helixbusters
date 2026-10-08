"""BLISS 5-prime end counting and UMI deduplication on coordinate-sorted BAMs.

Coordinates describe the terminal aligned base, not an interbase cut position.
Each input BAM must represent one biological library. No spatial merging is done.
"""

from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
import heapq
import json
import os
import tempfile

import pysam


@dataclass(frozen=True)
class UmiFamily:
    representative: str
    members: tuple[str, ...]
    read_count: int


def group_umis(counts, method="exact", max_distance=1):
    """Group ACGT UMIs at ONE chromosome/position/strand.

    Directional edges follow count(A) >= 2 * count(B) - 1. Descending abundance
    and lexical ties make representatives reproducible. Distance is Hamming;
    transitive paths can connect UMIs farther apart than max_distance.
    """
    if method not in {"exact", "directional"}:
        raise ValueError("method must be 'exact' or 'directional'")
    if type(max_distance) is not int or not 0 <= max_distance <= 2:
        raise ValueError("max_distance must be 0, 1 or 2")
    lengths = set()
    for umi, count in counts.items():
        if not isinstance(umi, str) or not umi or set(umi) - set("ACGT"):
            raise ValueError("UMIs must be nonempty ACGT strings")
        if type(count) is not int or count < 1:
            raise ValueError("UMI counts must be positive integers")
        lengths.add(len(umi))
    if len(lengths) > 1:
        raise ValueError("UMIs at a locus must have the same length")
    ordered = sorted(counts, key=lambda umi: (-counts[umi], umi))
    if method == "exact" or max_distance == 0 or len(counts) < 2:
        return [UmiFamily(umi, (umi,), counts[umi]) for umi in ordered]

    def neighbors(umi):
        # Enumerate substitutions instead of comparing every pair at dense sites.
        for distance in range(1, min(max_distance, len(umi)) + 1):
            for positions in combinations(range(len(umi)), distance):
                choices = [tuple(base for base in "ACGT" if base != umi[p])
                           for p in positions]
                for replacements in product(*choices):
                    sequence = list(umi)
                    for position, base in zip(positions, replacements):
                        sequence[position] = base
                    candidate = "".join(sequence)
                    if candidate in counts:
                        yield candidate

    assigned = set()
    families = []
    for root in ordered:
        if root in assigned:
            continue
        # Traverse the full directed component, including previously assigned
        # nodes; only newly assigned nodes belong to this family.
        visited = {root}
        stack = [root]
        while stack:
            current = stack.pop()
            for candidate in neighbors(current):
                if (candidate not in visited and
                        counts[current] >= 2 * counts[candidate] - 1):
                    visited.add(candidate)
                    stack.append(candidate)
        members = tuple(sorted(visited - assigned))
        assigned.update(members)
        families.append(UmiFamily(root, members, sum(counts[u] for u in members)))
    return families


def five_prime_position(read, policy="strict"):
    """Return a zero-based terminal-base coordinate, or None if ambiguous.

    strict requires an aligned base (M, = or X) at the biological 5-prime CIGAR
    end. aligned uses the alignment boundary even if that end is clipped.
    Neither policy extrapolates soft clips into an unobserved genomic sequence.
    Spliced alignments are excluded in both policies for this genomic assay.
    """
    if policy not in {"strict", "aligned"}:
        raise ValueError("five_prime_policy must be 'strict' or 'aligned'")
    if read.is_unmapped or not read.cigartuples or read.reference_start < 0:
        return None
    if any(op == 3 for op, length in read.cigartuples):
        return None
    if read.reference_end is None or read.reference_end <= read.reference_start:
        return None
    terminal_op = read.cigartuples[-1 if read.is_reverse else 0][0]
    if policy == "strict" and terminal_op not in {0, 7, 8}:
        return None
    return read.reference_end - 1 if read.is_reverse else read.reference_start


def _read_umi(read, umi_length):
    # RX is authoritative when present; never silently replace an invalid tag
    # with a potentially unrelated read-name suffix.
    if read.has_tag("RX"):
        umi = read.get_tag("RX")
    else:
        name = read.query_name or ""
        if "_" not in name:
            return None, "missing_umi"
        umi = name.rsplit("_", 1)[1]
    if not isinstance(umi, str) or len(umi) != umi_length or set(umi) - set("ACGT"):
        return None, "invalid_umi"
    return umi, None


def deduplicate_bam(
    input_bam, output_file_umi_pcr, output_file_umi_count, *,
    method="exact", umi_length=8, max_distance=1, min_mapq=0,
    five_prime_policy="strict", read_selection="single-end",
    output_molecules=None, output_sites=None, output_qc=None,
):
    """Write UMI families and end counts; return QC and parameter provenance.

    Existing positional outputs remain six-column family TSV and four-column
    BED (counts SUMMED across strands). Optional sites TSV retains strand and
    molecule BED6 emits one row per inferred molecule for downstream callers.
    No BAM index is required, but actual coordinate order is validated.
    """
    group_umis({}, method, max_distance)  # Validate even for empty BAMs.
    if type(umi_length) is not int or umi_length < 1:
        raise ValueError("umi_length must be a positive integer")
    if type(min_mapq) is not int or not 0 <= min_mapq <= 254:
        raise ValueError("min_mapq must be between 0 and 254")
    if five_prime_policy not in {"strict", "aligned"}:
        raise ValueError("five_prime_policy must be 'strict' or 'aligned'")
    if read_selection not in {"single-end", "read1", "read2"}:
        raise ValueError("read_selection must be 'single-end', 'read1' or 'read2'")
    if output_file_umi_pcr is None or output_file_umi_count is None:
        raise ValueError("Both family and count output paths are required")

    destinations = {"families": output_file_umi_pcr, "counts": output_file_umi_count,
                    "molecules": output_molecules, "sites": output_sites, "qc": output_qc}
    destinations = {name: Path(path) for name, path in destinations.items() if path is not None}
    paths = [Path(input_bam).resolve()] + [path.resolve() for path in destinations.values()]
    if len(paths) != len(set(paths)):
        raise ValueError("Input and output paths must all be distinct")

    stats = {key: 0 for key in ("input_reads", "accepted_reads", "exact_molecules",
                               "deduplicated_molecules", "strand_loci", "positions")}
    skipped = Counter()
    pending = {}
    positions = []
    temporaries = {}
    handles = {}
    try:
        with ExitStack() as stack:
            bam = stack.enter_context(pysam.AlignmentFile(str(input_bam), "rb"))
            samples = {rg["SM"] for rg in bam.header.to_dict().get("RG", []) if rg.get("SM")}
            if len(samples) > 1:
                raise ValueError("BAM contains multiple SM samples; deduplicate each library separately")
            for name, path in destinations.items():
                path.parent.mkdir(parents=True, exist_ok=True)
                handle = stack.enter_context(tempfile.NamedTemporaryFile(
                    mode="w", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False))
                handles[name] = handle
                temporaries[name] = Path(handle.name)
            if "sites" in handles:
                handles["sites"].write("chrom\tstart\tend\tstrand\tmolecules\treads\n")

            def flush(chrom, before=None):
                while positions and (before is None or positions[0] < before):
                    position = heapq.heappop(positions)
                    strand_counts = pending.pop(position)
                    total = 0
                    for strand, counts in sorted(strand_counts.items()):
                        families = group_umis(counts, method, max_distance)
                        stats["strand_loci"] += 1
                        stats["exact_molecules"] += len(counts)
                        stats["deduplicated_molecules"] += len(families)
                        total += len(families)
                        for family in families:
                            handles["families"].write(
                                f"{chrom}\t{position}\t{position + 1}\t{strand}\t"
                                f"{family.representative}\t{family.read_count}\n")
                            if "molecules" in handles:
                                handles["molecules"].write(
                                    f"{chrom}\t{position}\t{position + 1}\t"
                                    f"{family.representative}\t0\t{strand}\n")
                        if "sites" in handles:
                            handles["sites"].write(
                                f"{chrom}\t{position}\t{position + 1}\t{strand}\t"
                                f"{len(families)}\t{sum(counts.values())}\n")
                    handles["counts"].write(f"{chrom}\t{position}\t{position + 1}\t{total}\n")
                    stats["positions"] += 1

            previous = None
            current_reference = None
            unplaced_seen = False
            for read in bam.fetch(until_eof=True):
                stats["input_reads"] += 1
                # Coordinate order applies to placed unmapped records too.
                if read.reference_id >= 0 and read.reference_start >= 0:
                    key = (read.reference_id, read.reference_start)
                    if unplaced_seen or (previous is not None and key < previous):
                        raise ValueError("Input BAM must be coordinate-sorted")
                    if current_reference is not None and read.reference_id != current_reference:
                        flush(bam.get_reference_name(current_reference))
                    current_reference = read.reference_id
                    previous = key
                    flush(bam.get_reference_name(current_reference), read.reference_start)
                else:
                    unplaced_seen = True

                reason = None
                if read.is_unmapped:
                    reason = "unmapped"
                elif read.is_secondary:
                    reason = "secondary"
                elif read.is_supplementary:
                    reason = "supplementary"
                elif read.is_qcfail:
                    reason = "qc_fail"
                elif read.mapping_quality == 255:
                    reason = "unknown_mapq"
                elif read.mapping_quality < min_mapq:
                    reason = "low_mapq"
                if reason:
                    skipped[reason] += 1
                    continue

                if read_selection == "single-end" and read.is_paired:
                    raise ValueError("Paired reads require explicit read_selection='read1' or 'read2'")
                if read_selection != "single-end":
                    if not read.is_paired or read.is_read1 == read.is_read2:
                        raise ValueError("Mate selection requires valid paired-read flags")
                    if (read_selection == "read1" and not read.is_read1 or
                            read_selection == "read2" and not read.is_read2):
                        skipped["other_mate"] += 1
                        continue

                position = five_prime_position(read, five_prime_policy)
                if position is None:
                    skipped["ambiguous_five_prime"] += 1
                    continue
                if position >= bam.lengths[read.reference_id]:
                    skipped["outside_reference"] += 1
                    continue
                umi, reason = _read_umi(read, umi_length)
                if reason:
                    skipped[reason] += 1
                    continue
                # Coordinate-only duplicate flags do not override UMI evidence.
                if position not in pending:
                    pending[position] = {}
                    heapq.heappush(positions, position)
                strand = "-" if read.is_reverse else "+"
                pending[position].setdefault(strand, Counter())[umi] += 1
                stats["accepted_reads"] += 1

            if current_reference is not None:
                flush(bam.get_reference_name(current_reference))
            if not stats["accepted_reads"] and (skipped["missing_umi"] or skipped["invalid_umi"]):
                raise ValueError("No usable UMIs found; check RX tags, read-name suffixes and umi_length")
            stats["duplicate_reads"] = stats["accepted_reads"] - stats["deduplicated_molecules"]
            stats["umi_groups_merged"] = stats["exact_molecules"] - stats["deduplicated_molecules"]
            stats["skipped"] = dict(sorted(skipped.items()))
            stats["parameters"] = {
                "method": method, "umi_length": umi_length,
                "max_distance": max_distance if method == "directional" else 0,
                "min_mapq": min_mapq, "five_prime_policy": five_prime_policy,
                "read_selection": read_selection,
                "coordinates": "0-based terminal aligned base; half-open BED intervals",
            }
            if "qc" in handles:
                json.dump(stats, handles["qc"], indent=2)
                handles["qc"].write("\n")
        # Publish only after a complete successful read and close of the BAM.
        for name, path in destinations.items():
            os.replace(temporaries[name], path)
        return stats
    finally:
        for path in temporaries.values():
            path.unlink(missing_ok=True)
