#!/usr/bin/env python3
"""Inspect original FASTQ layouts without trimming or modifying input files."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import sys


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_manifest(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "group", "replicate", "barcode", "fastq"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Manifest requires tab-separated sample, group, replicate, barcode, fastq columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Manifest has no samples")
    seen = set()
    for row in rows:
        if any(not row.get(key) for key in required):
            raise ValueError("Manifest contains a missing required value")
        if row["sample"] in seen:
            raise ValueError(f"Duplicate sample: {row['sample']}")
        seen.add(row["sample"])
        row["barcode"] = row["barcode"].upper()
        if set(row["barcode"]) - set("ACGT"):
            raise ValueError(f"Invalid barcode for {row['sample']}")
        fastq = Path(row["fastq"]).expanduser()
        row["fastq"] = str((path.parent / fastq).resolve())
    return rows


def fastq_sequences(path, limit):
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        for number in range(1, limit + 1):
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().rstrip("\r\n").upper()
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if (not header.startswith("@") or not plus.startswith("+") or
                    not sequence or len(sequence) != len(quality) or
                    set(sequence) - set("ACGTN")):
                raise ValueError(f"Malformed four-line FASTQ record {number} in {path}")
            yield sequence


def inspect_fastq(path, barcodes, limit, scan_bases, prefix_bases, offset, umi_length):
    lengths, prefixes = Counter(), Counter()
    composition = [Counter() for _ in range(scan_bases)]
    patterns = [(barcode, orientation, pattern) for barcode in sorted(set(barcodes))
                for orientation, pattern in (("forward", barcode),
                                             ("reverse_complement", reverse_complement(barcode)))]
    positions = {pair[:2]: Counter() for pair in patterns}
    matched_reads = Counter()
    offset_matches = {barcode: Counter() for barcode in set(barcodes)}
    umi_counts = Counter()
    n = 0
    for sequence in fastq_sequences(path, limit):
        n += 1
        lengths[len(sequence)] += 1
        prefixes[sequence[:prefix_bases]] += 1
        for i, base in enumerate(sequence[:scan_bases]):
            composition[i][base] += 1
        window = sequence[:scan_bases]
        for barcode, orientation, pattern in patterns:
            # Count each position per read; retain all occurrences, not just the first.
            hits = []
            position = window.find(pattern)
            while position >= 0:
                hits.append(position)
                position = window.find(pattern, position + 1)
            if hits:
                matched_reads[(barcode, orientation)] += 1
                positions[(barcode, orientation)].update(hits)
        for barcode in offset_matches:
            observed = sequence[offset:offset + len(barcode)]
            if len(observed) == len(barcode):
                offset_matches[barcode]["eligible_reads"] += 1
                distance = sum(a != b for a, b in zip(observed, barcode))
                offset_matches[barcode]["exact_reads"] += distance == 0
                offset_matches[barcode]["within_one_mismatch_reads"] += distance <= 1
        if len(sequence) >= umi_length:
            umi_counts[sequence[:umi_length]] += 1
    if n == 0:
        raise ValueError(f"Empty FASTQ: {path}")
    barcode_results = []
    for barcode, orientation, _ in patterns:
        key = (barcode, orientation)
        barcode_results.append({
            "barcode": barcode, "orientation": orientation,
            "matching_reads": matched_reads[key],
            "matching_pct": 100 * matched_reads[key] / n,
            "positions": [{"offset_0based": i, "reads": count, "pct_all_reads": 100 * count / n}
                          for i, count in sorted(positions[key].items())],
        })
    return {
        "fastq": str(path), "compressed_or_plain_file_bytes": path.stat().st_size,
        "reads_examined": n, "read_length_histogram": dict(sorted(lengths.items())),
        "top_prefixes": [{"sequence": seq, "reads": count, "pct": 100 * count / n}
                         for seq, count in prefixes.most_common(10)],
        "base_composition": [{"offset_0based": i, "reads_covering_position": sum(counts.values()),
                              "counts": dict(counts),
                              "entropy_bits": -sum((c / sum(counts.values())) * math.log2(c / sum(counts.values()))
                                                   for c in counts.values())}
                             for i, counts in enumerate(composition) if counts],
        "barcode_search": barcode_results,
        "forward_barcode_at_tested_offset": {barcode: dict(counts) for barcode, counts in offset_matches.items()},
        "first_bases_candidate_umi": {"length": umi_length, "distinct_sequences": len(umi_counts),
                                      "top_sequences": umi_counts.most_common(10)},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New JSON report; existing files are refused")
    parser.add_argument("--max-reads", type=int, default=100000, help="First reads per unique FASTQ (default: 100000)")
    parser.add_argument("--scan-bases", type=int, default=80, help="Search within this many initial bases")
    parser.add_argument("--prefix-bases", type=int, default=20)
    parser.add_argument("--barcode-offset", type=int, default=8, help="Zero-based offset to test, not an assumed layout")
    parser.add_argument("--umi-length", type=int, default=8, help="Initial candidate UMI length to describe")
    args = parser.parse_args()
    if min(args.max_reads, args.scan_bases, args.prefix_bases, args.umi_length) < 1 or args.barcode_offset < 0:
        parser.error("Lengths and max-reads must be positive; barcode-offset must be nonnegative")
    if args.output.exists():
        parser.error(f"Output already exists: {args.output}")
    rows = read_manifest(args.manifest.expanduser().resolve())
    by_path = defaultdict(list)
    for row in rows:
        by_path[row["fastq"]].append(row["sample"])
    print(f"Inspecting {len(by_path)} unique FASTQs for {len(rows)} samples, up to {args.max_reads:,} reads each.", flush=True)
    duplicates = {path: samples for path, samples in by_path.items() if len(samples) > 1}
    for path, samples in duplicates.items():
        print(f"WARNING: shared FASTQ for {', '.join(samples)}: {path}", file=sys.stderr)
    results = []
    for path, samples in by_path.items():
        print(f"Reading {', '.join(samples)}", flush=True)
        result = inspect_fastq(Path(path), [row["barcode"] for row in rows], args.max_reads,
                               args.scan_bases, args.prefix_bases, args.barcode_offset, args.umi_length)
        result["samples_using_path"] = samples
        results.append(result)
        for row in (row for row in rows if row["fastq"] == path):
            exact = result["forward_barcode_at_tested_offset"][row["barcode"]].get("exact_reads", 0)
            barcode_hits = [hit for hit in result["barcode_search"] if hit["barcode"] == row["barcode"]]
            peaks = [(position["reads"], hit["orientation"], position["offset_0based"])
                     for hit in barcode_hits for position in hit["positions"]]
            peak = max(peaks, default=None)
            print(f"  {row['sample']}: n={result['reads_examined']:,}; barcode={row['barcode']}; "
                  f"exact at offset {args.barcode_offset}={100 * exact / result['reads_examined']:.3f}%; "
                  f"strongest exact barcode peak (reads, orientation, offset)={peak}", flush=True)
    report = {
        "sampling": "First N records per unique input file; not a random sample or a full-library count",
        "interpretation": "Exact short-motif matches can occur by chance. Positions and diversity do not prove UMI identity.",
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "manifest": rows, "shared_fastq_paths": duplicates, "files": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(f"Report: {args.output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, EOFError, UnicodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
