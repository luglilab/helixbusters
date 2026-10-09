"""Mapping QC summaries and sparse BLISS end tracks for biological replicates."""

from contextlib import ExitStack
import csv
import heapq
import itertools
from importlib import metadata
import json
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys

import pysam


def validate_label(value):
    """Labels become directory names and Nextflow script arguments."""
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"Invalid sample/group/replicate label: {value!r}; use letters, digits, _, . or -")
    return value


def write_json(path, data):
    with Path(path).open("x") as handle:
        handle.write(json.dumps(data, indent=2, allow_nan=False) + "\n")


def check_environment(reference_config, genome, aligner, check_macs3=False):
    """Fail before mapping if worker tools or selected references are missing."""
    from helixbusters.genomes import load_reference_config, make_genome_filter
    from helixbusters.mapping import validate_index

    tools = (aligner, "samtools", "umi_tools", "bamCoverage", "multiqc")
    if check_macs3:
        tools += ('macs3',)
    missing = [name for name in tools if shutil.which(name) is None]
    if missing:
        raise ValueError("Missing worker executables: " + ", ".join(missing) +
                         "; activate/update the helixbusters environment")
    import pyBigWig  # Check the binary extension, not only distribution metadata.
    reference = load_reference_config(reference_config, genome, aligner)
    validate_index(reference["genome_index"], aligner)
    genome_filter = make_genome_filter(reference["genome"], reference["blacklist_bed"],
                                     reference["blacklist_genome"])
    versions = {name: metadata.version(name) for name in ("pysam", "deepTools", "multiqc")}
    versions["pyBigWig"] = pyBigWig.__version__
    versions["samtools"] = subprocess.check_output(["samtools", "--version"], text=True)
    write_json("environment.json", {"versions": versions, "reference": reference,
                                    "python_version": sys.version,
                                    "genome_filter": genome_filter.metadata(),
                                    "executables": {name: shutil.which(name) for name in tools}})


def metrics(mapping, dedup):
    """Recompute proportions from counts; denominators are explicit."""
    primary = mapping["primary_records"]
    retained = mapping["retained_records"]
    accepted = dedup["accepted_reads"]
    result = {key: mapping[key] for key in (
        "alignment_records", "primary_records", "mapped_primary_records", "retained_records")}
    reasons = {"unmapped", "secondary", "supplementary", "qc_fail", "unknown_mapq",
               "low_mapq", "mitochondrial", "noncanonical", "blacklist"}
    reasons.update(mapping["excluded_by_reason"])
    result.update({f"excluded_{key}": mapping["excluded_by_reason"].get(key, 0)
                   for key in sorted(reasons)})
    result.update({key: dedup[key] for key in (
        "accepted_reads", "deduplicated_molecules", "duplicate_reads")})
    result["mapping_rate_pct"] = 100 * mapping["mapped_primary_records"] / primary if primary else None
    result["retained_of_primary_pct"] = 100 * retained / primary if primary else None
    result["umi_duplicate_pct"] = 100 * dedup["duplicate_reads"] / accepted if accepted else None
    result["dedup_accepted_of_retained_pct"] = 100 * accepted / retained if retained else None
    result["ambiguous_five_prime_reads"] = dedup.get("skipped", {}).get("ambiguous_five_prime", 0)
    return result


def write_table(path, rows):
    columns = sorted(set().union(*(row.keys() for row in rows.values())))
    with Path(path).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample"] + columns, delimiter="\t")
        writer.writeheader()
        for sample, row in sorted(rows.items()):
            writer.writerow({"sample": sample, **row})


def bam_stats(bam, prefix):
    for command in ("flagstat", "stats", "idxstats"):
        with Path(f"{prefix}.{command}.txt").open("x") as output:
            subprocess.run(["samtools", command, str(bam)], stdout=output, check=True)


def bam_header(bam):
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        return list(zip(handle.references, handle.lengths))


def bed_rows(handle, header):
    """Validate sorted 1-bp end counts without loading a genome-sized array."""
    rank = {chrom: i for i, (chrom, _) in enumerate(header)}
    lengths = dict(header)
    previous = None
    for number, line in enumerate(handle, 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"Invalid counts BED row {number}")
        chrom, start, end, count = fields
        start, end, count = int(start), int(end), int(count)
        if chrom not in rank or start < 0 or end != start + 1 or end > lengths[chrom] or count < 1:
            raise ValueError(f"Invalid end count on BED row {number}")
        key = (rank[chrom], start)
        if previous is not None and key <= previous:
            raise ValueError("End counts must be unique and sorted in BAM header order")
        previous = key
        yield key, chrom, start, end, count


def merge_end_counts(paths, totals, header):
    """Sum independently deduplicated libraries and average individual CPMs.

    The pooled CPM denominator is the sum of molecules. An equal-weight mean
    is undefined if any biological replicate has no molecules.
    """
    if len(paths) != len(totals) or not paths:
        raise ValueError("One molecule total is required per nonempty list of libraries")
    if any(type(total) is not int or total < 0 for total in totals):
        raise ValueError("Molecule totals must be nonnegative integers")
    pooled_total = sum(totals)
    mean_available = all(totals)
    with ExitStack() as stack:
        streams = []
        for path, total in zip(paths, totals):
            def stream(path=path, total=total):
                observed = 0
                handle = stack.enter_context(Path(path).open())
                for key, chrom, start, end, count in bed_rows(handle, header):
                    observed += count
                    yield key, chrom, start, end, count, count * 1e6 / total if total else 0.0
                if observed != total:
                    raise ValueError(f"BED molecule sum {observed} differs from QC total {total}: {path}")
            streams.append(stream())
        for _, rows in itertools.groupby(heapq.merge(*streams), key=lambda row: row[0]):
            rows = list(rows)
            _, chrom, start, end, _, _ = rows[0]
            raw = sum(row[4] for row in rows)
            pooled = raw * 1e6 / pooled_total if pooled_total else 0.0
            mean = sum(row[5] for row in rows) / len(totals) if mean_available else None
            yield chrom, start, end, raw, pooled, mean


def end_tracks(paths, totals, header, prefix, mean=False):
    """Write sparse 1-bp tracks, raw pooled BED and track provenance."""
    import pyBigWig

    fields = [("raw", 3), ("CPM", 4)]
    if mean and all(totals):
        fields.append(("mean.CPM", 5))
    destinations = [Path(f"{prefix}.ends.{suffix}.bw") for suffix, _ in fields]
    destinations += [Path(f"{prefix}.ends.counts.bed"), Path(f"{prefix}.ends.provenance.json")]
    if any(path.exists() for path in destinations):
        raise ValueError(f"End-track outputs already exist for {prefix}; use a new output prefix")
    writers = []
    try:
        for suffix, column in fields:
            writer = pyBigWig.open(f"{prefix}.ends.{suffix}.bw", "w")
            writers.append((writer, column))
            writer.addHeader(header)
        batch = []

        def flush():
            if not batch:
                return
            for writer, column in writers:
                writer.addEntries([row[0] for row in batch], [row[1] for row in batch],
                                  ends=[row[2] for row in batch],
                                  values=[float(row[column]) for row in batch])
            batch.clear()

        with Path(f"{prefix}.ends.counts.bed").open("x") as bed:
            for row in merge_end_counts(paths, totals, header):
                bed.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n")
                batch.append(row)
                if len(batch) == 10000:
                    flush()
            flush()
    finally:
        for writer, _ in writers:
            writer.close()
    write_json(f"{prefix}.ends.provenance.json", {
        "inputs": [str(path) for path in paths], "molecule_totals": totals,
        "coordinate_system": "0-based half-open", "bin_size": 1,
        "signal": "UMI-deduplicated 5-prime ends, strands summed",
        "pooled_CPM": "sum of end counts / sum of library molecule totals * 1e6",
        "mean_CPM": "equal-weight mean of individually normalized replicates" if mean else None,
        "mean_CPM_available": bool(mean and all(totals)),
        "zero_molecule_libraries": sum(total == 0 for total in totals),
        "pyBigWig_version": pyBigWig.__version__,
    })


def coverage_track(bam, prefix, threads, bin_size):
    """Mapping diagnostic coverage includes PCR duplicates; no read extension."""
    if any(Path(f"{prefix}.{suffix}").exists() for suffix in ("coverage.CPM.bw", "coverage.provenance.json")):
        raise ValueError(f"Coverage outputs already exist for {prefix}; use a new output prefix")
    command = ["bamCoverage", "--bam", str(bam), "--outFileName", f"{prefix}.coverage.CPM.bw",
               "--binSize", str(bin_size), "--normalizeUsing", "CPM", "--exactScaling",
               "--numberOfProcessors", str(threads)]
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        nonempty = next(handle.fetch(until_eof=True), None) is not None
    if nonempty:
        subprocess.run(command, check=True)
    else:
        import pyBigWig
        writer = pyBigWig.open(f"{prefix}.coverage.CPM.bw", "w")
        try:
            writer.addHeader(bam_header(bam))
        finally:
            writer.close()
    write_json(f"{prefix}.coverage.provenance.json", {
        "command": command, "signal": "filtered alignment coverage before UMI deduplication",
        "empty_bam": not nonempty,
        "normalization": "CPM per mapped alignment; relative coverage, not absolute break burden",
        "bamCoverage_version": subprocess.check_output(["bamCoverage", "--version"], text=True).strip(),
    })


def sample_qc(sample, group, replicate, mapping_path, dedup_path, all_bam, filtered_bam,
              counts, threads=2, bin_size=50):
    for label in (sample, group, replicate):
        validate_label(label)
    mapping = json.loads(Path(mapping_path).read_text())
    dedup = json.loads(Path(dedup_path).read_text())
    if mapping["sample"] != sample or dedup["input_reads"] != mapping["retained_records"]:
        raise ValueError("Mapping and deduplication QC do not match this sample")
    for stage, bam in (("all", all_bam), ("filtered", filtered_bam)):
        bam_stats(bam, f"sample__{sample}.{stage}")
    summary = {"sample": sample, "group": group, "replicate": replicate,
               "metrics": metrics(mapping, dedup), "mapping": mapping, "deduplication": dedup}
    write_json(f"{sample}.summary.json", summary)
    write_json(f"{sample}.chrom.sizes.json", bam_header(filtered_bam))
    write_table(f"{sample}.qc.tsv", {sample: {"group": group, "replicate": replicate, **summary["metrics"]}})
    end_tracks([counts], [dedup["deduplicated_molecules"]], bam_header(filtered_bam), sample)
    coverage_track(filtered_bam, sample, threads, bin_size)


def group_qc(group, summaries, counts, headers, bams=None, threads=2, bin_size=50):
    validate_label(group)
    # All three lists are ordered by sample by the Nextflow caller.
    records = [json.loads(Path(path).read_text()) for path in summaries]
    sizes = [json.loads(Path(path).read_text()) for path in headers]
    if not records or len(records) != len(counts) or len(records) != len(sizes):
        raise ValueError("Group inputs must have matching lengths")
    if any(record["group"] != group for record in records):
        raise ValueError("Condition metadata mismatch")
    if len({record["sample"] for record in records}) != len(records):
        raise ValueError("Duplicate samples in condition")
    if len({record["replicate"] for record in records}) != len(records):
        raise ValueError("One library per biological replicate is required within a condition")
    if any(header != sizes[0] for header in sizes):
        raise ValueError("Replicate BAM headers have different contigs/order/lengths")
    mapping = {key: sum(record["mapping"][key] for record in records) for key in (
        "alignment_records", "primary_records", "mapped_primary_records", "retained_records")}
    reasons = set().union(*(record["mapping"]["excluded_by_reason"] for record in records))
    mapping["excluded_by_reason"] = {
        key: sum(record["mapping"]["excluded_by_reason"].get(key, 0) for record in records)
        for key in sorted(reasons)}
    dedup = {key: sum(record["deduplication"][key] for record in records) for key in (
        "accepted_reads", "deduplicated_molecules", "duplicate_reads")}
    dedup["skipped"] = {"ambiguous_five_prime": sum(
        record["deduplication"].get("skipped", {}).get("ambiguous_five_prime", 0) for record in records)}
    # Refuse pools mixing mapping or UMI rules even when contig headers match.
    if any(record["mapping"]["parameters"]["aligner"] != records[0]["mapping"]["parameters"]["aligner"] or
           record["mapping"]["parameters"]["min_mapq"] != records[0]["mapping"]["parameters"]["min_mapq"] or
           record["deduplication"]["parameters"] != records[0]["deduplication"]["parameters"] or
           record["mapping"]["genome_filter"] != records[0]["mapping"]["genome_filter"]
           for record in records):
        raise ValueError("Replicates have incompatible filtering/deduplication parameters")
    result = metrics(mapping, dedup)
    result["biological_replicates"] = len(records)
    for key in ("mapping_rate_pct", "retained_of_primary_pct", "umi_duplicate_pct"):
        values = [record["metrics"][key] for record in records if record["metrics"][key] is not None]
        result[f"replicate_mean_{key}"] = statistics.mean(values) if values else None
        result[f"replicate_sd_{key}"] = statistics.stdev(values) if len(values) > 1 else None
        result[f"replicate_n_{key}"] = len(values)
    summary = {"group": group, "samples": [record["sample"] for record in records],
               "replicates": [record["replicate"] for record in records], "metrics": result,
               "mapping_parameters": records[0]["mapping"]["parameters"],
               "deduplication_parameters": records[0]["deduplication"]["parameters"],
               "genome_filter": records[0]["mapping"]["genome_filter"],
               "note": "Pooled counts after independent library deduplication; no inference from pooled replicates"}
    write_json(f"{group}.condition.summary.json", summary)
    write_table(f"{group}.condition.qc.tsv", {group: result})
    write_table(f"{group}.replicates.qc.tsv", {
        record["sample"]: {"replicate": record["replicate"], **record["metrics"]} for record in records})
    end_tracks(counts, [record["deduplication"]["deduplicated_molecules"] for record in records],
               [tuple(row) for row in sizes[0]], group, mean=True)

    if bams is not None:
        if len(bams) != len(records):
            raise ValueError("One filtered BAM per replicate is required")
        merged = f"{group}.condition.filtered.bam"
        command = ["samtools", "merge", "-@", str(max(0, threads - 1)), "-o", merged]
        subprocess.run(command + [str(path) for path in bams], check=True)
        subprocess.run(["samtools", "index", merged], check=True)
        bam_stats(merged, f"condition__{group}.filtered")
        coverage_track(merged, f"{group}.condition.filtered", threads, bin_size)


def report_content(samples, conditions, genome=None):
    """MultiQC self-contained JSON tables with sample and pooled sections."""
    multiqc_config = {
        'module_order': [
            'custom_content',
            {'samtools': {'name': 'SingleReplicate — Samtools',
                         'anchor': 'samtools_single_replicate',
                         'path_filters': ['*/sample__*.txt']}},
            {'samtools': {'name': 'MergedReplicate — Samtools',
                         'anchor': 'samtools_merged_replicate',
                         'path_filters': ['*/condition__*.txt']}},
        ],
    }
    if genome is not None:
        from helixbusters.genomes import canonical_lengths
        canonical = canonical_lengths(genome)
        allowed = set(canonical) | {chrom[3:] for chrom in canonical}
        observed = set()
        for path in Path('.').glob('*.idxstats.txt'):
            with path.open() as handle:
                for line in handle:
                    if line.strip():
                        observed.add(line.split('\t')[0])
        multiqc_config['samtools_idxstats_ignore'] = sorted(observed - allowed)
    write_json('helixbusters_multiqc_config.json', multiqc_config)
    description = ("Mapping rate uses initial primary reads; retention uses the same denominator. "
                   "Filter exclusions are sequential and mutually exclusive, not independent overlap fractions. "
                   "UMI duplication uses reads accepted for deduplication. Coverage tracks include PCR duplicates.")
    for section, paths in (("samples", samples), ("conditions", conditions)):
        records = [json.loads(Path(path).read_text()) for path in paths]
        rows = {(record["sample"] if section == "samples" else record["group"]):
                dict(record["metrics"]) for record in records}
        if section == "samples":
            for record in records:
                rows[record["sample"]].update(group=record["group"], replicate=record["replicate"])
        write_table(f"helixbusters_{section}.tsv", rows)
        # MultiQC omits undefined numeric values; TSV/JSON retain the missingness.
        data = {name: {key: value for key, value in row.items() if value is not None}
                for name, row in rows.items()}
        write_json(f"helixbusters_{section}_mqc.json", {
            "id": f"helixbusters_{section}", "section_name": f"Helixbusters {section}",
            "description": description + (" Condition percentages are ratios of pooled counts; replicate means and sample SDs are descriptive." if section == "conditions" else ""),
            "plot_type": "table", "pconfig": {"id": f"helixbusters_{section}_table",
                                                    "title": f"Helixbusters {section} QC"}, "data": data,
        })
