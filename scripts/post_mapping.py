#!/usr/bin/env python3
"""Nextflow bridge for per-library QC, condition aggregation and MultiQC."""

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.reporting import check_environment, group_qc, report_content, sample_qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="stage", required=True)
    check = commands.add_parser("check")
    check.add_argument("--reference-config", required=True)
    check.add_argument("--genome", required=True)
    check.add_argument("--aligner", choices=("bwa", "bowtie2"), required=True)
    sample = commands.add_parser("sample")
    for option in ("sample", "group", "replicate", "mapping", "dedup", "all-bam", "filtered-bam", "counts"):
        sample.add_argument(f"--{option}", required=True)
    sample.add_argument("--threads", type=int, default=2)
    sample.add_argument("--bin-size", type=int, default=50)
    group = commands.add_parser("group")
    group.add_argument("--group", required=True)
    for option in ("summaries", "counts", "headers"):
        group.add_argument(f"--{option}", nargs="+", required=True)
    group.add_argument("--bams", nargs="+", required=True)
    group.add_argument("--threads", type=int, default=2)
    group.add_argument("--bin-size", type=int, default=50)
    report = commands.add_parser("report")
    report.add_argument("--samples", nargs="+", required=True)
    report.add_argument("--conditions", nargs="+", required=True)
    report.add_argument("--genome", help="Restrict MultiQC chromosome plots to canonical nuclear contigs")
    args = parser.parse_args()
    if args.stage in {"sample", "group"} and (args.threads < 1 or args.bin_size < 1):
        parser.error("threads and bin-size must be positive")
    if args.stage == "check":
        check_environment(args.reference_config, args.genome, args.aligner)
    elif args.stage == "sample":
        sample_qc(args.sample, args.group, args.replicate, args.mapping, args.dedup,
                  args.all_bam, args.filtered_bam, args.counts, args.threads, args.bin_size)
    elif args.stage == "group":
        group_qc(args.group, args.summaries, args.counts, args.headers, args.bams, args.threads, args.bin_size)
    else:
        report_content(args.samples, args.conditions, args.genome)


if __name__ == "__main__":
    main()
