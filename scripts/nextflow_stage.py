#!/usr/bin/env python3
"""Small command-line bridge from Nextflow tasks to Helixbusters APIs."""

import argparse
import json
from pathlib import Path
import sys

# Nextflow runs this file from an isolated task directory. Make the project
# package importable from the checkout even when it was not pip-installed.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from helixbusters.deduplication import deduplicate_bam
from helixbusters.genomes import load_reference_config
from helixbusters.mapping import map_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="stage", required=True)
    mapping = commands.add_parser("map")
    mapping.add_argument("--sample", required=True)
    mapping.add_argument("--read1", required=True)
    mapping.add_argument("--genome", required=True)
    mapping.add_argument("--reference-config", required=True)
    mapping.add_argument("--aligner", choices=("bwa", "bowtie2"), default="bwa")
    mapping.add_argument("--mapq", type=int, default=20)
    mapping.add_argument("--threads", type=int, default=8)
    mapping.add_argument("--sort-threads", type=int, default=1)
    mapping.add_argument("--sort-memory", default="768M")
    mapping.add_argument("--outdir", required=True)

    dedup = commands.add_parser("dedup")
    dedup.add_argument("--sample", required=True)
    dedup.add_argument("--bam", required=True)
    dedup.add_argument("--method", choices=("exact", "directional"), default="directional")
    dedup.add_argument("--umi-length", type=int, default=8)
    dedup.add_argument("--mapq", type=int, default=20)
    dedup.add_argument("--outdir", required=True)

    args = parser.parse_args()
    output = Path(args.outdir)
    output.mkdir(parents=True, exist_ok=True)
    if args.stage == "map":
        reference = load_reference_config(args.reference_config, args.genome, args.aligner)
        result = map_sample(
            args.sample, args.read1, reference["genome_index"], output,
            aligner=args.aligner, min_mapq=args.mapq, threads=args.threads,
            sort_threads=args.sort_threads, sort_memory=args.sort_memory,
            genome=reference["genome"], blacklist_bed=reference["blacklist_bed"],
            blacklist_genome=reference["blacklist_genome"],
        )
        print(json.dumps(result, indent=2))
    else:
        prefix = output / args.sample
        result = deduplicate_bam(
            args.bam, Path(f"{prefix}.families.tsv"), Path(f"{prefix}.counts.bed"),
            method=args.method, umi_length=args.umi_length, min_mapq=args.mapq,
            read_selection="single-end", output_molecules=Path(f"{prefix}.molecules.bed"),
            output_sites=Path(f"{prefix}.sites.tsv"), output_qc=Path(f"{prefix}.dedup.json"),
        )
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
