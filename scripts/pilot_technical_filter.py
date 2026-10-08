#!/usr/bin/env python3
"""Compare identical prepared-read subsets before and after motif exclusion."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import subprocess
import sys

from inspect_technical_prefix import prefix_distance
from inspect_five_prime import inspect

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.genomes import load_reference_config
from helixbusters.mapping import map_sample, validate_index
from helixbusters.deduplication import deduplicate_bam


def subsets(source, output, limit, motif):
    """Write paired subsets without changing identifiers or qualities."""
    baseline, filtered = output / 'baseline.fastq.gz', output / 'exclude_2_errors.fastq.gz'
    counts = Counter({'examined': 0, 'excluded': 0, 'retained': 0})
    with gzip.open(source, 'rt') as inp, gzip.open(baseline, 'xt') as a, gzip.open(filtered, 'xt') as b:
        for _ in range(limit):
            record = [inp.readline() for _ in range(4)]
            if not record[0]:
                break
            seq = record[1].rstrip('\r\n').upper()
            if (not record[0].startswith('@') or not record[2].startswith('+') or
                    not seq or len(seq) != len(record[3].rstrip('\r\n'))):
                raise ValueError('Malformed prepared FASTQ')
            counts['examined'] += 1
            a.writelines(record)
            if prefix_distance(seq, motif) <= 2:
                counts['excluded'] += 1
            else:
                b.writelines(record)
                counts['retained'] += 1
    if not counts['retained']:
        raise ValueError('No reads retained in pilot')
    return dict(counts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fastq', type=Path, required=True, help='Existing prepared .umi.fastq.gz')
    p.add_argument('--reference-config', required=True)
    p.add_argument('--genome', default='hg38')
    p.add_argument('--outdir', type=Path, required=True, help='New directory; existing directories are refused')
    p.add_argument('--max-reads', type=int, default=100000)
    p.add_argument('--map-threads', type=int, default=2)
    p.add_argument('--sort-threads', type=int, default=1)
    p.add_argument('--umi-length', type=int, default=8)
    p.add_argument('--motif', default='CCCTATAGTGAGTCGTAT')
    args = p.parse_args()
    if min(args.max_reads, args.map_threads, args.umi_length) < 1 or args.sort_threads < 0:
        p.error('Invalid read limit or thread/UMI count')
    if len(args.motif) < 12 or set(args.motif) - set('ACGT'):
        p.error('Motif must contain at least 12 A/C/G/T bases')
    if args.outdir.exists():
        p.error('Output directory already exists; choose a new one')
    reference = load_reference_config(args.reference_config, args.genome, 'bwa')
    validate_index(reference['genome_index'], 'bwa')
    args.outdir.mkdir(parents=True)
    counts = subsets(args.fastq, args.outdir, args.max_reads, args.motif)
    rows = []
    for branch in ('baseline', 'exclude_2_errors'):
        print(f'Mapping pilot branch: {branch}', flush=True)
        folder = args.outdir / branch
        folder.mkdir()
        mapping = map_sample(branch, args.outdir / f'{branch}.fastq.gz', reference['genome_index'], folder,
                             aligner='bwa', min_mapq=20, threads=args.map_threads,
                             sort_threads=args.sort_threads, genome=reference['genome'],
                             blacklist_bed=reference['blacklist_bed'], blacklist_genome=reference['blacklist_genome'])
        mapped_qc = json.loads(Path(mapping['MappingQC']).read_text())
        bam = mapping['BamFilteredPath']
        diagnosis = inspect(bam, args.max_reads)
        (folder / 'five_prime.json').write_text(json.dumps(diagnosis, indent=2) + '\n')
        dedup = deduplicate_bam(bam, folder / f'{branch}.families.tsv', folder / f'{branch}.counts.bed',
                               method='directional', umi_length=args.umi_length, min_mapq=20,
                               five_prime_policy='strict', output_qc=folder / f'{branch}.dedup.json')
        primary = mapped_qc['primary_records']
        retained = mapped_qc['retained_records']
        rows.append({'branch': branch, 'input_subset_reads': primary,
                     'mapped_primary_reads': mapped_qc['mapped_primary_records'],
                     'mapping_pct': 100 * mapped_qc['mapped_primary_records'] / primary if primary else None,
                     'filtered_reads': retained, 'accepted_strict_reads': dedup['accepted_reads'],
                     'ambiguous_five_prime_reads': dedup['skipped'].get('ambiguous_five_prime', 0),
                     'strict_acceptance_pct': 100 * dedup['accepted_reads'] / retained if retained else None,
                     'molecules': dedup['deduplicated_molecules'],
                     'molecules_per_100k_baseline_reads': dedup['deduplicated_molecules'] * 100000 / counts['examined']})
    report = {'subset': counts, 'comparisons': rows, 'reference': reference,
              'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              'python': sys.version, 'samtools': subprocess.check_output(['samtools', '--version'], text=True),
              'interpretation': 'First N prepared reads, identical starting subset. Exclusion only; no clipping rescue. Retention percentages and absolute molecule yield must both be considered. Lower clipping alone does not validate DSB specificity.'}
    with (args.outdir / 'pilot_comparison.json').open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(report['comparisons'], indent=2))


if __name__ == '__main__':
    main()
