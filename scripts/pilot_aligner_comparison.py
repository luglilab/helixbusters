#!/usr/bin/env python3
"""Compare BLISS aligners on the same reproducible prepared-read reservoir."""
import argparse
import csv
import gzip
from itertools import combinations
import json
from pathlib import Path
import random
import shutil
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.deduplication import deduplicate_bam, five_prime_position
from helixbusters.genomes import load_reference_config
from helixbusters.mapping import map_sample, validate_index
from helixbusters.reporting import write_json
from inspect_five_prime import inspect


def reservoir(source, output, limit, seed):
    """Uniform reservoir over all prepared FASTQ records, without replacement."""
    if limit < 1:
        raise ValueError('Read limit must be positive')
    rng = random.Random(seed)
    records = []
    total = 0
    opener = gzip.open if str(source).endswith('.gz') else open
    with opener(source, 'rt') as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                break
            if (not record[0].startswith('@') or not record[2].startswith('+') or
                    not record[1].strip() or len(record[1].rstrip()) != len(record[3].rstrip())):
                raise ValueError(f'Malformed FASTQ record {total + 1}')
            total += 1
            if total <= limit:
                records.append(record)
            else:
                index = rng.randrange(total)
                if index < limit:
                    records[index] = record
    if not records:
        raise ValueError('Empty prepared FASTQ')
    identifiers = [record[0].split()[0] for record in records]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError('Sampled FASTQ identifiers are not unique')
    with Path(output).open('xb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as handle:
            for record in records:
                handle.write(''.join(record).encode('ascii'))
    return {'source_reads': total, 'sampled_reads': len(records), 'seed': seed,
            'sampling': 'Uniform reservoir of prepared FASTQ records; all branches use the same file'}


def accepted_reads(path):
    """Read accepted primary endpoints from an already filtered single-end BAM."""
    result = {}
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            if read.is_paired:
                raise ValueError('Pilot requires single-end prepared libraries')
            position = five_prime_position(read, 'strict')
            if position is None:
                continue
            if read.query_name in result:
                raise ValueError('Duplicate accepted primary read identifier')
            sequence = read.get_forward_sequence() or ''
            # Flags are diagnostic only: exact motif presence does not classify
            # a molecule as an artifact and does not remove it from counting.
            result[read.query_name] = {
                'key': (read.reference_name, position, '-' if read.is_reverse else '+'),
                'mapq': read.mapping_quality, 'cigar': read.cigarstring,
                'nm': read.get_tag('NM') if read.has_tag('NM') else None,
                'technical_motif_first80': any(motif in sequence[:80] for motif in
                    ('TAATACGACTCACTATAGGG', 'CCCTATAGTGAGTCGTATTA')),
                'sequence_prefix': sequence[:50],
            }
    return result


def compare_reads(left, right):
    shared = left.keys() & right.keys()
    identical = sum(left[name]['key'] == right[name]['key'] for name in shared)
    return {'accepted_both': len(shared), 'same_five_prime_coordinate_and_strand': identical,
            'discordant_coordinate_or_strand': len(shared) - identical,
            'accepted_only_left': len(left.keys() - right.keys()),
            'accepted_only_right': len(right.keys() - left.keys()),
            'coordinate_concordance_pct_of_shared': 100 * identical / len(shared) if shared else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared-dir', type=Path, required=True)
    parser.add_argument('--samples', nargs='+', default=['HD1_ACUTE', 'HD3_CHRONIC'])
    parser.add_argument('--reference-config', required=True)
    parser.add_argument('--genome', default='hg38')
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--max-reads', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=1729)
    parser.add_argument('--map-threads', type=int, default=6)
    parser.add_argument('--sort-threads', type=int, default=1)
    parser.add_argument('--mapq', type=int, default=20)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error('Output directory exists; select a new directory')
    if args.max_reads < 1 or args.map_threads < 1 or args.sort_threads < 0 or not 0 <= args.mapq <= 254:
        parser.error('Invalid read limit, thread count or MAPQ')
    if len(set(args.samples)) != len(args.samples):
        parser.error('Sample names must be unique')
    from helixbusters.reporting import validate_label
    for sample in args.samples:
        validate_label(sample)
        if not (args.prepared_dir / f'{sample}.umi.fastq.gz').is_file():
            parser.error(f'Missing prepared FASTQ: {sample}.umi.fastq.gz')
    for tool in ('bwa', 'bowtie2', 'samtools'):
        if shutil.which(tool) is None:
            parser.error(f'Missing executable: {tool}')
    references = {aligner: load_reference_config(args.reference_config, args.genome, aligner)
                  for aligner in ('bwa', 'bowtie2')}
    for aligner, reference in references.items():
        validate_index(reference['genome_index'], aligner)
    args.outdir.mkdir(parents=True)
    report = {'parameters': {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
              'samples': {}, 'interpretation': 'Exploratory mapper comparison, not biological inference. MAPQ scores are mapper-dependent. End-to-end acceptance alone does not validate DSB coordinates. Library complexity on the sampled reads is not full-library complexity.'}
    rows = []
    for sample in args.samples:
        folder = args.outdir / sample
        folder.mkdir()
        fastq = folder / 'shared_subset.fastq.gz'
        sampling = reservoir(args.prepared_dir / f'{sample}.umi.fastq.gz', fastq, args.max_reads, args.seed)
        branches, reads, headers = {}, {}, []
        for branch, aligner, mode in [('bwa', 'bwa', 'end-to-end'),
                                      ('bowtie2_end_to_end', 'bowtie2', 'end-to-end'),
                                      ('bowtie2_local', 'bowtie2', 'local')]:
            print(f'Mapping {sample}: {branch}', flush=True)
            output = folder / branch
            reference = references[aligner]
            mapped = map_sample(sample, fastq, reference['genome_index'], output,
                aligner=aligner, bowtie2_mode=mode, min_mapq=args.mapq,
                threads=args.map_threads, sort_threads=args.sort_threads,
                genome=reference['genome'], blacklist_bed=reference['blacklist_bed'],
                blacklist_genome=reference['blacklist_genome'])
            qc = json.loads(Path(mapped['MappingQC']).read_text())
            if qc['primary_records'] != sampling['sampled_reads']:
                raise ValueError('Mapper primary-record count differs from sampled input reads')
            with pysam.AlignmentFile(mapped['BamAllPath'], 'rb') as bam:
                headers.append(list(zip(bam.references, bam.lengths)))
            if headers[-1] != headers[0]:
                raise ValueError('BWA and Bowtie2 reference headers differ; coordinate comparison is invalid')
            dedup = deduplicate_bam(mapped['BamFilteredPath'], output / f'{sample}.families.tsv',
                output / f'{sample}.counts.bed', method='directional', umi_length=8,
                min_mapq=args.mapq, five_prime_policy='strict',
                output_qc=output / f'{sample}.dedup.json')
            write_json(output / 'five_prime.json', inspect(mapped['BamFilteredPath'], args.max_reads))
            reads[branch] = accepted_reads(mapped['BamFilteredPath'])
            record = {'sample': sample, 'branch': branch, 'input_reads': sampling['sampled_reads'],
                      'mapped_primary': qc['mapped_primary_records'], 'filtered_reads': qc['retained_records'],
                      'strict_accepted_reads': dedup['accepted_reads'], 'molecules': dedup['deduplicated_molecules'],
                      'accepted_with_exact_T7_motif_first80': sum(r['technical_motif_first80'] for r in reads[branch].values())}
            rows.append(record)
            branches[branch] = {'metrics': record, 'mapping': qc, 'deduplication': dedup}
        pairs = {f'{a}_vs_{b}': compare_reads(reads[a], reads[b]) for a, b in combinations(reads, 2)}
        with (folder / 'bowtie2_recovered_reads.tsv').open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['branch', 'read', 'chrom', 'position', 'strand', 'mapq', 'cigar', 'NM', 'exact_T7_motif_first80', 'sequence_prefix'])
            for branch in ('bowtie2_end_to_end', 'bowtie2_local'):
                for name in sorted(reads[branch].keys() - reads['bwa'].keys()):
                    read = reads[branch][name]
                    writer.writerow([branch, name, *read['key'], read['mapq'], read['cigar'], read['nm'], read['technical_motif_first80'], read['sequence_prefix']])
        report['samples'][sample] = {'sampling': sampling, 'branches': branches, 'pairwise': pairs}
    with (args.outdir / 'aligner_comparison.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.outdir / 'aligner_comparison.json', report)
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
