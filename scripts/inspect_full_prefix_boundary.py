#!/usr/bin/env python3
"""Diagnose residual T7 bases and compare exact 18/20-base prefix removal."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.deduplication import five_prime_position, deduplicate_bam
from helixbusters.genomes import load_reference_config, GenomeFilter
from helixbusters.mapping import map_sample, validate_index
from helixbusters.technical import T7_FORWARD, T7_REVERSE, full_t7_rescue_reason
from pilot_aligner_comparison import accepted_reads, compare_reads

def end_mismatches(read, length=10):
    """MD-backed reference comparison in biological 5-prime orientation."""
    if not read.has_tag('MD'):
        return [None] * length
    aligned = {q: base for q, r, base in read.get_aligned_pairs(with_seq=True)
               if q is not None and r is not None and base is not None}
    sequence = read.query_sequence or ''
    result = []
    for i in range(length):
        q = len(sequence) - 1 - i if read.is_reverse else i
        result.append(sequence[q].upper() != aligned[q].upper()
                      if 0 <= q < len(sequence) and q in aligned else None)
    return result


def diagnose_bam(path, sample, branch, writer):
    groups = {}
    prefixes = {}
    with pysam.AlignmentFile(path, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_paired:
                raise ValueError('Single-end BLISS reads required')
            strict = five_prime_position(read, 'strict')
            group = 'strict_accepted' if strict is not None else 'ambiguous_5prime'
            counts = groups.setdefault(group, Counter())
            sequence = read.get_forward_sequence() or ''
            counts['reads'] += 1
            counts['starts_TA'] += sequence.startswith('TA')
            counts['forward_T7_prefix12_in_first60'] += T7_FORWARD[:12] in sequence[:60]
            counts['reverse_T7_prefix12_in_first60'] += T7_REVERSE[:12] in sequence[:60]
            prefixes.setdefault(group, Counter())[sequence[:30]] += 1
            if strict is None:
                continue
            mismatch = end_mismatches(read)
            for i, value in enumerate(mismatch):
                if value is not None:
                    counts[f'position_{i + 1}_examined'] += 1
                    counts[f'position_{i + 1}_mismatch'] += value
            counts['first_two_any_mismatch'] += any(x is True for x in mismatch[:2])
            writer.writerow([sample, branch, read.query_name, read.reference_name, strict,
                             '-' if read.is_reverse else '+', read.mapping_quality,
                             read.get_tag('NM') if read.has_tag('NM') else None,
                             sequence[:40], *mismatch])
    return {group: {'counts': dict(counts), 'top_prefixes': prefixes[group].most_common(20)}
            for group, counts in groups.items()}


def prepare_comparison(source, outdir):
    """Exact full reverse T7 only; preserve names/UMIs, exclude both T7 repeats."""
    counts = Counter()
    paths = {branch: outdir / f'{branch}.fastq.gz' for branch in ('trim18', 'trim20')}
    handles = {branch: gzip.open(path, 'xt') for branch, path in paths.items()}
    names = set()
    try:
        with gzip.open(source, 'rt', encoding='ascii') as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip('\r\n').upper()
                plus = handle.readline(); quality = handle.readline().rstrip('\r\n')
                if not header.startswith('@') or not plus.startswith('+') or not sequence or len(sequence) != len(quality) or set(sequence) - set('ACGTN'):
                    raise ValueError('Malformed prefix-pilot FASTQ')
                name = header.split()[0]
                if name in names:
                    raise ValueError('Duplicate read identifier')
                names.add(name)
                umi = name.rsplit('_', 1)[-1]
                if len(umi) != 8 or set(umi) - set('ACGT'):
                    raise ValueError('Expected authoritative 8-base UMI suffix')
                counts['input_reads'] += 1
                reason = full_t7_rescue_reason(sequence, quality)
                if reason:
                    counts[reason] += 1
                else:
                    counts['comparison_reads'] += 1
                    for branch, length in (('trim18', 18), ('trim20', 20)):
                        handles[branch].write(f'{name}\n{sequence[length:]}\n+\n{quality[length:]}\n')
    finally:
        for handle in handles.values():
            handle.close()
    if sum(value for key, value in counts.items() if key != 'input_reads') != counts['input_reads']:
        raise ValueError('Comparison cohort classification does not conserve reads')
    return dict(counts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit-dir', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--reference-config', type=Path, help='Enable matched BWA 18/20-base trimming pilot')
    p.add_argument('--map-threads', type=int, default=6)
    p.add_argument('--sort-threads', type=int, default=1)
    args = p.parse_args()
    if args.outdir.exists() or args.map_threads < 1 or args.sort_threads < 0:
        p.error('Use a new output directory and valid CPU limits')
    source = args.audit_dir.resolve()
    audit = json.loads((source / 'audit.json').read_text())
    if audit['parameters']['motif'] != T7_REVERSE[:18]:
        p.error('This comparison requires the original 18-base T7 audit')
    reference = None
    if args.reference_config:
        reference = load_reference_config(args.reference_config, audit['parameters']['genome'], 'bwa')
        validate_index(reference['genome_index'], 'bwa')
        if reference['genome']:
            masking = GenomeFilter(reference['genome'], reference['blacklist_bed'], reference['blacklist_genome'])
            for sample, item in audit['samples'].items():
                original = item['mapping']['excluded_exact_trim']['mapping']['genome_filter']
                if original['blacklist_sha256'] != masking.sha256:
                    p.error(f'Blacklist changed since the original audit: {sample}')
    args.outdir.mkdir(parents=True)
    report = {'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              'technical_sequences': {'reverse_T7_20': T7_REVERSE, 'forward_T7_20': T7_FORWARD},
              'reference': reference, 'samples': {},
              'interpretation': 'Diagnostic and experimental matched-cohort comparison only. Exact prefix evidence does not independently validate physical DSB identity or the true insert boundary. Both T7 orientations are screened in the first 60 tail bases; other technical structures can remain. MD comparison is strand-aware; unmapped query bases and missing MD are not counted as matches. No production settings changed.'}
    summary = []
    with (args.outdir / 'strict_reads.boundary.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'branch', 'read', 'chrom', 'endpoint', 'strand', 'mapq', 'NM', 'sequence_prefix',
                         *[f'mismatch_position_{i}' for i in range(1, 11)]])
        for sample in audit['samples']:
            print(f'{sample}: diagnosing residual T7 sequence', flush=True)
            folder = args.outdir / sample; folder.mkdir()
            item = {'diagnostics': {}, 'mapping': {}}
            for branch in ('production_retained', 'excluded_exact_trim'):
                path = source / sample / branch / f'{sample}.q20.bam'
                if path.exists():
                    item['diagnostics'][branch] = diagnose_bam(path, sample, branch, writer)
            item['cohort'] = prepare_comparison(source / sample / 'excluded_exact_untrimmed.fastq.gz', folder)
            size = item['cohort'].get('comparison_reads', 0)
            endpoints = {}
            for branch in ('trim18', 'trim20'):
                if not reference or not size:
                    continue
                print(f'{sample}: mapping {branch}, {size:,} matched reads', flush=True)
                dest = folder / branch
                mapped = map_sample(sample, folder / f'{branch}.fastq.gz', reference['genome_index'], dest,
                                    aligner='bwa', min_mapq=20, threads=args.map_threads, sort_threads=args.sort_threads,
                                    genome=reference['genome'], blacklist_bed=reference['blacklist_bed'], blacklist_genome=reference['blacklist_genome'])
                qc = json.loads(Path(mapped['MappingQC']).read_text())
                if qc['primary_records'] != size:
                    raise ValueError('Mapping input differs from matched cohort')
                dedup = deduplicate_bam(mapped['BamFilteredPath'], dest / f'{sample}.families.tsv', dest / f'{sample}.counts.bed',
                                        method='directional', umi_length=8, min_mapq=20, five_prime_policy='strict', output_qc=dest / f'{sample}.dedup.json')
                item['mapping'][branch] = {'mapping': qc, 'deduplication': dedup,
                                          'diagnostics': diagnose_bam(Path(mapped['BamFilteredPath']), sample, branch, writer)}
                endpoints[branch] = accepted_reads(mapped['BamFilteredPath'])
                summary.append({'sample': sample, 'branch': branch, 'input_reads': size,
                                'accepted_reads': dedup['accepted_reads'], 'molecules': dedup['deduplicated_molecules']})
            if endpoints:
                item['endpoint_comparison'] = compare_reads(endpoints['trim18'], endpoints['trim20'])
                with (folder / 'shared_endpoints.tsv').open('x') as comparison:
                    w = csv.writer(comparison, delimiter='\t')
                    w.writerow(['read', 'chrom18', 'endpoint18', 'strand18', 'chrom20', 'endpoint20', 'strand20', 'delta20_minus18'])
                    for name in sorted(endpoints['trim18'].keys() & endpoints['trim20'].keys()):
                        a, b = endpoints['trim18'][name]['key'], endpoints['trim20'][name]['key']
                        delta = b[1] - a[1] if a[0] == b[0] and a[2] == b[2] else None
                        w.writerow([name, *a, *b, delta])
            report['samples'][sample] = item
            (args.outdir / 'boundary.json').write_text(json.dumps(report, indent=2) + '\n')
    if summary:
        with (args.outdir / 'boundary.metrics.tsv').open('x') as handle:
            w = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter='\t'); w.writeheader(); w.writerows(summary)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, EOFError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
