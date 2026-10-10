#!/usr/bin/env python3
"""Experimental partial-prefix correction on identical filtered-BAM read cohorts."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.deduplication import _read_umi, deduplicate_bam
from helixbusters.genomes import load_reference_config
from helixbusters.mapping import map_sample, validate_index
from audit_bliss_recovery import reservoir_bam
from pilot_aligner_comparison import accepted_reads, compare_reads


def candidate_clip(read, motif):
    """Conservative pilot rule; no genomic boundary validation is implied."""
    cigar = read.cigartuples or []
    if not cigar or read.is_unmapped or any(op == 3 for op, _ in cigar):
        return 0
    operation, length = cigar[-1 if read.is_reverse else 0]
    if operation != 4 or not 12 <= length < len(motif):
        return 0
    sequence = read.get_forward_sequence() or ''
    quality = read.get_forward_qualities()
    if sequence[:length] != motif[:length] or len(sequence) - length < 40 or quality is None:
        return 0
    if np.mean(quality[:length]) < 30 or min(quality[length:length + 5]) < 20:
        return 0
    return length


def write_cohort(reads, folder, motif):
    """Preserve original orientation and validated UMI; alter candidates only."""
    changed, names = {}, set()
    with gzip.open(folder / 'baseline.fastq.gz', 'xt') as baseline, gzip.open(folder / 'candidate_trim.fastq.gz', 'xt') as corrected:
        for read in reads:
            if read.is_paired:
                raise ValueError('This pilot requires single-end libraries')
            umi, reason = _read_umi(read, 8)
            if reason:
                raise ValueError('Sample contains missing/invalid UMIs; correction cannot safely proceed')
            name = read.query_name
            if name.rsplit('_', 1)[-1] != umi:
                name += '_' + umi  # Preserve authoritative RX when a suffix is absent/different.
            if name in names:
                raise ValueError('Duplicate read identifiers in pilot cohort')
            names.add(name)
            sequence, quality = read.get_forward_sequence(), read.get_forward_qualities()
            if not sequence or quality is None:
                raise ValueError('Pilot FASTQ reconstruction needs sequence and base qualities')
            qualities = ''.join(chr(q + 33) for q in quality)
            trim = candidate_clip(read, motif)
            baseline.write(f'@{name}\n{sequence}\n+\n{qualities}\n')
            corrected.write(f'@{name}\n{sequence[trim:]}\n+\n{qualities[trim:]}\n')
            if trim:
                changed[name] = {'removed_bases': trim, 'source_aligned_endpoint': (
                    read.reference_name, read.reference_end - 1 if read.is_reverse else read.reference_start,
                    '-' if read.is_reverse else '+')}
    return changed


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--single-replicate-dir', type=Path, required=True)
    p.add_argument('--samples', nargs='+', required=True)
    p.add_argument('--reference-config', type=Path, required=True)
    p.add_argument('--genome', default='hg38')
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--max-reads', type=int, default=50000)
    p.add_argument('--seed', type=int, default=1729)
    p.add_argument('--motif', default='CCCTATAGTGAGTCGTAT')
    p.add_argument('--map-threads', type=int, default=6)
    p.add_argument('--sort-threads', type=int, default=1)
    args = p.parse_args()
    if args.outdir.exists() or args.max_reads < 1 or len(set(args.samples)) != len(args.samples):
        p.error('Use a new output directory, positive cohort size and unique samples')
    if len(args.motif) < 13 or set(args.motif) - set('ACGT'):
        p.error('Motif must contain at least 13 A/C/G/T bases')
    reference = load_reference_config(args.reference_config, args.genome, 'bwa')
    validate_index(reference['genome_index'], 'bwa')
    args.outdir.mkdir(parents=True)
    report = {'sampling': 'Uniform primary-mapped read reservoir from an already filtered BAM; conditional on original mapping, not a whole-library yield estimate',
              'candidate_rule': '5-prime soft clip of 12..motif_length-1 bases exactly matching the motif prefix; clip mean Q>=30; next five bases Q>=20; remaining insert >=40; no spliced CIGAR',
              'interpretation': 'Pilot only. Remapping concordance does not independently establish true DSB coordinates. Primary pipeline is unchanged. UMI and strict endpoint policies remain enabled.',
              'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}, 'samples': {}}
    rows = []
    for sample in args.samples:
        source = args.single_replicate_dir / sample
        bams = list((source / 'mapping').glob(f'{sample}.q*.bam'))
        if len(bams) != 1:
            raise ValueError(f'Expected one filtered BAM for {sample}')
        reads, seen = reservoir_bam(bams[0], args.max_reads, args.seed)
        if not reads:
            raise ValueError('Empty pilot cohort')
        folder = args.outdir / sample; folder.mkdir()
        changed = write_cohort(reads, folder, args.motif)
        results, endpoints = {}, {}
        for branch in ('baseline', 'candidate_trim'):
            print(f'{sample}: {branch}; {len(reads):,} reads, {len(changed):,} correction candidates', flush=True)
            output = folder / branch
            mapped = map_sample(sample, folder / f'{branch}.fastq.gz', reference['genome_index'], output,
                                aligner='bwa', min_mapq=20, threads=args.map_threads, sort_threads=args.sort_threads,
                                genome=reference['genome'], blacklist_bed=reference['blacklist_bed'], blacklist_genome=reference['blacklist_genome'])
            mapping_qc = json.loads(Path(mapped['MappingQC']).read_text())
            if mapping_qc['primary_records'] != len(reads):
                raise ValueError('Remapping primary count differs from the shared cohort size')
            dedup = deduplicate_bam(mapped['BamFilteredPath'], output / f'{sample}.families.tsv', output / f'{sample}.counts.bed',
                                   method='directional', min_mapq=20, five_prime_policy='strict',
                                   output_qc=output / f'{sample}.dedup.json')
            endpoints[branch] = accepted_reads(mapped['BamFilteredPath'])
            results[branch] = dedup
            rows.append({'sample': sample, 'branch': branch, 'cohort_reads': len(reads), 'candidate_reads': len(changed),
                         'strict_accepted_reads': dedup['accepted_reads'], 'molecules': dedup['deduplicated_molecules']})
        shared = compare_reads(endpoints['baseline'], endpoints['candidate_trim'])
        newly_accepted = endpoints['candidate_trim'].keys() - endpoints['baseline'].keys()
        unchanged_controls = endpoints['baseline'].keys() - changed.keys()
        control_agreement = sum(name in endpoints['candidate_trim'] and endpoints['baseline'][name]['key'] == endpoints['candidate_trim'][name]['key']
                                for name in unchanged_controls)
        source_agreement = sum(name in changed and tuple(changed[name]['source_aligned_endpoint']) == endpoints['candidate_trim'][name]['key']
                               for name in newly_accepted)
        with (folder / 'newly_accepted_reads.tsv').open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['read', 'correction_candidate', 'removed_bases', 'chrom', 'endpoint', 'strand', 'mapq', 'NM', 'source_aligned_endpoint_agrees'])
            for name in sorted(newly_accepted):
                record = endpoints['candidate_trim'][name]
                writer.writerow([name, name in changed, changed.get(name, {}).get('removed_bases', 0), *record['key'], record['mapq'], record['nm'],
                                 name in changed and tuple(changed[name]['source_aligned_endpoint']) == record['key']])
        report['samples'][sample] = {'source_filtered_reads': seen, 'cohort_reads': len(reads), 'candidate_reads': len(changed),
                                     'baseline': results['baseline'], 'candidate_trim': results['candidate_trim'],
                                     'shared_endpoint_comparison': shared, 'newly_strict_accepted': len(newly_accepted),
                                     'new_endpoint_matches_source_aligned_boundary': source_agreement,
                                     'unchanged_accepted_controls': len(unchanged_controls), 'unchanged_controls_same_endpoint': control_agreement}
    (args.outdir / 'pilot.json').write_text(json.dumps(report, indent=2) + '\n')
    with (args.outdir / 'pilot.metrics.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t'); writer.writeheader(); writer.writerows(rows)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
