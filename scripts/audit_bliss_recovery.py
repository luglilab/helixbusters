#!/usr/bin/env python3
"""Read-only 5-prime clipping and observed UMI-family depth diagnostics."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import random
import sys

import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.deduplication import five_prime_position
from helixbusters.technical import prefix_distance


def anchored_match(sequence, motif):
    """Longest exact read-prefix match to any contiguous motif substring."""
    best = 0
    for start in range(len(motif)):
        length = 0
        for a, b in zip(sequence, motif[start:]):
            if a != b:
                break
            length += 1
        best = max(best, length)
    return best


def reservoir_bam(path, limit, seed):
    """Uniform primary mapped read sample; scan all records once."""
    rng = random.Random(seed)
    sampled = []
    seen = 0
    with pysam.AlignmentFile(str(path), 'rb') as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            seen += 1
            if len(sampled) < limit:
                sampled.append(read)
            else:
                index = rng.randrange(seen)
                if index < limit:
                    sampled[index] = read
    return sampled, seen


def read_features(read, motif):
    sequence = read.get_forward_sequence() or ''
    qualities = read.get_forward_qualities()
    cigar = read.cigartuples or []
    operation, length = cigar[-1 if read.is_reverse else 0] if cigar else (-1, 0)
    accepted = five_prime_position(read) is not None
    reason = 'accepted' if accepted else ('spliced_N' if any(op == 3 for op, _ in cigar)
                                           else f'five_prime_{"MIDNSHP=X"[operation] if 0 <= operation <= 8 else "missing"}')
    clipped = sequence[:length] if operation == 4 else ''
    observed_first10 = float(np.mean(qualities[:10])) if qualities is not None and len(qualities) else None
    clip_quality = float(np.mean(qualities[:length])) if operation == 4 and qualities is not None and length else None
    return {'read_name': read.query_name, 'reason': reason, 'strand': '-' if read.is_reverse else '+',
            'clip_length': length if operation in (4, 5) else 0,
            'soft_clip_length': length if operation == 4 else 0,
            'read_length': len(sequence), 'mapq': read.mapping_quality,
            'first10_mean_quality': observed_first10, 'clip_mean_quality': clip_quality,
            'first10_has_N': 'N' in sequence[:10],
            'prefix_motif_exact_bases': anchored_match(sequence, motif),
            'clip_motif_exact_bases': anchored_match(clipped, motif),
            'full_prefix_motif_within2_edits': prefix_distance(sequence, motif) <= 2,
            'sequence_prefix': sequence[:40], 'clip_prefix': clipped[:40]}


def observed_complexity(path, accepted, molecules):
    """Bernoulli thinning of fixed observed families; no unseen-library estimate."""
    sizes = Counter()
    with Path(path).open() as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if len(row) != 6 or int(row[5]) < 1:
                raise ValueError('Invalid molecular family table')
            sizes[int(row[5])] += 1
    if sum(sizes.values()) != molecules or sum(k * v for k, v in sizes.items()) != accepted:
        raise ValueError('Family sizes do not conserve accepted reads and deduplicated molecules')
    curve = [{'fraction': float(p), 'expected_accepted_reads': float(p * accepted),
              'expected_observed_families': float(sum(n * (1 - (1 - p) ** k) for k, n in sizes.items()))}
             for p in np.linspace(0, 1, 21)]
    return sizes, curve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--single-replicate-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--samples', nargs='+')
    parser.add_argument('--max-reads', type=int, default=50000)
    parser.add_argument('--seed', type=int, default=1729)
    parser.add_argument('--motif', default='CCCTATAGTGAGTCGTAT')
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error('Choose a new output directory; original results are protected')
    if args.max_reads < 1 or len(args.motif) < 12 or set(args.motif) - set('ACGT'):
        parser.error('Positive sample size and A/C/G/T motif of at least 12 bases required')
    folders = sorted(p for p in args.single_replicate_dir.iterdir() if p.is_dir()
                     and (args.samples is None or p.name in args.samples))
    if not folders or args.samples is not None and {p.name for p in folders} != set(args.samples):
        parser.error('Requested sample directories are missing')
    args.outdir.mkdir(parents=True)
    summaries, groups, curves, full = [], [], [], {}
    for folder in folders:
        sample = folder.name
        print(f'Auditing {sample}: one BAM scan; <= {args.max_reads:,} uniformly sampled reads', flush=True)
        qc = json.loads((folder / 'qc' / f'{sample}.summary.json').read_text())
        dedup = json.loads((folder / 'deduplication' / f'{sample}.dedup.json').read_text())
        paths = list((folder / 'mapping').glob(f'{sample}.q*.bam'))
        if len(paths) != 1:
            raise ValueError(f'Expected one filtered BAM for {sample}')
        reads, seen = reservoir_bam(paths[0], args.max_reads, args.seed)
        if seen != dedup['input_reads']:
            raise ValueError('BAM and deduplication input counts differ; use matching outputs')
        table = pd.DataFrame([read_features(read, args.motif) for read in reads])
        if table.empty:
            raise ValueError('No primary mapped alignments to diagnose')
        table.to_csv(args.outdir / f'{sample}.sampled_reads.tsv', sep='\t', index=False)
        hist = Counter(table.loc[table.reason == 'five_prime_S', 'soft_clip_length'].astype(int))
        for reason, subset in table.groupby('reason'):
            first_quality = subset.first10_mean_quality.dropna()
            clip_quality = subset.clip_mean_quality.dropna()
            groups.append({'sample': sample, 'reason': reason, 'sampled_reads': len(subset),
                           'reads_with_base_quality': len(first_quality),
                           'median_first10_quality': first_quality.median() if len(first_quality) else None,
                           'median_clip_quality': clip_quality.median() if len(clip_quality) else None,
                           'first10_mean_quality_below20_pct': 100 * first_quality.lt(20).mean() if len(first_quality) else None,
                           'first10_has_N_pct': 100 * subset.first10_has_N.mean(),
                           'prefix_exact_motif_segment_ge8_pct': 100 * subset.prefix_motif_exact_bases.ge(8).mean(),
                           'clip_exact_motif_segment_ge8_pct': 100 * subset.clip_motif_exact_bases.ge(8).mean(),
                           'full_prefix_motif_within2_edits_pct': 100 * subset.full_prefix_motif_within2_edits.mean()})
        families, curve = observed_complexity(folder / 'deduplication' / f'{sample}.families.tsv',
                                               dedup['accepted_reads'], dedup['deduplicated_molecules'])
        curves.extend({'sample': sample, **row} for row in curve)
        accepted, molecules = dedup['accepted_reads'], dedup['deduplicated_molecules']
        if accepted < 1 or molecules < 1:
            raise ValueError('No strictly accepted molecules for depth diagnostics')
        record = {'sample': sample, 'filtered_reads': seen, 'accepted_reads': accepted,
                  'molecules': molecules, 'ambiguous_reads': dedup['skipped'].get('ambiguous_five_prime', 0),
                  'accepted_pct_of_filtered': 100 * accepted / seen,
                  'molecules_pct_of_filtered': 100 * molecules / seen,
                  'duplicate_pct_of_accepted': 100 * (accepted - molecules) / accepted,
                  'singleton_family_pct': 100 * families[1] / molecules,
                  'sampled_reads': len(table), 'sampled_soft_clip_pct': 100 * table.reason.eq('five_prime_S').mean(),
                  'soft_clip_median_bases': float(table.loc[table.reason == 'five_prime_S', 'soft_clip_length'].median()) if hist else None,
                  'soft_clip_1_to5_pct_of_soft_clips': float(100 * table.loc[table.reason == 'five_prime_S', 'soft_clip_length'].le(5).mean()) if hist else None}
        summaries.append(record)
        full[sample] = {'summary': record, 'clip_length_histogram': {str(k): v for k, v in sorted(hist.items())},
                        'top_clip_prefixes': Counter(table.loc[table.reason == 'five_prime_S', 'clip_prefix']).most_common(30),
                        'family_size_histogram': {str(k): v for k, v in sorted(families.items())},
                        'mapping_parameters': qc.get('mapping', {}).get('parameters', {})}
    pd.DataFrame(summaries).to_csv(args.outdir / 'yield.samples.tsv', sep='\t', index=False)
    pd.DataFrame(groups).to_csv(args.outdir / 'clipping.groups.tsv', sep='\t', index=False)
    pd.DataFrame(curves).to_csv(args.outdir / 'observed_depth_curves.tsv', sep='\t', index=False)
    (args.outdir / 'audit.json').write_text(json.dumps({'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        'sampling': 'Uniform reservoir over all primary mapped reads in each filtered BAM; first N records are not used. Soft-clipped bases and qualities are restored to original sequencing orientation; hard-clipped bases are absent.',
        'motif_matching': 'Heuristic exact read-prefix match to any contiguous substring of the supplied motif; >=8 bases is diagnostic, not evidence of an insert boundary. Accepted reads are a control.',
        'complexity': 'Expected observed fixed families under Bernoulli thinning of accepted reads. Directional clusters are held fixed, not re-inferred; no extrapolation beyond observed depth.',
        'versions': {'pysam': pysam.__version__, 'numpy': np.__version__, 'pandas': pd.__version__},
        'samples': full}, indent=2) + '\n')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    data = pd.DataFrame(summaries)
    x = np.arange(len(data))
    axes[0].bar(x, data.filtered_reads, label='Filtered mapping reads', color='#dddddd')
    axes[0].bar(x, data.accepted_reads, label='Strict accepted reads', color='#edac58')
    axes[0].bar(x, data.molecules, label='Deduplicated molecules', color='#527a9f')
    axes[0].set_xticks(x, data['sample'], rotation=45, ha='right'); axes[0].set_ylabel('Records / molecules'); axes[0].legend()
    for sample, rows in pd.DataFrame(curves).groupby('sample'):
        axes[1].plot(rows.expected_accepted_reads, rows.expected_observed_families, label=sample)
    axes[1].set(xlabel='Expected accepted reads under thinning', ylabel='Expected observed fixed UMI families', title='Within observed depth only')
    axes[1].legend(fontsize=8); fig.tight_layout(); fig.savefig(args.outdir / 'usable_depth.pdf'); fig.savefig(args.outdir / 'usable_depth.png', dpi=160); plt.close(fig)
    print(args.outdir)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
