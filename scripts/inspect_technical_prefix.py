#!/usr/bin/env python3
"""Count anchored technical-prefix variants in prepared FASTQ, without trimming."""
import argparse
from collections import Counter
import json
from pathlib import Path

from inspect_fastq_layout import fastq_sequences


def prefix_distance(sequence, motif, maximum=2):
    """Bounded Levenshtein distance to a prefix, with no free leading bases."""
    # At most maximum insertions can extend a matching prefix.
    sequence = sequence[:len(motif) + maximum]
    previous = list(range(len(sequence) + 1))
    for i, base in enumerate(motif, 1):
        current = [i]
        for j, observed in enumerate(sequence, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1,
                               previous[j - 1] + (base != observed)))
        previous = current
    candidates = previous[max(0, len(motif) - maximum):]
    return min(candidates, default=maximum + 1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fastq', type=Path, required=True, help='Prepared FASTQ after UMI/barcode removal')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-reads', type=int, default=100000)
    p.add_argument('--motif', default='CCCTATAGTGAGTCGTAT')
    args = p.parse_args()
    if args.max_reads < 1 or len(args.motif) < 12 or set(args.motif) - set('ACGT'):
        p.error('Positive max-reads and an A/C/G/T motif of at least 12 bases are required')
    if args.output.exists():
        p.error('Output exists; choose a new filename')
    bins = Counter({str(i): 0 for i in range(3)})
    bins['greater_than_2'] = 0
    examples = {str(i): Counter() for i in range(3)}
    examined = 0
    for seq in fastq_sequences(args.fastq, args.max_reads):
        examined += 1
        distance = prefix_distance(seq, args.motif)
        category = str(distance) if distance <= 2 else 'greater_than_2'
        bins[category] += 1
        if distance <= 2:
            examples[category][seq[:30]] += 1
    report = {
        'fastq': str(args.fastq), 'motif': args.motif, 'reads_examined': examined,
        'sampling': 'First N prepared FASTQ records; not a random sample',
        'method': 'Anchored Levenshtein distance to full motif; substitutions, insertions and deletions cost one each',
        'exclusive_distance_counts': dict(bins),
        'cumulative_matching_pct': {str(i): 100 * sum(bins[str(j)] for j in range(i + 1)) / examined
                                    if examined else None for i in range(3)},
        'frequent_matching_prefixes': {key: value.most_common(20) for key, value in examples.items()},
        'interpretation': 'Motif similarity does not establish adapter-dimer identity or genomic insert boundaries. Current prepared FASTQ already excludes exact motif matches.',
    }
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps({key: report[key] for key in ('reads_examined', 'exclusive_distance_counts', 'cumulative_matching_pct')}, indent=2))


if __name__ == '__main__':
    main()
