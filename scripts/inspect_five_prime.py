#!/usr/bin/env python3
"""Diagnose biological 5-prime CIGAR ends in a bounded, read-only BAM scan."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.deduplication import five_prime_position


def inspect(bam, limit):
    reasons, operations, lengths, clips, starts, cigars, strands = (Counter() for _ in range(7))
    scanned = examined = 0
    with pysam.AlignmentFile(str(bam), 'rb') as handle:
        for read in handle.fetch(until_eof=True):
            scanned += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            examined += 1
            strand = 'reverse' if read.is_reverse else 'forward'
            strands[strand] += 1
            cigar = read.cigartuples or []
            operation, length = cigar[-1 if read.is_reverse else 0] if cigar else (-1, 0)
            label = 'MIDNSHP=X'[operation] if 0 <= operation <= 8 else 'missing'
            operations[f'{strand}:{label}'] += 1
            accepted = five_prime_position(read) is not None
            reason = 'accepted_strict' if accepted else ('spliced_N' if any(op == 3 for op, _ in cigar) else f'five_prime_{label}')
            reasons[reason] += 1
            if not accepted:
                cigars[read.cigarstring or 'missing'] += 1
                sequence = read.get_forward_sequence() or ''
                starts[sequence[:30]] += 1
                if operation in {4, 5}:
                    lengths[f'{label}:{length}'] += 1
                if operation == 4 and sequence:
                    clips[sequence[:length][:30]] += 1
            if examined >= limit:
                break
    return {
        'bam': str(bam), 'sampling': 'First N primary mapped alignments in BAM order; not a random sample',
        'limit': limit, 'records_scanned': scanned, 'primary_mapped_examined': examined,
        'exclusive_reasons': dict(reasons),
        'ambiguous_pct': 100 * (examined - reasons['accepted_strict']) / examined if examined else None,
        'strand_counts': dict(strands), 'five_prime_operations_by_strand': dict(operations),
        'clipping_length_histogram': dict(lengths),
        'top_soft_clip_prefixes_biological_orientation': clips.most_common(30),
        'top_ambiguous_read_prefixes_biological_orientation': starts.most_common(30),
        'top_ambiguous_cigars': cigars.most_common(30),
        'notes': 'Reverse alignments use the last CIGAR operation; sequences are restored to original read orientation. Hard-clipped bases are absent from BAM. Counts diagnose strict 5-prime acceptance, not UMI validity.',
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bam', type=Path, required=True)
    p.add_argument('--max-reads', type=int, default=100000)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.max_reads < 1:
        p.error('--max-reads must be positive')
    if args.output.exists():
        p.error('Output exists; choose a new filename')
    result = inspect(args.bam, args.max_reads)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ('primary_mapped_examined', 'ambiguous_pct', 'exclusive_reasons')}, indent=2))


if __name__ == '__main__':
    main()
