#!/usr/bin/env python3
"""Extract anchored BLISS UMIs and report exclusive preparation outcomes."""

import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.technical import prefix_distance


def prepare(reads, sample, barcode, orientation, umi_length, min_length, technical_prefix, outdir, technical_prefix_max_errors=0):
    if not barcode or set(barcode) - set('ACGT'):
        raise ValueError('Barcode must contain only A/C/G/T')
    if umi_length < 1 or min_length < 1:
        raise ValueError('UMI and minimum insert lengths must be positive')
    if orientation not in {'forward', 'reverse_complement'}:
        raise ValueError('Invalid barcode orientation')
    if technical_prefix and (len(technical_prefix) < 12 or set(technical_prefix) - set('ACGT')):
        raise ValueError('Technical prefix must contain at least 12 A/C/G/T bases, or be empty')
    if technical_prefix_max_errors not in {0, 1, 2} or (technical_prefix_max_errors and not technical_prefix):
        raise ValueError('Technical-prefix errors must be 0, 1 or 2; a motif is required for tolerant matching')
    observed_barcode = barcode if orientation == 'forward' else barcode.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fastq = outdir / f'{sample}.umi.fastq.gz'
    qc = outdir / f'{sample}.preparation.json'
    custom = outdir / f'{sample}.preparation_mqc.json'
    if any(p.exists() for p in (fastq, qc, custom)):
        raise ValueError('Preparation outputs already exist; use a new directory')
    counts = Counter({key: 0 for key in ('input_reads', 'barcode_mismatch_reads', 'invalid_umi_reads',
                                       'technical_prefix_reads', 'short_insert_reads', 'output_reads')})
    opener = gzip.open if str(reads).lower().endswith('.gz') else open
    end = umi_length + len(barcode)
    with opener(reads, 'rt', encoding='ascii') as source, fastq.open('xb') as raw:
        # Fixed gzip timestamp makes the compressed output reproducible.
        with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as destination:
            while True:
                header = source.readline()
                if not header:
                    break
                header = header.rstrip('\r\n')
                sequence = source.readline().rstrip('\r\n').upper()
                plus = source.readline().rstrip('\r\n')
                quality = source.readline().rstrip('\r\n')
                counts['input_reads'] += 1
                if (not header.startswith('@') or not plus.startswith('+') or not sequence or
                        len(sequence) != len(quality) or set(sequence) - set('ACGTN')):
                    raise ValueError(f'Malformed FASTQ record {counts["input_reads"]}')
                if sequence[umi_length:end] != observed_barcode:
                    counts['barcode_mismatch_reads'] += 1
                    continue
                umi = sequence[:umi_length]
                if len(umi) != umi_length or set(umi) - set('ACGT'):
                    counts['invalid_umi_reads'] += 1
                    continue
                insert = sequence[end:]
                if technical_prefix and prefix_distance(insert, technical_prefix, technical_prefix_max_errors) <= technical_prefix_max_errors:
                    counts['technical_prefix_reads'] += 1
                    continue
                if len(insert) < min_length:
                    counts['short_insert_reads'] += 1
                    continue
                identifier, *description = header.split(maxsplit=1)
                new_header = identifier + '_' + umi + (' ' + description[0] if description else '')
                destination.write(f'{new_header}\n{insert}\n+\n{quality[end:]}\n'.encode('ascii'))
                counts['output_reads'] += 1
    parameters = {'barcode': barcode, 'observed_barcode': observed_barcode,
                  'barcode_orientation': orientation, 'umi_length': umi_length,
                  'minimum_insert_length': min_length, 'excluded_exact_technical_prefix': technical_prefix if technical_prefix_max_errors == 0 else None,
                  'technical_prefix': technical_prefix, 'technical_prefix_max_errors': technical_prefix_max_errors,
                  'barcode_match': 'exact, anchored immediately after UMI; no mismatch rescue',
                  'technical_filter': 'exclude entire read by anchored motif edit distance; substitutions/insertions/deletions; no fixed trimming'}
    data = dict(counts)
    data['preparation_retained_pct'] = 100 * counts['output_reads'] / counts['input_reads'] if counts['input_reads'] else None
    qc.write_text(json.dumps({'sample': sample, 'counts': dict(counts), 'parameters': parameters}, indent=2) + '\n')
    custom.write_text(json.dumps({
        'id': 'helixbusters_preparation', 'section_name': 'Helixbusters read preparation',
        'description': f'Counts start from original FASTQ records. Exclusions are exclusive: barcode, UMI, anchored technical prefix (up to {technical_prefix_max_errors} edits), insert length. The technical-prefix category does not establish adapter-dimer identity.',
        'plot_type': 'table', 'pconfig': {'id': 'helixbusters_preparation_table', 'title': 'Original reads and BLISS preparation'},
        'data': {sample: {k: v for k, v in data.items() if v is not None}},
    }, indent=2) + '\n')
    print(json.dumps({'sample': sample, **data, 'parameters': parameters}, indent=2))
    if not counts['output_reads']:
        raise ValueError('No reads retained; inspect preparation JSON and verify the library layout')
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--barcode', required=True)
    parser.add_argument('--barcode-orientation', choices=('forward', 'reverse_complement'), required=True)
    parser.add_argument('--umi-length', type=int, default=8)
    parser.add_argument('--minimum-insert-length', type=int, default=20)
    parser.add_argument('--technical-prefix', default='')
    parser.add_argument('--technical-prefix-max-errors', type=int, choices=(0, 1, 2), default=0)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()
    if not args.sample or any(c in args.sample for c in '/\\\t\r\n') or args.sample in {'.', '..'}:
        parser.error('Invalid sample identifier')
    prepare(args.reads, args.sample, args.barcode.upper(), args.barcode_orientation, args.umi_length,
            args.minimum_insert_length, args.technical_prefix.upper(), args.outdir, args.technical_prefix_max_errors)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, EOFError, UnicodeError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
