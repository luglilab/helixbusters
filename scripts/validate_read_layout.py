#!/usr/bin/env python3
"""Check configured inline barcode orientation and offset before preparation."""
import argparse
import json
from pathlib import Path
import sys
from inspect_fastq_layout import inspect_fastq


def validate(result, barcode, orientation, offset, minimum_fraction):
    counts = {hit['orientation']: sum(p['reads'] for p in hit['positions'] if p['offset_0based'] == offset)
              for hit in result['barcode_search'] if hit['barcode'] == barcode}
    expected = counts.get(orientation, 0)
    alternative = counts.get('reverse_complement' if orientation == 'forward' else 'forward', 0)
    fraction = expected / result['reads_examined']
    passed = fraction >= minimum_fraction and expected > 0 and expected >= alternative
    return {'reads_examined': result['reads_examined'], 'configured_orientation': orientation,
            'barcode_offset': offset, 'configured_exact_pct': 100 * fraction,
            'alternative_exact_pct': 100 * alternative / result['reads_examined'],
            'minimum_fraction': minimum_fraction, 'passed': passed}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reads', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--barcode', required=True)
    p.add_argument('--orientation', choices=['forward', 'reverse_complement'], required=True)
    p.add_argument('--umi-length', type=int, default=8)
    p.add_argument('--max-reads', type=int, default=10000)
    p.add_argument('--minimum-fraction', type=float, default=.1)
    args = p.parse_args()
    if args.max_reads < 1 or args.umi_length < 1 or not 0 < args.minimum_fraction <= 1:
        p.error('Read/UMI lengths must be positive and fraction in (0,1]')
    result = inspect_fastq(args.reads, [args.barcode], args.max_reads,
                           max(80, args.umi_length + len(args.barcode)), 40, args.umi_length, args.umi_length)
    check = validate(result, args.barcode, args.orientation, args.umi_length, args.minimum_fraction)
    report = {'sample': args.sample, 'check': check, 'layout': result,
              'interpretation': 'First-read screen only; barcode matching does not prove biological identity or UMI validity. No orientation correction is performed.'}
    with Path(f'{args.sample}.layout.json').open('x') as handle:
        json.dump(report, handle, indent=2)
    with Path(f'{args.sample}.layout_mqc.json').open('x') as handle:
        json.dump({'id': 'helixbusters_layout', 'section_name': 'BLISS inline barcode layout',
                   'description': report['interpretation'], 'plot_type': 'table',
                   'data': {args.sample: check}}, handle, indent=2)
    if not check['passed']:
        raise SystemExit(f"Barcode layout check failed for {args.sample}: {check}. Review the layout report before changing metadata or thresholds.")
    print(f'{args.sample}: barcode layout passed ({check["configured_exact_pct"]:.2f}% exact)')


if __name__ == '__main__':
    main()
