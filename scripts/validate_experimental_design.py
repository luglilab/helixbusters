#!/usr/bin/env python3
"""Write validated experimental-design metadata and a MultiQC table."""
import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.design import validate_design


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--design', choices=('paired', 'unpaired', 'unspecified'), default='unspecified')
    p.add_argument('--purpose', choices=('biological', 'input_titration'), default='biological')
    args = p.parse_args()
    with args.manifest.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'sample', 'group', 'replicate'}.issubset(reader.fieldnames or []):
            p.error('Manifest requires sample, group and replicate columns')
        report = validate_design(list(reader), args.design)
    report['analysis_purpose'] = args.purpose
    if args.purpose == 'input_titration':
        members = report['sample_metadata']
        if args.design != 'unspecified' or report['declared_donors'] != 1 or any(not r['donor'] for r in members) or len(members) != len(report['conditions']):
            p.error('Input titration requires one declared donor and one library per input level, with unspecified design')
        report['notes'] = 'Technical input titration from one donor. Input levels are not independent biological replicates; no biological differential inference.'
    with Path('design.summary.json').open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    with Path('design.metadata.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=('sample', 'group', 'replicate', 'donor'), delimiter='\t')
        writer.writeheader()
        writer.writerows(report['sample_metadata'])
    with Path('design.metadata_mqc.json').open('x') as handle:
        json.dump({'id': 'helixbusters_design', 'section_name': 'Helixbusters experimental design',
                   'description': f"Declared design: {args.design}. {report['status']}. {report['notes']}",
                   'plot_type': 'table', 'pconfig': {'id': 'helixbusters_design_table', 'title': 'Experimental units'},
                   'data': {row['sample']: {'condition': row['group'], 'replicate': row['replicate'],
                                            'donor': row.get('donor') or 'not declared', 'design': args.design,
                                            'analysis_purpose': args.purpose}
                            for row in report['sample_metadata']}}, handle, indent=2)
    print(f"Validated {args.design} design: {report['samples']} samples, {len(report['conditions'])} conditions")


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
