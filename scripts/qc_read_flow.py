#!/usr/bin/env python3
"""Report sequential read losses with explicit raw-read denominators."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.reporting import write_json, write_table


def read_flow(sample, preparation):
    counts = preparation['counts']; metrics = sample['metrics']; total = counts['input_reads']
    if total <= 0 or preparation['sample'] != sample['sample']:
        raise ValueError('Invalid input depth or inconsistent sample identity')
    if counts['output_reads'] != metrics['primary_records']:
        raise ValueError('Prepared reads and primary alignment counts differ')
    losses = {key: counts.get(key, 0) for key in ('barcode_mismatch_reads', 'invalid_umi_reads', 'technical_prefix_reads', 'short_insert_reads')}
    losses.update({key: metrics.get(key, 0) for key in ('excluded_unmapped', 'excluded_low_mapq', 'excluded_blacklist', 'excluded_mitochondrial', 'excluded_noncanonical', 'excluded_qc_fail', 'excluded_unknown_mapq', 'ambiguous_five_prime_reads', 'duplicate_reads')})
    record = {'group': sample['group'], 'raw_reads': total, 'prepared_reads': counts['output_reads'],
              'usable_molecules': metrics['deduplicated_molecules'],
              'prepared_pct_of_raw': 100 * counts['output_reads'] / total,
              'usable_molecules_pct_of_raw': 100 * metrics['deduplicated_molecules'] / total}
    record.update({f'{key}_pct_of_raw': 100 * value / total for key, value in losses.items()})
    dominant = max(losses, key=losses.get)
    record['largest_recorded_loss'] = dominant
    record['largest_recorded_loss_pct_of_raw'] = 100 * losses[dominant] / total
    retained = metrics['retained_records']
    record['ambiguous_5prime_pct_of_retained'] = 100 * metrics['ambiguous_five_prime_reads'] / retained if retained else None
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samples', nargs='+', required=True)
    p.add_argument('--preparations', nargs='+', required=True)
    p.add_argument('--purpose', choices=['biological', 'input_titration'], default='biological')
    args = p.parse_args()
    preparations = {}
    for path in args.preparations:
        if path.endswith('.preparation.json'):
            value = json.loads(Path(path).read_text())
            if value['sample'] in preparations:
                raise ValueError('Duplicate preparation summary')
            preparations[value['sample']] = value
    rows = {}
    for path in args.samples:
        sample = json.loads(Path(path).read_text())
        if sample['sample'] in rows:
            raise ValueError('Duplicate sample summary')
        rows[sample['sample']] = read_flow(sample, preparations[sample['sample']])
    write_table('helixbusters_read_flow.tsv', rows)
    write_json('helixbusters_read_flow_mqc.json', {
        'id': 'helixbusters_read_flow', 'section_name': 'BLISS read retention and losses',
        'description': 'Loss fractions use original raw reads; ambiguous 5prime also reports the retained-alignment denominator. Supplementary alignment records are not counted as independent input reads. Largest loss is descriptive. Analysis purpose: ' + args.purpose,
        'plot_type': 'table', 'data': rows})


if __name__ == '__main__':
    main()
