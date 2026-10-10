#!/usr/bin/env python3
"""Paired categorical time-course analysis of relative DSB counts."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.design import validate_design
from helixbusters.differential import integer_counts, sha256, paired_effects
from helixbusters.reporting import validate_label


def validate_timecourse(metadata, times):
    if len(times) < 3 or len(set(times)) != len(times):
        raise ValueError('Specify at least three distinct time points, baseline first')
    for label in times:
        validate_label(label)
    if metadata['sample'].duplicated().any() or set(metadata.group) != set(times):
        raise ValueError('Time points must match all groups exactly; sample names must be unique')
    validate_design(metadata.to_dict('records'), 'paired')
    if metadata.donor.nunique() < 3:
        raise ValueError('Time-course inference requires at least three biological donors')
    return metadata.sort_values('sample', kind='stable').reset_index(drop=True)


def run(analysis, output, times, minimum_count=5, minimum_samples=2, fdr=0.05):
    if output.exists():
        raise ValueError('Choose a new output directory')
    if minimum_count < 1 or minimum_samples < 2 or (not 0 < fdr < 1):
        raise ValueError('Invalid abundance or FDR thresholds')
    metadata = validate_timecourse(pd.read_csv(analysis / 'analysis.samples.tsv', sep='\t', dtype=str, keep_default_na=False), times)
    names = metadata['sample'].tolist()
    if minimum_samples > len(names):
        raise ValueError('Abundance filter exceeds sample count')
    source = json.loads((analysis / 'analysis.summary.json').read_text())
    totals = np.array([source['samples'][n]['molecules'] for n in names], dtype=np.int64)
    if (totals <= 0).any():
        raise ValueError('Zero-depth libraries cannot enter inference')
    metadata['total_molecules'] = totals
    paths = sorted(analysis.glob('windows_*.counts.tsv'))
    paths += [p for p in sorted((analysis / 'Annotation/GeneSignal').glob('genes.*.counts.tsv')) if p.name in ('genes.promoter.counts.tsv', 'genes.gene_body.counts.tsv')]
    if not paths:
        raise ValueError('No count matrices found')
    script = Path(__file__).with_name('fit_dsb_deseq2.R')
    output.mkdir(parents=True)
    metadata.to_csv(output / 'samples.tsv', sep='\t', index=False)
    report = {
        'model': '~ donor + condition (categorical time)',
        'global_test': 'LRT versus ~ donor',
        'timepoints': times,
        'baseline': times[0],
        'minimum_count': minimum_count,
        'minimum_samples': minimum_samples,
        'fdr': fdr,
        'multiplicity': (
            'Global LRT: BH within family. Wald baseline contrasts: joint BH '
            'across features and contrasts within family, without selection '
            'by LRT. No joint FDR across widths or gene contexts.'
        ),
        'interpretation': (
            'Relative signal, not absolute DSB per cell. Without time-matched '
            'unstimulated controls, stimulation and culture time are not separable.'
        ),
        'input_hashes': {
            'metadata': sha256(analysis / 'analysis.samples.tsv'),
            'summary': sha256(analysis / 'analysis.summary.json'),
            'model': sha256(script),
            'reporter': sha256(__file__),
        },
        'families': {},
    }
    mqc = {}
    for path in paths:
        family = path.name.removesuffix('.counts.tsv')
        folder = output / family
        folder.mkdir()
        table = pd.read_csv(path, sep='\t', keep_default_na=False)
        raw = integer_counts(table, names)
        keys = ['chrom', 'gene_id'] if 'gene_id' in table else ['region']
        if table.duplicated(keys).any() or (raw.sum(axis=0) > totals).any():
            raise ValueError(f'Duplicate features or invalid depth: {path}')
        if path.name.startswith('windows_') and (not np.array_equal(raw.sum(axis=0), totals)):
            raise ValueError('Window counts do not conserve library depths')
        table.insert(0, 'feature_id', [f'f{i}' for i in range(len(table))])
        eligible = (raw >= minimum_count).sum(axis=1) >= minimum_samples
        table['passed_abundance_filter'] = eligible
        table.to_csv(folder / 'feature_universe.tsv', sep='\t', index=False)
        record = {'input_sha256': sha256(path), 'features': len(table), 'eligible': int(eligible.sum())}
        if eligible.sum() < 20:
            record.update(status='skipped', reason='Fewer than 20 abundance-eligible features')
        else:
            table.loc[eligible, ['feature_id', *names]].to_csv(folder / 'model_counts.tsv', sep='\t', index=False)
            with (folder / 'model.log').open('x') as log:
                subprocess.run(['Rscript', '--vanilla', str(script), str(folder / 'model_counts.tsv'), str(output / 'samples.tsv'), str(folder), 'timecourse', ','.join(times[1:]), times[0]], stdout=log, stderr=subprocess.STDOUT, check=True)
            for normalization in ('poscounts', 'library_total'):
                for contrast in ['global', *[f'{t}_vs_{times[0]}' for t in times[1:]]]:
                    filename = f'{normalization}.{contrast}.tsv'
                    fitted = pd.read_csv(folder / filename, sep='\t')
                    result = table.merge(fitted, on='feature_id', how='left', validate='one_to_one')
                    result['significant'] = result.padj.le(fdr)
                    if contrast != 'global':
                        time = next((t for t in times[1:] if contrast == f'{t}_vs_{times[0]}'))
                        donors, ratios, median, agree, _, stable = paired_effects(raw / totals[None, :] * 1000000.0, metadata, time, times[0])
                        result['paired_median_log2_CPM_ratio'] = median
                        result['all_donors_same_direction'] = agree
                        result['leave_one_donor_out_direction_stable'] = stable
                        for i, donor in enumerate(donors):
                            result[f'{donor}_log2_CPM_ratio'] = ratios[:, i]
                    elif normalization == 'poscounts':
                        record['global_significant'] = int(result.significant.sum())
                    result.to_csv(folder / filename, sep='\t', index=False)
            record['status'] = 'fitted'
        report['families'][family] = record
        mqc[family] = {'eligible': record['eligible'], 'status': record['status'], 'global_significant': record.get('global_significant', 0)}
    (output / 'timecourse.summary.json').write_text(json.dumps(report, indent=2))
    (output / 'differential_mqc.json').write_text(json.dumps({'id': 'helixbusters_timecourse', 'section_name': 'Paired DSB time course', 'description': report['multiplicity'], 'plot_type': 'table', 'data': mqc}, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--timepoints', nargs='+', required=True, help='Categorical levels, baseline first')
    parser.add_argument('--minimum-count', type=int, default=5)
    parser.add_argument('--minimum-samples', type=int, default=2)
    parser.add_argument('--fdr', type=float, default=0.05)
    args = parser.parse_args()
    try:
        run(args.analysis_dir, args.outdir, args.timepoints, args.minimum_count, args.minimum_samples, args.fdr)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
