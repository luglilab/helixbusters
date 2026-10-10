#!/usr/bin/env python3
"""Exploratory sample PCA of genomic-window counts, with depth sensitivity."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def pca(matrix):
    centered = matrix - matrix.mean(axis=0)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    variance = singular ** 2
    if variance.sum() == 0:
        raise ValueError('No variation remains for PCA')
    return u[:, :2] * singular[:2], variance[:2] / variance.sum(), vt[:2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--windows', nargs='+', type=int, default=[10000, 50000, 100000])
    parser.add_argument('--top-features', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=1729)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error('Choose a new output directory; existing results are protected')
    if args.top_features < 2:
        parser.error('--top-features must be >=2')
    metadata = pd.read_csv(args.analysis_dir / 'analysis.samples.tsv', sep='\t', keep_default_na=False)
    analysis_purpose = json.loads((args.analysis_dir / 'analysis.summary.json').read_text()).get('analysis_purpose', 'biological') if (args.analysis_dir / 'analysis.summary.json').exists() else 'biological'
    if metadata['sample'].duplicated().any() or len(metadata) < 1:
        raise ValueError('Sample metadata must be nonempty and unique')
    names = metadata['sample'].tolist()
    args.outdir.mkdir(parents=True)
    groups = metadata['group'].tolist()
    palette = dict(zip(sorted(set(groups)), plt.get_cmap('tab10').colors))
    fig, axes = plt.subplots(2, len(args.windows), figsize=(6 * len(args.windows), 10), squeeze=False)
    summary = {'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
               'analysis_purpose': analysis_purpose,
               'interpretation': 'Exploratory PCA, not a condition test. Donor pairing is not inferred. All-zero bins are omitted upstream. CPM measures relative genomic allocation, not absolute DSB burden.',
               'windows': {}, 'versions': {'numpy': np.__version__, 'pandas': pd.__version__, 'matplotlib': matplotlib.__version__}}
    coordinates = []
    multiqc = {}
    expected_totals = None
    for column, width in enumerate(args.windows):
        table = pd.read_csv(args.analysis_dir / f'windows_{width}.counts.tsv', sep='\t')
        counts = table[names].to_numpy(dtype=float).T
        if not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all():
            raise ValueError('Counts must be finite nonnegative integers')
        totals = counts.sum(axis=1).astype(np.int64)
        if expected_totals is not None and not np.array_equal(totals, expected_totals):
            raise ValueError('Window matrices have different library totals')
        expected_totals = totals
        reason = None
        if len(names) < 3:
            reason = 'Fewer than three samples'
        elif (totals <= 0).any():
            reason = 'At least one sample has no molecules'
        if reason:
            summary['windows'][str(width)] = {'status': 'skipped', 'reason': reason}
            multiqc[f'{width}_bp'] = {'status': reason}
            for ax in axes[:, column]:
                ax.text(.5, .5, f'{width} bp: PCA skipped\n{reason}', ha='center', transform=ax.transAxes)
                ax.set_axis_off()
            continue
        # Condition-independent abundance filter; no per-feature z scaling.
        eligible = (counts.sum(axis=0) >= 5) & ((counts > 0).sum(axis=0) >= 2)
        indices = np.flatnonzero(eligible)
        logcpm = np.log2(1 + counts / totals[:, None] * 1e6)
        variance = logcpm[:, indices].var(axis=0)
        indices = indices[np.argsort(-variance, kind='stable')[:args.top_features]]
        if len(indices) < 2:
            reason = 'Fewer than two eligible variable windows'
            summary['windows'][str(width)] = {'status': 'skipped', 'reason': reason}
            multiqc[f'{width}_bp'] = {'status': reason}
            for ax in axes[:, column]:
                ax.text(.5, .5, f'{width} bp: PCA skipped\n{reason}', ha='center', transform=ax.transAxes)
                ax.set_axis_off()
            continue
        if not np.any(variance > 0):
            reason = 'No variation in eligible windows'
            summary['windows'][str(width)] = {'status': 'skipped', 'reason': reason}
            multiqc[f'{width}_bp'] = {'status': reason}
            for ax in axes[:, column]:
                ax.text(.5, .5, f'{width} bp: PCA skipped\n{reason}', ha='center', transform=ax.transAxes)
                ax.set_axis_off()
            continue
        target = int(totals.min())
        rng = np.random.default_rng(args.seed)
        rare = np.vstack([rng.multivariate_hypergeometric(row.astype(np.int64), target) for row in counts])
        transformed_rare = np.log2(1 + rare / target * 1e6)
        record = {'library_molecules': dict(zip(names, map(int, totals))),
                  'eligible_windows': int(eligible.sum()), 'selected_windows': len(indices),
                  'abundance_filter': 'pooled count >=5 and nonzero in >=2 samples',
                  'transformation': 'log2(1 + CPM), feature-centered, no feature scaling',
                  'equal_depth_molecules': target, 'panels': {}}
        selected = table.iloc[indices][['region', 'chrom', 'start', 'end']].copy()
        for row, (label, data) in enumerate([('logCPM', logcpm), ('equal_depth_logCPM', transformed_rare)]):
            if not np.any(data[:, indices].var(axis=0) > 0):
                record['panels'][label] = {'status': 'skipped', 'reason': 'No variation after equal-depth sampling'}
                axes[row, column].set_axis_off()
                continue
            scores, explained, loadings = pca(data[:, indices])
            correlations = [float(np.corrcoef(scores[:, pc], np.log10(totals))[0, 1])
                            if np.std(np.log10(totals)) > 0 and np.std(scores[:, pc]) > 0 else None
                            for pc in range(2)]
            record['panels'][label] = {'explained_variance_ratio': explained.tolist(),
                                      'correlation_PC_with_original_log10_library_size': correlations}
            multiqc[f'{width}_bp_{label}'] = {'PC1_variance_pct': float(explained[0] * 100),
                'analysis_purpose': analysis_purpose,
                'minimum_library_molecules': int(totals.min()), 'maximum_library_molecules': int(totals.max()),
                'library_depth_fold_range': float(totals.max() / totals.min()),
                'depth_association_flag': any(abs(r) >= .7 for r in correlations if r is not None),
                'PC2_variance_pct': float(explained[1] * 100), 'PC1_depth_correlation': correlations[0],
                'PC2_depth_correlation': correlations[1], 'selected_windows': len(indices),
                'sampling_depth': target if row else 'full library'}
            ax = axes[row, column]
            for i, sample in enumerate(names):
                ax.scatter(*scores[i], color=palette[groups[i]], s=65)
                offset = (5, -14) if sample.startswith('HD1_') else (5, 10)
                ax.annotate(sample, scores[i], xytext=offset, textcoords='offset points', fontsize=9)
                coordinates.append({'window_bp': width, 'analysis': label, 'sample': sample,
                                    'group': groups[i], 'replicate': metadata.iloc[i]['replicate'],
                                    'donor': metadata.iloc[i]['donor'], 'library_molecules': int(totals[i]),
                                    'PC1': scores[i, 0], 'PC2': scores[i, 1]})
            ax.set_xlabel(f'PC1 ({explained[0]:.1%})')
            ax.set_ylabel(f'PC2 ({explained[1]:.1%})')
            ax.set_title(f'{width // 1000} kb | {label}\n{len(indices):,} variable windows')
            ax.margins(x=.22, y=.15)
            ax.spines[['top', 'right']].set_visible(False)
            for pc in range(2):
                selected[f'{label}_PC{pc+1}_loading'] = loadings[pc]
        selected.to_csv(args.outdir / f'windows_{width}.features.tsv', sep='\t', index=False)
        summary['windows'][str(width)] = record
    handles = [plt.Line2D([], [], marker='o', linestyle='', color=c, label=g) for g,c in palette.items()]
    fig.legend(handles=handles, loc='upper center', ncol=len(handles), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(args.outdir / 'windows_PCA.png', dpi=180)
    fig.savefig(args.outdir / 'windows_PCA.pdf')
    plt.close(fig)
    pd.DataFrame(coordinates, columns=['window_bp', 'analysis', 'sample', 'group', 'replicate',
                                     'donor', 'library_molecules', 'PC1', 'PC2']).to_csv(
        args.outdir / 'pca.coordinates.tsv', sep='\t', index=False)
    (args.outdir / 'pca.summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    (args.outdir / 'pca_mqc.json').write_text(json.dumps({
        'id': 'helixbusters_pca', 'section_name': 'Helixbusters genomic-window PCA',
        'description': 'Exploratory PCA; full-library and equal-depth panels: Analysis/PCA. Depth flag means absolute PC/depth correlation >=0.7, a descriptive diagnostic, not a test. Analysis purpose: ' + analysis_purpose,
        'plot_type': 'table', 'pconfig': {'id': 'helixbusters_pca_table', 'title': 'PCA variance and depth diagnostics'},
        'data': multiqc}, indent=2, allow_nan=False) + '\n')
    (args.outdir / 'README.md').write_text('''# Exploratory genomic-window PCA

One point per library. For input_titration, input levels come from one donor and
are not biological replicates. Columns: requested window sizes. Upper row:
log2(1+CPM) of full-library counts. Lower row: equal-depth molecule sampling
without replacement, then the same transformation. Both use the same windows
selected by full-library variance after a condition-independent abundance filter.
Features are centered but not scaled. Rarefaction is a sensitivity check, not a
replacement for count-based differential modeling, and one random seed is not
a stability analysis. PCA signs are arbitrary; explained variance uses all PCs
of the selected feature matrix. Library totals and PC-depth correlations are in
pca.summary.json. Files list exact features, loadings and sample coordinates.

Sparse low-count windows and log transformation may retain depth effects.
Grouping in a six-sample PCA is descriptive; it does not establish a treatment
effect. Donor metadata are preserved as supplied, without inferring pairing.
Absolute DSB load cannot be inferred from CPM or these PCA coordinates.
''')
    print(json.dumps(summary['windows'], indent=2))

if __name__ == '__main__':
    main()
