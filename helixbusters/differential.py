"""Count-based condition comparisons and descriptive robustness diagnostics."""
from pathlib import Path
import hashlib
import json
import subprocess

import numpy as np
import pandas as pd

from helixbusters.design import validate_design
from helixbusters.reporting import validate_label


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def contrast_metadata(metadata, design, numerator, denominator):
    """Never infer donor identities or choose a reference condition implicitly."""
    if design not in ('paired', 'unpaired'):
        raise ValueError('Differential analysis requires paired or unpaired design')
    for label in (numerator, denominator):
        validate_label(label)
    if numerator == denominator:
        raise ValueError('Contrast conditions must differ')
    if metadata['sample'].duplicated().any():
        raise ValueError('Duplicate sample metadata')
    metadata = metadata.copy()
    metadata['donor'] = metadata['donor'].astype(str).str.strip()
    validate_design(metadata.to_dict('records'), design)
    subset = metadata[metadata.group.isin([numerator, denominator])].copy()
    sizes = subset.groupby('group').size()
    if any(sizes.get(group, 0) < 3 for group in (numerator, denominator)):
        raise ValueError('Differential analysis requires >=3 biological samples per contrast condition')
    # Validate the selected contrast too, including donor matching.
    validate_design(subset.to_dict('records'), design)
    return subset.sort_values('sample', kind='stable').reset_index(drop=True)


def integer_counts(table, samples):
    raw = table[samples].to_numpy(dtype=float)
    if not np.isfinite(raw).all() or (raw < 0).any() or not np.equal(raw, np.floor(raw)).all():
        raise ValueError('Counts must be finite nonnegative integers')
    if (raw > np.iinfo(np.int32).max).any():
        raise ValueError('A feature count exceeds the DESeq2 integer limit')
    return raw.astype(np.int64)


def paired_effects(values, metadata, numerator, denominator):
    """Descriptive log ratios with a declared pseudocount, never a test."""
    pairs = []
    donors = sorted(metadata.donor.unique())
    for donor in donors:
        num = metadata.index[(metadata.donor == donor) & (metadata.group == numerator)][0]
        den = metadata.index[(metadata.donor == donor) & (metadata.group == denominator)][0]
        pairs.append(np.log2((values[:, num] + 1) / (values[:, den] + 1)))
    ratios = np.column_stack(pairs)
    median = np.median(ratios, axis=1)
    same_direction = ((ratios > 0).all(axis=1) | (ratios < 0).all(axis=1))
    omitted = np.column_stack([np.median(np.delete(ratios, i, axis=1), axis=1)
                               for i in range(len(donors))])
    stable = (np.sign(omitted) == np.sign(median[:, None])).all(axis=1) & (median != 0)
    return donors, ratios, median, same_direction, omitted, stable


def depth_sensitivity(raw, totals, groups, iterations, seed, minimum_molecules=2, minimum_replicates=2):
    """Sample disjoint molecular categories and an unassigned residual category."""
    if (raw.sum(axis=0) > totals).any() or (totals <= 0).any():
        raise ValueError('Feature counts exceed total retained molecules, or a library is empty')
    target = int(totals.min())
    rng = np.random.default_rng(seed)
    hits = {group: np.zeros(len(raw), dtype=np.int64) for group in groups}
    numbers = {group: [] for group in groups}
    categories = [np.append(raw[:, i], totals[i] - raw[:, i].sum()) for i in range(raw.shape[1])]
    for _ in range(iterations):
        sampled = np.column_stack([rng.multivariate_hypergeometric(row, target)[:-1] for row in categories])
        for group, columns in groups.items():
            support = (sampled[:, columns] >= minimum_molecules).sum(axis=1) >= minimum_replicates
            hits[group] += support
            numbers[group].append(int(support.sum()))
    return {group: value / iterations for group, value in hits.items()}, numbers, target


def run_comparison(analysis_dir, outdir, metadata, design, numerator, denominator,
                   minimum_count=5, minimum_samples=2, fdr=0.05, iterations=50, seed=1729):
    """One explicit contrast, independent feature-family FDR, read-only inputs."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    analysis_dir, outdir = Path(analysis_dir), Path(outdir)
    if outdir.exists():
        raise ValueError('Choose a new output directory; existing results are protected')
    if minimum_count < 1 or minimum_samples < 2 or iterations < 2 or not 0 < fdr < 1:
        raise ValueError('Invalid abundance, FDR or resampling settings')
    source_metadata = pd.read_csv(analysis_dir / 'analysis.samples.tsv', sep='\t', keep_default_na=False, dtype=str)
    if set(metadata['sample']) != set(source_metadata['sample']):
        raise ValueError('Metadata must identify exactly the source analysis samples')
    for field in ('group', 'replicate'):
        original = source_metadata.set_index('sample')[field].sort_index()
        supplied = metadata.set_index('sample')[field].sort_index()
        if not original.equals(supplied):
            raise ValueError(f'Metadata override changes source {field}')
    declared = source_metadata.set_index('sample').donor.str.strip()
    supplied_donors = metadata.set_index('sample').donor.str.strip()
    if not supplied_donors.loc[declared[declared != ''].index].equals(declared[declared != '']):
        raise ValueError('Metadata override changes an already declared donor')
    metadata = contrast_metadata(metadata, design, numerator, denominator)
    names = metadata['sample'].tolist()
    if minimum_samples > len(names):
        raise ValueError('Abundance filter requires too many samples')
    analysis_summary = json.loads((analysis_dir / 'analysis.summary.json').read_text())
    totals = np.array([analysis_summary['samples'][name]['molecules'] for name in names], dtype=np.int64)
    if (totals <= 0).any():
        raise ValueError('Differential analysis cannot use zero-depth samples')
    metadata['total_molecules'] = totals
    paths = sorted(analysis_dir.glob('windows_*.counts.tsv'))
    paths += [p for p in sorted((analysis_dir / 'Annotation/GeneSignal').glob('genes.*.counts.tsv'))
              if p.name in ('genes.promoter.counts.tsv', 'genes.gene_body.counts.tsv')]
    if not paths:
        raise ValueError('No window or promoter/gene-body count matrices were found')
    rscript = Path(__file__).resolve().parents[1] / 'scripts/fit_dsb_deseq2.R'
    subprocess.run(['Rscript', '--vanilla', '-e',
                    'if (!requireNamespace("DESeq2", quietly=TRUE)) stop("Install DESeq2 in the analysis environment")'],
                   check=True)
    outdir.mkdir(parents=True)
    metadata.to_csv(outdir / 'samples.tsv', sep='\t', index=False)
    groups = {group: np.flatnonzero(metadata.group.to_numpy() == group)
              for group in (numerator, denominator)}
    summary = {'design': design, 'model_formula': '~ donor + condition' if design == 'paired' else '~ condition',
               'contrast': {'numerator': numerator, 'denominator': denominator},
               'parameters': {'minimum_count': minimum_count, 'minimum_samples': minimum_samples,
                              'fdr': fdr, 'iterations': iterations, 'seed': seed},
               'normalization': 'Primary: DESeq2 poscounts per feature family. Sensitivity: total retained molecules. Both measure relative signal; no absolute DSB calibration.',
               'multiplicity': 'BH within each feature family; no joint FDR across widths or gene contexts. Peaks are excluded because discovery uses the same samples.',
               'robustness': 'Equal-depth support and leave-one-donor-out CPM medians are descriptive, not additional significance tests. Leave-one-donor-out does not refit the model.',
               'versions': {'numpy': np.__version__, 'pandas': pd.__version__, 'matplotlib': matplotlib.__version__},
               'inputs': {'analysis_summary_sha256': sha256(analysis_dir / 'analysis.summary.json'),
                          'source_metadata_sha256': sha256(analysis_dir / 'analysis.samples.tsv'),
                          'model_script_sha256': sha256(rscript),
                          'reporting_module_sha256': sha256(__file__)},
               'families': {}}
    multiqc = {}
    for path in paths:
        family = path.name.removesuffix('.counts.tsv')
        folder = outdir / family
        folder.mkdir()
        table = pd.read_csv(path, sep='\t', keep_default_na=False)
        raw = integer_counts(table, names)
        if path.name.startswith('windows_') and not np.array_equal(raw.sum(axis=0), totals):
            raise ValueError('Window counts do not conserve total retained molecules')
        if (raw.sum(axis=0) > totals).any():
            raise ValueError('Gene counts exceed total retained molecules')
        identifiers = ['chrom', 'gene_id'] if 'gene_id' in table else ['region']
        if table.duplicated(identifiers).any():
            raise ValueError(f'Duplicate feature IDs in {path}')
        # Stable numeric keys avoid interpreting versioned gene IDs or chromosome aliases.
        table.insert(0, 'feature_id', [f'f{i}' for i in range(len(table))])
        eligible = (raw >= minimum_count).sum(axis=1) >= minimum_samples
        features = table.drop(columns=names).copy()
        features['passed_abundance_filter'] = eligible
        features.to_csv(folder / 'feature_universe.tsv', sep='\t', index=False)
        frequencies, counts_per_draw, target = depth_sensitivity(raw, totals, groups, iterations, seed)
        result = table.copy()
        for group, columns in groups.items():
            result[f'{group}_mean_CPM'] = (raw[:, columns] / totals[columns] * 1e6).mean(axis=1)
            result[f'{group}_equal_depth_support_frequency'] = frequencies[group]
        cpm = raw / totals[None, :] * 1e6
        if design == 'paired':
            donors, ratios, median, concordant, omitted, stable = paired_effects(cpm, metadata, numerator, denominator)
            result['paired_median_log2_CPM_ratio'] = median
            result['all_donors_same_direction'] = concordant
            result['leave_one_donor_out_direction_stable'] = stable
            for i, donor in enumerate(donors):
                result[f'{donor}_log2_CPM_ratio'] = ratios[:, i]
                result[f'omit_{donor}_median_log2_CPM_ratio'] = omitted[:, i]
        pd.DataFrame({'iteration': np.arange(iterations),
                      **{f'{g}_supported_features': counts_per_draw[g] for g in groups}}).to_csv(
            folder / 'equal_depth_draws.tsv', sep='\t', index=False)
        record = {'input_sha256': sha256(path), 'features': len(raw), 'eligible_features': int(eligible.sum()),
                  'equal_depth_molecules': target, 'support_rule': '>=2 molecules in >=2 biological samples within condition'}
        if eligible.sum() < 20:
            record.update(status='skipped', reason='Fewer than 20 abundance-eligible features for dispersion trend estimation')
            result.to_csv(folder / 'robustness.tsv', sep='\t', index=False)
            multiqc[family] = {'status': record['reason'], 'eligible': int(eligible.sum())}
            summary['families'][family] = record
            continue
        result.loc[eligible, ['feature_id', *names]].to_csv(folder / 'model_counts.tsv', sep='\t', index=False)
        with (folder / 'model.log').open('x') as log:
            subprocess.run(['Rscript', '--vanilla', str(rscript), str(folder / 'model_counts.tsv'),
                            str(outdir / 'samples.tsv'), str(folder), design, numerator, denominator],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        for mode, prefix in [('poscounts', ''), ('library_total', 'library_total_')]:
            fitted = pd.read_csv(folder / f'{mode}.results.tsv', sep='\t').set_index('feature_id')
            for column in fitted:
                result[f'{prefix}{column}'] = result.feature_id.map(fitted[column])
            factors = pd.read_csv(folder / f'{mode}.size_factors.tsv', sep='\t').set_index('sample').loc[names]
            record[f'{mode}_size_factors'] = factors.size_factor.to_dict()
        result['passed_abundance_filter'] = eligible
        result['test_status'] = np.where(~eligible, 'below_abundance_filter',
                                         np.where(result.pvalue.notna(), 'tested', 'unavailable_or_Cooks_filtered'))
        failed = eligible & result.beta_converged.eq(False)
        result.loc[failed, 'test_status'] = 'coefficient_fit_failed'
        result['significant'] = result.padj.le(fdr)
        result['normalization_direction_agrees'] = (np.sign(result.log2FoldChange) ==
                                                   np.sign(result.library_total_log2FoldChange)) & result.log2FoldChange.notna()
        result['significant_both_normalizations'] = result.significant & result.library_total_padj.le(fdr)
        if 'effective_assignable_bp' in result:
            for group in groups:
                result[f'{group}_mean_CPM_per_effective_kb'] = np.divide(
                    result[f'{group}_mean_CPM'].to_numpy() * 1000, result.effective_assignable_bp.to_numpy(),
                    out=np.full(len(result), np.nan), where=result.effective_assignable_bp.to_numpy() > 0)
        result = result.sort_values(['padj', 'pvalue', 'feature_id'], na_position='last', kind='stable')
        result.to_csv(folder / 'results.tsv', sep='\t', index=False)
        result[result.significant].to_csv(folder / 'significant.tsv', sep='\t', index=False)
        for group, positive in [(numerator, True), (denominator, False)]:
            direction = result.log2FoldChange.gt(0) if positive else result.log2FoldChange.lt(0)
            result[result.significant & direction].to_csv(folder / f'{group}.higher_relative_signal.tsv', sep='\t', index=False)
        record.update(status='fitted', tested_features=int(result.pvalue.notna().sum()),
                      significant=int(result.significant.sum()),
                      significant_both_normalizations=int(result.significant_both_normalizations.sum()),
                      model_warning_file='model_warnings.txt', log2FoldChange='unshrunk DESeq2 maximum-likelihood estimate')
        if design == 'paired':
            record['all_donors_same_direction'] = int(result.all_donors_same_direction.sum())
            record['leave_one_donor_out_direction_stable'] = int(result.leave_one_donor_out_direction_stable.sum())
        fig, panels = plt.subplots(2, 2, figsize=(11, 8))
        axes = panels.ravel()
        valid = result.log2FoldChange.notna()
        axes[0].scatter(np.log10(result.loc[valid, 'baseMean'] + 1), result.loc[valid, 'log2FoldChange'],
                        c=np.where(result.loc[valid, 'significant'], '#bc3838', '#888888'), s=6, alpha=.5)
        axes[0].axhline(0, color='black', lw=.7)
        axes[0].set(xlabel='log10(1 + mean normalized count)', ylabel=f'log2({numerator}/{denominator})', title='Model effect (unshrunk)')
        axes[1].scatter(result.log2FoldChange, result.library_total_log2FoldChange, s=6, alpha=.4)
        axes[1].axhline(0, color='black', lw=.7); axes[1].axvline(0, color='black', lw=.7)
        axes[1].set(xlabel='DESeq2 poscounts log2 fold change', ylabel='Library-total log2 fold change', title='Normalization sensitivity')
        for i, (group, columns) in enumerate(groups.items()):
            axes[2].scatter(np.full(iterations, i + .12), counts_per_draw[group], alpha=.35, s=10)
            original = int(((raw[:, columns] >= 2).sum(axis=1) >= 2).sum())
            axes[2].scatter(i - .12, original, marker='s', color='black', s=35)
        axes[2].set_xticks(range(len(groups)), list(groups))
        axes[2].set(ylabel='Replicate-supported features', title=f'Full library (square) / equal depth ({target:,})')
        if design == 'paired':
            for donor in donors:
                axes[3].scatter(result.paired_median_log2_CPM_ratio,
                                result[f'omit_{donor}_median_log2_CPM_ratio'], s=4, alpha=.15, label=f'omit {donor}')
            axes[3].axhline(0, color='black', lw=.7); axes[3].axvline(0, color='black', lw=.7)
            axes[3].set(xlabel='Median log2 CPM ratio, all donors', ylabel='Median ratio after omitting one donor', title='Descriptive donor influence')
            axes[3].legend(fontsize=8)
        else:
            axes[3].set_axis_off()
        fig.suptitle(f'{family}: {numerator} / {denominator}')
        fig.tight_layout(); fig.savefig(folder / 'diagnostics.pdf'); fig.savefig(folder / 'diagnostics.png', dpi=150)
        plt.close(fig)
        summary['families'][family] = record
        multiqc[family] = {key: record[key] for key in ('status', 'eligible_features', 'tested_features', 'significant', 'significant_both_normalizations', 'equal_depth_molecules')}
    (outdir / 'differential.summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    (outdir / 'differential_mqc.json').write_text(json.dumps({
        'id': 'helixbusters_differential', 'section_name': 'Helixbusters relative DSB signal comparisons',
        'description': f'{numerator} / {denominator}; {summary["model_formula"]}. BH within feature family. Relative signal, not absolute DSB burden. Robustness diagnostics: Analysis/Differential.',
        'plot_type': 'table', 'pconfig': {'id': 'helixbusters_differential_table', 'title': 'Differential signal and normalization sensitivity'},
        'data': multiqc}, indent=2, allow_nan=False) + '\n')
    (outdir / 'README.md').write_text('''# Relative DSB signal comparison

Positive log2FoldChange means higher relative signal in the numerator condition.
results.tsv retains the complete feature universe, original integer counts,
abundance eligibility, primary DESeq2 Wald results, dispersion/Cook's diagnostics,
and library-total normalization sensitivity. Missing tests are not negative results.
significant.tsv applies BH FDR within that feature family. Condition higher-signal
lists contain significant features with the corresponding effect direction.
FDR is not pooled across window widths or promoter/gene-body summaries. Choose a
primary family before interpreting other families as sensitivity analyses.
MACS consensus and combined gene counts are excluded from these tests.

Primary normalization estimates DESeq2 poscounts size factors within each tested
family; it assumes adequate stable relative signal. Library-total sensitivity
uses all retained deduplicated molecules, including intergenic/unassigned signal.
Neither estimates absolute DSB burden per cell. Inspect normalization disagreement,
dispersion plots and model warnings before interpreting candidates. Effect sizes
are unshrunk; low-count estimates can be unstable. Failed convergence is not tested.

Equal-depth support frequency uses repeated sampling without replacement down to
the smallest library, including a residual category for molecules outside each
family. Support means >=2 molecules in >=2 samples within condition. Frequencies
are depth sensitivity diagnostics, not probabilities of biological reproducibility.
Paired CPM ratios use log2((numerator CPM + 1)/(denominator CPM + 1)). Leave-one-donor-
out medians do not refit the model and have no p-values. No such pairing diagnostics
are generated for independent samples. Robustness columns do not redefine FDR.
''')
    return summary
