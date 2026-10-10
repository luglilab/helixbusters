#!/usr/bin/env python3
"""Annotate exact molecular DSB ends and exploratory regions with a build-matched GTF."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.annotation import FEATURES, GTFAnnotation
from helixbusters.genomes import canonical_lengths, normalize_genome, GenomeFilter
from helixbusters.regions import load_sites
from helixbusters.reporting import validate_label
from helixbusters.gene_signal import write_gene_signal


def write_regions(annotation, regions, destination):
    with destination.open('x') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['region', 'chrom', 'start', 'end', 'dominant_feature', 'mixed_features',
                         *[f'{feature}_bp' for feature in FEATURES], 'gene_ids', 'gene_names', 'gene_biotypes'])
        for region, chrom, start, end in regions:
            dominant, coverage, genes = annotation.annotate(chrom, start, end)
            writer.writerow([region, chrom, start, end, dominant, sum(value > 0 for value in coverage.values()) > 1,
                             *[coverage[feature] for feature in FEATURES], *annotation.gene_fields(chrom, genes)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', nargs='+', required=True)
    parser.add_argument('--counts', nargs='+', required=True)
    parser.add_argument('--headers', nargs='+', required=True)
    parser.add_argument('--design-file', required=True)
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--gtf-genome', required=True)
    parser.add_argument('--genome', required=True)
    parser.add_argument('--promoter-upstream', type=int, default=2000)
    parser.add_argument('--promoter-downstream', type=int, default=500)
    parser.add_argument('--gene-min-reps', type=int, default=2)
    parser.add_argument('--gene-min-molecules', type=int, default=2)
    masking = parser.add_mutually_exclusive_group()
    masking.add_argument('--environment-file', help='Pipeline environment.json containing the mapping blacklist and checksum')
    masking.add_argument('--blacklist', help='Build-matched BED for a standalone annotation run')
    parser.add_argument('--analysis-dir', type=Path, default=Path('.'))
    parser.add_argument('--peak-files', nargs='+')
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error('Choose a new output directory; existing results are protected')
    if args.gene_min_reps < 1 or args.gene_min_molecules < 1:
        parser.error('Gene support thresholds must be positive integers')
    if normalize_genome(args.genome) != normalize_genome(args.gtf_genome):
        parser.error('GTF genome must match the mapping genome')
    if len(set(args.samples)) != len(args.samples) or not len(args.samples) == len(args.counts) == len(args.headers):
        parser.error('Sample names must be unique and match count/header lists')
    for sample in args.samples:
        validate_label(sample)
    if args.peak_files and len(args.peak_files) != len(args.samples):
        parser.error('Peak files must match the sample list')
    design = json.loads(Path(args.design_file).read_text())
    metadata = {row['sample']: row for row in design['sample_metadata']}
    if set(metadata) != set(args.samples):
        parser.error('Metadata differs from sample list')
    header, sites = load_sites(args.counts, args.headers)
    lengths = canonical_lengths(args.genome)
    for chrom, length in header:
        canonical = 'chr' + chrom.removeprefix('chr')
        if canonical in lengths and length != lengths[canonical]:
            parser.error(f'Reference chromosome length disagrees with {args.genome}: {chrom}')
    annotation = GTFAnnotation(args.gtf, header, args.promoter_upstream, args.promoter_downstream, args.genome)
    genome_filter = None
    if args.environment_file:
        environment = json.loads(Path(args.environment_file).read_text())
        reference = environment['reference']
        if normalize_genome(reference['genome']) != normalize_genome(args.genome):
            parser.error('Mapping environment genome differs from annotation genome')
        genome_filter = GenomeFilter(args.genome, reference['blacklist_bed'], reference['blacklist_genome'])
        if genome_filter.sha256 != environment['genome_filter']['blacklist_sha256']:
            parser.error('Blacklist changed since mapping QC; use the original mapping blacklist')
    elif args.blacklist:
        genome_filter = GenomeFilter(args.genome, args.blacklist, args.genome)
    annotated_chromosomes = {chrom for chrom, gene in annotation.genes}
    occupied = {chrom for sample in sites for chrom, rows in sample.items() if rows}
    if occupied - annotated_chromosomes:
        parser.error('GTF has no gene annotations for occupied chromosomes: ' + ', '.join(sorted(occupied - annotated_chromosomes)))
    args.outdir.mkdir(parents=True)
    rows = []
    total_expected = 0
    for sample, sample_sites in zip(args.samples, sites, strict=True):
        totals = Counter(dict.fromkeys(FEATURES, 0))
        with (args.outdir / f'{sample}.dsb_sites.annotations.tsv').open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['chrom', 'start', 'end', 'molecules', 'feature', 'gene_ids', 'gene_names', 'gene_biotypes'])
            for chrom, positions in sample_sites.items():
                for start, count in positions:
                    category, _, genes = annotation.annotate(chrom, start, start + 1)
                    totals[category] += count
                    writer.writerow([chrom, start, start + 1, count, category, *annotation.gene_fields(chrom, genes)])
        total = sum(count for positions in sample_sites.values() for _, count in positions)
        if sum(totals.values()) != total:
            raise ValueError('Annotation did not conserve molecule totals')
        total_expected += total
        for category in FEATURES:
            rows.append({'sample': sample, **{key: metadata[sample].get(key, '') for key in ('group', 'replicate', 'donor')},
                         'feature': category, 'molecules': totals[category], 'total_molecules': total,
                         'percentage': 100 * totals[category] / total if total else None})
    samples = pd.DataFrame(rows)
    samples.to_csv(args.outdir / 'dsb_feature_distribution.samples.tsv', sep='\t', index=False)
    conditions = samples.groupby(['group', 'feature'], sort=True).agg(
        mean_percentage=('percentage', 'mean'), sd_percentage=('percentage', 'std'),
        biological_samples=('sample', 'size'), nonempty_samples=('percentage', 'count'),
        pooled_molecules=('molecules', 'sum')).reset_index()
    conditions.to_csv(args.outdir / 'dsb_feature_distribution.conditions.tsv', sep='\t', index=False)
    body = samples[samples.feature.isin(['exon', 'intron'])].groupby('sample', sort=False).agg(
        genebody_molecules=('molecules', 'sum'), genebody_percentage=('percentage', lambda x: x.sum(min_count=1))).reset_index()
    body.to_csv(args.outdir / 'dsb_genebody.samples.tsv', sep='\t', index=False)
    gene_summary = write_gene_signal(annotation, args.samples, sites, metadata, args.outdir / 'GeneSignal',
                                    args.gene_min_reps, args.gene_min_molecules, genome_filter=genome_filter)
    (args.outdir / 'annotation_gene_signal_mqc.json').write_text(json.dumps({
        'id': 'helixbusters_gene_signal', 'section_name': 'Helixbusters per-condition gene candidates',
        'description': 'Descriptive rankings. Analysis purpose: ' + design.get('analysis_purpose', 'biological') + '. For input_titration, support refers to a single library per input level, not biological replication. Ambiguous gene assignments excluded; no differential significance test. Lists: Analysis/Annotation/GeneSignal.',
        'plot_type': 'table', 'pconfig': {'id': 'helixbusters_gene_signal_table', 'title': 'Technical input-level gene signal' if design.get('analysis_purpose') == 'input_titration' else 'Replicate-supported gene signal'},
        'data': gene_summary}, indent=2) + '\n')
    for path in sorted(args.analysis_dir.glob('*.counts.tsv')):
        if not (path.name.startswith('windows_') or path.name == 'peaks_consensus.counts.tsv'):
            continue
        def regions():
            with path.open() as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    yield row['region'], row['chrom'], int(row['start']), int(row['end'])
        write_regions(annotation, regions(), args.outdir / path.name.replace('.counts.tsv', '.annotations.tsv'))
    for sample, path in zip(args.samples, args.peak_files or [], strict=bool(args.peak_files)):
        def peaks():
            with Path(path).open() as handle:
                for index, line in enumerate(handle, 1):
                    if line.strip() and not line.startswith('#'):
                        fields = line.split()
                        yield fields[3] if len(fields) > 3 else f'peak_{index}', fields[0], int(fields[1]), int(fields[2])
        write_regions(annotation, peaks(), args.outdir / f'{sample}.peaks.annotations.tsv')
    for path in sorted(args.analysis_dir.glob('MergedReplicate/*/peaks/*.consensus.bed')):
        def consensus():
            with path.open() as handle:
                for line in handle:
                    fields = line.split()
                    yield fields[3], fields[0], int(fields[1]), int(fields[2])
        write_regions(annotation, consensus(), args.outdir / f'{path.stem}.annotations.tsv')
    colors = dict(zip(FEATURES, ['#e69f00', '#0072b2', '#56b4e9', '#999999']))
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    base = np.zeros(len(args.samples))
    for feature in FEATURES:
        values = samples[samples.feature == feature].set_index('sample').loc[args.samples, 'percentage'].fillna(0).to_numpy()
        axes[0].bar(np.arange(len(args.samples)), values, bottom=base, color=colors[feature], label=feature)
        base += values
    axes[0].set_xticks(np.arange(len(args.samples)), args.samples, rotation=40, ha='right')
    axes[0].set_title('DSB molecule allocation per sample')
    groups = sorted(samples['group'].unique())
    for g, group in enumerate(groups):
        for f, feature in enumerate(FEATURES):
            subset = samples[(samples.group == group) & (samples.feature == feature)]
            x = g * (len(FEATURES) + 1) + f
            axes[1].bar(x, subset.percentage.mean(), color=colors[feature], alpha=.65)
            values = subset.percentage.dropna().to_numpy()
            offsets = np.linspace(-.12, .12, len(values)) if len(values) > 1 else np.zeros(len(values))
            axes[1].scatter(x + offsets, values, s=25, color='black', zorder=3)
    axes[1].set_xticks([g * 5 + 1.5 for g in range(len(groups))], groups)
    axes[1].set_title('Condition means with biological sample points')
    for ax in axes:
        ax.set_ylabel('Percent of retained DSB molecules')
        ax.set_ylim(0, 100)
        ax.spines[['top', 'right']].set_visible(False)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, .93))
    fig.savefig(args.outdir / 'dsb_feature_distribution.png', dpi=180)
    fig.savefig(args.outdir / 'dsb_feature_distribution.pdf')
    plt.close(fig)
    provenance = {**annotation.provenance, 'design': design['design'], 'total_molecules': total_expected,
        'sample_order': args.samples, 'denominator': 'retained deduplicated DSB molecules, counted once per sample',
        'condition_aggregation': 'equal-weight mean of nonempty sample percentages; points represent biological samples; no significance test',
        'enrichment': 'not computed; percentages are not enrichment over genomic opportunity or mappability',
        'versions': {'numpy': np.__version__, 'pandas': pd.__version__, 'matplotlib': matplotlib.__version__}}
    (args.outdir / 'annotation.provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    for label, data in [('samples', {sample: {feature: rows.loc[rows.feature == feature, 'percentage'].iloc[0]
                            for feature in FEATURES} for sample, rows in samples.groupby('sample')}),
                        ('conditions', {group: {feature: rows.loc[rows.feature == feature, 'mean_percentage'].iloc[0]
                            for feature in FEATURES} for group, rows in conditions.groupby('group')})]:
        # JSON null represents empty-library percentages rather than a false zero.
        data = {key: {feature: float(value) if pd.notna(value) else None for feature, value in values.items()}
                for key, values in data.items()}
        (args.outdir / f'annotation_{label}_mqc.json').write_text(json.dumps({
            'id': f'helixbusters_annotation_{label}', 'section_name': f'Helixbusters DSB features: {label}',
            'description': 'Exclusive promoter > exon > intron > intergenic; molecule-weighted sample percentages. Conditions are equal-weight sample means. Not enrichment. Figures and gene associations: Analysis/Annotation.',
            'plot_type': 'bargraph', 'pconfig': {'id': f'helixbusters_annotation_{label}_plot',
                'title': f'DSB feature allocation: {label}', 'ylab': 'Molecules (%)', 'stacking': 'normal'},
            'data': data}, indent=2, allow_nan=False) + '\n')
    print(f'Annotated {total_expected:,} molecules across {len(args.samples)} samples; outputs: {args.outdir}')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
