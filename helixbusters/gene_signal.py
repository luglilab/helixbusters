"""Replicate-aware descriptive gene signal lists from uniquely assigned DSB ends."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import numpy as np
import pandas as pd
from helixbusters.reporting import validate_label

CONTEXTS = ('promoter', 'gene_body', 'combined')


def write_gene_signal(annotation, samples, sites, metadata, destination, minimum_replicates=2, minimum_molecules=2):
    if minimum_replicates < 1 or minimum_molecules < 1:
        raise ValueError('Gene support thresholds must be positive integers')
    destination = Path(destination)
    destination.mkdir()
    genes = sorted(annotation.genes)
    lengths = defaultdict(Counter)
    for chrom, segments in annotation.segments.items():
        for segment, candidates in zip(segments, annotation.gene_candidates[chrom], strict=True):
            if len(candidates) == 1:
                start, end, category, _ = segment
                key = (chrom, candidates[0])
                context = 'promoter' if category == 'promoter' else 'gene_body'
                lengths[key][context] += end - start
                lengths[key]['combined'] += end - start
    counts = {context: defaultdict(Counter) for context in CONTEXTS}
    audits = []
    with (destination / 'ambiguous_gene_sites.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'chrom', 'start', 'end', 'molecules', 'feature', 'candidate_gene_ids'])
        for sample, sample_sites in zip(samples, sites, strict=True):
            audit = Counter(unique_gene_molecules=0, ambiguous_gene_molecules=0, intergenic_molecules=0)
            for chrom, rows in sample_sites.items():
                for position, molecules in rows:
                    category, candidates = annotation.site_gene_candidates(chrom, position)
                    if not candidates:
                        audit['intergenic_molecules'] += molecules
                    elif len(candidates) > 1:
                        audit['ambiguous_gene_molecules'] += molecules
                        writer.writerow([sample, chrom, position, position + 1, molecules, category, ';'.join(candidates)])
                    else:
                        audit['unique_gene_molecules'] += molecules
                        key = (chrom, candidates[0])
                        context = 'promoter' if category == 'promoter' else 'gene_body'
                        counts[context][key][sample] += molecules
                        counts['combined'][key][sample] += molecules
            total = sum(count for rows in sample_sites.values() for _, count in rows)
            if sum(audit.values()) != total:
                raise ValueError('Gene assignment does not conserve sample molecules')
            audits.append({'sample': sample, 'total_molecules': total, **audit})
    pd.DataFrame(audits).to_csv(destination / 'gene_assignment.samples.tsv', sep='\t', index=False)
    totals = np.array([row['total_molecules'] for row in audits], dtype=float)
    info = pd.DataFrame([{'chrom': chrom, 'gene_id': gene_id, 'gene_name': annotation.genes[(chrom, gene_id)]['name'],
        'gene_biotype': annotation.genes[(chrom, gene_id)]['biotype'], 'strand': annotation.genes[(chrom, gene_id)]['strand'],
        'gene_start': annotation.genes[(chrom, gene_id)]['start'], 'gene_end': annotation.genes[(chrom, gene_id)]['end']}
        for chrom, gene_id in genes])
    groups = defaultdict(list)
    for i, sample in enumerate(samples):
        validate_label(metadata[sample]['group'])
        groups[metadata[sample]['group']].append(i)
    summary = {}
    for context in CONTEXTS:
        raw = np.array([[counts[context][gene].get(sample, 0) for sample in samples] for gene in genes], dtype=np.int64)
        cpm = np.divide(raw * 1e6, totals[None, :], out=np.full(raw.shape, np.nan), where=totals[None, :] > 0)
        bp = np.array([lengths[gene][context] for gene in genes], dtype=np.int64)
        features = info.copy()
        features['uniquely_assignable_bp'] = bp
        pd.concat([features, pd.DataFrame(raw, columns=samples)], axis=1).to_csv(
            destination / f'genes.{context}.counts.tsv', sep='\t', index=False)
        pd.concat([features, pd.DataFrame(cpm, columns=samples)], axis=1).to_csv(
            destination / f'genes.{context}.CPM.tsv', sep='\t', index=False)
        for group, columns in sorted(groups.items()):
            nonempty = [i for i in columns if totals[i] > 0]
            table = features.copy()
            subset = raw[:, columns]
            table['condition'] = group
            table['biological_samples'] = len(columns)
            table['nonempty_samples'] = len(nonempty)
            table['replicates_with_signal'] = (subset > 0).sum(axis=1)
            table['replicates_meeting_min_molecules'] = (subset >= minimum_molecules).sum(axis=1)
            table['pooled_molecules'] = subset.sum(axis=1)
            if nonempty:
                normalized = cpm[:, nonempty]
                table['mean_CPM'] = normalized.mean(axis=1)
                table['median_CPM'] = np.median(normalized, axis=1)
                table['sd_CPM'] = normalized.std(axis=1, ddof=1) if len(nonempty) > 1 else np.nan
            else:
                table['mean_CPM'] = table['median_CPM'] = table['sd_CPM'] = np.nan
            table['mean_CPM_per_kb'] = np.divide(table['mean_CPM'].to_numpy() * 1000, bp,
                out=np.full(len(bp), np.nan), where=bp > 0)
            table['supported_candidate'] = table['replicates_meeting_min_molecules'] >= minimum_replicates
            for index in columns:
                table[f'{samples[index]}_molecules'] = raw[:, index]
                table[f'{samples[index]}_CPM'] = cpm[:, index]
            # Rank within each condition without comparing or selecting by another condition.
            table = table[table.pooled_molecules > 0].sort_values(
                ['supported_candidate', 'median_CPM', 'mean_CPM', 'gene_id', 'chrom'],
                ascending=[False, False, False, True, True], kind='stable')
            table.insert(0, 'rank', np.arange(1, len(table) + 1))
            table.to_csv(destination / f'{group}.{context}.ranked_genes.tsv', sep='\t', index=False)
            candidates = table[table.supported_candidate]
            candidates.to_csv(destination / f'{group}.{context}.candidate_genes.tsv', sep='\t', index=False)
            summary[f'{group}_{context}'] = {'genes_with_signal': len(table), 'supported_candidates': len(candidates),
                'biological_samples': len(columns), 'minimum_replicates': minimum_replicates,
                'minimum_molecules_per_replicate': minimum_molecules}
    provenance = {'status': 'descriptive per-condition gene ranking; no differential model or significance test',
        'gene_assignment': 'global promoter > exon > intron priority; genes overlapping the highest-priority feature only; multiple candidates excluded from integer gene counts and reported separately',
        'contexts': {'promoter': 'promoter-assigned sites', 'gene_body': 'exon/intron-assigned sites excluding all promoters',
                     'combined': 'promoter plus gene_body; alternative summaries, not independent datasets'},
        'normalization': 'CPM denominator is all retained deduplicated molecules per sample, including ambiguous and intergenic molecules; zero-depth samples have undefined CPM',
        'ranking': 'supported candidates first, then median sample CPM, mean sample CPM and stable gene ID tie-break; equal sample weighting',
        'lengths': 'uniquely assignable feature bp from GTF partition; not corrected for blacklist or mappability; CPM per kb is descriptive density, not enrichment',
        'gtf_sha256': annotation.provenance['gtf_sha256'], 'sample_metadata': [metadata[sample] for sample in samples],
        'conditions': summary}
    (destination / 'gene_signal.provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    (destination / 'README.md').write_text('''# Per-condition gene signal candidates

Use CONDITION.promoter.candidate_genes.tsv and CONDITION.gene_body.candidate_genes.tsv
as reproducible descriptive candidates. ranked_genes.tsv retains all genes with any
signal and includes support flags. combined is an additional promoter-plus-body
summary. Wide integer counts/CPM matrices retain all annotated genes, including zeros.
All files retain gene IDs, names, biotypes, coordinates, assignable length, library
normalization and sample-level support. Ambiguous gene assignments are audited and
excluded from integer gene counts. Assignment does not establish functional targets.

Candidates meet configured minimum molecules in configured biological replicates.
This is an exploratory support filter, not statistical significance. A condition
with fewer available replicates than required produces an empty candidate list;
it does not silently relax the threshold. CPM ranks relative molecular allocation,
not absolute DSB burden. Long genes have more opportunity to accumulate DSBs;
CPM per kb provides descriptive density but lacks mappability/blacklist correction.
Differential comparisons need biological design, replicate-level integer counts,
appropriate normalization, an abundance filter and multiple-testing correction.
''')
    return summary
