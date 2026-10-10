"""Check biological units, depth conservation and a real paired DESeq2 contrast."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from helixbusters.differential import contrast_metadata, depth_sensitivity, integer_counts, paired_effects, run_comparison


def metadata():
    return pd.DataFrame([{'sample': f'{donor}_{group}', 'group': group, 'replicate': donor,
                          'donor': donor} for donor in ('D1', 'D2', 'D3') for group in ('A', 'B')])


class TestDifferential(unittest.TestCase):
    def test_explicit_design_complete_pairs_and_minimum_replicates(self):
        rows = metadata()
        self.assertEqual(len(contrast_metadata(rows, 'paired', 'B', 'A')), 6)
        for design in ('unspecified', 'unpaired'):
            with self.assertRaises(ValueError):
                contrast_metadata(rows, design, 'B', 'A')
        independent = rows.copy()
        independent['donor'] = independent['sample']
        self.assertEqual(len(contrast_metadata(independent, 'unpaired', 'B', 'A')), 6)
        rows.loc[0, 'donor'] = ''
        with self.assertRaises(ValueError):
            contrast_metadata(rows, 'paired', 'B', 'A')
        with self.assertRaises(ValueError):
            contrast_metadata(metadata().iloc[:4], 'paired', 'B', 'A')
        with self.assertRaises(ValueError):
            contrast_metadata(metadata(), 'paired', 'A', 'A')

    def test_integer_counts_and_equal_depth_residual(self):
        with self.assertRaises(ValueError):
            integer_counts(pd.DataFrame({'a': [1.1]}), ['a'])
        raw = np.array([[4, 40], [6, 60]])
        totals = np.array([10, 100])
        a, draws, target = depth_sensitivity(raw, totals, {'A': [0, 1]}, 10, 42)
        b, _, _ = depth_sensitivity(raw, totals, {'A': [0, 1]}, 10, 42)
        np.testing.assert_array_equal(a['A'], b['A'])
        self.assertEqual(target, 10)
        # When both features have support, the total matched-depth molecule count
        # makes supporting all features impossible if the threshold exceeds depth.
        residual = np.array([[0, 0]])
        frequencies, _, _ = depth_sensitivity(residual, totals, {'A': [0, 1]}, 10, 42)
        np.testing.assert_array_equal(frequencies['A'], [0])
        with self.assertRaises(ValueError):
            depth_sensitivity(raw, np.array([1, 10]), {'A': [0, 1]}, 10, 42)

    def test_paired_direction_and_leave_one_donor_influence(self):
        values = np.array([[10, 20, 10, 20, 10, 20], [10, 11, 10, 9, 10, 100]])
        donors, ratios, median, concordant, omitted, stable = paired_effects(values, metadata(), 'B', 'A')
        self.assertEqual(donors, ['D1', 'D2', 'D3'])
        self.assertTrue(concordant[0]); self.assertTrue(stable[0])
        self.assertFalse(concordant[1]); self.assertFalse(stable[1])
        self.assertTrue((ratios[0] > 0).all())

    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_real_models_direction_metadata_and_protected_outputs(self):
        check = subprocess.run(['Rscript', '--vanilla', '-e',
                                'quit(status=if (requireNamespace("DESeq2", quietly=TRUE)) 0 else 1)'],
                               capture_output=True, timeout=90)
        if check.returncode:
            self.skipTest('DESeq2 unavailable')
        rng = np.random.default_rng(19)
        rows = metadata()
        n = 250
        mean = rng.uniform(30, 150, size=(n, 1)) * np.array([[.8, .8, 1, 1, 1.2, 1.2]])
        mean[:15, [1, 3, 5]] *= 8
        mean[15:30, [1, 3, 5]] *= .125
        raw = rng.negative_binomial(30, 30 / (30 + mean)).astype(np.int64)
        # A filtered low-abundance row remains in the complete reported universe.
        raw[-1] = 0
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); analysis = root / 'Analysis'; analysis.mkdir()
            source_rows = rows.copy(); source_rows['donor'] = ''
            source_rows.to_csv(analysis / 'analysis.samples.tsv', sep='\t', index=False)
            rows.to_csv(root / 'confirmed_metadata.tsv', sep='\t', index=False)
            totals = raw.sum(axis=0)
            (analysis / 'analysis.summary.json').write_text(json.dumps({'samples': {
                name: {'molecules': int(depth)} for name, depth in zip(rows['sample'], totals)}}))
            table = pd.DataFrame({'region': [f'r{i}' for i in range(n)], 'chrom': 'chr1',
                                  'start': np.arange(n) * 10000, 'end': (np.arange(n) + 1) * 10000})
            table = pd.concat([table, pd.DataFrame(raw, columns=rows['sample'])], axis=1)
            table.to_csv(analysis / 'windows_10000.counts.tsv', sep='\t', index=False)
            script = Path(__file__).parents[1] / 'scripts/differential_dsb.py'
            command = subprocess.run([sys.executable, str(script), '--analysis-dir', str(analysis),
                                      '--outdir', str(root / 'output'), '--metadata', str(root / 'confirmed_metadata.tsv'),
                                      '--design', 'paired', '--contrast', 'B', 'A', '--iterations', '3'],
                                     capture_output=True, text=True, timeout=120)
            self.assertEqual(command.returncode, 0, command.stdout + command.stderr)
            summary = json.loads((root / 'output/differential.summary.json').read_text())
            result = pd.read_csv(root / 'output/windows_10000/results.tsv', sep='\t').set_index('region')
            self.assertEqual(summary['model_formula'], '~ donor + condition')
            self.assertGreater(result.loc[[f'r{i}' for i in range(15)], 'log2FoldChange'].median(), 2)
            self.assertLess(result.loc[[f'r{i}' for i in range(15, 30)], 'log2FoldChange'].median(), -2)
            self.assertTrue(pd.isna(result.loc[f'r{n-1}', 'pvalue']))
            self.assertEqual(result.loc[f'r{n-1}', 'test_status'], 'below_abundance_filter')
            self.assertEqual(len(result), n)
            np.testing.assert_array_equal(result.loc[table.region, rows['sample']].to_numpy(), raw)
            self.assertIn('library_total_padj', result)
            self.assertIn('omit_D1_median_log2_CPM_ratio', result)
            self.assertTrue((root / 'output/windows_10000/dispersion_fit.pdf').is_file())
            # Exercise independent samples and gene-only input without inventing pairing.
            gene_analysis = root / 'GeneAnalysis'; gene_analysis.mkdir()
            source_rows.to_csv(gene_analysis / 'analysis.samples.tsv', sep='\t', index=False)
            shutil.copy2(analysis / 'analysis.summary.json', gene_analysis / 'analysis.summary.json')
            gene_folder = gene_analysis / 'Annotation/GeneSignal'; gene_folder.mkdir(parents=True)
            genes = pd.DataFrame({'chrom': 'chr1', 'gene_id': [f'g{i}' for i in range(n)],
                                  'gene_name': [f'gene{i}' for i in range(n)],
                                  'effective_assignable_bp': np.arange(1, n + 1) * 1000})
            pd.concat([genes, pd.DataFrame(raw, columns=rows['sample'])], axis=1).to_csv(
                gene_folder / 'genes.gene_body.counts.tsv', sep='\t', index=False)
            independent = rows.copy(); independent['donor'] = independent['sample']
            independent.to_csv(root / 'independent.tsv', sep='\t', index=False)
            command = subprocess.run([sys.executable, str(script), '--analysis-dir', str(gene_analysis),
                                      '--outdir', str(root / 'independent'), '--metadata', str(root / 'independent.tsv'),
                                      '--design', 'unpaired', '--contrast', 'B', 'A', '--iterations', '3'],
                                     capture_output=True, text=True, timeout=120)
            self.assertEqual(command.returncode, 0, command.stdout + command.stderr)
            independent_result = pd.read_csv(root / 'independent/genes.gene_body/results.tsv', sep='\t')
            self.assertNotIn('all_donors_same_direction', independent_result)
            np.testing.assert_allclose(independent_result.B_mean_CPM_per_effective_kb,
                                       independent_result.B_mean_CPM * 1000 / independent_result.effective_assignable_bp)
            original = (analysis / 'windows_10000.counts.tsv').read_bytes()
            with self.assertRaises(ValueError):
                run_comparison(analysis, root / 'output', rows, 'paired', 'B', 'A', iterations=3)
            self.assertEqual((analysis / 'windows_10000.counts.tsv').read_bytes(), original)
            changed = rows.copy(); changed.loc[0, 'group'] = 'B'
            with self.assertRaises(ValueError):
                run_comparison(analysis, root / 'bad_metadata', changed, 'paired', 'B', 'A', iterations=3)
            rows.to_csv(analysis / 'analysis.samples.tsv', sep='\t', index=False)
            changed = rows.copy(); changed.loc[changed.donor == 'D1', 'donor'] = 'different'
            with self.assertRaises(ValueError):
                run_comparison(analysis, root / 'bad_donors', changed, 'paired', 'B', 'A', iterations=3)


if __name__ == '__main__':
    unittest.main()
