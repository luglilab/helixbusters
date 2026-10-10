"""Verify PCA sample identities, deterministic depth sensitivity and skips."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).parents[1]

class TestWindowPCA(unittest.TestCase):
    def run_pca(self, directory, counts, label):
        root = Path(directory)
        inputs = root / 'inputs'
        inputs.mkdir(exist_ok=True)
        names = [f's{i}' for i in range(counts.shape[0])]
        pd.DataFrame({'sample': names, 'group': ['A' if i % 2 else 'B' for i in range(len(names))],
                      'replicate': names, 'donor': [''] * len(names)}).to_csv(inputs / 'analysis.samples.tsv', sep='\t', index=False)
        table = pd.DataFrame(counts.T, columns=names)
        table.insert(0, 'end', np.arange(counts.shape[1]) * 10000 + 10000)
        table.insert(0, 'start', np.arange(counts.shape[1]) * 10000)
        table.insert(0, 'chrom', 'chr1')
        table.insert(0, 'region', [f'r{i}' for i in range(counts.shape[1])])
        table.to_csv(inputs / 'windows_10000.counts.tsv', sep='\t', index=False)
        out = root / label
        result = subprocess.run([sys.executable, str(ROOT / 'scripts/pca_windows.py'),
            '--analysis-dir', str(inputs), '--outdir', str(out), '--windows', '10000'],
            capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        return out, json.loads((out / 'pca.summary.json').read_text())

    def test_samples_normalization_and_reproducibility(self):
        counts = np.random.default_rng(4).poisson(np.array([1, 2, 4, 6])[:, None], size=(4, 30))
        with tempfile.TemporaryDirectory() as directory:
            a, summary = self.run_pca(directory, counts, 'first')
            b, _ = self.run_pca(directory, counts, 'second')
            self.assertEqual((a / 'pca.coordinates.tsv').read_bytes(), (b / 'pca.coordinates.tsv').read_bytes())
            coords = pd.read_csv(a / 'pca.coordinates.tsv', sep='\t')
            self.assertEqual(len(coords), 8)
            self.assertEqual(set(coords['sample']), {'s0', 's1', 's2', 's3'})
            window = summary['windows']['10000']
            self.assertEqual(window['equal_depth_molecules'], int(counts.sum(axis=1).min()))
            self.assertEqual(list(window['library_molecules'].values()), counts.sum(axis=1).tolist())
            self.assertGreater((a / 'windows_PCA.pdf').stat().st_size, 100)
            self.assertEqual(json.loads((a / 'pca_mqc.json').read_text())['plot_type'], 'table')
            diagnostics = json.loads((a / 'pca_mqc.json').read_text())['data']['10000_bp_logCPM']
            self.assertEqual(diagnostics['minimum_library_molecules'], int(counts.sum(axis=1).min()))
            self.assertEqual(diagnostics['maximum_library_molecules'], int(counts.sum(axis=1).max()))
            self.assertIsInstance(diagnostics['depth_association_flag'], bool)

    def test_insufficient_samples_and_empty_library_are_reported(self):
        for counts, reason in [(np.ones((2, 5), dtype=int), 'Fewer than three samples'),
                               (np.zeros((3, 5), dtype=int), 'At least one sample has no molecules')]:
            with self.subTest(reason=reason), tempfile.TemporaryDirectory() as directory:
                out, summary = self.run_pca(directory, counts, 'skipped')
                self.assertEqual(summary['windows']['10000']['reason'], reason)
                self.assertEqual(len(pd.read_csv(out / 'pca.coordinates.tsv', sep='\t')), 0)

    def test_no_variation_is_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            _, summary = self.run_pca(directory, np.full((3, 5), 3), 'skipped')
            self.assertEqual(summary['windows']['10000']['status'], 'skipped')

if __name__ == '__main__':
    unittest.main()
