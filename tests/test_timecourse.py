"""Verify paired time metadata and a real all-time-point DESeq2 fit."""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
import pandas as pd
spec = importlib.util.spec_from_file_location('timecourse', Path(__file__).resolve().parents[1] / 'scripts/timecourse_dsb.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
TIMES = ['Baseline', 'T2H', 'T6H', 'T24H']

def metadata():
    return pd.DataFrame([dict(sample=f'{d}_{t}', donor=d, replicate=d, group=t) for d in ['D1', 'D2', 'D3'] for t in TIMES])

class TestTimecourse(unittest.TestCase):

    def test_metadata_rejects_incomplete_or_unreplicated_design(self):
        self.assertEqual(len(module.validate_timecourse(metadata(), TIMES)), 12)
        for rows, times in [(metadata().iloc[1:], TIMES), (metadata().iloc[:8], TIMES), (metadata(), TIMES[:3]), (metadata(), ['Baseline', 'T2H', 'T2H'])]:
            with self.assertRaises(ValueError):
                module.validate_timecourse(rows, times)

    def test_real_fit_joint_correction_and_global_effect_semantics(self):
        check = subprocess.run(['Rscript', '--vanilla', '-e', 'quit(status=if(requireNamespace("DESeq2",quietly=TRUE)) 0 else 1)'], capture_output=True)
        if check.returncode:
            self.skipTest('DESeq2 unavailable')
        rows = metadata().sort_values('sample').reset_index(drop=True)
        rng = np.random.default_rng(1729)
        mean = rng.uniform(60, 200, (150, 1)) * np.ones((150, 12))
        late = np.flatnonzero(rows.group.to_numpy() == 'T24H')
        mean[:15, late] *= 8
        counts = rng.negative_binomial(50, 50 / (50 + mean))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            analysis = root / 'analysis'
            analysis.mkdir()
            rows.to_csv(analysis / 'analysis.samples.tsv', sep='\t', index=False)
            totals = counts.sum(axis=0)
            (analysis / 'analysis.summary.json').write_text(json.dumps({'samples': {n: {'molecules': int(t)} for n, t in zip(rows['sample'], totals)}}))
            table = pd.DataFrame(counts, columns=rows['sample'])
            table.insert(0, 'region', [f'r{i}' for i in range(150)])
            table.to_csv(analysis / 'windows_10000.counts.tsv', sep='\t', index=False)
            report = module.run(analysis, root / 'output', TIMES)
            folder = root / 'output/windows_10000'
            self.assertEqual(report['families']['windows_10000']['status'], 'fitted')
            global_table = pd.read_csv(folder / 'poscounts.global.tsv', sep='\t')
            self.assertNotIn('log2FoldChange', global_table)
            late_table = pd.read_csv(folder / 'poscounts.T24H_vs_Baseline.tsv', sep='\t')
            self.assertGreater(late_table.iloc[:15].log2FoldChange.median(), 1)
            self.assertGreater(global_table.iloc[:15].significant.sum(), 5)
            self.assertTrue((late_table.padj >= late_table.padj_within_contrast - 1e-12).all())
            self.assertIn('D1_log2_CPM_ratio', late_table)
            self.assertTrue((folder / 'library_total.global.tsv').exists())
            self.assertEqual(len(pd.read_csv(folder / 'poscounts.size_factors.tsv', sep='\t')), 12)
