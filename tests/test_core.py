import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from helixbusters.core import Helixbusters
from tests.bam_helpers import read, write_bam


class TestHelixbusters(unittest.TestCase):
    def test_samples_are_independent_and_options_reach_deduplication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hb = Helixbusters('unused.csv', 'human', 1, 'unused_index')
            hb.infofile = pd.DataFrame([
                {'Sample': sample, 'OutputPath': str(root / sample),
                 'BamFilteredPath': str(write_bam(root / f'{sample}.bam',
                                                 [read()] * 5 + [read(umi='AAAAAAAC')]))}
                for sample in ('acute', 'chronic')
            ])
            hb.generate_umi_output_for_samples(method='directional')
            for _, row in hb.infofile.iterrows():
                qc = json.loads(Path(row['Deduplication_QC']).read_text())
                self.assertEqual(qc['deduplicated_molecules'], 1)
                self.assertEqual(qc['accepted_reads'], 6)
                self.assertEqual(qc['parameters']['method'], 'directional')
                self.assertTrue(Path(row['Sites_Output']).is_file())
                self.assertTrue(Path(row['Molecules_Output']).is_file())

    def test_missing_bam_is_an_error(self):
        hb = Helixbusters('unused.csv', 'human', 1, 'unused_index')
        hb.infofile = pd.DataFrame([{'Sample': 'sample', 'OutputPath': '.', 'BamFilteredPath': None}])
        with self.assertRaises(FileNotFoundError):
            hb.generate_umi_output_for_samples()

    def test_requires_sample_information(self):
        hb = Helixbusters('unused.csv', 'human', 1, 'unused_index')
        with self.assertRaises(ValueError):
            hb.generate_umi_output_for_samples()
