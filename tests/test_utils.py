from pathlib import Path
import tempfile
import unittest

from helixbusters.utils import read_excel_column, process_bam_and_generate_umi_outputs
from tests.bam_helpers import read, write_bam


class TestUtils(unittest.TestCase):
    def test_single_end_samplesheet(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'samples.csv'
            path.write_text('Sample,Replicate,Group,PathReadForward,SampleBarcodeForward\n'
                            'S1,R1,Acute,/reads/S1.fastq.gz,CGTGTGAG\n')
            frame, modality = read_excel_column(path)
            self.assertEqual(modality, 'single-end')
            self.assertEqual(frame.loc[0, 'Sample'], 'S1')

    def test_existing_three_argument_bam_api(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam = write_bam(root / 'input.bam', [read(flag=16), read(flag=16)])
            stats = process_bam_and_generate_umi_outputs(bam, root / 'families.tsv', root / 'counts.bed')
            self.assertEqual(stats['deduplicated_molecules'], 1)
            self.assertEqual((root / 'counts.bed').read_text(), 'chr1\t119\t120\t1\n')
