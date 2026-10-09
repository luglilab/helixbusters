import gzip
import importlib.util
import json
import shutil
import subprocess
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location('prepare_bliss_reads', Path(__file__).parents[1] / 'scripts/prepare_bliss_reads.py')
preparation = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(preparation)


class TestBlissPreparation(unittest.TestCase):
    def test_tolerant_prefix_excludes_substitution_insertion_and_deletion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            motif = 'CCCTATAGTGAGTCGTAT'
            genomic = 'GATTACA' * 5
            inserts = [motif + genomic, 'AC' + motif[2:] + genomic,
                       motif[:8] + 'A' + motif[8:] + genomic,
                       motif[:8] + motif[9:] + genomic, genomic]
            source = root / 'input.fastq'
            with source.open('w') as handle:
                for i, insert in enumerate(inserts):
                    seq = 'AACCGGTTCGTGTGAG' + insert
                    handle.write(f'@r{i}\n{seq}\n+\n{"I" * len(seq)}\n')
            counts = preparation.prepare(source, 'test', 'CTCACACG', 'reverse_complement',
                                         8, 20, motif, root / 'out', 2)
            self.assertEqual(counts['technical_prefix_reads'], 4)
            self.assertEqual(counts['output_reads'], 1)
            with gzip.open(root / 'out/test.umi.fastq.gz', 'rt') as handle:
                self.assertEqual(handle.read().splitlines()[1], genomic)

    @unittest.skipUnless(shutil.which('multiqc'), 'MultiQC is not on PATH')
    def test_multiqc_combines_preparation_reports_from_two_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / 'input.fastq'
            seq = 'AACCGGTTCGTGTGAG' + 'GATTACA' * 5
            path.write_text(f'@r\n{seq}\n+\n{"I" * len(seq)}\n')
            for sample in ('sample1', 'sample2'):
                preparation.prepare(path, sample, 'CTCACACG', 'reverse_complement', 8, 20, '', root)
            result = subprocess.run(['multiqc', '.', '--filename', 'multiqc_report.html', '--outdir', '.',
                                     '--data-dir', '--cl-config', 'data_dir_name: multiqc_data',
                                     '--no-version-check'], cwd=root, capture_output=True, text=True, timeout=90)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('Found 2 samples', result.stderr)
            self.assertTrue((root / 'multiqc_data').is_dir())
            self.assertIn('Helixbusters read preparation', (root / 'multiqc_report.html').read_text())

    def test_orientation_trimming_exclusions_and_original_preservation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / 'original.fastq.gz'
            umi, barcode, observed = 'AACCGGTT', 'CTCACACG', 'CGTGTGAG'
            genomic = 'GATTACAGATTACAGATTACA'
            technical = 'CCCTATAGTGAGTCGTAT'
            seqs = [umi + observed + genomic, umi + barcode + genomic,
                    'NACCGGTT' + observed + genomic, umi + observed + technical + genomic,
                    umi + observed + 'ACG', umi + observed + genomic + technical]
            with gzip.open(path, 'wt') as handle:
                for i, seq in enumerate(seqs):
                    handle.write(f'@read{i} description\n{seq}\n+\n{"I" * len(seq)}\n')
            original = path.read_bytes()
            counts = preparation.prepare(path, 'sample', barcode, 'reverse_complement', 8, 20, technical, root / 'out')
            self.assertEqual(counts, {'input_reads': 6, 'barcode_mismatch_reads': 1, 'invalid_umi_reads': 1,
                                     'technical_prefix_reads': 1, 'short_insert_reads': 1, 'output_reads': 2})
            with gzip.open(root / 'out/sample.umi.fastq.gz', 'rt') as handle:
                lines = handle.read().splitlines()
            self.assertEqual(lines[:4], ['@read0_AACCGGTT description', genomic, '+', 'I' * len(genomic)])
            self.assertEqual(lines[5], genomic + technical)  # Internal motif is not clipped or rejected.
            self.assertEqual(path.read_bytes(), original)
            qc = json.loads((root / 'out/sample.preparation.json').read_text())
            self.assertEqual(qc['parameters']['observed_barcode'], observed)
            with self.assertRaisesRegex(ValueError, 'already exist'):
                preparation.prepare(path, 'sample', barcode, 'reverse_complement', 8, 20, technical, root / 'out')

    def test_zero_retention_reports_failure_and_malformed_record_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / 'input.fastq'
            path.write_text('@r\nAAAAAAAAACGT\n+\nIIIIIIIIIIII\n')
            with self.assertRaisesRegex(ValueError, 'No reads retained'):
                preparation.prepare(path, 's', 'CTCACACG', 'forward', 8, 20, '', root / 'zero')
            self.assertTrue((root / 'zero/s.preparation.json').exists())
            path.write_text('@r\nACGT\n+\nIII\n')
            with self.assertRaisesRegex(ValueError, 'Malformed'):
                preparation.prepare(path, 's', 'ACGT', 'forward', 8, 20, '', root / 'bad')
