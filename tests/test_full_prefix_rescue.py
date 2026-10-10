"""Preserve production preparation and prove shared molecules are jointly deduplicated."""
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from prepare_bliss_reads import prepare
from run_full_prefix_rescue import prepare_joint
import run_full_prefix_rescue
from helixbusters.technical import T7_FORWARD, T7_REVERSE


class TestFullPrefixRescue(unittest.TestCase):
    def test_baseline_is_identical_to_existing_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root / 'raw.fastq.gz'
            barcode = 'ACGTACGT'; dna = 'G' * 60
            records = [('keep', 'AAAAAAAA' + barcode + dna), ('rescue', 'CCCCCCCC' + barcode + T7_REVERSE + dna),
                       ('repeat', 'AAAAAAAA' + barcode + T7_REVERSE + T7_FORWARD + dna),
                       ('barcode', 'AAAAAAAA' + 'T' * 8 + dna), ('umi', 'NAAAAAAA' + barcode + dna),
                       ('short', 'AAAAAAAA' + barcode + 'G'),
                       ('edited', 'AAAAAAAA' + barcode + 'A' + T7_REVERSE[1:] + dna)]
            with gzip.open(source, 'wt') as handle:
                for name, seq in records:
                    handle.write(f'@{name} comment\n{seq}\n+\n' + 'I' * len(seq) + '\n')
            golden = root / 'golden'; candidate = root / 'candidate'; candidate.mkdir()
            production = prepare(source, 'toy', barcode, 'reverse_complement', 8, 20, T7_REVERSE[:18], golden, 2)
            result = prepare_joint(source, 'toy', barcode, 'reverse_complement', candidate)
            self.assertEqual(result['production_counts'], production)
            self.assertEqual(result['rescued_reads'], 1)
            self.assertEqual(result['rescue_exclusions']['residual_T7_in_either_orientation'], 1)
            with gzip.open(golden / 'toy.umi.fastq.gz', 'rt') as a, gzip.open(candidate / 'baseline.fastq.gz', 'rt') as b:
                self.assertEqual(a.read(), b.read())
            with gzip.open(candidate / 'candidate_trim.fastq.gz', 'rt') as handle:
                text = handle.read(); self.assertIn('@rescue_CCCCCCCC comment\n' + dna, text)

    @unittest.skipUnless(all(shutil.which(x) for x in ('bwa', 'samtools')), 'Mapping tools unavailable')
    def test_joint_deduplication_does_not_add_same_umi_molecule_twice(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rng = random.Random(47); genome = ''.join(rng.choices('ACGT', k=12000))
            fasta = root / 'ref.fa'; fasta.write_text('>chr1\n' + genome + '\n')
            subprocess.run(['bwa', 'index', str(fasta)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            rows = []
            for group in ('ACUTE', 'CHRONIC'):
                source = root / (group + '.fastq.gz'); insert = genome[100:250]
                with gzip.open(source, 'wt') as handle:
                    for name, body in [('original', insert), ('rescued', T7_REVERSE + insert)]:
                        seq = 'AAAAAAAAACGTACGT' + body
                        handle.write(f'@{name}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
                rows.append((group, str(source)))
            manifest = root / 'samples.tsv'
            manifest.write_text('sample\tgroup\treplicate\tbarcode\tfastq\n' + ''.join(f'{g}\t{g}\tR1\tACGTACGT\t{s}\n' for g, s in rows))
            metadata = root / 'metadata.tsv'; metadata.write_text('sample\tgroup\treplicate\tdonor\nACUTE\tACUTE\tR1\tHD1\nCHRONIC\tCHRONIC\tR1\tHD1\n')
            reference = {'genome': None, 'genome_index': str(fasta), 'blacklist_bed': None, 'blacklist_genome': None}
            argv = ['run_full_prefix_rescue.py', '--manifest', str(manifest), '--metadata', str(metadata),
                    '--reference-config', str(root / 'ref.json'), '--outdir', str(root / 'output'),
                    '--barcode-orientation', 'forward', '--map-threads', '1', '--sort-threads', '0']
            with patch.object(sys, 'argv', argv), patch.object(run_full_prefix_rescue, 'load_reference_config', return_value=reference):
                run_full_prefix_rescue.main()
            report = json.loads((root / 'output/pilot.json').read_text())
            for group, _ in rows:
                sample = report['samples'][group]
                self.assertEqual(sample['mapping']['baseline']['deduplication']['deduplicated_molecules'], 1)
                combined = sample['mapping']['candidate_trim']['deduplication']
                self.assertEqual(combined['accepted_reads'], 2)
                self.assertEqual(combined['deduplicated_molecules'], 1)
                self.assertEqual(sample['endpoint_comparison']['same_five_prime_coordinate_and_strand'], 1)
                self.assertEqual(sample['newly_accepted_source_classes'], {'exact20_rescue': 1})
