"""Check reproducible sampling and coordinate comparison, including real toy mapping."""
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

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from pilot_aligner_comparison import compare_reads, reservoir
import pilot_aligner_comparison


class TestAlignerPilot(unittest.TestCase):
    def test_reservoir_is_deterministic_and_preserves_records(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input.fastq'
            source.write_text(''.join(f'@read{i}_AAAAAAAA\nACGT\n+\nIIII\n' for i in range(100)))
            reports = [reservoir(source, root / f'{i}.gz', 10, 1729) for i in range(2)]
            self.assertEqual(reports[0], reports[1])
            self.assertEqual((root / '0.gz').read_bytes(), (root / '1.gz').read_bytes())
            self.assertEqual(reports[0]['source_reads'], 100)
            with gzip.open(root / '0.gz', 'rt') as handle:
                lines = handle.read().splitlines()
            self.assertEqual(len(lines), 40)
            self.assertTrue(all(lines[i + 1:i + 4] == ['ACGT', '+', 'IIII'] for i in range(0, 40, 4)))

    def test_concordance_distinguishes_coordinate_and_strand_changes(self):
        a = {'same': {'key': ('chr1', 10, '+')}, 'changed': {'key': ('chr1', 20, '+')}, 'a_only': {'key': ('chr1', 1, '+')}}
        b = {'same': {'key': ('chr1', 10, '+')}, 'changed': {'key': ('chr1', 20, '-')}, 'b_only': {'key': ('chr2', 1, '+')}}
        result = compare_reads(a, b)
        self.assertEqual(result['same_five_prime_coordinate_and_strand'], 1)
        self.assertEqual(result['discordant_coordinate_or_strand'], 1)
        self.assertEqual(result['accepted_only_left'], 1)
        self.assertEqual(result['accepted_only_right'], 1)

    @unittest.skipUnless(all(shutil.which(name) for name in ('bwa', 'bowtie2', 'bowtie2-build', 'samtools')), 'Mapping tools unavailable')
    def test_three_branches_map_identical_toy_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = random.Random(17)
            genome = ''.join(rng.choices('ACGT', k=6000))
            fasta = root / 'genome.fa'
            fasta.write_text('>chr1\n' + genome + '\n')
            subprocess.run(['bwa', 'index', str(fasta)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            bt = root / 'bowtie'
            subprocess.run(['bowtie2-build', str(fasta), str(bt)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            blacklist = root / 'blacklist.bed'
            blacklist.write_text('chr1\t0\t1\n')
            reference = root / 'references.json'
            reference.write_text(json.dumps({'schema_version': 1, 'references': {'hg38': {
                'bwa_index': str(fasta), 'bowtie2_index': str(bt), 'blacklist': {'genome': 'hg38', 'path': str(blacklist)}}}}))
            prepared = root / 'prepared'
            prepared.mkdir()
            with gzip.open(prepared / 'toy.umi.fastq.gz', 'wt') as handle:
                for i in range(30):
                    start = 100 + i * 120
                    seq = genome[start:start + 150]
                    handle.write(f'@read{i}_AAAAAAAA\n{seq}\n+\n' + 'I' * len(seq) + '\n')
            argv = [str(ROOT / 'scripts/pilot_aligner_comparison.py'),
                '--prepared-dir', str(prepared), '--samples', 'toy', '--reference-config', str(reference),
                '--outdir', str(root / 'results'), '--max-reads', '20', '--map-threads', '1', '--sort-threads', '0']
            # An unlabeled toy genome avoids claiming this 6 kb reference is hg38.
            # Production reference loading and hg38 length validation stay enabled.
            def toy_reference(config, genome_name, aligner):
                return {'genome': None, 'genome_index': str(fasta if aligner == 'bwa' else bt),
                        'blacklist_bed': None, 'blacklist_genome': None}
            with patch.object(sys, 'argv', argv), patch.object(
                    pilot_aligner_comparison, 'load_reference_config', side_effect=toy_reference):
                pilot_aligner_comparison.main()
            report = json.loads((root / 'results/aligner_comparison.json').read_text())['samples']['toy']
            for branch in report['branches'].values():
                self.assertEqual(branch['metrics']['input_reads'], 20)
                self.assertEqual(branch['metrics']['molecules'], 20)
            for pair in report['pairwise'].values():
                self.assertEqual(pair['same_five_prime_coordinate_and_strand'], 20)


if __name__ == '__main__':
    unittest.main()
