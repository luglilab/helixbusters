"""Validate T7 orientation, physical trimming length and strand-aware MD diagnostics."""
import gzip
from pathlib import Path
import sys
import random
import json
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pysam

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from inspect_full_prefix_boundary import end_mismatches, prepare_comparison, T7_FORWARD, T7_REVERSE
import inspect_full_prefix_boundary


class TestFullPrefixBoundary(unittest.TestCase):
    @unittest.skipUnless(all(shutil.which(x) for x in ('bwa', 'samtools')), 'Mapping tools unavailable')
    def test_real_twenty_base_remapping_recovers_known_boundaries_on_both_strands(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rng = random.Random(47)
            genome = ''.join(rng.choices('ACGT', k=12000))
            fasta = root / 'reference.fa'; fasta.write_text('>chr1\n' + genome + '\n')
            subprocess.run(['bwa', 'index', str(fasta)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            source = root / 'audit'; folder = source / 'toy'; folder.mkdir(parents=True)
            (source / 'audit.json').write_text(json.dumps({'parameters': {'motif': T7_REVERSE[:18], 'genome': 'hg38'}, 'samples': {'toy': {}}}))
            with gzip.open(folder / 'excluded_exact_untrimmed.fastq.gz', 'wt') as handle:
                for i in range(10):
                    genomic = genome[100 + i * 400:250 + i * 400]
                    if i % 2:
                        genomic = genomic.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                    seq = T7_REVERSE + genomic
                    handle.write(f'@r{i}_AAAAAAAA\n{seq}\n+\n' + 'I' * len(seq) + '\n')
            reference = {'genome': None, 'genome_index': str(fasta), 'blacklist_bed': None, 'blacklist_genome': None}
            argv = ['inspect_full_prefix_boundary.py', '--audit-dir', str(source), '--outdir', str(root / 'output'),
                    '--reference-config', str(root / 'references.json'), '--map-threads', '1', '--sort-threads', '0']
            with patch.object(sys, 'argv', argv), patch.object(inspect_full_prefix_boundary, 'load_reference_config', return_value=reference):
                inspect_full_prefix_boundary.main()
            result = json.loads((root / 'output/boundary.json').read_text())['samples']['toy']
            self.assertEqual(result['cohort']['comparison_reads'], 10)
            self.assertEqual(result['mapping']['trim20']['deduplication']['deduplicated_molecules'], 10)
            accepted = inspect_full_prefix_boundary.accepted_reads(root / 'output/toy/trim20/toy.q20.bam')
            for i in range(10):
                expected = 249 + i * 400 if i % 2 else 100 + i * 400
                self.assertEqual(accepted[f'r{i}_AAAAAAAA']['key'], ('chr1', expected, '-' if i % 2 else '+'))
    def test_known_T7_reverse_complement_is_twenty_bases(self):
        self.assertEqual(len(T7_REVERSE), 20)
        self.assertEqual(T7_REVERSE, T7_FORWARD.translate(str.maketrans('ACGT', 'TGCA'))[::-1])

    def test_MD_mismatch_uses_biological_five_prime_for_both_strands(self):
        for reverse in (False, True):
            a = pysam.AlignedSegment()
            a.query_sequence = ('T' if reverse else 'A') * 10
            a.flag = 16 if reverse else 0; a.reference_start = 100; a.cigartuples = [(0, 10)]
            a.set_tag('MD', '9G0' if reverse else '0C9')
            self.assertEqual(end_mismatches(a), [True] + [False] * 9)
            a.set_tag('MD', None)
            self.assertEqual(end_mismatches(a), [None] * 10)

    def test_comparison_preserves_names_umi_and_excludes_both_repeat_orientations(self):
        records = [T7_REVERSE + 'G' * 60,
                   T7_REVERSE + T7_FORWARD + 'G' * 60,
                   T7_REVERSE + T7_REVERSE + 'G' * 60,
                   T7_REVERSE[:18] + 'CC' + 'G' * 60,
                   T7_REVERSE + 'G' * 20]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root / 'source.gz'
            with gzip.open(source, 'wt') as handle:
                for i, seq in enumerate(records):
                    handle.write(f'@r{i}_ACGTACGT\n{seq}\n+\n' + 'I' * len(seq) + '\n')
            counts = prepare_comparison(source, root)
            self.assertEqual(counts, {'input_reads': 5, 'comparison_reads': 1,
                                     'residual_T7_in_either_orientation': 2,
                                     'not_exact_full20_prefix': 1, 'short_tail': 1})
            for branch, expected in [('trim18', 'TA' + 'G' * 60), ('trim20', 'G' * 60)]:
                with gzip.open(root / f'{branch}.fastq.gz', 'rt') as handle:
                    self.assertEqual(handle.read(), '@r0_ACGTACGT\n' + expected + '\n+\n' + 'I' * len(expected) + '\n')


if __name__ == '__main__':
    unittest.main()
