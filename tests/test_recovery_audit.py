"""Validate clipping orientation, uniform reservoir bounds and family conservation."""
from pathlib import Path
import sys
import gzip
import json
import random
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from audit_bliss_recovery import anchored_match, read_features, reservoir_bam, observed_complexity
from pilot_clipped_prefix import candidate_clip, write_cohort
import pilot_clipped_prefix


class TestRecoveryAudit(unittest.TestCase):
    @unittest.skipUnless(all(shutil.which(name) for name in ('bwa', 'samtools')), 'Mapping tools unavailable')
    def test_real_remapping_pilot_preserves_known_synthetic_boundaries(self):
        motif = 'CCCTATAGTGAGTCGTAT'
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = random.Random(47)
            sequence = ''.join(rng.choices('ACGT', k=12000))
            fasta = root / 'reference.fa'; fasta.write_text('>chr1\n' + sequence + '\n')
            subprocess.run(['bwa', 'index', str(fasta)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            folder = root / 'SingleReplicate/toy/mapping'; folder.mkdir(parents=True)
            with pysam.AlignmentFile(folder / 'toy.q20.bam', 'wb', header={'SQ': [{'SN': 'chr1', 'LN': len(sequence)}]}) as bam:
                for i in range(20):
                    start = 100 + i * 400; reverse = i % 2 == 1
                    genomic = sequence[start:start + 150]
                    original_genomic = genomic.translate(str.maketrans('ACGT', 'TGCA'))[::-1] if reverse else genomic
                    clip = 12 if i >= 10 else 0
                    original = motif[:clip] + original_genomic
                    read = pysam.AlignedSegment(bam.header)
                    read.query_name = f'r{i}_AAAAAAAA'; read.query_sequence = original.translate(str.maketrans('ACGT', 'TGCA'))[::-1] if reverse else original
                    read.query_qualities = [35] * len(original); read.reference_id = 0; read.reference_start = start
                    read.flag = 16 if reverse else 0; read.mapping_quality = 60
                    read.cigartuples = ([(0, 150), (4, clip)] if reverse else [(4, clip), (0, 150)]) if clip else [(0, 150)]
                    bam.write(read)
            argv = ['pilot_clipped_prefix.py', '--single-replicate-dir', str(root / 'SingleReplicate'),
                    '--samples', 'toy', '--reference-config', str(root / 'references.json'),
                    '--outdir', str(root / 'pilot'), '--all-reads', '--max-reads', '1', '--map-threads', '1', '--sort-threads', '0']
            reference = {'genome': None, 'genome_index': str(fasta), 'blacklist_bed': None, 'blacklist_genome': None}
            with patch.object(sys, 'argv', argv), patch.object(pilot_clipped_prefix, 'load_reference_config', return_value=reference):
                pilot_clipped_prefix.main()
            report = json.loads((root / 'pilot/pilot.json').read_text())['samples']['toy']
            self.assertEqual(report['cohort_reads'], 20)
            self.assertEqual(report['candidate_reads'], 10)
            self.assertEqual(report['candidate_trim']['deduplicated_molecules'], 20)
            self.assertEqual(report['new_endpoint_matches_source_aligned_boundary'], report['newly_strict_accepted'])
            self.assertEqual(report['unchanged_accepted_controls'], report['unchanged_controls_same_endpoint'])
    def test_pilot_changes_exact_technical_clip_only_and_preserves_rx(self):
        motif = 'CCCTATAGTGAGTCGTAT'
        for reverse in (False, True):
            header = pysam.AlignmentHeader.from_references(['chr1'], [1000])
            read = pysam.AlignedSegment(header)
            read.query_name = 'read_CCCCCCCC'; read.set_tag('RX', 'AAAAAAAA')
            sequence = motif[:12] + 'G' * 50
            read.query_sequence = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1] if reverse else sequence
            read.query_qualities = [35] * len(sequence)
            read.flag = 16 if reverse else 0
            read.reference_id = 0; read.reference_start = 100
            read.cigartuples = [(0, 50), (4, 12)] if reverse else [(4, 12), (0, 50)]
            self.assertEqual(candidate_clip(read, motif), 12)
            # A genomic motif prefix without 5-prime clipping must remain untouched.
            original = read.cigartuples; read.cigartuples = [(0, 62)]
            self.assertEqual(candidate_clip(read, motif), 0)
            read.cigartuples = original
            with tempfile.TemporaryDirectory() as tmp:
                changed = write_cohort([read], Path(tmp), motif)
                self.assertEqual(len(changed), 1)
                with gzip.open(Path(tmp) / 'candidate_trim.fastq.gz', 'rt') as handle:
                    lines = handle.read().splitlines()
                self.assertEqual(lines[0], '@read_CCCCCCCC_AAAAAAAA')
                self.assertEqual(lines[1], 'G' * 50)
                self.assertEqual(len(lines[3]), 50)
            quality = [5] * 12 + [35] * 50
            read.query_qualities = quality[::-1] if reverse else quality
            self.assertEqual(candidate_clip(read, motif), 0)
    def test_forward_reverse_clips_and_qualities_are_in_read_orientation(self):
        motif = 'CCCTATAGTGAGTCGTAT'
        sequence = motif[:8] + 'G' * 20
        qualities = [5] * 8 + [35] * 20
        for reverse in (False, True):
            read = pysam.AlignedSegment()
            read.query_name = 'test_AAAAAAAA'
            read.flag = 16 if reverse else 0
            read.reference_id = 0; read.reference_start = 100; read.mapping_quality = 60
            read.query_sequence = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1] if reverse else sequence
            read.query_qualities = qualities[::-1] if reverse else qualities
            read.cigartuples = [(0, 20), (4, 8)] if reverse else [(4, 8), (0, 20)]
            info = read_features(read, motif)
            self.assertEqual(info['reason'], 'five_prime_S')
            self.assertEqual(info['clip_prefix'], motif[:8])
            self.assertEqual(info['clip_mean_quality'], 5)
            self.assertEqual(info['clip_motif_exact_bases'], 8)
            self.assertEqual(info['first10_mean_quality'], 11)
        read.cigartuples = [(0, 28)]
        self.assertEqual(read_features(read, motif)['reason'], 'accepted')

    def test_exact_motif_segments_do_not_skip_leading_bases(self):
        motif = 'CCCTATAGTGAGTCGTAT'
        self.assertEqual(anchored_match('TAGTGAGTCGTATGG', motif), 13)
        self.assertLess(anchored_match('N' + motif, motif), 8)

    def test_observed_family_curve_conserves_reads_and_molecules(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'families.tsv'
            p.write_text('chr1\t1\t2\t+\tAAAAAAAA\t1\nchr1\t2\t3\t+\tCCCCCCCC\t3\n')
            sizes, curve = observed_complexity(p, 4, 2)
            self.assertEqual(curve[0]['expected_observed_families'], 0)
            self.assertEqual(curve[-1]['expected_observed_families'], 2)
            self.assertAlmostEqual(curve[10]['expected_observed_families'], .5 + .875)
            with self.assertRaises(ValueError):
                observed_complexity(p, 5, 2)

    def test_reservoir_is_deterministic_and_scans_all_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'reads.bam'
            with pysam.AlignmentFile(p, 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 1000}]}) as handle:
                for i in range(100):
                    read = pysam.AlignedSegment()
                    read.query_name = f'r{i}'; read.query_sequence = 'A' * 20
                    read.reference_id = 0; read.reference_start = i; read.cigartuples = [(0, 20)]
                    handle.write(read)
            a, seen = reservoir_bam(p, 10, 1729); b, _ = reservoir_bam(p, 10, 1729)
            self.assertEqual(seen, 100)
            self.assertEqual([r.query_name for r in a], [r.query_name for r in b])
            self.assertEqual(len(a), 10)
            self.assertTrue(any(r.reference_start > 50 for r in a))


if __name__ == '__main__':
    unittest.main()
