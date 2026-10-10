"""Prove exclusive raw-read fates and conservative experimental trimming."""
import gzip
import json
from pathlib import Path
import sys
import random
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from audit_original_bliss import audit_cohort
import audit_original_bliss
import analyze_prefix_sensitivity
from helixbusters.genomes import GenomeFilter
import pysam


class TestOriginalAudit(unittest.TestCase):
    @unittest.skipUnless(all(shutil.which(x) for x in ('bwa', 'samtools')), 'Mapping tools unavailable')
    def test_original_cohort_real_mapping_recovers_known_exact_prefix_inserts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rng = random.Random(47)
            genome = ''.join(rng.choices('ACGT', k=12000))
            fasta = root / 'reference.fa'; fasta.write_text('>chr1\n' + genome + '\n')
            subprocess.run(['bwa', 'index', str(fasta)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            source = root / 'original.fastq.gz'
            with gzip.open(source, 'wt') as handle:
                for i in range(20):
                    insert = genome[100 + i * 400:250 + i * 400]
                    if i >= 10:
                        insert = 'CCCTATAGTGAGTCGTAT' + insert
                    seq = 'AAAAAAAAACGTACGT' + insert
                    handle.write(f'@r{i}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
            manifest = root / 'manifest.tsv'
            manifest.write_text('sample\tgroup\treplicate\tbarcode\tfastq\ntoy\tA\tR1\tACGTACGT\t' + str(source) + '\n')
            argv = ['audit_original_bliss.py', '--manifest', str(manifest), '--outdir', str(root / 'audit'),
                    '--barcode-orientation', 'forward', '--reference-config', str(root / 'reference.json'),
                    '--max-reads', '100', '--map-threads', '1', '--sort-threads', '0']
            reference = {'genome': None, 'genome_index': str(fasta), 'blacklist_bed': None, 'blacklist_genome': None}
            with patch.object(sys, 'argv', argv), patch.object(audit_original_bliss, 'load_reference_config', return_value=reference):
                audit_original_bliss.main()
            report = json.loads((root / 'audit/audit.json').read_text())['samples']['toy']
            self.assertEqual(report['sampled_reads'], 20)
            self.assertEqual(report['mapping']['production_retained']['deduplication']['deduplicated_molecules'], 10)
            self.assertEqual(report['mapping']['excluded_exact_trim']['deduplication']['deduplicated_molecules'], 10)
    def test_exclusive_fates_preserve_umi_and_exclude_repeated_or_uncertain_prefixes(self):
        motif = 'CCCTATAGTGAGTCGTAT'; barcode = 'ACGTACGT'; genomic = 'G' * 60
        records = [('mismatch', 'AAAAAAAA' + 'T' * 8 + genomic),
                   ('invalid', 'NAAAAAAA' + barcode + genomic),
                   ('retained', 'AAAAAAAA' + barcode + genomic),
                   ('short', 'AAAAAAAA' + barcode + 'G'),
                   ('exact', 'CCCCCCCC' + barcode + motif + genomic),
                   ('repeat', 'AAAAAAAA' + barcode + motif + motif + genomic),
                   ('edited', 'AAAAAAAA' + barcode + 'A' + motif[1:] + genomic)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root / 'raw.gz'
            with gzip.open(source, 'wt') as handle:
                for name, seq in records:
                    handle.write(f'@{name}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
            result = audit_cohort(source, root, barcode, [barcode, 'TTTTTTTT'], motif)
            self.assertEqual(result['sampled_reads'], 7)
            self.assertEqual(result['exclusive_outcomes'], {'barcode_mismatch': 1, 'invalid_umi': 1,
                              'technical_prefix': 3, 'short_insert': 1, 'production_retained': 1})
            self.assertEqual(result['technical_prefix_classes']['exact_single_prefix_pilot'], 1)
            self.assertEqual(result['technical_prefix_classes']['exact_prefix_repeated_motif_first60'], 1)
            self.assertEqual(result['technical_prefix_classes']['edited_prefix_boundary_uncertain'], 1)
            self.assertEqual(result['barcode_diagnostics']['exact_other_manifest_barcode'], 1)
            with gzip.open(root / 'excluded_exact_trim.fastq.gz', 'rt') as handle:
                self.assertEqual(handle.read(), '@exact_CCCCCCCC\n' + genomic + '\n+\n' + 'I' * 60 + '\n')

    def test_low_quality_boundary_and_duplicate_identifiers_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root / 'raw.gz'; motif = 'CCCTATAGTGAGTCGTAT'
            seq = 'A' * 8 + 'ACGTACGT' + motif + 'G' * 60
            with gzip.open(source, 'wt') as handle:
                handle.write('@r\n' + seq + '\n+\n' + 'I' * 34 + '!' + 'I' * 59 + '\n')
            result = audit_cohort(source, root, 'ACGTACGT', ['ACGTACGT'], motif)
            self.assertEqual(result['technical_prefix_classes']['exact_prefix_low_quality_boundary'], 1)
            other = root / 'other'; other.mkdir()
            with gzip.open(source, 'at') as handle:
                handle.write('@r\n' + seq + '\n+\n' + 'I' * len(seq) + '\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                audit_cohort(source, other, 'ACGTACGT', ['ACGTACGT'], motif)

    def test_both_sensitivity_branches_are_annotated_with_explicit_paired_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root / 'source'; source.mkdir()
            samples = ['HD1_ACUTE', 'HD1_CHRONIC']
            (source / 'pilot.json').write_text(json.dumps({'parameters': {'all_reads': True, 'genome': 'hg38'},
                                                        'samples': {s: {} for s in samples}}))
            blacklist = root / 'blacklist.bed'; blacklist.write_text('chr1\t20000\t20100\n')
            sha = GenomeFilter('hg38', str(blacklist), 'hg38').sha256
            for sample in samples:
                for branch in ('baseline', 'candidate_trim'):
                    folder = source / sample / branch; folder.mkdir(parents=True)
                    (folder / f'{sample}.mapping.json').write_text(json.dumps({'genome_filter': {'blacklist_sha256': sha}}))
                    (folder / f'{sample}.counts.bed').write_text('chr1\t5000\t5001\t3\nchr1\t11000\t11001\t2\n')
                    (folder / f'{sample}.families.tsv').write_text('')
                    with pysam.AlignmentFile(folder / f'{sample}.q20.bam', 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 248956422}]}) as bam:
                        pass
            metadata = root / 'metadata.tsv'
            metadata.write_text('sample\tgroup\treplicate\tdonor\nHD1_ACUTE\tACUTE\tR1\tHD1\nHD1_CHRONIC\tCHRONIC\tR1\tHD1\n')
            gtf = root / 'genes.gtf'
            attrs = 'gene_id "G"; gene_name "G"; gene_type "protein_coding";'
            gtf.write_text(f'chr1\ttest\tgene\t1001\t9000\t.\t+\t.\t{attrs}\nchr1\ttest\ttranscript\t1001\t9000\t.\t+\t.\t{attrs}\nchr1\ttest\texon\t1001\t2000\t.\t+\t.\t{attrs}\n')
            reference = {'genome': 'hg38', 'blacklist_bed': str(blacklist), 'blacklist_genome': 'hg38'}
            argv = ['analyze_prefix_sensitivity.py', '--sensitivity-dir', str(source), '--metadata', str(metadata),
                    '--gtf', str(gtf), '--reference-config', str(root / 'references.json'), '--outdir', str(root / 'output')]
            with patch.object(sys, 'argv', argv), patch.object(analyze_prefix_sensitivity, 'load_reference_config', return_value=reference):
                analyze_prefix_sensitivity.main()
            for branch in ('baseline', 'candidate_trim'):
                folder = root / 'output' / branch
                self.assertTrue((folder / 'windows_100000.counts.tsv').is_file())
                self.assertTrue((folder / 'Annotation/dsb_feature_distribution.samples.tsv').is_file())
                self.assertIn('HD1', (folder / 'analysis.samples.tsv').read_text())
            self.assertFalse((root / 'output/baseline/Differential').exists())


if __name__ == '__main__':
    unittest.main()
