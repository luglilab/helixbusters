"""Verify molecule conservation, interval boundaries and replicate support."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from helixbusters.regions import consensus_regions, load_sites, union_regions, window_regions, write_matrix


class TestRegions(unittest.TestCase):
    def test_incompatible_references_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts = root / 'counts.bed'
            counts.write_text('chr1\t0\t1\t1\n')
            headers = [root / 'a.json', root / 'b.json']
            headers[0].write_text(json.dumps([['chr1', 100]]))
            headers[1].write_text(json.dumps([['chr1', 101]]))
            with self.assertRaisesRegex(ValueError, 'incompatible'):
                load_sites([counts, counts], headers)

    def test_windows_conserve_molecules_and_clip_last_bin(self):
        sites = [{'chr1': [(0, 2), (9, 3), (10, 1), (24, 4)]}, {'chr1': [(10, 7)]}]
        regions = window_regions(sites, [('chr1', 25)], 10)
        self.assertEqual(regions, [('chr1', 0, 10), ('chr1', 10, 20), ('chr1', 20, 25)])
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'windows'
            write_matrix(prefix, regions, ['a', 'b'], sites)
            with Path(f'{prefix}.counts.tsv').open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual([int(row['a']) for row in rows], [5, 1, 4])
            self.assertEqual([int(row['b']) for row in rows], [0, 7, 0])

    def test_consensus_uses_distinct_replicates_and_exact_overlap(self):
        # Duplicate overlapping peaks from one replicate cannot supply support.
        peaks = {'a': [('chr1', 0, 10), ('chr1', 1, 9)],
                 'b': [('chr1', 8, 20)], 'c': [('chr1', 18, 30)]}
        self.assertEqual(consensus_regions(peaks, 2),
                         [('chr1', 8, 10, ('a', 'b')), ('chr1', 18, 20, ('b', 'c'))])
        self.assertEqual(consensus_regions(peaks, 3), [])
        with self.assertRaises(ValueError):
            consensus_regions(peaks, 4)

    def test_union_and_zero_count_regions(self):
        header = [('chr1', 100), ('chr2', 100)]
        self.assertEqual(union_regions([('chr1', 10, 20), ('chr1', 18, 30), ('chr2', 0, 10)], header),
                         [('chr1', 10, 30), ('chr2', 0, 10)])
        with self.assertRaises(ValueError):
            union_regions([('chr1', 99, 101)], header)

    def test_windows_cli_metadata_and_reference_mismatch(self):
        script = Path(__file__).parents[1] / 'scripts' / 'analyze_regions.py'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'header.json').write_text(json.dumps([['chr1', 25]]))
            (root / 'counts.bed').write_text('chr1\t9\t10\t2\nchr1\t10\t11\t1\n')
            (root / 'molecules.bed').write_text('')
            metadata = {'sample': 'a', 'group': 'treated', 'replicate': 'R1', 'donor': ''}
            (root / 'design.json').write_text(json.dumps({'design': 'unpaired', 'model_formula': '~ condition', 'sample_metadata': [metadata]}))
            command = [sys.executable, str(script), '--samples', 'a', '--counts', 'counts.bed',
                       '--headers', 'header.json', '--molecules', 'molecules.bed', '--design-file', 'design.json', '--windows', '10']
            result = subprocess.run(command, cwd=root, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((root / 'analysis.summary.json').read_text())
            self.assertEqual(report['windows']['10']['regions'], 2)
            self.assertIn('no differential model', report['status'])
            # Original inputs remain unchanged and outputs cannot be overwritten.
            self.assertEqual((root / 'counts.bed').read_text(), 'chr1\t9\t10\t2\nchr1\t10\t11\t1\n')
            again = subprocess.run(command, cwd=root, capture_output=True, text=True)
            self.assertNotEqual(again.returncode, 0)

    def test_peak_cli_nolambda(self):
        self.test_peak_cli_consensus_matrix_with_fake_macs3(nolambda=True)

    def test_peak_cli_consensus_matrix_with_fake_macs3(self, nolambda=False):
        """Exercise orchestration, not MACS3 statistical peak detection."""
        script = Path(__file__).parents[1] / 'scripts' / 'analyze_regions.py'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake = root / 'macs3'
            fake.write_text(f'#!{sys.executable}\n' + '''import sys
from pathlib import Path
args = sys.argv[1:]
if args == ['--version']:
    print('MACS3 test double')
else:
    assert args[0] == 'callpeak'
    assert args[args.index('--keep-dup') + 1] == 'all'
    assert args[args.index('--shift') + 1] == '-50'
    assert '--nomodel' in args and '-c' not in args
    sample = args[args.index('-n') + 1]
    folder = Path(args[args.index('--outdir') + 1])
    (folder / 'arguments.json').write_text(__import__('json').dumps(args))
    interval = (0, 20) if sample == 'a' else (10, 30)
    (folder / f'{sample}_peaks.narrowPeak').write_text(f'chr1\\t{interval[0]}\\t{interval[1]}\\tp\\t0\\t.\\t1\\t2\\t3\\t1\\n')
    (folder / f'{sample}_peaks.xls').write_text('test peak output')
    (folder / f'{sample}_summits.bed').write_text('chr1\\t10\\t11\\tsummit\\n')
''')
            fake.chmod(0o755)
            (root / 'header.json').write_text(json.dumps([['chr1', 100]]))
            metadata = []
            for sample, position in [('a', 10), ('b', 20)]:
                (root / f'{sample}.counts.bed').write_text(f'chr1\t{position}\t{position+1}\t1\n')
                (root / f'{sample}.molecules.bed').write_text(f'chr1\t{position}\t{position+1}\tACGT\t0\t+\n')
                metadata.append({'sample': sample, 'group': 'treated', 'replicate': sample, 'donor': ''})
            (root / 'design.json').write_text(json.dumps({'design': 'unpaired', 'model_formula': '~ condition', 'sample_metadata': metadata}))
            command = [sys.executable, str(script), '--samples', 'a', 'b', '--counts', 'a.counts.bed', 'b.counts.bed',
                       '--headers', 'header.json', 'header.json', '--molecules', 'a.molecules.bed', 'b.molecules.bed',
                       '--design-file', 'design.json', '--peaks', '--min-reps-consensus', '2']
            if nolambda:
                command.append('--nolambda')
            result = subprocess.run(command, cwd=root, env={**os.environ, 'PATH': f'{root}:{os.environ.get("PATH", "")}'}, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            for sample in ('a', 'b'):
                arguments = json.loads((root / f'SingleReplicate/{sample}/peaks/arguments.json').read_text())
                self.assertEqual('--nolambda' in arguments, nolambda)
                provenance = json.loads((root / f'SingleReplicate/{sample}/peaks/provenance.json').read_text())
                self.assertEqual('--nolambda' in provenance['command'], nolambda)
            summary = json.loads((root / 'analysis.summary.json').read_text())
            self.assertEqual(summary['peak_parameters']['nolambda'], nolambda)
            self.assertEqual((root / 'MergedReplicate/treated/peaks/treated.consensus.bed').read_text(),
                             'chr1\t10\t20\ttreated_1\t2\ta,b\n')
            with (root / 'peaks_consensus.counts.tsv').open() as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual((row['a'], row['b']), ('1', '0'))
            # Dedicated tasks must match the previous direct invocation and
            # aggregation must consume their outputs without calling MACS again.
            cli = script.with_name('call_bliss_peaks.py')
            peak_files, provenance_files = [], []
            for sample in ('a', 'b'):
                task = root / f'task_{sample}'
                task.mkdir()
                call = [sys.executable, str(cli), '--sample', sample,
                        '--counts', str(root / f'{sample}.counts.bed'), '--header', str(root / 'header.json'),
                        '--molecules', str(root / f'{sample}.molecules.bed'), '--effective-genome-size', '2913022398']
                if nolambda:
                    call.append('--nolambda')
                result = subprocess.run(call, cwd=task, env={**os.environ, 'PATH': f'{root}:{os.environ.get("PATH", "")}'}, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                peak_files.append(str(task / f'{sample}_peaks.narrowPeak'))
                provenance_files.append(str(task / f'{sample}.provenance.json'))
                self.assertTrue((task / f'{sample}.versions.json').is_file())
            aggregate = root / 'aggregate'
            aggregate.mkdir()
            external = [sys.executable, str(script), '--samples', 'a', 'b',
                        '--counts', str(root / 'a.counts.bed'), str(root / 'b.counts.bed'),
                        '--headers', str(root / 'header.json'), str(root / 'header.json'),
                        '--molecules', str(root / 'a.molecules.bed'), str(root / 'b.molecules.bed'),
                        '--design-file', str(root / 'design.json'), '--peaks', '--min-reps-consensus', '2',
                        '--peak-files', *peak_files, '--peak-provenance', *provenance_files]
            if nolambda:
                external.append('--nolambda')
            result = subprocess.run(external, cwd=aggregate, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse((aggregate / 'SingleReplicate').exists())
            self.assertEqual((aggregate / 'peaks_consensus.counts.tsv').read_text(), (root / 'peaks_consensus.counts.tsv').read_text())
            self.assertEqual(json.loads((aggregate / 'analysis.summary.json').read_text()), summary)
