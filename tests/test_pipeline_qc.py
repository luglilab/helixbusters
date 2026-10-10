"""Regression checks for layout gates, read-flow denominators and titration labels."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from validate_read_layout import validate
from qc_read_flow import read_flow


class TestPipelineQC(unittest.TestCase):
    def test_layout_cli_stops_wrong_orientation_before_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);reads=root/'reads.fastq'
            sequence='AACCTTGGCGTGTGAG'+'GATTACA'*6
            reads.write_text(''.join(f'@r{i}\n{sequence}\n+\n'+ 'I'*len(sequence)+'\n' for i in range(20)))
            for orientation,passed in [('forward',True),('reverse_complement',False)]:
                output=root/orientation;output.mkdir()
                result=subprocess.run([sys.executable,str(ROOT/'scripts/validate_read_layout.py'),
                    '--reads',str(reads),'--sample','test','--barcode','CGTGTGAG','--orientation',orientation],
                    cwd=output,capture_output=True,text=True)
                self.assertEqual(result.returncode==0,passed,result.stderr)
                self.assertEqual(json.loads((output/'test.layout.json').read_text())['check']['passed'],passed)

    def test_layout_orientation_and_minimum_signal(self):
        result={'reads_examined':100,'barcode_search':[{'barcode':'AACCGGTA','orientation':orientation,'positions':[{'offset_0based':8,'reads':n}]} for orientation,n in [('forward',80),('reverse_complement',2)]]}
        self.assertTrue(validate(result,'AACCGGTA','forward',8,.1)['passed'])
        self.assertFalse(validate(result,'AACCGGTA','reverse_complement',8,.1)['passed'])
        self.assertFalse(validate(result,'AACCGGTA','forward',9,.1)['passed'])
        self.assertFalse(validate(result,'AACCGGTA','forward',8,.9)['passed'])

    def test_flow_uses_raw_denominators_without_counting_alignment_records(self):
        prep={'sample':'A','counts':{'input_reads':100,'output_reads':80,'barcode_mismatch_reads':20}}
        sample={'sample':'A','group':'G','metrics':{'primary_records':80,'alignment_records':120,'retained_records':20,'deduplicated_molecules':8,'ambiguous_five_prime_reads':10,'duplicate_reads':2,'excluded_low_mapq':50}}
        flow=read_flow(sample,prep)
        self.assertEqual(flow['usable_molecules_pct_of_raw'],8)
        self.assertEqual(flow['excluded_low_mapq_pct_of_raw'],50)
        self.assertEqual(flow['ambiguous_5prime_pct_of_retained'],50)
        self.assertEqual(flow['largest_recorded_loss'],'excluded_low_mapq')
        sample['metrics']['primary_records']=79
        with self.assertRaises(ValueError):read_flow(sample,prep)

    def test_input_titration_requires_one_explicit_donor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for donors,passed in [('D1,D1',True),('D1,D2',False)]:
                a,b=donors.split(',');manifest=root/'manifest.tsv'
                manifest.write_text(f'sample\tgroup\treplicate\tdonor\nA\t1M\tR1\t{a}\nB\t50k\tR1\t{b}\n')
                output=root/donors;output.mkdir()
                result=subprocess.run([sys.executable,str(ROOT/'scripts/validate_experimental_design.py'),'--manifest',str(manifest),'--purpose','input_titration'],cwd=output,capture_output=True,text=True)
                self.assertEqual(result.returncode==0,passed,result.stderr)
                if passed:
                    report=json.loads((output/'design.summary.json').read_text())
                    self.assertEqual(report['analysis_purpose'],'input_titration')
                    self.assertIn('not independent biological replicates',report['notes'])
