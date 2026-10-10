"""Exact fixed-bin membership, biological replication and depth sensitivity."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from window_upset import run


class TestWindowUpSet(unittest.TestCase):
    def test_input_titration_uses_technical_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.setup_input(root)
            metadata=pd.read_csv(root/'analysis.samples.tsv',sep='\t')
            metadata['group']=['1M','500k','100k','25k'];metadata['donor']='D1';metadata['replicate']='R1'
            metadata.to_csv(root/'analysis.samples.tsv',sep='\t',index=False)
            source=json.loads((root/'analysis.summary.json').read_text());source['analysis_purpose']='input_titration'
            (root/'analysis.summary.json').write_text(json.dumps(source))
            result=run(root,root/'overlap',[10],minimum_replicates=1,iterations=3)
            self.assertIn('not biological replicates',result['interpretation'])
            self.assertIn('technical', (root/'overlap/README.md').read_text().lower())

    def setup_input(self,root):
        names=['a1','a2','b1','b2']
        pd.DataFrame({'sample':names,'group':['A','A','B','B'],'replicate':['R1','R2','R1','R2'],
                      'donor':['D1','D2','D1','D2']}).to_csv(root/'analysis.samples.tsv',sep='\t',index=False)
        pd.DataFrame({'region':['r1','r2','r3'],'chrom':['chr1']*3,'start':[0,10,20],'end':[10,20,25],
                      'a1':[5,0,5],'a2':[5,0,5],'b1':[0,5,5],'b2':[0,5,5]}).to_csv(root/'windows_10.counts.tsv',sep='\t',index=False)
        (root/'analysis.summary.json').write_text(json.dumps({'samples':{n:{'molecules':10} for n in names}}))

    def test_fixed_bins_are_not_merged_and_shared_membership_is_exact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.setup_input(root)
            original=(root/'windows_10.counts.tsv').read_bytes()
            result=run(root,root/'overlap',[10],iterations=3)
            records=result['windows']['10']['full_depth']
            self.assertEqual({r['conditions']:r['windows'] for r in records},{'A':1,'B':1,'A,B':1})
            self.assertEqual({r['conditions']:r['covered_bp'] for r in records},{'A':10,'B':10,'A,B':5})
            self.assertEqual(records,result['windows']['10']['equal_depth'])
            membership=pd.read_csv(root/'overlap/windows_10/window_membership.tsv',sep='\t')
            self.assertEqual(membership.A_supported.tolist(),[True,False,True])
            self.assertTrue((root/'overlap/windows_10/full_depth/conditions_upset.pdf').is_file())
            self.assertEqual(original,(root/'windows_10.counts.tsv').read_bytes())
            with self.assertRaisesRegex(ValueError,'new output'):
                run(root,root/'overlap',[10],iterations=3)

    def test_duplicate_biological_units_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.setup_input(root)
            path=root/'analysis.samples.tsv';data=pd.read_csv(path,sep='\t');data.loc[1,'donor']='D1';data.to_csv(path,sep='\t',index=False)
            with self.assertRaisesRegex(ValueError,'Duplicate donor'):
                run(root,root/'overlap',[10],iterations=3)

    def test_empty_support_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.setup_input(root)
            result=run(root,root/'overlap',[10],minimum_molecules=6,iterations=3)
            self.assertTrue(all(r['windows']==0 for r in result['windows']['10']['full_depth']))
