"""Verify that BAM sensitivity preserves existing molecular units and strict ends."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import pysam

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_bam_peak_qvalues import representative_bam, filter_unique_primary
from helixbusters.deduplication import deduplicate_bam


class TestMolecularBam(unittest.TestCase):
    def test_clipping_branch_excludes_reported_multimappers_and_regroups_umis(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); bam=root/'source.bam'
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}) as output:
                specs=[('clip_AAAAAAAA',100,'5S45M',None,60),
                       ('strict_AAAAAAAA',100,'50M',None,60),
                       ('XA_CCCCCCCC',100,'50M',('XA','chr1,+201,50M,1;'),60),
                       ('NH_CCCCCCCC',100,'50M',('NH',2),60),
                       ('low_CCCCCCCC',100,'50M',None,0),
                       ('reverse_GGGGGGGG',200,'45M5S',None,60)]
                for name,start,cigar,tag,mapq in specs:
                    read=pysam.AlignedSegment();read.query_name=name;read.query_sequence='A'*50
                    read.reference_id=0;read.reference_start=start;read.cigarstring=cigar
                    read.flag=16 if name.startswith('reverse') else 0;read.mapping_quality=mapq
                    if tag: read.set_tag(*tag)
                    output.write(read)
            unique=root/'unique.bam'
            filtering=filter_unique_primary(bam,unique,20)
            self.assertEqual(filtering['reported_multimapping'],2)
            self.assertEqual(filtering['low_or_unknown_mapq'],1)
            self.assertEqual(filtering['retained_ambiguous_five_prime'],2)
            results={}
            for policy in ('strict','aligned'):
                results[policy]=deduplicate_bam(unique,root/f'{policy}.families',root/f'{policy}.counts',
                    method='directional',min_mapq=20,five_prime_policy=policy)
            self.assertEqual(results['strict']['deduplicated_molecules'],1)
            self.assertEqual(results['aligned']['accepted_reads'],3)
            self.assertEqual(results['aligned']['deduplicated_molecules'],2)
            target=root/'molecular.bam'
            self.assertEqual(representative_bam(unique,root/'aligned.families',target,20,'aligned'),2)
            with pysam.AlignmentFile(target,'rb') as handle:
                self.assertEqual([r.query_name for r in handle.fetch()],['clip_AAAAAAAA','reverse_GGGGGGGG'])

    def test_one_representative_per_family_both_strands(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); bam=root/'source.bam'
            header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}
            with pysam.AlignmentFile(bam,'wb',header=header) as output:
                specs=[('clipped_AAAAAAAA',100,False,'5S45M',False),
                       ('root_AAAAAAAA',100,False,'50M',True),
                       ('duplicate_AAAAAAAA',100,False,'50M',False),
                       ('member_AAAAAAAC',100,False,'50M',False),
                       ('different_CCCCCCCC',100,False,'50M',False),
                       ('reverse_GGGGGGGG',200,True,'50M',False),
                       ('reverse_dup_GGGGGGGG',200,True,'50M',True)]
                for name,start,reverse,cigar,duplicate in specs:
                    read=pysam.AlignedSegment(); read.query_name=name; read.query_sequence='A'*50
                    read.reference_id=0; read.reference_start=start; read.mapping_quality=60
                    read.flag=16 if reverse else 0; read.is_duplicate=duplicate; read.cigarstring=cigar
                    output.write(read)
            family=root/'families.tsv'
            family.write_text('chr1\t100\t101\t+\tAAAAAAAA\t3\nchr1\t100\t101\t+\tCCCCCCCC\t1\nchr1\t249\t250\t-\tGGGGGGGG\t2\n')
            target=root/'molecular.bam'
            self.assertEqual(representative_bam(bam,family,target),3)
            with pysam.AlignmentFile(target,'rb') as handle:
                reads=list(handle.fetch())
            self.assertEqual([r.query_name for r in reads],['root_AAAAAAAA','different_CCCCCCCC','reverse_GGGGGGGG'])
            self.assertFalse(any(r.is_duplicate for r in reads))
            self.assertTrue((root/'molecular.bam.bai').is_file())
            with self.assertRaises(FileExistsError):
                representative_bam(bam,family,target)

    def test_missing_family_fails_instead_of_changing_unit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); bam=root/'source.bam'
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':100}]}) as handle:
                pass
            family=root/'families.tsv'; family.write_text('chr1\t10\t11\t+\tAAAAAAAA\t1\n')
            with self.assertRaisesRegex(ValueError,'no matching representative'):
                representative_bam(bam,family,root/'molecular.bam')
