"""Regression proof for GTF coordinates, strand priority and molecule allocation."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from helixbusters.annotation import GTFAnnotation
from helixbusters.gene_signal import write_gene_signal
from helixbusters.genomes import GenomeFilter

ROOT = Path(__file__).parents[1]

def fixture_gtf(path, chromosome='chr1'):
    rows = []
    for gene, strand, start, end, exons in [('A', '+', 100, 300, [(100, 120), (200, 220)]),
                                           ('B', '-', 400, 600, [(400, 430), (580, 600)])]:
        attrs = f'gene_id "{gene}"; gene_name "gene{gene}"; gene_type "protein_coding";'
        for feature, intervals in [('gene', [(start, end)]), ('transcript', [(start, end)]), ('exon', exons)]:
            for left, right in intervals:
                rows.append(f'{chromosome}\ttest\t{feature}\t{left+1}\t{right}\t.\t{strand}\t.\t{attrs}\n')
    path.write_text(''.join(rows))

class TestAnnotation(unittest.TestCase):
    def test_density_subtracts_union_blacklist_preserves_counts_and_changes_rank(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture_gtf(root / 'genes.gtf')
            # Overlapping rows must be subtracted once; chr aliases must match.
            (root / 'blacklist.bed').write_text('1\t130\t210\nchr1\t180\t290\nchr1\t80\t100\n')
            mask = GenomeFilter('hg38', root / 'blacklist.bed', 'hg38')
            annotation = GTFAnnotation(root / 'genes.gtf', [('chr1', 1000)], 20, 5)
            metadata = {s: {'sample':s, 'group':'condition', 'replicate':s, 'donor':''} for s in ['s1','s2']}
            sites = [{'chr1':[(105,2),(400,3)]}] * 2
            write_gene_signal(annotation, ['s1','s2'], sites, metadata, root/'masked', genome_filter=mask)
            write_gene_signal(annotation, ['s1','s2'], sites, metadata, root/'unmasked')
            def read(path):
                with path.open() as handle:
                    return list(csv.DictReader(handle, delimiter='\t'))
            raw = read(root/'masked/condition.gene_body.candidate_genes.tsv')
            density = read(root/'masked/condition.gene_body.density_candidate_genes.tsv')
            self.assertEqual([r['gene_id'] for r in raw], ['B','A'])
            self.assertEqual([r['gene_id'] for r in density], ['A','B'])
            gene = density[0]
            self.assertEqual(int(gene['uniquely_assignable_bp']), 195)
            self.assertEqual(int(gene['effective_assignable_bp']), 35)
            self.assertEqual(int(gene['excluded_assignable_bp']), 160)
            self.assertAlmostEqual(float(gene['median_CPM_per_effective_kb']), (2/5)*1e6*1000/35)
            for context in ('promoter','gene_body','combined'):
                masked = read(root/f'masked/genes.{context}.counts.tsv')
                unmasked = read(root/f'unmasked/genes.{context}.counts.tsv')
                self.assertEqual([(r['gene_id'],r['s1'],r['s2']) for r in masked],
                                 [(r['gene_id'],r['s1'],r['s2']) for r in unmasked])
            provenance = json.loads((root/'masked/gene_signal.provenance.json').read_text())
            self.assertEqual(provenance['blacklist_correction']['blacklist_sha256'], mask.sha256)
            with self.assertRaisesRegex(ValueError, 'blacklisted'):
                write_gene_signal(annotation,['s1'],[{'chr1':[(130,2)]}],{'s1':metadata['s1']},root/'invalid',genome_filter=mask)

    def test_half_open_promoters_both_strands_and_partition(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'genes.gtf'
            fixture_gtf(path, '1')
            index = GTFAnnotation(path, [('chr1', 1000)], 20, 5)
            for site, feature in [(79, 'intergenic'), (80, 'promoter'), (104, 'promoter'),
                                  (105, 'exon'), (120, 'intron'), (299, 'intron'), (300, 'intergenic'),
                                  (594, 'exon'), (595, 'promoter'), (619, 'promoter'), (620, 'intergenic')]:
                self.assertEqual(index.annotate('chr1', site, site+1)[0], feature)
            dominant, covered, genes = index.annotate('chr1', 70, 130)
            self.assertEqual(dominant, 'promoter')
            self.assertEqual(covered, {'promoter':25, 'exon':15, 'intron':10, 'intergenic':10})
            self.assertEqual(genes, ['A'])

    def test_alternative_transcript_tss_and_assembly_header(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'genes.gtf'
            fixture_gtf(path)
            with path.open('a') as handle:
                handle.write('chr1\ttest\ttranscript\t201\t300\t.\t+\t.\tgene_id "A";\n')
            index = GTFAnnotation(path, [('chr1', 1000)], 20, 5)
            self.assertEqual(index.annotate('chr1', 190, 191)[0], 'promoter')
            path.write_text('#!genome-build GRCh37\n' + path.read_text())
            with self.assertRaisesRegex(ValueError, 'genome build'):
                GTFAnnotation(path, [('chr1', 1000)], genome='hg38')

    def test_cli_conserves_molecules_and_exports_regions_and_multiqc(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture_gtf(root / 'genes.gtf')
            (root / 'header.json').write_text(json.dumps([['chr1', 248956422]]))
            (root / 'blacklist.bed').write_text('chr1\t900\t910\n')
            mask = GenomeFilter('hg38', root / 'blacklist.bed', 'hg38')
            (root / 'environment.json').write_text(json.dumps({'reference': {'genome':'hg38',
                'blacklist_bed':str(root/'blacklist.bed'), 'blacklist_genome':'hg38'},
                'genome_filter':mask.metadata()}))
            (root / 'a.bed').write_text(''.join(f'chr1\t{s}\t{s+1}\t{c}\n' for s,c in [(79,2),(80,3),(105,4),(120,5),(595,6),(620,7)]))
            (root / 'b.bed').write_text('chr1\t80\t81\t7\nchr1\t105\t106\t2\n')
            (root / 'design.json').write_text(json.dumps({'design':'unspecified', 'sample_metadata':[
                {'sample':'a','group':'A','replicate':'R1','donor':''},
                {'sample':'b','group':'A','replicate':'R2','donor':''}]}))
            (root / 'windows_10000.counts.tsv').write_text('region\tchrom\tstart\tend\ta\tb\nr1\tchr1\t70\t130\t14\t9\n')
            result = subprocess.run([sys.executable, str(ROOT / 'scripts/annotate_regions.py'),
                '--samples','a','b','--counts',str(root/'a.bed'),str(root/'b.bed'),
                '--headers',str(root/'header.json'),str(root/'header.json'), '--design-file',str(root/'design.json'),
                '--gtf',str(root/'genes.gtf'),'--gtf-genome','hg38','--genome','hg38',
                '--environment-file',str(root/'environment.json'),
                '--promoter-upstream','20','--promoter-downstream','5','--analysis-dir',str(root),
                '--outdir',str(root/'Annotation')],capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stderr)
            provenance = json.loads((root/'Annotation/GeneSignal/gene_signal.provenance.json').read_text())
            self.assertEqual(provenance['blacklist_correction']['blacklist_sha256'], mask.sha256)
            with (root/'Annotation/dsb_feature_distribution.samples.tsv').open() as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(sum(int(r['molecules']) for r in rows if r['sample']=='a'),27)
            self.assertEqual(sum(float(r['percentage']) for r in rows if r['sample']=='a'),100)
            with (root/'Annotation/dsb_feature_distribution.conditions.tsv').open() as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            promoter=next(r for r in rows if r['feature']=='promoter')
            self.assertAlmostEqual(float(promoter['mean_percentage']), (100*9/27+100*7/9)/2)
            self.assertTrue((root/'Annotation/windows_10000.annotations.tsv').is_file())
            self.assertGreater((root/'Annotation/dsb_feature_distribution.pdf').stat().st_size,100)
            with (root/'Annotation/GeneSignal/A.promoter.candidate_genes.tsv').open() as handle:
                genes = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual([row['gene_id'] for row in genes], ['A'])
            self.assertEqual(int(genes[0]['a_molecules']), 3)
            self.assertEqual(int(genes[0]['b_molecules']), 7)
            self.assertEqual(int(genes[0]['replicates_meeting_min_molecules']), 2)
            self.assertAlmostEqual(float(genes[0]['mean_CPM']), (3/27 + 7/9) * 1e6 / 2)
            with (root/'Annotation/GeneSignal/A.gene_body.candidate_genes.tsv').open() as handle:
                genes = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(int(genes[0]['a_molecules']), 9)
            self.assertEqual(int(genes[0]['b_molecules']), 2)
            if shutil.which('multiqc'):
                result=subprocess.run(['multiqc',str(root/'Annotation'),'--outdir',str(root/'MultiQC')],capture_output=True,text=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertIn('helixbusters_annotation_samples',result.stderr)
                self.assertIn('helixbusters_annotation_conditions',result.stderr)
                self.assertIn('helixbusters_gene_signal',result.stderr)

    def test_ambiguous_genes_are_not_double_counted_and_support_is_not_relaxed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture_gtf(root / 'genes.gtf')
            with (root / 'genes.gtf').open('a') as handle:
                for feature, start, end in [('gene',100,300), ('transcript',100,300), ('exon',100,120)]:
                    handle.write(f'chr1\ttest\t{feature}\t{start+1}\t{end}\t.\t+\t.\tgene_id "C";\n')
            annotation = GTFAnnotation(root / 'genes.gtf', [('chr1',1000)],20,5)
            metadata = {'s': {'sample':'s','group':'single','replicate':'R1','donor':''}}
            sites = [{'chr1':[(79,2),(80,3),(595,5)]}]
            result = write_gene_signal(annotation,['s'],sites,metadata,root/'GeneSignal',2,2)
            with (root/'GeneSignal/gene_assignment.samples.tsv').open() as handle:
                audit = next(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(int(audit['unique_gene_molecules']),5)
            self.assertEqual(int(audit['ambiguous_gene_molecules']),3)
            self.assertEqual(int(audit['intergenic_molecules']),2)
            self.assertEqual(result['single_promoter']['supported_candidates'],0)
            with (root/'GeneSignal/single.promoter.ranked_genes.tsv').open() as handle:
                ranked=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r['gene_id'] for r in ranked],['B'])
            self.assertEqual(int(ranked[0]['pooled_molecules']),5)

if __name__ == '__main__':
    unittest.main()
