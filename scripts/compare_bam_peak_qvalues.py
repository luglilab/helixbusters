#!/usr/bin/env python3
"""Test strand-directed MACS3 BAM extension using one read per existing UMI family."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.deduplication import five_prime_position, _read_umi, deduplicate_bam


def representative_bam(source, families, destination, min_mapq=20, five_prime_policy='strict'):
    """Select one valid read carrying each existing directional-family root UMI."""
    destination=Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    pending=set()
    with Path(families).open() as handle:
        for line in handle:
            chrom,start,end,strand,umi,count=line.rstrip().split('\t')
            key=(chrom,int(start),strand,umi)
            if key in pending or int(end)!=int(start)+1:
                raise ValueError('Invalid or duplicate UMI family row')
            pending.add(key)
    expected=len(pending)
    with pysam.AlignmentFile(source,'rb') as original, pysam.AlignmentFile(destination,'wb',template=original) as output:
        for read in original.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if read.is_paired:
                raise ValueError('This pilot requires single-end BLISS libraries')
            if read.mapping_quality==255 or read.mapping_quality<min_mapq:
                continue
            position=five_prime_position(read,five_prime_policy)
            umi,reason=_read_umi(read,8)
            if position is None or reason:
                continue
            key=(read.reference_name,position,'-' if read.is_reverse else '+',umi)
            if key in pending:
                # Coordinate-only duplicate flags must not hide distinct UMI molecules.
                read.is_duplicate=False
                output.write(read); pending.remove(key)
    if pending:
        raise ValueError(f'{len(pending)} UMI families have no matching representative read')
    pysam.index(str(destination))
    return expected


def filter_unique_primary(source, destination, min_mapq):
    """Retain filtered primary alignments without reported alternative placements.

    XA and NH detect reported alternatives; their absence cannot prove perfect
    genomic uniqueness. MAPQ remains the primary alignment-confidence filter.
    """
    destination=Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    counts=Counter()
    with pysam.AlignmentFile(source,'rb') as original, pysam.AlignmentFile(destination,'wb',template=original) as output:
        for read in original.fetch(until_eof=True):
            counts['input_reads']+=1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                reason='nonprimary_or_unusable'
            elif read.mapping_quality==255 or read.mapping_quality<min_mapq:
                reason='low_or_unknown_mapq'
            elif ((read.has_tag('NH') and read.get_tag('NH')>1) or
                  (read.has_tag('XA') and bool(read.get_tag('XA')))):
                reason='reported_multimapping'
            else:
                if read.is_paired:
                    raise ValueError('Single-end BLISS input required')
                output.write(read); counts['retained_reads']+=1
                if five_prime_position(read,'strict') is None and five_prime_position(read,'aligned') is not None:
                    counts['retained_ambiguous_five_prime']+=1
                continue
            counts[reason]+=1
    return dict(counts)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--accept-five-prime-clipping',action='store_true',
                   help='MACS3-only branch: exclude reported multimappers and regroup UMIs at aligned boundaries')
    args=p.parse_args()
    if shutil.which('macs3') is None:
        p.error('Activate an environment containing MACS3')
    source=args.source.resolve(); out=args.outdir.resolve()
    if out.exists():
        p.error('Choose a new output directory')
    design=source/'MultiQC/design.summary.json'
    metadata=json.loads(design.read_text())
    summary=json.loads((source/'Analysis/analysis.summary.json').read_text())
    samples=[r['sample'] for r in metadata['sample_metadata']]
    paths=[]
    for sample in samples:
        base=source/'SingleReplicate'/sample
        bams=list((base/'mapping').glob('*.q*.bam'))
        family=base/'deduplication'/f'{sample}.families.tsv'
        counts=base/'deduplication'/f'{sample}.counts.bed'
        molecules=base/'deduplication'/f'{sample}.molecules.bed'
        if len(bams)!=1 or not all(x.is_file() for x in (family,counts,molecules)):
            p.error('Missing or ambiguous molecular inputs for '+sample)
        qc=json.loads((base/'qc'/f'{sample}.summary.json').read_text())
        if qc['deduplication']['parameters']['five_prime_policy']!='strict':
            p.error('Source must use strict 5-prime acceptance')
        if args.accept_five_prime_clipping:
            genome=qc['mapping']['genome_filter']
            if (qc['mapping']['parameters']['min_mapq']<20 or
                    not genome.get('exclude_mitochondrial') or not genome.get('blacklist_sha256') or
                    not genome.get('canonical_nuclear_chromosomes')):
                p.error('Clipping branch requires the canonical, mitochondrial/blacklist filtered MAPQ>=20 source')
        paths.append((bams[0],family,counts,molecules,qc['mapping']['parameters']['min_mapq']))
    out.mkdir(parents=True)
    branch_qc={}
    policy='aligned' if args.accept_five_prime_clipping else 'strict'
    if args.accept_five_prime_clipping:
        new_paths=[]
        for sample,(bam,_,_,_,mapq) in zip(samples,paths,strict=True):
            stage=out/'molecular_preparation'/sample; stage.mkdir(parents=True)
            unique=stage/f'{sample}.unique_primary.bam'
            filtering=filter_unique_primary(bam,unique,mapq)
            controls={}
            for branch in ('strict','aligned'):
                dest=stage/branch; dest.mkdir()
                family=dest/f'{sample}.families.tsv'; counts=dest/f'{sample}.counts.bed'
                molecules=dest/f'{sample}.molecules.bed'
                controls[branch]=deduplicate_bam(unique,family,counts,method='directional',
                    umi_length=8,min_mapq=mapq,five_prime_policy=branch,output_molecules=molecules,
                    output_qc=dest/f'{sample}.deduplication.json')
                if branch=='aligned':
                    new_paths.append((unique,family,counts,molecules,mapq))
            branch_qc[sample]={'original_strict_molecules':summary['samples'][sample]['molecules'],
                              'filtering':filtering,'unique_strict':controls['strict'],
                              'unique_aligned':controls['aligned']}
            summary['samples'][sample]['molecules']=controls['aligned']['deduplicated_molecules']
        paths=new_paths
        (out/'clipping_branch.QC.json').write_text(json.dumps(branch_qc,indent=2)+'\n')
    bam_paths=[]; headers=[]
    for sample,(bam,family,_,_,mapq) in zip(samples,paths,strict=True):
        target=out/f'{sample}.molecular.bam'
        print(f'{sample}: selecting one read per existing UMI family',flush=True)
        n=representative_bam(bam,family,target,mapq,policy)
        if n!=summary['samples'][sample]['molecules']:
            raise ValueError('Molecular BAM count disagrees with source summary')
        with pysam.AlignmentFile(target,'rb') as handle:
            header=list(zip(handle.references,handle.lengths))
        h=out/f'{sample}.header.json'; h.write_text(json.dumps(header)+'\n')
        headers.append(str(h)); bam_paths.append(target)
    version=subprocess.check_output(['macs3','--version'],text=True).strip()
    all_commands=[]
    for q in (0.05,0.10):
        folder=out/f'q_{q:.2f}'; folder.mkdir()
        peak_files=[]; provenances=[]
        for sample,bam in zip(samples,bam_paths,strict=True):
            task=folder/'SingleReplicate'/sample/'peaks'; task.mkdir(parents=True)
            command=['macs3','callpeak','-t',str(bam),'-f','BAM','-g','hs','-n',sample,
                     '--outdir',str(task),'--nomodel','--extsize','80','--keep-dup','all','-q',str(q)]
            with (task/'macs3.log').open('x') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            prov={'sample':sample,'molecules':summary['samples'][sample]['molecules'],
                  'width':80,'qvalue':q,'effective_genome_size':2913022398,'nolambda':False,
                  'shift':0,'format':'BAM','version':version,'command':command,
                  'five_prime_policy':policy,
                  'unit':'one read per directional UMI family; '+policy+' boundary'}
            path=task/f'{sample}.provenance.json'; path.write_text(json.dumps(prov,indent=2)+'\n')
            peak_files.append(str(task/f'{sample}_peaks.narrowPeak')); provenances.append(str(path))
            all_commands.append(command)
        analysis=folder/'Analysis'; analysis.mkdir()
        command=[sys.executable,str(Path(__file__).with_name('analyze_regions.py')),
                 '--samples',*samples,'--counts',*[str(x[2]) for x in paths],'--headers',*headers,
                 '--molecules',*[str(x[3]) for x in paths],'--design-file',str(design),'--peaks',
                 '--min-reps-consensus','2','--peak-width','80','--peak-qvalue',str(q),
                 '--effective-genome-size','2913022398','--peak-files',*peak_files,
                 '--peak-provenance',*provenances]
        subprocess.run(command,cwd=analysis,check=True)
    (out/'provenance.json').write_text(json.dumps({'source':str(source),'commands':all_commands,
        'five_prime_policy':policy,
        'interpretation':'exploratory sensitivity: BAM, extsize80, shift0, local background; not differential testing',
        'scope':'MACS3 only; source gene/window counts and models are unchanged'},indent=2)+'\n')


if __name__=='__main__':
    main()
