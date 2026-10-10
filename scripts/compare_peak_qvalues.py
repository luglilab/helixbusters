#!/usr/bin/env python3
"""Recall peaks at prespecified thresholds from existing molecular libraries."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import pysam


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--qvalues', nargs='+', type=float, default=[0.05,0.10])
    args=p.parse_args()
    if not args.qvalues or len(set(args.qvalues))!=len(args.qvalues) or any(not 0<q<1 for q in args.qvalues):
        p.error('Specify distinct qvalues between 0 and 1')
    if shutil.which('macs3') is None:
        p.error('MACS3 must be available in the active environment')
    source=args.source.resolve(); out=args.outdir.resolve()
    if out.exists():
        p.error('Output already exists; choose a new directory')
    design=source/'MultiQC/design.summary.json'
    metadata=json.loads(design.read_text())
    settings=json.loads((source/'Analysis/analysis.summary.json').read_text())['peak_parameters']
    samples=[row['sample'] for row in metadata['sample_metadata']]
    paths=[]
    for sample in samples:
        folder=source/'SingleReplicate'/sample/'deduplication'
        counts=folder/f'{sample}.counts.bed'; molecules=folder/f'{sample}.molecules.bed'
        bams=list((source/'SingleReplicate'/sample/'mapping').glob('*.q*.bam'))
        if not counts.is_file() or not molecules.is_file() or len(bams)!=1:
            p.error('Missing molecular inputs or ambiguous filtered BAM for '+sample)
        paths.append((counts,molecules,bams[0]))
    out.mkdir(parents=True)
    headers=[]
    for sample,(_,_,bam) in zip(samples,paths,strict=True):
        with pysam.AlignmentFile(bam,'rb') as handle:
            header=list(zip(handle.references,handle.lengths))
        path=out/f'{sample}.header.json'; path.write_text(json.dumps(header)+'\n'); headers.append(str(path))
    script=Path(__file__).with_name('analyze_regions.py')
    commands=[]
    for q in args.qvalues:
        folder=out/f'q_{q:g}'; folder.mkdir()
        command=[sys.executable,str(script),'--samples',*samples,'--counts',*[str(x[0]) for x in paths],
                 '--headers',*headers,'--molecules',*[str(x[1]) for x in paths],
                 '--design-file',str(design),'--peaks','--min-reps-consensus','2',
                 '--peak-width',str(settings['width']),'--peak-qvalue',str(q),
                 '--effective-genome-size',str(settings['effective_genome_size'])]
        if settings['nolambda']:
            command.append('--nolambda')
        commands.append(command)
        print(f'MACS3: q={q}; {len(samples)} sequential libraries; consensus >=2 replicates',flush=True)
        subprocess.run(command,cwd=folder,check=True)
    (out/'comparison.provenance.json').write_text(json.dumps({'source':str(source),'commands':commands,
        'qvalues':args.qvalues,'min_reps_consensus':2,'interpretation':'prespecified exploratory threshold sensitivity; not differential testing'},indent=2)+'\n')


if __name__=='__main__':
    main()
