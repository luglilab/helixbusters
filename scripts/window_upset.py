#!/usr/bin/env python3
"""Overlap of replicate-supported fixed windows, with equal-depth sensitivity."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from helixbusters.differential import integer_counts, depth_sensitivity, sha256


def write_intersections(table, membership, folder, title):
    """Count original bins once each; adjacent bins are never merged."""
    groups=list(membership)
    matrix=np.column_stack([membership[g] for g in groups])
    counts=Counter(); bases=Counter()
    for i,row in enumerate(matrix):
        members=tuple(g for g,keep in zip(groups,row,strict=True) if keep)
        if members:
            counts[members]+=1
            bases[members]+=int(table.iloc[i]['end'])-int(table.iloc[i]['start'])
    combinations=sorted(set(counts)|{(g,) for g in groups}|{tuple(groups)},key=lambda c:(len(c),c))
    rows=[{'conditions':','.join(c),'windows':counts[c],'covered_bp':bases[c]} for c in combinations]
    pd.DataFrame(rows).to_csv(folder/'intersections.tsv',sep='\t',index=False)
    sets=pd.DataFrame([{'condition':g,'windows':int(membership[g].sum()),
                       'covered_bp':int((table.end-table.start)[membership[g]].sum())} for g in groups])
    sets.to_csv(folder/'set_sizes.tsv',sep='\t',index=False)
    fig,(bar,dots)=plt.subplots(2,1,figsize=(max(5,len(combinations)*1.1),4.5),sharex=True,
                               gridspec_kw={'height_ratios':[3,1]},layout='constrained')
    bar.bar(range(len(combinations)),[counts[c] for c in combinations],color='#4879ae')
    bar.set_ylabel('Supported genomic windows');bar.set_title(title)
    for x,c in enumerate(combinations):
        bar.text(x,counts[c],f'{counts[c]:,}',ha='center',va='bottom',fontsize=8)
        dots.scatter([x]*len(groups),range(len(groups)),color='#dddddd',s=35)
        indices=[i for i,g in enumerate(groups) if g in c]
        dots.plot([x]*len(indices),indices,'o-',color='black')
    dots.set_yticks(range(len(groups)),groups);dots.set_xticks([])
    dots.set_xlim(-.6,max(.6,len(combinations)-.4));dots.set_ylim(-.5,max(.5,len(groups)-.5))
    fig.savefig(folder/'conditions_upset.png',dpi=180);fig.savefig(folder/'conditions_upset.pdf');plt.close(fig)
    return rows


def run(analysis_dir,outdir,widths,minimum_molecules=5,minimum_replicates=2,
        iterations=50,seed=1729,minimum_frequency=.8):
    analysis_dir,outdir=Path(analysis_dir),Path(outdir)
    if outdir.exists():raise ValueError('Choose a new output directory')
    if minimum_molecules<1 or minimum_replicates<1 or iterations<2 or not 0<minimum_frequency<=1:
        raise ValueError('Invalid support or resampling settings')
    if not widths or any(width<1 for width in widths):raise ValueError('Window widths must be positive')
    metadata=pd.read_csv(analysis_dir/'analysis.samples.tsv',sep='\t',dtype=str,keep_default_na=False)
    names=metadata['sample'].tolist()
    if not names or len(set(names))!=len(names):raise ValueError('Sample IDs must be nonempty and unique')
    if metadata.duplicated(['group','replicate']).any():raise ValueError('One library per biological replicate required')
    for group,rows in metadata.groupby('group'):
        donors=rows.get('donor',pd.Series(dtype=str))
        donors=donors[donors!='']
        if donors.duplicated().any():raise ValueError('Duplicate donor within condition')
    groups={g:np.flatnonzero(metadata.group.to_numpy()==g) for g in sorted(metadata.group.unique())}
    if any(len(cols)<minimum_replicates for cols in groups.values()):
        raise ValueError('Support threshold exceeds biological replicate count of a condition')
    source=json.loads((analysis_dir/'analysis.summary.json').read_text())
    purpose=source.get('analysis_purpose','biological')
    support_label='technical input-level support' if purpose=='input_titration' else 'replicate support'
    totals=np.array([source['samples'][n]['molecules'] for n in names],dtype=np.int64)
    outdir.mkdir(parents=True)
    summary={'interpretation':'Descriptive overlap of replicate-supported signal; no enrichment test or differential claim.',
             'analysis_purpose':purpose,
             'parameters':{'minimum_molecules':minimum_molecules,'minimum_replicates':minimum_replicates,
                           'iterations':iterations,'seed':seed,'equal_depth_minimum_frequency':minimum_frequency},
             'library_molecules':dict(zip(names,map(int,totals))),
             'metadata_sha256':sha256(analysis_dir/'analysis.samples.tsv'),'windows':{}}
    if purpose=='input_titration':
        summary['interpretation']='Descriptive technical input-level overlap from one donor. Input levels are not biological replicates; no biological consensus or differential inference.'
    mqc={}
    for width in sorted(set(widths)):
        path=analysis_dir/f'windows_{width}.counts.tsv'
        table=pd.read_csv(path,sep='\t');raw=integer_counts(table,names)
        if not np.array_equal(raw.sum(axis=0),totals):raise ValueError('Window counts do not conserve molecular totals')
        if table.duplicated(['chrom','start','end']).any():raise ValueError('Duplicate genomic window')
        if ((table.end<=table.start)|(table.start<0)|(table.start%width!=0)|((table.end-table.start)>width)).any():
            raise ValueError('Invalid fixed window coordinates')
        folder=outdir/f'windows_{width}';folder.mkdir()
        output=table.copy()
        membership={g:(raw[:,cols]>=minimum_molecules).sum(axis=1)>=minimum_replicates for g,cols in groups.items()}
        for i,name in enumerate(names):
            output[name+'_CPM']=raw[:,i]/totals[i]*1e6 if totals[i] else 0
        for g,cols in groups.items():
            output[g+'_supporting_replicates']=(raw[:,cols]>=minimum_molecules).sum(axis=1)
            output[g+'_supported']=membership[g]
        full=folder/'full_depth';full.mkdir()
        full_rows=write_intersections(table,membership,full,f'{width/1000:g} kb | full-depth {support_label}')
        record={'input_sha256':sha256(path),'observed_windows':len(table),'full_depth':full_rows}
        for row in full_rows:mqc[f'{width}_full_'+row['conditions']]={'windows':row['windows'],'covered_bp':row['covered_bp']}
        if (totals>0).all():
            freq,draws,target=depth_sensitivity(raw,totals,groups,iterations,seed,minimum_molecules,minimum_replicates)
            stable={g:membership[g] & (freq[g]>=minimum_frequency) for g in groups}
            for g in groups:
                output[g+'_equal_depth_support_frequency']=freq[g]
                output[g+'_equal_depth_stable']=stable[g]
            equal=folder/'equal_depth';equal.mkdir()
            record['equal_depth']=write_intersections(table,stable,equal,f'{width/1000:g} kb | equal-depth {support_label}')
            record['equal_depth_molecules']=int(target)
            pd.DataFrame(draws).to_csv(equal/'support_draws.tsv',sep='\t',index=False)
            for row in record['equal_depth']:mqc[f'{width}_equal_'+row['conditions']]={'windows':row['windows'],'covered_bp':row['covered_bp']}
        else:record['equal_depth_status']='skipped: a sample has zero molecules'
        output.to_csv(folder/'window_membership.tsv',sep='\t',index=False)
        summary['windows'][str(width)]=record
    (outdir/'window_upset.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (outdir/'window_upset_mqc.json').write_text(json.dumps({'id':'helixbusters_window_overlap',
        'section_name':'Helixbusters supported-window overlap','description':summary['interpretation'],
        'plot_type':'table','pconfig':{'id':'helixbusters_window_overlap_table','title':'Window overlap and depth sensitivity'},'data':mqc},indent=2)+'\n')
    (outdir/'README.md').write_text(summary['interpretation'] + '\n' + 'Fixed windows are matched by chromosome/start/end and never merged. Full-depth support requires the stated molecule threshold in the stated number of libraries within each group (biological replicates only for biological analyses). Equal-depth support requires that rule in at least 80% of 50 draws by default. All libraries are sampled to the smallest library, without replacement; seed and settings are recorded. Neither support nor resampling establishes statistical enrichment or a condition difference. Unique means supported in only one condition under this rule, not absence of reads in the other condition.\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-dir',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--windows',type=int,nargs='+',default=[10000,50000,100000])
    p.add_argument('--minimum-molecules',type=int,default=5);p.add_argument('--minimum-replicates',type=int,default=2)
    p.add_argument('--iterations',type=int,default=50);p.add_argument('--seed',type=int,default=1729)
    p.add_argument('--minimum-frequency',type=float,default=.8)
    args=p.parse_args()
    run(args.analysis_dir,args.outdir,args.windows,args.minimum_molecules,args.minimum_replicates,
        args.iterations,args.seed,args.minimum_frequency)
