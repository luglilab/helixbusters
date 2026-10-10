#!/usr/bin/env python3
"""Audit EXP1 descriptive gene lists with matched-depth support sensitivity."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--analysis-dir',type=Path,required=True)
parser.add_argument('--outdir',type=Path,required=True)
parser.add_argument('--iterations',type=int,default=50)
args=parser.parse_args()
if args.outdir.exists(): parser.error('Choose a new output directory')
if args.iterations<2: parser.error('Need at least two iterations')
args.outdir.mkdir(parents=True)
p=args.analysis_dir/'Annotation/GeneSignal'
provenance=json.loads((p/'gene_signal.provenance.json').read_text())
metadata=provenance['sample_metadata']; names=[r['sample'] for r in metadata]
audit=pd.read_csv(p/'gene_assignment.samples.tsv',sep='\t').set_index('sample').loc[names]
totals=audit.total_molecules.to_numpy(dtype=np.int64); target=int(totals.min())
groups={g:[i for i,r in enumerate(metadata) if r['group']==g] for g in sorted({r['group'] for r in metadata})}
summary={'matched_depth':target,'iterations':args.iterations,'seed':1729,'contexts':{},
 'interpretation':'Repeated molecule subsampling without replacement. Frequencies describe sensitivity to depth, not p-values, biological reproducibility probabilities or independent experiments. Other/unassigned molecules retain the original library denominator. GTF opportunity is not corrected for blacklist or mappability.'}
fig,axes=plt.subplots(1,2,figsize=(12,5))
for ax,context in zip(axes,('promoter','gene_body')):
 table=pd.read_csv(p/f'genes.{context}.counts.tsv',sep='\t')
 raw=table[names].to_numpy(dtype=np.int64)
 if (raw.sum(axis=0)>totals).any(): raise ValueError('Gene counts exceed retained library molecules')
 simulations={g:[] for g in groups}; stability={g:np.zeros(len(table),dtype=np.int64) for g in groups}
 rng=np.random.default_rng(1729)
 for iteration in range(args.iterations):
  sampled=np.zeros_like(raw)
  for i in range(len(names)):
   categories=np.append(raw[:,i],totals[i]-raw[:,i].sum())
   sampled[:,i]=rng.multivariate_hypergeometric(categories,target)[:-1]
  for g,columns in groups.items():
   settings=provenance['conditions'][f'{g}_{context}']
   supported=(sampled[:,columns]>=settings['minimum_molecules_per_replicate']).sum(axis=1)>=settings['minimum_replicates']
   stability[g]+=supported; simulations[g].append(int(supported.sum()))
 records=[]
 context_summary={}
 for index,(g,columns) in enumerate(groups.items()):
  settings=provenance['conditions'][f'{g}_{context}']
  observed=settings['supported_candidates']; values=np.array(simulations[g])
  context_summary[g]={'original_candidates':observed,'matched_depth_mean':float(values.mean()),
    'matched_depth_min':int(values.min()),'matched_depth_max':int(values.max()),
    'matched_depth_q025':float(np.quantile(values,.025)),'matched_depth_q975':float(np.quantile(values,.975))}
  ax.scatter(np.full(len(values),index+.15),values,alpha=.3,s=18)
  ax.scatter(index-.15,observed,marker='s',s=70,color='black')
  for iteration,n in enumerate(values): records.append({'condition':g,'iteration':iteration,'supported_candidates':n,'molecules_per_sample':target})
  features=table[['chrom','gene_id','gene_name','gene_biotype','uniquely_assignable_bp']].copy()
  features['support_frequency_equal_depth']=stability[g]/args.iterations
  features['original_supported']=(raw[:,columns]>=settings['minimum_molecules_per_replicate']).sum(axis=1)>=settings['minimum_replicates']
  cpm=raw[:,columns]/totals[columns][None,:]*1e6
  features['mean_CPM']=cpm.mean(axis=1); features['median_CPM']=np.median(cpm,axis=1)
  features.sort_values(['support_frequency_equal_depth','median_CPM'],ascending=False).to_csv(args.outdir/f'{g}.{context}.depth_sensitivity.tsv',sep='\t',index=False)
  candidates=pd.read_csv(p/f'{g}.{context}.candidate_genes.tsv',sep='\t')
  context_summary[g]['candidate_length_signal_spearman']=float(candidates.uniquely_assignable_bp.rank().corr(candidates.mean_CPM.rank()))
  context_summary[g]['three_replicate_support']=int((candidates.replicates_meeting_min_molecules==len(columns)).sum())
  context_summary[g]['top_genes']=candidates[['gene_name','uniquely_assignable_bp','pooled_molecules','mean_CPM','replicates_meeting_min_molecules']].head(10).to_dict('records')
  context_summary[g]['stable_candidates_frequency_ge080']=int(((stability[g]/args.iterations)>=.8).sum())
 pd.DataFrame(records).to_csv(args.outdir/f'{context}.equal_depth_candidate_counts.tsv',sep='\t',index=False)
 ax.set_xticks(range(len(groups)),list(groups));ax.set_title(context);ax.set_ylabel('Replicate-supported genes')
 ax.spines[['top','right']].set_visible(False)
 ax.scatter([],[],marker='s',color='black',label='Original libraries')
 ax.scatter([],[],alpha=.5,label=f'Equal depth: {target:,} molecules/sample')
 ax.legend(fontsize=8)
 summary['contexts'][context]=context_summary
fig.tight_layout();fig.savefig(args.outdir/'gene_support_depth_sensitivity.png',dpi=180);fig.savefig(args.outdir/'gene_support_depth_sensitivity.pdf');plt.close(fig)
summary['feature_percentages']=pd.read_csv(args.analysis_dir/'Annotation/dsb_feature_distribution.conditions.tsv',sep='\t')[['group','feature','mean_percentage','sd_percentage']].to_dict('records')
summary['assignment_audit']=audit.reset_index().to_dict('records')
(args.outdir/'review.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:{g:{key:value for key,value in r.items() if key!='top_genes'} for g,r in v.items()} for k,v in summary['contexts'].items()},indent=2))
