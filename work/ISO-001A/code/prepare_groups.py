"""Label-blind group features and dependence graph; no locked E4 read permission."""
from group_features import *
from collections import defaultdict
import sys,importlib.util
start=time.monotonic();reg=json.loads((W/'registry.json').read_text());extra=json.loads((W/'extra_registry.json').read_text())
split=pd.read_csv(I/'split',sep='\t').set_index('promoter_group_id');mem=pd.read_csv(I/'membership',sep='\t');c2g=mem.set_index('candidate_id').promoter_group_id
core=pd.read_parquet(I/'core').set_index('candidate_id');c=pd.read_csv(I/'candidate',sep='\t').set_index('candidate_id').loc[c2g.index]
x=group_design(core,mem);base=read_base(mem);meta=split.loc[x.index].copy();meta['candidate_count']=mem.groupby('promoter_group_id').size();meta['monoexon_fraction']=x.struct_monoexon;meta['known_fraction']=x.annotation_known_exact
par={g:g for g in split.index}
def find(a):
 while par[a]!=a:par[a]=par[par[a]];a=par[a]
 return a
def union(a,b):
 a,b=find(a),find(b)
 if a!=b:par[max(a,b)]=min(a,b)
def union_groups(groups):
 it=iter(groups);a=next(it,None)
 if a is not None:
  for b in it:union(a,b)
for comp,d in split.groupby('connected_component_id'):union_groups(d.index)
# Exact-member TSS opportunities; all are frozen coordinates, never CAGE-derived.
raw=pd.read_csv(I/'raw',sep='\t',usecols=['raw_model_key','caller','replicate','original_transcript_id','exact_structure_hash','TSS','strand','chromosome'])
cw=pd.read_csv(I/'crosswalk',sep='\t',usecols=['exact_structure_hash','cluster50_sha256']);cw['candidate_id']='ISOCRED50_'+cw.cluster50_sha256
raw=raw.merge(cw[['exact_structure_hash','candidate_id']],on='exact_structure_hash',validate='many_to_one');raw['promoter_group_id']=raw.candidate_id.map(c2g);raw=raw.dropna(subset=['promoter_group_id']);raw['position0']=raw.TSS-1
positions=raw[['promoter_group_id','chromosome','strand','position0']].drop_duplicates()
potential_edges=0
for ctx,d in positions.groupby(['chromosome','strand']):
 rows=d.sort_values('position0');active=[]
 for g,pos in zip(rows.promoter_group_id,rows.position0):
  active=[(pg,pp) for pg,pp in active if pp>=pos-100]
  for pg,pp in active:
   if pg!=g:union(pg,g);potential_edges+=1
  active.append((g,pos))
posg=positions.groupby('promoter_group_id').position0
meta['exact_TSS_count']=posg.nunique();meta['TSS_span']=posg.max()-posg.min()
def opportunity(s):
 a=sorted(set(s));return 101+sum(min(101,b-a) for a,b in zip(a,a[1:])) if a else 0
meta['window_opportunity_bp']=posg.apply(opportunity)
f=pd.DataFrame(index=x.index);rdetails=[]
for rep in ['R1','R2']:
 d=pq.ParquetFile(I/('read_'+rep)).read().to_pandas();d['promoter_group_id']=d.candidate_id.map(c2g);d=d.dropna(subset=['promoter_group_id'])
 for key,z in d[['read_key','promoter_group_id']].drop_duplicates().groupby('read_key',sort=False):
  if len(z)>1:union_groups(z.promoter_group_id)
 # Each physical read ID contributes once per group, even when many callers/members share it.
 unique=d[['read_key','promoter_group_id']].drop_duplicates();degree=unique.groupby('read_key').size();unique['shared']=unique.read_key.map(degree)>1
 f['read_union_'+rep]=unique.groupby('promoter_group_id').size().reindex(x.index,fill_value=0)
 f['read_crossgroup_fraction_'+rep]=unique.groupby('promoter_group_id').shared.mean().reindex(x.index)
 f['read_candidate_sum_'+rep]=d.groupby('promoter_group_id').size().reindex(x.index,fill_value=0)
 d['source_count']=d.source_caller_mask.map(int.bit_count)
 f['read_source_sum_'+rep]=d.groupby('promoter_group_id').source_count.sum().reindex(x.index,fill_value=0)
 for caller,bit in [('TALON',1),('FLAIR',2),('IsoQuant',4)]:
  dc=d[(d.source_caller_mask&bit)>0][['read_key','promoter_group_id']].drop_duplicates()
  f['read_'+caller+'_'+rep]=dc.groupby('promoter_group_id').size().reindex(x.index,fill_value=0)
 rdetails.append(dict(replicate=rep,candidate_links=len(d),unique_group_links=len(unique),caller_counted_edges=int(d.source_count.sum()),shared_read_ids=int((degree>1).sum())))
 print('reads completed',rep,flush=True)
keys=['contig','strand','donor_boundary','acceptor_boundary'];j=[]
for b in extra['bindings']:
 if b['input_id'].startswith('jlink_'):
  d=pq.ParquetFile(I/b['input_id']).read().to_pandas();d['promoter_group_id']=d.candidate_id.map(c2g);j.append(d.dropna(subset=['promoter_group_id']))
j=pd.concat(j,ignore_index=True);ju=j[['promoter_group_id']+keys].drop_duplicates()
for _,d in ju.groupby(keys,sort=False):
 if len(d)>1:union_groups(d.promoter_group_id)
f['junction_unique_count']=ju.groupby('promoter_group_id').size().reindex(x.index,fill_value=0)
f['junction_candidate_count']=j.groupby('promoter_group_id').size().reindex(x.index,fill_value=0)
for rep in ['R1','R2']:
 counts=pq.ParquetFile(I/('jcount_'+rep)).read().to_pandas()
 q=ju.merge(counts[keys+['unique_read_support']],on=keys,how='left',validate='many_to_one');q.unique_read_support=q.unique_read_support.fillna(0)
 f['junction_median_'+rep]=q.groupby('promoter_group_id').unique_read_support.median().reindex(x.index)
 f['junction_fraction_'+rep]=q.assign(supported=q.unique_read_support>0).groupby('promoter_group_id').supported.mean().reindex(x.index)
meta['dependency_component_id']=[find(g) for g in meta.index]
finaldeps=set(meta.loc[meta.split=='FINAL_BLIND_TEST','dependency_component_id']);meta['touches_final_dependency']=meta.dependency_component_id.isin(finaldeps)
meta['group_degree']=meta.groupby('dependency_component_id').dependency_component_id.transform('size')
v=pd.read_csv(I/'visible',sep='\t').set_index('promoter_group_id');meta['assay_eligibility']=v.assay_eligibility;meta['mapping_ambiguity']=v.mapping_ambiguity
# Bind and invoke the installed OFFICIAL SQANTI3 rules implementation, no ML training labels.
share=R/'env/ISO-000C-R1/sqanti3/share/sqanti3-6.0.2-1';sys.path.insert(0,str(share))
from src.utilities.filter.sqanti3_rules_filter import read_json_rules,apply_rules
rulesfile=share/'src/utilities/filter/filter_default.json';rules=read_json_rules(rulesfile);sqs=[];sq_audit=[]
for b in reg['bindings']:
 if not b['input_id'].startswith('sqanti_'):continue
 path=b['path'];caller=path.split('SQANTI3_')[1].split('_')[0];rep=path.split('SQANTI3_')[1].split('_')[1].split('/')[0]
 cols=['isoform','structural_category','exons','perc_A_downstream_TTS','RTS_stage','all_canonical','min_cov']
 d=pd.read_csv(I/b['input_id'],sep='\t',usecols=cols);u=d.drop(columns='isoform').drop_duplicates().copy()
 u['SQANTI3_RULE']=u.apply(lambda row:float(apply_rules(row,False,rules)=='Isoform'),axis=1)
 d=d.merge(u,on=cols[1:],validate='many_to_one');d['raw_model_key']=caller+'|'+rep+'|'+d.isoform
 q=raw[['raw_model_key','candidate_id','promoter_group_id']].merge(d[['raw_model_key','SQANTI3_RULE']],on='raw_model_key',validate='one_to_one')
 sqs.append(q);sq_audit.append(dict(caller=caller,replicate=rep,classification_rows=len(d),mapped_model_rows=len(q),rules_sha256=sha(rulesfile),code_sha256=sha(share/'src/utilities/filter/sqanti3_rules_filter.py'),ML_training_labels='NONE_RULES_MODE',E4_features_read=0))
sq=pd.concat(sqs,ignore_index=True).groupby('candidate_id').SQANTI3_RULE.mean();sqg=pd.DataFrame({'g':c2g,'s':sq.reindex(c2g.index)}).groupby('g').s.agg(['mean','count'])
base['SQANTI3_RULE']=sqg['mean'];meta['SQANTI3_candidate_coverage']=sqg['count']/meta.candidate_count
# PSQAN-style normalization uses neutral E1 support and frozen candidates, with NO new filtering.
# Ambiguous gene assignments are unscored, not set to a low score.
ps=pd.DataFrame(index=c.index);ps['gene']=c.primary_gene_id.where(c.gene_assignment_status=='KNOWN_UNIQUE')
for rep in ['r1','r2']:
 z=core.loc[c.index,'e1_unique_read_ids_'+rep].astype(float);den=z.groupby(ps.gene).transform('sum');ps[rep]=100*z/den.replace(0,np.nan)
ps['score']=ps[['r1','r2']].mean(axis=1,skipna=False);ps['g']=c2g;psg=ps.groupby('g').score.agg(['mean','count']);base['PSQAN_STYLE']=psg['mean'];meta['PSQAN_candidate_coverage']=psg['count']/meta.candidate_count
for caller in ['TALON','FLAIR','IsoQuant']:
 cs=raw[raw.caller==caller].groupby('promoter_group_id').candidate_id.nunique();meta['caller_'+caller+'_candidate_fraction']=cs.reindex(meta.index,fill_value=0)/meta.candidate_count
meta['read_support_union']=f.read_union_R1+f.read_union_R2
meta['read_duplication_ratio']=(f.read_source_sum_R1+f.read_source_sum_R2)/(f.read_union_R1+f.read_union_R2).replace(0,np.nan)
meta['junction_duplication_ratio']=f.junction_candidate_count/f.junction_unique_count.replace(0,np.nan)
dev=meta.index[meta.split=='DEVELOPMENT_POOL'];assert len(dev)==46408
parq('data/group_core_label_blind.parquet',x.loc[dev].reset_index());parq('data/group_evidence_label_blind.parquet',f.loc[dev].reset_index());parq('data/group_baselines_label_blind.parquet',base.loc[dev].reset_index());parq('data/group_metadata_label_blind.parquet',meta.loc[dev].reset_index())
# Dependency IDs/folds contain no evidence outcomes. Final groups themselves are not exported.
tsv('manifests/SQANTI3_rule_score_audit.tsv',sq_audit);tsv('manifests/evidence_deduplication_summary.tsv',rdetails)
js('receipts/label_blind_preparation.json',dict(development_groups=len(dev),frozen_final_groups=int((meta.split=='FINAL_BLIND_TEST').sum()),dev_touching_final_dependencies=int(meta.loc[dev].touches_final_dependency.sum()),dependency_components=meta.loc[dev].dependency_component_id.nunique(),potential_window_edges=potential_edges,max_dependency_groups=int(meta.group_degree.max()),all_shared_metadata_groups=int((v.mapping_ambiguity=='ALL_COMPATIBLE_CTSS_PROMOTER_GROUP_SHARED').sum()),ambiguous_5279_scope='historical outcome_status subset of 5912 all-shared metadata groups; no label-based reclassification',core_group_columns=x.columns.tolist(),dedup_features=f.columns.tolist()))
receipt('prepare_groups',start,development_groups=len(dev),E4_reads=0)
