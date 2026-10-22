"""Predeclared final statistics. Input is a frozen prediction table and exported final targets."""
from common import *
from sklearn.metrics import roc_auc_score,brier_score_loss
MAIN=['B0','B2','SOURCE_AWARE','READ_UNION','RANDOM_FOREST','ISOCRED_MAP']
OPTIONAL=['ELASTIC_NET','PSQAN_STYLE','SQANTI3_RULE']
CONTRASTS=[('ISOCRED_MAP','RANDOM_FOREST'),('RANDOM_FOREST','READ_UNION')]
def enrichment_at(y,s,weights=None,alpha=.1):
 y=np.asarray(y);s=np.asarray(s);w=np.ones(len(y)) if weights is None else np.asarray(weights)
 ok=np.isfinite(y)&np.isfinite(s)&(w>0);y=y[ok];s=s[ok];w=w[ok]
 if not len(y):return [np.nan]*4
 prev=w@y/w.sum()
 if prev==0:return [np.nan]*4
 order=np.argsort(-s,kind='stable');k=np.searchsorted(np.cumsum(w[order]),alpha*w.sum(),side='left');threshold=s[order[min(k,len(y)-1)]]
 top=s>=threshold;ppv=w[top]@y[top]/w[top].sum()
 return [ppv/prev,ppv,w[top].sum()/w.sum(),threshold]
def metric(y,s,alpha=.1):
 y=np.asarray(y);s=np.asarray(s);ok=np.isfinite(s)&np.isfinite(y)
 y=y[ok];s=s[ok];v=enrichment_at(y,s,alpha=alpha)
 return dict(n=len(y),positive=int(y.sum()),prevalence=float(y.mean()) if len(y) else np.nan,theoretical_enrichment_ceiling=float(1/y.mean()) if len(y) and y.mean()>0 else np.nan,enrichment=v[0],top_positive_rate=v[1],top_actual_coverage=v[2],threshold=v[3],AUC=float(roc_auc_score(y,s)) if len(np.unique(y))==2 else np.nan)
def strata(meta):
 return {
 'read_support':pd.cut(meta.read_support_union,[-1,0,4,19,99,np.inf],labels=['0','1-4','5-19','20-99','>=100']).astype(str),
 'annotation':np.where(meta.known_fraction==0,'novel_only',np.where(meta.known_fraction==1,'known_only','mixed')),
 'promoter_sharing':pd.cut(meta.candidate_count,[0,1,4,np.inf],labels=['1','2-4','>=5']).astype(str)
 }
def evaluate(pred,labels,meta,destination):
 assert not destination.exists();destination.mkdir()
 join=pred.join(labels[['primary_target']],validate='one_to_one')
 eligible=join.evaluation_eligible
 assert join.loc[eligible,'primary_target'].notna().all()
 common=eligible & join.primary_common_mask
 assert (common==eligible & join[MAIN].notna().all(axis=1)).all()
 d=join.loc[common];y=d.primary_target.to_numpy(dtype=int)
 comps,codes=np.unique(d.dependency_component_id,return_inverse=True)
 assert len(comps)>=20 and len(np.unique(y))==2
 rows=[];curve=[]
 for method in MAIN:
  rows.append(dict(analysis_set='PRIMARY_COMMON',method=method,eligible_denominator=int(eligible.sum()),scoring_coverage=float(common.sum()/eligible.sum()),independent_components=len(comps),**metric(y,d[method])))
  for alpha in [.05,.1,.2]:curve.append(dict(analysis_set='PRIMARY_COMMON',method=method,nominal_coverage=alpha,**metric(y,d[method],alpha)))
 for method in MAIN+OPTIONAL:
  use=eligible & np.isfinite(join[method])
  rows.append(dict(analysis_set='METHOD_SCORABLE_ELIGIBLE',method=method,eligible_denominator=int(eligible.sum()),scoring_coverage=float(use.sum()/eligible.sum()),independent_components=int(join.loc[use,'dependency_component_id'].nunique()),**metric(join.loc[use,'primary_target'],join.loc[use,method])))
 paired=[]
 for other in OPTIONAL:
  use=common & np.isfinite(join[other])
  if not use.any():continue
  for method in ['RANDOM_FOREST','ISOCRED_MAP',other]:paired.append(dict(analysis_set='COMMON_WITH_'+other,method=method,**metric(join.loc[use,'primary_target'],join.loc[use,method])))
 # A single dependency resample index per draw for all primary methods/coverages.
 rng=np.random.default_rng(SEED);draws=[];index_hash=hashlib.sha256()
 for b in range(2000):
  chosen=rng.integers(0,len(comps),len(comps),dtype=np.int64);index_hash.update(chosen.tobytes())
  weights=np.bincount(chosen,minlength=len(comps))[codes]
  row={'bootstrap':b}
  for method in MAIN:
   for alpha in [.05,.1,.2]:row[method+'__'+str(alpha)]=enrichment_at(y,d[method].to_numpy(),weights,alpha)[0]
  for a,z in CONTRASTS:row[a+'__minus__'+z]=row[a+'__0.1']-row[z+'__0.1']
  draws.append(row)
 boot=pd.DataFrame(draws);assert np.isfinite(boot.to_numpy()).all()
 intervals=[]
 for method in MAIN:
  for alpha in [.05,.1,.2]:
   values=boot[method+'__'+str(alpha)];point=metric(y,d[method],alpha)['enrichment']
   intervals.append(dict(kind='METHOD_COVERAGE',name=method,nominal_coverage=alpha,point=point,lower95=float(values.quantile(.025)),upper95=float(values.quantile(.975)),lower97_5=float(values.quantile(.0125)),upper97_5=float(values.quantile(.9875))))
 for a,z in CONTRASTS:
  values=boot[a+'__minus__'+z];point=metric(y,d[a])['enrichment']-metric(y,d[z])['enrichment']
  intervals.append(dict(kind='PRIMARY_CONTRAST',name=a+' minus '+z,nominal_coverage=.1,point=point,lower95=float(values.quantile(.025)),upper95=float(values.quantile(.975)),lower97_5=float(values.quantile(.0125)),upper97_5=float(values.quantile(.9875))))
 strat=[]
 sm=meta.loc[d.index]
 for axis,groups in strata(sm).items():
  groups=pd.Series(groups,index=sm.index)
  for level in sorted(groups.unique()):
   part=d.loc[groups==level];nc=part.dependency_component_id.nunique();pos=int(part.primary_target.sum())
   interpretable=len(part)>=200 and pos>=20 and nc>=20
   for method in MAIN:strat.append(dict(axis=axis,level=level,method=method,independent_components=nc,descriptive_interpretability_threshold_met=interpretable,**metric(part.primary_target,part[method])))
 tables={'metrics.tsv':rows,'coverage_curves.tsv':curve,'partial_common_support.tsv':paired,'cluster_intervals.tsv':intervals,'stratified_metrics.tsv':strat}
 for name,table in tables.items():pd.DataFrame(table).to_csv(destination/name,sep='\t',index=False)
 pq.write_table(pa.Table.from_pandas(boot,preserve_index=False),destination/'cluster_bootstrap.parquet',compression='zstd')
 info=dict(evidence_status='FINAL_HELDOUT_WITHIN_HEPG2',membership=len(join),eligible=int(eligible.sum()),primary_common=int(common.sum()),primary_positive=int(y.sum()),independent_components=len(comps),prevalence=float(y.mean()),theoretical_enrichment_ceiling=float(1/y.mean()),bootstrap_draws=2000,bootstrap_index_stream_sha256=index_hash.hexdigest(),seed=SEED,primary_contrasts=[a+' minus '+z for a,z in CONTRASTS],interval_policy='95% descriptive and nominal 97.5% Bonferroni-family intervals for the two primary contrasts; no p-values or significance-driven selection',no_final_refit=True)
 (destination/'summary.json').write_text(json.dumps(info,indent=2)+'\n')
 return info
