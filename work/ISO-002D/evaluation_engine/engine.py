"""ISO-002D memory-only paired evaluation engine. No file/network/model I/O."""
from __future__ import annotations
import hashlib,json,math
from dataclasses import dataclass
import numpy as np
METHODS=('RF10','EN10','ANNOTATION_ONLY','JUNCTION_MEDIAN_SUPPORT','WINDOW_OPPORTUNITY_ONLY')
SUPERVISED=frozenset(('RF10','EN10'))
CONTRASTS=(('RF10','EN10'),('RF10','ANNOTATION_ONLY'),('RF10','JUNCTION_MEDIAN_SUPPORT'),
 ('RF10','WINDOW_OPPORTUNITY_ONLY'),('EN10','ANNOTATION_ONLY'),('EN10','JUNCTION_MEDIAN_SUPPORT'),('EN10','WINDOW_OPPORTUNITY_ONLY'))
FRACTIONS=(.05,.1,.2)
SEED=20260903
BOOTSTRAP_DRAWS=2000
PROBABILITY_CLIP_EPSILON=np.finfo(np.float64).eps
EXTERNAL_ENDPOINTS=('primary','technical_pair_1','technical_pair_2','technical_pair_3','secondary_position_concordant_5bp')
STRATA=frozenset(('annotation_score','candidate_count','window_opportunity_bp','monoexon_fraction','junction_support_sum'))
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
 if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
 if isinstance(x,(bool,np.bool_)):return bool(x)
 if isinstance(x,(int,np.integer)):return int(x)
 return x
def token_hash(values):
 return hashlib.sha256(json.dumps(list(values),ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
def strict_bool(value,name):
 if isinstance(value,(bool,np.bool_)):return bool(value)
 if isinstance(value,(int,np.integer)) and value in (0,1):return bool(value)
 if isinstance(value,str) and value.lower() in ('true','false'):return value.lower()=='true'
 raise ValueError('INVALID_BOOLEAN:'+name)
def indexed(rows,role):
 out={}
 for r in rows:
  gid=r.get('promoter_group_id')
  if not isinstance(gid,str) or not gid or gid.strip()!=gid:raise ValueError('INVALID_GROUP_ID:'+role)
  if gid in out:raise ValueError('DUPLICATE_GROUP_ID:'+role+':'+gid)
  out[gid]=r
 return out
@dataclass
class Prepared:
 ids:list[str]
 eligible:np.ndarray
 labels:dict[str,np.ndarray]
 scores:dict[str,np.ndarray]
 components:list[str]
 mono:np.ndarray
 folds:list[str]|None
 case:str
 label_columns:dict[str,str]
def prepare_records(label_rows,prediction_rows,metadata_rows,*,case,label_columns,eligible_column='primary_eligible',
 component_column='dependency_component_id',mono_column='junction_not_applicable',fold_column=None):
 """Require 1:1 complete joins. The primary eligible mother set is never rebuilt per endpoint."""
 if case not in ('development','external'):raise ValueError('INVALID_CASE')
 labels=indexed(label_rows,'labels');pred=indexed(prediction_rows,'predictions');meta=indexed(metadata_rows,'metadata')
 if set(labels)!=set(pred) or set(labels)!=set(meta):raise ValueError('NONIDENTICAL_GROUP_ID_SETS')
 if set(label_columns)!=(set(EXTERNAL_ENDPOINTS) if case=='external' else {'primary'}):raise ValueError('WRONG_FIXED_ENDPOINT_SET')
 ids=sorted(labels)
 eligible=np.array([strict_bool(labels[g][eligible_column],eligible_column) for g in ids],dtype=bool)
 mono=np.array([strict_bool(meta[g][mono_column],mono_column) for g in ids],dtype=bool)
 y={}
 for ep,col in label_columns.items():
  vals=np.array([np.nan if labels[g].get(col) is None else float(labels[g][col]) for g in ids])
  if np.isinf(vals).any() or not set(vals[np.isfinite(vals)])<={0.,1.}:raise ValueError('NONBINARY_LABEL:'+ep)
  if not np.isfinite(vals[eligible]).all():raise ValueError('MISSING_LABEL_IN_PRIMARY_ELIGIBLE_SET:'+ep)
  y[ep]=vals
 if case=='external' and np.any(y['secondary_position_concordant_5bp'][eligible]>y['primary'][eligible]):
  raise ValueError('SECONDARY5BP_NOT_SUBSET_OF_PRIMARY')
 scores={}
 for method in METHODS:
  if any(method not in pred[g] for g in ids):raise ValueError('MISSING_METHOD_COLUMN:'+method)
  v=np.array([np.nan if pred[g].get(method) is None else float(pred[g][method]) for g in ids],dtype=np.float64)
  if method in SUPERVISED and np.any((v[np.isfinite(v)]<0)|(v[np.isfinite(v)]>1)):raise ValueError('INVALID_PROBABILITY:'+method)
  scores[method]=v
 comps=[meta[g].get(component_column) for g in ids]
 if any(not isinstance(c,str) or not c or c.strip()!=c for c in comps):raise ValueError('INVALID_COMPONENT_ID')
 folds=None
 if case=='development':
  if not fold_column:raise ValueError('DEVELOPMENT_FOLD_COLUMN_REQUIRED')
  folds=[str(meta[g][fold_column]) for g in ids]
  if any(f in ('None','nan','NA','') for f in folds):raise ValueError('MISSING_FOLD')
 return Prepared(ids,eligible,y,scores,comps,mono,folds,case,dict(label_columns))
class RankPlan:
 """One fixed finite-score population. Sorting/tie buckets are label-independent."""
 def __init__(self,scores,mask):
  self.indices=np.flatnonzero(mask);self.s=np.asarray(scores,dtype=np.float64)[self.indices]
  if not np.isfinite(self.s).all():raise ValueError('RANK_PLAN_REQUIRES_FINITE_SCORES')
  self.order=np.argsort(self.s,kind='stable')
  self.sorted_s=self.s[self.order]
  self.starts=np.flatnonzero(np.r_[True,self.sorted_s[1:]!=self.sorted_s[:-1]]) if len(self.s) else np.array([],dtype=int)
 def compute(self,labels,all_weights,denominator_weights,probability=False):
  y=np.asarray(labels,dtype=np.float64)[self.indices];w=np.asarray(all_weights,dtype=np.float64)[self.indices]
  if np.any(~np.isfinite(w)) or np.any(w<0):raise ValueError('INVALID_WEIGHTS')
  active=w>0
  if not np.isfinite(y[active]).all() or not set(y[active])<={0.,1.}:raise ValueError('INVALID_ACTIVE_LABELS')
  total=float(w.sum());den=float(denominator_weights)
  r={'active_unique_n':int(active.sum()),'weighted_n':total,'scoring_coverage':total/den if den>0 else np.nan}
  p=float(np.dot(w[active],y[active])/total) if total>0 else np.nan
  r.update(prevalence=p,prevalence_ceiling=1/p if p>0 else np.nan,AUC=np.nan)
  if total>0:
   yy=np.nan_to_num(y,nan=0.0);sw=w[self.order]
   pos=np.add.reduceat((w*yy)[self.order],self.starts);mass=np.add.reduceat(sw,self.starts);neg=mass-pos
   P=float(pos.sum());N=float(neg.sum())
   if P>0 and N>0:r['AUC']=float(np.dot(pos,np.cumsum(neg)-neg+.5*neg)/(P*N))
   desc_mass=mass[::-1];desc_pos=pos[::-1];cumulative=np.cumsum(desc_mass)
  for f in FRACTIONS:
   key=format(f,'.2g')
   if total>0:
    j=min(int(np.searchsorted(cumulative,f*total,side='left')),len(cumulative)-1)
    selected=float(cumulative[j]);pt=float(np.cumsum(desc_pos)[j]/selected)
    cut=float(self.sorted_s[self.starts][::-1][j]);cov=selected/total
    r.update({f'enrichment@{key}':pt/p if p>0 else np.nan,f'PPV@{key}':pt,
     f'actual_selection_coverage@{key}':cov,f'threshold@{key}':cut,
     f'actual_coverage_enrichment_ceiling@{key}':min(1/p,1/cov) if p>0 else np.nan})
   else:
    for k in ['enrichment','PPV','actual_selection_coverage','threshold','actual_coverage_enrichment_ceiling']:r[f'{k}@{key}']=np.nan
  if probability:
   if np.any((self.s<0)|(self.s>1)):raise ValueError('PROBABILITY_OUTSIDE_0_1')
   if total>0:
    q=np.clip(self.s,PROBABILITY_CLIP_EPSILON,1-PROBABILITY_CLIP_EPSILON)
    r['Brier']=float(np.dot(w[active],(y[active]-self.s[active])**2)/total)
    r['log_loss']=float(-np.dot(w[active],y[active]*np.log(q[active])+(1-y[active])*np.log1p(-q[active]))/total)
   else:r.update(Brier=np.nan,log_loss=np.nan)
  return r
def weighted_metrics(y,s,w=None,probability=False):
 """Pure single-score adapter; no score flipping. Nonfinite scores excluded with explicit coverage."""
 y=np.asarray(y,dtype=np.float64);s=np.asarray(s,dtype=np.float64)
 w=np.ones(len(y)) if w is None else np.asarray(w,dtype=np.float64)
 if y.shape!=s.shape or y.shape!=w.shape or y.ndim!=1:raise ValueError('SHAPE')
 if np.any(np.isinf(y)) or not set(y[np.isfinite(y)])<={0.,1.}:raise ValueError('NONBINARY')
 if np.any(~np.isfinite(w)) or np.any(w<0):raise ValueError('INVALID_WEIGHTS')
 den=np.isfinite(y)&(w>0);mask=den&np.isfinite(s)
 return RankPlan(s,mask).compute(y,w,float(w[den].sum()),probability)
def interval(samples,confidence=.95):
 values=np.asarray(samples,dtype=np.float64);ok=values[np.isfinite(values)];alpha=1-confidence
 if not len(ok):return dict(status='UNDEFINED_NO_VALID_DRAWS',low=None,high=None,confidence=confidence,
  total_draws=len(values),valid_draws=0,invalid_draws=len(values),lower_tail_expected_draws=0.0)
 low,high=np.quantile(ok,[alpha/2,1-alpha/2],method='linear')
 return dict(status='DEFINED',low=float(low),high=float(high),confidence=confidence,total_draws=len(values),
  valid_draws=len(ok),invalid_draws=len(values)-len(ok),quantile_method='linear',
  lower_quantile=alpha/2,upper_quantile=1-alpha/2,lower_tail_expected_draws=len(ok)*alpha/2,
  empirical_probability_grid=1/len(ok),order_statistic_reference_tail_scale=1/(len(ok)+1),
  tail_scale_interpretation='REFERENCE_SCALE_ONLY_NOT_INTERPOLATED_QUANTILE_ATTAINABILITY_LIMIT',
  coverage_claim='NOMINAL_PERCENTILE_BOOTSTRAP_APPROXIMATION_NOT_EXACT_FWER_GUARANTEE')
def component_draws(components,eligible,draws=BOOTSTRAP_DRAWS,seed=SEED):
 """Shared mother-population component draws; each member keeps component multiplicity."""
 eligible=np.asarray(eligible,dtype=bool);unique=sorted({components[i] for i in np.flatnonzero(eligible)})
 lookup={c:i for i,c in enumerate(unique)};index=np.array([lookup.get(c,-1) for c in components],dtype=np.int64)
 rng=np.random.Generator(np.random.PCG64(seed))
 for b in range(draws):
  counts=np.bincount(rng.integers(0,len(unique),size=len(unique)),minlength=len(unique)) if unique else np.array([],dtype=np.int64)
  w=np.zeros(len(components),dtype=np.float64)
  w[eligible]=counts[index[eligible]]
  yield b,counts,w
def population_info(data,mask,mother):
 ids=[data.ids[i] for i in np.flatnonzero(mask)]
 return dict(population_id=data.case+':'+token_hash(ids),n=len(ids),component_n=len({data.components[i] for i in np.flatnonzero(mask)}),
  eligible_mother_n=int(mother.sum()),eligible_mother_population_id=data.case+':'+token_hash(data.ids[i] for i in np.flatnonzero(mother)),
  excluded_nonfinite_n=int((mother&~mask).sum()),excluded_mono_n=int((mother&~mask&data.mono).sum()),
  included_mono_n=int((mask&data.mono).sum()),scoring_coverage=len(ids)/int(mother.sum()) if mother.any() else np.nan)
def eligible_endpoint_summary(data,endpoint,mask):
 n=int(mask.sum());positive=int(data.labels[endpoint][mask].sum())
 p=positive/n if n else np.nan
 return dict(n=n,positives=positive,prevalence=p,prevalence_ceiling=1/p if p>0 else np.nan,
  prevalence_status='DEFINED' if n else 'UNDEFINED_EMPTY_ELIGIBLE_POPULATION',
  ceiling_status='DEFINED' if p>0 else 'UNDEFINED_EMPTY_OR_ZERO_POSITIVES')
def status_map(point):
 out={}
 for k,v in point.items():
  if isinstance(v,(float,np.floating)) and not np.isfinite(v):
   if not point['weighted_n']:out[k]='UNDEFINED_EMPTY_SCORING_SUPPORT'
   elif k=='AUC':out[k]='UNDEFINED_SINGLE_CLASS'
   elif k.startswith('enrichment') or 'ceiling' in k:out[k]='UNDEFINED_ZERO_POSITIVES'
   else:out[k]='UNDEFINED'
 return out
def evaluate_records(label_rows,prediction_rows,metadata_rows,*,case,label_columns,contrasts,
 eligible_column='primary_eligible',component_column='dependency_component_id',mono_column='junction_not_applicable',
 fold_column=None,descriptive_strata=None,return_bootstrap_samples=False):
 """Memory-only API. Real execution is exclusively the gated parent adapter's responsibility."""
 if tuple(tuple(c) for c in contrasts)!=CONTRASTS:raise ValueError('CONTRAST_LIST_NOT_THE_FROZEN_7_ORDER')
 data=prepare_records(label_rows,prediction_rows,metadata_rows,case=case,label_columns=label_columns,
  eligible_column=eligible_column,component_column=component_column,mono_column=mono_column,fold_column=fold_column)
 scopes=[('pooled',data.eligible.copy(),True)]
 if data.folds is not None:
  for fold in sorted(set(data.folds)):scopes.append(('fold:'+fold,data.eligible&(np.asarray(data.folds)==fold),True))
 for name,mapping in (descriptive_strata or {}).items():
  if name not in STRATA:raise ValueError('UNDECLARED_STRATUM')
  if set(mapping)!=set(data.ids):raise ValueError('STRATUM_ID_SET_MISMATCH')
  values=[str(mapping[g]) if mapping[g] is not None else 'NA' for g in data.ids]
  for level in sorted(set(values)):scopes.append(('stratum:'+name+':'+level,data.eligible&(np.asarray(values)==level),False))
 plans={};planmask={};native=[];pairs=[];point={};arrays={}
 ones=np.ones(len(data.ids))
 bootkeys=['prevalence','scoring_coverage','AUC']+[f'{metric}@{format(f,".2g")}' for f in FRACTIONS for metric in ['enrichment','PPV','actual_selection_coverage']]+['Brier','log_loss']
 def plan(scope,ep,method,mask,do_boot):
  pid=(scope,ep,method,token_hash(np.flatnonzero(mask).tolist()))
  if pid not in plans:
   plans[pid]=RankPlan(data.scores[method],mask);planmask[pid]=mask
   mother=next(m for s,m,b in scopes if s==scope)
   point[pid]=plans[pid].compute(data.labels[ep],ones,float(mother.sum()),method in SUPERVISED)
   if do_boot:arrays[pid]={k:np.full(BOOTSTRAP_DRAWS,np.nan) for k in bootkeys if k in point[pid]}
  return pid
 for scope,mother,do_boot in scopes:
  for ep in label_columns:
   for method in METHODS:
    mask=mother&np.isfinite(data.scores[method]);pid=plan(scope,ep,method,mask,do_boot)
    native.append(dict(scope=scope,endpoint=ep,method=method,pid=pid,population=population_info(data,mask,mother),bootstrap=do_boot))
   for left,right in CONTRASTS:
    mask=mother&np.isfinite(data.scores[left])&np.isfinite(data.scores[right])
    lp=plan(scope,ep,left,mask,do_boot);rp=plan(scope,ep,right,mask,do_boot)
    pairs.append(dict(scope=scope,endpoint=ep,left=left,right=right,lpid=lp,rpid=rp,
     population=population_info(data,mask,mother),bootstrap=do_boot))
 drawhash=hashlib.sha256();denominators={scope:mother for scope,mother,boot in scopes if boot}
 for b,counts,w in component_draws(data.components,data.eligible):
  drawhash.update(np.asarray(counts,dtype='<i8').tobytes())
  den={s:float(w[m].sum()) for s,m in denominators.items()}
  for pid,metrics in arrays.items():
   scope,ep,method,_=pid
   r=plans[pid].compute(data.labels[ep],w,den[scope],method in SUPERVISED)
   for k in metrics:metrics[k][b]=r[k]
 results=[];contrasts_out=[]
 for row in native:
  pid=row.pop('pid');p=point[pid]
  row.update(metrics=p,undefined=status_map(p),CI95={k:interval(v) for k,v in arrays.get(pid,{}).items()},
   inference_role='DESCRIPTIVE_STRATUM_NO_INTERVAL' if not row['bootstrap'] else 'NOMINAL_DESCRIPTIVE',
   scoring_population='METHOD_OWN_FINITE_SUPPORT_WITH_COMMON_PRIMARY_ELIGIBILITY')
  results.append(row)
 for row in pairs:
  lp=row.pop('lpid');rp=row.pop('rpid');lpoint=point[lp];rpoint=point[rp]
  common=[k for k in bootkeys if k in lpoint and k in rpoint and k not in ('prevalence','scoring_coverage')]
  for metric in common:
   delta=lpoint[metric]-rpoint[metric]
   out={k:v for k,v in row.items() if k!='bootstrap'}
   out.update(metric=metric,left_estimate=lpoint[metric],right_estimate=rpoint[metric],difference_left_minus_right=delta,
    left_metrics=lpoint,right_metrics=rpoint,
    favorable_direction='NEGATIVE' if metric in ('Brier','log_loss') else 'POSITIVE' if metric.startswith(('enrichment','PPV')) or metric=='AUC' else 'NOT_AN_ACCURACY_DIRECTION',
    undefined_reason=None if np.isfinite(delta) else 'PAIRED_METRIC_UNDEFINED_EMPTY_SINGLE_CLASS_OR_ZERO_POSITIVES')
   if row['bootstrap']:
    v=arrays[lp][metric]-arrays[rp][metric];out['CI95']=interval(v)
    if row['scope']=='pooled' and row['endpoint']=='primary' and metric=='enrichment@0.1':
     out['Bonferroni_CI']=interval(v,confidence=1-.05/len(CONTRASTS))
     out['Bonferroni_family_id']=case+':pooled:primary:enrichment@0.1:7_prespecified_contrasts'
     out['family_contrasts']=7
    if return_bootstrap_samples:out['paired_bootstrap_differences']=v
   else:out.update(CI95=None,inference_role='DESCRIPTIVE_STRATUM_NO_INTERVAL')
   contrasts_out.append(out)
 unique=sorted({data.components[i] for i in np.flatnonzero(data.eligible)})
 out=dict(schema_version='1.0',case=case,methods=list(METHODS),contrasts=[list(c) for c in CONTRASTS],label_columns=data.label_columns,
  eligible_endpoint_summaries={scope:{ep:eligible_endpoint_summary(data,ep,mother) for ep in label_columns}
   for scope,mother,do_boot in scopes},
  primary_eligible_population=population_info(data,data.eligible,data.eligible),native_method_results=results,paired_contrast_results=contrasts_out,
  bootstrap=dict(draws=BOOTSTRAP_DRAWS,seed=SEED,generator='numpy.Generator(PCG64)',component_order='LEXICOGRAPHIC',
   component_count=len(unique),component_order_sha256=token_hash(unique),component_count_draw_stream_sha256=drawhash.hexdigest(),
   resampling_population='COMMON_PRIMARY_ELIGIBLE_MOTHER_COMPONENTS;FIXED_PAIR_OR_METHOD_FINITE_MASK_APPLIED_TO_SHARED_DRAW',
   weight_rule='EACH_GROUP_WEIGHT_EQUALS_NUMBER_OF_TIMES_ITS_COMPONENT_WAS_DRAWN;NEVER_COLLAPSE_REPEATS_TO_BOOLEAN',
   biological_interpretation='DEPENDENCY_COMPONENT_RESAMPLING;NOT_BIOLOGICAL_DONOR_OR_LIBRARY_REPLICATION'),
  interpretation=dict(primary='PROMOTER_CTSS_SUPPORT_NOT_FUNCTION',external='BOTH2_OF2_BIOSAMPLE_POOLS',development='INHERITED2_OF3',
   cross_system='ABSOLUTE_ENRICHMENT_NOT_COMPARABLE_AS_IDENTICAL_ENDPOINT_SCALE',
   finite_support='CONTRASTS_WITH_DIFFERENT_POPULATION_IDS_ARE_DIFFERENT_FINITE_SUPPORT_ESTIMANDS;NO_CROSS_POPULATION_WINNER',
   intervals='95_PERCENT_NOMINAL;BONFERRONI_MARGINAL99.2857142857_ONLY_7_POOLED_PRIMARY_ENRICHMENT_CONTRASTS_PER_DATASET',
   multiplicity='NO_CROSS_DATASET_OR_ALL_ENDPOINT_FWER_CLAIM;FINITE_BOOTSTRAP_TAIL_MONTE_CARLO_LIMITATION',
   calibration='BRIER_LOGLOSS_SUPERVISED_ONLY;NO_DIAGNOSTIC_PROBABILITY_INTERPRETATION',score_flipping=False,
   sensitivity='ALL3_TECHNICAL_PAIRS_AND5BP_REPORTED_ON_FIXED_PRIMARY_MASK;NO_BEST_PAIR_SELECTION'),
  runtime=dict(numpy_version=np.__version__,probability_clip_epsilon=float(PROBABILITY_CLIP_EPSILON)),
  no_input_file_or_network_io_in_engine=True)
 return clean(out)
