from common import *
from metrics import enrichment,metrics
from scipy.special import expit
from scipy.stats import norm
from sklearn.metrics import roc_auc_score,log_loss,brier_score_loss
import joblib
start=time.monotonic();assert not (O/'results/ISO-001A_metrics.tsv').exists()
d=pd.read_parquet(O/'results/ISO-001A_development_OOF.parquet').set_index('promoter_group_id');meta=pd.read_parquet(O/'data/group_metadata_label_blind.parquet').set_index('promoter_group_id').loc[d.index];labels=pd.read_parquet(O/'data/development_E4_v1.parquet').set_index('promoter_group_id').loc[d.index];y=d.target.to_numpy();s2=d.S2_target.to_numpy();comp=d.dependency_component_id.to_numpy();clusters,ci=np.unique(comp,return_inverse=True);K=len(clusters)
assert len(d)==len(set(d.index)) and labels.evaluation_eligible.all() and np.array_equal(labels.primary_target,y)
models=[c for c in d if c not in ['target','S2_target','inner_fold','dependency_component_id']];legit=['B0','B2','SOURCE_AWARE','READ_UNION','ELASTIC_NET','RANDOM_FOREST'];prob=['ELASTIC_NET','RANDOM_FOREST','ISOCRED_MAP','ISOCRED_NO_ANNOTATION','ABL_CANDIDATE_REPLICATION','ABL_SOURCE_CUMULATIVE','ABL_MISSING_ZERO','CALLER_TALON','CALLER_FLAIR','CALLER_IsoQuant'];mm=[];ff=[]
for m in models:
 for ep,yy in [('PRIMARY',y),('S2_SAME_ELIGIBILITY',s2)]:
  r=metrics(yy,d[m].to_numpy());ok=np.isfinite(d[m]);r.update(model=m,endpoint=ep,clusters=len(np.unique(comp[ok])))
  if m in prob and ep=='PRIMARY':r['Brier']=brier_score_loss(yy,d[m]);r['log_loss']=log_loss(yy,np.clip(d[m],1e-12,1-1e-12))
  mm.append(r)
 den=num=0
 for f,z in d.groupby('inner_fold'):
  ok=np.isfinite(z[m]);z=z.loc[ok];auc=roc_auc_score(z.target,z[m]);n1=z.target.sum();n0=len(z)-n1;num+=auc*n1*n0;den+=n1*n0
  ff.append(dict(model=m,fold=f,conditional_AUC=auc,groups=len(z),prevalence=z.target.mean(),**{k:v for k,v in metrics(z.target.to_numpy(),z[m].to_numpy()).items() if k not in ['n','prevalence']}))
 ff.append(dict(model=m,fold='PAIR_WEIGHTED',conditional_AUC=num/den))
ms=pd.DataFrame(mm);primary=ms[ms.endpoint=='PRIMARY'].set_index('model');strongest=primary.loc[legit].enrichment.idxmax();point_delta=primary.at['ISOCRED_MAP','enrichment']-primary.at[strongest,'enrichment']
class FastLift:
 def __init__(self,y,s,subset=None):
  keep=np.isfinite(s) if subset is None else np.isfinite(s)&subset;ii=np.flatnonzero(keep);o=np.argsort(-s[ii],kind='stable');self.i=ii[o];self.y=y[self.i];v=s[self.i]
  # End of each exact tie run, so tied groups are all selected in bootstrap draws.
  self.end=np.empty(len(v),dtype=int);ends=np.r_[np.flatnonzero(v[:-1]!=v[1:]),len(v)-1];st=0
  for en in ends:self.end[st:en+1]=en;st=en+1
 def __call__(self,w):
  a=w[self.i];cs=np.cumsum(a);cp=np.cumsum(a*self.y)
  if len(cs)==0 or cs[-1]==0 or cp[-1]==0:return np.nan
  k=self.end[min(np.searchsorted(cs,.1*cs[-1]),len(cs)-1)];return (cp[k]/cs[k])/(cp[-1]/cs[-1])
fast={m:FastLift(y,d[m].to_numpy()) for m in models}
for m in models:assert abs(fast[m](np.ones(len(d)))-primary.at[m,'enrichment'])<1e-10
part={m:np.isfinite(d[m].to_numpy()) for m in ['PSQAN_STYLE','SQANTI3_RULE']};paired={m:(FastLift(y,d.ISOCRED_MAP.to_numpy(),ok),FastLift(y,d[m].to_numpy(),ok)) for m,ok in part.items()}
readcat=pd.cut(meta.read_support_union,[-1,0,4,19,99,np.inf],labels=['0','1-4','5-19','20-99','>=100']).astype(str)
known=np.select([meta.known_fraction==0,meta.known_fraction==1],['novel_only','known_only'],default='mixed')
share=pd.cut(meta.candidate_count,[0,1,4,np.inf],labels=['1','2-4','>=5']).astype(str)
op=pd.cut(meta.window_opportunity_bp,[0,101,125,np.inf],labels=['101','102-125','>125']).astype(str)
mono=np.select([meta.monoexon_fraction==0,meta.monoexon_fraction==1],['none','all'],default='mixed')
strata={'read_support':readcat.to_numpy(),'known_structure':known,'promoter_sharing':share.to_numpy(),'window_opportunity':op.to_numpy(),'monoexon':mono}
for c in ['TALON','FLAIR','IsoQuant']:
 v=meta['caller_'+c+'_candidate_fraction'];strata['caller_'+c]=np.select([v==0,v==1],['absent','all'],default='partial')
strata['joint_read_annotation_sharing']=np.array([a+'|'+b+'|'+c for a,b,c in zip(readcat,known,share)])
strows=[];control=[]
for axis,values in strata.items():
 for level in sorted(set(values)):
  keep=values==level;n=keep.sum();pos=int(y[keep].sum());nc=len(np.unique(comp[keep]));interpretable=n>=200 and pos>=20 and nc>=20
  if n==0 or pos==0 or pos==n:continue
  for m in models:
   if np.isfinite(d.loc[keep,m]).sum()<2:continue
   r=metrics(y[keep],d.loc[keep,m].to_numpy());strows.append(dict(axis=axis,level=level,model=m,positive_groups=pos,clusters=nc,interpretable=interpretable,**r))
  if axis=='joint_read_annotation_sharing' and interpretable:
   control.append((n,FastLift(y,d.ISOCRED_MAP.to_numpy(),keep),FastLift(y,d[strongest].to_numpy(),keep)))
control_n=sum(z[0] for z in control);point_control=sum(n*(a(np.ones(len(d)))-b(np.ones(len(d)))) for n,a,b in control)/control_n if control_n else np.nan
br={m:(d[m].to_numpy()-y)**2 for m in prob};boot=[];rng=np.random.default_rng(SEED)
for b in range(2000):
 cw=np.bincount(rng.integers(K,size=K),minlength=K);w=cw[ci].astype(float);vals={m:f(w) for m,f in fast.items()};best=max(vals[m] for m in legit)
 row=dict(iteration=b,ISOCRED_MAP=vals['ISOCRED_MAP'],strongest_baseline=best,delta_vs_strongest=vals['ISOCRED_MAP']-best)
 for m in models:row['lift__'+m]=vals[m]
 for m,(a,bb) in paired.items():row['delta_common__'+m]=a(w)-bb(w)
 for m in ['ABL_CANDIDATE_REPLICATION','ABL_SOURCE_CUMULATIVE','ABL_MISSING_ZERO']:row['Brier_gain__'+m]=float(w@(br[m]-br['ISOCRED_MAP'])/w.sum())
 row['controlled_delta']=sum(n*(a(w)-bb(w)) for n,a,bb in control)/control_n if control_n else np.nan
 boot.append(row)
 if b%500==0:print('paired cluster bootstrap',b,flush=True)
boot=pd.DataFrame(boot);summ=[]
for m in models:
 lo,hi=boot['lift__'+m].quantile([.025,.975]);summ.append(dict(comparison=m,point=primary.at[m,'enrichment'],CI_low=lo,CI_high=hi,estimand='TOP_DECILE_ENRICHMENT'))
lo,hi=boot.delta_vs_strongest.quantile([.025,.975]);summ.append(dict(comparison='ISOCRED_MINUS_STRONGEST',point=point_delta,CI_low=lo,CI_high=hi,estimand='PAIRED_DELTA_MAX_BASELINE_IN_EACH_BOOTSTRAP'))
for m,(a,b) in paired.items():
 lo,hi=boot['delta_common__'+m].quantile([.025,.975]);summ.append(dict(comparison='ISOCRED_MINUS_'+m,point=a(np.ones(len(d)))-b(np.ones(len(d))),CI_low=lo,CI_high=hi,estimand='PAIRED_DELTA_COMMON_SUPPORT',common_groups=int(part[m].sum())))
for m in ['ABL_CANDIDATE_REPLICATION','ABL_SOURCE_CUMULATIVE','ABL_MISSING_ZERO']:
 lo,hi=boot['Brier_gain__'+m].quantile([.025,.975]);summ.append(dict(comparison='Brier_GAIN_VS_'+m,point=float((br[m]-br['ISOCRED_MAP']).mean()),CI_low=lo,CI_high=hi,estimand='CONTROL_BRIER_MINUS_ISOCRED_BRIER'))
lo,hi=boot.controlled_delta.quantile([.025,.975]);summ.append(dict(comparison='READ_ANNOTATION_SHARING_CONTROLLED',point=point_control,CI_low=lo,CI_high=hi,estimand='N_WEIGHTED_WITHIN_JOINT_STRATUM_DELTA',common_groups=control_n))
tsv('results/ISO-001A_metrics.tsv',mm);tsv('results/ISO-001A_fold_metrics.tsv',ff);tsv('results/ISO-001A_cluster_CI.tsv',summ);tsv('results/ISO-001A_stratified_metrics.tsv',strows);parq('results/ISO-001A_cluster_bootstrap.parquet',boot)
js('results/scientific_gate_numbers.json',dict(groups=len(d),positives=int(y.sum()),S2_positives=int(s2.sum()),independent_components=K,strongest_same_information_baseline=strongest,IsoCred_enrichment=primary.at['ISOCRED_MAP','enrichment'],strongest_baseline_enrichment=primary.at[strongest,'enrichment'],delta=point_delta,controlled_delta=point_control,controlled_groups=control_n,comparisons=summ,primary_metrics=primary.reset_index().to_dict('records')))
# Fixed prior predictive spread for unseen components; NOT posterior credible intervals.
ev=pd.read_parquet(O/'data/group_evidence_label_blind.parquet').set_index('promoter_group_id').loc[d.index];base=pd.read_parquet(O/'data/group_baselines_label_blind.parquet').set_index('promoter_group_id').loc[d.index];x=pd.DataFrame(index=d.index)
for rep in ['R1','R2']:
 x['LR_read_'+rep]=np.log1p(ev['read_union_'+rep]);x['LR_junction_'+rep]=np.log1p(ev['junction_median_'+rep]);x['LR_junction_fraction_'+rep]=ev['junction_fraction_'+rep];x['resolution_read_shared_'+rep]=ev['read_crossgroup_fraction_'+rep]
x['annotation']=base.B2/2;x['resolution_candidate_count']=np.log1p(meta.candidate_count);x['resolution_TSS_count']=np.log1p(meta.exact_TSS_count);x['resolution_window']=np.log1p(meta.window_opportunity_bp);x['resolution_monoexon']=meta.monoexon_fraction;x['status_junction_NOT_APPLICABLE']=(ev.junction_unique_count==0).astype(float)
unc=d[['ISOCRED_MAP']].copy()
for f in range(5):
 a=joblib.load(O/f'models/ISOCRED_MAP_fold{f}.joblib');ix=d.inner_fold==f;assert list(x.columns)==a['features'];xx=a['preprocessor'].transform(x.loc[ix]);m=a['model'];pv=m.predict_proba(xx,comp[ix]);assert np.max(abs(pv-d.loc[ix,'ISOCRED_MAP'].to_numpy()))<1e-12
 z=m.theta[0]+xx@m.theta[1:1+m.p];unc.loc[ix,'unseen_component_prior_p10']=expit(z+.5*norm.ppf(.1));unc.loc[ix,'unseen_component_prior_p90']=expit(z+.5*norm.ppf(.9))
unc['uncertainty_semantics']='UNSEEN_COMPONENT_PRIOR_ONLY_NOT_CALIBRATED_CREDIBLE_INTERVAL';parq('results/ISO-001A_group_event_predictions_with_limited_uncertainty.parquet',unc.reset_index())
receipt('evaluate_comparison',start,evaluation_groups=len(d),bootstrap_components=K,bootstrap_draws=2000,OOF_model_reload_max_difference_bound=1e-12,primary_comparison_unit='PROMOTER_GROUP',no_final_or_external_metrics=True)
print(json.dumps(dict(strongest=strongest,delta=point_delta,controlled_delta=point_control),default=str),flush=True)
