from group_features import *
from iso000g_model_utils import weighted_ecdf
from metrics import enrichment,metrics
from minimal_map import HierarchicalMAP
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
import warnings,joblib,sys
start=time.monotonic();contract=O/'contracts/ISO-001A_real_data_comparison_contract.yaml';freeze=json.loads((O/'receipts/protocol_freeze.json').read_text());assert sha(contract)==freeze['contract_sha256']
label=pd.read_parquet(O/'data/development_E4_v1.parquet').set_index('promoter_group_id');meta=pd.read_parquet(O/'data/group_metadata_label_blind.parquet').set_index('promoter_group_id')
core=pd.read_parquet(O/'data/group_core_label_blind.parquet').set_index('promoter_group_id');ev=pd.read_parquet(O/'data/group_evidence_label_blind.parquet').set_index('promoter_group_id');base=pd.read_parquet(O/'data/group_baselines_label_blind.parquet').set_index('promoter_group_id')
assert label.index.equals(meta.index) and not label.index.duplicated().any()
eligible=label.evaluation_eligible;g=meta.index[eligible];meta=meta.loc[g];ev=ev.loc[g];base=base.loc[g];core=core.loc[g];y=label.loc[g,'primary_target'].to_numpy(dtype=int);fold=meta.inner_fold.to_numpy();comp=meta.dependency_component_id.to_numpy()
assert (meta.split=='DEVELOPMENT_POOL').all() and not meta.touches_final_dependency.any()
def minimal(mode='proper'):
 x=pd.DataFrame(index=g)
 for rep in ['R1','R2']:
  count=ev['read_source_sum_'+rep] if mode=='cumulative' else ev['read_union_'+rep]
  if mode.startswith('caller_'):count=ev['read_'+mode.split('_',1)[1]+'_'+rep]
  x['LR_read_'+rep]=np.log1p(count);x['LR_junction_'+rep]=np.log1p(ev['junction_median_'+rep]);x['LR_junction_fraction_'+rep]=ev['junction_fraction_'+rep]
  x['resolution_read_shared_'+rep]=ev['read_crossgroup_fraction_'+rep]
 x['annotation']=base.B2/2
 x['resolution_candidate_count']=np.log1p(meta.candidate_count);x['resolution_TSS_count']=np.log1p(meta.exact_TSS_count);x['resolution_window']=np.log1p(meta.window_opportunity_bp)
 x['resolution_monoexon']=meta.monoexon_fraction
 x['status_junction_NOT_APPLICABLE']=(ev.junction_unique_count==0).astype(float)
 if mode=='no_annotation':x=x.drop(columns='annotation')
 if mode=='zero':x=x.drop(columns='status_junction_NOT_APPLICABLE').fillna(0)
 return x
mi=minimal();X=pd.concat([core,mi.add_prefix('group_')],axis=1)
for s in X.columns:assert not any(t in s.lower() for t in ['cage','ctss','outcome','anchor','target','assay_eligibility','mapping_ambiguity'])
en=[dict(C=C,l1_ratio=l) for C in [.01,.1,1.,10.] for l in [0.,.25,.5,.75,1.]]
rf=[dict(n_estimators=300,max_depth=d,min_samples_leaf=l,max_features=m) for d in [6,12,None] for l in [10,50,200] for m in ['sqrt',.5]]
def make(name,params):
 return LogisticRegression(**params,solver='saga',penalty='elasticnet',max_iter=5000,random_state=SEED,n_jobs=1,class_weight=None) if name=='ELASTIC_NET' else RandomForestClassifier(**params,bootstrap=True,class_weight=None,random_state=SEED,n_jobs=1)
def fit_model(name,params,a,y0,b):
 with warnings.catch_warnings(record=True) as ww:
  warnings.simplefilter('always');m=make(name,params).fit(a,y0)
 return m,m.predict_proba(b)[:,list(m.classes_).index(1)],sum(issubclass(v.category,ConvergenceWarning) for v in ww)
def safe_train(mask,valid):return mask&~np.isin(comp,np.unique(comp[valid]))
OOF=pd.DataFrame(index=g);OOF['target']=y;OOF['S2_target']=label.loc[g,'S2_target'];OOF['inner_fold']=fold;OOF['dependency_component_id']=comp
for b in ['B0','B2','SQANTI3_RULE','PSQAN_STYLE']:OOF[b]=base[b]
OOF.loc[meta.SQANTI3_candidate_coverage<1,'SQANTI3_RULE']=np.nan
OOF.loc[meta.PSQAN_candidate_coverage<1,'PSQAN_STYLE']=np.nan
tuning=[];fits=[];foldaudit=[]
(O/'models').mkdir(exist_ok=True)
for f in range(5):
 va=fold==f;tr=safe_train(fold!=f,va)
 assert not set(comp[tr])&set(comp[va])
 assert tr.sum()>1000 and len(np.unique(y[tr]))==2
 foldaudit.append(dict(outer_fold=f,train=int(tr.sum()),validation=int(va.sum()),dependency_purged=int(((fold!=f)&~tr).sum()),train_positive_rate=float(y[tr].mean()),validation_positive_rate=float(y[va].mean())))
 q=[]
 for vals in [mi[['LR_read_R1','LR_read_R2']].mean(axis=1).to_numpy(),mi[['LR_junction_R1','LR_junction_R2']].mean(axis=1).to_numpy(),base.B2.to_numpy()]:
  q.append(weighted_ecdf(vals[tr],np.ones(tr.sum()),vals[va]))
 q=np.column_stack(q);lr=np.nanmean(q[:,:2],axis=1);OOF.loc[va,'SOURCE_AWARE']=np.mean(np.column_stack([lr,q[:,2]]),axis=1);OOF.loc[va,'READ_UNION']=q[:,0]
 for name,grid in [('ELASTIC_NET',en),('RANDOM_FOREST',rf)]:
  # Two grouped inner folds, deterministic on merged dependence component ID.
  inner=np.array([int(hnum(c+'|ISO001A_INNER|'+str(SEED))>=.5) for c in comp]);cache=[]
  for inf in [0,1]:
   iv=tr&(inner==inf);it=safe_train(tr&(inner!=inf),iv);pp=prep();cache.append((it,iv,pp.fit_transform(X.loc[it]),pp.transform(X.loc[iv])))
  scores=[]
  for k,params in enumerate(grid):
   ep=[];aw=0
   for inf,(it,iv,a,b) in enumerate(cache):
    st=time.monotonic();m,pv,nw=fit_model(name,params,a,y[it],b);r=metrics(y[iv],pv);aw+=nw;ep.append(r['enrichment']);tuning.append(dict(outer_fold=f,model=name,config_index=k,params=json.dumps(params,sort_keys=True),inner_fold=inf,convergence_warnings=nw,wall_seconds=time.monotonic()-st,**r))
   scores.append(np.mean(ep) if aw==0 else -np.inf)
  assert np.isfinite(scores).any(),'No converged baseline configuration'
  k=int(np.argmax(scores));params=grid[k];pp=prep();a=pp.fit_transform(X.loc[tr]);b=pp.transform(X.loc[va]);st=time.monotonic();m,pv,nw=fit_model(name,params,a,y[tr],b);assert nw==0,'Chosen baseline failed convergence'
  OOF.loc[va,name]=pv;fits.append(dict(fold=f,model=name,params=params,selected_inner_enrichment=float(scores[k]),wall_seconds=time.monotonic()-st,convergence_warnings=nw))
  joblib.dump(dict(model=m,preprocessor=pp,features=X.columns.tolist(),params=params),O/f'models/{name}_fold{f}.joblib',compress=3)
  print(json.dumps(dict(event='baseline_outer_completed',fold=f,model=name,params=params,OOF_metrics=metrics(y[va],pv)),default=str),flush=True)
 # Finite preregistered model/ablation set; no tuning against outer labels.
 for mode,name in [('proper','ISOCRED_MAP'),('replicated','ABL_CANDIDATE_REPLICATION'),('cumulative','ABL_SOURCE_CUMULATIVE'),('zero','ABL_MISSING_ZERO'),('no_annotation','ISOCRED_NO_ANNOTATION'),('caller_TALON','CALLER_TALON'),('caller_FLAIR','CALLER_FLAIR'),('caller_IsoQuant','CALLER_IsoQuant')]:
  z=minimal(mode);pp=prep(mode=='zero');a=pp.fit_transform(z.loc[tr]);b=pp.transform(z.loc[va]);names=pp.get_feature_names_out();src=np.array([0 if 'LR_' in s else 1 if 'annotation' in s else 2 for s in names]);weights=meta.candidate_count.to_numpy()[tr] if mode=='replicated' else None
  st=time.monotonic();m=HierarchicalMAP().fit(a,y[tr],comp[tr],weight=weights,source=src);pv=m.predict_proba(b,comp[va]);OOF.loc[va,name]=pv
  fits.append(dict(fold=f,model=name,params=dict(coefficient_within_source_sd=1,source_mean_sd=1,component_sd=.5,intercept_sd=5,quadrature_nodes=15),train_loss=m.loss,iterations=m.iterations,max_gradient=m.max_gradient,wall_seconds=time.monotonic()-st,one_likelihood_per_group=(mode!='replicated')))
  joblib.dump(dict(model=m,preprocessor=pp,features=z.columns.tolist()),O/f'models/{name}_fold{f}.joblib',compress=3)
 parq(f'results/OOF_checkpoint_fold{f}.parquet',OOF.reset_index());tsv(f'results/tuning_through_fold{f}.tsv',tuning);print('all outer models complete',f,flush=True)
parq('results/ISO-001A_development_OOF.parquet',OOF.reset_index());tsv('results/ISO-001A_grouped_inner_tuning.tsv',tuning);tsv('results/ISO-001A_outer_fold_audit.tsv',foldaudit);js('results/ISO-001A_model_parameters.json',fits)
receipt('real_comparison',start,evaluation_groups=len(g),development_groups_exported=len(label),OOF_models=len(OOF.columns)-4,outer_folds=5,inner_folds=2,baseline_grid=en+rf,contract_sha256=sha(contract),E4_visibility='DEVELOPMENT_ONLY',forbidden_E4_reads=0)
