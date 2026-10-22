"""One predeclared development-only final selection/refit. Never reads final E4."""
from final_inference import *
from metrics import metrics
import joblib
start=time.monotonic()
protocol=R/'contracts/ISO-002B_final_evaluation_protocol.yaml'
fr=json.loads((R/'logs/ISO-002B_protocol_freeze.json').read_text())
assert sha(protocol)==fr['contract_sha256']
for n,h in fr['code_sha256'].items():assert sha(W/'code'/n)==h,n
assert not (O/'final_predictions/ISO-002B_final_predictions.parquet').exists()
d=load_features();meta=d['metadata'];base=d['baselines'];ev=d['evidence'];core=d['core']
labels=pd.read_parquet(R/'work/ISO-001A/output/data/development_E4_v1.parquet').set_index('promoter_group_id')
train=labels.index[labels.evaluation_eligible]
test=meta.index[meta.split=='FINAL_BLIND_TEST']
assert len(train)==30487 and len(test)==11602
assert (meta.loc[train,'split']=='DEVELOPMENT_POOL').all()
assert not meta.loc[train,'touches_final_dependency'].any()
comp=meta.loc[train,'dependency_component_id'].to_numpy()
test_comp=meta.loc[test,'dependency_component_id'].to_numpy()
assert not set(comp).intersection(test_comp)
mi=minimal(ev,base,meta);X=pd.concat([core,mi.add_prefix('group_')],axis=1)
for name in X.columns:assert not any(s in name.lower() for s in ['cage','ctss','target','outcome','anchor','assay_eligibility','mapping_ambiguity'])
y=labels.loc[train,'primary_target'].to_numpy(dtype=int)
assert set(y)=={0,1}
out=meta.loc[test,['split','dependency_component_id']].copy()
mask=pd.read_parquet(O/'data/final_eligibility_label_blind.parquet').set_index('promoter_group_id')
out['evaluation_eligible']=mask.loc[test,'evaluation_eligible']
for name in ['B0','B2','SQANTI3_RULE','PSQAN_STYLE']:out[name]=base.loc[test,name]
out.loc[meta.loc[test,'SQANTI3_candidate_coverage']<1,'SQANTI3_RULE']=np.nan
out.loc[meta.loc[test,'PSQAN_candidate_coverage']<1,'PSQAN_STYLE']=np.nan
for n,v in source_scores(mi.loc[train],base.loc[train],mi.loc[test],base.loc[test]).items():out[n]=v
reference=pd.DataFrame({'read':mi.loc[train,['LR_read_R1','LR_read_R2']].mean(axis=1),'junction':mi.loc[train,['LR_junction_R1','LR_junction_R2']].mean(axis=1),'annotation':base.loc[train,'B2']})
parq('models/source_aware_training_reference.parquet',reference.reset_index())
inner=np.array([int(hnum(c+'|ISO001A_INNER|'+str(SEED))>=.5) for c in comp])
cache=[]
for inf in [0,1]:
 it=inner!=inf;iv=inner==inf
 assert not set(comp[it]).intersection(comp[iv])
 pp=prep()
 cache.append((it,iv,pp.fit_transform(X.loc[train[it]]),pp.transform(X.loc[train[iv]])))
tuning=[];parameters=[]
for name,grid in [('ELASTIC_NET',EN_GRID),('RANDOM_FOREST',RF_GRID)]:
 scores=[]
 for k,params in enumerate(grid):
  ep=[];warning_count=0
  for inf,(it,iv,a,b) in enumerate(cache):
   begin=time.monotonic();model,warnings_count=fit_model(name,params,a,y[it]);pred=prob(model,b)
   result=metrics(y[iv],pred);ep.append(result['enrichment']);warning_count+=warnings_count
   tuning.append(dict(model=name,config_index=k,params=json.dumps(params,sort_keys=True),inner_fold=inf,train_groups=int(it.sum()),validation_groups=int(iv.sum()),convergence_warnings=warnings_count,wall_seconds=time.monotonic()-begin,**result))
  scores.append(float(np.mean(ep)) if warning_count==0 else -np.inf)
  tsv('results/final_training_inner_tuning.tsv',tuning)
  print(json.dumps({'event':'fixed_inner_configuration_attempted','model':name,'config_index':k,'inner_enrichment':float(np.mean(ep)),'converged':warning_count==0}),flush=True)
 if not np.isfinite(scores).any():
  assert name=='ELASTIC_NET','No valid RF configuration'
  record=dict(model=name,status='NO_CONVERGED_INNER_CONFIGURATION',scoring='UNAVAILABLE_NO_IMPUTATION')
  out[name]=np.nan;parameters.append(record);js('receipts/EN_final_training_unavailable.json',record);continue
 k=int(np.argmax(scores));params=grid[k]
 selected=dict(model=name,config_index=k,params=params,selected_inner_enrichment=scores[k],final_attempt_limit=1)
 js('receipts/'+name+'_selection_before_final_fit.json',selected)
 pp=prep();a=pp.fit_transform(X.loc[train]);b=pp.transform(X.loc[test])
 begin=time.monotonic();model,nw=fit_model(name,params,a,y)
 record=dict(**selected,convergence_warnings=nw,wall_seconds=time.monotonic()-begin,train_groups=len(train),features=X.columns.tolist())
 payload=dict(model=model,preprocessor=pp,features=X.columns.tolist(),params=params,training_component_ids=np.unique(comp),training_group_ids=train.to_numpy())
 if nw:
  assert name=='ELASTIC_NET','RF refit cannot be accepted'
  record.update(status='FAILED_FINAL_NONCONVERGENCE',scoring='UNAVAILABLE_NO_IMPUTATION',max_iter=5000)
  joblib.dump(payload,O/'models/FAILED_NOT_SCORED_ELASTIC_NET_final.joblib',compress=3)
  out[name]=np.nan;js('receipts/EN_final_training_unavailable.json',record)
 else:
  record.update(status='ACCEPTED')
  out[name]=prob(model,b);assert np.isfinite(out[name]).all()
  joblib.dump(payload,O/f'models/{name}_final.joblib',compress=3)
 parameters.append(record);js('results/final_model_parameters.json',parameters)
 print(json.dumps({'event':'one_final_refit_attempt','model':name,'status':record['status'],'params':params}),flush=True)
pp=prep();a=pp.fit_transform(mi.loc[train]);b=pp.transform(mi.loc[test])
names=pp.get_feature_names_out();source=np.array([0 if 'LR_' in n else 1 if 'annotation' in n else 2 for n in names])
begin=time.monotonic();model=HierarchicalMAP().fit(a,y,comp,source=source)
out['ISOCRED_MAP']=model.predict_proba(b,test_comp)
joblib.dump(dict(model=model,preprocessor=pp,features=mi.columns.tolist(),training_component_ids=np.unique(comp),training_group_ids=train.to_numpy()),O/'models/ISOCRED_MAP_final.joblib',compress=3)
parameters.append(dict(model='ISOCRED_MAP',status='ACCEPTED',fit='ONE_GROUP_BERNOULLI_MAP',train_groups=len(train),iterations=model.iterations,max_gradient=model.max_gradient,loss=model.loss,wall_seconds=time.monotonic()-begin,posterior_calibration_claim=False))
js('results/final_model_parameters.json',parameters)
for n in MAIN+OPTIONAL:out[n+'_scorable']=np.isfinite(out[n])
out['primary_common_mask']=out.evaluation_eligible & out[[n+'_scorable' for n in MAIN]].all(axis=1)
assert out.loc[out.primary_common_mask,MAIN].notna().all().all()
out['evidence_status']='FINAL_PREDICTIONS_FROZEN_BEFORE_LABEL_UNLOCK'
p=parq('final_predictions/ISO-002B_final_predictions.parquet',out.reset_index())
tsv('final_predictions/ISO-002B_final_prediction_coverage.tsv',[dict(method=n,membership=len(out),eligible=int(out.evaluation_eligible.sum()),scorable_eligible=int((out.evaluation_eligible & out[n+'_scorable']).sum()),primary_common=int(out.primary_common_mask.sum())) for n in MAIN+OPTIONAL])
receipt('final_training',start,train_groups=len(train),final_predictions=len(out),primary_common=int(out.primary_common_mask.sum()),final_label_reads=0,protocol_sha256=sha(protocol),prediction_sha256=sha(p),inner_fits=len(tuning),final_refit_attempts=sum(x.get('status')!='NO_CONVERGED_INNER_CONFIGURATION' for x in parameters),no_post_final_selection=True)
