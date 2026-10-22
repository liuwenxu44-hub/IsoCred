"""Finite grouped OOF supplement. Only development targets are readable."""
from common import *
from metrics import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
import warnings,joblib,sys
contract=validate_contract();started=time.monotonic();d=load_development();meta=d['metadata'];mi=d['minimal']
name=sys.argv[1];assert name in contract['inputs']
X=mi[contract['inputs'][name]];y=d['labels'].primary_target.to_numpy(dtype=int)
assert set(y)=={0,1}
comp=meta.dependency_component_id.to_numpy();fold=meta.inner_fold.to_numpy()
inner=np.array([int(hnum(c+'|ISO001A_INNER|'+str(SEED))>=.5) for c in comp])
spec=contract['models'][name if name in ['RF14','EN14'] else 'L2_controls'];grid=spec['grid']
def fit(params,a,y0):
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always')
  if name=='RF14':m=RandomForestClassifier(**params,bootstrap=True,class_weight=None,random_state=SEED,n_jobs=1)
  else:m=LogisticRegression(**params,solver=spec['solver'],penalty=spec['penalty'],max_iter=spec['max_iter'],tol=spec['tol'],random_state=SEED,n_jobs=1,class_weight=None)
  m.fit(a,y0)
 return m,[str(w.message) for w in caught if issubclass(w.category,ConvergenceWarning)]
def predict(m,a):
 assert list(m.classes_)==[0,1]
 return m.predict_proba(a)[:,list(m.classes_).index(1)]
OOF=pd.DataFrame({'promoter_group_id':X.index,name:np.nan})
tuning=[];fits=[]
for f in range(5):
 va=fold==f;tr=(fold!=f)&~np.isin(comp,np.unique(comp[va]))
 assert not set(comp[tr]).intersection(comp[va])
 cache=[]
 for inf in [0,1]:
  it=tr&(inner!=inf);iv=tr&(inner==inf)
  assert not set(comp[it]).intersection(comp[iv])
  pp=prep();a=pp.fit_transform(X.loc[it]);b=pp.transform(X.loc[iv]);cache.append((it,iv,a,b))
 scores=[]
 for k,params in enumerate(grid):
  es=[];nw=0
  for inf,(it,iv,a,b) in enumerate(cache):
   st=time.monotonic();m,ww=fit(params,a,y[it]);pred=predict(m,b);r=metrics(y[iv],pred)
   es.append(r['enrichment']);nw+=len(ww)
   tuning.append(dict(model=name,outer_fold=f,inner_fold=inf,config_index=k,params=json.dumps(params,sort_keys=True),convergence_warnings=len(ww),warning_messages=' | '.join(ww),train_n=int(it.sum()),validation_n=int(iv.sum()),wall_seconds=time.monotonic()-st,**r))
  scores.append(float(np.mean(es)) if nw==0 else -np.inf)
 valid=np.isfinite(scores)
 if not valid.any():
  fits.append(dict(model=name,outer_fold=f,status='NO_CONVERGED_INNER_CONFIG_NO_SCORE',features=X.columns.tolist()))
 else:
  k=int(np.argmax(scores));params=grid[k]
  js(O/'receipts'/f'{name}_fold{f}_selection.json',dict(config_index=k,params=params,inner_enrichment=scores[k],outer_refit_attempts_allowed=1))
  pp=prep();a=pp.fit_transform(X.loc[tr]);b=pp.transform(X.loc[va])
  if name in ['RF14','EN14']:
   inherited=joblib.load(R/f'work/ISO-001A/output/models/ISOCRED_MAP_fold{f}.joblib')
   assert inherited['features']==X.columns.tolist()
   np.testing.assert_array_equal(b,inherited['preprocessor'].transform(X.loc[va]))
  st=time.monotonic();m,ww=fit(params,a,y[tr])
  payload=dict(model=m,preprocessor=pp,features=X.columns.tolist(),params=params,training_group_ids=X.index[tr].to_numpy(),training_component_ids=np.unique(comp[tr]),validation_group_ids=X.index[va].to_numpy(),contract_sha256=sha(R/'contracts/ISO-002C_phase1_comparison_contract.yaml'),status='ACCEPTED' if not ww else 'FAILED_NONCONVERGENCE_NO_SCORE',evidence_status=EVIDENCE)
  dest=O/'models'/f'{name}_fold{f}.joblib';assert not dest.exists();joblib.dump(payload,dest,compress=3)
  if not ww:
   pv=predict(m,b);assert np.isfinite(pv).all();OOF.loc[va,name]=pv
  fits.append(dict(model=name,outer_fold=f,status=payload['status'],params=params,warning_messages=ww,train_n=int(tr.sum()),validation_n=int(va.sum()),selected_inner_enrichment=scores[k],transformed_columns=pp.get_feature_names_out().tolist(),features=X.columns.tolist(),wall_seconds=time.monotonic()-st))
 parq(O/'results'/f'{name}_OOF_through_fold{f}.parquet',OOF)
 print(json.dumps(dict(event='outer_fold_completed',model=name,fold=f,status=fits[-1]['status']),default=str),flush=True)
parq(O/'results'/f'{name}_OOF.parquet',OOF)
pd.DataFrame(tuning).to_csv(O/'results'/f'{name}_inner_tuning.tsv',sep='\t',index=False)
js(O/'results'/f'{name}_fit_parameters.json',fits)
js(O/'receipts'/f'{name}_training_complete.json',dict(status='COMPLETE_WITH_FAILURES_RETAINED' if OOF[name].isna().any() else 'COMPLETE',wall_seconds=time.monotonic()-started,cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,inner_fits=len(tuning),outer_refit_attempts=sum('params' in v for v in fits),scored=int(OOF[name].notna().sum()),eligible=len(OOF),new_final_predictions=0,new_external_outcome_reads=0,contract_sha256=sha(R/'contracts/ISO-002C_phase1_comparison_contract.yaml')))
