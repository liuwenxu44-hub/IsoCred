"""Only authorized common10 development fitting; requires completed P1 and frozen protocol."""
from common import *
from metrics import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
import warnings,joblib,sys,time
def fit_model(name,spec,params,a,y):
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always')
  if name=='RF10':
   m=RandomForestClassifier(**params,bootstrap=True,class_weight=None,random_state=SEED,n_jobs=1)
  else:
   m=LogisticRegression(**params,solver=spec['solver'],penalty=spec['penalty'],max_iter=spec['max_iter'],tol=spec['tol'],random_state=SEED,n_jobs=1,class_weight=None)
  m.fit(a,y)
 return m,[str(w.message) for w in caught if issubclass(w.category,ConvergenceWarning)],[str(w.message) for w in caught]
def positive(m,a):
 assert list(m.classes_)==[0,1]
 p=m.predict_proba(a)[:,list(m.classes_).index(1)]
 assert np.isfinite(p).all() and ((p>=0)&(p<=1)).all()
 return p
def main(name):
 assert name in ['RF10','EN10']
 contract=validate_contract();started=time.monotonic();d=load_development()
 X=d['features'];meta=d['metadata'];y=d['labels'].primary_target.to_numpy(dtype=int)
 assert set(y)=={0,1};assert X.columns.tolist()==contract['features']
 comp=meta.dependency_component_id.to_numpy();fold=meta.inner_fold.to_numpy()
 assert set(fold)==set(range(5))
 inner=np.array([int(hnum(c+'|ISO001A_INNER|'+str(SEED))>=.5) for c in comp])
 spec=contract['models'][name];grid=spec['grid']
 assert len(grid)==(18 if name=='RF10' else 20)
 OOF=pd.DataFrame({'promoter_group_id':X.index,name:np.nan})
 tuning=[];fits=[]
 for f in range(5):
  va=fold==f;tr=(fold!=f)&~np.isin(comp,np.unique(comp[va]))
  assert not set(comp[tr]).intersection(comp[va])
  assert not set(X.index[tr]).intersection(X.index[va])
  cache=[]
  for inf in [0,1]:
   it=tr&(inner!=inf);iv=tr&(inner==inf)
   assert not set(comp[it]).intersection(comp[iv])
   assert set(y[it])=={0,1} and set(y[iv])=={0,1}
   pp=prep();a=pp.fit_transform(X.loc[it]);b=pp.transform(X.loc[iv])
   cache.append((it,iv,a,b))
  scores=[]
  for k,params in enumerate(grid):
   es=[];nw=0
   for inf,(it,iv,a,b) in enumerate(cache):
    st=time.monotonic();m,ww,allwarnings=fit_model(name,spec,params,a,y[it]);pred=positive(m,b)
    r=metrics(y[iv],pred,probability=True);es.append(r['enrichment']);nw+=len(ww)
    row=dict(model=name,outer_fold=f,inner_fold=inf,config_index=k,params=json.dumps(params,sort_keys=True),
     convergence_warnings=len(ww),warning_messages=' | '.join(allwarnings),
     train_n=int(it.sum()),validation_n=int(iv.sum()),train_prevalence=float(y[it].mean()),
     validation_prevalence=float(y[iv].mean()),wall_seconds=time.monotonic()-st,**r)
    tuning.append(row)
    js(W/'training'/f'{name}_outer{f}_config{k}_inner{inf}.json',row)
   scores.append(float(np.mean(es)) if nw==0 else -np.inf)
  if not np.isfinite(scores).any():
   fits.append(dict(model=name,outer_fold=f,status='NO_CONVERGED_INNER_CONFIG_NO_SCORE',features=FEATURES))
  else:
   k=int(np.argmax(scores));params=grid[k]
   js(W/'training'/f'{name}_fold{f}_selection.json',dict(config_index=k,params=params,inner_enrichment=scores[k],outer_refit_attempts_allowed=1))
   pp=prep();a=pp.fit_transform(X.loc[tr]);b=pp.transform(X.loc[va])
   st=time.monotonic();m,ww,allwarnings=fit_model(name,spec,params,a,y[tr])
   payload=dict(model=m,preprocessor=pp,features=FEATURES,params=params,model_family=name,
    training_group_ids=X.index[tr].to_numpy(),training_component_ids=np.unique(comp[tr]),
    validation_group_ids=X.index[va].to_numpy(),outer_fold=f,seed=SEED,
    contract_sha256=sha(R/'contracts/ISO-002D_common10_computation_contract.json'),
    status='ACCEPTED' if not ww else 'FAILED_NONCONVERGENCE_NO_SCORE',evidence_status=EVIDENCE)
   dest=O/'models'/f'{name}_fold{f}.joblib';assert not dest.exists();joblib.dump(payload,dest,compress=3)
   if not ww:OOF.loc[va,name]=positive(m,b)
   fits.append(dict(model=name,outer_fold=f,status=payload['status'],params=params,
    warning_messages=allwarnings,train_n=int(tr.sum()),validation_n=int(va.sum()),
    selected_inner_enrichment=scores[k],transformed_columns=pp.get_feature_names_out().tolist(),
    features=FEATURES,wall_seconds=time.monotonic()-st))
  parq(O/'results'/f'{name}_OOF_through_fold{f}.parquet',OOF)
  print(json.dumps(dict(event='outer_fold_completed',model=name,fold=f,status=fits[-1]['status'])),flush=True)
 parq(O/'results'/f'{name}_OOF.parquet',OOF)
 dest=O/'results'/f'{name}_inner_tuning.tsv';assert not dest.exists()
 pd.DataFrame(tuning).to_csv(dest,sep='\t',index=False)
 js(O/'results'/f'{name}_fit_parameters.json',fits)
 js(W/'receipts'/f'{name}_training_complete.json',dict(
  status='COMPLETE_WITH_FAILURES_RETAINED' if OOF[name].isna().any() else 'COMPLETE',
  wall_seconds=time.monotonic()-started,
  cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
  peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,inner_fits=len(tuning),
  outer_refit_attempts=sum('params' in v for v in fits),scored=int(OOF[name].notna().sum()),
  eligible=len(OOF),new_final_predictions=0,new_external_outcome_reads=0,
  contract_sha256=sha(R/'contracts/ISO-002D_common10_computation_contract.json')))
if __name__=='__main__':main(sys.argv[1])
