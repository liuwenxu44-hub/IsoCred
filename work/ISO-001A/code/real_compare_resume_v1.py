"""Versioned I/O-only recovery of frozen real_compare.py after recorded SIGPIPE."""
from common import *
import difflib,joblib
src=W/'code/real_compare.py'
frozen=json.loads((O/'receipts/protocol_freeze.json').read_text())
original=src.read_text()
def replace_once(s,a,b):
 assert s.count(a)==1,(a,s.count(a))
 return s.replace(a,b)
text=replace_once(original,"tuning=[];fits=[];foldaudit=[]",'''tuning=pd.read_csv(O/'results/tuning_through_fold0.tsv',sep='\\t').to_dict('records');fits=[];foldaudit=[]
old_checkpoint=pd.read_parquet(O/'results/OOF_checkpoint_fold0.parquet').set_index('promoter_group_id')
assert old_checkpoint.index.equals(OOF.index)
OOF=old_checkpoint.copy()
for name in ['ELASTIC_NET','RANDOM_FOREST','ISOCRED_MAP','ABL_CANDIDATE_REPLICATION','ABL_SOURCE_CUMULATIVE','ABL_MISSING_ZERO','ISOCRED_NO_ANNOTATION','CALLER_TALON','CALLER_FLAIR','CALLER_IsoQuant']:
 saved=joblib.load(O/f'models/{name}_fold0.joblib');m=saved['model'];r=dict(fold=0,model=name,recovered_after_SIGPIPE=True,wall_seconds=None,wall_seconds_status='IN_MEMORY_TIMING_LOST_NOT_RECONSTRUCTED')
 if name in ['ELASTIC_NET','RANDOM_FOREST']:
  r['params']=saved['params'];rows=[t for t in tuning if t['model']==name and json.loads(t['params'])==saved['params']];assert len(rows)==2
  r['selected_inner_enrichment']=float(np.mean([t['enrichment'] for t in rows]));r['convergence_warnings']=sum(t['convergence_warnings'] for t in rows)
 else:r.update(train_loss=m.loss,iterations=m.iterations,max_gradient=m.max_gradient,one_likelihood_per_group=(name!='ABL_CANDIDATE_REPLICATION'),params=dict(coefficient_within_source_sd=1,source_mean_sd=1,component_sd=.5,intercept_sd=5,quadrature_nodes=15))
 fits.append(r)''')
text=replace_once(text," q=[]", " if f==0:continue\n q=[]")
text=replace_once(text,"  joblib.dump(dict(model=m,preprocessor=pp,features=X.columns.tolist(),params=params),O/f'models/{name}_fold{f}.joblib',compress=3)",'''  path=O/f'models/{name}_fold{f}.joblib'
  if path.exists():
   assert f==1 and name=='ELASTIC_NET'
   prior=joblib.load(path);assert prior['params']==params and prior['features']==list(X.columns)
   previous=prior['model'].predict_proba(prior['preprocessor'].transform(X.loc[va]))[:,list(prior['model'].classes_).index(1)]
   delta=float(np.max(abs(previous-pv)));assert delta<1e-12
   js('receipts/recovered_EN_fold1_identity.json',dict(params=params,max_prediction_difference=delta,old_model_sha256=sha(path),old_model_preserved=True,inner_grid_repeated_only_to_recover_lost_records=True))
  else:joblib.dump(dict(model=m,preprocessor=pp,features=X.columns.tolist(),params=params),path,compress=3)
  tsv(f'results/recovery_tuning_fold{f}_{name}.tsv',tuning)
  js(f'results/recovery_fits_fold{f}_{name}.json',fits)''')
assert not (O/'receipts/real_compare_resume_v1_transform.json').exists()
js('receipts/real_compare_resume_v1_transform.json',dict(original_code_sha256=sha(src),resume_code_sha256=sha(Path(__file__)),transformed_code_sha256=hashlib.sha256(text.encode()).hexdigest(),diff=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),fromfile='frozen_real_compare.py',tofile='runtime_recovery_v1')),scientific_contract_unchanged=True,seed_unchanged=True,old_fold0_reused=True))
exec(compile(text,str(src)+'[runtime_recovery_v1]','exec'),globals())
