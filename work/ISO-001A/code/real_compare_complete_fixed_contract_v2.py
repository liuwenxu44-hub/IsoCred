"""Complete fixed comparators, preserving EN nonconvergence as unavailable."""
from common import *
import difflib
src=W/'code/real_compare.py';original=src.read_text();fr=json.loads((O/'receipts/protocol_freeze.json').read_text());assert sha(src)==fr['science_code']['real_compare.py']
def sub(s,a,b):
 assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)
text=sub(original,"tuning=[];fits=[];foldaudit=[]",'''tuning=pd.read_csv(O/'results/tuning_through_fold2.tsv',sep='\\t').to_dict('records');fits=json.loads((O/'results/recovery_fits_fold2_RANDOM_FOREST.json').read_text());foldaudit=[]
old_checkpoint=pd.read_parquet(O/'results/OOF_checkpoint_fold2.parquet').set_index('promoter_group_id');assert old_checkpoint.index.equals(OOF.index);OOF=old_checkpoint.copy()
for name in ['ISOCRED_MAP','ABL_CANDIDATE_REPLICATION','ABL_SOURCE_CUMULATIVE','ABL_MISSING_ZERO','ISOCRED_NO_ANNOTATION','CALLER_TALON','CALLER_FLAIR','CALLER_IsoQuant']:
 m=joblib.load(O/f'models/{name}_fold2.joblib')['model'];assert not any(r['fold']==2 and r['model']==name for r in fits)
 fits.append(dict(fold=2,model=name,params=dict(coefficient_within_source_sd=1,source_mean_sd=1,component_sd=.5,intercept_sd=5,quadrature_nodes=15),train_loss=m.loss,iterations=m.iterations,max_gradient=m.max_gradient,one_likelihood_per_group=name!='ABL_CANDIDATE_REPLICATION',wall_seconds=None,wall_seconds_status='LOST_IN_MEMORY_TIMING_AFTER_FAILURE_NOT_INVENTED'))
unavailable=[]''')
text=sub(text," q=[]"," if f<3:continue\n q=[]")
needle=";assert nw==0,'Chosen baseline failed convergence'"
text=sub(text,needle,'''\n  tsv(f'results/completion_v2_tuning_fold{f}_{name}.tsv',tuning)
  if nw:
   assert name=='ELASTIC_NET'
   repeat,again,nw2=fit_model(name,params,a,y[tr],b)
   issue=dict(fold=f,model=name,params=params,selected_inner_enrichment=float(scores[k]),status='UNAVAILABLE_OUTER_NONCONVERGENCE',convergence_warnings=nw,repeat_convergence_warnings=nw2,first_n_iter=m.n_iter_.tolist(),repeat_n_iter=repeat.n_iter_.tolist(),repeat_prediction_difference=float(np.max(abs(pv-again))),max_iter_unchanged=5000,seed_unchanged=SEED,outer_validation_predictions_not_accepted=True,alternative_configuration_not_selected=True,wall_seconds=time.monotonic()-st,one_likelihood_per_group=True)
   js(f'receipts/EN_outer_nonconvergence_fold{f}.json',issue)
   joblib.dump(dict(model=m,preprocessor=pp,features=X.columns.tolist(),params=params,acceptance='FAILED_NOT_SCORED'),O/f'models/FAILED_NOT_SCORED_{name}_fold{f}.joblib',compress=3)
   OOF.loc[va,name]=np.nan;fits.append(issue);unavailable.append(issue)
   print(json.dumps(dict(event='FROZEN_MAX_ITER_NONCONVERGENCE_REPRODUCED',**issue)),flush=True)
   js(f'results/completion_v2_fits_fold{f}_{name}.json',fits)
   continue''')
text=sub(text,"  joblib.dump(dict(model=m,preprocessor=pp,features=X.columns.tolist(),params=params),O/f'models/{name}_fold{f}.joblib',compress=3)","  assert not (O/f'models/{name}_fold{f}.joblib').exists()\n  joblib.dump(dict(model=m,preprocessor=pp,features=X.columns.tolist(),params=params),O/f'models/{name}_fold{f}.joblib',compress=3)\n  js(f'results/completion_v2_fits_fold{f}_{name}.json',fits)")
text=sub(text,"receipt('real_comparison',start,", "js('receipts/mandatory_comparator_nonconvergence.json',dict(unavailable=unavailable,complete_required_comparison=not bool(unavailable),max_iter_unchanged=5000,GO_prohibited_if_incomplete=True))\nreceipt('real_comparison',start,")
js('receipts/completion_v2_transform.json',dict(original_code_sha256=sha(src),adapter_sha256=sha(Path(__file__)),diff=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),fromfile='original_frozen_driver',tofile='failure_handling_v2')),scientific_unit_metric_seed_grid_unchanged=True,allowed_output='complete all remaining registered models; nonconverged EN gets no OOF score and cannot support GO'))
exec(compile(text,str(src)+'[fixed_contract_completion_v2]','exec'),globals())
