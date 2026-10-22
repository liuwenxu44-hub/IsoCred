"""Fresh model inference rebuild, without fitting or final labels."""
from final_inference import *
import joblib,sys
start=time.monotonic();name=sys.argv[1];assert name in ['A','B']
destination=O/'results'/('final_predictions_'+name)
assert not destination.exists();destination.mkdir()
d=load_features();meta=d['metadata'];base=d['baselines'];ev=d['evidence'];core=d['core']
test=meta.index[meta.split=='FINAL_BLIND_TEST'];mi=minimal(ev,base,meta);X=pd.concat([core,mi.add_prefix('group_')],axis=1)
out=meta.loc[test,['split','dependency_component_id']].copy()
mask=pd.read_parquet(O/'data/final_eligibility_label_blind.parquet').set_index('promoter_group_id');out['evaluation_eligible']=mask.loc[test,'evaluation_eligible']
for n in ['B0','B2','SQANTI3_RULE','PSQAN_STYLE']:out[n]=base.loc[test,n]
out.loc[meta.loc[test,'SQANTI3_candidate_coverage']<1,'SQANTI3_RULE']=np.nan
out.loc[meta.loc[test,'PSQAN_candidate_coverage']<1,'PSQAN_STYLE']=np.nan
ref=pd.read_parquet(O/'models/source_aware_training_reference.parquet').set_index('promoter_group_id')
q=[]
for name_col,v in [('read',mi.loc[test,['LR_read_R1','LR_read_R2']].mean(axis=1)),('junction',mi.loc[test,['LR_junction_R1','LR_junction_R2']].mean(axis=1)),('annotation',base.loc[test,'B2'])]:
 q.append(weighted_ecdf(ref[name_col].to_numpy(),np.ones(len(ref)),v.to_numpy()))
q=np.column_stack(q)
with warnings.catch_warnings():
 warnings.simplefilter('ignore',RuntimeWarning);lr=np.nanmean(q[:,:2],axis=1)
out['SOURCE_AWARE']=np.mean(np.column_stack([lr,q[:,2]]),axis=1);out['READ_UNION']=q[:,0]
params=json.loads((O/'results/final_model_parameters.json').read_text());status={x['model']:x['status'] for x in params}
for n in ['ELASTIC_NET','RANDOM_FOREST','ISOCRED_MAP']:
 if status[n]!='ACCEPTED':
  assert n=='ELASTIC_NET' and status[n] in ['FAILED_FINAL_NONCONVERGENCE','NO_CONVERGED_INNER_CONFIGURATION']
  out[n]=np.nan;continue
 p=joblib.load(O/f'models/{n}_final.joblib')
 assert not set(p['training_component_ids']).intersection(out.dependency_component_id)
 z=mi if n=='ISOCRED_MAP' else X
 b=p['preprocessor'].transform(z.loc[test,p['features']])
 out[n]=p['model'].predict_proba(b,out.dependency_component_id.to_numpy()) if n=='ISOCRED_MAP' else prob(p['model'],b)
for n in MAIN+OPTIONAL:out[n+'_scorable']=np.isfinite(out[n])
out['primary_common_mask']=out.evaluation_eligible & out[[n+'_scorable' for n in MAIN]].all(axis=1)
out['evidence_status']='FINAL_PREDICTIONS_FROZEN_BEFORE_LABEL_UNLOCK'
p=destination/'predictions.parquet';pq.write_table(pa.Table.from_pandas(out.reset_index(),preserve_index=False),p,compression='zstd')
old=pd.read_parquet(O/'final_predictions/ISO-002B_final_predictions.parquet').set_index('promoter_group_id')
pd.testing.assert_frame_equal(out,old,check_exact=True)
result=dict(status='PASS_EXACT_INFERENCE_REBUILD',root=name,rows=len(out),sha256=sha(p),original_sha256=sha(O/'final_predictions/ISO-002B_final_predictions.parquet'),no_fitting=True,final_label_reads=0)
(destination/'rebuild_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
receipt('prediction_rebuild_'+name,start,**result)
