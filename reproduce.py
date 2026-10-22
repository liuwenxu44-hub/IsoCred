#!/usr/bin/env python3
"""Processed-table verification and fit-disabled checkpoint replay; no raw-data I/O."""
import argparse, gzip, hashlib, importlib.util, json, os, sys
from pathlib import Path
sys.dont_write_bytecode=True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
 os.environ[key]='1'
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
threadpool_limits(1)
ROOT=Path(__file__).resolve().parent
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def parquet(path): return pd.read_parquet(ROOT/path).set_index('promoter_group_id')
def same(a,b): np.testing.assert_allclose(a,b,rtol=0,atol=1e-12,equal_nan=True)
def load_py(name,path):
 spec=importlib.util.spec_from_file_location(name,ROOT/path);module=importlib.util.module_from_spec(spec)
 sys.modules[name]=module;spec.loader.exec_module(module);return module
engine=load_py('public_common10_engine','work/ISO-002D/evaluation_engine/engine.py')
def metric(y,s): return engine.weighted_metrics(np.asarray(y,dtype=float),np.asarray(s,dtype=float))
def verify():
 provenance=json.loads((ROOT/'provenance/source_to_release.json').read_text())
 for row in provenance['sources']:
  assert sha(ROOT/row['path'])==row['release_sha256'],row['path']
 manifest_path=ROOT/'provenance/release_manifest.json'
 manifest_count=None
 if manifest_path.exists():
  manifest=json.loads(manifest_path.read_text())
  for row in manifest['files']:
   path=ROOT/row['path'];assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],row['path']
  manifest_count=len(manifest['files'])
 checks=[]
 # Final heldout points: no bootstrap or fitted-model loading.
 final=parquet('output/ISO-002B/final_predictions/ISO-002B_final_predictions.parquet')
 labels=parquet('output/ISO-002B/final_metrics/final_E4_eligible_targets.parquet')
 ref=pd.read_csv(ROOT/'output/ISO-002B/final_metrics/evaluation/metrics.tsv',sep='\t')
 for _,row in ref.iterrows():
  mask=final.evaluation_eligible.astype(bool)&final[row.method].notna()
  if row.analysis_set=='PRIMARY_COMMON':mask&=final.primary_common_mask
  result=metric(labels.loc[mask,'primary_target'],final.loc[mask,row.method])
  same(result['enrichment@0.1'],row.enrichment);same(result['AUC'],row.AUC)
  same(result['actual_selection_coverage@0.1'],row.top_actual_coverage)
  checks.append(dict(case='HepG2_final',method=row.method,scope=row.analysis_set,n=int(mask.sum()),point_metrics='MATCH_1e-12'))
 # Post-holdout development points and failed EN71 fold remain missing.
 old=parquet('work/ISO-001A/output/results/ISO-001A_development_OOF.parquet')
 aliases={'RF71':'RANDOM_FOREST','EN71':'ELASTIC_NET','MAP14':'ISOCRED_MAP'}
 table=pd.read_csv(ROOT/'work/ISO-002C/evaluation/results_A/model_metrics.tsv',sep='\t')
 for _,row in table[table.cohort=='ELIGIBLE_DEVELOPMENT'].iterrows():
  name=row.model
  scores=old[aliases.get(name,name)] if name in aliases or name in old else parquet(f'output/ISO-002C/results/{name}_OOF.parquet')[name].reindex(old.index)
  result=metric(old.target,scores)
  same(result['enrichment@0.1'],row.full_ties_10_enrichment);same(result['AUC'],row.AUC)
  assert int(np.isfinite(scores).sum())==int(row.scored_groups)
  checks.append(dict(case='HepG2_post_holdout_development_only',method=name,n=int(np.isfinite(scores).sum()),point_metrics='MATCH_1e-12'))
 assert old.ELASTIC_NET.isna().sum()==6264 and old.loc[old.inner_fold==3,'ELASTIC_NET'].isna().all()
 # External group endpoint already unlocked in the historical analysis.
 pred=parquet('output/ISO-002D/predictions/WTC11_all_five_methods.parquet')
 lab=pd.read_csv(ROOT/'work/ISO-002D/endpoint_extractor/runs/WTC11_once_01/promoter_group_endpoints.tsv.gz',sep='\t').set_index('promoter_group_id').reindex(pred.index)
 eligible=lab.assay_eligible.astype(str).str.lower().isin(('true','1'))
 ledger=pd.read_csv(ROOT/'work/ISO-002D/evaluation_adapter/run/external/scientific_numbers.tsv.gz',sep='\t',low_memory=False)
 subset=ledger[(ledger.kind=='native')&(ledger.scope=='pooled')&(ledger.endpoint=='primary')]
 assert int(eligible.sum())==34150
 for name in engine.METHODS:
  result=metric(lab.loc[eligible,'primary_label'],pred.loc[eligible,name])
  for key in ('AUC','enrichment@0.1','actual_selection_coverage@0.1'):
   rows=subset[(subset.method==name)&(subset.metric==key)]
   assert len(rows)==1,(name,key,len(rows));same(result[key],rows.estimate_full_precision.iloc[0])
  checks.append(dict(case='WTC11_common10',method=name,n=result['active_unique_n'],point_metrics='MATCH_1e-12'))
 rf=metric(lab.loc[eligible,'primary_label'],pred.loc[eligible,'RF10'])['enrichment@0.1']
 en=metric(lab.loc[eligible,'primary_label'],pred.loc[eligible,'EN10'])['enrichment@0.1']
 assert rf-en<0
 return dict(status='PASS_PROCESSED_POINT_VERIFICATION',hash_bound_source_files=len(provenance['sources']),public_manifest_files_verified=manifest_count,checks=checks,EN71_failed_fold3_retained=True,WTC11_RF_minus_EN_enrichment=rf-en,bootstrap_rerun=False,model_fits=0,raw_data_read=False)
def replay():
 import joblib
 from sklearn.linear_model import LogisticRegression
 from sklearn.ensemble import RandomForestClassifier
 from sklearn.impute import SimpleImputer
 from sklearn.preprocessing import RobustScaler
 def no_fit(*args,**kwargs):raise RuntimeError('Model/preprocessing fitting disabled during public replay')
 for cls in (LogisticRegression,RandomForestClassifier,SimpleImputer,RobustScaler):cls.fit=no_fit
 # Use the original mathematical/helper code, with only its root made relative.
 os.chdir(ROOT);sys.path.insert(0,str(ROOT/'work/ISO-002B/code'))
 import final_inference as fi
 fi.HierarchicalMAP.fit=no_fit
 data=fi.load_features();meta=data['metadata'];base=data['baselines']
 ids=meta.index[meta.split=='FINAL_BLIND_TEST']
 mi=fi.minimal(data['evidence'],base,meta);full=pd.concat([data['core'],mi.add_prefix('group_')],axis=1)
 reference=parquet('output/ISO-002B/final_predictions/ISO-002B_final_predictions.parquet')
 checks=[]
 final_status={p['model']:p['status'] for p in json.loads((ROOT/'output/ISO-002B/results/final_model_parameters.json').read_text())}
 for name in ('ELASTIC_NET','RANDOM_FOREST','ISOCRED_MAP'):
  payload=joblib.load(ROOT/f'output/ISO-002B/models/{name}_final.joblib')
  assert final_status[name]=='ACCEPTED'
  assert not set(payload['training_component_ids']).intersection(meta.loc[ids,'dependency_component_id'])
  x=mi if name=='ISOCRED_MAP' else full
  transformed=payload['preprocessor'].transform(x.loc[ids,payload['features']])
  score=payload['model'].predict_proba(transformed,meta.loc[ids,'dependency_component_id'].to_numpy()) if name=='ISOCRED_MAP' else fi.prob(payload['model'],transformed)
  np.testing.assert_array_equal(score,reference.loc[ids,name])
  checks.append(dict(case='HepG2_final',method=name,rows=len(ids),exact_equal=True))
 train=parquet('output/ISO-002B/models/source_aware_training_reference.parquet')
 q=[]
 for key,value in [('read',mi.loc[ids,['LR_read_R1','LR_read_R2']].mean(axis=1)),('junction',mi.loc[ids,['LR_junction_R1','LR_junction_R2']].mean(axis=1)),('annotation',base.loc[ids,'B2'])]:
  q.append(fi.weighted_ecdf(train[key].to_numpy(),np.ones(len(train)),value.to_numpy()))
 q=np.column_stack(q)
 with np.errstate(invalid='ignore'):
  source=np.mean(np.column_stack([np.nanmean(q[:,:2],axis=1),q[:,2]]),axis=1)
 for name,score in [('SOURCE_AWARE',source),('READ_UNION',q[:,0])]+[(n,base.loc[ids,n].where(meta.loc[ids,'SQANTI3_candidate_coverage' if n=='SQANTI3_RULE' else 'PSQAN_candidate_coverage']>=1).to_numpy() if n in ('SQANTI3_RULE','PSQAN_STYLE') else base.loc[ids,n].to_numpy()) for n in ('B0','B2','SQANTI3_RULE','PSQAN_STYLE')]:
  np.testing.assert_array_equal(score,reference.loc[ids,name]);checks.append(dict(case='HepG2_final',method=name,rows=len(ids),exact_equal=True))
 # Complete matched-input OOF replay using the shipped training-only preprocessors.
 features=parquet('output/ISO-002C/data/minimal_development_features.parquet')
 folds=parquet('output/ISO-002C/data/development_fold_metadata.parquet')
 for name in ('RF14','EN14','L2_READ','L2_READ_ANNOTATION','L2_READ_ANNOTATION_OPPORTUNITY'):
  score=pd.Series(np.nan,index=features.index)
  for fold in range(5):
   p=joblib.load(ROOT/f'output/ISO-002C/models/{name}_fold{fold}.joblib');use=folds.inner_fold==fold
   assert not set(p['training_component_ids']).intersection(folds.loc[use,'dependency_component_id'])
   if p['status']!='ACCEPTED':continue
   score.loc[use]=p['model'].predict_proba(p['preprocessor'].transform(features.loc[use,p['features']]))[:,1]
  reference=parquet(f'output/ISO-002C/results/{name}_OOF.parquet')[name].reindex(features.index)
  np.testing.assert_array_equal(score,reference);checks.append(dict(case='HepG2_post_holdout_development_only',method=name,rows=len(score),exact_equal=True))
 # WTC11 scores are the mean of each family's five fixed HepG2 models.
 features=parquet('output/ISO-002D/data/WTC11_common10.parquet')
 reference=parquet('output/ISO-002D/predictions/WTC11_folds_and_mean.parquet')
 for family in ('RF10','EN10'):
  values=[]
  for fold in range(5):
   p=joblib.load(ROOT/f'output/ISO-002D/models/{family}_fold{fold}.joblib');assert p['status']=='ACCEPTED'
   score=p['model'].predict_proba(p['preprocessor'].transform(features[p['features']]))[:,1]
   np.testing.assert_array_equal(score,reference.loc[features.index,f'{family}_fold{fold}']);values.append(score)
  mean=np.mean(np.stack(values),axis=0);np.testing.assert_array_equal(mean,reference.loc[features.index,family])
  checks.append(dict(case='WTC11_common10',method=family,rows=len(mean),folds=5,exact_equal=True))
 return dict(status='PASS_EXACT_CHECKPOINT_REPLAY',checks=checks,model_fits=0,raw_data_read=False,label_reads_during_replay=False,bootstrap_rerun=False)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('verify','replay'));p.add_argument('--output',type=Path)
 args=p.parse_args();destination=args.output or ROOT/'replay-output'/args.mode
 assert not destination.exists(),'Choose a new output directory; existing results are never overwritten'
 destination.mkdir(parents=True)
 result=verify() if args.mode=='verify' else replay()
 (destination/'check.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
