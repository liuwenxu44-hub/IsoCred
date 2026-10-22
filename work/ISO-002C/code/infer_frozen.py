"""Minimal label-free scoring interface. Does not train, tune or calibrate."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd,joblib
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from threadpoolctl import threadpool_limits
threadpool_limits(1)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def forbidden_fit(*args,**kwargs):raise RuntimeError('Fitting is forbidden in frozen inference')
LogisticRegression.fit=forbidden_fit
RandomForestClassifier.fit=forbidden_fit
def infer(model_path,model_sha,features_path,metadata_path,output):
 assert sha(model_path)==model_sha,'Model hash mismatch'
 output=Path(output);assert not output.exists()
 model=joblib.load(model_path);assert model['status']=='ACCEPTED','Failed model is unscorable'
 x=pd.read_parquet(features_path).set_index('promoter_group_id')
 meta=pd.read_parquet(metadata_path).set_index('promoter_group_id')
 assert not x.index.duplicated().any() and x.index.equals(meta.index)
 assert set(model['features']).issubset(x.columns),'Required input columns absent: never fill them with zero'
 assert not set(meta.dependency_component_id).intersection(model['training_component_ids']),'Training dependency overlap'
 assert not set(x.index).intersection(model['training_group_ids']),'Training group overlap'
 a=model['preprocessor'].transform(x[model['features']])
 assert list(model['model'].classes_)==[0,1]
 pred=model['model'].predict_proba(a)[:,1]
 assert np.isfinite(pred).all()
 out=pd.DataFrame({'promoter_group_id':x.index,'support_event_score':pred,'model_sha256':model_sha,'meaning':'promoter_CAGE_support_event_not_function'})
 out.to_parquet(output,index=False)
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser()
 for n in ['model','model-sha256','features','metadata','output']:ap.add_argument('--'+n,required=True)
 a=ap.parse_args();infer(a.model,a.model_sha256,a.features,a.metadata,a.output)
