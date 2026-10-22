"""Reusable ISO-002D functions; no data is opened on import."""
import os,json,hashlib,resource
from pathlib import Path
import numpy as np,pandas as pd,pyarrow as pa,pyarrow.parquet as pq
from threadpoolctl import threadpool_limits
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler
from sklearn.pipeline import make_pipeline
threadpool_limits(1);pa.set_cpu_count(1);pa.set_io_thread_count(1)
R=Path('.');W=R/'work/ISO-002D';O=R/'output/ISO-002D'
SEED=20260903
FEATURES=['LR_junction_R1','LR_junction_fraction_R1','LR_junction_R2','LR_junction_fraction_R2','annotation','resolution_candidate_count','resolution_TSS_count','resolution_window','resolution_monoexon','status_junction_NOT_APPLICABLE']
EVIDENCE='POST_HOC_AFTER_HISTORICAL_FINAL_EXPOSURE_DEVELOPMENT_ONLY'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def js(p,x):
 with Path(p).open('x') as f:json.dump(x,f,indent=2,ensure_ascii=False,default=str);f.write('\n')
def parq(p,x):
 p=Path(p);assert not p.exists(),str(p)
 pq.write_table(pa.Table.from_pandas(x,preserve_index=False),p,compression='zstd')
def prep():return make_pipeline(SimpleImputer(strategy='median',fill_value=0,add_indicator=True,keep_empty_features=True),RobustScaler())
def hnum(s):return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8],'big')/2**64
def common10(ev,base,meta):
 x=pd.DataFrame(index=meta.index)
 assert ev.index.equals(meta.index) and base.index.equals(meta.index)
 count=ev['junction_unique_count'].to_numpy(dtype=float)
 assert np.isfinite(count).all() and (count>=0).all() and (count==np.floor(count)).all()
 assert np.array_equal(count==0,meta['monoexon_fraction'].to_numpy(dtype=float)==1)
 for rep in ['R1','R2']:
  x['LR_junction_'+rep]=np.log1p(ev['junction_median_'+rep])
  x['LR_junction_fraction_'+rep]=ev['junction_fraction_'+rep]
 x['annotation']=base['B2']/2
 x['resolution_candidate_count']=np.log1p(meta['candidate_count'])
 x['resolution_TSS_count']=np.log1p(meta['exact_TSS_count'])
 x['resolution_window']=np.log1p(meta['window_opportunity_bp'])
 x['resolution_monoexon']=meta['monoexon_fraction']
 x['status_junction_NOT_APPLICABLE']=(ev['junction_unique_count']==0).astype(float)
 assert x.columns.tolist()==FEATURES
 validate_features(x)
 return x
def validate_features(x):
 assert x.columns.tolist()==FEATURES and not x.index.duplicated().any()
 assert not np.isinf(x.to_numpy(dtype=float)).any()
 na=x.status_junction_NOT_APPLICABLE.to_numpy()==1
 assert set(x.status_junction_NOT_APPLICABLE.unique())<={0.,1.}
 for rep in ['R1','R2']:
  a=x['LR_junction_'+rep].to_numpy();f=x['LR_junction_fraction_'+rep].to_numpy()
  assert np.isnan(a[na]).all() and np.isnan(f[na]).all()
  assert np.isfinite(a[~na]).all() and (a[~na]>=0).all()
  assert np.isfinite(f[~na]).all() and ((f[~na]>=0)&(f[~na]<=1)).all()
 for c in FEATURES[4:]:assert np.isfinite(x[c]).all()
 for c in ['annotation','resolution_monoexon']:assert x[c].between(0,1).all()
 assert (x.resolution_candidate_count>=np.log(2)).all() and (x.resolution_TSS_count>=np.log(2)).all()
 assert (x.resolution_window>=np.log(102)-1e-12).all()
def diagnostics(x):
 return pd.DataFrame({'ANNOTATION_ONLY':x.annotation,'JUNCTION_MEDIAN_SUPPORT':x[['LR_junction_R1','LR_junction_R2']].mean(axis=1,skipna=False),'WINDOW_OPPORTUNITY_ONLY':x.resolution_window},index=x.index)
def load_development():
 d={n:pd.read_parquet(R/'work/ISO-001A/output/data'/('group_'+n+'_label_blind.parquet')).set_index('promoter_group_id') for n in ['metadata','evidence','baselines']}
 labels=pd.read_parquet(R/'work/ISO-001A/output/data/development_E4_v1.parquet').set_index('promoter_group_id')
 for v in d.values():assert v.index.equals(labels.index) and not v.index.duplicated().any()
 ids=labels.index[labels.evaluation_eligible]
 assert len(labels)==46408 and len(ids)==30487
 d={k:v.loc[ids].copy() for k,v in d.items()};labels=labels.loc[ids]
 assert (d['metadata'].split=='DEVELOPMENT_POOL').all()
 assert not d['metadata'].touches_final_dependency.any()
 x=common10(d['evidence'],d['baselines'],d['metadata'])
 d.update(labels=labels,features=x);return d
def validate_contract():
 p=R/'contracts/ISO-002D_common10_computation_contract.json'
 fr=json.loads((W/'receipts/computation_freeze.json').read_text())
 assert sha(p)==fr['contract_sha256']
 for k in ['code_sha256','training_input_sha256','gate_sha256']:
  for rel,h in fr[k].items():assert sha(R/rel)==h,rel
 c=json.loads(p.read_text())
 assert c['features']==FEATURES and c['seed']==SEED
 assert c['training']['eligible_groups']==30487 and c['input_gate_status']=='PASS'
 return c
