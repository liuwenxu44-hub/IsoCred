"""Phase-1 helpers. All data and outputs are project-contained."""
import os,time,json,hashlib,resource
from pathlib import Path
import numpy as np,pandas as pd,pyarrow as pa,pyarrow.parquet as pq
from threadpoolctl import threadpool_limits
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler
from sklearn.pipeline import make_pipeline
threadpool_limits(1);pa.set_cpu_count(1);pa.set_io_thread_count(1)
R=Path('.');W=R/'work/ISO-002C';O=R/'output/ISO-002C';SEED=20260903
EVIDENCE='POST_HOC_SECONDARY_AFTER_FINAL_EXPOSURE_DEVELOPMENT_GROUPED_OOF'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def js(p,x):
 p=Path(p)
 with p.open('x') as f:json.dump(x,f,indent=2,ensure_ascii=False,default=str);f.write('\n')
 return p
def parq(p,x):
 p=Path(p);assert not p.exists(),str(p)
 pq.write_table(pa.Table.from_pandas(x,preserve_index=False),p,compression='zstd');return p
def prep():return make_pipeline(SimpleImputer(strategy='median',fill_value=0,add_indicator=True,keep_empty_features=True),RobustScaler())
def hnum(s):return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8],'big')/2**64
def minimal(ev,base,meta):
 x=pd.DataFrame(index=meta.index)
 for rep in ['R1','R2']:
  x['LR_read_'+rep]=np.log1p(ev['read_union_'+rep])
  x['LR_junction_'+rep]=np.log1p(ev['junction_median_'+rep])
  x['LR_junction_fraction_'+rep]=ev['junction_fraction_'+rep]
  x['resolution_read_shared_'+rep]=ev['read_crossgroup_fraction_'+rep]
 x['annotation']=base.B2/2
 x['resolution_candidate_count']=np.log1p(meta.candidate_count)
 x['resolution_TSS_count']=np.log1p(meta.exact_TSS_count)
 x['resolution_window']=np.log1p(meta.window_opportunity_bp)
 x['resolution_monoexon']=meta.monoexon_fraction
 x['status_junction_NOT_APPLICABLE']=(ev.junction_unique_count==0).astype(float)
 return x
def load_development():
 d={n:pd.read_parquet(R/'work/ISO-001A/output/data'/('group_'+n+'_label_blind.parquet')).set_index('promoter_group_id') for n in ['core','metadata','evidence','baselines']}
 labels=pd.read_parquet(R/'work/ISO-001A/output/data/development_E4_v1.parquet').set_index('promoter_group_id')
 for v in d.values():assert v.index.equals(labels.index) and not v.index.duplicated().any()
 ids=labels.index[labels.evaluation_eligible]
 assert len(labels)==46408 and len(ids)==30487
 d={k:v.loc[ids].copy() for k,v in d.items()};labels=labels.loc[ids]
 assert (d['metadata'].split=='DEVELOPMENT_POOL').all()
 assert not d['metadata'].touches_final_dependency.any()
 x=minimal(d['evidence'],d['baselines'],d['metadata'])
 d.update(labels=labels,minimal=x);return d
def validate_contract():
 p=R/'contracts/ISO-002C_phase1_comparison_contract.yaml'
 fr=json.loads((W/'receipts/protocol_freeze.json').read_text())
 assert sha(p)==fr['contract_sha256']
 for rel,h in fr['code_sha256'].items():assert sha(R/rel)==h,rel
 for rel,h in fr['training_input_sha256'].items():assert sha(R/rel)==h,rel
 return json.loads(p.read_text())
