import os,time,json,hashlib,resource
from pathlib import Path
import numpy as np,pandas as pd,pyarrow as pa,pyarrow.parquet as pq
from threadpoolctl import threadpool_limits,threadpool_info
threadpool_limits(1);pa.set_cpu_count(1);pa.set_io_thread_count(1)
R=Path('.');W=R/'work/ISO-002B';I=W/'inputs';O=R/'output/ISO-002B';SEED=20260903
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def js(p,x):
 p=O/p;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,default=str,sort_keys=True)+'\n');return p
def report(name,s):
 p=O/'reports'/('ISO-002B_'+name+'.md');p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s.rstrip()+'\n')
def tsv(p,x):
 p=O/p;p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(x).to_csv(p,sep='\t',index=False);return p
def parq(p,x):
 p=O/p;p.parent.mkdir(parents=True,exist_ok=True);pq.write_table(pa.Table.from_pandas(x,preserve_index=False),p,compression='zstd');return p
def hnum(s):return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8],'big')/2**64
def ids(s):return [] if pd.isna(s) or str(s) in ['NA','nan','','.','NOT_APPLICABLE'] else str(s).replace(',',';').split(';')
def receipt(name,start,**kw):
 u=resource.getrusage(resource.RUSAGE_SELF);x=dict(kw,wall_seconds=time.monotonic()-start,cpu_user_seconds=u.ru_utime,cpu_system_seconds=u.ru_stime,RSS_bytes=u.ru_maxrss*1024,pools=threadpool_info(),arrow_CPU=pa.cpu_count(),arrow_IO=pa.io_thread_count());js('receipts/'+name+'.json',x);print(json.dumps(x,default=str),flush=True)
