"""Close OOF prediction generation, bind evaluation, replay without any refit."""
from pathlib import Path
import os,json,subprocess,hashlib,time,sys
R=Path('.');W=R/'work/ISO-002C';O=R/'output/ISO-002C';E=W/'evaluation'
PY=R/'env/ISO-000C-R1/sqanti3/bin/python'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','GOTO_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','POLARS_MAX_THREADS']:env[k]='1'
env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE',PYTHONHASHSEED='20260903',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(W/'tmp'))
assert json.loads((W/'receipts/pipeline_complete.json').read_text())['all_completed']
assert json.loads((W/'review/stat_review.json').read_text())
names=['RF14','EN14','L2_READ','L2_READ_ANNOTATION','L2_READ_ANNOTATION_OPPORTUNITY'];files=[];ninner=nouter=0
for name in names:
 r=json.loads((O/'receipts'/f'{name}_training_complete.json').read_text());ninner+=r['inner_fits'];nouter+=r['outer_refit_attempts']
 for folder,tail,role in [('results','_OOF.parquet','SUPPLEMENTARY_OOF'),('results','_inner_tuning.tsv','FINITE_TUNING'),('results','_fit_parameters.json','FIT_PARAMETERS'),('receipts','_training_complete.json','TRAINING_RECEIPT')]:
  files.append((role,str((O/folder/(name+tail)).relative_to(R))))
 pars=json.loads((O/'results'/f'{name}_fit_parameters.json').read_text())
 for f in pars:
  if f['status'] in ['ACCEPTED','FAILED_NONCONVERGENCE_NO_SCORE']:
   p=O/'models'/f"{name}_fold{f['outer_fold']}.joblib";assert p.is_file()
   files.append(('ACCEPTED_MODEL' if f['status']=='ACCEPTED' else 'FAILED_MODEL_RETAINED',str(p.relative_to(R))))
assert ninner==500 and nouter<=25,(ninner,nouter)
common=['python3',str(W/'code/freeze_manifest.py')]
for state in ['candidate','frozen']:
 p=O/'manifests'/f'OOF_{state}.json'
 cmd=common+['create','--root',str(R),'--output',str(p),'--analysis-id','ISO-002C_SUPPLEMENTARY_OOF_PREDICTIONS','--status',state,'--note','Post-hoc development OOF; not new blind final/external predictions']
 for role,rel in files:cmd+=['--file',role+'='+rel]
 subprocess.run(cmd,check=True,stdout=subprocess.PIPE,env=env)
 c=subprocess.run(common+['check','--root',str(R),'--manifest',str(p),'--strict'],check=True,stdout=subprocess.PIPE,env=env,text=True)
 assert json.loads(c.stdout)['status']=='PASS'
 with (O/'manifests'/f'OOF_{state}_check.json').open('x') as f:f.write(c.stdout)
with (W/'receipts/evaluation_runtime_freeze.json').open('x') as f:json.dump(dict(files={n:sha(W/'code'/n) for n in ['evaluation_guard.py','postprocess.py','infer_frozen.py','replay_check.py']},new_input_generations_frozen=True,inner_fits=ninner,outer_refits=nouter),f,indent=2)
subprocess.run([str(PY),str(E/'evaluate_phase1.py'),'bind-inputs'],env=env,check=True)
for run in ['A','B']:
 label='metrics_'+run
 cmd=['python3',str(W/'code/monitor.py'),'--label',label,'--output-root',str(W),'--hard-cap','8','--interval','.5','python3',str(W/'code/evaluation_guard.py'),label,str(E/('results_'+run))]
 with (W/'logs'/f'{label}_console.log').open('xb') as f:result=subprocess.run(cmd,env=env,cwd=W,stdout=f,stderr=subprocess.STDOUT)
 assert result.returncode==0,dict(run=run,returncode=result.returncode)
 print('METRICS_COMPLETE',run,flush=True)
a=json.loads((E/'results_A/evaluation_complete.json').read_text())
checks=[]
for entry in a['files']:
 p=R/entry['path'];q=E/'results_B'/p.name
 assert sha(p)==sha(q),p.name
 checks.append(dict(file=p.name,A_sha256=sha(p),B_sha256=sha(q),exact=True))
with (W/'receipts/independent_statistics_replay.json').open('x') as f:json.dump(dict(status='PASS',checks=checks),f,indent=2)
with (W/'logs/replay_console.log').open('xb') as f:
 subprocess.run([str(PY),str(W/'code/replay_check.py')],env=env,check=True,stdout=f,stderr=subprocess.STDOUT)
with (W/'receipts/postprocess_complete.json').open('x') as f:json.dump(dict(status='PASS_SCORES_AND_STATISTICS_REPLAYED',new_external_outcome_reads=0,new_final_evaluation=False),f,indent=2)
print('POSTPROCESS_COMPLETE',flush=True)
