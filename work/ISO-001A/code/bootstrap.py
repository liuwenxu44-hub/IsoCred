import csv,hashlib,json,os,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
R=Path('.');W=R/'work/ISO-001A'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
assert not W.exists() and not (R/'.iso001a_writer_lock').exists()
for d in ['reports','contracts','logs','manifests']:
 assert not list((R/d).glob('ISO-001A*'))
roots={'G-R1':('logs/ISO-000G-R1_freeze_manifest.json','75b4338cec96fa4a452ed085a05dc77532b157b83381120e3718ec8b67869050'),'G':('logs/ISO-000G_freeze_manifest.json','b8dfa6f471ed74335aa22c560c073a3ea1c904b70ce8ad5f0fa8aacc8b4dfd9e'),'F':('logs/ISO-000F_freeze_manifest.json','f580c23e7e11688879b21332a29235dc3504ace6d18561d3842a322063100c3c'),'R2B-R1':('logs/ISO-000D-R2B-R1_freeze_manifest.json','647dd707d59d8831bb3693b955005a8b37c966378bbe59a150c4cf24f485e671')}
reg={};man={}
for name,(path,h) in roots.items():
 assert sha(R/path)==h
 man[name]=json.loads((R/path).read_text())
 entries=man[name].get('files')
 if entries is None:
  release=man[name].get('release_manifest',man[name].get('release_manifest_path'))
  assert sha(R/release)==man[name]['release_manifest_sha256']
  entries=list(csv.DictReader((R/release).open(),delimiter='\t'))
 for x in entries:
  if x['path'] in reg:assert reg[x['path']]['sha256']==x['sha256'],x['path']
  reg[x['path']]=dict(x,binding_manifest=path)
required={'candidate':'manifests/ISO-000D-R2B-R1_CANDIDATE_UNIVERSE_FROZEN.tsv','core':'data/features/ISO-000G-R1/ISO-000G-R1_candidate_core_features.parquet','split':'manifests/ISO-000G_blind_split_manifest.tsv','membership':'manifests/ISO-000G_candidate_promoter_group_membership.tsv','holdout':'manifests/ISO-000D-R2B-R1_candidate_holdout_map.tsv','raw':'manifests/ISO-000D-R2B-R1_raw_models.tsv','crosswalk':'manifests/ISO-000D-R2B-R1_exact_cluster_crosswalk.tsv','visible':'manifests/ISO-000F_development_visible_promoter_groups.tsv','old_oof':'data/canary/ISO-000G-R1/synthetic_OOF.parquet','old_fixed':'data/features/ISO-000G/fixed_baseline_scores.parquet','model_utils':'scripts/iso000g_model_utils.py','grid':'contracts/ISO-000G_supervised_baseline_grid.yaml','null_contract':'contracts/ISO-000G_synthetic_canary_contract.yaml','locked_registry':'manifests/ISO-000F_E4_locked_assets.tsv'}
for name in ['B0_read_fixed','B1_junction_fixed','B2_annotation_fixed']:required[name]=f'data/baselines/ISO-000G-R1/{name}.parquet'
for i,(p,x) in enumerate(reg.items()):
 if p.endswith('_classification.txt') and 'SQANTI3_' in p:required[f'sqanti_{i}']=p
for key,expected in [('candidate','31ed901689f0ce70039aa1dd6a80cdb2d5a3001441c3a3550947c0eb56cf6c9d'),('core','17dae2cff24022c9146bc6d0193c25d4f71eb2a53db2908d6c69937b31279d67'),('split','eb216b7d30b4dde15df87a2304a9a3e9891eef3da1a2e41c42269e20d527c23a')]:assert reg[required[key]]['sha256']==expected
checked=[]
for k,p in required.items():
 x=reg[p];assert sha(R/p)==x['sha256'],p;checked.append(dict(input_id=k,**x))
for p,h in [('PROJECT_PLAN_IsoCred.md','7f3619608206c1db5d3100fa975212bd2b9c712c229e666b236fe3146dd3c789'),('upstream_sources/Top5_Q1_algorithm_research_plans_V4_final_freeze_specification.docx','7370cfbc3f01ef9ead47284747f4bbe32a1d8f997308ec4235d3c3a6b08ef355')]:assert sha(R/p)==h
(R/'.iso001a_writer_lock').mkdir();W.mkdir()
for d in ['code','inputs','output','metadata','logs','receipts','tmp']: (W/d).mkdir()
for x in checked:(W/'inputs'/x['input_id']).symlink_to(R/x['path'])
(W/'registry.json').write_text(json.dumps({'bindings':checked,'historical_manifests':roots,'locked_assets_metadata_only':list(csv.DictReader((R/required['locked_registry']).open(),delimiter='\t'))},indent=2)+'\n')
now=datetime.now(timezone.utc).isoformat();owner={'stage':'ISO-001A_Q1_SCIENTIFIC_VALUE_GATE','started_utc':now,'maximum_workdays':10,'technical_repair_workdays':2,'authority':'user message 2026-09-07','G_R2_completed':False}
(R/'.iso001a_writer_lock/owner.json').write_text(json.dumps(owner,indent=2)+'\n')
(W/'receipts/bootstrap.json').write_text(json.dumps(dict(owner,status='PASS',inputs=len(checked),host=os.uname().nodename,process_snapshot=subprocess.check_output(['ps','-u','lwx','-o','pid,ppid,pgid,etime,nlwp,args'],text=True)),indent=2)+'\n')
print(json.dumps(dict(status='PASS',inputs=len(checked),sqanti_files=sum(x['input_id'].startswith('sqanti_') for x in checked))))
