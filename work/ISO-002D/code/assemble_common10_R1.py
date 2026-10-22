"""Assemble common10 only from closed P1 branches; no outcome or model access."""
from common import *
from operational_recovery_R1 import validate_recovered_guard
from collections import defaultdict,Counter
from datetime import datetime,timezone
import sys,time,json
KEYS=['contig','strand','donor_boundary','acceptor_boundary']
S=W/'structure';N=W/'neutral'
def read_tsv(path):return pd.read_csv(path,sep='\t',keep_default_na=False)
def assemble_values(groups,junctions,counts):
 ids=groups.index
 assert not ids.duplicated().any() and not junctions.duplicated(['promoter_group_id']+KEYS).any()
 assert set(junctions.promoter_group_id)<=set(ids)
 jn=junctions.groupby('promoter_group_id').size().reindex(ids,fill_value=0)
 np.testing.assert_array_equal(jn.to_numpy(),groups.junction_unique_count.to_numpy())
 ev=pd.DataFrame({'junction_unique_count':jn},index=ids)
 for rep in ['R1','R2']:
  cnt=counts[rep]
  assert not cnt.duplicated(KEYS).any()
  a=cnt.unique_read_support.to_numpy(dtype=float)
  assert np.isfinite(a).all() and (a>=1).all() and (a==np.floor(a)).all()
  q=junctions.merge(cnt[KEYS+['unique_read_support']],on=KEYS,how='left',validate='many_to_one')
  q['unique_read_support']=q.unique_read_support.fillna(0)
  ev['junction_median_'+rep]=q.groupby('promoter_group_id').unique_read_support.median().reindex(ids)
  ev['junction_fraction_'+rep]=q.assign(supported=q.unique_read_support>0).groupby('promoter_group_id').supported.mean().reindex(ids)
 base=pd.DataFrame({'B2':groups.B2_annotation_fixed_group_mean},index=ids)
 return common10(ev,base,groups),ev
class DSU:
 def __init__(self,ids):self.parent={g:g for g in ids};self.merges=0
 def find(self,g):
  p=self.parent
  while p[g]!=g:p[g]=p[p[g]];g=p[g]
  return g
 def union(self,a,b):
  a,b=self.find(a),self.find(b)
  if a!=b:self.parent[max(a,b)]=min(a,b);self.merges+=1
 def join(self,ids):
  ids=list(ids)
  for g in ids[1:]:self.union(ids[0],g)
def record(path):return dict(path=str(path.relative_to(R)),bytes=path.stat().st_size,sha256=sha(path))
def main():
 start=time.monotonic()
 recovery_entry=validate_recovered_guard()
 requirement=W/'receipts/P1_operational_acceptance_requirement_R1.json'
 req=json.loads(requirement.read_text())
 guard_freeze=R/req['guard_freeze'];guard_complete=R/req['guard_completion']
 gate=json.loads(guard_complete.read_text());assert gate['status']=='PASS'
 deadline=datetime.fromisoformat(req['original_input_deadline_utc']).timestamp()
 assert datetime.fromisoformat(gate['finished_utc']).timestamp()<=deadline and time.time()<=deadline
 gf=json.loads(guard_freeze.read_text());assert sha(R/gf['guard_path'])==gf['guard_sha256']
 sm=S/'branch_frozen.json';sc=S/'branch_frozen_check.json'
 nm=N/'receipts/branch_closed_inventory.json';nc=N/'receipts/branch_inventory_check.json'
 ss=json.loads(sc.read_text());nn=json.loads(nc.read_text())
 assert ss['status']=='PASS' and ss['failed_file_count']==0 and ss['manifest_sha256']==sha(sm)
 assert nn['status']=='PASS' and nn['inventory_sha256']==sha(nm)
 sb={v['path']:v for v in json.loads(sm.read_text())['files']}
 nb={v['path']:v for v in json.loads(nm.read_text())['files']}
 sourcefiles=[recovery_entry]+[record(p) for p in [requirement,guard_freeze,guard_complete]]
 def bound(path,bindings):
  rel=str(path.relative_to(R));v=bindings[rel]
  assert path.stat().st_size==v['bytes'] and sha(path)==v['sha256'],rel
  sourcefiles.append(record(path));return path
 # Every consumed table/observation is checked against its one authoritative branch.
 sr={name:bound(S/'build1'/name,sb) for name in [
  'promoter_groups.tsv','common10_structure_transformed6.tsv','group_junction_union.tsv',
  'metadata_dependency_components.tsv','group_exact_TSS_positions.tsv','candidate_group_membership.tsv']}
 counts={};obs={}
 for rep in ['R1','R2']:
  for name in ['measurement_complete.json','independent_replay.json']:
   rp=bound(N/'measurements'/rep/name,nb)
   assert json.loads(rp.read_text())['status']=='PASS'
  cp=bound(N/'measurements'/rep/'junction/junction_counts.parquet',nb)
  counts[rep]=pd.read_parquet(cp)
  prefix=str((N/'measurements'/rep/'junction/observations').relative_to(R))+'/'
  obs[rep]=[bound(R/p,nb) for p in sorted(nb) if p.startswith(prefix) and p.endswith('.parquet')]
  assert obs[rep]
 for name in ['measurement_complete.json','acquisition_complete.json']:
  rp=bound(N/'receipts'/name,nb);assert json.loads(rp.read_text())['status']=='PASS'
 assert json.loads((N/'receipts/measurement_complete.json').read_text())['zero_support_interpretation_authorized']
 for p in [sm,sc,nm,nc,R/'contracts/ISO-002D_common10_external_scope.yaml',
  W/'independence/peer_adjudication.json',W/'code/common.py',Path(__file__)]:sourcefiles.append(record(p))
 js(W/'receipts/assembly_input_code_freeze.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
  files=sourcefiles,real_outcomes_read=0,model_fits=0,
  measurement='Original group union unique-junction median and supported fraction; zero only after complete neutral measurement',
  graph='Start structural gene/block/junction/TSS200 components, then join observed shared-read junction components across all contigs per biological replicate. Bootstrap only; no feature/candidate alteration.'))
 groups=read_tsv(sr['promoter_groups.tsv']).set_index('promoter_group_id')
 assert (groups.global_holdout==False).all() and (groups.exposure_exclusion==False).all()
 junctions=read_tsv(sr['group_junction_union.tsv'])
 x,ev=assemble_values(groups,junctions,counts)
 structural=read_tsv(sr['common10_structure_transformed6.tsv']).set_index('promoter_group_id')
 pd.testing.assert_index_equal(x.index,structural.index)
 np.testing.assert_allclose(x[FEATURES[4:]].to_numpy(),structural[FEATURES[4:]].to_numpy(),rtol=1e-14,atol=1e-14)
 d=DSU(groups.index)
 meta=read_tsv(sr['metadata_dependency_components.tsv']).set_index('promoter_group_id')
 assert set(meta.index)==set(groups.index)
 for _,part in meta.groupby('metadata_dependency_component_id'):d.join(part.index)
 initial_components=len({d.find(g) for g in groups.index})
 jrepresentative={}
 for key,part in junctions.groupby(KEYS,sort=False):
  d.join(part.promoter_group_id);jrepresentative[key]=part.promoter_group_id.iloc[0]
 read_audit=[]
 for rep in ['R1','R2']:
  first={};read_components=0;linked_observations=0;before=d.merges
  for path in obs[rep]:
   ob=pq.ParquetFile(path).read(columns=['read_key']+KEYS).to_pandas()
   assert ob.read_key.str.startswith(rep+'|').all()
   for rid,ctg,strand,donor,acceptor in ob.itertuples(index=False,name=None):
    g=jrepresentative.get((ctg,strand,donor,acceptor))
    if g is None:continue
    linked_observations+=1
    if rid in first:
     if d.find(first[rid])!=d.find(g):read_components+=1;d.union(first[rid],g)
    else:first[rid]=g
  read_audit.append(dict(replicate=rep,observed_reads_supporting_group_junctions=len(first),
   supporting_junction_observations=linked_observations,component_union_events=d.merges-before,
   intercomponent_read_events=read_components,readmap_not_used=True))
  print(json.dumps(read_audit[-1]),flush=True)
 components=defaultdict(list)
 for g in groups.index:components[d.find(g)].append(g)
 compids={}
 for members in components.values():
  cid='WTC11_FULLDEP_'+hashlib.sha256(json.dumps(sorted(members),separators=(',',':')).encode()).hexdigest()[:24]
  for g in members:compids[g]=cid
 comp=pd.DataFrame({'promoter_group_id':groups.index,'dependency_component_id':[compids[g] for g in groups.index]})
 comp['component_group_count']=comp.groupby('dependency_component_id').promoter_group_id.transform('size')
 assert time.time()<=deadline,'P1_INPUT_DEADLINE_BEFORE_OUTPUT'
 # All structurally legal groups remain, regardless of future assay eligibility.
 parq(O/'data/WTC11_common10.parquet',x.reset_index())
 parq(O/'data/WTC11_group_evidence_label_blind.parquet',ev.reset_index())
 parq(O/'data/WTC11_group_metadata_label_blind.parquet',groups.reset_index())
 parq(O/'data/WTC11_dependency_components.parquet',comp)
 parq(O/'data/WTC11_simple_diagnostics.parquet',diagnostics(x).reset_index())
 outputfiles=[record(O/'data'/n) for n in ['WTC11_common10.parquet','WTC11_group_evidence_label_blind.parquet','WTC11_group_metadata_label_blind.parquet','WTC11_dependency_components.parquet','WTC11_simple_diagnostics.parquet']]
 assert time.time()<=deadline,'P1_INPUT_DEADLINE_AFTER_OUTPUT'
 js(W/'receipts/assembly_complete.json',dict(status='PASS_PENDING_PARENT_P1_ACCEPTANCE',
  finished_utc=datetime.now(timezone.utc).isoformat(),input_deadline_utc=req['original_input_deadline_utc'],
  features=FEATURES,structurally_legal_groups=len(x),scored_input_groups=len(x),
  monoexon_only_groups=int(x.status_junction_NOT_APPLICABLE.sum()),missing_values=x.isna().sum().to_dict(),
  junction_union_links=len(junctions),initial_metadata_components=initial_components,
  final_dependency_components=len(components),largest_component=max(map(len,components.values())),
  component_size_quantiles=np.quantile(list(map(len,components.values())),[0,.25,.5,.75,.95,.99,1]).tolist(),
  read_dependency_audit=read_audit,files=outputfiles,model_fits=0,outcomes_read=0,
  scientific_claim='COMMON10_MEASUREMENT_COMPLETED_NOT_EXTERNAL_VALIDATION',
  wall_seconds=time.monotonic()-start,
  CPU_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
  peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
if __name__=='__main__':main()
