"""Independent E4 extractor: immediately skip non-development IDs before using values."""
from common import *
from collections import defaultdict,Counter
from datetime import datetime,timezone
import csv,gzip
start=time.monotonic();contract=O/'contracts/ISO-001A_real_data_comparison_contract.yaml';fr=json.loads((O/'receipts/protocol_freeze.json').read_text());assert sha(contract)==fr['contract_sha256'];assert sha(W/'code/extract_development.py')==fr['science_code']['extract_development.py']
assert not (O/'data/development_E4_v1.parquet').exists()
meta=pd.read_parquet(O/'data/group_metadata_label_blind.parquet').set_index('promoter_group_id');assert len(meta)==46408 and (meta.split=='DEVELOPMENT_POOL').all();dev=set(meta.index)
reg=json.loads((W/'registry.json').read_text());locked={Path(r['asset_path']).name:r for r in reg['locked_assets_metadata_only']}
outcome=locked['promoter_group_outcomes.tsv.gz'];links=locked['CTSS_to_promoter_group_links.tsv.gz']
for a in [outcome,links]:assert sha(R/a['asset_path'])==a['sha256']
rows=[]
with gzip.open(R/outcome['asset_path'],'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  g=r['promoter_group_id']
  if g not in dev:continue
  assert r['holdout_status']=='NON_HOLDOUT' and r['exposure_exclusion_status']=='NOT_EXCLUDED'
  m=meta.loc[g];assert (r['assay_eligibility']=='true')==bool(m.assay_eligibility) and r['mapping_ambiguity']==m.mapping_ambiguity
  reasons=[]
  if not m.assay_eligibility:reasons.append('ASSAY_INELIGIBLE')
  if m.mapping_ambiguity=='ALL_COMPATIBLE_CTSS_PROMOTER_GROUP_SHARED':reasons.append('ALL_SHARED_UNRESOLVED')
  if m.touches_final_dependency:reasons.append('DEPENDENCE_TOUCHES_FINAL')
  eligible=not reasons
  primary=int(r['primary_anchor_status']=='E4_REPRODUCIBLE_PROMOTER_SUPPORT');s2=int(r['S2_POSITION_CONCORDANT_5BP']=='true')
  if eligible:
   assert primary==int(int(r['replicate_detection_count'])>=2)
   assert s2<=primary
  rows.append(dict(promoter_group_id=g,evaluation_eligible=eligible,exclusion_reason=';'.join(reasons) or 'ELIGIBLE',primary_target=primary if eligible else np.nan,S2_target=s2 if eligible else np.nan,legacy_outcome_status=r['outcome_status'],legacy_primary_target=primary,legacy_S2_target=s2))
table=pd.DataFrame(rows).set_index('promoter_group_id').loc[meta.index];assert len(table)==46408 and not table.index.duplicated().any()
# Actual shared CTSS edges must be contained by the label-blind dependence graph.
obs={};nlinks=0;cross=0;bad=0
with gzip.open(R/links['asset_path'],'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  g=r['promoter_group_id']
  if g not in dev:continue
  if r['within_50bp']!='true':continue
  nlinks+=1;c=meta.at[g,'dependency_component_id'];oid=r['ctss_observation_id']
  if oid in obs:
   pc,pg=obs[oid]
   if pg!=g:cross+=1;bad+=int(pc!=c)
  else:obs[oid]=(c,g)
assert bad==0,('CTSS dependence missing from frozen graph',bad)
p=parq('data/development_E4_v1.parquet',table.reset_index());p.chmod(0o444)
audit=dict(created_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=sha(contract),extractor_sha256=sha(W/'code/extract_development.py'),source_hashes={Path(a['asset_path']).name:a['sha256'] for a in [outcome,links]},membership_groups=len(table),eligible_groups=int(table.evaluation_eligible.sum()),eligible_primary_positive=int(table.primary_target.sum()),eligible_S2_positive=int(table.S2_target.sum()),exclusion_patterns=table.exclusion_reason.value_counts().to_dict(),legacy_development_ambiguous=int((table.legacy_outcome_status=='E4_AMBIGUOUS_MAPPING').sum()),development_primary_CTSS_links=nlinks,development_shared_CTSS_edges=cross,CTSS_edges_crossing_dependency_components=bad,output_sha256=sha(p),E4_development_visibility='UNLOCKED_DEVELOPMENT_ONLY',FINAL_BLIND_TEST='STILL_LOCKED_NONDEVELOPMENT_ROWS_SKIPPED',global_external_holdout='LOCKED',exposure_exclusions='LOCKED',Davies_A375_outcomes='LOCKED',full_E4_access_consumer='INDEPENDENT_EXTRACTOR_ONLY')
js('receipts/development_unlock_receipt.json',audit);print(json.dumps(audit,default=str),flush=True);receipt('extract_development',start,output_sha256=sha(p),development_groups=len(table),final_outcome_rows_exported=0)
