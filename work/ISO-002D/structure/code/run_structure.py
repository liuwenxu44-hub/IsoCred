#!/usr/bin/env python3
"""Target-blind fixed GTF structural build, with preexecution code freeze."""
import argparse,csv,gzip,json,math,os,resource,subprocess,sys,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from bootstrap_qc import ROOT,W,CONTRACT,CONTRACT_SHA,REF,GTF,sha,record,js
from structure_lib import *
def tsv(path,rows,fields=None):
 rows=list(rows);fields=fields or list(rows[0])
 def enc(x):
  if isinstance(x,bool):return str(x).lower()
  if isinstance(x,(list,tuple,set,dict)):return json.dumps(sorted(x) if isinstance(x,set) else x,sort_keys=True,separators=(',',':'))
  if x is None:return 'NA'
  if isinstance(x,float):return format(x,'.17g')
  return x
 with path.open('x',newline='') as f:
  dw=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');dw.writeheader()
  dw.writerows({k:enc(r.get(k)) for k in fields} for r in rows)
def main():
 a=argparse.ArgumentParser();a.add_argument('--generation',default='build1');args=a.parse_args()
 assert args.generation in ['build1','replay1']
 out=W/args.generation;assert not out.exists(),out
 assert sha(CONTRACT)==CONTRACT_SHA
 assert json.loads((W/'receipts/QC_exit.json').read_text())['exit_code']==0
 qce=json.loads((W/'receipts/QC_execution_contract.json').read_text())
 fz=json.loads((W/'receipts/QC_execution_contract_freeze.json').read_text())
 assert sha(W/'receipts/QC_execution_contract.json')==fz['sha256']
 # Original code is hash checked before any named safe AST definition is used on data.
 for f in qce['files']:
  if f['role'] in ['ORIGINAL_SOURCE_SEMANTICS','FIXED_INPUT']:
   assert sha(ROOT/f['path'])==f['sha256'],f['path']
 cls=W/'qc/WTC11_FIXED_classification.txt';corr=W/'qc/WTC11_FIXED_corrected.gtf'
 tests=subprocess.run([sys.executable,str(W/'code/test_structure_semantics.py')],capture_output=True,text=True,check=True)
 tr=json.loads(tests.stdout);assert tr['status']=='PASS'
 out.mkdir()
 js(out/'synthetic_semantic_tests.json',tr)
 inputfiles=[record(CONTRACT,'FROZEN_SCOPE'),record(W/'receipts/QC_execution_contract.json','QC_INPUT_CODE_FREEZE'),
  record(cls,'STRUCTURAL_CLASSIFICATION_NO_LABELS'),record(corr,'CORRECTED_COORDINATE_CROSSWALK_ONLY'),
  record(W/'receipts/QC_exit.json','COMPLETED_QC_RECEIPT')]
 inputfiles += [record(W/'code'/n,'BUILD_CODE') for n in ['run_structure.py','structure_lib.py','test_structure_semantics.py','bootstrap_qc.py']]
 js(out/'preexecution_contract.json',dict(stage_contract_sha256=CONTRACT_SHA,created_utc=datetime.now(timezone.utc).isoformat(),
  files=inputfiles,command=[sys.executable,str(Path(__file__)),'--generation',args.generation],
  parent_QC_fixed_inputs=qce['files'],
  semantics=dict(coordinate_authority='ARCHIVE_GTF_1_BASED_CLOSED_TRANSCRIPT_ORDER',gene_functions='ORIGINAL_AST_NAMED_DEFS',
  candidate='EXACT_DEDUP_THEN_ORIGINAL_DETERMINISTIC_COMPLETE_LINK_50BP_MONO_RO_0.80',
  promoter='ORIGINAL_PROMOTER_GROUPS_50BP_COMPLETE_TSS_ENVELOPE',
  windows='ALL_EXACT_MEMBER_TSS_MINUS1_UNION_INCLUSIVE_PLUS_MINUS50',
  annotation='MEAN_OVER_CANDIDATES_OF_PRIORITY_SCORE_EXACT2_CHAIN1_OTHER0_DIVIDED2',
  junctions='ALL_GROUP_MEMBER_DISTINCT_KEYS;PLUS(a,b-1);MINUS(a-1,b)',
  global_holdout='ORIGINAL_850_GENE_SYMBOL_AND_10KB_GENOMIC_BLOCK_ANTISENSE_UNION',
  exposure='ORIGINAL_11_EXPOSURE_GENE_SYMBOL_EXCLUSION',
  coding_field='UNMEASURED_PLACEHOLDER_ORF_DISABLED;NOT_READ_INTO_FEATURES_OR_EXCLUSIONS',
  input_domain='FIXED_TALON4_LAPA_POOLED3_CANDIDATES;NOT_ORIGINAL3CALLER',
  CTSS_read=False,model_training=False,external_preprocessing_fit=False,shared_read_graph=False)))
 js(out/'preexecution_contract_freeze.json',record(out/'preexecution_contract.json','PREEXECUTION_FREEZE'))
 start=time.monotonic();u=resource.getrusage(resource.RUSAGE_SELF)
 ns=reference_index(raw_namespace(),REF)
 print('REFERENCE_INDEX_READY',flush=True)
 tx,malformed=ns['parse_gtf'](GTF);assert not malformed and len(tx)==94486
 roles={r['contig_name']:r['role'] for r in csv.DictReader((ROOT/'manifests/ISO-000C-R2_reference_contigs.tsv').open(),delimiter='\t')}
 human={k for k,v in roles.items() if v=='HUMAN_PRIMARY_ANALYSIS_SEQUENCE'}
 observed=Counter();maprows=0
 with gzip.open(ROOT/'work/ISO-002C-R1/archive/long_only_cDNA_PacBio/read_model_map.tsv.gz','rt') as f:
  reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==['read_id','transcript_id']
  for r in reader:observed[r['transcript_id']]+=1;maprows+=1
 assert set(observed)==set(tx) and maprows==5830665
 sq_byid,corr_id,corr_hash=ns['parse_sqanti'](cls,corr)
 assert len(sq_byid)==94250 and all(len(x)==1 for x in sq_byid.values())
 raw=[];cross=[];reasoncounts=Counter();qcounts=Counter();classcounts=Counter()
 for tid,t in sorted(tx.items()):
  rawkey='ARCHIVAL_TALON4_LAPA|POOL3|'+tid
  chrom=next(iter(t['chroms'])) if len(t['chroms'])==1 else 'MULTICONTIG'
  strand=next(iter(t['strands'])) if len(t['strands'])==1 else 'MULTISTRAND'
  reasons=[];c=None
  if len(t['chroms'])!=1 or len(t['strands'])!=1:reasons.append('explicit_fusion_chimera')
  else:
   try:c=canonical_structure(chrom,strand,t['exons'])
   except ValueError as e:reasons.append(str(e))
  if c and chrom not in human:reasons.append('chrEBV' if chrom=='chrEBV' else 'nonhuman_contig')
  if not observed[tid]:reasons.append('annotation_only_no_observation')
  qrows=sq_byid.get(tid,[]);qrow=qrows[0] if len(qrows)==1 else {}
  qstatus='EXACT_ID_ONE_TO_ONE' if len(qrows)==1 else ('EXACT_ID_MULTIPLE_QC_ROWS' if len(qrows)>1 else 'NO_EXACT_ID_MATCH')
  coordinate_status='NO_CORRECTED_GTF_RECORD';corrids=[]
  if c:
   if tid in corr_id:
    corrids=[tid];coordinate_status='UNCHANGED' if corr_id[tid]==c['exact_structure_hash'] else 'QC_NORMALIZED_DIFFERS_RAW_RETAINED'
   elif corr_hash.get(c['exact_structure_hash']):
    corrids=sorted(corr_hash[c['exact_structure_hash']]);coordinate_status='STRUCTURE_MATCH_DIFFERENT_QC_ID'
    if not qrows and len(corrids)==1 and len(sq_byid.get(corrids[0],[]))==1:
     qrow=sq_byid[corrids[0]][0];qstatus='STRUCTURE_CROSSWALK_ONE_TO_ONE'
    elif len(corrids)>1:qstatus='STRUCTURE_CROSSWALK_AMBIGUOUS'
  if chrom in human:
   assert qstatus in ['EXACT_ID_ONE_TO_ONE','STRUCTURE_CROSSWALK_ONE_TO_ONE'],(tid,qstatus)
   assert qrow.get('structural_category') and qrow.get('associated_gene') is not None,tid
  else:qstatus='NOT_RUN_TECHNICAL_NONHUMAN'
  if str(qrow.get('structural_category','')).lower()=='fusion':reasons.append('explicit_fusion_chimera')
  compat=set()
  for value in (t['attrs'].get('gene_id'),t['attrs'].get('gene_name'),qrow.get('associated_gene')):compat.update(ns['normalize_gene'](value))
  if c:
   span=(min(a for a,b in c['ordered_exons']),max(b for a,b in c['ordered_exons']))
   overlaps={g['gene_id'] for g in ns['gene_intervals'][(chrom,strand)] if g['start']<=span[1] and g['end']>=span[0]}
  else:span=(0,0);overlaps=set()
  if len(compat)==1:lstatus='KNOWN_UNIQUE_ASSOCIATION'
  elif len(compat)>1:lstatus='AMBIGUOUS'
  elif len(overlaps)==1:compat=overlaps;lstatus='KNOWN_UNIQUE_OVERLAP'
  elif len(overlaps)>1:compat=overlaps;lstatus='AMBIGUOUS'
  else:lstatus='NOVEL'
  rel=ns['gencode_relation'](chrom,strand,c,lstatus) if c else 'NA'
  reasons=sorted(set(reasons));reasoncounts.update(reasons);qcounts[qstatus]+=1
  if qrow:classcounts[qrow['structural_category']]+=1
  raw.append(dict(raw_model_key=rawkey,original_transcript_id=tid,original_gene_id=t['attrs'].get('gene_id','NA'),
   raw_gtf_path=str(GTF.relative_to(ROOT)),raw_gtf_sha256='66b12bb0e58285716c3b19d2437206043bd2fa8e564b43c550d71a69cdfc96d7',
   chromosome=chrom,strand=strand,contig_role=roles.get(chrom,'CONTIG_NOT_IN_ORIGINAL_REFERENCE'),
   ordered_exon_chain=c['ordered_exons'] if c else [],junction_chain=c['ordered_internal_junctions'] if c else [],
   TSS=c['TSS'] if c else None,TES=c['TES'] if c else None,exon_count=len(c['ordered_exons']) if c else 0,
   transcript_length=sum(b-a+1 for a,b in c['ordered_exons']) if c else 0,
   exact_structure_hash=c['exact_structure_hash'] if c else 'NA',genomic_span=span,
   compatible_gene_ids=sorted(compat),gene_assignment_status=lstatus,GENCODE_structural_relation=rel,
   SQANTI3_structural_category=qrow.get('structural_category','NA'),SQANTI3_crosswalk_status=qstatus,
   observation_status='OBSERVED_TRANSCRIPT_MODEL' if observed[tid] else 'ANNOTATION_ONLY',
   archived_assigned_read_count_provenance_only=observed[tid],technical_eligibility=not reasons,technical_exclusion_reason=reasons))
  cross.append(dict(raw_model_key=rawkey,original_transcript_id=tid,crosswalk_status=qstatus,
   QC_record_ids=corrids or ([tid] if qrows else []),QC_record_count=len(qrows),
   raw_exact_structure_hash=c['exact_structure_hash'] if c else 'NA',QC_coordinate_status=coordinate_status,
   SQANTI3_structural_category=qrow.get('structural_category','NA'),SQANTI3_subcategory=qrow.get('subcategory','NA'),
   SQANTI3_associated_gene=qrow.get('associated_gene','NA'),SQANTI3_associated_transcript=qrow.get('associated_transcript','NA'),
   candidate_coordinate_source='RAW_ARCHIVED_GTF_ONLY'))
 print('RAW_CLASSIFICATION_CROSSWALK_READY',flush=True)
 exact=exacts_from_raw(raw);cand,cw=candidates_from_exacts(exact)
 hold=list(csv.DictReader((ROOT/'manifests/ISO-000C_global_holdout_reconciliation.tsv').open(),delimiter='\t'))
 ctx=original_holdout_blocks(REF,hold,ns['attrs']);holds=candidate_holdouts(cand,raw,ctx)
 groups,positions,junctions,members=groups_and_features(cand,exact,cw,holds)
 for name,rows in [('raw_models.tsv',raw),('sqanti_crosswalk.tsv',cross),('exact_structures.tsv',exact),
  ('candidates_50bp.tsv',cand),('exact_candidate_crosswalk.tsv',cw),('candidate_holdout_map.tsv',holds),
  ('reference_genomic_blocks.tsv',sorted(ctx['blocks'],key=lambda x:x['block_id'])),
  ('promoter_groups.tsv',groups),('group_exact_TSS_positions.tsv',positions),('group_junction_union.tsv',junctions),
  ('candidate_group_membership.tsv',members)]:tsv(out/name,rows)
 base=[{k:r[k] for k in ['promoter_group_id','B2_annotation_fixed_group_mean','candidate_count','exact_TSS_count',
  'window_opportunity_bp','monoexon_fraction','junction_unique_count']} for r in groups]
 tsv(out/'common10_structure_primitives.tsv',base)
 transformed=[dict(promoter_group_id=r['promoter_group_id'],annotation=r['B2_annotation_fixed_group_mean']/2,
  resolution_candidate_count=math.log1p(r['candidate_count']),resolution_TSS_count=math.log1p(r['exact_TSS_count']),
  resolution_window=math.log1p(r['window_opportunity_bp']),resolution_monoexon=r['monoexon_fraction'],
  status_junction_NOT_APPLICABLE=float(r['junction_unique_count']==0)) for r in groups]
 tsv(out/'common10_structure_transformed6.tsv',transformed)
 # Metadata-only conservative dependence; no imported HepG2 graph or CTSS reads.
 par={g['promoter_group_id']:g['promoter_group_id'] for g in groups}
 def find(g):
  while par[g]!=g:par[g]=par[par[g]];g=par[g]
  return g
 def union(a,b):
  a,b=find(a),find(b)
  if a!=b:par[max(a,b)]=min(a,b)
 def join(gs):
  gs=list(gs)
  for b in gs[1:]:union(gs[0],b)
 byg=defaultdict(set);byb=defaultdict(set);byj=defaultdict(set);bycs=defaultdict(list)
 for r in groups:
  for g in r['gene_ids']:byg[g].add(r['promoter_group_id'])
  for b in r['block_ids']:byb[b].add(r['promoter_group_id'])
 for gs in list(byg.values())+list(byb.values()):join(gs)
 for r in junctions:byj[tuple(r[k] for k in ['contig','strand','donor_boundary','acceptor_boundary'])].add(r['promoter_group_id'])
 for gs in byj.values():join(gs)
 for r in positions:bycs[(r['contig'],r['strand'])].append((r['position0'],r['promoter_group_id']))
 edge_count=0
 for vals in bycs.values():
  active=[]
  for pos,g in sorted(vals):
   active=[(p,h) for p,h in active if p>=pos-200]
   for p,h in active:
    if h!=g:union(g,h);edge_count+=1
   active.append((pos,g))
 components=defaultdict(list)
 for g in par:components[find(g)].append(g)
 rows=[]
 for gs in components.values():
  cid='WTC11_METADEP_'+hjson({'groups':sorted(gs)})[:24]
  rows.extend(dict(promoter_group_id=g,metadata_dependency_component_id=cid,component_group_count=len(gs)) for g in sorted(gs))
 tsv(out/'metadata_dependency_components.tsv',sorted(rows,key=lambda x:x['promoter_group_id']))
 warnings=[line.rstrip('\n') for line in (W/'logs/QC.stderr.log').read_text().splitlines() if 'WARNING' in line]
 summary=dict(status='PASS_WITH_DOCUMENTED_QC_WARNINGS_AND_GENERATOR_DOMAIN_LIMITATIONS',raw_transcripts=len(raw),
  human_QC_rows=len(sq_byid),nonhuman_control_exclusions=reasoncounts['nonhuman_contig'],
  technically_eligible_raw=sum(r['technical_eligibility'] for r in raw),technical_exclusion_reasons=dict(reasoncounts),
  raw_classification_counts=dict(classcounts),crosswalk_counts=dict(qcounts),
  corrected_coordinate_status_counts=dict(Counter(r['QC_coordinate_status'] for r in cross)),
  exact_structures=len(exact),candidate50_count=len(cand),candidate_global_holdout=sum(r['global_holdout'] for r in holds),
  candidate_exposure_excluded=sum(r['permanent_exposure_exclusion'] for r in holds),analysis_candidate_count=len(members),
  promoter_group_count=len(groups),group_TSS_positions=len(positions),group_junction_links=len(junctions),
  monoexon_only_groups=sum(r['junction_unique_count']==0 for r in groups),metadata_dependency_components=len(components),
  metadata_TSS_200bp_pair_edges=edge_count,metadata_component_largest=max(map(len,components.values())),
  window_graph_policy='SAME_CONTIG_STRAND_ALL_MEMBER_TSS_DISTANCE_LE200_FOR_BOOTSTRAP_ONLY;PROMOTER_GROUPING_UNCHANGED',
  warnings=warnings,coding_semantics='ORF_NOT_RUN_NONCODING_IS_UNMEASURED_PLACEHOLDER;NOT_USED',
  outcomes_read=0,new_models=0,discovery_rerun=0,old_assets_modified=False)
 js(out/'build_summary.json',summary)
 u1=resource.getrusage(resource.RUSAGE_SELF)
 js(out/'builder_resource.json',dict(wall_seconds=time.monotonic()-start,
  CPU_seconds=u1.ru_utime+u1.ru_stime-u.ru_utime-u.ru_stime,peak_RSS_bytes=u1.ru_maxrss*1024))
 outputs=[record(p,'STRUCTURE_OUTPUT') for p in sorted(out.iterdir()) if p.is_file()]
 js(out/'closed_output_inventory.json',dict(status='CLOSED_GENERATION_OUTPUTS',generation=args.generation,files=outputs,
  inventory_excludes_its_own_hash=True,analysis_proof='STRUCTURE_AND_COMMON10_6_ONLY_NOT_BAM_COMPLETION'))
 print(json.dumps(summary,sort_keys=True),flush=True)
if __name__=='__main__':main()
