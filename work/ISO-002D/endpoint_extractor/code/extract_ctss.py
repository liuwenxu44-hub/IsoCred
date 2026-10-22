#!/usr/bin/env python3
"""Outcome-isolated ISO-002D endpoint extractor. Import is inert; real CLI requires P4 authorization."""
import argparse,bisect,csv,gzip,hashlib,io,json,os,re,sys,time,traceback
from collections import Counter,defaultdict
from datetime import datetime,timezone
from decimal import Decimal,localcontext,MAX_EMAX,MIN_EMIN
from pathlib import Path

ROOT=Path('.')
BASE=ROOT/'work/ISO-002D/endpoint_extractor'
SCOPE_SHA='360b541a032319be77829b8483b68290d4d47e6340c12ba6ed470442d8bf2289'
ENDPOINT_SHA='83561a30b925227efabaf9dbc9eb3ae4c0fd7baee810229f2b72782f1a90b5f9'
LIBS=('GSM5625727','GSM5625728','GSM5625729','GSM5625730')
NAMES=('GSM5625727_WTC11-1_rep1.ctss.bed.gz','GSM5625728_WTC11-1_rep2.ctss.bed.gz','GSM5625729_WTC11-1_rep3.ctss.bed.gz','GSM5625730_WTC11-3_rep1.ctss.bed.gz')
URLS={a:'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5625nnn/'+a+'/suppl/'+n for a,n in zip(LIBS,NAMES)}
ALLOWED={f'chr{i}' for i in range(1,23)}|{'chrX','chrY','chrM'}
MAPPING_FIELDS=['site_id','contig','strand','position0','compatible_group_count100','promoter_group_shared','bio1_recorded','bio3_recorded','bio1_pooled_signal','bio3_pooled_signal']+[a+'_signal' for a in LIBS]+['promoter_group_id','minimum_absolute_distance','within_primary50','dependency_component_id']
META_ROLES=('scope','endpoint','computation_contract','P1_acceptance','candidate_independence','model_manifest','prediction_manifest','prediction_replay','download_authorization','download_manifest','extractor_freeze')
ASSERTIONS=('candidate_independence_accepted','P1_complete','model_freeze_passed','all_predictions_frozen','fit_disabled_replay_passed','exact_four_downloads_completed','extractor_independently_accepted','no_further_fits','one_time_extraction','not_training_process')
STRUCT_PATHS={
 'groups':'work/ISO-002D/structure/build1/promoter_groups.tsv',
 'positions':'work/ISO-002D/structure/build1/group_exact_TSS_positions.tsv',
 'components':'output/ISO-002D/data/WTC11_dependency_components.parquet'}

class EndpointError(RuntimeError):pass
def require(ok,message):
 if not ok:raise EndpointError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def hjson(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def record(path,root=ROOT):
 path=Path(path);return dict(path=str(path.relative_to(root)),bytes=path.stat().st_size,sha256=sha(path))
def write_json(path,x):
 with Path(path).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def decimal_add(a,b):
 """Exact nonnegative decimal addition; no floating-point cancellation/underflow."""
 if not a:return b
 if not b:return a
 at,bt=a.as_tuple(),b.as_tuple();low=min(at.exponent,bt.exponent)
 with localcontext() as ctx:
  ctx.prec=max(len(at.digits)+at.exponent-low,len(bt.digits)+bt.exponent-low)+1
  ctx.Emax=MAX_EMAX;ctx.Emin=MIN_EMIN
  value=a+b
 require(value.is_finite() and value>=0,'NONFINITE_POOLED_SIGNAL')
 return value
def pooled(values):
 out=Decimal(0)
 for x in values:
  if x is not None:out=decimal_add(out,x)
 return out
def integer(value,name):
 require(isinstance(value,(str,int)) and not isinstance(value,bool) and re.fullmatch(r'[0-9]+',str(value)) is not None,name)
 return int(value)
def boolean(value,name):
 if isinstance(value,bool):return value
 require(value in ('true','false','True','False'),name)
 return value.lower()=='true'
def parse_bed(lines):
 """BED6 or wider; only valid data rows count as recorded sites. Headers are ignored."""
 for number,line in enumerate(lines,1):
  if not line.strip() or line.startswith('#') or line.startswith('track ') or line.startswith('browser '):continue
  z=line.rstrip('\r\n').split('\t')
  require(len(z)>=6,f'BED_FIELD_COUNT:{number}')
  require(bool(z[0]) and not any(c.isspace() for c in z[0]),f'BED_CONTIG:{number}')
  start=integer(z[1],f'BED_START:{number}');end=integer(z[2],f'BED_END:{number}')
  require(end==start+1,f'BED_NOT_ONE_BASE:{number}')
  require(z[5] in ('+','-'),f'BED_STRAND:{number}')
  try:signal=Decimal(z[4])
  except Exception as e:raise EndpointError(f'BED_SIGNAL:{number}') from e
  require(signal.is_finite() and signal>=0,f'BED_SIGNAL_NONFINITE_OR_NEGATIVE:{number}')
  yield z[0],z[5],start,signal
def pair_within(a,b,window=5):
 a=sorted(a);b=sorted(b);i=j=0
 while i<len(a) and j<len(b):
  delta=a[i]-b[j]
  if abs(delta)<=window:return True
  if delta<0:i+=1
  else:j+=1
 return False
def validate_groups(group_rows,position_rows,component_rows):
 groups={}
 for r in group_rows:
  gid=r['promoter_group_id'];require(isinstance(gid,str) and gid and gid not in groups,'DUPLICATE_OR_EMPTY_GROUP_ID')
  require(r['strand'] in ('+','-') and bool(r['contig']),'INVALID_GROUP_CONTEXT')
  require(not boolean(r['global_holdout'],'GLOBAL_FLAG') and not boolean(r['exposure_exclusion'],'EXPOSURE_FLAG'),'FORBIDDEN_GLOBAL_OR_EXPOSURE_GROUP')
  if 'technical_eligibility' in r:require(boolean(r['technical_eligibility'],'TECHNICAL_FLAG'),'TECHNICALLY_INELIGIBLE_GROUP')
  groups[gid]=dict(r)
 require(bool(groups),'EMPTY_STRUCTURAL_GROUP_UNIVERSE')
 positions=defaultdict(set);seen=set()
 for r in position_rows:
  gid=r['promoter_group_id'];require(gid in groups,'TSS_UNKNOWN_GROUP')
  require((r['contig'],r['strand'])==(groups[gid]['contig'],groups[gid]['strand']),'TSS_CONTEXT_MISMATCH')
  p=integer(r['position0'],'INVALID_TSS_0BASED');key=(gid,p)
  require(key not in seen,'DUPLICATE_GROUP_TSS');seen.add(key);positions[gid].add(p)
 require(set(positions)==set(groups),'GROUP_WITHOUT_ALLMEMBER_TSS')
 for gid,ps in positions.items():
  require(max(ps)-min(ps)<=50,'PROMOTER_TSS_ENVELOPE_EXCEEDS_50')
  g=groups[gid]
  for field,expected in [('exact_TSS_count',len(ps)),('TSS_min_0based',min(ps)),('TSS_max_0based',max(ps)),('TSS_span',max(ps)-min(ps))]:
   if field in g:require(integer(g[field],field)==expected,'TSS_METADATA_MISMATCH:'+field)
 components={};sizes={}
 for r in component_rows:
  gid=r['promoter_group_id'];cid=r['dependency_component_id']
  require(gid in groups and gid not in components and isinstance(cid,str) and bool(cid),'INVALID_COMPONENT_PARTITION')
  components[gid]=cid
  if 'component_group_count' in r:sizes[gid]=integer(r['component_group_count'],'COMPONENT_SIZE')
 require(set(components)==set(groups),'INCOMPLETE_FROZEN_COMPONENT_PARTITION')
 counts=Counter(components.values())
 for gid,size in sizes.items():require(size==counts[components[gid]],'COMPONENT_SIZE_MISMATCH')
 return groups,positions,components

class Extractor:
 def __init__(self,group_rows,position_rows,component_rows):
  self.groups,self.positions,self.components=validate_groups(group_rows,position_rows,component_rows)
  self.index=defaultdict(list)
  for gid,ps in self.positions.items():
   g=self.groups[gid]
   self.index[(g['contig'],g['strand'])].extend((p,gid) for p in ps)
  self.coords={}
  for key,arr in self.index.items():arr.sort();self.coords[key]=[p for p,g in arr]
  self.sites={};self.contigs={a:set() for a in LIBS};self.library_stats={a:Counter() for a in LIBS}
  self.fed=set()
 def hits(self,chrom,strand,p):
  key=(chrom,strand);coords=self.coords.get(key,[]);arr=self.index.get(key,[])
  lo=bisect.bisect_left(coords,p-100);hi=bisect.bisect_right(coords,p+100);hits={}
  for tss,gid in arr[lo:hi]:hits[gid]=min(hits.get(gid,101),abs(p-tss))
  return tuple(sorted(hits.items()))
 def feed(self,accession,lines):
  require(accession in LIBS and accession not in self.fed,'LIBRARY_DUPLICATE_OR_UNAUTHORIZED')
  self.fed.add(accession);lib=LIBS.index(accession);stats=self.library_stats[accession]
  for chrom,strand,p,signal in parse_bed(lines):
   stats['valid_rows']+=1;stats['zero_signal_rows']+=int(signal==0);self.contigs[accession].add(chrom)
   key=(chrom,strand,p)
   if key in self.sites:hits=self.sites[key]['hits']
   else:hits=self.hits(chrom,strand,p)
   if not hits:stats['unmapped_rows_within100']+=1;continue
   stats['mapped_rows_within100']+=1
   if len(hits)>1:
    require(len({self.components[g] for g,d in hits})==1,'SHARED_CTSS_CROSSES_FROZEN_COMPONENT:'+str(key))
   site=self.sites.setdefault(key,dict(hits=hits,signals=[None]*4))
   if site['signals'][lib] is not None:stats['duplicate_mapped_coordinate_rows']+=1
   site['signals'][lib]=decimal_add(site['signals'][lib] or Decimal(0),signal)
 def finish(self):
  require(self.fed==set(LIBS),'INCOMPLETE_EXACT_FOUR_LIBRARIES')
  states={g:dict(shared=0,unique=0,min_distance=None,pos=[set() for _ in LIBS],signals=[Decimal(0) for _ in LIBS],mapped=0) for g in self.groups}
  shared_sites=shared_pairs=0;pool_mapped_sites=[0,0]
  for key,site in self.sites.items():
   chrom,strand,p=key;hits=site['hits'];signals=site['signals'];shared=len(hits)>1
   shared_sites+=int(shared);shared_pairs+=len(hits)*(len(hits)-1)//2 if shared else 0
   pool_mapped_sites[0]+=int(any(x is not None for x in signals[:3]));pool_mapped_sites[1]+=int(signals[3] is not None)
   for gid,distance in hits:
    s=states[gid];s['shared' if shared else 'unique']+=1;s['mapped']+=1
    s['min_distance']=distance if s['min_distance'] is None else min(s['min_distance'],distance)
    if distance<=50:
     for k,value in enumerate(signals):
      if value is not None and value>0:s['pos'][k].add(p);s['signals'][k]=decimal_add(s['signals'][k],value)
  pool_contigs=[set.union(*(self.contigs[a] for a in LIBS[:3])),self.contigs[LIBS[3]]]
  common=pool_contigs[0]&pool_contigs[1];rows=[];eligible=[]
  for gid in sorted(self.groups):
   g=self.groups[gid];s=states[gid];bio1=set.union(*s['pos'][:3]);bio3=s['pos'][3]
   primary=bool(bio1 and bio3);concordant=pair_within(bio1,bio3,5)
   require(not concordant or primary,'SECONDARY_NOT_PRIMARY_SUBSET')
   ambiguity='ALL_COMPATIBLE_CTSS_PROMOTER_GROUP_SHARED' if s['shared']>0 and s['unique']==0 else 'PARTLY_PROMOTER_GROUP_SHARED' if s['shared'] else 'NO_PROMOTER_GROUP_SHARING'
   reasons=[]
   if g['contig'] not in ALLOWED:reasons.append('CONTIG_NOT_ALLOWED')
   if g['contig'] not in pool_contigs[0]:reasons.append('CONTIG_ABSENT_BIO1_RECORDED_ROWS')
   if g['contig'] not in pool_contigs[1]:reasons.append('CONTIG_ABSENT_BIO3_RECORDED_ROWS')
   if ambiguity=='ALL_COMPATIBLE_CTSS_PROMOTER_GROUP_SHARED':reasons.append(ambiguity)
   ok=not reasons
   if ok:eligible.append(gid)
   row=dict(promoter_group_id=gid,contig=g['contig'],strand=g['strand'],dependency_component_id=self.components[gid],
    exact_TSS_count=len(self.positions[gid]),TSS_min_0based=min(self.positions[gid]),TSS_max_0based=max(self.positions[gid]),
    assay_eligible=ok,exclusion_reasons=reasons,mapping_ambiguity=ambiguity,shared_compatible_sites100=s['shared'],unique_compatible_sites100=s['unique'],
    recorded_compatible_sites100=s['mapped'],minimum_absolute_distance=s['min_distance'],
    bio1_detected50=bool(bio1),bio3_detected50=bool(bio3),biological_detection_count=int(bool(bio1))+int(bool(bio3)),
    bio1_signal_sum50=str(pooled(s['signals'][:3])),bio3_signal_sum50=str(s['signals'][3]),
    primary_detected_2of2=primary,primary_label=int(primary) if ok else None,
    secondary_position_concordant_5bp=concordant,secondary_5bp_label=int(concordant) if ok else None,
    outcome_status='ASSAY_INELIGIBLE' if not ok else 'REPRODUCIBLE_PROMOTER_CTSS_SUPPORT' if primary else 'ASSAY_ELIGIBLE_NOT_REPRODUCIBLY_DETECTED')
   for k,a in enumerate(LIBS[:3]):
    positive=bool(s['pos'][k] and bio3)
    row[a+'_vs_GSM5625730_detected']=positive
    row[a+'_vs_GSM5625730_eligible']=ok
    row[a+'_vs_GSM5625730_label']=int(positive) if ok else None
   rows.append(row)
  masksha=hjson(eligible)
  summary=dict(status='PASS',unit='promoter_group_id',groups=len(rows),eligible_groups=len(eligible),
   primary_positive_eligible=sum(r['primary_label']==1 for r in rows),secondary_5bp_positive_eligible=sum(r['secondary_5bp_label']==1 for r in rows),
   exclusion_reason_counts=dict(Counter(z for r in rows for z in r['exclusion_reasons'])),
   primary_mask_sha256=masksha,technical_pair_mask_sha256={a+'_vs_GSM5625730':masksha for a in LIBS[:3]},
   mapped_unique_coordinate_strand_sites100=len(self.sites),pooled_mapped_unique_sites100=dict(bio1=pool_mapped_sites[0],bio3=pool_mapped_sites[1]),
   shared_recorded_sites100=shared_sites,actual_shared_pair_site_checks=shared_pairs,shared_dependency_violations=0,
   recorded_contigs_by_library={a:sorted(self.contigs[a]) for a in LIBS},
   recorded_contigs_by_pool=dict(bio1=sorted(pool_contigs[0]),bio3=sorted(pool_contigs[1])),common_recorded_contigs=sorted(common),
   library_row_statistics={a:dict(self.library_stats[a]) for a in LIBS},mapping_unit='unique recorded coordinate/strand over all four libraries; rows outside100 counted separately',
   biological_pools=2,donors=1,technical_pairs_all_reported=True,all_pairs_same_primary_mask=True,
   recorded_zero_sites_in_ambiguity=True,secondary5_is_subset_primary=True,candidate_label_expansion=False,
   interpretation='CTSS promoter support event; eligible negative is not biological absence or transcript nonfunction')
  return rows,summary
 def mapping_rows(self):
  for (chrom,strand,p),site in sorted(self.sites.items()):
   signals=site['signals'];hits=site['hits'];base=dict(site_id='CTSS_'+hjson([chrom,strand,p])[:24],contig=chrom,strand=strand,position0=p,
    compatible_group_count100=len(hits),promoter_group_shared=len(hits)>1,bio1_recorded=any(v is not None for v in signals[:3]),bio3_recorded=signals[3] is not None,
    bio1_pooled_signal=str(pooled(signals[:3])),bio3_pooled_signal=str(signals[3]) if signals[3] is not None else None)
   for k,a in enumerate(LIBS):base[a+'_signal']=str(signals[k]) if signals[k] is not None else None
   for gid,d in hits:yield dict(base,promoter_group_id=gid,minimum_absolute_distance=d,within_primary50=d<=50,dependency_component_id=self.components[gid])

def project_path(text,root=ROOT):
 require(isinstance(text,str) and text and not Path(text).is_absolute() and '..' not in Path(text).parts,'NONPROJECT_RELATIVE_PATH')
 path=root/text
 require(path.resolve().is_relative_to(root.resolve()),'PROJECT_PATH_ESCAPE')
 cursor=path
 while cursor!=root:
  require(not cursor.is_symlink(),'SYMLINK_INPUT_OR_OUTPUT')
  cursor=cursor.parent
 return path
def verify_record(binding,root=ROOT):
 require(isinstance(binding,dict) and set(('path','bytes','sha256'))<=set(binding),'INCOMPLETE_FILE_BINDING')
 require(re.fullmatch(r'[0-9a-f]{64}',binding['sha256']) is not None and isinstance(binding['bytes'],int) and binding['bytes']>=0,'INVALID_FILE_BINDING')
 path=project_path(binding['path'],root)
 require(path.is_file() and path.stat().st_size==binding['bytes'],'BOUND_FILE_SIZE:'+binding['path'])
 require(sha(path)==binding['sha256'],'BOUND_FILE_SHA:'+binding['path'])
 return path
def validate_unlock_shape(u):
 require(u.get('schema_version')=='1.0' and u.get('type')=='ISO002D_P4_INDEPENDENT_CTSS_EXTRACTION_UNLOCK','UNLOCK_TYPE')
 require(u.get('status')=='PASS' and u.get('authorized_action')=='ONE_TIME_CTSS_EXTRACTION','UNLOCK_NOT_AUTHORIZED')
 require(u.get('independent_authority') is True and isinstance(u.get('auditor_id'),str) and bool(u['auditor_id']),'INDEPENDENT_UNLOCK_AUTHORITY_REQUIRED')
 require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}',u.get('unlock_id','')) is not None,'INVALID_UNLOCK_ID')
 require(set(u.get('bindings',{}))==set(META_ROLES),'UNLOCK_METADATA_BINDINGS_INCOMPLETE')
 require(set(u.get('structural_inputs',{}))==set(STRUCT_PATHS),'UNLOCK_STRUCTURAL_BINDINGS_INCOMPLETE')
 require(all(u.get('assertions',{}).get(k) is True for k in ASSERTIONS),'UNLOCK_GATE_ASSERTION_NOT_PASSED')
 require(u['bindings']['scope']['sha256']==SCOPE_SHA and u['bindings']['endpoint']['sha256']==ENDPOINT_SHA,'SCOPE_ENDPOINT_DRIFT')
 require(u['bindings']['scope']['path']=='contracts/ISO-002D_common10_external_scope.yaml' and u['bindings']['endpoint']['path']=='work/ISO-002D/endpoint/endpoint_contract.json','WRONG_CONTRACT_PATH')
 for role,path in STRUCT_PATHS.items():require(u['structural_inputs'][role]['path']==path,'STRUCTURAL_ROLE_PATH:'+role)
 files=u.get('ctss_inputs',[])
 require([x.get('accession') for x in files]==list(LIBS),'EXACT_FOUR_ORDERED_CTSS_REQUIRED')
 for x,name in zip(files,NAMES):
  require(x.get('url')==URLS[x['accession']] and Path(x.get('path','')).name==name,'CTSS_SOURCE_NOT_ALLOWLISTED')
  require(isinstance(x.get('bytes'),int) and x['bytes']>0 and re.fullmatch(r'[0-9a-f]{64}',x.get('sha256','')) is not None,'CTSS_EXACT_CHECKSUM_REQUIRED')
 require(u.get('output_directory')=='work/ISO-002D/endpoint_extractor/runs/'+u['unlock_id'],'OUTPUT_DIRECTORY_NOT_EXCLUSIVE')
 return u
def prepare_authorized_run(unlock_path,expected_sha):
 """Validate independent receipt/metadata first; no CTSS bytes are opened here."""
 require(re.fullmatch(r'[0-9a-f]{64}',expected_sha) is not None,'UNLOCK_EXTERNAL_SHA_REQUIRED')
 path=project_path(unlock_path);require(sha(path)==expected_sha,'UNLOCK_RECEIPT_SHA')
 u=validate_unlock_shape(json.loads(path.read_text()));bound={}
 for role,binding in u['bindings'].items():
  require(Path(binding['path']).suffix.lower() in ('.json','.yaml','.yml'),'METADATA_ONLY_BINDING_REQUIRED:'+role)
  bound[role]=verify_record(binding)
 for role,binding in u['structural_inputs'].items():bound[role]=verify_record(binding)
 endpoint=json.loads(bound['endpoint'].read_text())
 require(endpoint['primary']['window_bp']==50 and endpoint['eligibility']['mapping_enumeration_bp']==100,'ENDPOINT_WINDOW_DRIFT')
 require(endpoint['primary']['bio1']==list(LIBS[:3]) and endpoint['primary']['bio3']==list(LIBS[3:]),'BIOLOGICAL_POOL_DRIFT')
 freeze=json.loads(bound['extractor_freeze'].read_text())
 code_entries=[f for f in freeze['files'] if f['path']==str(Path(__file__).relative_to(ROOT))]
 require(len(code_entries)==1,'EXTRACTOR_CODE_NOT_FROZEN');verify_record(code_entries[0])
 download=json.loads(bound['download_manifest'].read_text())
 require(download.get('status')=='PASS' and download.get('files')==u['ctss_inputs'],'DOWNLOAD_MANIFEST_BINDING_MISMATCH')
 # A single fixed lock across all unlock IDs prevents repeated real extraction.
 out=project_path(u['output_directory'])
 require(not out.exists(),'OUTPUT_ALREADY_EXISTS')
 require(not (BASE/'receipts/P4_consumption.json').exists(),'P4_ALREADY_CONSUMED_NO_AUTOMATIC_RETRY')
 return u,bound,out
def table_rows(path):
 if Path(path).suffix=='.parquet':
  import pyarrow.parquet as pq
  for batch in pq.ParquetFile(path).iter_batches(batch_size=65536,use_threads=False):
   yield from batch.to_pylist()
 else:
  with Path(path).open(newline='') as f:yield from csv.DictReader(f,delimiter='\t')
def tsv(path,rows,fields=None):
 iterator=iter(rows);first=next(iterator,None)
 if fields is None:
  require(first is not None,'EMPTY_OUTPUT_REQUIRES_DECLARED_SCHEMA')
  fields=list(first)
 raw=Path(path).open('xb')
 zipped=str(path).endswith('.gz')
 gz=gzip.GzipFile(fileobj=raw,filename='',mode='wb',mtime=0) if zipped else raw
 text=io.TextIOWrapper(gz,encoding='utf-8',newline='')
 def enc(v):
  if v is None:return 'NA'
  if isinstance(v,bool):return str(v).lower()
  if isinstance(v,(list,dict)):return json.dumps(v,sort_keys=True,separators=(',',':'))
  return v
 try:
  writer=csv.DictWriter(text,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  if first is not None:writer.writerow({k:enc(v) for k,v in first.items()})
  for row in iterator:writer.writerow({k:enc(v) for k,v in row.items()})
 finally:
  text.close()
  if not raw.closed:raw.close()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--unlock-receipt',required=True,help='Project-relative independent P4 authorization JSON')
 p.add_argument('--unlock-sha256',required=True,help='Externally frozen receipt SHA256; no inference or autoauthorization')
 args=p.parse_args();u,bound,out=prepare_authorized_run(args.unlock_receipt,args.unlock_sha256)
 started=time.monotonic()
 write_json(BASE/'receipts/P4_consumption.json',dict(status='CONSUMED_BEFORE_CTSS_READ_NO_AUTOMATIC_RETRY',
  unlock_id=u['unlock_id'],unlock_path=args.unlock_receipt,unlock_sha256=args.unlock_sha256,started_utc=datetime.now(timezone.utc).isoformat(),
  extractor_sha256=sha(Path(__file__)),no_fit=True,no_network=True))
 try:
  out.parent.mkdir(exist_ok=True);out.mkdir(exist_ok=False)
  write_json(out/'execution_binding.json',dict(unlock_sha256=args.unlock_sha256,bindings=u['bindings'],structural_inputs=u['structural_inputs'],ctss_inputs=u['ctss_inputs']))
  engine=Extractor(table_rows(bound['groups']),table_rows(bound['positions']),table_rows(bound['components']))
  for f in u['ctss_inputs']:
   path=verify_record(f)
   with gzip.open(path,'rt',encoding='utf-8',newline='') as lines:engine.feed(f['accession'],lines)
  rows,summary=engine.finish()
  tsv(out/'promoter_group_endpoints.tsv',rows)
  tsv(out/'CTSS_to_group_links.tsv.gz',engine.mapping_rows(),MAPPING_FIELDS)
  write_json(out/'mapping_summary.json',summary)
  result=dict(status='PASS_EXTRACTION_ONLY_NO_FITS',unlock_sha256=args.unlock_sha256,scope_sha256=SCOPE_SHA,endpoint_sha256=ENDPOINT_SHA,
   group_label_copies_to_candidates=0,model_reads=0,model_fits=0,network_calls=0,real_CTSS_files_read=4,
   wall_seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat(),
   files=[record(out/name) for name in ['execution_binding.json','promoter_group_endpoints.tsv','CTSS_to_group_links.tsv.gz','mapping_summary.json']])
  write_json(out/'extraction_complete.json',result);print(json.dumps(result,sort_keys=True))
 except BaseException as error:
  failure=dict(status='FAIL_TARGET_EXTRACTION_CONSUMED_STOP',exception=repr(error),traceback=traceback.format_exc(),
   unlock_sha256=args.unlock_sha256,automatic_retry_forbidden=True,wall_seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat())
  write_json(BASE/'receipts/P4_extraction_failure.json',failure);raise
if __name__=='__main__':main()
