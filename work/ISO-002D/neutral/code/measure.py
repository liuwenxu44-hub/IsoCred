"""ISO-002D full neutral measurement. Source-bound, no caller map or outcome inputs."""
from pathlib import Path
import collections,csv,datetime,gzip,hashlib,json,os,resource,subprocess,time,traceback
import numpy as np,pandas as pd,pyarrow as pa,pyarrow.parquet as pq,pysam
R=Path('.');N=R/'work/ISO-002D/neutral';SAM=R/'env/ISO-000C-R1/core/bin/samtools'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def js(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,sort_keys=True)
def identity(p):return dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=sha(p))
def put(p,frame):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);assert not p.exists()
 pq.write_table(pa.Table.from_pandas(frame,preserve_index=False),p,compression='zstd')
def readpq(p):return pq.ParquetFile(p).read().to_pandas()
def events(a):
 ts=a.get_tag('ts') if a.has_tag('ts') else 'MISSING'
 assert ts in ('+','-','MISSING'),'INVALID_TS'
 strand='-' if a.is_reverse else '+'
 if ts=='-':strand='+' if strand=='-' else '-'
 pos=a.reference_start;j=[]
 for op,length in a.cigartuples or []:
  if op==3:j.append((pos,pos+length) if strand=='+' else (pos+length,pos))
  if op in (0,2,3,7,8):pos+=length
 assert pos==a.reference_end,'CIGAR_REFERENCE_END'
 return (a.reference_id,a.reference_start,a.reference_end,strand),strand,j
def audit_and_filter(source,raw,rep,expected_raw_count,contigs,human,outdir):
 """Full source audit; temporary neutral BAM is not promoted until every gate passes."""
 outdir.mkdir(parents=True,exist_ok=False)
 idmap={}
 with gzip.open(raw,'rt') as f:
  while True:
   h=f.readline()
   if not h:break
   seq=f.readline();plus=f.readline();qual=f.readline()
   assert h.startswith('@') and plus.startswith('+') and qual.endswith('\n'),'RAW_FASTQ_FORMAT'
   assert len(seq.rstrip('\r\n'))==len(qual.rstrip('\r\n')),'RAW_FASTQ_LENGTH'
   rid=h[1:].split()[0];assert rid not in idmap,'RAW_DUPLICATE_ID';idmap[rid]=len(idmap)
 assert len(idmap)==expected_raw_count,'RAW_RECORD_COUNT'
 seen=bytearray(len(idmap));counts=collections.Counter();flags=collections.Counter();tags=collections.Counter();neutral_flags=collections.Counter()
 p=outdir/'neutral_unsorted.bam'
 with pysam.AlignmentFile(str(source),'rb') as bam:
  head=bam.header.to_dict()
  assert len(bam.references)==len(contigs) and len(set(bam.references))==len(contigs),'REF_COUNT'
  assert dict(zip(bam.references,bam.lengths))==contigs,'REF_DICT'
  assert any('-uf' in x.get('CL','') for x in head.get('PG',[])),'MISSING_UF'
  with pysam.AlignmentFile(str(p),'wb',header=bam.header) as dest:
   for a in bam.fetch(until_eof=True):
    counts['bam_rows']+=1;flags[str(a.flag)]+=1
    rid=a.query_name
    assert rid in idmap,('BAM_QNAME_OUTSIDE_EXACT_RAW',rid)
    seen[idmap[rid]]=1
    counts['secondary_rows']+=int(a.is_secondary);counts['supplementary_rows']+=int(a.is_supplementary)
    if a.is_unmapped:counts['unmapped_rows']+=1;continue
    counts['mapped_rows']+=1
    ts=a.get_tag('ts') if a.has_tag('ts') else 'MISSING';tags[ts]+=1
    assert ts in ('+','-','MISSING'),'INVALID_TS'
    if a.mapping_quality<20:counts['mapped_below_Q20']+=1;continue
    counts['mapped_Q20_rows']+=1
    if a.reference_name not in human:counts['Q20_nonhuman_rows']+=1;continue
    counts['neutral_rows']+=1;neutral_flags[str(a.flag)]+=1;dest.write(a)
 counts['raw_unique_QNAMEs']=len(idmap);counts['bam_unique_QNAMEs']=sum(seen);counts['raw_QNAMEs_missing_from_BAM']=len(seen)-sum(seen)
 v=dict(replicate=rep,counts=dict(counts),flags=dict(flags),neutral_flags=dict(neutral_flags),mapped_ts_tags=dict(tags),header=head,source_QNAME_set_equals_raw=(counts['raw_QNAMEs_missing_from_BAM']==0),readmap_used=False,supplementary_retained=True)
 js(outdir/'full_source_audit.json',v)
 assert counts['secondary_rows']==0,'UNEXPECTED_SECONDARY_REVIEW_REQUIRED'
 assert counts['raw_QNAMEs_missing_from_BAM']==0,'RAW_QNAME_ABSENT_IN_FULL_BAM'
 assert counts['neutral_rows']>0,'EMPTY_NEUTRAL'
 del idmap,seen
 return p,v
def extract(neutral,rep,outdir):
 outdir.mkdir(parents=True,exist_ok=False);temps=[];first={};amb=set();stats=collections.Counter();tags=collections.Counter()
 current=None;seen=set();lastkey=None;completed=set()
 def flush(contig):
  if not seen:return
  d=pd.DataFrame(sorted(seen),columns=['read_key','strand','donor_boundary','acceptor_boundary'])
  d.insert(0,'replicate',rep);d.insert(2,'contig',contig)
  p=outdir/(contig+'.pending.parquet');put(p,d);temps.append((contig,p))
 with pysam.AlignmentFile(str(neutral),'rb') as bam:
  for a in bam.fetch(until_eof=True):
   assert not a.is_unmapped and not a.is_secondary and a.mapping_quality>=20,'NEUTRAL_FILTER'
   order=(a.reference_id,a.reference_start)
   assert lastkey is None or order>=lastkey,'COORDINATE_ORDER';lastkey=order
   if a.reference_name!=current:
    if current is not None:flush(current);completed.add(current)
    assert a.reference_name not in completed,'CONTIG_REENTRY'
    current=a.reference_name;seen=set()
   stats['bam_rows']+=1;rid=rep+'|'+a.query_name
   span,strand,junctions=events(a)
   tags[a.get_tag('ts') if a.has_tag('ts') else 'MISSING']+=1
   if rid in first and first[rid]!=span:amb.add(rid)
   else:first[rid]=span
   for d,ac in junctions:
    k=(rid,strand,d,ac);stats['raw_junction_observations']+=1
    stats['duplicate_junction_observations_removed']+=int(k in seen);seen.add(k)
  if current is not None:flush(current)
 for counter_key in ['raw_junction_observations','unique_junction_observations','duplicate_junction_observations_removed']:
  stats.setdefault(counter_key,0)
 allcounts=[]
 for contig,p in temps:
  d=readpq(p);d['read_alignment_ambiguous']=d.read_key.isin(amb)
  dest=outdir/'observations'/contig/'part.parquet';put(dest,d)
  c=d.groupby(['contig','strand','donor_boundary','acceptor_boundary'],sort=True).agg(unique_read_support=('read_key','size'),ambiguous_read_support=('read_alignment_ambiguous','sum')).reset_index()
  c['replicate']=rep;allcounts.append(c);stats['unique_junction_observations']+=len(d)
  p.unlink()
 c=pd.concat(allcounts,ignore_index=True) if allcounts else pd.DataFrame(columns=['contig','strand','donor_boundary','acceptor_boundary','unique_read_support','ambiguous_read_support','replicate'])
 put(outdir/'junction_counts.parquet',c)
 put(outdir/'HUMAN_unique_read_ids.parquet',pd.DataFrame({'replicate':rep,'read_key':sorted(first),'read_alignment_ambiguous':[rid in amb for rid in sorted(first)]}))
 assert stats['raw_junction_observations']-stats['duplicate_junction_observations_removed']==stats['unique_junction_observations']
 result=dict(replicate=rep,**dict(stats),unique_mapped_read_ids=len(first),ambiguous_alignment_read_ids=len(amb),unique_junctions=len(c),ts_tag_counts=dict(tags),status='PASS',no_outcome_inputs=True)
 js(outdir/'denominators.json',result)
 return result
def main():
 S=json.loads((N/'specification.json').read_text());A=json.loads((N/'receipts/acquisition_complete.json').read_text());assert A['status']=='PASS'
 assert sha(R/S['contract']['path'])==S['contract']['sha256']
 tests=json.loads((N/'tests/synthetic_result.json').read_text());assert tests['status']=='PASS'
 runtime=json.loads((N/'receipts/runner_code_freeze.json').read_text())
 for frozen in runtime['files']:assert sha(R/frozen['path'])==frozen['sha256'],('RUNTIME_CODE_DRIFT',frozen['path'])
 codefiles=[N/'code/measure.py',N/'code/replay.py',N/'code/synthetic_tests.py',N/'code/supervise.py']
 inputs=[]
 for x in A['files']:
  assert sha(R/x['path'])==x['sha256'],'BAM_HASH_DRIFT'
  assert sha(R/x['raw_path'])==x['raw_sha256'],'FASTQ_HASH_DRIFT'
  inputs.extend([identity(R/x['path']),identity(R/x['raw_path'])])
 for p,h in S['inputs'].items():assert sha(R/p)==h,(p,'INPUT_HASH_DRIFT')
 js(N/'receipts/measurement_code_input_freeze.json',dict(before_full_BAM_content=True,contract=S['contract'],specification=identity(N/'specification.json'),inputs=inputs,code=[identity(p) for p in codefiles],synthetic_tests=identity(N/'tests/synthetic_result.json'),environment={'python':os.sys.version,'pysam':pysam.__version__,'pandas':pd.__version__,'pyarrow':pa.__version__,'numpy':np.__version__},readmap_used=False,no_outcome_inputs=True))
 rows=list(csv.DictReader((R/'manifests/ISO-000C-R2_reference_contigs.tsv').open(),delimiter='\t'))
 contigs={x['contig_name']:int(x['length']) for x in rows};human={x['contig_name'] for x in rows if x['role']=='HUMAN_PRIMARY_ANALYSIS_SEQUENCE'}
 assert len(contigs)==195 and len(human)==194
 totals=[]
 for x in A['files']:
  rep=x['replicate'];source=R/x['path'];repdir=N/'measurements'/rep
  start=time.monotonic()
  subprocess.run([str(SAM),'quickcheck','-v',str(source)],check=True)
  unsorted,audit=audit_and_filter(source,R/x['raw_path'],rep,x['raw_records'],contigs,human,N/'tmp'/rep)
  repdir.mkdir(exist_ok=False);neutral=repdir/'neutral.HUMAN_Q20.bam'
  command=[str(SAM),'sort','-@','3','-m','2G','-T',str(N/'tmp'/rep/'sort'),'-o',str(neutral),str(unsorted)]
  js(repdir/'sort_command.json',dict(command=command))
  subprocess.run(command,check=True)
  subprocess.run([str(SAM),'index','-@','3',str(neutral)],check=True)
  subprocess.run([str(SAM),'quickcheck','-v',str(neutral)],check=True)
  result=extract(neutral,rep,repdir/'junction')
  assert result['bam_rows']==audit['counts']['neutral_rows'],'SORT_FILTER_ROW_ACCOUNTING'
  result.update(neutral=identity(neutral),index=identity(Path(str(neutral)+'.bai')),whole_source_audit=identity(N/'tmp'/rep/'full_source_audit.json'),source_acquisition_sha256=x['sha256'],wall_seconds=time.monotonic()-start)
  js(repdir/'measurement_complete.json',result)
  subprocess.run([os.sys.executable,'-B',str(N/'code/replay.py'),rep],check=True)
  totals.append(result);print(json.dumps(result),flush=True)
 js(N/'receipts/measurement_complete.json',dict(status='PASS',replicates=totals,no_outcomes=True,processing_domain=S['processing_domain'],all_libraries_complete=True,zero_support_interpretation_authorized=True))
if __name__=='__main__':
 try:main()
 except BaseException as e:
  js(N/'receipts/measurement_failure.json',dict(status='FAIL',exception=repr(e),traceback=traceback.format_exc(),time=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  raise
