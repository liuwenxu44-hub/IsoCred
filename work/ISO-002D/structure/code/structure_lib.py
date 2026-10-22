"""Pure fixed-structure routines. Original functions loaded by named AST definitions only."""
import ast,builtins,csv,gzip,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path('.')
def functions(path,names,extra=None):
 source=Path(path).read_text();tree=ast.parse(source)
 nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in names]
 assert {n.name for n in nodes}==set(names)
 ns=dict(hashlib=hashlib,json=json,defaultdict=defaultdict,Counter=Counter,csv=csv,re=re,Path=Path)
 if extra:ns.update(extra)
 exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns)
 return ns
ID=functions(ROOT/'scripts/iso000d_r2b_identity.py',
 ['canonical_structure','reciprocal_overlap','can_add_complete_link','deterministic_clusters','cluster_sha256'])
PROM=functions(ROOT/'scripts/iso000f_build_promoter_group_outcomes.py',['stable_id','promoter_groups'])
canonical_structure=ID['canonical_structure']
def hjson(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def text_open(path,*args,**kwargs):
 if str(path).endswith('.gz'):
  mode=args[0] if args else kwargs.pop('mode','rt')
  if mode=='r':mode='rt'
  return gzip.open(path,mode,**kwargs)
 return builtins.open(path,*args,**kwargs)
def raw_namespace():
 return functions(ROOT/'scripts/iso000d_r2b_r1_build_raw_and_crosswalk.py',
 ['attrs','parse_gtf','normalize_gene','gencode_relation','parse_sqanti'],
 dict(canonical_structure=canonical_structure,open=text_open))
def reference_index(ns,path):
 ref,malformed=ns['parse_gtf'](path);assert not malformed
 genes={};gene_by_base={};gene_by_symbol=defaultdict(set);ref_exact=defaultdict(set)
 ref_chain=defaultdict(set);ref_junctions=defaultdict(set)
 for tid,t in ref.items():
  if len(t['chroms'])!=1 or len(t['strands'])!=1:continue
  chrom=next(iter(t['chroms']));strand=next(iter(t['strands']))
  try:c=canonical_structure(chrom,strand,t['exons'])
  except ValueError:continue
  gid=t['attrs'].get('gene_id','');gname=t['attrs'].get('gene_name','NA')
  if not gid:continue
  base=gid.split('.')[0];gene_by_base[base]=gid;gene_by_symbol[gname].add(gid)
  lo=min(a for a,b in c['ordered_exons']);hi=max(b for a,b in c['ordered_exons'])
  z=genes.setdefault(gid,dict(gene_id=gid,gene_name=gname,chromosome=chrom,strand=strand,start=lo,end=hi))
  z['start']=min(z['start'],lo);z['end']=max(z['end'],hi)
  ex=tuple(c['ordered_exons']);j=tuple(c['ordered_internal_junctions'])
  ref_exact[(chrom,strand,ex,c['TSS'],c['TES'])].add(gid)
  ref_chain[(chrom,strand,j)].add((c['TSS'],c['TES'],gid))
  ref_junctions[(chrom,strand)].update(j)
 gene_intervals=defaultdict(list)
 for g in genes.values():gene_intervals[(g['chromosome'],g['strand'])].append(g)
 for k in gene_intervals:gene_intervals[k].sort(key=lambda x:(x['start'],x['end'],x['gene_id']))
 ns.update(genes=genes,gene_by_base=gene_by_base,gene_by_symbol=gene_by_symbol,ref_exact=ref_exact,
  ref_chain=ref_chain,ref_junctions=ref_junctions,gene_intervals=gene_intervals)
 return ns
def exacts_from_raw(raw):
 byhash=defaultdict(list)
 for r in raw:
  if r['technical_eligibility']:byhash[r['exact_structure_hash']].append(r)
 exact=[]
 for h,m in sorted(byhash.items()):
  first=m[0];compat=sorted({g for r in m for g in r['compatible_gene_ids']})
  statuses={r['gene_assignment_status'] for r in m}
  if len(compat)==1 and 'AMBIGUOUS' not in statuses:
   locus=compat[0].split('.')[0];primary=compat[0];status='KNOWN_UNIQUE'
  elif compat:
   locus='AMBIGUOUS_LOCUS_'+hjson({'compatible_gene_ids':compat})[:24];primary='NA';status='AMBIGUOUS'
  else:locus='PENDING_NOVEL';primary='NA';status='NOVEL'
  exact.append(dict(exact_structure_hash=h,member_raw_model_keys=sorted(r['raw_model_key'] for r in m),
   chromosome=first['chromosome'],strand=first['strand'],ordered_exon_chain=first['ordered_exon_chain'],
   junction_chain=json.dumps(first['junction_chain'],separators=(',',':')),TSS=first['TSS'],TES=first['TES'],
   exon_count=first['exon_count'],genomic_span=first['genomic_span'],locus_group_id=locus,primary_gene_id=primary,
   all_compatible_gene_ids=compat,gene_assignment_status=status,
   GENCODE_structural_relations=sorted({r['GENCODE_structural_relation'] for r in m}),
   SQANTI3_classes=sorted({r['SQANTI3_structural_category'] for r in m}),monoexon=first['exon_count']==1))
 for x in exact:
  if x['locus_group_id']=='PENDING_NOVEL' and x['exon_count']>1:
   x['locus_group_id']='NOVEL_LOCUS_'+x['chromosome']+'_'+x['strand']+'_'+hjson(
    {'contig':x['chromosome'],'strand':x['strand'],'junction_chain':json.loads(x['junction_chain'])})[:24]
 pending=defaultdict(list)
 for x in exact:
  if x['locus_group_id']=='PENDING_NOVEL':
   a,b=x['genomic_span'];pending[(x['chromosome'],x['strand'])].append((a,b,x))
 for key,vals in pending.items():
  vals.sort(key=lambda z:(z[0],z[1],z[2]['exact_structure_hash']));groups=[];cur=[];end=-1
  for a,b,x in vals:
   if cur and a>end:groups.append(cur);cur=[];end=-1
   cur.append((a,b,x));end=max(end,b)
  if cur:groups.append(cur)
  for group in groups:
   lo=min(a for a,b,x in group);hi=max(b for a,b,x in group)
   lid='NOVEL_LOCUS_'+key[0]+'_'+key[1]+'_'+hjson({'contig':key[0],'strand':key[1],'interval':[lo,hi]})[:24]
   for a,b,x in group:x['locus_group_id']=lid
 return exact
def candidates_from_exacts(exact):
 clusters=ID['deterministic_clusters'](exact,50,.8);rows=[];xwalk=[]
 for c in sorted(clusters,key=lambda z:(z[0]['chromosome'],z[0]['strand'],min(x['exact_structure_hash'] for x in z))):
  hs=sorted(x['exact_structure_hash'] for x in c);own=ID['cluster_sha256'](50,c);cid='ISOCRED50_'+own
  compat=sorted({g for x in c for g in x['all_compatible_gene_ids']})
  prim=sorted({x['primary_gene_id'] for x in c if x['primary_gene_id']!='NA'});statuses={x['gene_assignment_status'] for x in c}
  status='AMBIGUOUS' if len(prim)>1 or 'AMBIGUOUS' in statuses else ('KNOWN_UNIQUE' if len(prim)==1 else 'NOVEL')
  loci=sorted({x['locus_group_id'] for x in c});assert len(loci)==1
  rows.append(dict(candidate_id=cid,cluster_sha256=own,tolerance_bp=50,member_exact_hashes=hs,
   member_raw_model_keys=sorted({r for x in c for r in x['member_raw_model_keys']}),
   representative_exact_structure_hash=min(hs),representative_is_file_encoding_only=True,
   contig=c[0]['chromosome'],strand=c[0]['strand'],junction_chain=c[0]['junction_chain'],
   endpoint_envelope=json.dumps({'TSS':[min(x['TSS'] for x in c),max(x['TSS'] for x in c)],
     'TES':[min(x['TES'] for x in c),max(x['TES'] for x in c)]},sort_keys=True,separators=(',',':')),
   locus_group=loci[0],primary_gene_id=prim[0] if len(prim)==1 else 'NA',
   all_compatible_gene_ids=compat,gene_assignment_status=status,monoexon=all(x['monoexon'] for x in c),
   GENCODE_structural_relations=sorted({z for x in c for z in x['GENCODE_structural_relations']}),
   candidate_generator='ARCHIVAL_TALON4_LAPA_LONG_ONLY_POOL3',technical_eligibility=True))
  xwalk.extend(dict(exact_structure_hash=h,candidate_id=cid) for h in hs)
 return rows,xwalk
def original_holdout_blocks(ref,hold,attrs):
 assert len(hold)==850 and sum(r['exposure_exclusion'].lower()=='true' for r in hold)==11
 holdids={r['gene_id'].split('.')[0] for r in hold};exids={r['gene_id'].split('.')[0] for r in hold if r['exposure_exclusion'].lower()=='true'}
 holdsyms={r['gene_symbol'] for r in hold};exsyms={r['gene_symbol'] for r in hold if r['exposure_exclusion'].lower()=='true'}
 gene={}
 with open(ref) as f:
  for line in f:
   if line.startswith('#'):continue
   z=line.rstrip('\n').split('\t')
   if len(z)!=9 or z[2]!='gene':continue
   a=attrs(z[8]);gid=a.get('gene_id')
   if gid:gene[gid.split('.')[0]]=dict(gene_id=gid,symbol=a.get('gene_name','NA'),chrom=z[0],strand=z[6],start=int(z[3]),end=int(z[4]))
 components=[];bycs=defaultdict(list)
 for g in gene.values():bycs[(g['chrom'],g['strand'])].append(g)
 for key,vals in sorted(bycs.items()):
  vals.sort(key=lambda g:(g['start'],g['end'],g['gene_id']));cur=[];end=-1
  for g in vals:
   a=max(1,g['start']-10000);b=g['end']+10000
   if cur and a>end:
    components.append(dict(chrom=key[0],strand=key[1],start=min(max(1,x['start']-10000) for x in cur),end=max(x['end']+10000 for x in cur),genes=cur));cur=[];end=-1
   cur.append(g);end=max(end,b)
  if cur:components.append(dict(chrom=key[0],strand=key[1],start=min(max(1,x['start']-10000) for x in cur),end=max(x['end']+10000 for x in cur),genes=cur))
 par=list(range(len(components)))
 def find(x):
  while par[x]!=x:par[x]=par[par[x]];x=par[x]
  return x
 def union(a,b):
  a,b=find(a),find(b)
  if a!=b:par[max(a,b)]=min(a,b)
 bych=defaultdict(list)
 for i,c in enumerate(components):bych[c['chrom']].append(i)
 for chrom,inds in bych.items():
  inds.sort(key=lambda i:(components[i]['start'],components[i]['end'],components[i]['strand']));active=[]
  for i in inds:
   active=[j for j in active if components[j]['end']>=components[i]['start']]
   for j in active:
    if components[j]['strand']!=components[i]['strand'] and components[j]['start']<=components[i]['end']:union(i,j)
   active.append(i)
 blocks=defaultdict(list)
 for i,c in enumerate(components):blocks[find(i)].append(c)
 final=[]
 for group in blocks.values():
  gs=sorted({g['gene_id'].split('.')[0]:g for c in group for g in c['genes']}.values(),key=lambda x:x['gene_id'])
  chrom=group[0]['chrom'];start=min(c['start'] for c in group);end=max(c['end'] for c in group)
  bid='GENOMIC_BLOCK_'+chrom+'_'+hjson({'chrom':chrom,'start':start,'end':end,'genes':[g['gene_id'].split('.')[0] for g in gs]})[:24]
  final.append(dict(block_id=bid,chrom=chrom,start=start,end=end,
   genes=[g['gene_id'] for g in gs],holdout=any(g['gene_id'].split('.')[0] in holdids or g['symbol'] in holdsyms for g in gs)))
 return dict(gene=gene,blocks=final,holdids=holdids,exids=exids,holdsyms=holdsyms,exsyms=exsyms)
def candidate_holdouts(candidates,raw,ctx):
 byraw={x['raw_model_key']:x for x in raw};bych=defaultdict(list)
 for b in ctx['blocks']:bych[b['chrom']].append(b)
 out=[]
 for x in candidates:
  bases={g.split('.')[0] for g in x['all_compatible_gene_ids']};syms={ctx['gene'].get(b,{}).get('symbol','NA') for b in bases}
  direct=sorted((bases&ctx['holdids'])|{s for s in syms if s in ctx['holdsyms']})
  exposure=sorted((bases&ctx['exids'])|{s for s in syms if s in ctx['exsyms']})
  spans=[byraw[k]['genomic_span'] for k in x['member_raw_model_keys']];lo=min(s[0] for s in spans);hi=max(s[1] for s in spans)
  overlap=sorted({b['block_id'] for b in bych[x['contig']] if b['start']<=hi and b['end']>=lo})
  holdblocks=sorted({b['block_id'] for b in bych[x['contig']] if b['holdout'] and b['start']<=hi and b['end']>=lo})
  held=bool(direct or holdblocks)
  out.append(dict(candidate_id=x['candidate_id'],all_compatible_gene_ids=x['all_compatible_gene_ids'],
   compatible_gene_symbols=sorted(syms),gene_assignment_status=x['gene_assignment_status'],matched_holdout_entries=direct,
   direct_gene_holdout=bool(direct),overlapping_frozen_genomic_blocks=overlap,matched_holdout_genomic_blocks=holdblocks,
   regional_isolation_holdout=bool(holdblocks),global_holdout=held,exposure_exclusion_entries=exposure,
   permanent_exposure_exclusion=bool(exposure),analysis_eligible=not held and not exposure,
   candidate_retained_in_universe=True))
 return out
def opportunity(positions):
 a=sorted(set(positions));return 101+sum(min(101,y-x) for x,y in zip(a,a[1:])) if a else 0
def junction_keys(candidate):
 pairs=json.loads(candidate['junction_chain']);s=candidate['strand']
 return [(candidate['contig'],s,int(a) if s=='+' else int(a)-1,int(b)-1 if s=='+' else int(b)) for a,b in pairs]
def groups_and_features(candidates,exact,xwalk,holds):
 hby={r['candidate_id']:r for r in holds};cby={c['candidate_id']:c for c in candidates}
 groups=PROM['promoter_groups']([c for c in candidates if hby[c['candidate_id']]['analysis_eligible']])
 c2g={c:g for g,cs in groups.items() for c in cs};e2c={r['exact_structure_hash']:r['candidate_id'] for r in xwalk}
 pos=defaultdict(set);hashes=defaultdict(set)
 for e in exact:
  g=c2g.get(e2c[e['exact_structure_hash']])
  if g:
   pos[g].add(e['TSS']-1);hashes[g].add(e['exact_structure_hash'])
 rows=[];positions=[];junctions=[];members=[]
 for gid in sorted(groups):
  cs=sorted(groups[gid]);first=cby[cs[0]];jkeys=set();scores=[];genes=set();blocks=set()
  for cid in cs:
   c=cby[cid];jkeys.update(junction_keys(c));rel=set(c['GENCODE_structural_relations'])
   scores.append(2 if 'KNOWN_EXACT' in rel else 1 if 'KNOWN_INTRON_CHAIN' in rel else 0)
   genes.update(g.split('.')[0] for g in c['all_compatible_gene_ids']);blocks.update(hby[cid]['overlapping_frozen_genomic_blocks'])
   members.append(dict(promoter_group_id=gid,candidate_id=cid))
  ps=sorted(pos[gid]);assert ps and max(ps)-min(ps)<=50
  rows.append(dict(promoter_group_id=gid,contig=first['contig'],strand=first['strand'],locus_group=first['locus_group'],
   candidate_count=len(cs),exact_structure_count=len(hashes[gid]),exact_TSS_count=len(ps),TSS_min_0based=min(ps),TSS_max_0based=max(ps),
   TSS_span=max(ps)-min(ps),window_opportunity_bp=opportunity(ps),
   monoexon_fraction=sum(cby[c]['monoexon'] for c in cs)/len(cs),
   B2_annotation_fixed_group_mean=sum(scores)/len(scores),junction_unique_count=len(jkeys),
   gene_ids=sorted(genes),block_ids=sorted(blocks),global_holdout=False,exposure_exclusion=False,
   study_domain='WTC11_ARCHIVAL_TALON4_LAPA_POOL3_BAM_BIO1_BIO3'))
  positions.extend(dict(promoter_group_id=gid,contig=first['contig'],strand=first['strand'],position0=p) for p in ps)
  junctions.extend(dict(promoter_group_id=gid,contig=k[0],strand=k[1],donor_boundary=k[2],acceptor_boundary=k[3]) for k in sorted(jkeys))
 return rows,positions,junctions,members
