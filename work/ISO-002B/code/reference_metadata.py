"""Small primary-reference metadata inventory only; no manuscript text or new study data."""
import json,hashlib,datetime,urllib.request,urllib.parse,time
from pathlib import Path
W=Path('./work/ISO-002B');D=W/'metadata/reference_identity';assert not D.exists();D.mkdir()
items={
'FLAIR':'10.1038/s41467-020-15171-6','IsoQuant':'10.1038/s41587-022-01565-y','SQANTI3':'10.1038/s41592-024-02229-2','SQANTI':'10.1101/gr.222976.117','TALON_preprint':'10.1101/672931','PSQAN':'10.1093/bioadv/vbaf293','LRGASP':'10.1038/s41592-024-02298-3','ENCODE2020':'10.1038/s41586-020-2493-4','ENCODEportal':'10.1093/nar/gkv1160','GENCODE':'10.1093/nar/gky955','FANTOM_promoter_atlas':'10.1038/nature13182','CAGEr':'10.1093/nar/gkv054','minimap2':'10.1093/bioinformatics/bty191','SAMtools':'10.1093/bioinformatics/btp352','RandomForest':'10.1023/A:1010933404324','ElasticNet':'10.1111/j.1467-9868.2005.00503.x','NestedCV_bias':'10.1186/1471-2105-7-91','Bootstrap':'10.1214/aos/1176344552','HierarchicalCV':'10.1111/ecog.02881','WeaklyInformativePriors':'10.1214/08-AOAS191'}
rows=[]
for name,doi in items.items():
 url='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='')
 row={'reference_id':name,'doi':doi,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  with urllib.request.urlopen(url,timeout=25) as response:raw=response.read();status=response.status
  p=D/(name+'.json');p.write_bytes(raw);m=json.loads(raw)['message'];assert m['DOI'].lower()==doi.lower()
  row.update(status='PRIMARY_METADATA_IDENTITY_FETCHED_NOT_CLAIM_REVIEWED',path=str(p.relative_to(W)),sha256=hashlib.sha256(raw).hexdigest(),title=m.get('title'),journal=m.get('container-title'),published=m.get('published'),authors=m.get('author'),type=m.get('type'))
 except Exception as e:row.update(status='UNVERIFIED',error=repr(e))
 rows.append(row);print(json.dumps({k:row[k] for k in ['reference_id','doi','status','title'] if k in row}),flush=True);time.sleep(.4)
(D/'reference_metadata_inventory.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False)+'\n')
