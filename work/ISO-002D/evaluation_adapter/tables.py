"""Label-preserving adapters and fixed descriptive bins; memory-only, no file I/O."""
import math
FEATURES=('LR_junction_R1','LR_junction_fraction_R1','LR_junction_R2','LR_junction_fraction_R2','annotation',
 'resolution_candidate_count','resolution_TSS_count','resolution_window','resolution_monoexon','status_junction_NOT_APPLICABLE')
METHODS=('RF10','EN10','ANNOTATION_ONLY','JUNCTION_MEDIAN_SUPPORT','WINDOW_OPPORTUNITY_ONLY')
LABEL_COLUMNS={'development':{'primary':'primary_target'},'external':{
 'primary':'primary_label','technical_pair_1':'GSM5625727_vs_GSM5625730_label',
 'technical_pair_2':'GSM5625728_vs_GSM5625730_label','technical_pair_3':'GSM5625729_vs_GSM5625730_label',
 'secondary_position_concordant_5bp':'secondary_5bp_label'}}
STRATUM_LEVELS={'annotation_score':['0','(0,1)','1'],'candidate_count':['1','2-4','>=5'],
 'window_opportunity_bp':['101','102-125','>125'],'monoexon_fraction':['none','mixed','all'],
 'junction_support_sum':['0','(0,5)','[5,20)','[20,100)','>=100','NOT_APPLICABLE']}
def require(v,name):
 if not v:raise ValueError(name)
def boolean(x):
 if isinstance(x,str):
  require(x.lower() in ('true','false'),'INVALID_BOOLEAN');return x.lower()=='true'
 require(x is not None and x in (True,False,0,1),'INVALID_BOOLEAN')
 return bool(x)
def numeric(x,nullable=False):
 if x is None or isinstance(x,str) and x in ('NA','NaN','nan',''):
  require(nullable,'MISSING_NUMERIC');return None
 v=float(x)
 if not math.isfinite(v):
  require(nullable and math.isnan(v),'NONFINITE_NUMERIC');return None
 return v
def index(rows,role):
 out={}
 for r in rows:
  g=r.get('promoter_group_id')
  require(isinstance(g,str) and bool(g) and g.strip()==g,'INVALID_GROUP_ID:'+role)
  require(g not in out,'DUPLICATE_GROUP_ID:'+role);out[g]=r
 return out
def close(a,b):
 return a is None and b is None or a is not None and b is not None and math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)
def positive_integer(x,name,minimum=1):
 v=numeric(x);require(v>=minimum and int(v)==v,'INVALID_INTEGER:'+name);return int(v)
def bins(annotation,candidate_count,window,mono,junction_sum):
 require(0<=annotation<=1 and 0<=mono<=1,'BIN_RANGE')
 require(candidate_count>=1 and int(candidate_count)==candidate_count,'CANDIDATE_BIN_RANGE')
 require(window>=101 and int(window)==window,'WINDOW_BIN_RANGE')
 a='0' if annotation==0 else '1' if annotation==1 else '(0,1)'
 c='1' if candidate_count==1 else '2-4' if candidate_count<=4 else '>=5'
 w='101' if window==101 else '102-125' if window<=125 else '>125'
 m='none' if mono==0 else 'all' if mono==1 else 'mixed'
 if junction_sum is None:j='NOT_APPLICABLE'
 else:
  require(junction_sum>=0 and math.isfinite(junction_sum),'JUNCTION_BIN_RANGE')
  j='0' if junction_sum==0 else '(0,5)' if junction_sum<5 else '[5,20)' if junction_sum<20 else '[20,100)' if junction_sum<100 else '>=100'
 return dict(zip(STRATUM_LEVELS,(a,c,w,m,j)))
def prepare(case,label_rows,prediction_rows,metadata_rows,evidence_rows,feature_rows,component_rows=None,*,expected_n=None,baselines=None):
 require(case in LABEL_COLUMNS,'INVALID_CASE')
 tables={k:index(v,k) for k,v in {'labels':label_rows,'predictions':prediction_rows,'metadata':metadata_rows,
  'evidence':evidence_rows,'features':feature_rows}.items()}
 if component_rows is not None:tables['components']=index(component_rows,'components')
 if baselines is not None:tables['baselines']=index(baselines,'baselines')
 ids=sorted(tables['labels'])
 require(all(set(t)==set(ids) for t in tables.values()),'NONIDENTICAL_GROUP_ID_SETS')
 if expected_n is not None:require(len(ids)==expected_n,'WRONG_FIXED_GROUP_COUNT')
 labels=[];predictions=[];metadata=[];strata={k:{} for k in STRATUM_LEVELS};folds=set()
 for g in ids:
  m=tables['metadata'][g];ev=tables['evidence'][g];x=tables['features'][g];p=tables['predictions'][g];y=tables['labels'][g]
  require(set(x)==set(FEATURES)|{'promoter_group_id'},'WRONG_EXACT_COMMON10_SCHEMA')
  require(set(p)==set(METHODS)|{'promoter_group_id'},'WRONG_EXACT_FIVE_METHOD_SCHEMA')
  nc=positive_integer(m['candidate_count'],'candidate_count');nt=positive_integer(m['exact_TSS_count'],'exact_TSS_count')
  window=positive_integer(m['window_opportunity_bp'],'window_opportunity_bp',101)
  mono=numeric(m['monoexon_fraction']);require(0<=mono<=1,'MONO_FRACTION_RANGE')
  nj=positive_integer(ev['junction_unique_count'],'junction_unique_count',0);is_mono=nj==0
  if 'junction_unique_count' in m:require(nj==positive_integer(m['junction_unique_count'],'metadata_junction_count',0),'METADATA_EVIDENCE_JUNCTION_COUNT_MISMATCH')
  require(is_mono==(mono==1),'MONO_JUNCTION_INCONSISTENCY')
  med=[];expected={}
  for rep in ('R1','R2'):
   v=numeric(ev['junction_median_'+rep],nullable=True);f=numeric(ev['junction_fraction_'+rep],nullable=True)
   if is_mono:require(v is None and f is None,'MONO_MEASUREMENT_NOT_NA')
   else:require(v is not None and v>=0 and f is not None and 0<=f<=1,'MISSING_OR_INVALID_JUNCTION_MEASUREMENT')
   med.append(v);expected['LR_junction_'+rep]=None if v is None else math.log1p(v);expected['LR_junction_fraction_'+rep]=f
  annotation=numeric(x['annotation']);require(0<=annotation<=1,'ANNOTATION_RANGE')
  b2=tables['baselines'][g]['B2'] if 'baselines' in tables else m.get('B2_annotation_fixed_group_mean')
  require(b2 is not None and close(annotation,numeric(b2)/2),'ANNOTATION_SOURCE_INCONSISTENCY')
  expected.update(annotation=annotation,resolution_candidate_count=math.log1p(nc),resolution_TSS_count=math.log1p(nt),
   resolution_window=math.log1p(window),resolution_monoexon=mono,status_junction_NOT_APPLICABLE=float(is_mono))
  require(all(close(numeric(x[k],nullable=True),v) for k,v in expected.items()),'COMMON10_RAW_PRIMITIVE_MISMATCH')
  pp={'promoter_group_id':g}
  for method in METHODS:
   v=numeric(p[method],nullable=method=='JUNCTION_MEDIAN_SUPPORT')
   if method in ('RF10','EN10'):require(v is not None and 0<=v<=1,'SUPERVISED_NOT_COMPLETE_PROBABILITY')
   pp[method]=v
  diag={'ANNOTATION_ONLY':annotation,'JUNCTION_MEDIAN_SUPPORT':None if is_mono else sum(math.log1p(v) for v in med)/2,
   'WINDOW_OPPORTUNITY_ONLY':math.log1p(window)}
  require(all(close(pp[k],v) for k,v in diag.items()),'FROZEN_DIAGNOSTIC_MISMATCH')
  eligible=boolean(y['evaluation_eligible'] if case=='development' else y['assay_eligible'])
  if case=='development':require(eligible,'DEVELOPMENT_LOADER_MUST_ALREADY_SELECT_FIXED_ELIGIBLE_IDS')
  yy={'promoter_group_id':g,'primary_eligible':eligible}
  for key,col in LABEL_COLUMNS[case].items():
   v=numeric(y[col],nullable=True);require(v is None or v in (0,1),'NONBINARY_LABEL')
   require(not eligible or v is not None,'MISSING_LABEL_IN_PRIMARY_ELIGIBLE_SET')
   yy[key]=None if v is None else int(v)
  if case=='external' and eligible:
   require(yy['secondary_position_concordant_5bp']<=yy['primary'],'SECONDARY_NOT_SUBSET_OF_PRIMARY')
   require(all(yy['technical_pair_'+str(k)]<=yy['primary'] for k in (1,2,3)),'TECHNICAL_PAIR_NOT_SUBSET_OF_POOLED_PRIMARY')
  cm=tables['components'][g] if 'components' in tables else m
  component=cm['dependency_component_id']
  require(isinstance(component,str) and component and component.strip()==component,'INVALID_COMPONENT_ID')
  mm=dict(promoter_group_id=g,dependency_component_id=component,junction_not_applicable=is_mono)
  if case=='development':
   fold=positive_integer(m['inner_fold'],'inner_fold',0);require(fold in range(5),'UNKNOWN_FROZEN_DEVELOPMENT_FOLD')
   mm['fold_id']=str(fold);folds.add(fold)
   require(m['split']=='DEVELOPMENT_POOL' and not boolean(m['touches_final_dependency']),'FORBIDDEN_DEVELOPMENT_DOMAIN')
  level=bins(annotation,nc,window,mono,None if is_mono else sum(med))
  for name,value in level.items():strata[name][g]=value
  labels.append(yy);predictions.append(pp);metadata.append(mm)
 if case=='development' and expected_n is not None:require(folds==set(range(5)),'MISSING_FROZEN_DEVELOPMENT_FOLD')
 return dict(labels=labels,predictions=predictions,metadata=metadata,label_columns={k:k for k in LABEL_COLUMNS[case]},
  strata=strata,case=case,fold_column='fold_id' if case=='development' else None)
