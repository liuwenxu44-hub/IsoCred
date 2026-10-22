"""Deterministic scientific-number ledger and descriptive report; no journal/winner claims."""
import csv,hashlib,json,math
from pathlib import Path
ADEQUACY=dict(groups=200,positives=20,negatives=20,components=20)
EVIDENCE={'development':'POST_HOC_AFTER_HISTORICAL_FINAL_EXPOSURE_DEVELOPMENT_ONLY',
 'external':'PRESPECIFIED_NEW_COMMON10_TRANSFER_UNDER_CELL_GENERATOR_ALIGNER_SHIFTS'}
def pointer(s):return str(s).replace('~','~0').replace('/','~1')
def finite(x):return isinstance(x,(int,float)) and math.isfinite(x)
def fmt(x):
 if x is None or isinstance(x,float) and not math.isfinite(x):return 'NA'
 return format(x,'.6g') if isinstance(x,(float,int)) else str(x)
def full(x):
 if x is None or isinstance(x,float) and not math.isfinite(x):return 'NA'
 if isinstance(x,float):return format(x,'.17g')
 if isinstance(x,(list,dict)):return json.dumps(x,sort_keys=True,separators=(',',':'))
 return str(x)
def adequacy(pop,metrics):
 n=int(pop['n']);p=metrics['prevalence'];positive=int(round(n*p)) if n and p is not None else 0
 if n and (p is None or not math.isclose(n*p,positive,abs_tol=1e-8)):raise ValueError('NONINTEGER_POINT_POSITIVE_COUNT')
 negative=n-positive;components=int(pop['component_n'])
 counts=dict(groups=n,positives=positive,negatives=negative,components=components)
 flags=[k.upper()+'_BELOW_'+str(v) for k,v in ADEQUACY.items() if counts[k]<v]
 if not n:flags.append('EMPTY_SCORING_SUPPORT')
 elif not positive or not negative:flags.append('SINGLE_CLASS_ENDPOINT')
 return dict(n=n,positives=positive,negatives=negative,component_n=components,
  adequacy_status='LOW_ADEQUACY_UNDERPOWERED_FLAG' if flags else 'ABOVE_FIXED_DESCRIPTIVE_ADEQUACY_THRESHOLDS',
  adequacy_flags=flags,adequacy_interpretation='PRESPECIFIED_DESCRIPTIVE_FLAGS_NOT_FORMAL_POWER_CALCULATION')
def build_ledgers(result):
 ledger=[];populations={}
 case=result['case']
 for kind,rows in [('native',result['native_method_results']),('paired',result['paired_contrast_results'])]:
  for i,row in enumerate(rows):
   pop=row['population'];point=row['metrics'] if kind=='native' else row['left_metrics']
   a=adequacy(pop,point);key=(row['scope'],row['endpoint'],pop['population_id'])
   descriptor=row['method'] if kind=='native' else row['left']+'-'+row['right']
   if key not in populations:
    populations[key]=dict(case=case,evidence_status=EVIDENCE[case],scope=row['scope'],endpoint=row['endpoint'],
     population_id=pop['population_id'],eligible_mother_n=pop['eligible_mother_n'],
     excluded_nonfinite_n=pop['excluded_nonfinite_n'],excluded_mono_n=pop['excluded_mono_n'],
     included_mono_n=pop['included_mono_n'],scoring_coverage=pop['scoring_coverage'],methods_or_contrasts=[],
     source_json_pointer=('/native_method_results/' if kind=='native' else '/paired_contrast_results/')+str(i),
     positive_count_rule='round(population.n * point_prevalence);integer_consistency_checked',**a)
   populations[key]['methods_or_contrasts'].append(descriptor)
   values=row['metrics'].items() if kind=='native' else [(row['metric'],row['difference_left_minus_right'])]
   for metric,value in values:
    ci=(row.get('CI95',{}).get(metric,{}) if kind=='native' else row.get('CI95')) or {}
    bonf=row.get('Bonferroni_CI') or {}
    source='/native_method_results/'+str(i)+'/metrics/'+pointer(metric) if kind=='native' else '/paired_contrast_results/'+str(i)+'/difference_left_minus_right'
    ident=[case,kind,row['scope'],row['endpoint'],descriptor,metric,pop['population_id']]
    ledger.append(dict(value_id='ISO002D_'+hashlib.sha256(json.dumps(ident,separators=(',',':')).encode()).hexdigest()[:24],
     case=case,evidence_status=EVIDENCE[case],kind=kind,scope=row['scope'],endpoint=row['endpoint'],
     method=row.get('method',''),left=row.get('left',''),right=row.get('right',''),metric=metric,
     population_id=pop['population_id'],n=a['n'],positives=a['positives'],negatives=a['negatives'],component_n=a['component_n'],
     eligible_mother_n=pop['eligible_mother_n'],scoring_coverage=pop['scoring_coverage'],
     excluded_mono_n=pop['excluded_mono_n'],included_mono_n=pop['included_mono_n'],estimate_full_precision=value,
     left_estimate=row.get('left_estimate'),right_estimate=row.get('right_estimate'),
     CI95_low=ci.get('low'),CI95_high=ci.get('high'),CI95_valid_draws=ci.get('valid_draws'),CI95_invalid_draws=ci.get('invalid_draws'),
     Bonferroni_low=bonf.get('low'),Bonferroni_high=bonf.get('high'),Bonferroni_confidence=bonf.get('confidence'),
     Bonferroni_valid_draws=bonf.get('valid_draws'),Bonferroni_tail_expected_draws=bonf.get('lower_tail_expected_draws'),
     Bonferroni_family_id=row.get('Bonferroni_family_id',''),
     favorable_direction=row.get('favorable_direction','NOT_APPLICABLE_TO_NATIVE_VALUE'),
     adequacy_status=a['adequacy_status'],adequacy_flags=a['adequacy_flags'],
     undefined_reason=row.get('undefined',{}).get(metric,'') if kind=='native' else row.get('undefined_reason'),
     source_file='evaluation.json',json_pointer=source,row_json_pointer=('/native_method_results/' if kind=='native' else '/paired_contrast_results/')+str(i),display_rule='six_significant_digits_in_MD;17_significant_digits_in_TSV;full_JSON_retained'))
 for i,r in enumerate(result.get('declared_stratum_inventory',[])):
  if r['engine_scope_exists']:continue
  key=(r['scope'],r['endpoint'],r['population_id'])
  populations[key]=dict(case=case,evidence_status=EVIDENCE[case],scope=r['scope'],endpoint=r['endpoint'],population_id=r['population_id'],
   eligible_mother_n=0,excluded_nonfinite_n=0,excluded_mono_n=0,included_mono_n=0,scoring_coverage=None,
   methods_or_contrasts=['ALL_FIXED_METHODS_AND_CONTRASTS_UNDEFINED_EMPTY_DECLARED_LEVEL'],
   source_json_pointer='/declared_stratum_inventory/'+str(i),positive_count_rule='explicit_empty_predeclared_level',
   **adequacy({'n':0,'component_n':0},{'prevalence':None}))
 for p in populations.values():p['methods_or_contrasts']=sorted(set(p['methods_or_contrasts']))
 return ledger,list(populations.values())
def resolve_pointer(obj,p):
 for key in p.strip('/').split('/'):
  key=key.replace('~1','/').replace('~0','~')
  obj=obj[int(key)] if isinstance(obj,list) else obj[key]
 return obj
def verify_ledger(result,ledger):
 seen=set()
 for r in ledger:
  if r['value_id'] in seen:raise ValueError('DUPLICATE_VALUE_ID')
  seen.add(r['value_id']);v=resolve_pointer(result,r['json_pointer'])
  if v!=r['estimate_full_precision']:raise ValueError('NUMBER_SOURCE_POINTER_MISMATCH')
  source=resolve_pointer(result,r['row_json_pointer']);native=r['kind']=='native';metric=r['metric']
  pop=source['population'];point=source['metrics'] if native else source['left_metrics'];a=adequacy(pop,point)
  ci=(source.get('CI95',{}).get(metric,{}) if native else source.get('CI95')) or {};bc=source.get('Bonferroni_CI') or {}
  expected=dict(n=a['n'],positives=a['positives'],negatives=a['negatives'],component_n=a['component_n'],
   eligible_mother_n=pop['eligible_mother_n'],scoring_coverage=pop['scoring_coverage'],excluded_mono_n=pop['excluded_mono_n'],included_mono_n=pop['included_mono_n'],
   left_estimate=source.get('left_estimate'),right_estimate=source.get('right_estimate'),
   CI95_low=ci.get('low'),CI95_high=ci.get('high'),CI95_valid_draws=ci.get('valid_draws'),CI95_invalid_draws=ci.get('invalid_draws'),
   Bonferroni_low=bc.get('low'),Bonferroni_high=bc.get('high'),Bonferroni_confidence=bc.get('confidence'),
   Bonferroni_valid_draws=bc.get('valid_draws'),Bonferroni_tail_expected_draws=bc.get('lower_tail_expected_draws'))
  for key,value in expected.items():
   if r[key]!=value:raise ValueError('NUMBER_AUXILIARY_FIELD_MISMATCH:'+key)
  if r['population_id']!=pop['population_id'] or r['adequacy_status']!=a['adequacy_status'] or r['adequacy_flags']!=a['adequacy_flags']:
   raise ValueError('NUMBER_POPULATION_OR_ADEQUACY_MISMATCH')
 return dict(status='PASS',values=len(ledger),all_full_precision_values_resolved=True,
  verified_fields='point_estimates;left/right;all_reported_CI_fields;valid/invalid_draws;population_counts;coverage;integer-derived_positive/negative_counts;adequacy')
def markdown(result,ledger,populations):
 case=result['case'];n=result['primary_eligible_population']['n'];boot=result['bootstrap']
 lines=['# Common10 CAGE-support comparison','',
  'Evidence: '+EVIDENCE[case]+'.','',
  'This comparison evaluates common10 promoter-CTSS-support predictions. It does not test the original RF71 model, establish transcript function, isolate a cell-type effect, or assess previously unseen genes. The negative result for the original MAP model is unchanged.','',
  'Development was post hoc after historical final exposure. External uses a new prespecified 2/2 endpoint; inherited development uses 2/3. Absolute enrichment across these two datasets is not treated as the same measurement scale.','',
  'All signs, endpoints and fixed strata are retained. Eligible negatives mean no reproducible support in this assay, not biological absence. Diagnostics are ranking scores, not calibrated probabilities.','',
  '## Fixed analysis and uncertainty','',
  'Primary-eligible groups: '+str(n)+'. Bootstrap: '+str(boot['draws'])+' shared dependency-component draws, seed '+str(boot['seed'])+', '+boot['generator']+'. Each sampled component retains its multiplicity for every member; components are not biological donors or libraries.','',
  'Each method uses its own finite support; each contrast uses its pair-common finite support. Junction diagnostic NA does not remove mono groups from RF10–EN10. Different population IDs are different estimands and cannot establish a cross-population winner.','',
  'All exact score-boundary ties are included. Requested selection fractions are 0.05, 0.1 and 0.2; actual coverage and thresholds are in scientific_numbers.tsv. No low-AUC score flipping or threshold tuning is performed.','',
  'Intervals are nominal 95% paired percentile bootstrap. Only each dataset’s seven pooled primary enrichment@0.1 contrasts additionally have nominal 99.2857142857% Bonferroni marginal intervals. There is no exact or cross-dataset/all-endpoint FWER guarantee. At 2000 draws the adjusted tail contains about 7.14 draws before undefined draws; all valid/undefined draw counts are retained. Nonsignificance is not equivalence.','',
  '## All pooled method estimates','',
  '| Endpoint | Method | n | Positive | Scoring coverage | Enrichment@0.1 | AUC | Brier | Log loss |',
  '|---|---|---:|---:|---:|---:|---:|---:|---:|']
 for r in result['native_method_results']:
  if r['scope']!='pooled':continue
  m=r['metrics'];a=adequacy(r['population'],m)
  lines.append('| '+' | '.join(map(fmt,[r['endpoint'],r['method'],a['n'],a['positives'],m['scoring_coverage'],m['enrichment@0.1'],m['AUC'],m.get('Brier'),m.get('log_loss')]))+' |')
 lines+=['','## All pooled primary and sensitivity signed contrasts','',
  'Each row is left minus right, recomputed on the population shown. A positive enrichment/AUC difference favors the left score within that population only. Brier/log-loss differences and all additional fixed-coverage results remain in the complete number ledger.','',
  '| Endpoint | Contrast | Metric | n | Population suffix | Difference | Nominal 95% CI | Valid draws | Bonferroni CI (primary only) |',
  '|---|---|---|---:|---|---:|---|---:|---|']
 for r in result['paired_contrast_results']:
  if r['scope']!='pooled' or r['metric'] not in ('enrichment@0.1','AUC'):continue
  ci=r.get('CI95') or {};bc=r.get('Bonferroni_CI') or {}
  ci_text='['+fmt(ci.get('low'))+', '+fmt(ci.get('high'))+']'
  bc_text='['+fmt(bc.get('low'))+', '+fmt(bc.get('high'))+']' if bc else 'not in family'
  lines.append('| '+' | '.join(map(fmt,[r['endpoint'],r['left']+' - '+r['right'],r['metric'],r['population']['n'],
   r['population']['population_id'][-12:],r['difference_left_minus_right'],ci_text,ci.get('valid_draws'),bc_text]))+' |')
 if case=='development':
  lines+=['','## Every frozen development fold','',
   '| Fold | Method | n | Enrichment@0.1 | AUC | Scoring coverage |','|---|---|---:|---:|---:|---:|']
  for r in result['native_method_results']:
   if r['scope'].startswith('fold:'):
    m=r['metrics'];lines.append('| '+' | '.join(map(fmt,[r['scope'],r['method'],r['population']['n'],m['enrichment@0.1'],m['AUC'],m['scoring_coverage']]))+' |')
 lines+=['','## Every fixed stratum and support-population adequacy','',
  'Prespecified descriptive flags: groups <200, positives <20, negatives <20, or components <20. “Underpowered” flags are not a formal power calculation and never delete or suppress a stratum. Annotation bins are annotation scores, not known/novel genes. Junction-support bins use raw R1+R2 median junction support, not E1 read union. NOT_APPLICABLE is retained.','',
  '| Scope | Endpoint | Population suffix | n | Positive | Negative | Components | Adequacy flags |',
  '|---|---|---|---:|---:|---:|---:|---|']
 for r in populations:
  flags=', '.join(r['adequacy_flags']) or 'above fixed thresholds'
  lines.append('| '+' | '.join(map(fmt,[r['scope'],r['endpoint'],r['population_id'][-12:],r['n'],r['positives'],r['negatives'],r['component_n'],flags]))+' |')
 lines+=['','## Complete result files','',
  'evaluation.json contains every method, fixed contrast, coverage fraction, endpoint, fold and stratum at full precision. scientific_numbers.tsv contains '+str(len(ledger))+' numeric values linked to their JSON source fields, including negative effects and undefined metrics; adequacy_flags.tsv lists every support population. Displayed finite numbers use six significant digits, while TSV values use 17 significant digits. Adequacy count rows retain JSON source pointers; positives are round(n × point prevalence), with integer consistency checked. Empty bins remain explicit n=0 rows. Low-adequacy strata, saturation and undefined endpoints are reported as observed.','']
 return '\n'.join(lines)
def write_tsv(path,rows):
 if not rows:raise ValueError('EMPTY_LEDGER')
 with Path(path).open('x',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader()
  for r in rows:writer.writerow({k:full(v) for k,v in r.items()})
def write_reports(output,result):
 output=Path(output);target=output/'evaluation.json'
 if target.exists():raise FileExistsError(target)
 raw=json.dumps(result,indent=2,allow_nan=False)+'\n';target.write_text(raw)
 digest=hashlib.sha256(raw.encode()).hexdigest()
 ledger,populations=build_ledgers(result);check=verify_ledger(result,ledger)
 for row in ledger:row['source_sha256']=digest
 write_tsv(output/'scientific_numbers.tsv',ledger);write_tsv(output/'adequacy_flags.tsv',populations)
 with (output/'scientific_report.md').open('x') as f:f.write(markdown(result,ledger,populations))
 with (output/'number_source_check.json').open('x') as f:json.dump(dict(check,source_sha256=digest),f,indent=2);f.write('\n')
 return dict(values=len(ledger),populations=len(populations),low_adequacy_populations=sum(r['adequacy_status'].startswith('LOW_') for r in populations),
  all_fixed_strata_retained=True,formal_power_analysis_performed=False,number_source_check=check)
