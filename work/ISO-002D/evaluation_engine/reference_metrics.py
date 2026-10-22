"""Independent slow reference: Python loops/pair enumeration, no production imports."""
import math
def reference_metrics(y,s,w=None,fraction=.1,probability=False):
 w=[1.0]*len(y) if w is None else [float(z) for z in w]
 eligible=[i for i in range(len(y)) if math.isfinite(y[i]) and w[i]>0]
 ix=[i for i in eligible if math.isfinite(s[i])]
 total=sum(w[i] for i in ix);den=sum(w[i] for i in eligible)
 p=sum(w[i]*y[i] for i in ix)/total if total else math.nan
 out=dict(weighted_n=total,scoring_coverage=total/den if den else math.nan,prevalence=p,AUC=math.nan)
 if total:
  order=sorted(ix,key=lambda i:-s[i]);cum=0
  for i in order:
   cum+=w[i]
   if cum>=fraction*total:
    threshold=s[i];break
  selected=[i for i in ix if s[i]>=threshold];sw=sum(w[i] for i in selected)
  ppv=sum(w[i]*y[i] for i in selected)/sw
  out.update(enrichment=ppv/p if p else math.nan,PPV=ppv,actual_selection_coverage=sw/total,threshold=threshold)
  pos=[i for i in ix if y[i]==1];neg=[i for i in ix if y[i]==0]
  pn=sum(w[i] for i in pos)*sum(w[i] for i in neg)
  if pn:
   out['AUC']=sum(w[i]*w[j]*(1.0 if s[i]>s[j] else .5 if s[i]==s[j] else 0.0) for i in pos for j in neg)/pn
 else:
  out.update(enrichment=math.nan,PPV=math.nan,actual_selection_coverage=math.nan,threshold=math.nan)
 if probability:
  eps=2.220446049250313e-16
  out['Brier']=sum(w[i]*(y[i]-s[i])**2 for i in ix)/total if total else math.nan
  out['log_loss']=-sum(w[i]*(y[i]*math.log(min(1-eps,max(eps,s[i])))+(1-y[i])*math.log(1-min(1-eps,max(eps,s[i])))) for i in ix)/total if total else math.nan
 return out
