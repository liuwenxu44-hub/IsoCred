"""Inherited exact-tie enrichment with explicit undefined metrics, no score flipping."""
import numpy as np
from sklearn.metrics import roc_auc_score,log_loss,brier_score_loss
def enrichment(y,s,w=None,fraction=.1):
 y=np.asarray(y,dtype=float);s=np.asarray(s,dtype=float)
 w=np.ones(len(y)) if w is None else np.asarray(w,dtype=float)
 assert y.shape==s.shape==w.shape and 0<fraction<=1
 assert not np.isinf(y).any() and set(y[np.isfinite(y)])<={0.,1.}
 keep=np.isfinite(s)&np.isfinite(y)&np.isfinite(w)&(w>0)
 y=y[keep];s=s[keep];w=w[keep]
 if not len(y):return (np.nan,np.nan,np.nan,np.nan)
 order=np.argsort(-s,kind='stable');cs=np.cumsum(w[order])
 k=np.searchsorted(cs,fraction*w.sum(),side='left')
 cut=s[order[min(k,len(y)-1)]];top=s>=cut
 p=w@y/w.sum();pt=w[top]@y[top]/w[top].sum()
 return (pt/p if p>0 else np.nan,pt,w[top].sum()/w.sum(),cut)
def metrics(y,s,w=None,probability=False):
 y=np.asarray(y,dtype=float);s=np.asarray(s,dtype=float)
 w=np.ones(len(y)) if w is None else np.asarray(w,dtype=float)
 assert y.shape==s.shape==w.shape
 assert not np.isinf(y).any() and set(y[np.isfinite(y)])<={0.,1.}
 den=np.isfinite(y)&np.isfinite(w)&(w>0)
 ok=np.isfinite(s)&den
 yy,ss,ww=y[ok],s[ok],w[ok]
 a,p,c,t=enrichment(yy,ss,ww)
 r=dict(n=int(ok.sum()),weighted_n=float(ww.sum()),scoring_coverage=float(w[ok].sum()/w[den].sum()) if den.any() else np.nan,
 active_row_fraction=float(ok.mean()) if len(ok) else np.nan,
 enrichment=a,top_positive_rate=p,top_actual_coverage=c,threshold=t,
 prevalence=float(ww@yy/ww.sum()) if len(yy) else np.nan,
 AUC=float(roc_auc_score(yy,ss,sample_weight=ww)) if len(np.unique(yy))==2 else np.nan,
 AUC_status='DEFINED' if len(np.unique(yy))==2 else 'UNDEFINED_SINGLE_CLASS_OR_EMPTY')
 if probability:
  assert ((ss>=0)&(ss<=1)).all()
  r.update(Brier=float(brier_score_loss(yy,ss,sample_weight=ww)) if len(yy) else np.nan,
   log_loss=float(log_loss(yy,ss,sample_weight=ww,labels=[0,1])) if len(yy) else np.nan)
 return r
