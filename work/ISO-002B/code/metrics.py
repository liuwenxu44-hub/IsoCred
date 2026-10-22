import numpy as np
from sklearn.metrics import roc_auc_score,log_loss,brier_score_loss
def enrichment(y,s,w=None):
 y=np.asarray(y);s=np.asarray(s);w=np.ones(len(y)) if w is None else np.asarray(w)
 keep=np.isfinite(s)&np.isfinite(y)&(w>0);y=y[keep];s=s[keep];w=w[keep]
 if not len(y) or (w@y)==0:return (np.nan,np.nan,np.nan,np.nan)
 order=np.argsort(-s,kind='stable');cs=np.cumsum(w[order]);k=np.searchsorted(cs,.1*w.sum(),side='left');cut=s[order[min(k,len(y)-1)]];top=s>=cut
 p=w@y/w.sum();pt=w[top]@y[top]/w[top].sum()
 return (pt/p,pt,w[top].sum()/w.sum(),cut)
def metrics(y,s):
 ok=np.isfinite(s);a,p,c,t=enrichment(np.asarray(y)[ok],np.asarray(s)[ok]);r=dict(n=int(ok.sum()),scoring_coverage=float(ok.mean()),enrichment=a,top_positive_rate=p,top_actual_coverage=c,threshold=t,prevalence=float(np.mean(np.asarray(y)[ok])),AUC=float(roc_auc_score(np.asarray(y)[ok],np.asarray(s)[ok])))
 return r
