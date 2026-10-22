"""Fixed inherited inference interface; contains no outcome reading."""
from group_features import *
from minimal_map import HierarchicalMAP
from iso000g_model_utils import weighted_ecdf
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
import warnings
EN_GRID=[dict(C=C,l1_ratio=l) for C in [.01,.1,1.,10.] for l in [0.,.25,.5,.75,1.]]
RF_GRID=[dict(n_estimators=300,max_depth=d,min_samples_leaf=l,max_features=m) for d in [6,12,None] for l in [10,50,200] for m in ['sqrt',.5]]
MAIN=['B0','B2','SOURCE_AWARE','READ_UNION','RANDOM_FOREST','ISOCRED_MAP']
OPTIONAL=['ELASTIC_NET','PSQAN_STYLE','SQANTI3_RULE']
def minimal(ev,base,meta):
 x=pd.DataFrame(index=meta.index)
 for rep in ['R1','R2']:
  x['LR_read_'+rep]=np.log1p(ev['read_union_'+rep])
  x['LR_junction_'+rep]=np.log1p(ev['junction_median_'+rep])
  x['LR_junction_fraction_'+rep]=ev['junction_fraction_'+rep]
  x['resolution_read_shared_'+rep]=ev['read_crossgroup_fraction_'+rep]
 x['annotation']=base.B2/2
 x['resolution_candidate_count']=np.log1p(meta.candidate_count)
 x['resolution_TSS_count']=np.log1p(meta.exact_TSS_count)
 x['resolution_window']=np.log1p(meta.window_opportunity_bp)
 x['resolution_monoexon']=meta.monoexon_fraction
 x['status_junction_NOT_APPLICABLE']=(ev.junction_unique_count==0).astype(float)
 return x
def make(name,params):
 if name=='ELASTIC_NET':return LogisticRegression(**params,solver='saga',penalty='elasticnet',max_iter=5000,random_state=SEED,n_jobs=1,class_weight=None)
 assert name=='RANDOM_FOREST'
 return RandomForestClassifier(**params,bootstrap=True,class_weight=None,random_state=SEED,n_jobs=1)
def fit_model(name,params,X,y):
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always')
  model=make(name,params).fit(X,y)
 return model,sum(issubclass(w.category,ConvergenceWarning) for w in caught)
def prob(model,X):
 return model.predict_proba(X)[:,list(model.classes_).index(1)]
def source_scores(mi_train,base_train,mi_final,base_final):
 q=[]
 for a,b in [(mi_train[['LR_read_R1','LR_read_R2']].mean(axis=1),mi_final[['LR_read_R1','LR_read_R2']].mean(axis=1)),(mi_train[['LR_junction_R1','LR_junction_R2']].mean(axis=1),mi_final[['LR_junction_R1','LR_junction_R2']].mean(axis=1)),(base_train.B2,base_final.B2)]:
  q.append(weighted_ecdf(a.to_numpy(),np.ones(len(a)),b.to_numpy()))
 q=np.column_stack(q)
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',RuntimeWarning)
  lr=np.nanmean(q[:,:2],axis=1)
 return dict(SOURCE_AWARE=np.mean(np.column_stack([lr,q[:,2]]),axis=1),READ_UNION=q[:,0])
def load_features():
 names=['core','evidence','baselines','metadata']
 return {n:pd.read_parquet(O/f'data/group_{n}_label_blind.parquet').set_index('promoter_group_id') for n in names}
