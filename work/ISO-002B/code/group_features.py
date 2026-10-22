from common import *
from iso000g_model_utils import design
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler
from sklearn.pipeline import make_pipeline
def candidate_design(core):
 x,co,bi,ca=design(core)
 # Missing numerical support and NOT_APPLICABLE remain distinct states.
 return pd.concat([x[co+bi],pd.get_dummies(x[ca],prefix=ca,dtype=float)],axis=1)
def group_design(core,membership):
 x=candidate_design(core.loc[membership.candidate_id].reset_index(drop=True))
 x.index=membership.promoter_group_id
 return x.groupby(level=0,sort=True).mean()
def prep(zero=False):
 return make_pipeline(SimpleImputer(strategy='constant' if zero else 'median',fill_value=0,add_indicator=not zero,keep_empty_features=True),RobustScaler())
def read_base(membership):
 data={}
 for short,alias,col in [('B0','B0_read_fixed','B0_READ_FIXED'),('B1','B1_junction_fixed','B1_JUNCTION_FIXED'),('B2','B2_annotation_fixed','B2_ANNOTATION_FIXED')]:
  b=pd.read_parquet(I/alias).set_index('candidate_id');data[short]=b.loc[membership.candidate_id,col].to_numpy()
 return pd.DataFrame(data,index=membership.promoter_group_id).groupby(level=0,sort=True).mean()
