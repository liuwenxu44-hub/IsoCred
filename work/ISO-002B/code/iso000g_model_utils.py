"""Frozen label-blind transformations. No fitting or data access at import."""
import numpy as np,pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler,OneHotEncoder
from sklearn.impute import SimpleImputer
def design(frame):
 X=frame.copy();cat=[c for c in X if c.endswith('_status')]
 binary=[c for c in X if c not in cat and (pd.api.types.is_bool_dtype(X[c]) or c.startswith(('annotation_','sqanti_','e1_observed')) or 'all_junctions_' in c or c in ['struct_monoexon','struct_gene_assignment_uncertain'])]
 continuous=[c for c in X if c not in cat+binary]
 for col in continuous+binary:
  X[col]=pd.to_numeric(X[col],errors='raise').astype(float)
  if col in continuous and any(s in col for s in ['unique_read_ids','unique_to_candidate','ambiguous_read_ids','_cpm_','unique_support','transcript_span','junction_count','exon_count']):X[col]=np.log1p(X[col])
 return X,continuous,binary,cat
def preprocessor(continuous,binary,cat):
 return ColumnTransformer([('continuous',Pipeline([('impute',SimpleImputer(strategy='median',add_indicator=True)),('scale',RobustScaler())]),continuous),('binary',SimpleImputer(strategy='most_frequent',add_indicator=True),binary),('categories',OneHotEncoder(handle_unknown='ignore',sparse_output=False),cat)],remainder='drop',sparse_threshold=0)
def weighted_ecdf(train,weight,values):
 keep=np.isfinite(train);t=train[keep];w=weight[keep]
 if not len(t):return np.full(len(values),np.nan)
 order=np.argsort(t,kind='stable');t=t[order];cum=np.cumsum(w[order]);ix=np.searchsorted(t,values,side='right')-1
 out=np.zeros(len(values));ok=ix>=0;out[ok]=cum[ix[ok]]/cum[-1];out[~np.isfinite(values)]=np.nan;return out
