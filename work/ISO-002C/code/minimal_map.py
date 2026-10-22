"""Finite hierarchical MAP model of a group-level support event (not function)."""
import numpy as np
from scipy.special import expit
from scipy.optimize import minimize
from numpy.polynomial.hermite import hermgauss
class HierarchicalMAP:
 def fit(self,X,y,components,weight=None,source=None):
  X=np.asarray(X);y=np.asarray(y);n,p=X.shape;weight=np.ones(n) if weight is None else np.asarray(weight)
  comps,ci=np.unique(components,return_inverse=True);nc=len(comps)
  # Exchangeable coefficients within prespecified evidence-source blocks.
  source=np.zeros(p,dtype=int) if source is None else np.asarray(source,dtype=int);S=int(source.max())+1
  # theta=[intercept,beta_p,source_means_S,training_component_intercepts_nc]
  def fun(t):
   b=t[1:1+p];mu=t[1+p:1+p+S];u=t[1+p+S:];z=t[0]+X@b+u[ci]
   loss=np.sum(weight*(np.logaddexp(0,z)-y*z));diff=b-mu[source]
   loss+=.5*np.dot(diff,diff)+.5*np.dot(mu,mu)+2*np.dot(u,u)+.5*t[0]**2/25
   e=weight*(expit(z)-y);g=np.zeros_like(t);g[0]=e.sum()+t[0]/25;g[1:1+p]=X.T@e+diff;g[1+p:1+p+S]=mu-np.bincount(source,weights=diff,minlength=S);g[1+p+S:]=np.bincount(ci,weights=e,minlength=nc)+4*u
   return loss,g
  initial=np.zeros(1+p+S+nc);pr=np.average(y,weights=weight);initial[0]=np.log(pr/(1-pr))
  fit=minimize(fun,initial,jac=True,method='L-BFGS-B',options={'maxiter':2000,'ftol':1e-10,'gtol':1e-5,'maxls':40})
  if not fit.success:raise RuntimeError('MAP optimization failed: '+str(fit.message))
  self.theta=fit.x;self.p=p;self.S=S;self.loss=float(fit.fun);self.iterations=int(fit.nit);self.max_gradient=float(abs(fit.jac).max());self.fit_components=set(comps)
  return self
 def predict_proba(self,X,components):
  assert not self.fit_components.intersection(components),'group/component leakage'
  z=self.theta[0]+np.asarray(X)@self.theta[1:1+self.p]
  # Held-out component is unobserved: integrate N(0,0.5^2) prior, not its E4 label.
  q,w=hermgauss(15);p=expit(z[:,None]+np.sqrt(2)*.5*q)@(w/np.sqrt(np.pi))
  return p
