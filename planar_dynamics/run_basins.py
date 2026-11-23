import sys,time; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from crk import *
def grid(N):
    xs=np.linspace(-1.6,1.6,N); vs=np.linspace(-1,1,N); X,V=np.meshgrid(xs,vs)
    return xs,vs,np.array([X.ravel(),V.ravel(),np.full(X.size,1e-4),np.zeros(X.size)])
xs,vs,z0=grid(200); res={}
cases=[('c0.1',0.1,0),('c0.2',0.2,0),('c0.3',0.3,0),('c0.4',0.4,0),('c0.5',0.5,0),('c1.0',1.0,0),('p0.5',0.5,1),('p1.49',1.49,1),('p1.8',1.8,1)]
for tag,A,pl in cases:
    z,_=integ(z0,A,40,planar=bool(pl)); res[tag]=np.sign(z[0]).reshape(200,200); print(tag,round((res[tag]>0).mean(),3),flush=True)
np.savez('results/basins200.npz',xs=xs,vs=vs,**res)
