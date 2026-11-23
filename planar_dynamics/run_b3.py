import sys,time; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from crk import *
N=300; xs=np.linspace(-1.6,1.6,N); vs=np.linspace(-1,1,N); X,V=np.meshgrid(xs,vs)
z0=np.array([X.ravel(),V.ravel(),np.full(X.size,1e-4),np.zeros(X.size)])
cases=[('c0.1',0.1,0),('c0.2',0.2,0),('c0.3',0.3,0),('c0.4',0.4,0),('p0.1',0.1,1),('p0.2',0.2,1),('p0.5',0.5,1),('p1.49',1.49,1)]
res={}
for tag,A,pl in cases:
    z,(rx,rv,rr)=integ(z0,A,50,planar=bool(pl),record=True)
    last=rx[:,-10:]
    O=np.zeros(N*N,int); O[(last<-0.5).all(1)]=-1; O[(last>0.5).all(1)]=1
    res[tag]=O.reshape(N,N); print(tag,'W-',round((O==-1).mean(),3),'W+',round((O==1).mean(),3),'other',round((O==0).mean(),3),flush=True)
np.savez('results/basins3.npz',xs=xs,vs=vs,**res)
