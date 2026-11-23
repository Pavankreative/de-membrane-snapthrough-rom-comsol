import sys,time; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from crk import *
# (a) transient-length check on 100x100 subgrid
N=100; xs=np.linspace(-1.6,1.6,N); vs=np.linspace(-1,1,N); X,V=np.meshgrid(xs,vs)
z0=np.array([X.ravel(),V.ravel(),np.full(X.size,1e-4),np.zeros(X.size)])
for tag,A,pl in [('c0.2',0.2,0),('c0.3',0.3,0),('p0.5',0.5,1)]:
    out=[]
    for nper in (50,200,400):
        z,(rx,rv,rr)=integ(z0,A,nper,planar=bool(pl),record=True); last=rx[:,-10:]
        O=np.zeros(N*N,int); O[(last<-0.5).all(1)]=-1; O[(last>0.5).all(1)]=1; out.append(round((O==0).mean(),3))
    print('other fraction vs run length (50,200,400 periods)',tag,out,flush=True)
# (b) uncertainty exponent from 300x300 3-outcome maps, small eps only
d=np.load('results/basins3.npz'); dx=d['xs'][1]-d['xs'][0]; ks=np.array([1,2,3,4])
print('uncertainty fraction f(k*dx), k=1..4, and slope')
for t in [t for t in d.files if t[0] in 'cp']:
    O=d[t]; f=np.array([0.5*((O[:,:-k]!=O[:,k:]).mean()+(O[:-k,:]!=O[k:,:]).mean()) for k in ks])
    print(t,np.round(f,4),round(np.polyfit(np.log(ks*dx),np.log(np.maximum(f,1e-9)),1)[0],2))
# (c) Lyapunov scan from 6 initial conditions
ics=[(-0.95,0),(-1.0,0.1),(-0.9,-0.1),(1.0,0.0),(0.95,0.05),(-1.1,0.0)]
As=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.7,1.0,1.3]
print('Lyapunov (full 4D, seeded rho=1e-4): A: mean max-LE over ICs, range')
LE={}
for A in As:
    v=[lyap([x,xd,1e-4,0],A,300,500) for x,xd in ics]; LE[A]=v; print(A,round(np.mean(v),4),round(min(v),4),round(max(v),4),flush=True)
print('planar rho=0 Lyapunov')
LP={}
for A in [0.1,0.2,0.3,0.5,0.8,1.0,1.3,1.49,1.6,1.8,2.0]:
    v=[lyap([x,xd,0,0],A,300,500,planar=True) for x,xd in ics]; LP[A]=v; print(A,round(np.mean(v),4),round(min(v),4),round(max(v),4),flush=True)
np.savez('results/lyap.npz',As=np.array(As),LE=np.array([LE[a] for a in As]),Ap=np.array(list(LP)),LP=np.array([LP[a] for a in LP]))
# (d) normal (rho) exponent along planar in-well trajectories
print('normal exponent of rho along planar trajectory from (-0.95,0)')
NL={}
for A in [0.1,0.2,0.3,0.35,0.4,0.45,0.5,0.7,1.0,1.3]:
    NL[A]=nlyap([-0.95,0,0,0],A,100,300); print(A,round(NL[A],4),flush=True)
np.savez('results/nlyap.npz',A=np.array(list(NL)),L=np.array(list(NL.values())))
