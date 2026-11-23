import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from crk import *
from scipy.integrate import solve_bvp
# (3) numerical M2 for the model with forcing f_n in the rho equation (gamma=3): scale per f_n^2
P_=dict(sig2=1.593**2,gam=3.0,w=1.72); s=np.sqrt(P_['sig2']); Lh=14.0
taus=np.linspace(-Lh,Lh,2001)
def M2_unit(tau0):
    # rho1'' = (sig2 - 2 gam sech^2) rho1 + cos(w(tau+tau0)),  rho1 -> 0 at +-inf (f_n = 1)
    f=lambda t,y:np.vstack([y[1],(P_['sig2']-2*P_['gam']/np.cosh(t)**2)*y[0]+np.cos(P_['w']*(t+tau0))])
    bc=lambda ya,yb:np.array([ya[1]-s*ya[0],yb[1]+s*yb[0]])
    t=np.linspace(-Lh,Lh,1201); sol=solve_bvp(f,bc,t,np.zeros((2,t.size)),tol=1e-8,max_nodes=40000)
    r=sol.sol(taus)[0]; rp=sol.sol(taus)[1]; x0sq=2/np.cosh(taus)**2
    return np.trapezoid(x0sq*r*rp,taus)   # = gamma^{-1} * M2/(eps^2 f_n^2)
vals=[M2_unit(t0) for t0 in np.linspace(0,2*np.pi/P_['w'],25)]
print('max |int x0^2 rho1 rho1\' dtau| per f_n^2 (gamma=3):',np.max(np.abs(vals)),' => max|M2|/(eps^2 f_n^2)=',3*np.max(np.abs(vals)))
M1=lambda w,f,g1,g3,d,d2:np.pi*w*(np.sqrt(2)*f/np.cosh(np.pi*w/2)+w*(g1-(w*w+4)/6*g3)/np.sinh(np.pi*w/2))
print('first-order amplitude pi*w[...] at reference params:',M1(1.72,0.12,0.5,0.08,0.08,0.03))
# (4) Poincare sections
for tag,A,pl in [('c0.30',0.30,0),('c0.40',0.40,0),('c0.50',0.50,0),('p1.49',1.49,1),('p0.60',0.6,1)]:
    z0=np.array([[-0.95],[0],[0.0 if pl else 1e-4],[0]]); _,(rx,rv,rr)=integ(z0,A,4200,record=True,planar=bool(pl))
    np.savez(f'results/poin_{tag}.npz',x=rx[0,200:],v=rv[0,200:],r=rr[0,200:])
    print(tag,'x range',rx[0,200:].min().round(3),rx[0,200:].max().round(3),'distinct pts approx',len(np.unique(np.round(rx[0,200:],3))))
