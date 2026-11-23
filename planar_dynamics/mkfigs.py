import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from crk import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({"font.family":"serif","font.size":8.5,"axes.grid":False})
cm=ListedColormap(['#2b6cb0','#111111','#d1453b'])      # -1 left well, 0 unsettled, +1 right well
xs=np.linspace(-1.6,1.6,200); vs=np.linspace(-1,1,200); xs3=np.linspace(-1.6,1.6,300); vs3=np.linspace(-1,1,300)
b3=np.load('results/basins3.npz')
panels=[('p, A=0.1  (A<A_c)',xs3,vs3,b3['p0.1']),('p, A=0.2  (A≈A_c)',xs,vs,np.load('b200_p0.2.npy')),('p, A=0.5',xs,vs,np.load('b200_p0.5.npy')),
        ('c, A=0  (unforced)',xs,vs,np.load('b200_c0.0.npy')),('c, A=0.2',xs,vs,np.load('b200_c0.2.npy')),('c, A=0.4',xs,vs,np.load('b200_c0.4.npy'))]
fig,ax=plt.subplots(2,3,figsize=(7.4,5.0),sharex=True,sharey=True)
for a,(t,x,v,O) in zip(ax.ravel(),panels):
    a.pcolormesh(x,v,O,cmap=cm,vmin=-1,vmax=1,shading='auto',rasterized=True)
    a.set_title(t.replace('p,','planar').replace('c,','coupled'),fontsize=8.5)
for a in ax[1]: a.set_xlabel("x(0)")
for a in ax[:,0]: a.set_ylabel("x'(0)")
fig.tight_layout(); fig.savefig('image_Figure5.png',dpi=250)
# Fig 6 Poincare
fig,ax=plt.subplots(2,3,figsize=(7.4,5.0))
def load(t): return np.load(f'results/poin_{t}.npz')
for a,(t,ttl) in zip(ax[0],[('c0.30','coupled A=0.30'),('c0.40','coupled A=0.40'),('c0.50','coupled A=0.50')]):
    d=load(t); a.plot(d['x'],d['v'],'.',ms=1.2,color='k'); a.set_title(ttl+'  (x, x\')',fontsize=8.5); a.set_xlabel('x'); a.set_ylabel("x'")
d=load('c0.50'); ax[1,0].plot(d['x'],d['r'],'.',ms=1.2,color='k'); ax[1,0].set_title('coupled A=0.50  (x, ρ)',fontsize=8.5); ax[1,0].set_xlabel('x'); ax[1,0].set_ylabel('ρ')
for a,(t,ttl) in zip(ax[1,1:],[('p0.60','planar A=0.60'),('p1.49','planar A=1.49')]):
    d=load(t); a.plot(d['x'],d['v'],'o',ms=4,color='k'); a.set_title(ttl+'  (x, x\')',fontsize=8.5); a.set_xlabel('x'); a.set_ylabel("x'")
fig.tight_layout(); fig.savefig('image_Figure6.png',dpi=250)
# Fig 7 exponents
Ag=np.array([0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.85,1.0,1.15,1.3])
Lm=np.array([nlyap([-0.95,0,0,0],A,200,300) for A in Ag]); Lp=np.array([nlyap([0.95,0,0,0],A,200,300) for A in Ag])
ly=np.load('results/lyap.npz'); sw=np.array([[0.05,0.99],[0.1,0.99],[0.15,0.98],[0.2,0.81],[0.25,0.62],[0.3,0.62],[0.4,0.49],[0.5,0.28],[0.7,0.22]])
fig,ax=plt.subplots(1,3,figsize=(7.6,2.9))
ax[0].plot(sw[:,0],sw[:,1],'o-',color='k'); ax[0].axvline(0.2044,ls=':',color='gray'); ax[0].set_xlabel('A'); ax[0].set_ylabel('uncertainty exponent α'); ax[0].set_title('planar basins',fontsize=8.5); ax[0].set_ylim(0,1.1)
ax[1].plot(Ag,Lm,'o-',label='left well (x≈−1)'); ax[1].plot(Ag,Lp,'s--',label='right well (x≈+1)'); ax[1].axhline(0,color='k',lw=.5)
for v,c in ((0.2044,'gray'),(0.385,'C3')): ax[1].axvline(v,ls=':',color=c)
ax[1].set_xlabel('A'); ax[1].set_ylabel('λ_ρ (normal exponent)'); ax[1].legend(fontsize=7); ax[1].set_title('azimuthal stability of in-well orbit',fontsize=8.5)
ax[2].plot(ly['As'],ly['LE'].mean(1),'o-',label='coupled (4D)'); ax[2].plot(ly['Ap'],ly['LP'].mean(1),'s--',label='planar (ρ≡0)')
ax[2].axhline(0,color='k',lw=.5); ax[2].axvline(0.2044,ls=':',color='gray'); ax[2].axvline(0.385,ls=':',color='C3')
ax[2].set_xlabel('A'); ax[2].set_ylabel('largest Lyapunov exponent'); ax[2].legend(fontsize=7); ax[2].set_title('attractor chaos',fontsize=8.5)
fig.tight_layout(); fig.savefig('image_Figure7.png',dpi=250)
print(np.round(Lm,4)); print(np.round(Lp,4))
