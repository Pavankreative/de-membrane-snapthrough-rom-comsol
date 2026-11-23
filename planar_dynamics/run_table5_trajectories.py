"""Table 5: DOP853 reference trajectories. Count at the manuscript sampling interval."""
import numpy as np
from scipy.integrate import solve_ivp
from common17 import DATA, T, write_rows

def trajectory(gamma,knn,amp,damped=True,unforced_time=None,seed=1e-4):
    delta,delta2=(.08,.03) if damped else (0.,0.)
    def rhs(t,y):
        x,v,r,u=y
        force=amp*(.12+(.50-.08*x*x)*x)*np.cos(1.72*t)
        return [v,x-x**3-gamma*x*r*r-delta*v-delta2*v*abs(v)+force,
                u,1.593**2*r-knn*r**3-gamma*x*x*r-delta*u-delta2*u*abs(u)]
    end=120*T if unforced_time is None else unforced_time
    t=np.linspace(0,end,24001) if unforced_time is None else np.r_[np.arange(0,end,T/200),end]
    sol=solve_ivp(rhs,(0,end),[-.95,0,seed,0],method='DOP853',t_eval=t,rtol=1e-9,atol=1e-11)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.t,sol.y

def run():
    rows=[]; traces={}
    for g,k,amps in [(.15,1.25,[0]),(3.,1.25,[.5,1.]),(3.,8.,[0.,.1,.2,.3,.4,.5,1.])]:
        for a in amps:
            for damped in ([True,False] if g==.15 else [True]):
                end=None if damped else 300.
                t,z=trajectory(g,k,a,damped,end)
                x,v,r,u=z; idx=np.flatnonzero(x[1:]*x[:-1]<0)
                rows.append(dict(gamma=g,Knn=k,A=a,damped=int(damped),end=t[-1],n_cross=len(idx),max_abs_rho=np.max(abs(r)),final_x=x[-1],final_rho=r[-1],first_cross_right_time=t[idx[0]+1] if len(idx) else '',tail_max_abs_x=np.max(abs(x[t>100*T])) if end is None else ''))
                traces[f'run{len(rows)-1}']=np.stack([t,x,r])
                if k==8 and a in (.3,.5): traces[f'velocity{len(rows)-1}']=v
    # Include the Figure 4 sweep; overlapping reference cases use the same solve.
    sweep=[]
    for model,seed in [('coupled',1e-4),('planar',0.)]:
        for amp in (.1,.2,.3,.4,.5,.7,1.,1.3,1.6):
            match=[r for r in rows if r['gamma']==3 and r['Knn']==8 and r['A']==amp]
            if match and model=='coupled':
                crossings=match[0]['n_cross']
            else:
                _,z=trajectory(3,8,amp,seed=seed)
                crossings=int(np.sum(z[0,1:]*z[0,:-1]<0))
            sweep.append(dict(model=model,A=amp,crossings=crossings))
    write_rows('figure4_sweep.csv',sweep)
    np.savez_compressed(DATA/'table5_traces.npz',**traces)
    write_rows('table5.csv',rows)

if __name__=='__main__':
    run()
