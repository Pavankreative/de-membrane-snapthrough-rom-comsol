"""Planar primary stable/unstable arcs by time-domain collocation.
At the turning section (v=0), integrate backward to the saddle for Wu
and forward to it for Ws. Finite-tail and phase refinements are archived.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq, minimize_scalar
from common17 import T, write_rows
from melnikov_envelope import melnikov_margin

OMEGA=1.72

def rhs(t, z, amp):
    x,v=z
    return np.array([v,x-x**3-.08*v-.03*v*abs(v)+amp*(.12+(.5-.08*x*x)*x)*np.cos(OMEGA*t)])

@lru_cache(maxsize=32)
def saddle(amp):
    t=np.linspace(0,T,65)
    orbit=solve_bvp(lambda t,z:rhs(t,z,amp),lambda a,b:a-b,t,np.zeros((2,t.size)),tol=1e-10)
    if not orbit.success or abs(orbit.y[0]).max()>.2:
        raise RuntimeError('central periodic orbit failed')
    def variational(t,flat):
        x,v=orbit.sol(t)
        jac=np.array([[0.,1.],[1.-3*x*x+amp*(.5-.24*x*x)*np.cos(OMEGA*t),-.08-.06*abs(v)]])
        return (jac @ flat.reshape(2,2)).ravel()
    flow=solve_ivp(variational,(0,T),np.eye(2).ravel(),rtol=1e-10,atol=1e-12)
    multipliers=np.sort(abs(np.linalg.eigvals(flow.y[:,-1].reshape(2,2))))
    if not flow.success or not multipliers[0]<1<multipliers[1]:
        raise RuntimeError('periodic orbit is not hyperbolic')
    return orbit.sol

def splitting(amp, side, phi, length=14., tol=2e-7):
    # t runs outwards from the primary turn on both half-orbits.
    # At the far end, fix the displacement to the saddle displacement.
    t=np.linspace(0,length,int(12*length)+1)
    phase_time=phi/OMEGA
    orbit=saddle(amp)
    def center(t):
        return orbit(np.mod(t,T))
    x0=side*np.sqrt(2)/np.cosh(t)
    v0=-x0*np.tanh(t)
    initial=np.vstack([center(phase_time-t)+[x0,-v0],center(phase_time+t)+[x0,v0]])
    far=center(phase_time+np.array([-length,length]))[0]
    def ode(t,y):
        return np.vstack([-rhs(phase_time-t,y[:2],amp),rhs(phase_time+t,y[2:],amp)])
    def bc(a,b):
        return np.array([a[1],a[3],b[0]-far[0],b[2]-far[1]])
    sol=solve_bvp(ode,bc,t,initial,tol=tol,max_nodes=6000)
    if not sol.success:
        raise RuntimeError(sol.message)
    xu,xs=sol.y[[0,2],0]
    if min(side*xu,side*xs)<1:
        raise RuntimeError('turning section selected a different branch')
    return xu,xs,side*(xu-xs)

def run():
    samples=[]; summary=[]; zeros=[]; tails=[]
    for name,side,amps in [('right',1,[.199,.203,.204,.2044,.205,.208,.216]),('left',-1,[.32,.33,.333,.3336,.3345,.338,.35])]:
        phases=np.arange(32)*(2*np.pi/32)
        maxima={16:[],32:[]}
        for amp in amps:
            vals=np.array([splitting(amp,side,phi) for phi in phases])
            for n in (16,32):
                ix=np.arange(0,32,32//n)
                maxima[n].append(vals[ix,2].max())
                for k in ix:
                    samples.append(dict(loop=name,A=amp,N_phase=n,phase_rad=phases[k],xu=vals[k,0],xs=vals[k,1],d=vals[k,2]))
            print(name,amp,maxima[32][-1],flush=True)
        analytic=melnikov_margin(OMEGA,side*.12,.5,.08,.08,.03)[1]
        row=dict(loop=name,first_order=analytic,S=(4*.08/3+8*np.sqrt(2)*.03/15)/analytic)
        for n in (16,32):
            maxima_n=np.array(maxima[n])
            above=np.flatnonzero(maxima_n>0)
            if len(above)==0 or above[0]==0: raise RuntimeError('onset bracket missing')
            j=above[0]
            lo,hi=amps[j-1:j+1]
            a=lo-maxima_n[j-1]*(hi-lo)/(maxima_n[j]-maxima_n[j-1])
            row.update({f'lo{n}':lo,f'hi{n}':hi,f'A{n}':a})
        def continuous_peak(a):
            result=minimize_scalar(lambda phi:-splitting(a,side,phi)[2],
                                   bounds=(np.pi/2-np.pi/16,np.pi/2+np.pi/16),
                                   method='bounded',options={'xatol':1e-8})
            if not result.success: raise RuntimeError('phase optimisation failed')
            return -result.fun
        row['phase_optimised']=brentq(continuous_peak,row['lo32'],row['hi32'],xtol=1e-10)
        row['percent']=100*(row['A32']/analytic-1)
        row['phase_shift']=row['A32']-row['A16']
        summary.append(row)
        # Refine both simple roots above onset. The finite-difference derivative
        # is evaluated twice; sign reversal alone cannot distinguish cubic contact.
        for n in (16,32):
            stride=32//n
            for k in range(0,32,stride):
                nxt=(k+stride)%32
                if vals[k,2]*vals[nxt,2]>=0: continue
                a=phases[k]; b=a+2*np.pi/n
                d=lambda phi:splitting(amps[-1],side,phi)[2]
                root=brentq(d,a,b,xtol=1e-9)
                slopes=[(d(root+h)-d(root-h))/(2*h) for h in (.002,.001)]
                if abs(slopes[1])<1e-5 or abs(slopes[0]/slopes[1]-1)>.01:
                    raise RuntimeError('phase derivative unresolved')
                zeros.append(dict(loop=name,A=amps[-1],N_phase=n,phi_lo=a,phi_hi=b,d_lo=d(a),d_hi=d(b),root=root,slope_h002=slopes[0],slope_h001=slopes[1]))
        for a in (row['lo32'],row['hi32']):
            phase=max((s for s in samples if s['loop']==name and s['A']==a and s['N_phase']==32),key=lambda s:s['d'])
            longer=splitting(a,side,phase['phase_rad'],18,2e-9)[2]
            tails.append(dict(loop=name,A=a,phase=phase['phase_rad'],base=phase['d'],refined=longer,absolute_change=abs(longer-phase['d'])))
    if len(zeros)!=8: raise RuntimeError('expected two primary phase roots per grid and loop')
    write_rows('manifold_splitting_table3.csv',samples)
    write_rows('table3.csv',summary)
    write_rows('manifold_phase_zeros.csv',zeros)
    write_rows('manifold_tail_check.csv',tails)

if __name__=='__main__':
    run()
