"""Recompute every Table 4 row and the stroboscopic points of Figure 8."""
import numpy as np
from common17 import DATA, PTR, library, write_rows

def run():
    lib=library()
    summaries=[]; raw=[]; sections={}
    for model,amps in [('planar',[.10,.25,.40,1.0,1.49]),('coupled',[.2044,.39,.50])]:
        planar=(model=='planar')
        for amp in amps:
            starts=np.array([[x,v,0 if planar else 1e-4,0] for x in np.linspace(-1.3,-.3,6) for v in np.linspace(-.4,.4,5)])
            vals=[lib.lyap(z.ctypes.data_as(PTR),amp,100,600,200,planar,1e-8) for z in starts]
            for i,(z,lam) in enumerate(zip(starts,vals)):
                raw.append(dict(model=model,A=amp,start=i,x=z[0],v=z[1],exponent=lam))
            summaries.append(dict(model=model,A=amp,median=np.median(vals),minimum=min(vals),maximum=max(vals),positive=int(np.sum(np.array(vals)>0))))
            z=starts[:8].copy()
            records=[np.empty((8,700)) for _ in range(3)]
            lib.integrate(z.ctypes.data_as(PTR),8,amp,700,200,planar,*(a.ctypes.data_as(PTR) for a in records))
            sections[f'{model}_{amp}']=np.stack([a[:,100:] for a in records[:2]])
    np.savez_compressed(DATA/'poincare_table4.npz',**sections)
    write_rows('table4.csv',summaries)
    write_rows('lyapunov_per_start.csv',raw)

if __name__=='__main__':
    run()
