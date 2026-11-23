"""Generate the Figure 5 first-hit classifications and Figure 6 fits."""
import ctypes as ct
import numpy as np
from common17 import DATA, PTR, library, write_rows

CASES = [(m,a) for m,amps in [('coupled',[.1,.2044,.30,.39,.50]),('planar',[.39,1.,1.49])] for a in amps]
OFFSETS = [1,2,4,8]

def run():
    lib=library()
    x=np.linspace(-1.5,0,200); v=np.linspace(-.75,.75,200)
    X,V=np.meshgrid(x,v)
    mask=(X<0)&(.5*V**2-.5*X**2+.25*X**4 < -.01)
    starts=[np.ascontiguousarray(a[mask]) for a in (X,V)]
    arrays={'x':x,'v':v,'eligible':mask}
    summary=[]; fits=[]
    for idx,(model,A) in enumerate(CASES):
        grids={}
        for steps in (100,200):
            grid=np.zeros_like(X,dtype=np.uint8)
            hits=np.empty(mask.sum(),dtype=np.uint8)
            lib.first_cross(*(a.ctypes.data_as(PTR) for a in starts),hits.size,A,model=='planar',steps,hits.ctypes.data_as(ct.POINTER(ct.c_uint8)))
            grid[mask]=hits
            grids[steps]=grid
            arrays[f'case{idx}_s{steps}']=grid
        labels=grids[200]
        fractions=[]
        for k in OFFSETS:
            pairs=mask[:,k:] & mask[:,:-k]
            eligible=int(pairs.sum())
            different=int(((labels[:,k:]!=labels[:,:-k]) & pairs).sum())
            frac=different/eligible
            fractions.append(frac)
            fits.append(dict(case=idx,model=model,A=A,pixels=k,epsilon=k*(x[1]-x[0]),pairs=eligible,disagree=different,fraction=frac))
        alpha=np.polyfit(np.log(OFFSETS),np.log(fractions),1)[0]
        summary.append(dict(case=idx,model=model,A=A,eligible=int(mask.sum()),escaped=int(labels.sum()),fraction=labels[mask].mean(),alpha=alpha,agreement=np.mean(grids[100][mask]==labels[mask])))
        print(summary[-1],flush=True)
    np.savez_compressed(DATA/'basin_labels.npz',**arrays)
    write_rows('basin_summary.csv',summary)
    write_rows('uncertainty_fits.csv',fits)

if __name__=='__main__':
    run()
