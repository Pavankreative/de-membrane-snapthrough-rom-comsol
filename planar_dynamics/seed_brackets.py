"""Local finite-horizon crossing brackets near the reported Table 2 transitions.
Bisection of two opposite labels is a local diagnostic; finite-time escape
is not assumed monotone across all amplitudes.
"""
import numpy as np
from run_table5_trajectories import trajectory
from common17 import write_rows

def run():
    rows=[]; summaries=[]
    for seed in (1e-2,1e-4,1e-8,1e-12,0.):
        def label(a):
            t,z=trajectory(3,8,a,seed=seed)
            hits=np.flatnonzero(z[0,1:]*z[0,:-1]<0)
            rows.append(dict(seed=seed,A=a,changes=len(hits),maximum_rho=max(abs(z[2]))))
            return len(hits)>0
        grid=([k/100 for k in range(30,71,2)] if seed else [k/100 for k in range(130,171,5)])
        left=grid[0]
        if label(left): raise RuntimeError('initial amplitude already crosses')
        for right in grid[1:]:
            if label(right): break
            left=right
        else: raise RuntimeError(f'unbracketed seed {seed}')
        while right-left>.001:
            mid=(left+right)/2
            if label(mid): right=mid
            else: left=mid
        summaries.append(dict(seed=seed,no_cross=left,cross=right,ratio_left=right/.33355766410496485,ratio_right=right/.20438657575275765))
    write_rows('seed_brackets.csv',summaries)
    write_rows('seed_trials.csv',rows)
if __name__=='__main__':
    run()
