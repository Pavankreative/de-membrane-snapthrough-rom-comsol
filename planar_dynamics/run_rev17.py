"""Execute the manuscript analysis from source. Run from any directory."""
from pathlib import Path
import os
import subprocess
import sys
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4',MPLBACKEND='Agg')
    for name in ('direct_manifold_onset','run_basin_20p_first_crossing','verify_lyap',
                 'run_table5_trajectories','seed_brackets','make_rev17_figures','verify_rev17'):
        subprocess.run([sys.executable,str(ROOT/(name+'.py'))],cwd=ROOT,env=env,check=True)
