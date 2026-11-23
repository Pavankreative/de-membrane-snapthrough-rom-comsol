"""Re-compute selected Table 4 entries and compare with results/lyap.npz and results/nlyap.npz.   python3 verify_lyap.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crk import *
here = os.path.dirname(os.path.abspath(__file__))
ics = [(-0.95,0),(-1.0,0.1),(-0.9,-0.1),(1.0,0.0),(0.95,0.05),(-1.1,0.0)]
L = np.load(os.path.join(here, 'results/lyap.npz')); N = np.load(os.path.join(here, 'results/nlyap.npz'))
for A in (0.2, 0.3, 0.5):
    i = list(L['As']).index(A); v = [lyap([x, xd, 1e-4, 0], A, 300, 500) for x, xd in ics]
    print(f'coupled A={A}: median {np.median(v):+.4f} (stored {np.median(L["LE"][i]):+.4f})  max|diff| {np.abs(np.array(v)-L["LE"][i]).max():.2e}')
for A in (0.3, 0.5):
    i = list(L['Ap']).index(A); v = [lyap([x, xd, 0, 0], A, 300, 500, planar=True) for x, xd in ics]
    print(f'planar  A={A}: median {np.median(v):+.4f} (stored {np.median(L["LP"][i]):+.4f})  max|diff| {np.abs(np.array(v)-L["LP"][i]).max():.2e}')
for A in (0.3, 0.4, 0.5):
    i = list(N['A']).index(A); v = nlyap([-0.95, 0, 0, 0], A, 100, 300)
    print(f'normal exponent A={A}: {v:+.5f} (stored {N["L"][i]:+.5f})')
