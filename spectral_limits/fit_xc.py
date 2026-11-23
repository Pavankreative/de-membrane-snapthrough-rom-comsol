"""Stable-side extrapolation of the zero-frequency crossings x_c (Table 7, 'COMSOL' column).
   python3 fit_xc.py data/raw/2D_pinned_1000_raw.txt [xmax=-2e-5] [nmodes=4] [rel_tol=5e-3]
Per xx <= xmax the distinct w2 values are sorted ascending (rank i = i-th axisymmetric mode; w2 = a_i (x_c,i - x)
is ordered by wavenumber on the stable side).  Each rank is fitted with a quadratic in xx; the real root nearest 0 is x_c.
All sampled xx are NEGATIVE: every reported crossing is an extrapolation, not a measured sign change."""
import sys
import numpy as np
from comsol_io import read_export

def crossings(path, xmax=-2e-5, nm=4, tol=5e-3, verbose=True):
    x, w = read_export(path)
    xs = sorted({a for a in x if a <= xmax})
    tracks = [[] for _ in range(nm)]
    for xx in xs:
        v = np.sort(w[x == xx]); cl = [[v[0]]]
        for a in v[1:]:
            (cl[-1].append(a) if abs(a - cl[-1][-1]) <= tol * max(abs(a), 1e-3) else cl.append([a]))
        if verbose: print(f'xx={xx:+.2e}  cluster sizes {[len(c) for c in cl[:nm + 2]]}')
        for i in range(min(nm, len(cl))): tracks[i].append((xx, np.mean(cl[i])))
    xc, xlin = [], []
    for t in tracks:
        X = np.array([a for a, _ in t]); W = np.array([b for _, b in t])
        c = np.polyfit(X, W, 2); r = np.roots(c); r = r[np.isreal(r)].real
        xc.append(r[np.argmin(abs(r))]); l = np.polyfit(X, W, 1); xlin.append(-l[1] / l[0])
    return np.array(xc), np.array(xlin)

if __name__ == '__main__':
    a = sys.argv
    xc, xl = crossings(a[1], float(a[2]) if len(a) > 2 else -2e-5, int(a[3]) if len(a) > 3 else 4,
                       float(a[4]) if len(a) > 4 else 5e-3)
    for i, (q, l) in enumerate(zip(xc, xl), 1): print(f'rank {i}: x_c = {q:.4e}   (linear-only {l:.4e})')
    print('differences to rank 1:', [f'{q - xc[0]:.4e}' for q in xc[1:]])
