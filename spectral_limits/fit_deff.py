"""Fit effective bending rigidity D_eff and offset delta to the stable-side axisymmetric spectra.
   python3 fit_deff.py            (reads data/raw/2D_pinned_1000_raw.txt and 2D_clamped_250_raw.txt)
Model: measured modes at xx equal plate-dispersion modes at x_true = xx - delta; first four distinct w2 per xx.
Manuscript values: pinned D = 6.738e-7, delta = -5.1e-6 ; clamped D = 6.719e-7, delta = +1.06e-6."""
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from comsol_io import read_export
from plate_dispersion import modes
H = Path(__file__).parent / 'data' / 'raw'

def first_four(x, w, xx, tol=5e-3):
    v = np.sort(w[np.isclose(x, xx, rtol=1e-9)])
    cl = [[v[0]]]
    for a in v[1:]:
        (cl[-1].append(a) if abs(a - cl[-1][-1]) <= tol * max(abs(a), 1e-3) else cl.append([a]))
    return np.array([np.mean(c) for c in cl[:4]])

def fit(rim, fname, xs):
    x, w = read_export(H / fname)
    data = {xx: first_four(x, w, xx) for xx in xs}
    def res(p):
        D, dl = p[0] * 1e-6, p[1] * 1e-5; out = []
        for xx, v in data.items():
            m = modes(xx - dl, D, rim, 0, 4, wmax=4500, ngrid=9000)
            if len(m) < 4: return np.ones(24) * 1e3
            out += list((m - v) / v)
        return np.array(out)
    s = least_squares(res, [0.7, 0.0], x_scale=[0.1, 1])
    return s.x[0] * 1e-6, s.x[1] * 1e-5, np.sqrt(np.mean(s.fun**2)), abs(s.fun).max()

if __name__ == '__main__':
    P = (-5e-5, -1e-4, -2e-4, -4e-4, -6e-4, -8e-4)
    C = (-2e-5, -5e-5, -1e-4, -2e-4, -4e-4, -8e-4)
    for rim, f, xs in (('pinned', '2D_pinned_1000_raw.txt', P), ('clamped', '2D_clamped_250_raw.txt', C)):
        D, dl, rms, mx = fit(rim, f, xs)
        print(f'{rim:8s} D = {D:.4e} N m   delta = {dl:+.3e}   rms rel resid {rms:.2e}  max {mx:.2e}')
