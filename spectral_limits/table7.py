"""Table 7: zero-frequency crossings x_c (voltage excess).  python3 table7.py  -> table7_crossings.csv
 - plate_case_A / plate_case_B : plate dispersion with D = 3.37e-8 and 0.98e-6 N m (assumed rigidities)
 - plate_fitted                : D_eff = 6.73e-7 N m (fitted to the axisymmetric stable-side spectra)
 - comsol_extrapolated         : quadratic extrapolation from stable-side COMSOL spectra (n = 0 only; rank 1 = (0,1), rank 2 = (0,2)); n>=1: not obtained"""
import csv
from pathlib import Path
from plate_dispersion import KR, crossing
from fit_xc import crossings
H = Path(__file__).parent
DEFF = 6.73e-7
sets = {'pinned': ('data/raw/2D_pinned_1000_raw.txt', 'pinned 1000 el'),
        'clamped': ('data/raw/2D_clamped_250_raw.txt', 'clamped 250 el')}
rows = []
for rim, (f, mesh) in sets.items():
    xc, xl = crossings(H / f, verbose=False)
    meas = {(0, 1): xc[0], (0, 2): xc[1]}; lin = {(0, 1): xl[0], (0, 2): xl[1]}
    for mode, kR in KR[rim].items():
        rows.append(dict(rim=rim, n=mode[0], m=mode[1], kR=kR,
                         plate_case_A=crossing(3.37e-8, kR), plate_case_B=crossing(0.98e-6, kR),
                         plate_fitted=crossing(DEFF, kR),
                         comsol_extrapolated=meas.get(mode, ''), comsol_linear_only=lin.get(mode, ''),
                         comsol_mesh=mesh if mode in meas else 'not obtained (n>=1)'))
out = H / 'table7_crossings.csv'
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, rows[0].keys()); w.writeheader()
    for r in rows: w.writerow({k: (f'{v:.6e}' if isinstance(v, float) else v) for k, v in r.items()})
for r in rows: print(r['rim'], (r['n'], r['m']), f"fitted {r['plate_fitted']:.6e}", 'comsol', r['comsol_extrapolated'] if r['comsol_extrapolated']=='' else f"{r['comsol_extrapolated']:.4e}")
