"""Fit p = a*|w| - b*|w|**3 to FE pressure-deflection data and report the fold and the cubic zero.

Usage:  python fit_cubic.py [csv] [w_max_mm ...]
Default: data/membrane_frozen_mean_control.csv, fit ranges 4, 6 and 10 mm.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parents[1]
path = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "data" / "membrane_frozen_mean_control.csv"
ranges = [float(x) for x in sys.argv[2:]] or [4.0, 6.0, 10.0]

with open(path) as f:
    r = list(csv.DictReader(f))
w = np.array([-float(x["dh_mm"]) for x in r])
p = np.array([float(x["pu_Pa"]) for x in r])
i = p.argmax()
print(f"{path.name}: {len(w)} points, FE maximum p = {p[i]:.2f} Pa at |w| = {w[i]:.2f} mm")
for wmax in ranges:
    m = (w > 0) & (w <= wmax)
    A = np.c_[w[m], w[m] ** 3]
    a, mb = np.linalg.lstsq(A, p[m], rcond=None)[0]
    b = -mb
    rms = math.sqrt(np.mean((A @ [a, mb] - p[m]) ** 2))
    msg = f"  fit to {wmax:4.1f} mm: a = {a:.3f} Pa/mm, b = {b:.4f} Pa/mm^3, rms = {rms:.3f} Pa"
    if b > 0:
        msg += f"; cubic fold at {math.sqrt(a / (3 * b)):.2f} mm, p_max = {2 * a * math.sqrt(a / (3 * b)) / 3:.1f} Pa, zero at {math.sqrt(a / b):.2f} mm"
    print(msg)
