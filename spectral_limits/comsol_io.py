"""Read COMSOL 'Global Evaluation' text exports and return (xx, w2).

w2 is recomputed from the complex Eigenfrequency column, not taken from the
exported w2 column:  lambda = 2*pi*i*f ;  w2 = Re[-lambda*(lambda + alpha)],  alpha = 14.4 1/s.
The w2 column exported by COMSOL in the raw files used a different expression and is NOT used.
Overdamped modes then appear as two entries with equal w2 (a real-root pair).
"""
import re
import numpy as np
ALPHA = 14.4
_num = re.compile(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')

def read_export(path):
    xs, ws = [], []
    for ln in open(path):
        if ln.startswith('%') or not ln.strip():
            continue
        p = ln.split()
        x = float(p[0])
        f = complex(p[1].replace('i', 'j'))
        lam = 2j * np.pi * f
        xs.append(x); ws.append((-lam * (lam + ALPHA)).real)
    return np.array(xs), np.array(ws)

def write_csv(path, xs, ws):
    with open(path, 'w') as fh:
        fh.write('xx,w2_per_s2\n')
        for a, b in zip(xs, ws):
            fh.write(f'{a:.6g},{b:.10g}\n')
