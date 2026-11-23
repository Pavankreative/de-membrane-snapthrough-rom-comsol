"""Shared paths, tabular output and compiler settings for revision 17."""
from pathlib import Path
import csv
import ctypes as ct
import subprocess
import numpy as np
ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
FIG = ROOT / 'figures'
DATA.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
T = 2*np.pi/1.72
PTR = ct.POINTER(ct.c_double)

def write_rows(name, rows):
    with (DATA/name).open('w', newline='') as f:
        out = csv.DictWriter(f, fieldnames=list(rows[0]))
        out.writeheader()
        out.writerows(rows)

def read_rows(name):
    with (DATA/name).open(newline='') as f:
        return list(csv.DictReader(f))

def library():
    target = ROOT/'librk17.so'
    source = ROOT/'rk.c'
    if not target.exists() or source.stat().st_mtime_ns > target.stat().st_mtime_ns:
        subprocess.run(['gcc','-O3','-ffp-contract=off','-shared','-fPIC','-fopenmp',str(source),'-o',str(target),'-lm'],check=True)
    lib = ct.CDLL(str(target))
    lib.lyap.argtypes = [PTR,ct.c_double,ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_double]
    lib.lyap.restype = ct.c_double
    lib.first_cross.argtypes = [PTR,PTR,ct.c_int,ct.c_double,ct.c_int,ct.c_int,ct.POINTER(ct.c_uint8)]
    lib.integrate.argtypes = [PTR,ct.c_int,ct.c_double,ct.c_int,ct.c_int,ct.c_int,PTR,PTR,PTR]
    return lib
