"""Exact free-vibration dispersion relation of a prestressed Kirchhoff-Love circular plate.

  D lap^2 w - T lap w + rho_h w_tt = 0 ,   T = N0*(-x)  (x = voltage excess (V/Vm)^2 - 1; x<0 stable side)
  Solutions J_n(k r) + I_n(b r); rim 'clamped' (w = w_r = 0) or 'pinned' (w = 0, M_r = 0, nu = 0.5).
Constants are those of the manuscript (lambda_p = 1.8, H0 = 0.5 mm, Rs = 90 mm, rho = 960 kg/m3).
"""
import numpy as np
from scipy.special import jv, ive
from scipy.optimize import brentq

N0 = 14.352607647580          # N/m, in-plane membrane force at lambda_p = 1.8
R = 0.09                      # m
RHO_H = 960 * 0.5e-3 / 1.8**2  # kg/m2, areal mass of the pre-stretched membrane
NU = 0.5

def _jd(n, z):  return 0.5 * (jv(n - 1, z) - jv(n + 1, z))
def _jdd(n, z): return 0.25 * (jv(n - 2, z) - 2 * jv(n, z) + jv(n + 2, z))
def _id(n, z):  return 0.5 * (ive(n - 1, z) + ive(n + 1, z))
def _idd(n, z): return 0.25 * (ive(n - 2, z) + 2 * ive(n, z) + ive(n + 2, z))

def det(w2, x, D, rim, n=0):
    T = -N0 * x
    d = np.sqrt(T * T + 4 * D * RHO_H * w2)
    s1, s2 = (T + d) / (2 * D), (T - d) / (2 * D)
    if s2 >= 0:
        return np.nan
    b, k = np.sqrt(s1), np.sqrt(-s2)
    kR, bR = k * R, b * R
    J, I = jv(n, kR), ive(n, bR)       # exponentially scaled I_n: common factor drops out of the sign
    if rim == 'clamped':
        return J * b * _id(n, bR) - I * k * _jd(n, kR)
    Jp, Jpp = k * _jd(n, kR), k * k * _jdd(n, kR)
    Ip, Ipp = b * _id(n, bR), b * b * _idd(n, bR)
    Mj = Jpp + NU * (Jp / R - n * n * J / R**2)
    Mi = Ipp + NU * (Ip / R - n * n * I / R**2)
    return J * Mi - I * Mj

def modes(x, D, rim, n=0, nm=4, wmax=3000.0, ngrid=6000):
    """Lowest nm squared frequencies w2 (1/s^2) of azimuthal order n."""
    g = np.linspace(1e-6, wmax, ngrid)
    v = [det(w, x, D, rim, n) for w in g]
    out = []
    for a, b_, u, vv in zip(g, g[1:], v, v[1:]):
        if np.isfinite(u) and np.isfinite(vv) and u * vv < 0:
            out.append(brentq(det, a, b_, args=(x, D, rim, n)))
            if len(out) >= nm:
                break
    return np.array(out)

# first zeros of the wavenumber condition k_R at zero tension and zero frequency
KR = {'pinned':  {(0, 1): 2.165871, (1, 1): 3.689191, (2, 1): 5.031408, (0, 2): 5.427433},
      'clamped': {(0, 1): 3.831706, (1, 1): 5.135622, (2, 1): 6.380162, (0, 2): 7.015587}}

def crossing(D, kR):
    """Zero-frequency crossing x_c = D (kR)^2 / (N0 R^2)."""
    return D * kR**2 / (N0 * R**2)
