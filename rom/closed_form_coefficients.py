"""Closed-form reduced-order coefficients (K1, K3) and the mean-deflection pressure law.

Implements the formulas of docs/theory.md section 3 for an equibiaxially pre-stretched
Gent membrane under a Maxwell stress, and maps them to the pressure law

    p(w) = a*|w| - b*|w|**3          [Pa, w in mm]

for the mean deflection w of the mode-1 shape J0(j01 r / R).

Usage:  python closed_form_coefficients.py [V0_volts]
"""
import math
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import j1, jn_zeros

EPS0 = 8.854e-12


def coefficients(mu=27460.0, Jm=50.0, lam=1.8, H0=0.5e-3, R=0.05, eps_r=4.7, V0=9000.0):
    eps = eps_r * EPS0
    j01 = jn_zeros(0, 1)[0]
    I0 = 2 * lam**2 + lam**-4 - 3
    s = I0 / Jm
    beta1 = 1 - lam**-6
    beta2 = lam**-8
    A0 = math.pi * R**2
    abar = 2 * j1(j01) / j01                      # mean of J0(j01 r/R) over the disc
    J2 = math.pi * j01**2 * j1(j01) ** 2          # int |grad Phi1|^2 dA0  (units of m^0 for this scaling)
    q = quad(lambda t: j1(j01 * t) ** 4 * t, 0, 1)[0]
    J4 = (j01 / R) ** 4 * 2 * math.pi * R**2 * q  # int |grad Phi1|^4 dA0
    G = beta2 / (1 - s) + beta1**2 / (2 * Jm * (1 - s) ** 2)

    # N = sigma0 * h, with h = H0/lam^2 and sigma0 the equibiaxial Cauchy stress of the pre-stretched Gent state
    h = H0 / lam**2
    sigma0 = mu * (lam**2 - lam**-4) / (1 - s)
    N = sigma0 * h
    Vc = math.sqrt(N * h / eps)

    KG = mu * H0 * beta1 * J2 / (1 - s)
    KE = eps * V0**2 * lam**2 * J2 / H0
    K1 = KE - KG
    K3_gent = 2 * mu * H0 * G * J4
    K3E = eps * V0**2 * (A0 * J4 - J2**2) / (2 * H0 * A0)

    scale_a = 1.0 / (A0 * abar**2) * 1e-3      # N/m (K1) -> Pa/mm for the mean deflection
    scale_b = 1.0 / (A0 * abar**4) * 1e-9
    out = dict(
        Vc=Vc, N=N, sigma0=sigma0, h=h, abar=abar, K1=K1, K3_gent=K3_gent, K3E=K3E,
        a_Pa_per_mm=K1 * scale_a,
        b_frozen_Pa_per_mm3=K3_gent * scale_b,          # Gent part only, i.e. without the Maxwell cubic term K3E
        b_full_Pa_per_mm3=(K3_gent + K3E) * scale_b,    # with the Maxwell cubic term K3E
    )
    return out


if __name__ == "__main__":
    V0 = float(sys.argv[1]) if len(sys.argv) > 1 else 9000.0
    c = coefficients(V0=V0)
    print(f"V0 = {V0:.0f} V   V0/Vc = {V0 / c['Vc']:.4f}")
    for k in ("Vc", "N", "sigma0", "h", "abar", "K1", "K3_gent", "K3E",
              "a_Pa_per_mm", "b_frozen_Pa_per_mm3", "b_full_Pa_per_mm3"):  # 'frozen' = without K3E
        print(f"  {k:22s} {c[k]:.6g}")
    a = c["a_Pa_per_mm"]
    for label, b in (("with K3E   ", c["b_full_Pa_per_mm3"]), ("Gent part  ", c["b_frozen_Pa_per_mm3"])):
        wf = math.sqrt(a / (3 * b))
        print(f"  cubic law ({label}): fold at |w| = {wf:.1f} mm, p_max = {2 * a * wf / 3:.0f} Pa; zero at {math.sqrt(a / b):.1f} mm")
