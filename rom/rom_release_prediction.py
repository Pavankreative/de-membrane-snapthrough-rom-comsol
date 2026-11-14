"""Single-mode ROM release test: a PREDICTION, not validated by FE.

Mean-deflection ROM at zero load:   M y'' + C y' = a*y - b*y**3 ,   y = |mean deflection| (> 0)

 a, b : FE linear slope and FE cubic coefficient of the frozen-thickness membrane at V0 = 9000 V
 M    : effective mass per reference area, from the FE first-mode growth rate at 8755 V
        (s1 = 174.5 1/s) and the linear slope a(8755 V) from the closed form scaling
 C/M  : Rayleigh mass-damping rate of the COMSOL models, 14.4 1/s

A time-domain FE run of the membrane model failed at 0.69 ms from short-wave noise growth
(docs/limitations.md), so these numbers have NOT been compared with FE.
"""
import numpy as np
from scipy.integrate import solve_ivp

Vc = 7295.57
a = 23.968e3          # Pa/m   (FE slope, 9000 V)
b = 0.134e9           # Pa/m^3 (FE cubic coefficient, 0-10 mm fit)
a8755 = a * ((8755 / Vc) ** 2 - 1) / ((9000 / Vc) ** 2 - 1)
M = a8755 / 174.5**2  # kg/m^2
c = 14.4              # 1/s


def rhs(t, z):
    y, v = z
    return [v, (a * y - b * y**3) / M - c * v]


def run(y0, tmax=0.06):
    sol = solve_ivp(rhs, [0, tmax], [y0, 0.0], max_step=2e-5, rtol=1e-10, atol=1e-13)
    return sol.t, sol.y[0]


if __name__ == "__main__":
    print(f"M = {M:.4f} kg/m^2, s0(9000 V) = {np.sqrt(a / M):.1f} 1/s, b/M = {b / M:.3g}")
    for y0 in (0.1e-3, 2e-3, 3e-3):
        t, y = run(y0)
        hit = {L: (t[np.argmax(y >= L * 1e-3)] * 1e3 if y.max() >= L * 1e-3 else None) for L in (1, 2, 5, 8, 10)}
        print(f"y0 = {y0 * 1e3:.1f} mm -> time (ms) to reach",
              {k: (round(v, 2) if v else None) for k, v in hit.items()},
              f"| peak {y.max() * 1e3:.1f} mm at {t[y.argmax()] * 1e3:.1f} ms")
