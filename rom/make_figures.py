"""Regenerate the figures in figures/ from the CSV files in data/.  Usage: python make_figures.py"""
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[1]
data, out = root / "data", root / "figures"
out.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 140})
C = {"fe": "#0072B2", "fit": "#D55E00", "def": "#009E73", "rom": "#7f7f7f", "ln": "#CC79A7"}


def col(name, key, sign=1.0):
    with open(data / name) as f:
        return np.array([sign * float(r[key]) for r in csv.DictReader(f)])


# FE data
w_fz, p_fz = col("membrane_frozen_mean_control.csv", "dh_mm", -1), col("membrane_frozen_mean_control.csv", "pu_Pa")
w_df, p_df = col("membrane_deformed_thickness_mean_control.csv", "dh_mm", -1), col("membrane_deformed_thickness_mean_control.csv", "pu_Pa")
wc_c, wc_m, p_c = (col("membrane_frozen_centre_control.csv", "dh_centre_mm", -1),
                   col("membrane_frozen_centre_control.csv", "mean_deflection_mm", -1),
                   col("membrane_frozen_centre_control.csv", "pu_Pa"))



# cubic fit to 0-10 mm (frozen)
m = (w_fz > 0) & (w_fz <= 10)
a, mb = np.linalg.lstsq(np.c_[w_fz[m], w_fz[m] ** 3], p_fz[m], rcond=None)[0]
b = -mb
a_rom, b_rom_full, b_rom_gent = 25.0493, 0.0106923, 0.00698391   # rom/closed_form_coefficients.py, V0 = 9000 V
wf = np.linspace(0, 16, 400)

# Figure 1: pressure law
fig, ax = plt.subplots(1, 2, figsize=(10, 4.0))
ax[0].plot(w_fz, p_fz, color=C["fe"], lw=2, label="FE, frozen thickness")
ax[0].plot(wf, a * wf - b * wf**3, "--", color=C["fit"], lw=1.6, label=f"cubic fit to 0-10 mm\n(a={a:.2f}, b={b:.3f})")
ax[0].plot(w_df, p_df, color=C["def"], lw=2, label="FE, deformed thickness")
ax[0].plot(w_df[-1], p_df[-1], "o", color=C["def"]); ax[0].annotate("solver stops\n2.2 mm, 46.8 Pa", (w_df[-1], p_df[-1]), (3.2, 22), arrowprops=dict(arrowstyle="->", color=C["def"]), color=C["def"])
ax[0].plot(w_fz[p_fz.argmax()], p_fz.max(), "o", color=C["fe"]); ax[0].annotate("fold\n122.8 Pa, 7.7 mm", (7.7, 122.8), (2.0, 126), color=C["fe"])
ax[0].set(xlabel="mean deflection |w| (mm)", ylabel="holding pressure p (Pa)", xlim=(0, 16), ylim=(0, 145),
          title="Membrane model, V0/Vc = 1.234")
ax[0].legend(frameon=False, loc="upper right", fontsize=8.5)
ax[1].plot(w_fz, p_fz, color=C["fe"], lw=2, label="FE, frozen thickness")
ax[1].plot(wf, a_rom * wf - b_rom_full * wf**3, color=C["rom"], lw=1.8, label="closed-form ROM, with K3E")
ax[1].plot(wf, a_rom * wf - b_rom_gent * wf**3, ":", color=C["rom"], lw=1.8, label="closed-form ROM, Gent part only")
ax[1].set(xlabel="mean deflection |w| (mm)", ylabel="p (Pa)", xlim=(0, 16), ylim=(0, 340), title="Closed-form ROM against FE")
ax[1].legend(frameon=False, loc="upper left", fontsize=8.5)
fig.tight_layout(); fig.savefig(out / "fig1_pressure_vs_deflection.png"); plt.close(fig)

# Figure 2: saturation under centre-deflection control
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(wc_c, wc_m, color=C["fe"], lw=2, label="mean deflection (deformed-area mean)")
ax[0].plot(wc_c, 0.5 * wc_c, "--", color=C["rom"], label="0.5 x centre (paraboloid)")
ax[0].set(xlabel="centre deflection (mm)", ylabel="mean deflection |w| (mm)", title="Mean saturates, apex keeps moving"); ax[0].legend(frameon=False, fontsize=8.5)
ax[1].plot(wc_c, p_c, color=C["fe"], lw=2)
ax[1].set(xlabel="centre deflection (mm)", ylabel="holding pressure (Pa)", ylim=(0, 135), title="Pressure never reaches zero (to 25 mm)")
ax[1].axhline(0, color="k", lw=0.8)
fig.tight_layout(); fig.savefig(out / "fig2_centre_control_saturation.png"); plt.close(fig)

# Figure 3: unstable ladder
n, fe, th = [], [], []
with open(data / "shell_unstable_growth_rates_V8755.csv") as f:
    for r in csv.DictReader(f):
        n.append(int(r["mode_n"])); fe.append(float(r["fe_growth_rate_per_s"])); th.append(float(r["membrane_theory_growth_rate_per_s"]))
n, fe, th = map(np.array, (n, fe, th))
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(n, th, "-", color=C["rom"], label="membrane theory  c j0n / Rs"); ax[0].plot(n, fe, "o", color=C["fe"], label="FE shell")
ax[0].set(xlabel="mode number n", ylabel="growth rate (1/s)", title="Unstable flat state, V0/Vc = 1.20"); ax[0].legend(frameon=False)
ax[1].bar(range(len(n)), (fe / th - 1) * 100, color=C["fe"]); ax[1].set_xticks(range(len(n)), n)
ax[1].set(xlabel="mode number n", ylabel="FE / theory - 1 (%)", title="Departure grows as k^2 (bending)")
fig.tight_layout(); fig.savefig(out / "fig3_unstable_ladder.png"); plt.close(fig)

# Figure 4: energy per area along the branch, U(w) = -int p dw (trapezoid for FE; analytic for the cubic)
y = np.linspace(0, 16, 800)
U_fit = -(0.5 * a * y**2 - 0.25 * b * y**4)        # U(y) = -int p dy  (p = holding pressure; U falls as y grows)
U_fe = -np.concatenate([[0], np.cumsum(0.5 * (p_fz[1:] + p_fz[:-1]) * np.diff(w_fz))])
fig, ax = plt.subplots(figsize=(5.6, 4.0))
ax.plot(w_fz, U_fe, color=C["fe"], lw=2, label="FE, frozen thickness (integrated)")
ax.plot(y, U_fit, "--", color=C["fit"], label="cubic ROM (extended)")
ax.axvline(math.sqrt(a / b), color=C["fit"], lw=0.8, ls=":"); ax.text(math.sqrt(a / b) + 0.1, -700, "cubic zero-load\npoint, 13.3 mm\n(extrapolation)", color=C["fit"], fontsize=8.5)
ax.axvspan(w_fz.max(), 16, color="0.9", zorder=0)
ax.text(w_fz.max() + 0.15, -150, "not reached by FE\n(branch ends at 11.75 mm,\np = 60 Pa, not zero)", fontsize=8, color="0.3")
ax.set(xlabel="mean deflection |w| (mm)", ylabel="U(w) per area (Pa mm)", title="Energy along the axisymmetric branch", xlim=(0, 16))
ax.legend(frameon=False, loc="lower left", fontsize=8.5)
fig.tight_layout(); fig.savefig(out / "fig4_energy_landscape.png"); plt.close(fig)
print("figures written to", out)
