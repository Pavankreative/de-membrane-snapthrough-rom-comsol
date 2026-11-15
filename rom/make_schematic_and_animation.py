"""Schematic figures and an animation built from the FE tables in data/.

fig5_model_schematic.png   : schematic of the model (not to scale)
fig6_workflow.png          : workflow and status of each step
branch_trace.gif           : the centre-deflection-controlled FE run, traced point by point.
                             Curves are the FE data. The dome profile is a schematic family
                             w = w0 (1 - (r/R)^m) with m chosen so that its area-weighted mean
                             equals the FE mean deflection; it is not an FE field.
Usage: python make_schematic_and_animation.py
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

root = Path(__file__).resolve().parents[1]
out = root / "figures"
plt.rcParams.update({"font.size": 10, "figure.dpi": 140})
BLUE, ORANGE, GREEN, GREY, RED = "#0072B2", "#D55E00", "#009E73", "#7f7f7f", "#b03a2e"


# ---------------------------------------------------------------- schematic
fig, ax = plt.subplots(figsize=(10, 4.2))
ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 42)
# left: pre-stretch
ax.text(4, 39, "1  Pre-stretch (rim displaced radially)", fontsize=10, weight="bold")
ax.plot([14, 36], [31, 31], color=BLUE, lw=5, solid_capstyle="butt")
ax.text(25, 34, "reference  R = 50 mm,  H0 = 0.5 mm", ha="center", fontsize=8.5)
ax.annotate("", xy=(25, 21), xytext=(25, 28), arrowprops=dict(arrowstyle="-|>", color=GREY))
ax.text(27, 24, "lambda_p = 1.8", fontsize=8.5, color=GREY)
ax.plot([6, 44], [16, 16], color=BLUE, lw=2, solid_capstyle="butt")
ax.text(25, 11.5, "stretched  Rs = 90 mm,  h = 0.154 mm", ha="center", fontsize=8.5)
# right: loaded membrane
ax.text(54, 39, "2  Voltage V0 and holding pressure p", fontsize=10, weight="bold")
x = np.linspace(56, 96, 100)
w = 9 * (1 - ((x - 76) / 20) ** 2)
ax.plot([56, 96], [28, 28], color=GREY, lw=1, ls=":")
ax.plot(x, 28 - w, color=ORANGE, lw=2.5)
for xx in (56, 96):
    ax.plot([xx, xx], [25, 31], color="k", lw=3)
ax.text(76, 31.5, "clamped rim", ha="center", fontsize=8.5)
for xx in np.linspace(62, 90, 5):
    ax.add_patch(FancyArrowPatch((xx, 23.5), (xx, 20.5), arrowstyle="-|>", mutation_scale=9, color=GREEN))
ax.text(76, 17, "transverse load p (holding pressure)", ha="center", color=GREEN, fontsize=8.5)
ax.text(76, 24.8, "deflection w (schematic)", ha="center", color=ORANGE, fontsize=8, alpha=0)
# bottom: statements
ax.text(4, 6, "Voltage adds an in-plane Maxwell stress  -eps V0^2 / h^2.  The net tension N_eff = N - eps V0^2 / h\n"
        "reaches zero at the critical voltage Vc = 7295.6 V; above it the flat state is unstable (V0 = 1.234 Vc here).",
        fontsize=8.8, va="center")
fig.suptitle("Model schematic (not to scale)", fontsize=11)
fig.tight_layout(); fig.savefig(out / "fig5_model_schematic.png"); plt.close(fig)


# ---------------------------------------------------------------- workflow
fig, ax = plt.subplots(figsize=(10, 4.6))
ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 5)


def box(x, y, w, h, text, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12", fc=color, ec="k", lw=0.8, alpha=0.9))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.8)


def arrow(p, q):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=12, color="k", lw=1))


ok, bad, todo = "#c8e6c9", "#ffe0b2", "#eeeeee"
box(0.2, 3.6, 2.6, 1.0, "Gent + Maxwell free energy\n(docs/theory.md)", ok)
box(3.5, 3.6, 2.9, 1.0, "Closed-form reduced model\nVc, K1, K3, Melnikov thresholds", ok)
box(7.2, 3.6, 2.6, 1.0, "ROM numerics\nthresholds, basins, Lyapunov", ok)
box(0.2, 1.9, 2.6, 1.0, "COMSOL membrane + shell\nVHB 4910 constants", ok)
box(3.5, 1.9, 2.9, 1.0, "Checks: Vc, stress, frequencies,\ngrowth rate  (all agree)", ok)
box(7.2, 1.9, 2.6, 1.0, "Pressure fold, K1\nagree;  K3 does not", bad)
box(0.2, 0.2, 2.6, 1.0, "Buckled wells (FE)\nnot found on the branch", bad)
box(3.5, 0.2, 2.9, 1.0, "Time-domain FE vs ROM\nnot possible (short-wave noise)", bad)
box(7.2, 0.2, 2.6, 1.0, "Azimuthal 3D model,\ncalibration: next", todo)
arrow((2.8, 4.1), (3.5, 4.1)); arrow((6.4, 4.1), (7.2, 4.1))
arrow((1.5, 3.6), (1.5, 2.9)); arrow((2.8, 2.4), (3.5, 2.4)); arrow((6.4, 2.4), (7.2, 2.4))
arrow((4.9, 3.6), (4.9, 2.9)); arrow((8.5, 1.9), (8.5, 1.2)); arrow((1.5, 1.9), (1.5, 1.2)); arrow((4.9, 1.9), (4.9, 1.2))
ax.text(0.2, 4.85, "green: done and agrees    orange: partial or does not agree    grey: planned", fontsize=8.5, color="k")
fig.tight_layout(); fig.savefig(out / "fig6_workflow.png"); plt.close(fig)

# ---------------------------------------------------------------- animation
with open(root / "data" / "membrane_frozen_centre_control.csv") as f:
    rows = list(csv.DictReader(f))
c = -np.array([float(r["dh_centre_mm"]) for r in rows])
m = -np.array([float(r["mean_deflection_mm"]) for r in rows])
p = np.array([float(r["pu_Pa"]) for r in rows])
idx = list(range(0, len(c), 4)) + [len(c) - 1]
idx = sorted(set(idx))
Rs = 90.0
fig = plt.figure(figsize=(10, 4.2))
g = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.25, 1])
axp = fig.add_subplot(g[:, 0]); axs = fig.add_subplot(g[0, 1]); axm = fig.add_subplot(g[1, 1])
axp.plot(c, p, color="0.85", lw=1.5)
(line,) = axp.plot([], [], color=BLUE, lw=2.2); (dot,) = axp.plot([], [], "o", color=ORANGE)
axp.set(xlim=(0, 26), ylim=(0, 135), xlabel="centre deflection (mm)", ylabel="holding pressure (Pa)",
        title="FE branch, membrane, frozen thickness, V0/Vc = 1.234")
axp.grid(alpha=0.25); axp.spines[["top", "right"]].set_visible(False)
txt = axp.text(0.98, 0.05, "", transform=axp.transAxes, ha="right", fontsize=9)
axs.set_aspect("equal"); axs.axis("off"); axs.set_xlim(-95, 95); axs.set_ylim(-30, 6)
r = np.linspace(-1, 1, 200)
(prof,) = axs.plot([], [], color=ORANGE, lw=2.4); axs.plot([-Rs, Rs], [0, 0], color="0.8", lw=1)
axs.set_title("schematic profile (not an FE field)", fontsize=9)
axm.plot(c, m, color="0.85", lw=1.5); (ml,) = axm.plot([], [], color=BLUE, lw=2)
axm.set(xlim=(0, 26), ylim=(0, 13), xlabel="centre deflection (mm)", ylabel="mean deflection (mm)")
axm.grid(alpha=0.25); axm.spines[["top", "right"]].set_visible(False)
fig.tight_layout()


def update(k):
    i = idx[k]
    line.set_data(c[: i + 1], p[: i + 1]); dot.set_data([c[i]], [p[i]])
    ratio = m[i] / c[i] if c[i] > 0 else 0.5
    ex = 2 * ratio / max(1 - ratio, 1e-6)
    prof.set_data(Rs * r, -c[i] * (1 - np.abs(r) ** ex))
    ml.set_data(c[: i + 1], m[: i + 1])
    txt.set_text(f"centre {c[i]:.1f} mm   mean {m[i]:.2f} mm   p = {p[i]:.1f} Pa")
    return line, dot, prof, ml, txt


ani = FuncAnimation(fig, update, frames=len(idx), interval=70, blit=False)
ani.save(out / "branch_trace.gif", writer=PillowWriter(fps=14))
plt.close(fig)
print("written:", [x.name for x in sorted(out.glob("*")) if x.name.startswith(("fig5", "fig6", "branch"))])
