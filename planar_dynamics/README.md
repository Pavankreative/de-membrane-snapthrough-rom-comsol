# Supplement: planar subsystem and coupled-model dynamics (Tables 3-5, Figs 5-7)

Reference parameters (nondimensional): delta = 0.08, delta2 = 0.03, f = 0.12, g1 = 0.50, g3 = 0.08, sigma_n = 1.593, omega = 1.72, gamma = 3.0, K_nn = 8.0.
Requires NumPy, SciPy, Matplotlib and gcc (`rk.c` compiles automatically into `librk.so` on first import of `crk.py`). Run from this folder.

| File | Purpose |
|---|---|
| `melnikov_envelope.py` | validity conditions and first-order Melnikov threshold; `python3 melnikov_envelope.py` gives A_c = 0.20439 for the reference case (Table 3, first-order column) |
| `rk.c`, `crk.py` | RK4 integrator and Benettin Lyapunov routines (OpenMP) |
| `run_more.py` | Table 4: Lyapunov scans, run-length check, uncertainty exponents -> `results/lyap.npz`, `nlyap.npz` |
| `run_b3.py`, `run_basins.py` | basin maps -> `results/basins3.npz` (Fig. 5) |
| `run_AN.py`, `run_AN2.py` | normal-instability threshold, M2 integral, Poincare sections -> `results/poin_*.npz` (Fig. 6) |
| `mkfigs.py` | Figs 5-7 from `results/` |
| `verify_lyap.py` | recomputes selected Table 4 entries; they match the stored results bit for bit |

## Not included
* The script of the direct (numerical) manifold computation behind the "direct onset" column of Table 3.
* The exact trajectory counts of Table 5 are unresolved (see manuscript audit).
