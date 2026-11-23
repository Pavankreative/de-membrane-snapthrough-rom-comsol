# Revision 17 dynamics supplement

This directory contains the planar and coupled oscillator computations for *Homoclinic splitting and the spectral limits of a reduced dielectric-elastomer membrane model*. The parameters are dimensionless. The archived trajectories are simulations of the reduced equations. COMSOL exports and the plate calculation are in the repository's separate `spectral_limits/` directory.

## Reproduce

The update ZIP overlays the repository at commit `66fb93af74fcfb68e0f6c59d0ac80257de4aa842`. All files needed for this dynamics workflow are also included in the ZIP. Extract to a writable directory and run the following from `planar_dynamics/`.

```sh
python -m pip install -r requirements-rev17.txt
python verify_rev17.py
python run_rev17.py
python verify_rev17.py /path/to/Membrane_Spectral_Limits_rev17_FINAL.docx
```

The first verification reads the archived results. The full run generates them again in place and verifies the SHA-256 manifest. Use an empty working copy to preserve the reference deposit. The optional DOCX argument compares every displayed body cell of Tables 3–5 with its CSV-derived formatting and compares embedded PNG bytes for Figures 1, 4, 5, 6, 7 and 8. Legacy exploration scripts elsewhere in the base repository use different horizons and event definitions.

## Execution environment

Reference runs used Linux x86-64, CPython 3.12.14, GCC 14.2.0, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8 (Agg, bundled DejaVu Serif). The wrapper uses one OpenBLAS thread and four OpenMP threads. GCC builds `rk.c` with `-O3 -ffp-contract=off -shared -fPIC -fopenmp` and links libm. The reference checksums identify exact bytes for this environment. Differences in compilers, math libraries and floating-point evaluations can change late-time counts in chaotic trajectories. Preserve the case settings when comparing results.

## Results and entry points

| Manuscript | Command or script | Saved evidence |
|---|---|---|
| Table 3, Figure 7 | `direct_manifold_onset.py` | `manifold_splitting_table3.csv` contains the two section coordinates and oriented splitting for each phase. `manifold_phase_zeros.csv` gives brackets and two slope estimates. `manifold_tail_check.csv` and `table3.csv` give refinement and onset results. |
| Figures 5–6 | `run_basin_20p_first_crossing.py` | `basin_labels.npz` contains coordinate axes, eligibility mask and eight binary maps at each step size. `basin_summary.csv` and `uncertainty_fits.csv` give counts and fit inputs. |
| Table 4, Figure 8 | `verify_lyap.py` | `lyapunov_per_start.csv`, `table4.csv`, `poincare_table4.npz`. |
| Table 5, Figures 1 and 4 | `run_table5_trajectories.py` | `table5.csv`, `table5_traces.npz`, `figure4_sweep.csv`. |
| Table 2 and seed comparison | `seed_brackets.py` | `seed_brackets.csv`, `seed_trials.csv`. |
| Plot all six figures | `make_rev17_figures.py` | PNGs in `figures/` read only archived CSV/NPZ files. |
| Published precision | `published_tables.py`, `verify_rev17.py` | Manuscript values are formatted once from full-precision CSVs. |

### Manifold convention

Collocation integrates the nonlinear stable and unstable primary arcs to the turning section at prescribed forcing phase. The central orbit satisfies periodic boundary conditions and its Floquet multipliers bracket unit modulus. The arc endpoints match the saddle displacement at finite time. The signed separation is `s*(xu-xs)`, with s = +1 right and −1 left. The 16-phase grid is the even-index subset of the 32-phase grid. A shared peak phase explains the identical sampled onsets. A continuous local phase maximisation provides a separate sensitivity check. Root slopes were tested above onset, where the two phase roots are distinct.

| Loop | S | First order | Direct 16 / 32 | Difference |
| --- | --- | --- | --- | --- |
| Right | 0.632596 | 0.204387 | 0.204608 / 0.204608 | +0.108% |
| Left | 0.387621 | 0.333558 | 0.333946 / 0.333946 | +0.116% |

### Basin classifications and fits

The grid comprises 200 endpoints-inclusive values per axis, x in [−1.5, 0] and velocity in [−0.75, 0.75]. Eligibility requires x < 0 and v²/2 − x²/2 + x⁴/4 < −0.01. At each accepted RK4 step the code tests x ≥ 0 and retains the first hit through 20 forcing periods. Phase starts at zero. The coupled seed is 0.0001 and the planar seed is zero. Each `caseI_sN` key corresponds to row I of the summary using N steps per period. The mask distinguishes ineligible pixels from eligible zero labels. The uncertainty fit counts horizontal pairs with both endpoints eligible at 1, 2, 4 and 8 pixel offsets. It uses ordinary least squares in log–log coordinates.

| Model | A | Escaped / eligible | Escape % | α | T/100 vs T/200 agreement % |
| --- | --- | --- | --- | --- | --- |
| coupled | 0.1 | 1068 / 21988 | 4.857 | 0.765491 | 100.000 |
| coupled | 0.2044 | 1296 / 21988 | 5.894 | 0.424051 | 100.000 |
| coupled | 0.3 | 3512 / 21988 | 15.972 | 0.098536 | 99.827 |
| coupled | 0.39 | 7555 / 21988 | 34.360 | 0.045731 | 99.072 |
| coupled | 0.5 | 12844 / 21988 | 58.414 | 0.066633 | 98.172 |
| planar | 0.39 | 616 / 21988 | 2.802 | 0.267793 | 100.000 |
| planar | 1.0 | 14381 / 21988 | 65.404 | 0.296226 | 99.973 |
| planar | 1.49 | 21238 / 21988 | 96.589 | 0.596754 | 99.995 |

### Finite-time Lyapunov exponents

Each row uses 30 starts, 100 transient and 600 measured periods, RK4 at T/200, and a two-trajectory separation of 10⁻⁸ renormalised each period. The first eight starts supply the stroboscopic sections.

| Model | A | Median | Range | Positive starts |
| --- | --- | --- | --- | --- |
| Planar | 0.1 | -0.0410 | -0.0418 to -0.0410 | 0/30 |
| Planar | 0.25 | -0.0424 | -0.0424 to -0.0424 | 0/30 |
| Planar | 0.4 | -0.0437 | -0.0462 to -0.0437 | 0/30 |
| Planar | 1.0 | -0.0561 | -0.0570 to -0.0475 | 0/30 |
| Planar | 1.49 | -0.0682 | -0.0682 to -0.0678 | 0/30 |
| Coupled | 0.2044 | -0.0420 | -0.0426 to -0.0217 | 0/30 |
| Coupled | 0.39 | +0.2002 | +0.1790 to +0.2187 | 30/30 |
| Coupled | 0.5 | +0.2191 | +0.2041 to +0.2318 | 30/30 |

### Reference trajectories

DOP853 uses relative tolerance 10⁻⁹ and absolute tolerance 10⁻¹¹. The count is the number of successive sampled displacement products strictly below zero. Except for the explicitly undamped control of duration 300, the sampling grid is 24,001 points over 120 periods. Array `runI` contains time, x and ρ for row I of `table5.csv` (zero-based). The two Figure 1 cases also store velocity. A finite-time crossing bracket depends on the coarse amplitude scan and subsequent local bisection.

| γ, Kₙₙ | Conditions | Outcome |
| --- | --- | --- |
| 0.15, 1.25 | 1, 2, 3 fail | Unforced, damped: final (x, ρ) = (-0.842, 1.395). Unforced, undamped: 74 sign changes in 300 time units. |
| 3.0, 1.25 | 2 fails | A = 0.5: 0 sign changes; max |ρ| = 0.210. A = 1.0: 223 sign changes; max |x| over the last 20 periods = 0.057, final ρ = 1.424. |
| 3.0, 8.0 | All pass | At A = 0, 0.1, 0.2 and 0.3: zero sign changes. At A = 0.4, 0.5 and 1.0: 2, 15, 91 sign changes, with max |ρ| = 0.677, 0.741, 0.900, respectively. |
