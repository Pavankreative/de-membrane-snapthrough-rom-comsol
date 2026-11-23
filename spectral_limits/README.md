# Supplement: flat-state spectrum, Table 7 and exclusion criterion

Companion to *Membrane_Spectral_Limits* (Proc. R. Soc. A submission). Python 3, NumPy, SciPy only.
Run every script from this folder.

| Script | Reproduces | Check value |
|---|---|---|
| `fit_deff.py` | D_eff and offset delta from the stable-side axisymmetric spectra | pinned 6.738e-7 N m, delta -5.1e-6; clamped 6.719e-7, +1.06e-6 |
| `fit_xc.py` | stable-side quadratic extrapolation of crossings (COMSOL column of Table 7) | pinned (0,1) 2.196e-5, (0,2) 1.659e-4; clamped 9.167e-5, 2.741e-4 |
| `table7.py` -> `table7_crossings.csv` | Table 7 (plate cases A, B, fitted D_eff = 6.73e-7; COMSOL column) | |
| `exclusion_criterion.py` | Sec. 5: Delta >= 0.19048 vs (1,1) crossings | factors 2,418 (pinned), 1,248 (clamped) |
| `plate_test.py` | shell-bending check at lambda = 1 | D = 1.148e-6 vs 1.144e-6 N m |
| `plate_dispersion.py`, `comsol_io.py` | exact prestressed plate dispersion (any n); COMSOL export reader | |

## Data (`data/`)
* `raw/` COMSOL 6.4 Global Evaluation exports, unedited: `2D_pinned_1000` (pinned, 1000 elements; used for pinned results), `2D_pinned_250`, `2D_clamped_250` (clamped results).
* `*_w2.csv` squared frequencies recomputed by `make_w2_csv.py` as w2 = Re[-lambda(lambda + 14.4)], lambda = 2 pi i f. The `w2` column inside the raw files used a different expression and is not used.
* `plate_test/plate_test.csv` result of `DE_plate_test.mph`.
* `diagnostic/` runs **not used for any reported number**: `3D_pinned_5mm` (four-fold symmetric mesh mixes n = 0 and n = 4, no crossing extracted), `2D_pinned_250_n1` and `2D_pinned_250_m0check` (azimuthal-extension runs; n = 1 results disagree with plate theory, unresolved).

## Limits
* All sampled xx are negative: every COMSOL crossing is an extrapolation from the stable side.
* (1,1) and (2,1) crossings were not obtained from COMSOL; the fitted-rigidity values are conditional plate-model predictions.
* The origin of D_eff = 0.59 x the unstretched-plate rigidity is not identified.
* The `.mph` models are not included here.
