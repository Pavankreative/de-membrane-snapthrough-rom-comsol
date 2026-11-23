# Execution and evidence record — 21 November 2026

Python 3.12 on Linux. All five public scripts completed with exit code 0 in an isolated copy of the repository, using `MPLBACKEND=Agg` for figures:

- `python rom/closed_form_coefficients.py 9000`
- `python rom/fit_cubic.py`
- `python rom/rom_release_prediction.py`
- `python rom/make_figures.py`
- `python rom/make_schematic_and_animation.py`

The last two regenerated PNG/GIF outputs; existing tracked artwork was not replaced merely because of rendering-library differences. This is a script execution check, not a COMSOL rerun or experimental validation.

Reproduced closed-form values: pre-stress 93004.9 Pa, critical voltage 7295.57 V, linear slope 25.0493 Pa/mm and full cubic coefficient 0.0106923 Pa/mm³. Fitting 219 archived FE points over 0–10 mm gives 23.908 Pa/mm, 0.1341 Pa/mm³ and RMS residual 0.075 Pa. The archived fold is 122.80 Pa at 7.70 mm.

The public COMSOL check table reports 41.858 Hz at zero voltage.

The Gent/Maxwell derivation and Melnikov criterion remain documented in the theory notes. The manifold/basin scripts underlying the reported Melnikov numerical checks are not public here. Axisymmetric FE cannot verify the azimuthal mode; nonlinear cubic disagreement, absent zero-pressure wells and the failed FE time-domain comparison remain open limitations. No manuscript acceptance or peer-review claim is made.
