# Dielectric-elastomer membrane snap-through: reduced-order model and COMSOL validation

A pre-stretched dielectric-elastomer (DE) membrane under voltage becomes unstable when the Maxwell stress cancels its in-plane tension. This project builds an analytical reduced-order model (ROM) of the snap-through, checks it against axisymmetric finite-element (FE) models in COMSOL Multiphysics 6.4, and records plainly where the reduced model holds and where it does not.

**Status:** research in progress, by one author, not peer reviewed. A manuscript is in preparation. All FE results use material constants from the published range for VHB 4910, not a fit to one film, and no experiments were done.

## What the project contains

* **Analytical part** (`docs/theory.md`): Gent free energy with a Maxwell stress; the critical voltage; a first-mode reduced potential $V(x)=-\tfrac12K_1x^2+\tfrac14K_3x^4$ with closed-form coefficients; an O(2)-symmetric extension with one azimuthal mode and its nine critical points; a Melnikov criterion for homoclinic tangling with two loop thresholds; where the one-mode reduction stops being valid.
* **Numerical part** (`docs/numerical_methods.md`): Melnikov quadrature, manifold computation, crossing thresholds, basin erosion and Lyapunov estimates on the reduced model, with numbers.
* **FE validation** (`docs/comsol_validation.md`): membrane and shell models, checks against closed forms, the spectrum of the unstable flat state, and pressure-deflection curves through a fold using a deflection-controlled continuation.
* **Code** (`rom/`): closed-form coefficients, cubic fits to the FE data, a ROM release-test prediction, and the figure script. Everything in `rom/` runs from the files in `data/`.

## Main results

| Result | Value | Status |
|---|---|---|
| Critical voltage $V_c$ (closed form) | 7295.6 V | reproduced in FE to 0.006% |
| Pre-stress, Maxwell tension reduction, frequency law, linear deflection response | 33 204 Pa at 5850 V; 41.858 Hz at 0 V | match closed forms to $\le0.11\%$ (stress) and $\le0.015\%$ (frequency) away from $V_c$ |
| First-mode growth rate of the unstable flat state, $V_0/V_c=1.20$ | 174.5 s⁻¹ | FE shell = membrane theory |
| Unstable flat state | 22 unstable pairs, growth rates to $4.7\times10^3$ s⁻¹ | the slow mode is not the fastest-growing direction |
| Linear stiffness $K_1$ of the ROM | within 4.5% of FE | agrees |
| Pressure fold, membrane, $V_0/V_c=1.234$ | 122.8 Pa at a mean deflection of 7.7 mm | FE, frozen thickness |
| Cubic coefficient $K_3$ of the ROM | closed form 3.6 to 83 times too small | **does not agree; must be calibrated** |
| Buckled wells of the cubic ROM | not found on the axisymmetric FE branch | **not confirmed** |
| Time-domain FE vs ROM | not possible; solver fails at 0.69 ms | **not done** |

![Pressure against deflection](figures/fig1_pressure_vs_deflection.png)

*Left: holding pressure against mean deflection in the membrane model. The frozen-thickness run passes the fold at 122.8 Pa and follows a cubic law to 11.7 mm; with the deformed thickness the solver stops at 2.2 mm and 46.8 Pa. Right: the closed-form ROM against FE. The slope is right, the cubic term is not.*

![Saturation under centre control](figures/fig2_centre_control_saturation.png)

*Beyond the fold the mean deflection saturates at about 11.7 mm while the apex keeps moving, and the pressure stays near 62 Pa. No zero-pressure equilibrium (no well) is found on this axisymmetric branch up to 25 mm.*

## What this does and does not show

It shows that the first-mode reduction is quantitatively right for the stiffness, the critical voltage and the first growth rate, and that the closed-form cubic coefficient is not. It also shows that above the critical voltage the flat state has many faster unstable directions than the one the ROM keeps.

It does **not** show the buckled wells, their stability, the time-domain escape under forcing, the azimuthal pass in FE, the coupled fluid-structure response, or a calibrated ROM. `docs/limitations.md` lists each with the evidence, the open questions and the next steps (a 3D or harmonic-expansion model of the azimuthal modes is the next piece).

## Repository layout

```
README.md
docs/
  theory.md              analytical model and criterion
  numerical_methods.md   ROM numerics, with results
  comsol_validation.md   FE models, checks, pressure-deflection results, practical notes
  limitations.md         what is not shown, open questions, next steps
data/                    tables exported from COMSOL (see data/README.md)
rom/
  closed_form_coefficients.py   K1, K3, a, b from the formulas of docs/theory.md
  fit_cubic.py                  fit p = a|w| - b|w|^3 to FE data, fold and zero
  rom_release_prediction.py     single-mode release test (prediction only, not validated)
  make_figures.py               regenerates figures/ from data/
figures/
models/                  COMSOL files (see models/README.md)
```

## Reproduce

```bash
pip install -r requirements.txt
python rom/closed_form_coefficients.py 9000    # V0 in volts; prints Vc, K1, K3, a, b
python rom/fit_cubic.py                        # fit to data/membrane_frozen_mean_control.csv
python rom/rom_release_prediction.py
python rom/make_figures.py
```

`closed_form_coefficients.py 9000` returns $V_c=7295.6$ V, $a=25.05$ Pa/mm and $b=0.0107$ (with the Maxwell cubic term) or $0.0070$ (without). `fit_cubic.py` returns $a=23.91$ Pa/mm and $b=0.134$ for the 0-10 mm fit, with an RMS residual of 0.075 Pa.

## Tools

COMSOL Multiphysics 6.4 (Membrane, Shell and Global ODE interfaces, eigenfrequency, stationary and time-dependent studies); Python with NumPy, SciPy and Matplotlib.

## Author and licence

Pavan Kumar Karsinapu, ORCID 0009-0006-8040-5244. Code and text under the MIT licence (`LICENSE`). Source data for material constants are to be cited in `docs/theory.md` (section 7) before this is relied on.
