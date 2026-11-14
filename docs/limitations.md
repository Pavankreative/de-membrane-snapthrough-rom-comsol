# Limitations, open questions and next steps

This file lists what the work does **not** show. Read it before relying on any number in this repository.

## Status of the claims

| Claim | Status | Evidence |
|---|---|---|
| Pre-stress, Maxwell tension reduction, frequency law, deflection response match closed forms | Verified in COMSOL | `data/membrane_checks_summary.csv` |
| Critical voltage $V_c=7295.6$ V reproduced | Verified in COMSOL (to 0.006%) | same file; the Maxwell stress is imposed in the same closed form, so this checks pre-stress, thickness bookkeeping and the solver, not the electromechanics |
| First-mode growth rate of the unstable flat state matches membrane theory | Verified (174.5 vs 174.5 1/s) | `data/shell_unstable_growth_rates_V8755.csv` |
| Linear stiffness $K_1$ of the ROM | Agrees within 4.5% | `rom/closed_form_coefficients.py` vs `data/membrane_frozen_mean_control.csv` |
| Pressure fold | Found: 122.8 Pa at 7.7 mm (membrane, 1.234 $V_c$); about 180 Pa at 13-14 mm (shell, 1.20 $V_c$) | `data/` |
| Closed-form cubic coefficient $K_3$ | **Wrong by a large factor** (3.6 to 83 times too small) | `rom/fit_cubic.py` |
| Buckled wells predicted by the cubic ROM | **Not confirmed** | below |
| Stability of the wells | **Not established** | below |
| Time-domain escape under forcing (FE) | **Not done** | below |
| Azimuthal pass (the $\rho$ coordinate of the ROM) | **Not tested** (models are axisymmetric) | |
| Coupled fluid-structure response | **Not done** | |
| Calibration of the ROM parameters against FE | **Not done** (the target matrix has no solver outputs) | |

## 1. No well on the axisymmetric branch

With a frozen thickness, the membrane branch passes the fold and then saturates: the mean deflection stops at about $11.7$ mm while the centre deflection goes on to $25$ mm, and the holding pressure falls to $62.5$ Pa and slows down. It never reaches zero (`figures/fig2_centre_control_saturation.png`). The cubic model puts a zero-load point at $13.3$ mm mean deflection; that is an extrapolation of a cubic fit to 0-10 mm and is not a point of the FE branch (`figures/fig4_energy_landscape.png`). Whether the true buckled states are non-axisymmetric or localised is open; that is the reason for the planned 3D (azimuthal) model.

## 2. The thickness law matters

* With the **deformed thickness** in the Maxwell stress, $\sigma_M=-\varepsilon V^2/h_\text{def}^2$, the membrane solver stops at a mean deflection of 2.2 mm and 46.8 Pa on both a 20-element and a 100-element mesh ($46.816$ and $46.830$ Pa at the last step). The probable cause is a local thinning limit, since the electric stress rises as $1/h_\text{def}^2$. This was not isolated.
* With a **frozen thickness** the branch passes the fold. All the pressure-law results above are for this assumption.
* The shell model freezes the thickness. Whether its layered shell thins with stretch was not verified; the thickness variable reports the input value and a section-force check was inconclusive.

## 3. Shell and membrane disagree on the cubic coefficient

FE cubic coefficients (frozen thickness): $0.0378$ Pa/mm³ for the shell at $V_0/V_c=1.20$ and $0.134$ for the membrane at $1.234$, a factor 3.5, although the closed form changes by under 2% between the two voltages. The cause (bending, the definition of the mean, the fit range) was not identified.

## 4. Short-wave instability and the time-domain run

For $V_0>V_c$ the net tension is negative and all axisymmetric modes are unstable (`docs/theory.md`, section 6). The growth rate rises with mode number: $190\ \mathrm{s^{-1}}$ for mode 1 at $V_0=9000$ V, but about $1.8\times10^3$ at $n=8$, $4.5\times10^3$ at $n=20$ and $2.3\times10^4$ at $n=100$ on a 100-element membrane mesh.

A time-domain run from rest at zero load (frozen thickness, 100 elements, start at a mean deflection of 0.1 mm) failed at $0.69$ ms with a jagged deformed shape and a von Mises stress of $3\times10^6$ N/m². A growth by a factor of $10^6$ in $0.7$ ms is of the order needed to turn solver tolerance into a visible signal. Mass-proportional damping cannot remove a negative-stiffness instability. **No FE-ROM dynamic comparison is reported.** The one-mode ROM describes the escape only if the short-wave band is suppressed or not seeded.

FE estimates for the shell, from the spectrum at $V_0/V_c=1.20$ (22 unstable pairs returned): an effective bending stiffness $D_\text{eff}\approx0.98\times10^{-6}$ N m, an unstable band ending near $n_c\approx73$ (wavelength about 2.5 mm) and a maximum growth rate of about $8\times10^3\ \mathrm{s^{-1}}$, about 48 times the first-mode rate. These estimates assume the shell does not thin.

## 5. Inputs

* The material constants ($\mu=27.46$ kPa, $J_m=50$, $\varepsilon_r=4.7$, $\rho=960$ kg/m³) are values in the published range for VHB 4910, not a fit to one film. The source for $\varepsilon_r=4.7$ is not verified.
* No experiments were done.

## 6. Method

* The Maxwell stress is imposed as an in-plane Cauchy stress; the electric field is not solved. The full-electrostatics route was not realised. This was also the case in an earlier ANSYS attempt (Student licence limits), whose placeholder constants are not used here.
* Meshes: shell with 1000 edge elements along the radius; membrane with 20 and 100 elements.

## Next steps

1. A 3D (or harmonic-expansion) membrane model to test the azimuthal pass and look for non-axisymmetric buckled states.
2. Calibration of $K_3$ and the other ROM parameters against FE, with the target matrix filled with computed values.
3. A dynamic test in a model where the short-wave band is suppressed, compared with `rom/rom_release_prediction.py`.
4. Aerodynamic loading as a prescribed projected pressure (not real CFD) after the structural part is settled.
