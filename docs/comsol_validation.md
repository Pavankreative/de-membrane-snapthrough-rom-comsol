# COMSOL finite-element validation

Two axisymmetric COMSOL Multiphysics 6.4 models of the pre-stretched dielectric-elastomer membrane were built to test the reduced model of `docs/theory.md`: a **Membrane**-interface model and a **Shell**-interface model (with bending). This note describes them, the checks, the results and the pitfalls. Data are in `data/`; figures in `figures/`. The COMSOL files are in `models/` (see its README).

## 1. Models

| Item | Value |
|---|---|
| Geometry | axisymmetric line, reference radius $R=50$ mm, reference thickness $H_0=0.5$ mm |
| Pre-stretch | equibiaxial $\lambda_p=1.8$ by a radial rim displacement; stretched radius $R_s=90$ mm, thickness $h=0.15432$ mm |
| Material | incompressible Gent, $\mu=27.46$ kPa, $J_m=50$, density $960$ kg/m³ |
| Permittivity | $\varepsilon_r=4.7$ (the electric field is not solved) |
| Maxwell stress | in-plane isotropic Cauchy stress $\sigma_M=-\varepsilon_0\varepsilon_rV_0^2/h^2$, with $h$ frozen (`hf`) or deformed (`mbrn.ddef`) |
| Damping | Rayleigh mass damping, $14.4\ \mathrm{s^{-1}}$ |
| Membrane mesh | 20 and 100 elements along the radius |
| Shell | layered hyperelastic material, 1000 edge elements along the radius ($0.09$ mm stretched), geometric nonlinearity in every step |

Derived values: $\sigma_0=93\,005$ Pa, $N=14.353$ N/m, $V_c=7295.6$ V. The voltages used are $V_0=8755$ V ($V_0/V_c=1.20$, shell) and $9000$ V ($1.234$, membrane).

The Maxwell stress is imposed through the Membrane interface's *External Stress*. The electrostatics interface was not coupled, so the checks below test the pre-stress, the thickness bookkeeping and the solver. They are not an independent test of the electromechanical coupling.

### Studies in the membrane file

The numbering follows the `.mph` file; the descriptions are summaries.

| Study | Purpose |
|---|---|
| 1 | Voltage sweep, stationary + eigenfrequency (frequency law) |
| 2 | Pressure sweep $p_0$, first pressure fold attempts |
| 3 | Stationary + time-dependent run (dynamic study) |
| 4 | Voltage sweep 5850-8750 V for the stress and deflection checks (flat constraint off) |
| 5 | Volume-controlled continuation (mean or centre deflection prescribed) |
| 6 | Pressure continuation from a state of Study 5 |
| 7 | Time-dependent release from rest at zero load (failed at 0.69 ms) |

## 2. Checks against closed forms (`data/membrane_checks_summary.csv`)

| Check | FE | Reference | Difference |
|---|---|---|---|
| In-plane stress, 5850-8750 V (13 voltages) | 33 204 Pa at 5850 V | closed form | $\le0.11\%$ (7350 V excepted) |
| Lowest frequency, 0 V | 41.858 Hz | 41.858 Hz | $<0.01\%$ |
| Lowest frequency, 2 kV | 40.254 Hz | 40.254 Hz | $<0.01\%$ |
| Frequency law $f_0\sqrt{1-(V/V_c)^2}$, 0-7 kV | 11.79 Hz at 7 kV | closed form | $\le0.015\%$ |
| Linear response to $p_0=0.1$ Pa | $-0.2885$ mm at 7350 V | $w_0/(1-(V/V_c)^2)$ | $\le0.06\%$ (7350 V: $-0.79\%$) |
| Effective critical voltage | 7295.1-7295.5 V | 7295.6 V | $\le0.006\%$ |

The 7350 V point lies 54 V above $V_c$, where $w$ is singular: its $0.79\%$ difference is a $0.006\%$ shift in $V_c$ magnified by the singularity. Frequencies are undamped (modulus of the complex eigenfrequency); the real part differs by up to $0.45\%$ at 7 kV because of damping. The frequency law holds for the lowest mode only, since the pure tension scaling is valid for membrane-dominated modes.

## 3. Unstable spectrum of the flat state (shell, $V_0/V_c=1.20$)

Linearising the stretched flat state at $V_0=8755$ V, the 30 eigenvalues nearest 1 Hz were computed; 22 real unstable pairs came out. `data/shell_unstable_growth_rates_V8755.csv` and `figures/fig3_unstable_ladder.png`: the slowest mode agrees with membrane theory ($174.5$ vs $174.5\ \mathrm{s^{-1}}$), the pairs are nearly equally spaced ($226$-$227\ \mathrm{s^{-1}}$ per mode against $228$ in theory), and the departure from theory grows as $k^2$ (as bending predicts): $-4.25\%$ at $n=22$. That fits an effective bending stiffness $D_\text{eff}=(0.91$-$1.03)\times10^{-6}$ N m over $n=10$-$22$.

## 4. Pressure-deflection response

### Volume (deflection) control

Pressure control cannot pass a fold. Instead the transverse pressure $p_u$ is made an unknown (a **Global ODE**) and the deflection is prescribed:

$$f=(1-sw)\,\Big(\frac{\int w\,dA}{\int dA}-dh\Big)\cdot10^{3}\ \mathrm{[1/m]}+sw\,(p_u-p_c),$$

with the Face Load set to $p_u$. (The factor $10^3$ was added to condition the equation; the solution does not depend on it.) With $sw=0$ the mean deflection equals $dh$; with $sw=1$ the pressure equals $p_c$ (used for the pressure-control continuation and the time-domain release). The mean is an integration operator over all boundaries in the spatial frame with *Compute integral in revolved geometry* on. For centre-deflection control a second operator on the axis point is used, with the revolved option **off** (otherwise the factor $2\pi r$ makes it zero).

### Shell, $V_0/V_c=1.20$ (`data/shell_pressure_vs_mean_deflection_V8755.csv`)

The pressure rises to about $180$ Pa near $13$-$14$ mm and then falls. A cubic fit to the points to 10 mm gives $a=20.2$ Pa/mm and $b=0.0378$ Pa/mm³. The shell continuation could not be pushed beyond about $16$-$18$ mm.

### Membrane, $V_0/V_c=1.234$ (`figures/fig1_pressure_vs_deflection.png`)

| $\lvert\bar w\rvert$ (mm) | 1 | 2 | 4 | 6 | 7.7 | 8 | 10 | 11 | 11.7 |
|---|---|---|---|---|---|---|---|---|---|
| $p$ (Pa), frozen thickness | 23.8 | 46.8 | 87.1 | 114.4 | **122.8** | 122.6 | 105.2 | 85.0 | 64.6 |

`python rom/fit_cubic.py` gives, for $|\bar w|\le10$ mm, $a=23.91$ Pa/mm, $b=0.134$ Pa/mm³, an RMS residual of $0.075$ Pa ($0.06\%$ of the fold), a cubic fold at $7.7$ mm and a cubic zero at $13.35$ mm. The FE linear slope is $23.968$ Pa/mm against $25.05$ in the closed form ($-4.3\%$ in $a$; $4.5\%$ as closed form over FE). The deformed-thickness run stops at $2.2$ mm and $46.8$ Pa.

Below the fold the branch continues to a mean deflection of $11.754$ mm ($60.3$ Pa). A pressure-controlled continuation from $71.5$ Pa then turned back in mean deflection at about $62$-$63$ Pa ($11.72$ mm at $63.5$ Pa, $11.68$ mm at $61.1$ Pa), which is why mean-deflection control stalls there. Controlling the centre deflection instead goes on to $25$ mm: the mean saturates at $11.73$ mm and the pressure is still $62.5$ Pa (`figures/fig2_centre_control_saturation.png`).

### Closed-form ROM against FE

| Quantity | Closed form | FE | FE / closed form |
|---|---|---|---|
| Membrane slope $a$ (Pa/mm) | 25.05 | 23.97 | 0.957 |
| Membrane $b$, frozen thickness (Pa/mm³) | 0.0107 (with $K_{3E}$); 0.0070 (without) | 0.134 | 12.5; 19 |
| Membrane $b$, deformed thickness | 0.0107; 0.0070 | 0.58 | 54; 83 |
| Shell slope $a$ at $1.20V_c$ | 21.13 | 20.22 | 0.957 |
| Shell $b$ at $1.20V_c$ | 0.0105; 0.0070 | 0.0378 | 3.6; 5.4 |
| Fold, membrane (mm; Pa) | 27.9; 467 or 34.6; 577 | 7.7; 122.8 | |

$K_1$ is reproduced, $K_3$ is not. The deformed-thickness $b$ is fitted to points up to 2.2 mm only. Which closed form corresponds to the frozen-thickness FE model was not established, so both are listed.

## 5. Practical notes

Things that cost time and may save yours:

* Gauss-point variables such as `mbrn.misesGp` cannot be evaluated at a point. Use a **Line Average** over all boundaries. The line average of an axisymmetric model is area-weighted.
* Eigenmode datasets give mode-shape stresses, not the base-state stress. Use the stationary study's dataset.
* A second prescribed displacement that holds the membrane flat ($z=0$ everywhere) must be switched off for deflection studies, or every deflection is zero.
* An empty but ticked *Auxiliary sweep* gives "Parameter names must not be empty". Lists in a *Specified combinations* sweep must have equal length, otherwise use *All combinations*.
* Pressure sweeps stall at folds; a deflection (volume) control through a Global ODE passes them.
* A mean deflection defined on the deformed area differs from the reference-area mean by a few per cent (3.4% at 2.2 mm for the deformed-thickness run).
* Time-domain integration of the unstable flat state blows up through short-wave noise (`docs/limitations.md`).
