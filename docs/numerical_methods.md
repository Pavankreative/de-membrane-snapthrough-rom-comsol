# Numerical methods and results of the reduced-order model

This note records what was computed numerically on the reduced model of `docs/theory.md`, with the numbers as obtained. The scripts that produced the ROM runs (manifold computation, basin maps, Lyapunov estimates) are **not yet in this repository**; they are listed under "To be added" at the end. Everything in `rom/` is reproducible from the data in `data/`.

## 1. Melnikov quadrature

The closed-form integrals of the Melnikov function were checked against adaptive quadrature for $\omega\in[0.3,3]$, with relative error below $10^{-8}$.

## 2. Stable and unstable manifolds

Direct numerical computation of the stable and unstable manifolds of the saddle reproduces both loop thresholds ($A_c^R=0.204$, $A_c^L=0.334$) to $0.1\%$, and the splitting amplitude to within $1.7\%$ up to twice the threshold.

## 3. Single-trajectory crossing tests

Starting inside the left well, $(x,x',\rho,\rho')=(-0.95,0,10^{-4},0)$, over 120 forcing periods, with a forcing scale $A$ multiplying every drive term:

| Case | Crossing threshold $A$ | $A/A_c^L$ | $A/A_c^R$ |
|---|---|---|---|
| Planar control, $\rho\equiv0$ | $\approx1.49$ | $\approx4.5$ | $\approx7.3$ |
| Coupled, $\gamma=3.0$, $K_{nn}=8.0$ (valid landscape) | $\approx0.39$ | $\approx1.17$ | $\approx1.9$ |

In the coupled case $\rho$ reaches about $0.75$ at crossing, which proceeds through the azimuthal pass (barrier $0.049$ against $0.25$ for the planar barrier). No crossing was seen for $A\le0.3$, below the left-loop onset. That tangle onset bounds crossing from below for this parameter set is observed, not proven.

**Seed dependence.** The first crossing in the coupled case moves with the seed $\rho(0)$: $A\approx0.361,\,0.383,\,0.405,\,0.463$ for $\rho(0)=10^{-2},\,10^{-4},\,10^{-8},\,10^{-12}$ (that is, $1.08$ to $1.39\,A_c^L$). The plane $\rho=0$ is transversally unstable, any seed lets the azimuthal mode grow, and the coupled threshold converges only logarithmically to the planar one. In a real membrane the mode is always seeded by imperfections.

## 4. Second-order shift $M_2$

The normal equation $\rho''=(\sigma_n^2-\nu(\nu+1)\,\mathrm{sech}^2\tau)\rho$ with $\nu(\nu+1)=2\gamma$ was solved numerically (decaying solutions integrated from $\pm14$, Wronskian constant to $10^{-10}$). Per unit $f_n^2$, $\max|M_2|=0.56$. The smallest forcing with a sign change of $M$ moves from $1.00A_c^R$ to $0.95A_c^R$ for $f_n=1$. It does not explain the drop from $\approx7.3A_c^R$ to $\approx1.9A_c^R$. $M_2$ for the left loop was not evaluated.

## 5. Basin erosion

Maps on a $200\times200$ grid of in-well starts with energy $H<-0.01$, seeded at $\rho(0)=10^{-4}$, and a horizon of 20 forcing periods. Labels agree to at least $99\%$ between step sizes $T/100$ and $T/200$. At 60 periods the agreement at $A=0.3$ falls to about $83\%$, so longer horizons are not converged and were not used.

| $A$ (coupled) | 0.10 | 0.2044 | 0.30 | 0.39 | 0.50 |
|---|---|---|---|---|---|
| Escaped fraction (%) | 4.9 | 5.9 | 16.0 | 34.7 | 58.3 |
| Uncertainty exponent $\alpha$ | 0.75 | 0.41 | 0.09 | 0.055 | 0.059 |

Planar case: escaped fraction $2.8,\,65.3,\,96.5\%$ and $\alpha=0.28,\,0.29,\,0.57$ at $A=0.39,\,1.0,\,1.49$. ($\alpha=1$ is a smooth boundary.)

## 6. Lyapunov exponents

The largest exponent is negative for the planar subsystem (tangling without sustained chaos) and about $+0.2$ for the coupled system at $A=0.39$ and $0.50$. These are finite-time estimates from one parameter set.

## 7. Single-mode ROM release test (prediction only)

`rom/rom_release_prediction.py` integrates the mean-deflection ROM at zero load,

$$M\ddot y+C\dot y=a\,y-b\,y^3,\qquad y=|\bar w| ,$$

with the FE slope $a=23.968$ Pa/mm and the FE cubic coefficient $b=0.134$ Pa/mm³ (frozen-thickness membrane, $V_0=9000$ V), the effective mass $M=0.664$ kg/m² obtained from the FE first-mode growth rate at 8755 V, and the mass-damping rate $14.4\ \mathrm{s^{-1}}$ of the FE model. From a start at $y_0=0.1$ mm it predicts $2.0$ mm at $20$ ms, $4.9$ mm at $25$ ms, $8$ mm at $27.8$ ms, and an overshoot to $18.4$ mm at $36$ ms.

**This has not been compared with FE.** A time-domain FE run failed within $0.7$ ms through growth of short-wave noise (`docs/limitations.md`). The prediction is kept for transparency and for a later comparison in a model where the short-wave band is suppressed.

## To be added

* Planar and coupled ROM integrators, manifold computation, Melnikov check
* Basin-map and Lyapunov scripts
* A fit-and-compare script for the full calibration protocol (target matrix) once FE data exist
