# Numerical methods and results of the reduced-order model

The [revision 17 dynamics record](../planar_dynamics/README.md) gives the executed methods, exact tabulated outputs, 20-period first-crossing definition, and reproducibility commands for Tables 3–5 and Figures 5–8. Local finite-time crossing brackets and all sampled amplitudes for Table 2 are in `planar_dynamics/data/seed_brackets.csv` and `seed_trials.csv`.

## 1. Melnikov quadrature

The closed-form integrals of the Melnikov function were checked against adaptive quadrature for $\omega\in[0.3,3]$, with relative error below $10^{-8}$.

## 2. Second-order shift $M_2$

The normal equation $\rho''=(\sigma_n^2-\nu(\nu+1)\,\mathrm{sech}^2\tau)\rho$ with $\nu(\nu+1)=2\gamma$ was solved numerically (decaying solutions integrated from $\pm14$, Wronskian constant to $10^{-10}$). Per unit $f_n^2$, $\max|M_2|=0.56$. The smallest forcing with a sign change of $M$ moves from $1.00A_c^R$ to $0.95A_c^R$ for $f_n=1$. It does not explain the drop from $\approx7.3A_c^R$ to $\approx1.9A_c^R$. $M_2$ for the left loop was not evaluated.

## 3. Single-mode ROM release test (prediction only)

`rom/rom_release_prediction.py` integrates the mean-deflection ROM at zero load,

$$M\ddot y+C\dot y=a\,y-b\,y^3,\qquad y=|\bar w| ,$$

with the FE slope $a=23.968$ Pa/mm and the FE cubic coefficient $b=0.134$ Pa/mm³ (frozen-thickness membrane, $V_0=9000$ V), the effective mass $M=0.664$ kg/m² obtained from the FE first-mode growth rate at 8755 V, and the mass-damping rate $14.4\ \mathrm{s^{-1}}$ of the FE model. From a start at $y_0=0.1$ mm it predicts $2.0$ mm at $20$ ms, $4.9$ mm at $25$ ms, $8$ mm at $27.8$ ms, and an overshoot to $18.4$ mm at $36$ ms.

**This has not been compared with FE.** A time-domain FE run failed within $0.7$ ms through growth of short-wave noise (`docs/limitations.md`). The prediction is kept for transparency and for a later comparison in a model where the short-wave band is suppressed.
