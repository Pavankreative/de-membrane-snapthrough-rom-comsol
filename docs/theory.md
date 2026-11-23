# Theory: reduced-order model of a pre-stretched dielectric-elastomer membrane

This note collects the analytical part of the project. Everything here is derived from the free energy of a Gent membrane with a Maxwell stress. Section 6 states where the single-mode reduction stops being valid. The numbers quoted are the reference values used in the FE comparison (`docs/comsol_validation.md`).

Notation: $R$ reference radius, $H_0$ reference thickness, $\lambda_p$ equibiaxial pre-stretch, $R_s=\lambda_pR$ stretched radius, $h=H_0/\lambda_p^2$ stretched thickness, $\mu$ shear modulus, $J_m$ Gent locking parameter, $\varepsilon=\varepsilon_0\varepsilon_r$ permittivity, $V$ voltage.

## 1. Constitutive model

Helmholtz free energy per unit reference volume, with $I_1=\operatorname{tr}\mathbf C$:

$$\Psi_{\text{Gent}}=-\frac{\mu J_m}{2}\ln\!\Big(1-\frac{I_1-3}{J_m}\Big),\qquad
\Psi_{\text{Maxwell}}=-\frac{\varepsilon}{2}\,J\,C^{-1}_{IJ}E_IE_J .$$

The Gent form bounds the admissible stretch through the singularity $(1-x_I)^{-2}$, $x_I=(I_1-3)/J_m$. The Maxwell stress of a field $V/h$ acts through the thickness and reduces the in-plane tension.

## 2. Pre-stretched flat state and the critical voltage

For an equibiaxially stretched incompressible Gent membrane, with $s=(I_0-3)/J_m$ and $I_0=2\lambda_p^2+\lambda_p^{-4}$:

$$\sigma_0=\frac{\mu\,(\lambda_p^2-\lambda_p^{-4})}{1-s},\qquad N=\sigma_0h .$$

The Maxwell stress reduces the net in-plane tension to

$$N_\text{eff}=N-\frac{\varepsilon V^2}{h}.$$

The reference voltage $V_c$ is where $N_\text{eff}$ vanishes:

$$V_c=\sqrt{\frac{Nh}{\varepsilon}} .$$

With the constants used here ($\mu=27.46$ kPa, $J_m=50$, $\lambda_p=1.8$, $\varepsilon_r=4.7$, $R=50$ mm, $H_0=0.5$ mm) this gives $\sigma_0=93\,005$ Pa, $N=14.353$ N/m and $V_c=7295.6$ V (reproduced by `rom/closed_form_coefficients.py`). For $V>V_c$ the net tension is negative and the flat state is unstable.

**Sensitivity.** $V_c$ varies by about $\pm7\%$ (7160 to 8240 V) across published ranges of $\mu$ and $J_m$ for VHB 4910, so voltages are quoted as ratios $V_0/V_c$.

## 3. Reduction to the first axisymmetric mode

Take the mode-1 shape $\Phi_1=J_0(j_{0,1}r/R)$ with amplitude $x=q_1$, where $j_{0,1}=2.4048$ is the first zero of $J_0$. Expanding the Gent energy in the slope $e=|\nabla w|^2$:

$$I_1-3=(I_0-3)+\beta_1e+\beta_2e^2+\cdots,\quad \beta_1=1-\lambda_p^{-6},\ \ \beta_2=\lambda_p^{-8},$$

$$\Psi_\text{Gent}\approx\Psi_0+\frac{\mu}{2}\frac{\beta_1}{1-s}\,e+\frac{\mu}{2}\,G(s)\,e^2,\qquad
G(s)=\frac{\beta_2}{1-s}+\frac{\beta_1^2}{2J_m(1-s)^2}.$$

The reduced potential in the first mode is

$$V(x)=-\tfrac12K_1x^2+\tfrac14K_3x^4,$$

$$K_1=K_E-K_G=K_G\Big[\Big(\frac{V_0}{V_c}\Big)^2-1\Big],\qquad
K_G=\frac{\mu H_0\beta_1J_2}{1-s},\qquad K_E=\frac{\varepsilon V_0^2\lambda_p^2J_2}{H_0},$$

$$K_3=2\mu H_0\,G(s)\,J_4+K_{3E},\qquad K_{3E}=\frac{\varepsilon V_0^2\,(A_0J_4-J_2^2)}{2H_0A_0}\ \ge0,$$

with $J_2=\int|\nabla\Phi_1|^2dA_0=\pi j_{0,1}^2J_1(j_{0,1})^2=4.897$ and $J_4=\int|\nabla\Phi_1|^4dA_0$, $A_0=\pi R^2$. For $V_0>V_c$ the coefficient $K_1>0$ and the flat state is a maximum of $V$ with two wells at $x=\pm\sqrt{K_1/K_3}$. The modal mass is $m_\text{eff}=\rho H_0\pi R^2J_1(j_{0,1})^2$ and $\omega_0^2=K_1/m_\text{eff}$.

### Mapping to a pressure law

For a uniform load $p$ per reference area the equilibrium of the mode gives, in terms of the mean deflection $\bar w=\bar a\,x$ with $\bar a=2J_1(j_{0,1})/j_{0,1}=0.4317$,

$$p=a\,|\bar w|-b\,|\bar w|^3,\qquad a=\frac{K_1}{A_0\bar a^2},\quad b=\frac{K_3}{A_0\bar a^4}.$$

The pressure here is the load that holds the deflection. Because the flat state has negative stiffness, it acts opposite to the deflection. The curve has a fold at $|\bar w|_f=\sqrt{a/(3b)}$ with $p_\max=\tfrac23a|\bar w|_f$ and a zero-load point at $|\bar w|=\sqrt{a/b}$ (the well of the cubic model).

Closed-form values at $V_0=9000$ V ($V_0/V_c=1.234$), from `rom/closed_form_coefficients.py`:

| Quantity | Value |
|---|---|
| $K_1$ | 36.67 |
| $K_3$ (Gent part) | 1906 |
| $K_{3E}$ | 1012 |
| $a$ | 25.05 Pa/mm |
| $b$ with $K_{3E}$ | 0.0107 Pa/mm³ |
| $b$ without $K_{3E}$ | 0.0070 Pa/mm³ |

The FE comparison in `docs/comsol_validation.md` shows that $a$ agrees to $4.5\%$ and that the closed-form $b$ is far too small. Use the FE-calibrated coefficient, not the closed form.

## 4. Normalised planar system and critical points

Rescale $x=x_sX$ with $x_s^2=K_1/K_3$ and time by $\omega_0^{-1}$. Retaining one azimuthal mode $\Psi_n\cos n\phi$ with amplitude $\rho$ (the O(2) invariant) and a viscous damping $\delta$, the autonomous part is the gradient of

$$V(x,\rho)=-\tfrac12x^2+\tfrac14x^4-\tfrac12\sigma_n^2\rho^2+\tfrac14K_{nn}\rho^4+\tfrac12\gamma\,x^2\rho^2,\qquad \sigma_n=\frac{j_{n,1}}{j_{0,1}}.$$

For the planar subsystem ($\rho=0$) the wells are $x=\pm1$ with $V=-\tfrac14$, joined to the saddle at $x=0$ by two homoclinic loops $x_0=\pm\sqrt2\,\mathrm{sech}\,\tau$. When all exist, the system has nine critical points: the saddle $(0,0)$, the wells $(\pm1,0)$, the azimuthal points $(0,\pm\rho^*)$ with $\rho^{*2}=\sigma_n^2/K_{nn}$, and four mixed points. The wells are the intended global minima if

$$\gamma>\sigma_n^2,\qquad K_{nn}>\sigma_n^4 .$$

The lowest route between the wells is the azimuthal pass $(0,\pm\rho^*)$ when $K_{nn}>\gamma\sigma_n^2$, with barrier $\Delta V_\rho=\tfrac14(1-\sigma_n^4/K_{nn})<\tfrac14$.

## 5. Melnikov criterion for homoclinic tangling

Forced planar subsystem with quadratic drag $\delta_2$, an aerodynamic drive $f\cos\omega\tau$ and electrical modulation $(g_1-g_3x^2)x\cos\omega\tau$:

$$x''=x-x^3+\epsilon\Big[-\delta x'-\delta_2|x'|x'+f\cos\omega\tau+(g_1-g_3x^2)\,x\cos\omega\tau\Big].$$

On the loop with sign $s=\pm1$, the Melnikov function is

$$\frac{M_s(\tau_0)}{\epsilon}=-\Big[\tfrac43\delta+\tfrac{8\sqrt2}{15}\delta_2\Big]+\pi\omega\Big[s\sqrt2\,f\,\mathrm{sech}\tfrac{\pi\omega}{2}+\frac{\omega\,g_\text{eff}}{\sinh(\pi\omega/2)}\Big]\sin\omega\tau_0,
\qquad g_\text{eff}=g_1-\frac{\omega^2+4}{6}g_3 .$$

It has a simple zero, and the stable and unstable manifolds intersect transversally, if

$$\Big|\pi\omega\Big[s\sqrt2\,f\,\mathrm{sech}\tfrac{\pi\omega}{2}+\frac{\omega\,g_\text{eff}}{\sinh(\pi\omega/2)}\Big]\Big|>\tfrac43\delta+\tfrac{8\sqrt2}{15}\delta_2 .$$

The tabulated integrals used are $\int x_0'^2d\tau=\tfrac43$, $\int|x_0'|^3d\tau=\tfrac{8\sqrt2}{15}$, $\int\mathrm{sech}\,\tau\cos\omega\tau\,d\tau=\pi\,\mathrm{sech}\tfrac{\pi\omega}{2}$ and the two higher $\mathrm{sech}$ integrals, which agree with adaptive quadrature to better than $10^{-8}$ for $\omega\in[0.3,3]$.

The electrical amplitude threshold, with $\kappa_E(\omega)=K_E/K_1-\tfrac{\omega^2+4}{6}K_{3E}/K_3$:

$$\Big(\frac{V_1}{V_0}\Big)_c=\frac{\Big[\tfrac43\delta+\tfrac{8\sqrt2}{15}\delta_2-s\sqrt2\,\pi\omega f\,\mathrm{sech}\tfrac{\pi\omega}{2}\Big]\sinh\tfrac{\pi\omega}{2}}{2\pi\omega^2\,\kappa_E(\omega)} .$$

At the reference parameters the two loops have different thresholds, because $f$ and $g_\text{eff}$ add on one loop and subtract on the other: $A_c^R=0.204$ ($s=+1$) and $A_c^L=0.334$ ($s=-1$).

### What the criterion is, and is not

* It is the sharp analytical boundary for the onset of a transverse homoclinic tangle (Smale horseshoe, fractal basin erosion), to first order in $\epsilon$. Corrections are $O(\epsilon^2)$.
* It is **not** a deterministic well-to-well escape boundary. Crossing needs more forcing, at the point where resonance layers overlap. In single-trajectory tests from the left well, crossing needed about $1.1$ to $1.4$ times $A_c^L$ in the coupled, seed-dependent case and about $4.5$ times in the planar case.
* It is meaningful only when the landscape conditions of section 4 hold.

A second-order shift $M_2$ from the azimuthal mode was evaluated numerically. It does not explain the lower coupled threshold, which is a property of the dynamics off the plane $\rho=0$ (energy carried through the azimuthal pass).

## 6. Where the one-mode reduction stops

For $V_0>V_c$ the net tension is negative, so every axisymmetric mode of the flat state is unstable, not only the first. For a membrane without bending the growth rate of mode $n$ is

$$s_n=c\,k_n,\quad k_n=\frac{j_{0,n}}{R_s},\quad c=\sqrt{\frac{|N_\text{eff}|}{\rho h}}$$

($c=6.53$ m/s and $s_1=174.5\ \mathrm{s^{-1}}$ at $V_0/V_c=1.20$). The growth rate rises linearly with $n$, so the first mode is the slowest unstable direction, not the fastest. Bending stiffness $D$ cuts the band off at $k_c=\sqrt{|N_\text{eff}|/D}$ (FE estimates in `docs/comsol_validation.md`). The ROM therefore describes the first-mode subsystem only. It is the dominant direction if the short-wave band is suppressed, for example by a stiff backing, or if the buckled states are stable against it. Neither was shown here.

## 7. References

Primary sources still to be added with checked bibliographic details (Gent 1996 for the strain-energy function, Guckenheimer and Holmes for Melnikov theory, and the VHB 4910 data sources for $\mu$, $J_m$, $\varepsilon_r$). The value $\varepsilon_r=4.7$ follows Keplinger et al. (2012); Molberg et al. (2009) report 4.65. $J_m=50$ is assumed and has no source.
