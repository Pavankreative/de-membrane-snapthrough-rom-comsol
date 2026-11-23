import numpy as np
from scipy.integrate import solve_ivp

def validity(gamma, K_nn, sigma_n):
    """Dimensionless (K1=K3=1, x_s=1) global well-validity checks of Sec. 3.2."""
    s2 = sigma_n**2
    c1 = gamma > s2                                   # local well stability: V_rr(+-1,0) = gamma - sigma_n^2 > 0
    c2 = K_nn > s2**2                                 # wells are global minima: -sigma_n^4/(4 K_nn) > -1/4
    c3 = True                                         # no mixed minimum: det H = 4 x^2 rho^2 (K_nn - gamma^2) with rho_m^2=(s2-gamma)/(K_nn-gamma^2)
    if abs(K_nn - gamma**2) > 1e-12:
        rho2 = (s2 - gamma) / (K_nn - gamma**2); x2 = 1.0 - gamma * rho2
        if rho2 > 0 and x2 > 0 and K_nn > gamma**2:
            c3 = False
    pass_is_saddle = K_nn > gamma * s2                # (0,rho*) is a saddle (not a minimum) iff K_nn > gamma sigma_n^2
    dV_rho = 0.25 * (1.0 - s2**2 / K_nn)              # lowest barrier through the azimuthal pass (valid if pass_is_saddle)
    return dict(c1=c1, c2=c2, c3=c3, valid=c1 and c2 and c3, pass_is_saddle=pass_is_saddle, dV_rho=dV_rho)

def melnikov_margin(omega, f, g1, g3, delta, delta_2):
    s2 = np.sqrt(2.0)
    g_eff = g1 - (omega**2 + 4.0) / 6.0 * g3
    lhs = np.pi * omega * (s2 * f / np.cosh(np.pi * omega / 2.0) + omega * g_eff / np.sinh(np.pi * omega / 2.0))
    rhs = 4.0 / 3.0 * delta + 8.0 * s2 / 15.0 * delta_2
    return abs(lhs) - rhs, rhs / abs(lhs)            # margin Delta and critical forcing scale A_c

def vratio_c(omega, f, delta, delta_2, KE_over_K1, K3E_over_K3):
    s2 = np.sqrt(2.0)
    kappa = KE_over_K1 - (omega**2 + 4.0) / 6.0 * K3E_over_K3
    num = (4.0/3.0*delta + 8.0*s2/15.0*delta_2 - s2*np.pi*omega*f/np.cosh(np.pi*omega/2.0)) * np.sinh(np.pi*omega/2.0)
    return num / (2.0*np.pi*omega**2*kappa)

def stabilized_hyperelastic_system(t, z, p):
    x, x_dot, rho, rho_dot = z
    forcing = (p['f_aero'] + (p['g1'] - p['g3']*x**2)*x) * np.cos(p['omega']*t)
    dx_dot = (x - x**3 - p['delta']*x_dot - p['delta_2']*x_dot*abs(x_dot) + forcing - p['gamma']*x*rho**2)
    drho_dot = (p['sigma_n']**2*rho - p['K_nn']*rho**3 - p['gamma']*x**2*rho
                - p['delta']*rho_dot - p['delta_2']*rho_dot*abs(rho_dot))
    return [x_dot, dx_dot, rho_dot, drho_dot]

def envelope_point(p, tol_periods=120):
    v = validity(p['gamma'], p['K_nn'], p['sigma_n'])
    if not v['valid']:
        return dict(status='REJECT_INVALID_LANDSCAPE', **v)
    margin, A_c = melnikov_margin(p['omega'], p['f_aero'], p['g1'], p['g3'], p['delta'], p['delta_2'])
    return dict(status='TANGLE' if margin > 0 else 'NO_TANGLE', margin=margin, A_c=A_c, **v)

if __name__ == '__main__':
    base = dict(delta=0.08, delta_2=0.03, f_aero=0.12, g1=0.50, g3=0.08, sigma_n=1.593, omega=1.72)
    for g, K in ((0.15, 1.25), (3.0, 1.25), (3.0, 8.0)):
        print(g, K, envelope_point(dict(base, gamma=g, K_nn=K)))
    print('Vratio_c(w=1)=', vratio_c(1.0, 0.01, 0.08, 0.02, 1.5, 0.3))
