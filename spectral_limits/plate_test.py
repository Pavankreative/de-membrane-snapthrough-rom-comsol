"""Shell-bending check at lambda = 1 (unstretched plate): centre deflection under uniform pressure, clamped rim.
Inputs are the COMSOL result of DE_plate_test.mph (Study 7); D is compared with E h^3 / [12 (1 - nu^2)], E = 3 mu."""
MU, H0, NU, RMEM = 27.46e3, 0.5e-3, 0.5, 0.05
Q = 5e-4                      # Pa
W_CENTRE = 42.528e-6          # m  (COMSOL, 3 mu E, clamped)
D_FE = Q * RMEM**4 / (64 * W_CENTRE)
D_TH = 3 * MU * H0**3 / (12 * (1 - NU**2))
print(f'D from COMSOL deflection = {D_FE:.4e} N m ; thin-plate E h^3/[12(1-nu^2)] = {D_TH:.4e} N m ; ratio {D_FE/D_TH:.4f}')
print(f'D_eff(6.73e-7) / unstretched plate = {6.73e-7/D_TH:.3f}')
