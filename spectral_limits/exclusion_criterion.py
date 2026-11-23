"""Exclusion criterion of Sec. 5:  g3/g1 = (K3E/K3) * Delta/(1+Delta) <= Delta/(1+Delta),  Delta = (V0/Vc)^2 - 1.
Reference coefficients g1 = 0.50, g3 = 0.08 need  Delta >= r/(1-r),  r = g3/g1.
Compared with the first competing (1,1) crossings from the fitted-rigidity plate model (D_eff = 6.73e-7)."""
from plate_dispersion import KR, crossing
G1, G3, DEFF = 0.50, 0.08, 6.73e-7
r = G3 / G1; delta_min = r / (1 - r)
print(f'g3/g1 = {r:.2f}  ->  Delta >= {delta_min:.5f}')
for rim in ('pinned', 'clamped'):
    xc = crossing(DEFF, KR[rim][(1, 1)])
    print(f'{rim:8s} (1,1) crossing {xc:.4e}   Delta_min / x_c = {delta_min / xc:,.0f}')
