"""Single formatting source for the displayed revision-17 table values."""
from common17 import read_rows

def table3():
    return [[r['loop'].capitalize(),f"{float(r['S']):.6f}",f"{float(r['first_order']):.6f}",
             f"{float(r['A16']):.6f} / {float(r['A32']):.6f}",f"{float(r['percent']):+.3f}%"] for r in read_rows('table3.csv')]

def table4():
    return [[r['model'].capitalize(),r['A'],f"{float(r['median']):+.4f}",
             f"{float(r['minimum']):+.4f} to {float(r['maximum']):+.4f}",f"{r['positive']}/30"] for r in read_rows('table4.csv')]

def table5():
    rows=read_rows('table5.csv')
    def row(g,k,a,d=1):
        return next(r for r in rows if (float(r['gamma']),float(r['Knn']),float(r['A']),int(r['damped']))==(g,k,a,d))
    damp=row(.15,1.25,0); cons=row(.15,1.25,0,0)
    half=row(3,1.25,.5); one=row(3,1.25,1.)
    cases=[row(3,8,a) for a in (.4,.5,1.)]
    assert all(int(row(3,8,a)['n_cross'])==0 for a in (0,.1,.2,.3))
    return [
        ['0.15, 1.25','1, 2, 3 fail',f"Unforced, damped: final (x, ρ) = ({float(damp['final_x']):.3f}, {float(damp['final_rho']):.3f}). Unforced, undamped: {cons['n_cross']} sign changes in 300 time units."],
        ['3.0, 1.25','2 fails',f"A = 0.5: {half['n_cross']} sign changes; max |ρ| = {float(half['max_abs_rho']):.3f}. A = 1.0: {one['n_cross']} sign changes; max |x| over the last 20 periods = {float(one['tail_max_abs_x']):.3f}, final ρ = {float(one['final_rho']):.3f}."],
        ['3.0, 8.0','All pass',"At A = 0, 0.1, 0.2 and 0.3: zero sign changes. At A = 0.4, 0.5 and 1.0: "+', '.join(r['n_cross'] for r in cases)+" sign changes, with max |ρ| = "+', '.join(f"{float(r['max_abs_rho']):.3f}" for r in cases)+", respectively."]
    ]
