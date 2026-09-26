"""LP relaxation + exact Farkas certificate for the pair-level models.
A model is (A, lo, hi, lb, ub): find x in the box [lb, ub] with lo <= A x <= hi.  Certificate y (rational, one per row):
  U(y) = sum_r (y_r hi_r if y_r > 0 else y_r lo_r)  >=  y^T A x  >=  L(y) = sum_j min(c_j lb_j, c_j ub_j),  c = A^T y;
so L(y) > U(y) proves infeasibility of the LP relaxation (hence of the integer problem).  y is obtained from the duals of the
phase-1 LP  min sum(s+ + s-)  s.t.  lo <= A x + s+ - s- <= hi, then rationalised and CHECKED in exact arithmetic."""
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog

def lp_infeasible_certificate(A, lo, hi, lb, ub, denom=10**6):
    m, n = A.shape
    # variables: x (n), sp (m), sm (m);  constraints split into A x + sp - sm <= hi and -(A x + sp - sm) <= -lo
    I = np.eye(m)
    Aub = np.vstack([np.hstack([A, I, -I]), np.hstack([-A, -I, I])])
    bub = np.concatenate([hi, -lo])
    c = np.concatenate([np.zeros(n), np.ones(2 * m)])
    bounds = [(lb[j], ub[j]) for j in range(n)] + [(0, None)] * (2 * m)
    res = linprog(c, A_ub=Aub, b_ub=bub, bounds=bounds, method='highs')
    if res.status != 0 or res.fun < 1e-7:
        return None                     # LP relaxation feasible (or solver trouble)
    dual = res.ineqlin.marginals      # <= 0 for <= rows (scipy convention)
    yhi, ylo = -dual[:m], -dual[m:]   # multipliers >= 0 on 'A x <= hi' and on '-A x <= -lo'
    y = [Fraction(round(a * denom), denom) - Fraction(round(b * denom), denom) for a, b in zip(yhi, ylo)]
    return y if check_certificate(A, lo, hi, lb, ub, y) else 'rounding failed'

def check_certificate(A, lo, hi, lb, ub, y):
    m, n = A.shape
    U = sum((yr * Fraction(int(hi[r])) if yr > 0 else yr * Fraction(int(lo[r]))) for r, yr in enumerate(y) if yr != 0)
    cvec = [Fraction(0)] * n
    for r, yr in enumerate(y):
        if yr == 0: continue
        row = A[r]
        for j in np.nonzero(row)[0]:
            cvec[j] += yr * int(row[j])
    L = sum(min(cj * Fraction(int(lb[j])), cj * Fraction(int(ub[j]))) for j, cj in enumerate(cvec) if cj != 0)
    return L > U

def lp_point(A, lo, hi, lb, ub):
    res = linprog(np.zeros(A.shape[1]), A_ub=np.vstack([A, -A]), b_ub=np.concatenate([hi, -lo]),
                  bounds=list(zip(lb, ub)), method='highs')
    return res.x if res.status == 0 else None

def branch_certify(A, lo, hi, lb, ub, binaries, depth=0, maxdepth=40, stats=None):
    """Branch on binary variables until every node's LP relaxation is infeasible with an exact certificate.
    Returns a tree {'fix': ..., 'cert': y} / {'var': j, 'zero': subtree, 'one': subtree}, or ('FEASIBLE', x) / None on failure."""
    if stats is not None: stats['nodes'] = stats.get('nodes', 0) + 1
    y = lp_infeasible_certificate(A, lo, hi, lb, ub)
    if isinstance(y, list):
        return {'cert': [str(v) for v in y if True]}
    x = lp_point(A, lo, hi, lb, ub)
    if x is None:
        return None
    frac = [(abs(x[j] - 0.5), j) for j in binaries if lb[j] < ub[j] and 1e-6 < x[j] < 1 - 1e-6]
    if not frac:
        # try integer rounding check
        return ('FEASIBLE?', None)
    if depth >= maxdepth:
        return None
    _, j = min(frac)
    out = {'var': int(j)}
    for val, name in ((0, 'zero'), (1, 'one')):
        lb2, ub2 = lb.copy(), ub.copy(); lb2[j] = ub2[j] = val
        sub = branch_certify(A, lo, hi, lb2, ub2, binaries, depth + 1, maxdepth, stats)
        if sub is None or isinstance(sub, tuple):
            return sub
        out[name] = sub
    return out

def check_tree(A, lo, hi, lb, ub, tree):
    if 'cert' in tree:
        return check_certificate(A, lo, hi, lb, ub, [Fraction(v) for v in tree['cert']])
    j = tree['var']
    for val, name in ((0, 'zero'), (1, 'one')):
        lb2, ub2 = lb.copy(), ub.copy(); lb2[j] = ub2[j] = val
        if not check_tree(A, lo, hi, lb2, ub2, tree[name]): return False
    return True
