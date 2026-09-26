"""Forced values of phi_W(P) (min/max over the pair-level feasible set) for a given AD configuration."""
import sys, json, numpy as np
from pair_level import *
from run_all import dd_matching
from scipy.optimize import milp, LinearConstraint, Bounds

def forced(In, dd, ddc, Ds, show=True):
    ddm = dd_matching(ddc)
    A, lo, hi, lb, ub, nvar, v, nsidx, orbits = build(In, dd, ddm, Ds)
    cons = [LinearConstraint(A, lo, hi)]
    out = {}
    for o in orbits:
        row = {}
        for P in PAIRS:
            j = v(o, P); c = np.zeros(nvar); c[j] = 1
            r1 = milp(c=c, constraints=cons, integrality=np.ones(nvar), bounds=Bounds(lb, ub))
            r2 = milp(c=-c, constraints=cons, integrality=np.ones(nvar), bounds=Bounds(lb, ub))
            row[P] = (round(r1.fun), round(-r2.fun))
        out[o] = row
    ns = {}
    for p, j in nsidx.items():
        c = np.zeros(nvar); c[j] = 1
        r1 = milp(c=c, constraints=cons, integrality=np.ones(nvar), bounds=Bounds(lb, ub))
        r2 = milp(c=-c, constraints=cons, integrality=np.ones(nvar), bounds=Bounds(lb, ub))
        ns[p] = (round(r1.fun), round(-r2.fun))
    if show:
        print('n_s ranges:', ns)
        for o in orbits:
            s = ' '.join('%d%d:%s' % (P[0], P[1], ('%d' % a) if a == b else '%d-%d' % (a, b)) for P, (a, b) in out[o].items() if b > 0)
            print('%d%d/%d  D=%s  %s' % (o[0][0], o[0][1], o[1], '', s))
    return out, ns

if __name__ == '__main__':
    reps = json.load(open('ad_reps_upto5.json'))
    rep = reps[int(sys.argv[1])]
    dd = [tuple(p) for p in rep['dd']]
    head = {pr(t, h): h for t, h in rep['arcs']}
    In = {c: [x for x in range(N) if x != c and head.get(pr(c, x)) == c] for c in range(N)}
    print('In', In, 'dd', dd)
    ddc = eval(sys.argv[2]); Ds = {k: [tuple(In[k])] for k in range(N) if len(In[k]) in (2, 3)}
    Ds.update(eval(sys.argv[3]))
    forced(In, dd, ddc, Ds)
