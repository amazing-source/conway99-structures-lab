"""Effect of (SEC) on step 2 for (0,6) types: secant iota-pairs (common label v) may only be assigned to sigma-pairs whose
root meets block(v) (valid when there is no dis D-square, e.g. in the Sigma = empty all-one-block family)."""
import json, sys, time, itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
import iota_design as D
for fn in sys.argv[1:]:
    cfg = json.load(open(fn)); mu, tau = cfg['mu'], cfg['tau']
    K2, K3 = D.K2K3(mu, tau)
    t0 = time.time()
    st, sol, edges = D.step1(mu, tau, 1, 70)
    if sol is None:
        print(fn, 'step 1', st); continue
    far = [z for z in itertools.combinations(D.L, 2) if z[0] // 2 != z[1] // 2]
    reps, seen = [], set()
    for z in far:
        if z in seen: continue
        sz = tuple(sorted((D.PI[z[0]], D.PI[z[1]]))); seen |= {z, sz}; reps.append((z, sz))
    for use_sec in (False, True):
        var = []
        for Zi in range(30):
            z, sz = reps[Zi]
            for Pj, (a, b) in enumerate(sol):
                for (e, f) in ((edges[a], edges[b]), (edges[b], edges[a])):
                    com = set(e) & set(f)
                    if use_sec and com and not any(v // 2 in {w // 2 for w in z} for v in com):
                        continue
                    var.append((Zi, Pj, e, f))
        byZ, byP, acc = {}, {}, {}
        for k, (Zi, Pj, e, f) in enumerate(var):
            byZ.setdefault(Zi, []).append(k); byP.setdefault(Pj, []).append(k)
            z, sz = reps[Zi]
            for a in z:
                for b in e: acc.setdefault((a, b), []).append(k)
            for a in sz:
                for b in f: acc.setdefault((a, b), []).append(k)
        rows = [(byZ.get(Zi, []), 1) for Zi in range(30)] + [(byP.get(Pj, []), 1) for Pj in range(30)] + \
               [(acc.get((a, b), []), int(K2[a, b])) for a in D.L for b in D.L]
        ri = [r for r, (c, _) in enumerate(rows) for _ in c]; ci = [x for c, _ in rows for x in c]
        A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(rows), len(var))); bb = np.array([t for _, t in rows], float)
        lp = linprog(np.zeros(len(var)), A_eq=A, b_eq=bb, bounds=(0, 1), method='highs')
        s = {0: 'feasible', 2: 'INFEASIBLE'}.get(lp.status, lp.status)
        ms = ''
        if lp.status == 0:
            ip = milp(np.zeros(len(var)), constraints=[LinearConstraint(A, bb, bb)], integrality=np.ones(len(var)), bounds=Bounds(0, 1), options={'time_limit': 40})
            ms = {0: 'FEASIBLE', 2: 'INFEASIBLE'}.get(ip.status, 'unknown(%d)' % ip.status)
        print('%s  SEC=%-5s step-2 vars %5d  LP %-10s MILP %s  %.0fs' % (fn, use_sec, len(var), s, ms, time.time() - t0), flush=True)
