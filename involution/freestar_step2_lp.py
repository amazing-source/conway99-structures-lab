"""Free-star check: for a saved random type (from random_types_twostep.py), find an iota (step 1) and test step 2 by LP
relaxation and by MILP (longer)."""
import json, sys, time, itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
import iota_design as D
cfg = json.load(open(sys.argv[1])); mu, tau = cfg['mu'], cfg['tau']
K2, K3 = D.K2K3(mu, tau)
t0 = time.time()
st, sol, edges = D.step1(mu, tau, 1, 100)
print('step 1', st, '%.0fs' % (time.time() - t0), flush=True)
if sol:
    far = [z for z in itertools.combinations(D.L, 2) if z[0] // 2 != z[1] // 2]
    reps, seen = [], set()
    for z in far:
        if z in seen: continue
        sz = tuple(sorted((D.PI[z[0]], D.PI[z[1]]))); seen |= {z, sz}; reps.append((z, sz))
    var = [(Zi, Pj, e, f) for Zi in range(30) for Pj, (i, j) in enumerate(sol) for (e, f) in ((edges[i], edges[j]), (edges[j], edges[i]))]
    R = {}
    byZ, byP, acc = {}, {}, {}
    for k, (Zi, Pj, e, f) in enumerate(var):
        byZ.setdefault(Zi, []).append(k); byP.setdefault(Pj, []).append(k)
        z, sz = reps[Zi]
        for a in z:
            for b in e: acc.setdefault((a, b), []).append(k)
        for a in sz:
            for b in f: acc.setdefault((a, b), []).append(k)
    rows = [(byZ[Zi], 1) for Zi in range(30)] + [(byP[Pj], 1) for Pj in range(30)] + \
           [(acc.get((a, b), []), int(K2[a, b])) for a in D.L for b in D.L]
    ri = [r for r, (c, _) in enumerate(rows) for _ in c]; ci = [x for c, _ in rows for x in c]
    A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(rows), len(var))); bb = np.array([t for _, t in rows], float)
    lp = linprog(np.zeros(len(var)), A_eq=A, b_eq=bb, bounds=(0, 1), method='highs')
    print('step 2 LP:', {0: 'feasible', 2: 'INFEASIBLE'}.get(lp.status, lp.status), '%.0fs' % (time.time() - t0), flush=True)
    ip = milp(np.zeros(len(var)), constraints=[LinearConstraint(A, bb, bb)], integrality=np.ones(len(var)), bounds=Bounds(0, 1),
              options={'time_limit': 150})
    print('step 2 MILP:', {0: 'FEASIBLE', 2: 'INFEASIBLE'}.get(ip.status, 'unknown(%d)' % ip.status), '%.0fs' % (time.time() - t0), flush=True)
