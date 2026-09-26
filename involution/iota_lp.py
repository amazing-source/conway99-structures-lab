"""For iota's found by step 1, test the LP relaxation of step 2 (is (ii) already LP-infeasible?), and save the iota's."""
import json, sys, time, itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
import iota_design as D
from fullblock_rigid import mates_from_D

cfg = json.load(open(sys.argv[1])); n = int(sys.argv[2]); tl = float(sys.argv[3])
m = mates_from_D([tuple(c) for c in cfg['cycles']]); mu = tau = m
found = []
out = []
for k in range(n):
    st, sol, edges = D.step1(mu, tau, 1, tl, forbid=found)
    if sol is None:
        print('step1', st); break
    found.append(sol)
    # rebuild step-2 matrices and run the LP relaxation
    K2, K3 = D.K2K3(mu, tau)
    far = [z for z in itertools.combinations(D.L, 2) if z[0] // 2 != z[1] // 2]
    reps, seen = [], set()
    for z in far:
        if z in seen: continue
        sz = tuple(sorted((D.PI[z[0]], D.PI[z[1]]))); seen |= {z, sz}; reps.append((z, sz))
    var = [(Zi, Pi_, e, f) for Zi in range(30) for Pi_, (i, j) in enumerate(sol)
           for (e, f) in ((edges[i], edges[j]), (edges[j], edges[i]))]
    rows = []; lo = []
    byZ = {}; byP = {}; acc = {}
    for q, (Zi, Pi_, e, f) in enumerate(var):
        byZ.setdefault(Zi, []).append(q); byP.setdefault(Pi_, []).append(q)
        z, sz = reps[Zi]
        for a in z:
            for b in e: acc.setdefault((a, b), []).append(q)
        for a in sz:
            for b in f: acc.setdefault((a, b), []).append(q)
    R = []
    for Zi in range(30): R.append((byZ[Zi], 1))
    for Pi_ in range(30): R.append((byP[Pi_], 1))
    for a in D.L:
        for b in D.L: R.append((acc.get((a, b), []), int(K2[a, b])))
    ri = [r for r, (cols, _) in enumerate(R) for c in cols]; ci = [c for cols, _ in R for c in cols]
    A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(R), len(var))); b = np.array([v for _, v in R], float)
    lp = linprog(np.zeros(len(var)), A_eq=A, b_eq=b, bounds=(0, 1), method='highs')
    print('iota #%d: LP relaxation of step 2:' % k, {0: 'feasible', 2: 'INFEASIBLE'}.get(lp.status, lp.status), flush=True)
    out.append([[list(edges[i]), list(edges[j])] for i, j in sol])
json.dump(out, open(sys.argv[1].replace('.json', '_iotas.json'), 'w'))
