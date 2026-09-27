"""Two-step evidence on C1's relevant tier-8 classes (profile (2,4)): step 1 (iota near-design, CP-SAT), then step 2 LP and
MILP.  Usage: python tier8_twostep.py idx1 idx2 ...  (indices into largeT_tier8_rel.json)"""
import json, sys, time, itertools
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
import iota_design as D
T = json.load(open('C:/home/researcher/Desktop/c1_recherche/data/largeT_tier8_rel.json'))
t0 = time.time()
for idx in map(int, sys.argv[1:]):
    mu, tau = T[idx]['mu'], T[idx]['tau']
    K2, K3 = D.K2K3(mu, tau)
    st, sol, edges = D.step1(mu, tau, 1, 70)
    if sol is None:
        print('class %d: step 1 %s' % (idx, st), '%.0fs' % (time.time() - t0), flush=True); continue
    far = [z for z in itertools.combinations(D.L, 2) if z[0] // 2 != z[1] // 2]
    reps, seen = [], set()
    for z in far:
        if z in seen: continue
        sz = tuple(sorted((D.PI[z[0]], D.PI[z[1]]))); seen |= {z, sz}; reps.append((z, sz))
    var = [(Zi, Pj, e, f) for Zi in range(30) for Pj, (i, j) in enumerate(sol) for (e, f) in ((edges[i], edges[j]), (edges[j], edges[i]))]
    byZ, byP, acc = {}, {}, {}
    for k, (Zi, Pj, e, f) in enumerate(var):
        byZ.setdefault(Zi, []).append(k); byP.setdefault(Pj, []).append(k)
        z, sz = reps[Zi]
        for a in z:
            for b in e: acc.setdefault((a, b), []).append(k)
        for a in sz:
            for b in f: acc.setdefault((a, b), []).append(k)
    rows = [(byZ[Zi], 1) for Zi in range(30)] + [(byP[Pj], 1) for Pj in range(30)] + [(acc.get((a, b), []), int(K2[a, b])) for a in D.L for b in D.L]
    ri = [r for r, (c, _) in enumerate(rows) for _ in c]; ci = [x for c, _ in rows for x in c]
    A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(rows), len(var))); bb = np.array([t for _, t in rows], float)
    lp = linprog(np.zeros(len(var)), A_eq=A, b_eq=bb, bounds=(0, 1), method='highs')
    lps = {0: 'feasible', 2: 'INFEASIBLE'}.get(lp.status, lp.status)
    st2 = D.step2(mu, tau, sol, edges, tl=40)[0] if lp.status == 0 else 'INFEASIBLE (LP)'
    print('class %d (pivots %s): step 1 FEASIBLE; step 2 LP %s; step 2 %s' % (idx, T[idx]['pivots'], lps, st2), '%.0fs' % (time.time() - t0), flush=True)
