"""Coarsest level: only the BLOCK PROFILE n_P (labels of each block in the quadruple) matters.
  sum_{P: p in x_P} n_{P,r} = T[p,r] = sum_{b in r} (K2[p+,b] + K2[p-,b])   for all blocks p, r.
Is the LP already infeasible there?  Also print the profile statistics of each near-design."""
import json, sys, itertools
from collections import Counter
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
from fullblock_rigid import mates_from_D
import iota_design as D

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
m = mates_from_D([tuple(c) for c in cfg['cycles']]); K2, K3 = D.K2K3(m, m)
BP = list(itertools.combinations(range(6), 2))
T = np.array([[sum(K2[2 * p, b] + K2[2 * p + 1, b] for b in (2 * r, 2 * r + 1)) for r in range(6)] for p in range(6)])
print('T =\n', T)
for k, pairs in enumerate(iotas):
    Q = [set(e) | set(f) for e, f in pairs]
    prof = [tuple(sum(1 for b in q if b // 2 == r) for r in range(6)) for q in Q]
    types = Counter(tuple(sorted([x for x in pp if x], reverse=True)) for pp in prof)
    nv = 15 * 30; idx = lambda i, j: i * 30 + j
    R = []
    for i in range(15): R.append(([idx(i, j) for j in range(30)], [1] * 30, 2))
    for j in range(30): R.append(([idx(i, j) for i in range(15)], [1] * 15, 1))
    for p in range(6):
        for r in range(6):
            cols, vals = [], []
            for i, bp in enumerate(BP):
                if p in bp:
                    for j in range(30):
                        if prof[j][r]:
                            cols.append(idx(i, j)); vals.append(prof[j][r])
            R.append((cols, vals, int(T[p, r])))
    ri = [rr for rr, (c, v, _) in enumerate(R) for _ in c]; ci = [x for c, _, _ in R for x in c]
    vv = [x for _, v, _ in R for x in v]
    A = sp.csr_matrix((np.array(vv, float), (ri, ci)), shape=(len(R), nv)); bb = np.array([t for _, _, t in R], float)
    lp = linprog(np.zeros(nv), A_eq=A, b_eq=bb, bounds=(0, 2), method='highs')
    print('iota #%d profile types' % k, dict(types), '| profile-level LP:', {0: 'feasible', 2: 'INFEASIBLE'}.get(lp.status, lp.status))
