"""Is step 2 already infeasible at the BLOCK level?  (ii-block): for each block p and label b,
   sum over the 10 sigma-pairs through p of [b in quadruple assigned] = K2[p+,b] + K2[p-,b]  (= 4 - 2[b in p u m(p)] at a full block).
Variables: w[bp, P] = number of the 2 sigma-pairs of block pair bp assigned to iota-pair P (0..2), row sums 2, column sums 1.
Tested as LP and MILP on the iota's saved by iota_lp.py."""
import json, sys, itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
from fullblock_rigid import mates_from_D
import iota_design as D

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
m = mates_from_D([tuple(c) for c in cfg['cycles']]); K2, K3 = D.K2K3(m, m)
BP = list(itertools.combinations(range(6), 2))
for k, pairs in enumerate(iotas):
    Q = [set(e) | set(f) for e, f in pairs]
    nv = 15 * 30
    idx = lambda i, j: i * 30 + j
    R = []
    for i in range(15): R.append(([idx(i, j) for j in range(30)], [1] * 30, 2))
    for j in range(30): R.append(([idx(i, j) for i in range(15)], [1] * 15, 1))
    for p in range(6):
        for b in range(12):
            cols = [idx(i, j) for i, bp in enumerate(BP) if p in bp for j in range(30) if b in Q[j]]
            R.append((cols, [1] * len(cols), int(K2[2 * p, b] + K2[2 * p + 1, b])))
    ri = [r for r, (c, v, _) in enumerate(R) for _ in c]; ci = [x for c, _, _ in R for x in c]
    A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(R), nv)); bb = np.array([t for _, _, t in R], float)
    lp = linprog(np.zeros(nv), A_eq=A, b_eq=bb, bounds=(0, 2), method='highs')
    ip = milp(np.zeros(nv), constraints=[LinearConstraint(A, bb, bb)], integrality=np.ones(nv), bounds=Bounds(0, 2),
              options={'time_limit': 30})
    print('iota #%d  block-level: LP' % k, {0: 'feasible', 2: 'INFEASIBLE'}.get(lp.status, lp.status),
          '| MILP', {0: 'feasible', 2: 'INFEASIBLE'}.get(ip.status, ip.status), flush=True)
