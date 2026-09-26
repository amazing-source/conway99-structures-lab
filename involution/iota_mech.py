"""Which sub-family of the block-level (ii) constraints already makes step 2 LP-infeasible (for saved iota's)?
 c(Q, r) = [r in assigned block pair of Q].  Families:
   OWN : for each label b, sum_{Q containing b} c(Q, blk b)  = 2
   MATE: for each label b, sum_{Q containing b} c(Q, blk m b) = 2
   REST: for each label b and other block r, sum_{Q containing b} c(Q, r) = 4"""
import json, sys, itertools
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
from fullblock_rigid import mates_from_D

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
m = mates_from_D([tuple(c) for c in cfg['cycles']])
BP = list(itertools.combinations(range(6), 2))
for k, pairs in enumerate(iotas):
    Q = [set(e) | set(f) for e, f in pairs]
    idx = lambda i, j: i * 30 + j
    base = []
    for i in range(15): base.append(([idx(i, j) for j in range(30)], 2))
    for j in range(30): base.append(([idx(i, j) for i in range(15)], 1))
    fam = {'OWN': [], 'MATE': [], 'REST': []}
    for b in range(12):
        for r in range(6):
            cols = [idx(i, j) for i, bp in enumerate(BP) if r in bp for j in range(30) if b in Q[j]]
            if r == b // 2: fam['OWN'].append((cols, 2))
            elif r == m[b] // 2: fam['MATE'].append((cols, 2))
            else: fam['REST'].append((cols, 4))
    res = {}
    for name in ('OWN', 'MATE', 'OWN+MATE', 'OWN+REST', 'MATE+REST', 'ALL'):
        R = base + sum((fam[f] for f in (['OWN', 'MATE', 'REST'] if name == 'ALL' else name.split('+'))), [])
        ri = [r for r, (c, _) in enumerate(R) for _ in c]; ci = [x for c, _ in R for x in c]
        A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(R), 450)); bb = np.array([t for _, t in R], float)
        lp = linprog(np.zeros(450), A_eq=A, b_eq=bb, bounds=(0, 2), method='highs')
        res[name] = {0: 'feas', 2: 'INF'}.get(lp.status, lp.status)
    print('iota #%d' % k, res)
