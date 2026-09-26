"""Farkas certificate of the block-level step-2 LP for a saved iota; print its support to understand the mechanism."""
import json, sys, itertools
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, '../ad')
from certify import certify, check
from fullblock_rigid import mates_from_D
import iota_design as D

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
m = mates_from_D([tuple(c) for c in cfg['cycles']]); K2, K3 = D.K2K3(m, m)
BP = list(itertools.combinations(range(6), 2))
pairs = iotas[k]
Q = [set(e) | set(f) for e, f in pairs]
nv = 15 * 30
idx = lambda i, j: i * 30 + j
R, names = [], []
for i in range(15): R.append(([idx(i, j) for j in range(30)], 2)); names.append('pair%s' % (BP[i],))
for j in range(30): R.append(([idx(i, j) for i in range(15)], 1)); names.append('quad%d%s' % (j, sorted(Q[j])))
for p in range(6):
    for b in range(12):
        cols = [idx(i, j) for i, bp in enumerate(BP) if p in bp for j in range(30) if b in Q[j]]
        R.append((cols, int(K2[2 * p, b] + K2[2 * p + 1, b]))); names.append('blk%d_lab%d' % (p, b))
ri = [r for r, (c, _) in enumerate(R) for _ in c]; ci = [x for c, _ in R for x in c]
A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(R), nv)); bb = np.array([t for _, t in R], float)
lb = np.zeros(nv); ub = np.full(nv, 2.0)
res = certify(A, bb, bb, lb, ub, maxnodes=1)
print('status', res['status'], res.get('level'))
if res['status'] == 'dead' and res['level'] == 'LP':
    y = res['tree']['cert']
    y = [Fraction(v) for v in y]
    print('check:', check(A, bb, bb, lb, ub, res['tree']))
    nz = [(names[r], y[r]) for r in range(len(R)) if y[r] != 0]
    print('nonzero multipliers:', len(nz))
    for n_, v in nz:
        print('  ', n_, v)
