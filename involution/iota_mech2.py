"""Finer mechanism search on a saved near-design: find a SMALL infeasible subsystem of the block-level (ii) constraints
(OWN/MATE/REST rows for specific labels b), by greedy deletion (deletion filter for an IIS on the constraint families per label)."""
import json, sys, itertools
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
from fullblock_rigid import mates_from_D

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
m = mates_from_D([tuple(c) for c in cfg['cycles']])
BP = list(itertools.combinations(range(6), 2))
pairs = iotas[k]
Q = [set(e) | set(f) for e, f in pairs]
idx = lambda i, j: i * 30 + j
base = []
for i in range(15): base.append(([idx(i, j) for j in range(30)], 2))
for j in range(30): base.append(([idx(i, j) for i in range(15)], 1))
rows = {}
for b in range(12):
    for r in range(6):
        cols = [idx(i, j) for i, bp in enumerate(BP) if r in bp for j in range(30) if b in Q[j]]
        t = 2 if r in (b // 2, m[b] // 2) else 4
        rows[(b, r)] = (cols, t)

def feas(keys):
    R = base + [rows[x] for x in keys]
    ri = [rr for rr, (c, _) in enumerate(R) for _ in c]; ci = [x for c, _ in R for x in c]
    A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(R), 450)); bb = np.array([t for _, t in R], float)
    return linprog(np.zeros(450), A_eq=A, b_eq=bb, bounds=(0, 2), method='highs').status == 0

keys = list(rows)
assert not feas(keys)
core = list(keys)
for x in list(keys):
    trial = [y for y in core if y != x]
    if not feas(trial):
        core = trial
print('irreducible infeasible set of (label b, block r) rows:', len(core))
for b, r in sorted(core):
    kind = 'OWN' if r == b // 2 else ('MATE' if r == m[b] // 2 else 'REST')
    print('  label', b, '(block %d)' % (b // 2), 'block r =', r, kind, 'target', rows[(b, r)][1])

# exact Farkas certificate on the IIS
import sys as _s
_s.path.insert(0, '../ad')
from certify import certify, check
from fractions import Fraction
R = base + [rows[x] for x in core]
names = ['pair%s=2' % (BP[i],) for i in range(15)] + ['quad%d%s=1' % (j, sorted(Q[j])) for j in range(30)] + \
        ['N(b=%d,r=%d)=%d' % (b, r, rows[(b, r)][1]) for b, r in core]
ri = [rr for rr, (c, _) in enumerate(R) for _ in c]; ci = [x for c, _ in R for x in c]
A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(R), 450)); bb = np.array([t for _, t in R], float)
res = certify(A, bb, bb, np.zeros(450), np.full(450, 2.0), maxnodes=1, denom=60)
print('certificate:', res['status'], res.get('level'))
if res['status'] == 'dead':
    y = [Fraction(v) for v in res['tree']['cert']]
    print('check', check(A, bb, bb, np.zeros(450), np.full(450, 2.0), res['tree']))
    for n_, v in zip(names, y):
        if v != 0:
            print('   %-45s %s' % (n_, v))
