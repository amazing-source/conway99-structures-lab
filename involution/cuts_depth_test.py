"""Does the lemma-cut version of the step-1' LP die at a smaller depth than the plain one (partial near-designs)?"""
import json, sys, time, random, itertools
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import iota_design as D
cfg = json.load(open(sys.argv[1])); base = sys.argv[2]
mu, tau = cfg['mu'], cfg['tau']


def load(fn):
    d = np.load(fn + '.npz'); A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape']))
    return A, d['lo'], d['hi'], json.load(open(fn + '_index.json'))
M = {'plain': load(base), 'cuts': load(base + '_cuts')}
idx = M['plain'][3]
edges = [tuple(e) for e in idx['edges (non-mu pairs)']]
pv = {tuple(p): n for n, p in enumerate(idx['iota_pairs'])}
found = []; rng = random.Random(0)
for k in range(2):
    st, sol, E2 = D.step1(mu, tau, 1, 60, forbid=found)
    if sol is None:
        print('step 1', st); break
    found.append(sol)
    assert [tuple(e) for e in E2] == edges
    order = list(range(30)); rng.shuffle(order)
    for name, (A, lo, hi, _) in M.items():
        eq = lo == hi
        res = []
        for kfix in (4, 6, 8, 10, 12):
            lb = np.zeros(A.shape[1])
            for q in order[:kfix]:
                i, j = sol[q]; lb[pv[(min(i, j), max(i, j))]] = 1
            r = linprog(np.zeros(A.shape[1]), A_eq=A[eq], b_eq=lo[eq], A_ub=A[~eq] if (~eq).any() else None,
                        b_ub=hi[~eq] if (~eq).any() else None, bounds=list(zip(lb, np.ones(A.shape[1]))), method='highs')
            res.append('%d:%s' % (kfix, 'F' if r.status == 0 else 'X'))
        print('iota #%d %-5s' % (k, name), ' '.join(res), flush=True)
