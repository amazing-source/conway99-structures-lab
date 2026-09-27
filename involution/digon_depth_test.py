import json, sys, random
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import iota_design as D
cfg = json.load(open('tier8_0.json')); mu, tau = cfg['mu'], cfg['tau']
def load(fn):
    d = np.load(fn + '.npz'); A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape']))
    return A, d['lo'], d['hi'], d['ub'], json.load(open(fn + '_index.json'))
M = {'plain': load('s1p_tier8_0'), 'digon': load('s1p_tier8_0_digon')}
idx = M['plain'][4]; edges = [tuple(e) for e in idx['edges (non-mu pairs)']]
pv = {tuple(p): n for n, p in enumerate(idx['iota_pairs'])}
r0 = linprog(np.zeros(M['digon'][0].shape[1]), A_eq=M['digon'][0], b_eq=M['digon'][1], bounds=list(zip(np.zeros(len(M['digon'][3])), M['digon'][3])), method='highs')
print('digon model root LP:', {0: 'feasible', 2: 'INFEASIBLE'}.get(r0.status, r0.status), flush=True)
found = []; rng = random.Random(0)
for k in range(3):
    st, sol, E2 = D.step1(mu, tau, 1, 60, forbid=found)
    if sol is None: print('step1', st); break
    found.append(sol); order = list(range(30)); rng.shuffle(order)
    for name, (A, lo, hi, ub, _) in M.items():
        out = []
        for kfix in (2, 4, 6, 8, 10, 12):
            lb = np.zeros(A.shape[1])
            for q in order[:kfix]:
                i, j = sol[q]; lb[pv[(min(i, j), max(i, j))]] = 1
            if (lb > ub).any():
                out.append('%d:X(fix)' % kfix); break
            r = linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=lo, bounds=list(zip(lb, ub)), method='highs')
            out.append('%d:%s' % (kfix, 'F' if r.status == 0 else 'X'))
            if r.status == 2: break
        print('iota #%d %-5s' % (k, name), ' '.join(out), flush=True)
