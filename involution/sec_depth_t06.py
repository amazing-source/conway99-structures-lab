"""SEC depth test (structures 27/09). Along true near-designs iota of rel 25 (t06_0): fix k iota-pairs (random order) and test the step-1' LP with and without SEC."""
import json, sys, random
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
sys.path.insert(0, 'C:/home/researcher/Desktop/c99_structures_lab/involution')
import iota_design as D
SP = ''
cfg = json.load(open('t06_0.json')); mu, tau = cfg['mu'], cfg['tau']
def load(fn):
    d = np.load(fn + '.npz'); A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape']))
    return A, d['lo'], d['ub'], json.load(open(fn + '_index.json'))
M = {'sec': load('s1p_t06_rel25'), 'plain': load('s1p_t06_rel25_plain')}
idx = M['sec'][3]; edges = [tuple(e) for e in idx['edges (non-mu pairs)']]
pv = {tuple(p): n for n, p in enumerate(idx['iota_pairs'])}
found = []; rng = random.Random(0)
for it in range(3):
    st, sol, E2 = D.step1(mu, tau, 1, 50, forbid=found)
    if sol is None:
        print('step 1', st); break
    found.append(sol); assert [tuple(e) for e in E2] == edges
    order = list(range(30)); rng.shuffle(order)
    for name, (A, lo, ub, _) in M.items():
        out = []
        for kfix in (0, 4, 8, 12, 16, 20, 25, 30):
            lb = np.zeros(A.shape[1])
            for q in order[:kfix]:
                i, j = sol[q]; lb[pv[(min(i, j), max(i, j))]] = 1
            if (lb > ub).any():
                out.append('%d:X(fix)' % kfix); break
            r = linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=lo, bounds=list(zip(lb, ub)), method='highs')
            out.append('%d:%s' % (kfix, 'F' if r.status == 0 else 'X'))
            if r.status == 2:
                break
        print('iota #%d %-5s' % (it, name), ' '.join(out), flush=True)
