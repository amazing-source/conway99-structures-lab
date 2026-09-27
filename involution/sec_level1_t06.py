"""Level-1 scan of the secant-first tree (structures 27/09): for each non-pivot label v = p deg of a (0,6) type and each
candidate secant iota-pair centred at v, fix it alone in the step-1' LP with SEC and report how many stay feasible."""
import json, sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
rel = int(sys.argv[1]) if len(sys.argv) > 1 else 25
base = 's1p_t06_rel%d' % rel
d = np.load(base + '.npz'); idx = json.load(open(base + '_index.json'))
A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape'])); lo, ub = d['lo'], d['ub']
mu, tau = idx['mu'], idx['tau']
edges = [tuple(e) for e in idx['edges (non-mu pairs)']]; pairs = [tuple(p) for p in idx['iota_pairs']]
nonpiv = [l for l in range(12) if mu[l] != tau[l]]
t0 = time.time()
for v in nonpiv:
    cands = [n for n, (i, j) in enumerate(pairs) if set(edges[i]) & set(edges[j]) == {v}]
    feas = []
    for n in cands:
        lb = np.zeros(A.shape[1]); lb[n] = 1
        r = linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=lo, bounds=list(zip(lb, ub)), method='highs')
        if r.status == 0:
            feas.append(n)
    print('rel %d  non-pivot %2d (block %d): %d candidate secant pairs, %d LP-feasible  (%.0fs)'
          % (rel, v, v // 2, len(cands), len(feas), time.time() - t0), flush=True)
