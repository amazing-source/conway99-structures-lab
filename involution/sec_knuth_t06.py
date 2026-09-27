"""Knuth-style size estimate of the secant-first LP tree for a (0,6) type (structures 27/09): along random paths, at each
level (one non-pivot label, in order) sample up to 12 candidate secant iota-pairs, estimate the number of LP-feasible
children as (fraction feasible) x (number of candidates compatible with the fixed ones), descend into a random feasible one.
Estimate = sum over levels of the product of the branching estimates. Budget ~270 s."""
import json, sys, time, random
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
rel = int(sys.argv[1]) if len(sys.argv) > 1 else 25
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
base = 's1p_t06_rel%d' % rel
d = np.load(base + '.npz'); idx = json.load(open(base + '_index.json'))
A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape'])); lo, ub = d['lo'], d['ub']
mu, tau = idx['mu'], idx['tau']
edges = [tuple(e) for e in idx['edges (non-mu pairs)']]; pairs = [tuple(p) for p in idx['iota_pairs']]
nonpiv = [l for l in range(12) if mu[l] != tau[l]]
cand = {v: [n for n, (i, j) in enumerate(pairs) if set(edges[i]) & set(edges[j]) == {v}] for v in nonpiv}
def feas(fix):
    lb = np.zeros(A.shape[1]); lb[list(fix)] = 1
    return linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=lo, bounds=list(zip(lb, ub)), method='highs').status == 0
def compatible(n, fixed):         # iota is a matching of the non-mu pairs: no non-mu pair used twice
    used = {e for m in fixed for e in pairs[m]}
    return not (set(pairs[n]) & used)
rng = random.Random(seed); t0 = time.time()
while time.time() - t0 < 200:
    fixed, prod, est, path = [], 1.0, 1.0, []
    for lev, v in enumerate(nonpiv):
        cs = [n for n in cand[v] if compatible(n, fixed)]
        rng.shuffle(cs)
        samp = cs[:12]
        ok = [n for n in samp if feas(fixed + [n])]
        b = len(cs) * len(ok) / max(len(samp), 1)
        path.append('%d:%d/%d of %d' % (lev + 1, len(ok), len(samp), len(cs)))
        prod *= b
        est += prod
        if not ok:
            break
        fixed.append(rng.choice(ok))
    print('rel %d path: %s | estimated tree nodes %.3g  (%.0fs)' % (rel, ' '.join(path), est, time.time() - t0), flush=True)
