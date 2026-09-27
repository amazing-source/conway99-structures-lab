"""Level-2 sample of the secant-first tree (structures 27/09): a few LP-feasible secants at the first non-pivot label, then
every candidate secant at the second non-pivot label; survival rate estimates the tree size.  Budget: stops at 270 s."""
import json, sys, time, random
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
rel = int(sys.argv[1]) if len(sys.argv) > 1 else 25
v1, v2 = int(sys.argv[2]), int(sys.argv[3])
base = 's1p_t06_rel%d' % rel
d = np.load(base + '.npz'); idx = json.load(open(base + '_index.json'))
A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape'])); lo, ub = d['lo'], d['ub']
edges = [tuple(e) for e in idx['edges (non-mu pairs)']]; pairs = [tuple(p) for p in idx['iota_pairs']]
cand = lambda v: [n for n, (i, j) in enumerate(pairs) if set(edges[i]) & set(edges[j]) == {v}]
def feas(fix):
    lb = np.zeros(A.shape[1]); lb[list(fix)] = 1
    return linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=lo, bounds=list(zip(lb, ub)), method='highs').status == 0
rng = random.Random(7); t0 = time.time()
c1 = cand(v1); rng.shuffle(c1)
tested = 0
for n1 in c1:
    if time.time() - t0 > 200:
        break
    if not feas([n1]):
        continue
    alive = 0; tot = 0
    for n2 in cand(v2):
        if time.time() - t0 > 270:
            break
        tot += 1; alive += feas([n1, n2])
    print('rel %d: secant %s at %d (feasible) -> at %d: %d / %d children feasible  (%.0fs)'
          % (rel, pairs[n1], v1, v2, alive, tot, time.time() - t0), flush=True)
    tested += 1
    if tested == 3:
        break
