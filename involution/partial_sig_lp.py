"""Partial-fixing test on the strengthened model (export_sigma_model): fix k iota-pairs of saved near-designs (p = 1)."""
import json, sys, time, random
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
base = sys.argv[1]
d = np.load(base + '.npz'); A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape']))
idx = json.load(open(base + '_index.json'))
ny = len(idx['y']); pk = [tuple(map(tuple, q)) for q in idx['p']]; pidx = {q: ny + i for i, q in enumerate(pk)}
eq = d['lo'] == d['hi']
Aeq, beq, Aub, bub = A[eq], d['lo'][eq], A[~eq], d['hi'][~eq]
iotas = json.load(open(sys.argv[2]))
rng = random.Random(0)
for k_iota, pairs in enumerate(iotas):
    qs = [tuple(sorted((tuple(e), tuple(f)))) for e, f in pairs]
    missing = [q for q in qs if q not in pidx]
    print('iota #%d: %d of its 30 pairs are excluded by the pointwise filters / lemma' % (k_iota, len(missing)), flush=True)
    order = list(range(30)); rng.shuffle(order)
    for kfix in (4, 6, 8, 10, 12):
        lb = np.zeros(A.shape[1]); ok = True
        for i in order[:kfix]:
            if qs[i] not in pidx:
                ok = False; break
            lb[pidx[qs[i]]] = 1
        if not ok:
            print('   %2d fixed: a fixed pair is already excluded' % kfix); break
        t0 = time.time()
        r = linprog(np.zeros(A.shape[1]), A_eq=Aeq, b_eq=beq, A_ub=Aub, b_ub=bub, bounds=list(zip(lb, np.ones(A.shape[1]))), method='highs')
        print('   %2d fixed: LP' % kfix, {0: 'feasible', 2: 'INFEASIBLE'}.get(r.status, r.status), '%.1fs' % (time.time() - t0), flush=True)
        if r.status == 2:
            break
