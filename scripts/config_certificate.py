"""Sparse Farkas certificate (L1-minimal) of the configuration LP with (P2), for a given decomposition: which rows
carry weight.  A reading aid for a hand proof."""
import sys
from collections import Counter, defaultdict
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack, hstack, identity, csr_matrix
sys.path.insert(0, 'scripts')
import config_model as cm


def certificate(dec, p2=True):
    A, lo, hi, cols = cm.build(dec, p2=p2)
    labels = list(cm.LABELS)
    nv = A.shape[1]
    eq = lo == hi
    ie, ii = np.flatnonzero(eq), np.flatnonzero(~eq)
    Aeq, beq, Ain, hin, lin_ = A[ie], lo[ie], A[ii], hi[ii], lo[ii]
    ne, ni = len(ie), len(ii)
    M = hstack([Aeq.T, -Aeq.T, Ain.T, -Ain.T, identity(nv)]).tocsr()
    nvar = 2 * ne + 2 * ni + nv
    A_ub = vstack([-M, csr_matrix(np.concatenate([beq, -beq, hin, -lin_, np.ones(nv)]).reshape(1, -1))])
    b_ub = np.concatenate([np.zeros(nv), [-1.0]])
    res = linprog(np.ones(nvar), A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method='highs', options={'time_limit': 120})
    if res.status != 0:
        return None, None, cols
    y = res.x
    w = defaultdict(float)
    for t in range(ne):
        w[labels[ie[t]]] += y[t] - y[ne + t]
    for t in range(ni):
        w[labels[ii[t]] + ' (<=hi)'] += y[2 * ne + t]
        w[labels[ii[t]] + ' (>=lo)'] -= y[2 * ne + ni + t]
    ub = {c: y[2 * ne + 2 * ni + c] for c in range(nv) if y[2 * ne + 2 * ni + c] > 1e-9}
    return {k: v for k, v in w.items() if abs(v) > 1e-9}, ub, cols


if __name__ == '__main__':
    from ppe_typed_lp import SAMPLES
    w, ub, cols = certificate(SAMPLES[int(sys.argv[1]) if len(sys.argv) > 1 else 0])
    if w is None:
        print('no certificate (LP feasible or time limit)'); sys.exit()
    fam = Counter(k.split()[0] for k in w)
    print('rows with weight by family:', dict(fam), '; upper-bound multipliers:', len(ub))
    for k, v in sorted(w.items(), key=lambda kv: (kv[0].split()[0], -abs(kv[1]))):
        if not k.startswith('O '):
            print(f'  {v:+.4f}  {k}')
