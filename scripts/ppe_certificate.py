"""Sparse Farkas certificate of the LP relaxation restricted to: the local system, the free 6-cycles, and a GIVEN list
of PPE / coupling rows (the minimal ones found by ppe_minimal_rows2.py).  A reading aid for the hand proofs, not a proof."""
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack, hstack, identity, csr_matrix

sys.path.insert(0, 'scripts')
import ppe_test
from ppe_lp_relax import REPS

LEMMA_A = ['coupling P=[0, 5] Q=[1, 2]', 'coupling P=[0, 6] Q=[1, 2]', 'coupling P=[0, 5] Q=[3, 4]',
           'coupling P=[0, 6] Q=[3, 4]', 'PPE Q=[1, 2]', 'PPE Q=[3, 4]']


def certificate(dec, keep):
    A, lo, hi, nv = ppe_test.build(dec)
    labels = list(ppe_test.LABELS)
    act = np.array([not l.startswith(('colsum', 'PPE', 'coupling')) or l in keep for l in labels])
    idx = np.flatnonzero(act)
    A, lo, hi = A[idx], lo[idx], hi[idx]
    labels = [labels[r] for r in idx]
    eq = lo == hi
    ie, ii = np.flatnonzero(eq), np.flatnonzero(~eq)
    Aeq, beq, Ain, hin, lin_ = A[ie], lo[ie], A[ii], hi[ii], lo[ii]
    ne, ni = len(ie), len(ii)
    M = hstack([Aeq.T, -Aeq.T, Ain.T, -Ain.T, identity(nv)]).tocsr()
    nvar = 2 * ne + 2 * ni + nv
    A_ub = vstack([-M, csr_matrix(np.concatenate([beq, -beq, hin, -lin_, np.ones(nv)]).reshape(1, -1))])
    b_ub = np.concatenate([np.zeros(nv), [-1.0]])
    res = linprog(np.ones(nvar), A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method='highs')
    if res.status != 0:
        return None
    y = res.x
    w = defaultdict(float)
    for t in range(ne):
        w[labels[ie[t]]] += y[t] - y[ne + t]
    for t in range(ni):
        w[labels[ii[t]] + ' (<=hi)'] += y[2 * ne + t]
        w[labels[ii[t]] + ' (>=lo)'] -= y[2 * ne + ni + t]
    ub = {ppe_test.VARS[t]: y[2 * ne + 2 * ni + t] for t in range(nv) if y[2 * ne + 2 * ni + t] > 1e-9}
    return res.fun, {k: v for k, v in w.items() if abs(v) > 1e-9}, ub


if __name__ == '__main__':
    w = int(sys.argv[1])
    extra = sys.argv[2:]
    out = certificate(REPS[w], LEMMA_A + extra)
    if out is None:
        print('restricted LP is feasible: no certificate')
        sys.exit()
    val, cert, ub = out
    print(f'class {w}: {REPS[w]}; L1 = {val:.3f}; {len(cert)} rows, {len(ub)} upper bounds')
    for k, v in sorted(cert.items(), key=lambda kv: kv[0]):
        print(f'  {v:+.4f}  {k}')
    for key, v in sorted(ub.items(), key=str):
        print(f'  ub {v:.4f}  {key}')
