"""Sparse Farkas certificate for the LP relaxation of the ppe_test model (a check that orients the hand proof).
Find y with y^T G >= 0 on the variables and y^T h <= -1, minimising the L1 norm of y, so that the certificate is
short enough to be read as a hand argument.  Prints the rows with nonzero weight, grouped by family."""
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack, hstack, identity, csr_matrix

sys.path.insert(0, 'scripts')
import ppe_test
from ppe_lp_relax import REPS


def farkas(dec):
    A, lo, hi, nv = ppe_test.build(dec)
    labels = list(ppe_test.LABELS)
    eq = lo == hi
    idx_eq = np.flatnonzero(eq)
    idx_in = np.flatnonzero(~eq)
    Aeq, beq = A[idx_eq], lo[idx_eq]
    Ain, hin, lin_ = A[idx_in], hi[idx_in], lo[idx_in]
    ne, ni = len(idx_eq), len(idx_in)
    M = hstack([Aeq.T, -Aeq.T, Ain.T, -Ain.T, identity(nv)]).tocsr()
    nvar = 2 * ne + 2 * ni + nv
    A_ub = vstack([-M, csr_matrix(np.concatenate([beq, -beq, hin, -lin_, np.ones(nv)]).reshape(1, -1))])
    b_ub = np.concatenate([np.zeros(nv), [-1.0]])
    res = linprog(np.ones(nvar), A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method='highs')
    y = res.x
    w = defaultdict(float)
    for t in range(ne):
        w[labels[idx_eq[t]]] += y[t] - y[ne + t]
    for t in range(ni):
        w[labels[idx_in[t]] + ' (<= hi)'] += y[2 * ne + t]
        w[labels[idx_in[t]] + ' (>= lo)'] -= y[2 * ne + ni + t]
    ub = {ppe_test.VARS[t]: y[2 * ne + 2 * ni + t] for t in range(nv) if y[2 * ne + 2 * ni + t] > 1e-9}
    return res.fun, {k: v for k, v in w.items() if abs(v) > 1e-9}, ub


if __name__ == '__main__':
    which = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    dec = REPS[which]
    val, w, ub = farkas(dec)
    print('decomposition', dec, ' L1 norm', round(val, 4))
    fam = Counter()
    for k in w:
        fam[k.split()[0]] += 1
    print('rows with nonzero weight by family:', dict(fam), ' upper-bound multipliers:', len(ub))
    for k, v in sorted(w.items(), key=lambda kv: (kv[0].split()[0], kv[0])):
        print(f'  {v:+.4f}  {k}')
    for key, v in sorted(ub.items(), key=str):
        print(f'  x<=1 weight {v:.4f} on {key}')
