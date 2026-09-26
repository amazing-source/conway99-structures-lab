"""Consistency check of the (P2) rows of config_model.py.  On any solution of the model WITHOUT (P2) (local systems,
cycles, PPE, coupling), the (P2) residuals must sum to ZERO over all unordered pairs {Q1,Q2}: the global identity
sum_{Q1 != Q2} sum_Z val val = sum_Z (144 - 14 - 2 d_Z) = 5376 uses only sum_Q val_Z(Q) = 12 and sum d_Z = 42.
A nonzero sum would reveal an implementation error; individual residuals show where (P2) is violated."""
import itertools
import sys

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

sys.path.insert(0, 'scripts')
import config_model as cm
from ppe_typed_lp import SAMPLES

dec = SAMPLES[int(sys.argv[1]) if len(sys.argv) > 1 else 0]
cols = cm.columns(dec)
A, lo, hi, _ = cm.build(dec, p2=False, cols=cols)
res = milp(c=np.zeros(A.shape[1]), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(A.shape[1]),
           bounds=Bounds(0, 1), options={'time_limit': 200})
print('IP without P2 status:', res.status, res.message[:60])
if res.status == 0:
    x = np.round(res.x).astype(int)
    chosen = [cols[c] for c in range(len(cols)) if x[c]]
    assert len(chosen) == 42
    for Z, p, vals in chosen:
        assert sum(vals.values()) == 12, (Z, sum(vals.values()))
    tot, res_list = 0, []
    for Q1, Q2 in itertools.combinations(cm.PAIRS, 2):
        if len(Q1 & Q2) == 2:
            continue
        lhs = sum(v[Q1] * v[Q2] for Z, p, v in chosen) + sum(v[Q2] for Z, p, v in chosen if frozenset(Z) == Q1)
        rhs = 16 if not (Q1 & Q2) else 12
        tot += lhs - rhs
        res_list.append(lhs - rhs)
    print(f'sum of (P2) residuals over the 210 unordered pairs: {tot}  (must be 0)')
    print(f'nonzero residuals: {sum(1 for r in res_list if r)} of 210; min {min(res_list)}, max {max(res_list)}')
