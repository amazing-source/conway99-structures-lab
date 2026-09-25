"""TEST (at most 5 minutes; orients the maths).  Add C1's off-diagonal (P2) identities to the orbit model
(local system + PPE + coupling), on the one-triangle samples the model left feasible.
Controls:
- the classes already killed stay infeasible (trivially);
- the model with P2 alone (without PPE and coupling) is run first."""
import sys
import time
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import vstack
sys.path.insert(0, 'scripts')
import ppe_test
from ppe_typed_lp import SAMPLES


def lp(dec, **kw):
    A, lo, hi, nv = ppe_test.build(dec, **kw)
    eq = lo == hi
    r = linprog(np.zeros(nv), A_ub=vstack([A[~eq], -A[~eq]]), b_ub=np.concatenate([hi[~eq], -lo[~eq]]),
                A_eq=A[eq], b_eq=lo[eq], bounds=(0, 1), method='highs')
    return r.status == 0, nv, A.shape[0]


def ip(dec, limit, **kw):
    A, lo, hi, nv = ppe_test.build(dec, **kw)
    r = milp(c=np.zeros(nv), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(nv), bounds=Bounds(0, 1),
             options={'time_limit': limit})
    return {0: True, 2: False}.get(r.status, None)


t0 = time.time()
for dec in SAMPLES:
    f1, nv, nr = lp(dec, ppe=False, coupling=False, p2=True)
    f2, nv2, nr2 = lp(dec, p2=True)
    print(tuple(sorted(len(c) for c in dec)), f'LP: P2 alone feasible {f1}; P2 + PPE + coupling feasible {f2} '
          f'({nv2} vars, {nr2} rows)  [{time.time() - t0:.0f} s]', flush=True)
for dec in SAMPLES[:2]:
    if time.time() - t0 > 170:
        break
    print(tuple(sorted(len(c) for c in dec)), 'MILP P2 + PPE + coupling:', ip(dec, 60, p2=True),
          f'[{time.time() - t0:.0f} s]', flush=True)
