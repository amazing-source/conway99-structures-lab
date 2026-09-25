"""TEST (seconds).  LP relaxation of the signed model on the three one-triangle samples of ppe_typed_test.py."""
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack
sys.path.insert(0, 'scripts')
import ppe_test

SAMPLES = [[(0, 1, 5, 3), (0, 2, 3, 6, 5), (0, 4, 5, 2, 6), (1, 2, 4), (1, 3, 4, 6)],
           [(0, 1, 2, 5), (0, 2, 3, 5, 1, 4), (0, 3, 6), (1, 3, 4, 6), (2, 4, 5, 6)],
           [(0, 1, 6, 5), (0, 2, 4, 3), (0, 4, 6), (1, 2, 6, 3, 5), (1, 3, 2, 5, 4)]]
for dec in SAMPLES:
    for typed in (False, True):
        A, lo, hi, nv = ppe_test.build(dec, typed=typed)
        eq = lo == hi
        res = linprog(np.zeros(nv), A_ub=vstack([A[~eq], -A[~eq]]), b_ub=np.concatenate([hi[~eq], -lo[~eq]]),
                      A_eq=A[eq], b_eq=lo[eq], bounds=(0, 1), method='highs')
        print(tuple(sorted(len(c) for c in dec)), 'typed' if typed else 'unsigned', 'LP relaxation feasible:',
              res.status == 0, flush=True)
