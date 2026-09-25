"""TEST (seconds).  Is the LP relaxation of the ppe_test model already infeasible for the 8 classes with two
intersecting triangles?  If so, a Farkas certificate (a linear combination of the constraints) exists, and it can be
read as a hand proof."""
import sys
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, 'scripts')
from ppe_test import build

REPS = [
    [(1, 2, 0), (1, 3, 2, 5, 0, 6, 4), (1, 5, 3, 6), (2, 6, 5, 4), (0, 3, 4)],
    [(3, 4, 0), (3, 1, 4, 6, 0, 5), (3, 6, 1, 5, 2), (4, 5, 6, 2), (0, 1, 2)],
    [(1, 2, 0), (1, 5, 2, 6), (1, 3, 5, 0, 6, 4), (2, 3, 6, 5, 4), (0, 3, 4)],
    [(1, 0, 2), (1, 5, 0, 6, 3), (1, 6, 5, 3, 2, 4), (0, 3, 4), (2, 5, 4, 6)],
    [(1, 0, 2), (1, 5, 0, 6, 3), (1, 6, 2, 3, 5, 4), (0, 3, 4), (2, 5, 6, 4)],
    [(1, 0, 2), (1, 5, 0, 6), (1, 3, 2, 5, 6, 4), (0, 3, 4), (2, 6, 3, 5, 4)],
    [(1, 0, 2), (1, 5, 0, 6), (1, 3, 5, 2, 6, 4), (0, 3, 4), (2, 3, 6, 5, 4)],
    [(3, 4, 0), (3, 5, 4, 6, 2), (3, 6, 0, 5, 1), (4, 2, 5, 6, 1), (0, 2, 1)],
]
if __name__ == '__main__':
  for dec in REPS:
    A, lo, hi, nv = build(dec)
    eq = lo == hi
    A_eq, b_eq = A[eq], lo[eq]
    from scipy.sparse import vstack
    Aub = vstack([A[~eq], -A[~eq]])
    bub = np.concatenate([hi[~eq], -lo[~eq]])
    res = linprog(np.zeros(nv), A_ub=Aub, b_ub=bub, A_eq=A_eq, b_eq=b_eq, bounds=(0, 1), method='highs')
    print(tuple(sorted(len(c) for c in dec)), 'LP relaxation status:', res.status, res.message[:60], flush=True)
