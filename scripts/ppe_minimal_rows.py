"""TEST (at most 5 minutes; it guides the hand proofs and proves nothing).  For each class with two intersecting
triangles: keep the local system and the free 6-cycles, and find a minimal set of PPE / coupling rows that keeps the LP
relaxation infeasible.  Method: deletion filter, one LP per row; column sums are dropped first."""
import sys
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack

sys.path.insert(0, 'scripts')
import ppe_test
from ppe_lp_relax import REPS


def lp_infeasible(A, lo, hi, nv, mask):
    Aa, la, ha = A[mask], lo[mask], hi[mask]
    eq = la == ha
    Aub = vstack([Aa[~eq], -Aa[~eq]])
    bub = np.concatenate([ha[~eq], -la[~eq]])
    res = linprog(np.zeros(nv), A_ub=Aub, b_ub=bub, A_eq=Aa[eq], b_eq=la[eq], bounds=(0, 1), method='highs')
    return res.status == 2


if __name__ == '__main__':
  t0 = time.time()
  which = [int(a) for a in sys.argv[1:]] or range(len(REPS))
  for w in which:
    dec = REPS[w]
    A, lo, hi, nv = ppe_test.build(dec)
    labels = list(ppe_test.LABELS)
    mask = np.ones(len(labels), dtype=bool)
    removable = [r for r, l in enumerate(labels) if l.startswith(('colsum', 'PPE', 'coupling'))]
    removable.sort(key=lambda r: (not labels[r].startswith('colsum'), labels[r]))
    assert lp_infeasible(A, lo, hi, nv, mask)
    for r in removable:
        mask[r] = False
        if not lp_infeasible(A, lo, hi, nv, mask):
            mask[r] = True
    kept = [labels[r] for r in removable if mask[r]]
    print(f'class {w} {tuple(sorted(len(c) for c in dec))} {dec}\n   minimal rows ({len(kept)}): {kept}  '
          f'[{time.time() - t0:.0f} s]', flush=True)
