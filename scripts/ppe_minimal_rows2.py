"""TEST (at most 5 minutes; it guides the hand proofs).  As ppe_minimal_rows.py, but the rows of Lemma A (proved by
hand for every class) are protected: coupling {1,2}-{0,5}, {1,2}-{0,6}, {3,4}-{0,5}, {3,4}-{0,6}, PPE at {1,2} and at
{3,4}.  The deletion filter then finds the extra rows each class needs."""
import sys
import time

import numpy as np

sys.path.insert(0, 'scripts')
import ppe_test
from ppe_lp_relax import REPS
from ppe_minimal_rows import lp_infeasible

PROTECT = {'coupling P=[0, 5] Q=[1, 2]', 'coupling P=[0, 6] Q=[1, 2]', 'coupling P=[0, 5] Q=[3, 4]',
           'coupling P=[0, 6] Q=[3, 4]', 'PPE Q=[1, 2]', 'PPE Q=[3, 4]'}
t0 = time.time()
for w in range(len(REPS)):
    dec = REPS[w]
    A, lo, hi, nv = ppe_test.build(dec)
    labels = list(ppe_test.LABELS)
    mask = np.array([not l.startswith(('colsum', 'PPE', 'coupling')) or l in PROTECT for l in labels])
    base = lp_infeasible(A, lo, hi, nv, mask)
    extra = []
    if not base:
        cand = [r for r, l in enumerate(labels) if l.startswith(('PPE', 'coupling')) and l not in PROTECT]
        for r in cand:
            mask[r] = True
        assert lp_infeasible(A, lo, hi, nv, mask)
        for r in cand:
            mask[r] = False
            if not lp_infeasible(A, lo, hi, nv, mask):
                mask[r] = True
        extra = [labels[r] for r in cand if mask[r]]
    print(f'class {w} {tuple(sorted(len(c) for c in dec))}: Lemma A rows alone infeasible: {base}; extra rows '
          f'({len(extra)}): {extra}  [{time.time() - t0:.0f} s]', flush=True)
