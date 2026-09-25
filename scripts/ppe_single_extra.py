"""TEST (short; guides the hand proofs).  For a class, protect the rows of Lemma A plus a given list, then report
which SINGLE extra row (coupling or PPE) already makes the LP relaxation infeasible."""
import sys

import numpy as np

sys.path.insert(0, 'scripts')
import ppe_test
from ppe_lp_relax import REPS
from ppe_minimal_rows import lp_infeasible
from ppe_certificate import LEMMA_A

w = int(sys.argv[1])
prot = set(LEMMA_A) | set(sys.argv[2:])
dec = REPS[w]
A, lo, hi, nv = ppe_test.build(dec)
labels = list(ppe_test.LABELS)
base = np.array([not l.startswith(('colsum', 'PPE', 'coupling')) or l in prot for l in labels])
print('protected rows alone infeasible:', lp_infeasible(A, lo, hi, nv, base))
hits = []
for r, l in enumerate(labels):
    if l.startswith(('PPE', 'coupling')) and l not in prot:
        m = base.copy()
        m[r] = True
        if lp_infeasible(A, lo, hi, nv, m):
            hits.append(l)
print(f'class {w}: single extra rows that suffice ({len(hits)}):', hits)
