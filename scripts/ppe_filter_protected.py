"""TEST (short; guides the hand proof).  Deletion filter with a protected set: which extra PPE / coupling rows does
the class need beyond the protected rows?"""
import sys

import numpy as np

sys.path.insert(0, 'scripts')
import ppe_test
from ppe_lp_relax import REPS
from ppe_minimal_rows import lp_infeasible
from ppe_certificate import LEMMA_A

w = int(sys.argv[1])
prot = set(LEMMA_A) | set(sys.argv[2:])
A, lo, hi, nv = ppe_test.build(REPS[w])
labels = list(ppe_test.LABELS)
mask = np.array([not l.startswith(('colsum', 'PPE', 'coupling')) or l.startswith(('PPE', 'coupling')) for l in labels])
assert lp_infeasible(A, lo, hi, nv, mask)
cand = [r for r, l in enumerate(labels) if l.startswith(('PPE', 'coupling')) and l not in prot]
cand.sort(key=lambda r: (labels[r].startswith('PPE'), labels[r]))      # try couplings first, keep PPE if possible
for r in cand:
    mask[r] = False
    if not lp_infeasible(A, lo, hi, nv, mask):
        mask[r] = True
print(f'class {w}: extra rows beyond the protected set ({sum(mask[r] for r in cand)}):',
      [labels[r] for r in cand if mask[r]])
