"""TEST (short; guides the hand proof of class b).  With a protected set of rows, try every PAIR of extra PPE rows and
report the pairs that make the LP relaxation infeasible."""
import itertools
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
base = np.array([not l.startswith(('colsum', 'PPE', 'coupling')) or l in prot for l in labels])
pool = [r for r, l in enumerate(labels) if l.startswith('PPE') and l not in prot]
hits = []
for r1, r2 in itertools.combinations(pool, 2):
    m = base.copy()
    m[r1] = m[r2] = True
    if lp_infeasible(A, lo, hi, nv, m):
        hits.append((labels[r1], labels[r2]))
print(f'class {w}: pairs of extra PPE rows that suffice ({len(hits)}):')
for h in hits:
    print('  ', h)
