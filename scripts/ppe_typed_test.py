"""TEST (at most 5 minutes; orients the maths, proves nothing).  Signed refinement of PPE (types of doubles) on the
one-triangle pure twisted decompositions that the unsigned model (local system + PPE + coupling) left feasible in
ppe_sample_rest.py (same seed, same samples).
Controls:
- the typed LOCAL model alone (no PPE, no coupling) must stay feasible;
- the unsigned model must be feasible again (reproduction)."""
import random
import sys
import time

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions
import ppe_test


def feas(dec, limit, **kw):
    A, lo, hi, nv = ppe_test.build(dec, **kw)
    res = milp(c=np.zeros(nv), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(nv), bounds=Bounds(0, 1),
               options={'time_limit': limit})
    return {0: True, 2: False}.get(res.status, None)


t0 = time.time()
by_type = {}
for dec in decompositions((1 << 21) - 1, []):
    if sum(1 for c in dec if len(c) == 3) <= 1:
        by_type.setdefault(tuple(sorted(len(c) for c in dec)), []).append([tuple(c) for c in dec])
random.seed(20260925)
samples = {typ: random.sample(decs, 6) for typ, decs in sorted(by_type.items(), key=lambda kv: (-kv[0].count(3), kv[0]))}
targets = [((3, 4, 4, 5, 5), 0), ((3, 4, 4, 4, 6), 5), ((3, 4, 4, 5, 5), 2), ((3, 4, 4, 5, 5), 5)]
for typ, idx in targets:
    dec = samples[typ][idx]
    if time.time() - t0 > 200:
        print('time budget reached, stop'); break
    loc = feas(dec, 20, ppe=False, coupling=False, typed=True)
    uns = feas(dec, 30)
    typd = feas(dec, 45, typed=True)
    print(f'{typ} sample {idx}: {dec}\n   typed local only: {loc} ; unsigned PPE+coupling: {uns} ; '
          f'SIGNED PPE+coupling: {typd}   [{time.time() - t0:.0f} s]', flush=True)
