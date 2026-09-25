"""TEST (at most 5 minutes; it orients the maths and proves nothing).  Model of ppe_test.py (PPE + coupling + local
system), run on a random sample of labelled pure twisted decompositions with at most one triangle, 6 per type."""
import random
import sys
import time

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions
from ppe_test import feasible

t0 = time.time()
by_type = {}
for dec in decompositions((1 << 21) - 1, []):
    nt = sum(1 for c in dec if len(c) == 3)
    if nt <= 1:
        by_type.setdefault(tuple(sorted(len(c) for c in dec)), []).append([tuple(c) for c in dec])
random.seed(20260925)
print(f'enumerated [{time.time() - t0:.0f} s]', flush=True)
for typ, decs in sorted(by_type.items(), key=lambda kv: (-kv[0].count(3), kv[0])):
    sample = random.sample(decs, 6)
    res = []
    for dec in sample:
        if time.time() - t0 > 260:
            res.append('skipped')
            continue
        res.append(feasible(dec))
    print(f'{typ}: {len(decs)} labelled; sample of 6, feasible with PPE + coupling: {res}  [{time.time() - t0:.0f} s]',
          flush=True)
