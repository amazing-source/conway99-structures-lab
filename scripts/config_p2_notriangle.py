"""TEST (at most 5 minutes).  Configuration LP with (P2) on random triangle-free pure twisted decompositions."""
import random, sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack
sys.path.insert(0, 'scripts')
from zeros_budget import decompositions
from config_model import columns, build

t0 = time.time()
by_type = {}
for dec in decompositions((1 << 21) - 1, []):
    if not any(len(c) == 3 for c in dec):
        by_type.setdefault(tuple(sorted(len(c) for c in dec)), []).append([tuple(c) for c in dec])
rng = random.Random(260926)
for typ in sorted(by_type):
    if time.time() - t0 > 230:
        print('time budget reached'); break
    dec = rng.choice(by_type[typ])
    cols = columns(dec)
    for p2 in (False, True):
        A, lo, hi, _ = build(dec, p2=p2, cols=cols)
        eq = lo == hi
        res = linprog(np.zeros(A.shape[1]), A_ub=vstack([A[~eq], -A[~eq]]), b_ub=np.concatenate([hi[~eq], -lo[~eq]]),
                      A_eq=A[eq], b_eq=lo[eq], bounds=(0, 1), method='highs', options={'time_limit': 25})
        st = {0: 'feasible', 2: 'INFEASIBLE'}.get(res.status, 'undecided')
        print(f'{typ} ({len(by_type[typ])} labelled), {len(cols)} columns, P2={p2}: LP {st}  [{time.time()-t0:.0f} s]',
              flush=True)
