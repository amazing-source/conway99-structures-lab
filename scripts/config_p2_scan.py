"""TEST (at most 5 minutes).  LP relaxation of the configuration model WITH C1's (P2), on the three one-triangle
samples and on random samples of every one-triangle and triangle-free type (seed fixed).  An LP-infeasible sample
has a Farkas certificate: a candidate hand proof."""
import random
import sys
import time

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions
from config_model import columns, build, lp_feasible
from ppe_typed_lp import SAMPLES

t0 = time.time()
for t, dec in enumerate(SAMPLES):
    cols = columns(dec)
    A, lo, hi, _ = build(dec, p2=True, cols=cols)
    print(f'sample {t}: LP with P2 feasible: {lp_feasible(A, lo, hi)}  [{time.time() - t0:.0f} s]', flush=True)
by_type = {}
for dec in decompositions((1 << 21) - 1, []):
    if sum(1 for c in dec if len(c) == 3) <= 1:
        by_type.setdefault(tuple(sorted(len(c) for c in dec)), []).append([tuple(c) for c in dec])
rng = random.Random(260926)
for typ in sorted(by_type, key=lambda t: (-t.count(3), t)):
    if time.time() - t0 > 250:
        print('time budget reached'); break
    dec = rng.choice(by_type[typ])
    cols = columns(dec)
    A, lo, hi, _ = build(dec, p2=True, cols=cols)
    print(f'{typ}: random sample, {len(cols)} columns: LP with P2 feasible: {lp_feasible(A, lo, hi)}  '
          f'[{time.time() - t0:.0f} s]', flush=True)
