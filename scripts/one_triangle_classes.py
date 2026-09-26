"""TEST (at most 5 minutes).  Isomorphism classes of the pure twisted decompositions with exactly one triangle
(canonical form over the 144 relabellings fixing the triangle setwise), and the LP relaxation of the configuration
model with (P2) on one representative per class (20 s limit per LP)."""
import itertools
import sys
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions, V
from config_model import columns, build


def canon1(dec):
    (T,) = [set(c) for c in dec if len(c) == 3]
    rest = sorted(set(V) - T)
    best = None
    for a in itertools.permutations(sorted(T)):
        for b in itertools.permutations(rest):
            mp = dict(zip(a + b, range(7)))
            form = tuple(sorted(tuple(sorted(tuple(sorted((mp[c[t]], mp[c[(t + 1) % len(c)]])))
                                             for t in range(len(c)))) for c in dec))
            if best is None or form < best[0]:
                best = (form, [tuple(mp[x] for x in c) for c in dec])
    return best


if __name__ == "__main__":
  t0 = time.time()
  classes = {}
  for dec in decompositions((1 << 21) - 1, []):
      if sum(1 for c in dec if len(c) == 3) != 1:
          continue
      form, rep = canon1(dec)
      typ = tuple(sorted(len(c) for c in dec))
      classes.setdefault((typ, form), [rep, 0])[1] += 1
  print('(subtour cuts on)', f'{len(classes)} classes of one-triangle decompositions ({sum(v[1] for v in classes.values())} labelled) '
        f'[{time.time() - t0:.0f} s]', flush=True)
  dead = alive = unk = 0
  for (typ, form), (rep, cnt) in sorted(classes.items()):
      if time.time() - t0 > 265:
          print('time budget reached'); break
      A, lo, hi, cols = build(rep, p2=True)
      eq = lo == hi
      res = linprog(np.zeros(A.shape[1]), A_ub=vstack([A[~eq], -A[~eq]]), b_ub=np.concatenate([hi[~eq], -lo[~eq]]),
                    A_eq=A[eq], b_eq=lo[eq], bounds=(0, 1), method='highs', options={'time_limit': 20})
      st = {0: 'LP feasible', 2: 'LP INFEASIBLE'}.get(res.status, 'undecided')
      dead += res.status == 2; alive += res.status == 0; unk += res.status not in (0, 2)
      print(f'  {typ} x{cnt}: {st}  [{time.time() - t0:.0f} s]', flush=True)
  print(f'summary: LP-infeasible {dead}, LP-feasible {alive}, undecided {unk}')
