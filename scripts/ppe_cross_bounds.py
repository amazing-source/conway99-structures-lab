"""TEST (seconds; guides the hand proof of class b).  Under the protected rows, LP bounds on the number of doubles
received by the four cross pairs {1,3}, {1,4}, {2,3}, {2,4}, and which orbits are forced to double there."""
import sys

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack

sys.path.insert(0, 'scripts')
import ppe_test
from ppe_lp_relax import REPS
from ppe_certificate import LEMMA_A

w = 2
prot = set(LEMMA_A) | {"coupling P=[1, 2] Q=[5, 6]", "coupling P=[3, 4] Q=[5, 6]", "PPE Q=[3, 5]", "PPE Q=[4, 6]",
                       "PPE Q=[3, 6]", "PPE Q=[4, 5]"}
A, lo, hi, nv = ppe_test.build(REPS[w])
labels = list(ppe_test.LABELS)
mask = np.array([not l.startswith(('colsum', 'PPE', 'coupling')) or l in prot for l in labels])
Aa, la, ha = A[mask], lo[mask], hi[mask]
eq = la == ha
Aub = vstack([Aa[~eq], -Aa[~eq]])
bub = np.concatenate([ha[~eq], -la[~eq]])
idx = {key: ix for ix, key in ppe_test.VARS.items()}
cross = [frozenset(q) for q in ((1, 3), (1, 4), (2, 3), (2, 4))]


def opt(cvec):
    r = linprog(cvec, A_ub=Aub, b_ub=bub, A_eq=Aa[eq], b_eq=la[eq], bounds=(0, 1), method='highs')
    return r.fun if r.status == 0 else None


tot = np.zeros(nv)
for key, ix in idx.items():
    if key[0] == 't' and key[3] == 2 and key[2] in cross:
        tot[ix] = 1
print('LP min / max of the total doubles at the four cross pairs:', round(opt(tot), 3), round(-opt(-tot), 3))
for Q in cross:
    c = np.zeros(nv)
    for key, ix in idx.items():
        if key[0] == 't' and key[3] == 2 and key[2] == Q:
            c[ix] = 1
    print(f'  pair {sorted(Q)}: LP min {round(opt(c), 3)}, max {round(-opt(-c), 3)}')
forced = []
for key, ix in idx.items():
    if key[0] == 't' and key[3] == 2 and key[2] in cross:
        c = np.zeros(nv)
        c[ix] = 1
        m = opt(c)
        if m is not None and m > 1e-6:
            forced.append((key[1], sorted(key[2]), round(m, 3)))
print('orbits with a positive LP lower bound on doubling at a cross pair:', forced)

per_orbit = {}
for key, ix in idx.items():
    if key[0] == 't' and key[3] == 2 and key[2] in cross:
        per_orbit.setdefault(key[1], []).append(ix)
print('orbits with a positive LP lower bound on their number of doubles at the cross pairs:')
for X, ixs in sorted(per_orbit.items()):
    c = np.zeros(nv)
    c[ixs] = 1
    m = opt(c)
    if m is not None and m > 1e-6:
        print(f'   orbit {X}: >= {round(m, 3)}  (cross pairs disjoint from it: '
              f'{[sorted(ppe_test.VARS[i][2]) for i in ixs]})')
