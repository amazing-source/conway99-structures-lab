import sys, time
sys.path.insert(0, 'scripts')
from ppe_test import feasible
reps = {
 '(3,3,4,4,7)': [(1, 2, 0), (1, 3, 2, 5, 0, 6, 4), (1, 5, 3, 6), (2, 6, 5, 4), (0, 3, 4)],
 '(3,3,4,5,6)a': [(3, 4, 0), (3, 1, 4, 6, 0, 5), (3, 6, 1, 5, 2), (4, 5, 6, 2), (0, 1, 2)],
 '(3,3,4,5,6)f': [(1, 0, 2), (1, 5, 0, 6), (1, 3, 5, 2, 6, 4), (0, 3, 4), (2, 3, 6, 5, 4)],
 '(3,3,5,5,5)': [(3, 4, 0), (3, 5, 4, 6, 2), (3, 6, 0, 5, 1), (4, 2, 5, 6, 1), (0, 2, 1)],
}
t0 = time.time()
for name, dec in reps.items():
    r = {}
    for label, kw in (('PPE doubles only', dict(ppe=True, colsum=False, coupling=False)),
                      ('column sums only', dict(ppe=False, colsum=True, coupling=False)),
                      ('coupling only', dict(ppe=False, colsum=False, coupling=True)),
                      ('PPE + colsum', dict(ppe=True, colsum=True, coupling=False)),
                      ('PPE + coupling', dict(ppe=True, colsum=False, coupling=True))):
        r[label] = feasible(dec, **kw)
    print(name, r, f'[{time.time()-t0:.0f} s]', flush=True)
