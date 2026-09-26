"""All-dd + near structures that pass every SINGLE-ORBIT check (per-orbit subsystem of the pair model feasible), found by
greedy repair; then the coupled pair model.  Is the coupling what kills?"""
import random, sys, time, numpy as np
from collections import Counter
from scipy.optimize import milp, LinearConstraint, Bounds
from general_model import build, feasible, N, PAIRS
from alldd import skeleton_from_pair_cycles, random_cycles
from alldd_near import rand_near
def orbit_ok(A, lo, hi, lb, ub, nvar, phi, i):
    cols = [phi[(i, P)] for P in PAIRS]
    A = A.tocsr()
    rows = sorted(set(A[:, cols].nonzero()[0]))
    # keep only rows that involve ONLY this orbit's phi (and n_s / indicator vars): drop symmetry rows (they involve other orbits)
    keep = []
    colset = set(cols)
    for r in rows:
        nz = A.indices[A.indptr[r]:A.indptr[r + 1]]
        if all((j in colset) or (j >= 42 * 21) for j in nz): keep.append(r)
    res = milp(np.zeros(nvar), constraints=[LinearConstraint(A[keep], lo[keep], hi[keep])], integrality=np.ones(nvar),
               bounds=Bounds(lb, ub), options={'time_limit': 10})
    return res.status == 0
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
stats = Counter(); t0 = time.time(); attempts = 0
target = int(sys.argv[2]) if len(sys.argv) > 2 else 4
while sum(stats.values()) < target and attempts < 400:
    attempts += 1
    orbs = skeleton_from_pair_cycles(random_cycles(rng))
    near = rand_near(rng, orbs)
    if near is None: continue
    A, lo, hi, lb, ub, nvar, phi, ns = build(orbs, dichotomy=False, near=near)
    if np.any(lb > ub): continue
    if not all(orbit_ok(A, lo, hi, lb, ub, nvar, phi, i) for i in range(42)): continue
    r = feasible(orbs, near=near, time_limit=120)
    stats['dead' if r is None else ('unknown' if isinstance(r, str) else 'ALIVE')] += 1
    print(stats, 'after', attempts, 'attempts', flush=True)
print('done', dict(stats), 'attempts', attempts, '%.0fs' % (time.time() - t0))
