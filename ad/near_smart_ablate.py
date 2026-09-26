import random, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from general_model import build, N, PAIRS
from alldd import skeleton_from_pair_cycles, random_cycles
from alldd_near import rand_near
from near_smart import orbit_ok
def solve(A, lo, hi, lb, ub, nvar, keep=None):
    A = A.tocsr()
    if keep is not None: A, lo, hi = A[keep], lo[keep], hi[keep]
    r = milp(np.zeros(nvar), constraints=[LinearConstraint(A, lo, hi)], integrality=np.ones(nvar), bounds=Bounds(lb, ub), options={'time_limit': 60})
    return {0: 'alive', 2: 'dead'}.get(r.status, 'unk')
rng = random.Random(3); done = 0
while done < 4:
    orbs = skeleton_from_pair_cycles(random_cycles(rng)); near = rand_near(rng, orbs)
    if near is None: continue
    A, lo, hi, lb, ub, nvar, phi, ns = build(orbs, dichotomy=False, near=near)
    if np.any(lb > ub) or not all(orbit_ok(A, lo, hi, lb, ub, nvar, phi, i) for i in range(42)): continue
    done += 1
    nosym = np.arange(A.shape[0] - 210)
    A2, lo2, hi2, lb2, ub2, nvar2, phi2, ns2 = build(orbs, dichotomy=True, near=near)
    print({'sym_nodich': solve(A, lo, hi, lb, ub, nvar), 'nosym_nodich': solve(A, lo, hi, lb, ub, nvar, nosym),
           'sym_dich': solve(A2, lo2, hi2, lb2, ub2, nvar2)})
