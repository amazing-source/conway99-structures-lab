import random, numpy as np
from general_model import build
from alldd import skeleton_from_pair_cycles, random_cycles
from alldd_near import rand_near
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp
rng = random.Random(11)
def solve(A, lo, hi, lb, ub, nvar, keep=None):
    if keep is not None:
        A = A[keep]; lo = lo[keep]; hi = hi[keep]
    r = milp(np.zeros(nvar), constraints=[LinearConstraint(A, lo, hi)], integrality=np.ones(nvar), bounds=Bounds(lb, ub), options={'time_limit': 60})
    return {0: 'alive', 2: 'dead'}.get(r.status, 'unk')
for trial in range(4):
    orbs = skeleton_from_pair_cycles(random_cycles(rng))
    near = rand_near(rng, orbs)
    res = {}
    A, lo, hi, lb, ub, nvar, phi, ns = build(orbs, dichotomy=False, near=near)
    res['nodich'] = solve(A, lo, hi, lb, ub, nvar)
    # drop symmetry rows: they are the rows with both +1 and -1 on phi of different orbits ... identify by construction order:
    # rebuild without balance too
    A2, lo2, hi2, lb2, ub2, nvar2, _, _ = build(orbs, dichotomy=False, balance=False, near=near)
    res['nodich_nobal'] = solve(A2, lo2, hi2, lb2, ub2, nvar2)
    # symmetry rows are the last C(21,2)=210 rows in the no-dichotomy build
    n = A2.shape[0]
    res['nodich_nobal_nosym'] = solve(A2, lo2, hi2, lb2, ub2, nvar2, keep=np.arange(n - 210))
    res['nodich_nosym'] = solve(A, lo, hi, lb, ub, nvar, keep=np.arange(A.shape[0] - 210))
    print(trial, res)
