"""Why does a random near structure kill all-dd at pair level?  Per-orbit feasibility (drop symmetry + dichotomy)."""
import random, numpy as np
from general_model import build, PAIRS
from alldd import skeleton_from_pair_cycles, random_cycles
from alldd_near import rand_near
from scipy.optimize import milp, LinearConstraint, Bounds
rng = random.Random(5)
for trial in range(3):
    orbs = skeleton_from_pair_cycles(random_cycles(rng))
    near = rand_near(rng, orbs)
    res = {}
    for dich in (False, True):
        A, lo, hi, lb, ub, nvar, phi, ns = build(orbs, dichotomy=dich, near=near)
        r = milp(np.zeros(nvar), constraints=[LinearConstraint(A, lo, hi)], integrality=np.ones(nvar), bounds=Bounds(lb, ub), options={'time_limit': 60})
        res['dich' if dich else 'nodich'] = r.status
    # no balance
    A, lo, hi, lb, ub, nvar, phi, ns = build(orbs, dichotomy=False, balance=False, near=near)
    r = milp(np.zeros(nvar), constraints=[LinearConstraint(A, lo, hi)], integrality=np.ones(nvar), bounds=Bounds(lb, ub), options={'time_limit': 60})
    res['nobal_nodich'] = r.status
    print(trial, res)
