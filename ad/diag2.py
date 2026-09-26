"""For each AD2 survivor: which single orbits W (not through the A-block) are infeasible ALONE, using at the A-block only
the union balance (sum_m phi(A m) = 4, i.e. no knowledge of the D-structure at A)?"""
import itertools, numpy as np
from pair_level import *
from scipy.optimize import milp, LinearConstraint, Bounds
C2 = r"C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/c2exp/logs/ad2_tail_survivors.txt"
reps = parse_survivors(C2)
dd = [(0, 1), (2, 3)]
ddm = {((0, 1), 0): ((2, 3), 0), ((2, 3), 0): ((0, 1), 0), ((0, 1), 1): ((2, 3), 1), ((2, 3), 1): ((0, 1), 1)}
for R in reps:
    In = R['In']
    base = {k: [tuple(In[k])] for k in range(N) if len(In[k]) in (2, 3)}
    Ds = dict(base); Ds[6] = [(0, 1, 2, 3, 4, 5)]      # a 6-cycle: at k = 6 its only component is the union
    A, lo, hi, lb, ub, nvar, v, nsidx, orbits = build(In, dd, ddm, Ds)
    Aper, loper, hiper = A[:-210], lo[:-210], hi[:-210]
    bad = []
    for o in orbits:
        if 6 in o[0]: continue
        cols = [orbits.index(o) * 21 + j for j in range(21)]
        rows = [i for i in range(len(Aper)) if np.any(Aper[i][cols] != 0)]
        r = milp(c=np.zeros(nvar), constraints=[LinearConstraint(Aper[rows], loper[rows], hiper[rows])], integrality=np.ones(nvar), bounds=Bounds(lb, ub))
        if r.status != 0: bad.append(o)
    print('rep', R['rep'], 'In', In, '\n   infeasible alone (A-free):', bad)
