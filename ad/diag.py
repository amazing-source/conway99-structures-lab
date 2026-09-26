import itertools, numpy as np
from pair_level import *
from scipy.optimize import milp, LinearConstraint, Bounds
C2 = r"C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/c2exp/logs/ad2_tail_survivors.txt"
reps = parse_survivors(C2)
dd = [(0, 1), (2, 3)]
ddm = {((0, 1), 0): ((2, 3), 0), ((2, 3), 0): ((0, 1), 0), ((0, 1), 1): ((2, 3), 1), ((2, 3), 1): ((0, 1), 1)}
R = [r for r in reps if r['rep'] == 130][0]; In = R['In']
base = {k: [tuple(In[k])] for k in range(N) if len(In[k]) in (2, 3)}
Ds = dict(base); Ds[6] = [(0, 1), (2, 3), (4, 5)]
A, lo, hi, lb, ub, nvar, v, nsidx, orbits = build(In, dd, ddm, Ds)
# rows: identify sym rows (last C(21,2)=210 rows)
nsym = 210
Aper, loper, hiper = A[:-nsym], lo[:-nsym], hi[:-nsym]
res = milp(c=np.zeros(nvar), constraints=[LinearConstraint(Aper, loper, hiper)], integrality=np.ones(nvar), bounds=Bounds(lb, ub))
print('without sym:', res.status, res.message)
# per orbit: which orbits individually infeasible (with ns free)
for o in orbits:
    idx = [orbits.index(o) * 21 + j for j in range(21)] + list(nsidx.values())
    rows = [i for i in range(len(Aper)) if np.any(Aper[i][[orbits.index(o) * 21 + j for j in range(21)]] != 0)]
    r = milp(c=np.zeros(nvar), constraints=[LinearConstraint(Aper[rows], loper[rows], hiper[rows])], integrality=np.ones(nvar), bounds=Bounds(lb, ub))
    if r.status != 0: print('orbit', o, 'infeasible alone')
