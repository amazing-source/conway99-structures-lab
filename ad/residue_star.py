"""AD residues (dd partners as variables, A-structure fixed) + (*): does (*) kill them?"""
import json, re, sys, time, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from pair_level import N, pr, two_regular_graphs, PAIRS
import pair_level_relaxed as PLR
from star_generic import augment
reps = json.load(open(sys.argv[1]))
logf = sys.argv[2]
lines = [l for l in open(logf) if re.match(r'rep \d+ n_dd \d+ free', l)]
for line in lines:
    m = re.match(r"rep (\d+) n_dd (\d+) free blocks \[(\d+)\] structures \d+ alive (\d+) unknown (\d+)", line)
    t, nd, A, al, un = map(int, m.groups())
    if al == 0 and un == 0: continue
    rep = reps[t]; dd = [tuple(p) for p in rep['dd']]
    head = {pr(a, h): h for a, h in rep['arcs']}
    In = {c: [x for x in range(N) if x != c and head.get(pr(c, x)) == c] for c in range(N)}
    alive = []; t0 = time.time()
    for cyc in two_regular_graphs(In[A]):
        pre = PLR.feasible_relaxed(In, dd, dichotomy=True, fixedDs={A: cyc}, time_limit=60)
        if pre is None: continue
        out = PLR.build_relaxed(In, dd, dichotomy=True, fixedDs={A: cyc})
        if out is None or isinstance(out, str): continue
        Am, lo, hi, lb, ub, nvar = out
        # phi indices in build_relaxed: phi = {(o,P): idx} created first in order orbits x PAIRS
        orbits = [(p, s) for p in PAIRS for s in (0, 1)]
        phi = {}; idx = 0
        for o in orbits:
            for P in PAIRS: phi[(o, P)] = idx; idx += 1
        A2, lo2, hi2, lb2, ub2, nv2 = augment(Am, lo, hi, lb, ub, nvar, phi, {o: o[0] for o in orbits})
        integ = np.concatenate([np.ones(nvar), np.zeros(nv2 - nvar)])
        res = milp(np.zeros(nv2), constraints=[LinearConstraint(A2, lo2, hi2)], integrality=integ, bounds=Bounds(lb2, ub2),
                   options={'time_limit': 120})
        st = {0: 'alive', 2: 'dead'}.get(res.status, 'unknown')
        if st != 'dead': alive.append((cyc, st))
    print('rep', t, 'n_dd', nd, 'with (*): A-structures not dead', len(alive), alive[:3], '%.0fs' % (time.time() - t0), flush=True)
