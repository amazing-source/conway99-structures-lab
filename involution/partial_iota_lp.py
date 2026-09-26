"""How many iota-pairs of a (dead) near-design must be fixed before the LP of the full sigma-pair model ((i),(ii),(iii),
(iv),(v) all linear) becomes infeasible?  Fixing {e,f}: every variable (Z, e, g) or (Z, g, e) with g != f is set to 0."""
import json, sys, time, random
import numpy as np
from scipy.optimize import linprog
from fullblock_rigid import mates_from_D, build

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
cycles = [tuple(c) for c in cfg['cycles']]
m = mates_from_D(cycles)
var, reps, A, lo = build(m, cycles, {})
print('sigma-pair model: vars', len(var), 'rows', A.shape[0], flush=True)
t0 = time.time()
r = linprog(np.zeros(len(var)), A_eq=A, b_eq=lo, bounds=(0, 1), method='highs')
print('root LP:', {0: 'feasible', 2: 'INFEASIBLE'}.get(r.status, r.status), '%.1fs' % (time.time() - t0), flush=True)
rng = random.Random(0)
for k_iota in range(min(2, len(iotas))):
    pairs = [(tuple(e), tuple(f)) for e, f in iotas[k_iota]]
    order = list(range(30)); rng.shuffle(order)
    for kfix in (3, 6, 10, 15, 20, 30):
        fixed = [pairs[i] for i in order[:kfix]]
        partner = {}
        for e, f in fixed:
            partner[e] = f; partner[f] = e
        ub = np.ones(len(var))
        for j, (Zi, e, e2) in enumerate(var):
            if (e in partner and partner[e] != e2) or (e2 in partner and partner[e2] != e):
                ub[j] = 0
        t0 = time.time()
        r = linprog(np.zeros(len(var)), A_eq=A, b_eq=lo, bounds=list(zip(np.zeros(len(var)), ub)), method='highs')
        print('iota #%d, %2d pairs fixed: LP' % (k_iota, kfix), {0: 'feasible', 2: 'INFEASIBLE'}.get(r.status, r.status),
              '%.1fs' % (time.time() - t0), flush=True)
