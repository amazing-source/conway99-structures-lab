"""Size of the configuration model: number of local solutions per (orbit, neighbour pair) on a one-triangle sample."""
import itertools, sys, time
sys.path.insert(0, 'scripts')
from zeros_budget import V
from local_system import frame, solve_core
from ppe_typed_lp import SAMPLES

t0 = time.time()
for dec in SAMPLES[:1]:
    tot, per = 0, []
    for j in V:
        for i in V:
            if i == j:
                continue
            cyc, k, R = frame(dec, j, i)
            tris = [set(c) for c in dec if len(c) == 3 and set(c) != set(cyc)]
            c = 0
            for pair in itertools.combinations([r for r in V if r not in (i, j)], 2):
                c += len(solve_core(tris, len(cyc) == 3, i, j, k, R, set(pair)))
            per.append(c); tot += c
    print(dec, 'total (orbit, pair, solution) triples:', tot, ' per orbit min/max:', min(per), max(per),
          f'[{time.time() - t0:.0f} s]')
