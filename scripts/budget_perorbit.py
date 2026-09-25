"""For the two budget classes (STAR (3^3,4^3) class A and PASCH): per-orbit lower bound on doubles that does NOT depend
on the cyclic order: m(i,j) = min over admissible neighbour pairs of the local minimum.  If sum m >= 43 the budget
argument is order-free (easy to write by hand).  Also prints, per orbit, the local minimum for each admissible pair."""
import sys, itertools
sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V

REPS = {
    'STAR-A (3^3,4^3)': [(0, 1, 2), (0, 3, 1, 4), (0, 5, 1, 6), (2, 3, 4), (2, 5, 6), (3, 5, 4, 6)],
    'PASCH (3^3,4^3)': [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 6), (1, 5, 6), (2, 3, 4, 5), (2, 4, 6)],
}
for name, dec in REPS.items():
    cyc_of_edge = {}
    for c in dec:
        for t in range(len(c)):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % len(c)]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]
    def other_nb(i, j_):
        c = cyc_of_edge[frozenset((i, j_))]; m = len(c); t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]
        return b if a == j_ else a
    total_m = 0
    print(f'\n=== {name}: {dec}')
    for j in V:
        others = [r for r in V if r != j]
        blk = 0
        for i in others:
            k = other_nb(i, j)
            Rx = [r for r in V if r not in (i, j, k)]
            vals = {}
            for l, r in itertools.combinations([o for o in others if o != i], 2):
                v = d3.local_min(d3.canon(tris, i, j, k, Rx, frozenset((l, r))))
                if v is not None:
                    vals[(l, r)] = v
            m = min(vals.values()) if vals else None
            blk += m
            print(f'  block {j} orbit {i} (k={k}): order-free min doubles {m}; per pair {vals}')
        total_m += blk
        print(f'  block {j}: sum of order-free minima {blk}')
    print(f'  TOTAL order-free lower bound {total_m}  (budget 42)')
