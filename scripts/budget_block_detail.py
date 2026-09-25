"""Detail for a budget class at one block: admissible pairs with their local minimum of doubles, the compatible cyclic
orders and their doubles sums."""
import sys, itertools
sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V, CYC6

def detail(name, dec, j):
    cyc_of_edge = {}
    for c in dec:
        for t in range(len(c)):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % len(c)]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]
    def other_nb(i, j_):
        c = cyc_of_edge[frozenset((i, j_))]; m = len(c); t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]
        return b if a == j_ else a
    others = [r for r in V if r != j]
    val = {}
    print(f'=== {name}, block {j}')
    for i in others:
        k = other_nb(i, j)
        Rx = [r for r in V if r not in (i, j, k)]
        adm = {}
        for l, r in itertools.combinations([o for o in others if o != i], 2):
            v = d3.local_min(d3.canon(tris, i, j, k, Rx, frozenset((l, r))))
            if v is not None:
                adm[(l, r)] = v
                val[(i, frozenset((l, r)))] = v
        print(f'  orbit {i} (k={k}, R={Rx}): admissible pairs with min doubles {adm}')
    res = []
    for order in CYC6:
        lab = [others[p] for p in order]
        tot = 0
        for t in range(6):
            key = (lab[t], frozenset((lab[t - 1], lab[(t + 1) % 6])))
            if key not in val:
                tot = None; break
            tot += val[key]
        if tot is not None:
            res.append((tot, lab))
    for tot, lab in sorted(res):
        print(f'  compatible order {lab}: doubles {tot}')

if __name__ == '__main__':
    starA = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 1, 6), (2, 3, 4), (2, 5, 6), (3, 5, 4, 6)]
    detail('STAR-A', starA, 6)
    detail('STAR-A', starA, 2)
