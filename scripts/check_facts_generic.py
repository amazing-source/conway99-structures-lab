"""Generic check of a hand proof at one block: (1) each declared-forbidden (orbit, pair) is infeasible for the local
system, unless marked as a 'consequence' (derived from another orbit's fact); (2) no 6-cycle survives the declared
restrictions alone."""
import sys
sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V, CYC6

def other_nb(dec, i, j):
    for c in dec:
        m = len(c)
        for t in range(m):
            if {c[t], c[(t + 1) % m]} == {i, j}:
                a, b = c[t - 1] if False else None, None
                pos = c.index(i)
                u, w = c[(pos - 1) % m], c[(pos + 1) % m]
                return w if u == j else u
    raise ValueError

def check(name, dec, j, forbidden, consequence=()):
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]
    others = [r for r in V if r != j]
    ok1 = True
    for i, pairs in forbidden.items():
        if i in consequence:
            continue
        k = other_nb(dec, i, j)
        Rx = [r for r in V if r not in (i, j, k)]
        for pr in pairs:
            v = d3.local_min(d3.canon(tris, i, j, k, Rx, frozenset(pr)))
            if v is not None:
                ok1 = False
                print('   NOT infeasible:', i, sorted(pr), v)
    cnt = 0
    for order in CYC6:
        lab = [others[p] for p in order]
        if all(frozenset((lab[t - 1], lab[(t + 1) % 6])) not in forbidden.get(lab[t], set()) for t in range(6)):
            cnt += 1
    print(f'{name}: (1) declared-forbidden pairs infeasible: {ok1}; (2) surviving 6-cycles: {cnt}')

F = lambda *ps: {frozenset(p) for p in ps}
# STAR class B, block 6
starB = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 1, 6), (2, 3, 5), (2, 4, 6), (3, 4, 5, 6)]
check('STAR-B block 6', starB, 6, {
    4: F((0, 1), (3, 5), (0, 2), (1, 2), (2, 3), (2, 5)),     # ML1: one in {0,1}, one in {3,5}; never 2 (k)
    0: F((2, 3), (2, 4), (2, 5)),
    1: F((2, 3), (2, 4), (2, 5)),
    3: F((2, 4), (4, 5)),
    5: F((2, 4), (3, 4)),
    2: F((2, 4)) and F((4, 0), (4, 1), (4, 3), (4, 5)),        # 2 never adjacent to its k = 4
})
# QUAD, block 1
quad = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 4, 6), (1, 5, 6), (2, 3, 6), (2, 4, 5)]
check('QUAD block 1', quad, 1, {
    0: F((3, 6), (4, 5), (5, 6), (2, 3), (2, 4), (2, 5), (2, 6)),   # ML1 + C2 ; C0: never 2
    2: F((5, 6), (0, 3), (0, 4), (0, 5), (0, 6)),                    # C2 ; C0: never 0
    5: F((0, 2), (2, 3), (6, 0), (6, 2), (6, 3), (6, 4)),            # C2 ; C8 ; C0: never 6
    6: F((0, 2), (2, 4), (5, 0), (5, 2), (5, 3), (5, 4)),            # C2 ; C8 ; C0: never 5
    3: F((0, 6),),
    4: F((0, 5),),
})
