"""Check of the hand proof for STAR (3^3,6,6), block 0 (representative
[(0,1,2),(0,3,1,4,5,6),(0,4,3,6,1,5),(2,3,5),(2,4,6)]):
 (1) every pair declared forbidden by facts F_a..F_e is infeasible for the local system (local_min is None);
 (2) with ONLY these restrictions (all other pairs allowed), no 6-cycle order exists."""
import sys, itertools
sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V, CYC6
dec = [(0, 1, 2), (0, 3, 1, 4, 5, 6), (0, 4, 3, 6, 1, 5), (2, 3, 5), (2, 4, 6)]
tris = [tuple(sorted(c)) for c in dec if len(c) == 3]
j = 0
k_of = {1: 2, 2: 1, 3: 1, 4: 3, 5: 1, 6: 5}
others = [1, 2, 3, 4, 5, 6]
P = lambda a, b: frozenset((a, b))
forbidden = {
    1: {P(2, x) for x in (3, 4, 5, 6)} | {P(3, 5), P(4, 6)},          # F_a
    2: {P(1, x) for x in (3, 4, 5, 6)},                                 # consequence of F_a (1 not adjacent to 2)
    4: {P(2, 1), P(2, 3), P(2, 5)},                                     # F_b
    6: {P(2, 1), P(2, 3), P(2, 5)},                                     # F_c
    3: {P(1, 5), P(1, 2)},                                              # F_d
    5: {P(1, 3), P(1, 2)},                                              # F_e
}
ok1 = True
for i, pairs in forbidden.items():
    k = k_of[i]
    Rx = [r for r in V if r not in (i, j, k)]
    for pr in pairs:
        if i in pr:
            continue
        v = d3.local_min(d3.canon(tris, i, j, k, Rx, pr))
        if i == 2:          # forbidden for 2 only as a consequence of F_a, not by its own local system
            continue
        if v is not None:
            ok1 = False
            print('NOT infeasible:', i, sorted(pr), v)
print('(1) all declared-forbidden pairs infeasible:', ok1)
cnt = 0
for order in CYC6:
    lab = [others[p] for p in order]
    good = True
    for t in range(6):
        i, l, r = lab[t], lab[t - 1], lab[(t + 1) % 6]
        if P(l, r) in forbidden.get(i, set()):
            good = False; break
    cnt += good
print('(2) 6-cycles compatible with F_a..F_e only:', cnt)
