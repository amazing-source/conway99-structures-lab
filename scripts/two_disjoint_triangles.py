"""Check (a check only; the proof is by hand) of the corollary of the Complement Lemma: in the pure twisted case, a
decomposition with two vertex-disjoint triangles T, T' is impossible.
Hand proof: let i be the 7th block.  The 3 cycles through i pair up its 6 neighbours; |T| = 3 is odd, so some cycle
passes s - i - t with s in T, t in T'.  At block j = s, the orbit of {i,s} paired via i has k = t and
R_x = (T - s) u (T' - t): C2 on T gives phi(ab) = 1 - N_ab, CL on T' gives phi(ab) = 3 - N_ab, where {a,b} = T - s.
Contradiction: the orbit is dead.
Checks: (1) the dead key is infeasible for every neighbour pair, with or without an own triangle; (2) every labelled
decomposition with exactly two triangles, disjoint, has such an x - i - y; their counts per type equal the counts
'impossible outright' of the exhaustive run (logs/doubles_budget3.log), so these are exactly the ones it killed."""
import itertools
import sys
from collections import Counter

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions, V
from local_system import solve_core, frame

# (1) the dead key, in the frame i = 0, j = 1, k = 2, R = {3,4,5,6}
for ab in itertools.combinations((3, 4, 5, 6), 2):
    cd = tuple(r for r in (3, 4, 5, 6) if r not in ab)
    for own in (False, True):
        for pair in itertools.combinations(range(2, 7), 2):
            assert not solve_core([{1, *ab}, {2, *cd}], own, 0, 1, 2, [3, 4, 5, 6], set(pair))
print('(1) triangles {j,a,b} and {k,c,d} covering R_x: no solution for any neighbour pair (with or without own triangle)')

_memo = {}


def dead(tris, own, i, j, k, R, pair):
    """solve_core, memoised on the key relabelled into the frame i=0, j=1, k=2, R -> 3..6 (sorted)."""
    mp = {i: 0, j: 1, k: 2}
    mp.update({r: 3 + t for t, r in enumerate(sorted(R))})
    key = (tuple(sorted(tuple(sorted(mp[v] for v in t)) for t in tris)), own, tuple(sorted(mp[v] for v in pair)))
    if key not in _memo:
        _memo[key] = not solve_core([set(t) for t in key[0]], own, 0, 1, 2, [3, 4, 5, 6], set(key[2]))
    return _memo[key]


# (2) all decompositions with exactly two triangles
count = Counter()
for dec in decompositions((1 << 21) - 1, []):
    tris = [set(c) for c in dec if len(c) == 3]
    if len(tris) != 2:
        continue
    typ = tuple(sorted(len(c) for c in dec))
    T, T2 = tris
    if T & T2:
        count[(typ, 'sharing')] += 1
        continue
    (i,) = set(V) - T - T2
    mixed = []
    for c in dec:
        if i in c:
            m = len(c)
            t = c.index(i)
            u, w = c[(t - 1) % m], c[(t + 1) % m]
            if (u in T) != (w in T):
                mixed.append((u, w) if u in T else (w, u))
    assert mixed, dec
    for x, y in mixed:
        for jj, kk in ((x, y), (y, x)):
            cyc, k, R = frame(dec, jj, i)
            assert k == kk
            others = [set(c) for c in dec if len(c) == 3 and set(c) != set(cyc)]
            for pair in itertools.combinations([r for r in V if r not in (i, jj)], 2):
                assert dead(others, len(cyc) == 3, i, jj, k, R, pair)
    count[(typ, 'disjoint')] += 1
for key, v in sorted(count.items()):
    print(f'(2) type {key[0]}, triangles {key[1]}: {v} labelled' + ('  -> dead orbit found and checked' if key[1] ==
                                                                  'disjoint' else ''))
logged = {(3, 3, 4, 4, 7): 0, (3, 3, 4, 5, 6): 5040, (3, 3, 5, 5, 5): 840}
for typ, n in logged.items():
    assert count[(typ, 'disjoint')] == n, (typ, count[(typ, 'disjoint')], n)
print(f'    ({len(_memo)} distinct frame keys solved)')
print('    disjoint counts = impossible-outright counts of the exhaustive run for every 2-triangle type')
