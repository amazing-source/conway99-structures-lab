"""C2 §17.3 / §18.3 in pure twisted configurations: each orbit X (pair ij, paired via i to {i,k}, free at j) gets a zero
at rs ⊂ R_x = [7] - {i,j,k} for each R_j-neighbour lying in a pair jr of a D-triangle {j,r,s} with r, s in R_x
(zeros as a SET: two neighbours giving the same rs count once).  doubles(x) >= omega(x) + zeros(x) >= zeros(x), and
16.1 gives sum over the 84 labels of doubles = 84 (n_s = 0 for 'un bloc' types), i.e. sum over the 42 orbits <= 42.
The R_j-neighbours of X are its two neighbours on the free 6-cycle at j, whose cyclic order is unknown: we MINIMISE the
zero count over all 60 cyclic orders, independently at each block (lower bound).  A type is killed if min > 42.
Only types with >= 4 triangles can exceed 42 (each block j gives <= 4 * #triangles through j)."""
import itertools, sys
sys.path.insert(0, 'scripts')
from k7_cycle_decomp import EDGES, EID, eid, cycles_through

V = list(range(7))
def decompositions(mask, acc):
    if mask == 0:
        yield list(acc); return
    e0 = (mask & -mask).bit_length() - 1
    for es, cyc in cycles_through(mask, e0):
        sub = mask
        for e in es: sub &= ~(1 << e)
        acc.append(cyc)
        yield from decompositions(sub, acc)
        acc.pop()

CYC6 = []
for perm in itertools.permutations(range(1, 6)):
    if perm[0] < perm[-1]:
        CYC6.append((0,) + perm)            # 60 cyclic orders of 6 positions

def block_min_zeros(dec, j):
    others = [r for r in V if r != j]
    # neighbours of each vertex in each cycle, and which cycle contains each edge
    cyc_of_edge = {}
    for c in dec:
        n = len(c)
        for t in range(n):
            a, b = c[t], c[(t + 1) % n]
            cyc_of_edge[frozenset((a, b))] = (c, t)
    def other_nb(i, j_):     # neighbour of i, in the cycle containing edge ij_, other than j_
        c, _ = cyc_of_edge[frozenset((i, j_))]
        n = len(c); t = c.index(i)
        a, b = c[(t - 1) % n], c[(t + 1) % n]
        return b if a == j_ else a
    def tri(r):              # third vertex if edge jr lies on a triangle
        c, _ = cyc_of_edge[frozenset((j, r))]
        return next(x for x in c if x not in (j, r)) if len(c) == 3 else None
    k_of = {i: other_nb(i, j) for i in others}
    zp = {}
    for i in others:
        R = set(V) - {i, j, k_of[i]}
        for r in others:
            s = tri(r)
            zp[(i, r)] = frozenset((r, s)) if (s is not None and r in R and s in R) else None
    # C2 §17 (rappels): if x's own cycle is a triangle {i,j,k}, phi_x(jk) = 0, so x is not R_j-adjacent to the pair jk
    forbid = {frozenset((i, k_of[i])) for i in others if tri(i) is not None}
    best = 99
    for order in CYC6:
        lab = [others[p] for p in order]
        if any(frozenset((lab[t], lab[(t + 1) % 6])) in forbid for t in range(6)):
            continue
        # exact count of a complete D-triangle {j,r,s} with i not in it: x has exactly one neighbour among its labels,
        # so x's two R_j-neighbours cannot both lie in pairs of the same such triangle
        bad = False
        for t in range(6):
            i, l, rr = lab[t], lab[t - 1], lab[(t + 1) % 6]
            if tri(l) is not None and tri(l) == rr and i not in (l, rr):
                bad = True; break
        if bad:
            continue
        tot = 0
        for t in range(6):
            i, l, rr = lab[t], lab[t - 1], lab[(t + 1) % 6]
            tot += len({z for z in (zp[(i, l)], zp[(i, rr)]) if z is not None})
        best = min(best, tot)
    return best

if __name__ == '__main__':
    results = {}
    for dec in decompositions((1 << 21) - 1, []):
        typ = tuple(sorted(len(c) for c in dec))
        if sum(1 for c in dec if len(c) == 3) < 4:
            continue
        z = sum(block_min_zeros(dec, j) for j in V)
        results.setdefault(typ, []).append(z)
    for typ, zs in sorted(results.items()):
        print(f'{typ}: {len(zs)} labelled decompositions; min over them of (min total zeros per orbit-sum) = {min(zs)}, '
              f'max = {max(zs)}; killed (all > 42): {min(zs) > 42}; count with >42: {sum(z > 42 for z in zs)}')
