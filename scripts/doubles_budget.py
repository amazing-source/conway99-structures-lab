"""Stronger version of zeros_budget.py: lower bound on doubles(x) from ALL exact triangle counts (C2 §16).
x in pair ij, paired via i to {i,k} (k = other neighbour of i on the cycle through ij), free at j; R_x = [7] - {i,j,k}.
Unknowns phi = phi_x(P) in {0,1,2} for the 6 pairs P in R_x, with sum = 6 + omega, omega >= 0 (C2 18.3).
Exact counts |N(x) ∩ Lab(tau)| for a complete D-triangle tau: 2 if x's pair in tau, 1 if it meets tau in one block,
3 if disjoint.  Facts used: x's only neighbours meeting block i are v, sigma v (in pair ik), so phi_x(ia) = 0 for a != k;
x's neighbours in pairs ja are exactly its R_j-neighbours there (one or none).
  (T0) tau inside R_x:                   phi over its 3 pairs sums to 3;
  (Ti) tau = {i,a,b}, j not in, a,b in R_x: phi(ab) = 1;
  (Tj) tau = {j,a,b}, i not in, a,b in R_x: phi(ab) = 1 - #(R_j-neighbours in ja, jb)  (that number must be <= 1);
  (Tj') tau = {j,k,b} (i not in): phi(kb) is outside R_x, but #R_j-nbrs in {jk, jb} <= 1;
  (a) own triangle {i,j,k}: no R_j-neighbour in jk.
doubles(x) >= min #{phi = 2} over feasible phi.  Budget (16.1, n_s = 0): sum over the 42 orbits of doubles <= 42.
Minimised over the free 6-cycle order at each block (independently); an order is discarded if some x is infeasible."""
import itertools, sys
sys.path.insert(0, 'scripts')
from zeros_budget import decompositions, CYC6, V

PAIRS_OF = {}
def pairs4(R):
    return [frozenset(p) for p in itertools.combinations(sorted(R), 2)]

def min_doubles(i, j, k, tris, nbrs):
    """tris: list of triangles (frozensets); nbrs: set of blocks r with an R_j-neighbour of x in pair jr."""
    R = set(V) - {i, j, k}
    P = pairs4(R)
    fixed, sums3 = {}, []
    for t in tris:
        if t <= R:
            sums3.append([frozenset(e) for e in itertools.combinations(sorted(t), 2)])
        elif i in t and j not in t:
            a, b = sorted(t - {i})
            if a in R and b in R:
                if fixed.get(frozenset((a, b)), 1) != 1: return None
                fixed[frozenset((a, b))] = 1
        elif j in t and i not in t:
            a, b = sorted(t - {j})
            n = (a in nbrs) + (b in nbrs)
            if n > 1: return None
            if a in R and b in R:
                val = 1 - n
                if fixed.get(frozenset((a, b)), val) != val: return None
                fixed[frozenset((a, b))] = val
    best = None
    free = [p for p in P if p not in fixed]
    for vals in itertools.product((0, 1, 2), repeat=len(free)):
        phi = dict(fixed); phi.update(zip(free, vals))
        if sum(phi.values()) < 6: continue
        if any(sum(phi[e] for e in s) != 3 for s in sums3): continue
        d = sum(1 for v in phi.values() if v == 2)
        if best is None or d < best: best = d
    return best

def analyse(dec):
    cyc_of_edge = {}
    for c in dec:
        n = len(c)
        for t in range(n):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % n]))] = c
    tris = [frozenset(c) for c in dec if len(c) == 3]
    def other_nb(i, j_):
        c = cyc_of_edge[frozenset((i, j_))]; n = len(c); t = c.index(i)
        a, b = c[(t - 1) % n], c[(t + 1) % n]
        return b if a == j_ else a
    total = 0
    for j in V:
        others = [r for r in V if r != j]
        k_of = {i: other_nb(i, j) for i in others}
        own_tri = {i: (frozenset((i, j, k_of[i])) in tris) for i in others}
        best = None
        for order in CYC6:
            lab = [others[p] for p in order]
            tot = 0
            for t in range(6):
                i, l, rr = lab[t], lab[t - 1], lab[(t + 1) % 6]
                if own_tri[i] and k_of[i] in (l, rr):
                    tot = None; break
                md = min_doubles(i, j, k_of[i], tris, {l, rr})
                if md is None:
                    tot = None; break
                tot += md
            if tot is not None and (best is None or tot < best):
                best = tot
        if best is None:
            return None          # no admissible order at block j: configuration impossible outright
        total += best
    return total

if __name__ == '__main__':
    mint = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    res = {}
    for dec in decompositions((1 << 21) - 1, []):
        nt = sum(1 for c in dec if len(c) == 3)
        if nt < mint: continue
        typ = tuple(sorted(len(c) for c in dec))
        res.setdefault(typ, []).append(analyse(dec))
    for typ, zs in sorted(res.items(), key=lambda kv: -kv[0].count(3)):
        imp = sum(z is None for z in zs); vals = [z for z in zs if z is not None]
        print(f'{typ}: {len(zs)} labelled; impossible outright {imp}; min doubles {min(vals) if vals else "-"}; '
              f'max {max(vals) if vals else "-"}; killed (all None or > 42): {all(z is None or z > 42 for z in zs)}')
