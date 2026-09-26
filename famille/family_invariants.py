"""Invariants of C1/C2 (lines, holonomy, prisms, Hamming points, p-ranks) computed on the known members m = 2, 11,
with the family formulas in m. Generic code for any srg(v, 2m, 1, 2) given by adjacency sets."""
import itertools, sys
from family_graphs import hamming_9, bvls_243

def lines_of(adj):
    L = set()
    for a in range(len(adj)):
        for b in adj[a]:
            if b > a:
                c = (adj[a] & adj[b]).pop()
                L.add(frozenset((a, b, c)))
    return [tuple(sorted(l)) for l in L]

def holonomy(adj, T):
    """sigma_T on the layer of T[0]: layer_j = N(t_j) - T ; cross matching layer_j -> layer_j' : p -> the unique
    neighbour of p in layer_j' (second common neighbour of p and t_j')."""
    t = T
    layer = [sorted(adj[x] - set(T)) for x in t]
    def cross(p, j2):
        c = [q for q in adj[p] if q in adj[t[j2]] and q not in T]
        assert len(c) == 1
        return c[0]
    sigma = {}
    for p in layer[0]:
        q = cross(p, 1); r = cross(q, 2); s = cross(r, 0)
        sigma[p] = s
    return sigma

def analyse(name, adj, m):
    n = len(adj)
    L = lines_of(adj)
    fixes = []
    for T in L:
        s = holonomy(adj, T)
        fixes.append(sum(1 for p in s if s[p] == p))
    P = sum(fixes) // 2
    Pmax = len(L) * (m - 1)
    # m(e): number of prisms in which the edge e is a rung = number of lines through each end ... computed directly:
    # a prism = lines T, T' with a perfect matching; rung (a, a') with a in T, a' in T'.
    rung = {}
    for T in L:
        s = holonomy(adj, T)
        for p in s:
            if s[p] == p:                                     # prism on T through label p: T' = {p, ...}
                pass
    # Hamming points: every pair of lines through x spans a complete 3x3 grid
    def is_hamming(x):
        lx = [T for T in L if x in T]
        for T1, T2 in itertools.combinations(lx, 2):
            a = [y for y in T1 if y != x]; b = [y for y in T2 if y != x]
            pts = {x, *a, *b}
            for ai in a:
                for bj in b:
                    c = (adj[ai] & adj[bj]) - {x}
                    assert len(c) == 1
                    pts |= c
            if len(pts) != 9:
                return False
            # the 9 points must induce K3 x K3 (degree 4 inside)
            if any(len(adj[p] & pts) != 4 for p in pts):
                return False
        return True
    ham = sum(1 for x in range(n) if is_hamming(x))
    print(f'{name}: v = {n}, m = {m}, lines = {len(L)} (formula v*m/3 = {n * m // 3}); holonomy fixed points per line: '
          f'min {min(fixes)}, max {max(fixes)} (untwisted value 2m-2 = {2 * m - 2}); P = {P} (P_max = L(m-1) = {Pmax}); '
          f'Hamming points: {ham}/{n}')
    return L

def rank_mod(M, p):
    M = [[x % p for x in row] for row in M]
    r = 0
    ncol = len(M[0])
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r

def smith_p_part(M, p):
    """p-adic valuations of the elementary divisors (Smith normal form over Z_(p)) of an integer matrix."""
    from fractions import Fraction
    M = [row[:] for row in M]
    n, mcols = len(M), len(M[0])
    def val(x):
        if x == 0:
            return None
        v = 0
        while x % p == 0:
            x //= p; v += 1
        return v
    vals = []
    r = 0
    rows = list(range(n)); cols = list(range(mcols))
    # work with Fractions whose denominators are prime to p
    F = [[Fraction(x) for x in row] for row in M]
    def fval(f):
        if f == 0:
            return None
        return val(f.numerator) - val(f.denominator)
    active_r, active_c = set(range(n)), set(range(mcols))
    while active_r and active_c:
        best = None
        for i in active_r:
            for j in active_c:
                if F[i][j] != 0:
                    v = fval(F[i][j])
                    if best is None or v < best[0]:
                        best = (v, i, j)
                        if v == 0:
                            break
            if best and best[0] == 0:
                break
        if best is None:
            break
        v, i, j = best
        vals.append(v)
        piv = F[i][j]
        for i2 in active_r:
            if i2 != i and F[i2][j] != 0:
                f = F[i2][j] / piv
                F[i2] = [a - f * b for a, b in zip(F[i2], F[i])]
        for j2 in active_c:
            if j2 != j:
                F[i][j2] = Fraction(0)
        active_r.discard(i); active_c.discard(j)
    zeros = min(n, mcols) - len(vals)
    from collections import Counter
    c = Counter(vals)
    return dict(sorted(c.items())), zeros

if __name__ == '__main__':
    for name, (V, adj), m, r, s in [('H(2,3)', hamming_9(), 2, 1, -2), ('BvLS', bvls_243(), 11, 4, -5)]:
        analyse(name, adj, m)
        n = len(adj)
        A = [[int(j in adj[i]) for j in range(n)] for i in range(n)]
        for th in (r, s):
            M = [[A[i][j] - th * (i == j) for j in range(n)] for i in range(n)]
            print(f'   rank_3(A - ({th})I) = {rank_mod(M, 3)}', end='')
        print()
        if n <= 243:
            for th in (r, s):
                M = [[A[i][j] - th * (i == j) for j in range(n)] for i in range(n)]
                d, z = smith_p_part(M, 3)
                print(f'   Smith 3-part of A - ({th})I: valuation -> count {d}, zero divisors {z}')
