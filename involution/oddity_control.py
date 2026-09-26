"""Positive control of the genus/oddity machinery on real srgs with the involution x -> -x (one fixed point):
H(2,3) (m = 2) and BvLS (m = 11).  Odd part: basis d_o = e_o - e_{-o}, form = standard/2 (Z^n unimodular, odd, characteristic
vector 1).  For the eigenlattice Lam_s (s = the negative eigenvalue) we check the oddity formula
    rank + p-excess(Lam_s) == |E_s 1|^2  (mod 8)
(2-adic oddity from the projected characteristic vector; p-adic Jordan by exact diagonalisation).  Own code, no sympy."""
import sys
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c99_structures_lab\famille')
from family_graphs import hamming_9, bvls_243

def odd_part(V, adj):
    idx = {v: i for i, v in enumerate(V)}
    neg = lambda v: tuple((-x) % 3 for x in v)
    reps, seen = [], set()
    for v in V:
        if any(v) and v not in seen:
            reps.append(v); seen |= {v, neg(v)}
    return [[int(idx[rb] in adj[idx[ra]]) - int(idx[neg(rb)] in adj[idx[ra]]) for rb in reps] for ra in reps]

def kernel_Z(M):
    """Z-basis (saturated) of {x in Z^n : M x = 0}, by unimodular column reduction."""
    m, n = len(M), len(M[0])
    A = [row[:] for row in M]                   # we column-reduce A, tracking U (n x n)
    U = [[int(i == j) for j in range(n)] for i in range(n)]
    def colop_add(dst, src, f):                 # col dst += f * col src
        for r in range(m): A[r][dst] += f * A[r][src]
        for r in range(n): U[r][dst] += f * U[r][src]
    def colswap(a, b):
        for r in range(m): A[r][a], A[r][b] = A[r][b], A[r][a]
        for r in range(n): U[r][a], U[r][b] = U[r][b], U[r][a]
    piv = 0
    for r in range(m):
        if piv >= n: break
        while True:
            nz = [c for c in range(piv, n) if A[r][c] != 0]
            if not nz: break
            c0 = min(nz, key=lambda c: abs(A[r][c]))
            colswap(piv, c0)
            done = True
            for c in range(piv + 1, n):
                if A[r][c] != 0:
                    q = A[r][c] // A[r][piv]
                    colop_add(c, piv, -q)
                    if A[r][c] != 0: done = False
            if done: break
        if any(A[r][c] != 0 for c in range(piv, piv + 1)):
            piv += 1
    return [[U[r][c] for r in range(n)] for c in range(piv, n)]

def solve(G, b):
    n = len(G)
    M = [[Fraction(x) for x in G[i]] + [Fraction(b[i])] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        inv = 1 / M[c][c]; M[c] = [x * inv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][n] for i in range(n)]

def det_frac(G):
    M = [[Fraction(x) for x in row] for row in G]; n = len(M); d = Fraction(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None: return Fraction(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            if f: M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d

def padic_jordan(G, p):
    M = [[Fraction(x) for x in row] for row in G]; n = len(M)
    def val(x):
        if x == 0: return 10 ** 9
        num, den, v = x.numerator, x.denominator, 0
        while num % p == 0: num //= p; v += 1
        while den % p == 0: den //= p; v -= 1
        return v
    active = list(range(n)); diag = []
    while active:
        vd, i = min((val(M[i][i]), i) for i in active)
        off = min(((val(M[a][b]), a, b) for a in active for b in active if a < b), default=(10 ** 9, None, None))
        if off[0] < vd:
            _, a, b = off
            for k in range(n): M[a][k] += M[b][k]
            for k in range(n): M[k][a] += M[k][b]
            continue
        piv = M[i][i]; diag.append(piv)
        for j in active:
            if j != i and M[j][i] != 0:
                f = M[j][i] / piv
                for k in range(n): M[j][k] -= f * M[i][k]
                for k in range(n): M[k][j] -= f * M[k][i]
        active.remove(i)
    comp = defaultdict(list)
    for a in diag:
        v = val(a); comp[v].append(a / Fraction(p) ** v)
    out = []
    for v in sorted(comp):
        prod = Fraction(1)
        for u in comp[v]: prod *= u
        num = prod.numerator * pow(prod.denominator, -1, p) % p
        out.append((v, len(comp[v]), 1 if pow(num, (p - 1) // 2, p) == 1 else -1))
    return out

def excess(jordan, p):
    return sum(d * (p ** v - 1) + (4 if v % 2 == 1 and leg == -1 else 0) for v, d, leg in jordan) % 8

def control(name, V, adj, s, p):
    S = odd_part(V, adj); n = len(S)
    B = kernel_Z([[S[i][j] - s * (i == j) for j in range(n)] for i in range(n)])
    assert all(sum(S[i][k] * b[k] for k in range(n)) == s * b[i] for b in B for i in range(n))
    G = [[sum(x * y for x, y in zip(a, b)) for b in B] for a in B]
    jord = padic_jordan(G, p)
    rhs = [sum(b) for b in B]                     # B . 1
    coef = solve(G, rhs)
    norm = sum(c * r for c, r in zip(coef, rhs))  # |E_s 1|^2 = (B1)^T G^-1 (B1)
    odd2 = (norm.numerator * pow(norm.denominator, -1, 8)) % 8
    lhs = (len(B) + excess(jord, p)) % 8
    print(f'{name}: odd dim {n}, rank Lam_{s} = {len(B)}, det = {det_frac(G)}, {p}-adic Jordan {jord}; '
          f'|E 1|^2 = {norm}; rank + excess = {lhs}; oddity_2 = {odd2}; formula holds: {lhs == odd2}', flush=True)

V, adj = hamming_9(); control('H(2,3), x -> -x', V, adj, -2, 3)
V, adj = bvls_243(); control('BvLS, x -> -x', V, adj, -5, 3)
