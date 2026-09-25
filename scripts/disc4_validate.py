"""S5 validation on real graphs in the NL / LS families (v = n^2, every prime of n ≡ 1 mod 4), where the criterion bites:
compute Lambda_theta exactly, check it is even with odd det, and that the parameter-level prediction
    D_theta ≡ P_theta (mod 4)  and  (-1)^{m/2} D_theta ≡ 1 (mod 4)
holds.  Negative control: a lattice with the same rank but det ≡ 3 (mod 4) is never produced."""
import itertools
import sys

sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c99_structures_lab\scripts')
import genus_tools as gt
from eigenlattice_check import int_kernel
import genus_scan as gs


def gf_p2(p, c):
    """GF(p^2) = F_p[x]/(x^2 - c); elements (a, b) = a + b x."""
    els = [(a, b) for a in range(p) for b in range(p)]
    def mul(u, w):
        return ((u[0] * w[0] + c * u[1] * w[1]) % p, (u[0] * w[1] + u[1] * w[0]) % p)
    return els, mul


def paley_p2(p):
    c = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
    els, mul = gf_p2(p, c)
    sq = {mul(u, u) for u in els if u != (0, 0)}
    idx = {u: i for i, u in enumerate(els)}
    n = len(els)
    A = [[0] * n for _ in range(n)]
    for u in els:
        for w in els:
            d = ((u[0] - w[0]) % p, (u[1] - w[1]) % p)
            if d in sq:
                A[idx[u]][idx[w]] = 1
    return A


def latin_square_graph(n, m):
    """LS_m(n) from Z_n: rows, columns and m - 2 'diagonal' classes y - a x = const (n prime)."""
    V = [(x, y) for x in range(n) for y in range(n)]
    def cls(u):
        return [('r', u[0]), ('c', u[1])] + [('d%d' % a, (u[1] - a * u[0]) % n) for a in range(1, m - 1)]
    C = [set(cls(u)) for u in V]
    return [[1 if (i != j and C[i] & C[j]) else 0 for j in range(len(V))] for i in range(len(V))]


def params(A):
    v = len(A); k = sum(A[0])
    adj = next(j for j in range(1, v) if A[0][j]); non = next(j for j in range(1, v) if not A[0][j])
    return v, k, sum(A[0][t] * A[adj][t] for t in range(v)), sum(A[0][t] * A[non][t] for t in range(v))


def check(name, A):
    v, k, l, mu = params(A)
    r, s, f, g = gs.feasible(v, k, l, mu)
    ok = True
    for th, oth, m in ((r, s, f), (s, r, g)):
        M = [[A[i][j] - (th if i == j else 0) for j in range(v)] for i in range(v)]
        B = int_kernel(M)
        G = [[sum(a * b for a, b in zip(u, w)) for w in B] for u in B]
        D = int(gt.det_frac(G))
        even = all(G[i][i] % 2 == 0 for i in range(len(G)))
        P = 1
        for q in gs.primes(v):
            if q == 2 or (r - s) % q == 0:
                continue
            if gs.vp(k - oth, q) < gs.vp(v, q):
                P *= q ** gs.vp(v, q)
        pred_ok = (D % 4 == P % 4) if all(q % 4 == 1 for q in gs.primes(r - s)) else None
        disc_ok = ((-1) ** (m // 2) * D) % 4 == 1 if D % 2 else None
        print(f'   {name} srg({v},{k},{l},{mu}) Lambda_{th}: rank {len(B)} (m={m}), even {even}, det = {gs.primes(D)}-number '
              f'≡ {D % 4} (mod 4), predicted P = {P}: {pred_ok}, (-1)^(m/2) det ≡ 1: {disc_ok}')
        ok &= even and (pred_ok is not False) and (disc_ok is not False)
    return ok


if __name__ == '__main__':
    allok = True
    allok &= check('Paley(25)', paley_p2(5))
    allok &= check('LS_3(5)', latin_square_graph(5, 3))
    allok &= check('LS_3(13)', latin_square_graph(13, 3))
    allok &= check('LS_4(13)', latin_square_graph(13, 4))
    allok &= check('Paley(169)', paley_p2(13))
    print('ALL CONSISTENT' if allok else '*** INCONSISTENT ***')
