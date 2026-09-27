"""Calibration of the 2-adic even/odd coupling (C1's idea (b), 27/09) on the two EXISTING graphs with an involution fixing one
vertex: Paley 9 = H(2,3) (m = 2, t = 3) and BvLS (m = 11, t = 9), sigma = -1.  Any invariant claimed to exclude m = 7 must
hold on these two unless it uses m-specific numbers.
For theta in {r, s} (the non-principal eigenvalues):
  Lambda   = Z^n kernel of (A - theta I)                                (full eigenlattice, standard form)
  K        = Z^{m^2} kernel of (B - theta I), B = odd part on orbits      (Lambda_odd = sqrt2 K)
  Lambda_e = integer kernel of the even action on (c0, c) (x and the orbit sums f_i), form c0^2 + 2|c|^2
Reports rank, det, parity, odd part of det mod 8, the index [Lambda : Lambda_e + Lambda_odd], and for K mod 2: the code
dimension and G2 = sum over the code of i^wt.  Usage: python calib_2adic.py"""
import sys, math
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, '../famille')
from family_graphs import hamming_9, bvls_243


def int_kernel(M):
    """Z-basis of {x in Z^n : M x = 0} for an integer matrix M (rows), via column-style elimination on [M^T | I]"""
    M = [list(map(int, r)) for r in M]
    nr, n = len(M), len(M[0])
    rows = [[M[i][j] for i in range(nr)] + [int(j == k) for k in range(n)] for j in range(n)]   # row j = (column j of M | e_j)
    piv_row = 0
    for c in range(nr):
        # gcd-reduce column c among rows >= piv_row
        while True:
            nz = [r for r in range(piv_row, n) if rows[r][c]]
            if len(nz) <= 1:
                break
            p = min(nz, key=lambda r: abs(rows[r][c]))
            for r in nz:
                if r != p:
                    q = rows[r][c] // rows[p][c]
                    rows[r] = [a - q * b for a, b in zip(rows[r], rows[p])]
        nz = [r for r in range(piv_row, n) if rows[r][c]]
        if nz:
            p = nz[0]
            rows[piv_row], rows[p] = rows[p], rows[piv_row]
            piv_row += 1
    return [r[nr:] for r in rows[piv_row:]]


def lll_rows(B):
    """size-reduce/LLL integer rows to keep numbers small (Gram-based float LLL on the standard form)"""
    B = [list(r) for r in B]
    n = len(B)
    changed = True
    while changed:
        changed = False
        G = np.array([[float(sum(a * b for a, b in zip(B[i], B[j]))) for j in range(n)] for i in range(n)])
        for i in range(n):
            for j in range(n):
                if i != j and G[j, j] > 0:
                    q = int(round(G[i, j] / G[j, j]))
                    if q and abs(G[i, j]) > G[j, j] / 2 + 1e-9:
                        B[i] = [a - q * b for a, b in zip(B[i], B[j])]
                        changed = True
                        break
            if changed:
                break
    return B


def det_exact(G):
    M = [list(map(int, r)) for r in G]
    n, sign, prev = len(M), 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k]), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def v2(x):
    x = abs(x); k = 0
    while x and x % 2 == 0:
        x //= 2; k += 1
    return k


def describe(name, G):
    d = det_exact(G)
    odd = any(G[i][i] % 2 for i in range(len(G)))
    k = v2(d)
    u = abs(d) // (2 ** k)
    return '%s: rank %d, det = 2^%d * %d (odd part mod 8 = %d), %s' % (name, len(G), k, u, u % 8, 'odd' if odd else 'even')


def code_info(Kb):
    """binary code spanned by the rows mod 2: dimension and G2 = sum i^wt over the code"""
    rows = []
    for r in Kb:
        v = 0
        for i, x in enumerate(r):
            if x % 2:
                v |= 1 << i
        for b in rows:
            v = min(v, v ^ b)
        if v:
            rows.append(v)
    dim = len(rows)
    if dim > 26:
        return dim, None
    cnt = [0, 0, 0, 0]
    words = [0]
    for b in rows:
        words += [w ^ b for w in words]
    for w in words:
        cnt[bin(w).count('1') % 4] += 1
    G2 = complex(cnt[0] - cnt[2], cnt[1] - cnt[3])
    return dim, (cnt, G2)


def run(name, V, adj):
    n = len(V)
    idx = {v: i for i, v in enumerate(V)}
    neg = [idx[tuple((-x) % 3 for x in v)] for v in V]
    A = np.zeros((n, n), dtype=np.int64)
    for i, nb in enumerate(adj):
        for j in nb:
            A[i, j] = 1
    x0 = idx[tuple([0] * len(V[0]))]
    k = int(A[x0].sum()); m = k // 2
    t = int(round(math.sqrt(8 * m - 7))); r, s = (t - 1) // 2, (-1 - t) // 2
    reps, seen = [], set()
    for v in range(n):
        if v != x0 and v not in seen:
            reps.append(v); seen |= {v, neg[v]}
    N = len(reps)
    B = np.array([[A[u, w] - A[u, neg[w]] for w in reps] for u in reps], dtype=np.int64)
    Q = np.array([[A[u, w] + A[u, neg[w]] for w in reps] for u in reps], dtype=np.int64)
    abar = np.array([1 if A[x0, u] else 0 for u in reps], dtype=np.int64)
    print('== %s: n = %d, m = %d, t = %d, eigenvalues %d, %d, %d; %d orbit pairs' % (name, n, m, t, k, r, s, N), flush=True)
    for th in (r, s):
        Lb = int_kernel((A - th * np.eye(n, dtype=np.int64)).tolist())
        Lb = lll_rows(Lb)
        GL = [[sum(a * b for a, b in zip(p, q)) for q in Lb] for p in Lb]
        Kb = lll_rows(int_kernel((B - th * np.eye(N, dtype=np.int64)).tolist()))
        GK = [[2 * sum(a * b for a, b in zip(p, q)) for q in Kb] for p in Kb]
        # even action on (c0, c): equations theta c0 - 2 abar.c = 0 and c0 abar + Q c - theta c = 0
        Me = [[th] + [-2 * int(a) for a in abar]] + [[int(abar[i])] + [int(Q[i, j]) - (th if i == j else 0) for j in range(N)] for i in range(N)]
        Eb = lll_rows(int_kernel(Me))
        GE = [[p[0] * q[0] + 2 * sum(a * b for a, b in zip(p[1:], q[1:])) for q in Eb] for p in Eb]
        dL, dK, dE = det_exact(GL), det_exact(GK), det_exact(GE)
        idx2 = F(dE * dK, dL)
        print('  theta = %d:' % th)
        print('    ' + describe('Lambda (full)', GL))
        print('    ' + describe('Lambda_odd = sqrt2 K', GK))
        print('    ' + describe('Lambda_even', GE))
        print('    index^2 [Lambda : even + odd]^2 = %s' % idx2)
        dim, ci = code_info(Kb)
        if ci:
            cnt, G2 = ci
            print('    K mod 2: dim %d, weights mod 4 counts %s, G2 = %s' % (dim, cnt, G2))
        else:
            print('    K mod 2: dim %d (G2 not enumerated)' % dim)
        sys.stdout.flush()


if __name__ == '__main__':
    which = sys.argv[1:] or ['paley9']
    if 'paley9' in which:
        run('Paley 9', *hamming_9())
    if 'bvls' in which:
        run('BvLS', *bvls_243())
