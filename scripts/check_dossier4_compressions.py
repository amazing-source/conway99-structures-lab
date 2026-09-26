"""Exact check (integers/fractions, no solver) of the P-compressions in Astra's dossier-4 report, sections 6.1-6.3.
P = 6I - 2Q - Sh + J on orbits; P^2 = 14P, so every principal submatrix of P is PSD.
Entries on U, W, A, B, V (U = {i,j}, W = D(U) = {k,l}, A = interior apex {i,k}, B = other apex, V = D(A)):
  Q_UW = 2, Q_UA = Q_WA = Q_UB = Q_WB = 1, Q_AB = 0, Q_AV = 2, Q_UV = Q_WV = 0 (Gram budget), Q_BV in {0, 1}.
For each non-PSD case, prints a short integer vector x with x^T P x < 0 (a hand-checkable certificate)."""
import itertools
from fractions import Fraction

def pmat(names, supp, Q):
    n = len(names)
    M = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a == b:
                M[a][b] = 5                                   # 6 - 0 - 2 + 1
            else:
                q = Q.get(frozenset((names[a], names[b])), 0)
                M[a][b] = -2 * q - len(supp[names[a]] & supp[names[b]]) + 1
    return M

def det(M):
    M = [[Fraction(v) for v in row] for row in M]
    n, d = len(M), Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d

def psd(M):
    n = len(M)
    for k in range(1, n + 1):
        for S in itertools.combinations(range(n), k):
            if det([[M[a][b] for b in S] for a in S]) < 0:
                return False
    return True

def cert(M):
    """Integer x with x^T M x < 0, via the Schur complement on a positive definite leading block (U, W, A)."""
    n = len(M)
    F = [[Fraction(v) for v in row] for row in M]
    m = 3                                                    # leading block U, W, A (checked PD)
    A11 = [row[:m] for row in F[:m]]
    # inverse of the 3 x 3 block by Gauss-Jordan
    aug = [A11[r][:] + [Fraction(int(r == c)) for c in range(m)] for r in range(m)]
    for c in range(m):
        p = next(r for r in range(c, m) if aug[r][c] != 0)
        aug[c], aug[p] = aug[p], aug[c]
        piv = aug[c][c]
        aug[c] = [v / piv for v in aug[c]]
        for r in range(m):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [a - f * b for a, b in zip(aug[r], aug[c])]
    inv = [row[m:] for row in aug]
    rest = list(range(m, n))
    X = [[sum(inv[a][c] * F[c][b] for c in range(m)) for b in rest] for a in range(m)]    # A11^-1 A12
    S = [[F[a][b] - sum(F[a][c] * X[c][rest.index(b)] for c in range(m)) for b in rest] for a in rest]
    k = len(rest)
    y = None
    for t in range(k):
        if S[t][t] < 0:
            y = [Fraction(int(s == t)) for s in range(k)]
            break
    if y is None and k == 2:
        s11, s12 = S[0][0], S[0][1]
        y = [-s12, s11]                                      # y^T S y = s11 * det S
    if y is None:
        return None
    x = [-sum(X[a][t] * y[t] for t in range(k)) for a in range(m)] + y
    from math import lcm, gcd
    L = 1
    for v in x:
        L = lcm(L, v.denominator)
    xi = [int(v * L) for v in x]
    g = 0
    for v in xi:
        g = gcd(g, abs(v))
    xi = [v // g for v in xi]
    val = sum(xi[a] * M[a][b] * xi[b] for a in range(n) for b in range(n))
    return (None, tuple(xi), val, S)

base = {'U': {'i', 'j'}, 'W': {'k', 'l'}, 'A': {'i', 'k'}}
Q0 = {frozenset(p): 2 for p in [('U', 'W')]}
Q0.update({frozenset(p): 1 for p in [('U', 'A'), ('W', 'A')]})
print('6.1  P[U,W,A] =', pmat(['U', 'W', 'A'], base, Q0), ' PSD:', psd(pmat(['U', 'W', 'A'], base, Q0)))
# 6.1: second apex B (Q_UB = Q_WB = 1, Q_AB = 0) with various supports
for sB in [{'j', 'p'}, {'l', 'p'}, {'i', 'p'}, {'k', 'p'}, {'j', 'l'}, {'i', 'l'}, {'j', 'k'}, {'p', 'q'}]:
    supp = dict(base, B=sB)
    Q = dict(Q0); Q[frozenset(('U', 'B'))] = 1; Q[frozenset(('W', 'B'))] = 1
    M = pmat(['U', 'W', 'A', 'B'], supp, Q)
    ok = psd(M)
    print(f'6.1  B = {sorted(sB)}: PSD {ok}', '' if ok else f'certificate {cert(M)[1:]}')
# 6.2-6.3: V = D(A), support disjoint from {i,k}
supp = dict(base, B={'p', 'q'})
Q = dict(Q0)
Q.update({frozenset(('U', 'B')): 1, frozenset(('W', 'B')): 1, frozenset(('A', 'V')): 2})
print('6.2  P[U,W,A,V] with V = {j,l}:', pmat(['U', 'W', 'A', 'V'], dict(supp, V={'j', 'l'}), Q))
for sV in itertools.combinations(['j', 'l', 'p', 'q', 'r'], 2):
    res = []
    for qbv in (0, 1):
        Qc = dict(Q); Qc[frozenset(('B', 'V'))] = qbv
        M = pmat(['U', 'W', 'A', 'B', 'V'], dict(supp, V=set(sV)), Qc)
        ok = psd(M)
        res.append(f'Q_BV={qbv}: ' + ('PSD' if ok else f'NOT PSD, x={cert(M)[1]} gives {cert(M)[2]}'))
    print('6.3  supp V =', ''.join(sV), ' | ', ' ; '.join(res))

# Cross-check of the two tables of NOTE_S23 (Schur residues, exact fractions)
def schur(M, m=3):
    F = [[Fraction(v) for v in row] for row in M]
    n = len(F)
    K = [row[:m] for row in F[:m]]
    dK = det(K)
    # adjugate-free inverse via Cramer on 3x3
    inv = [[None] * m for _ in range(m)]
    for r in range(m):
        for c in range(m):
            minor = [[K[a][b] for b in range(m) if b != r] for a in range(m) if a != c]
            inv[r][c] = (-1) ** (r + c) * det(minor) / dK
    rest = list(range(m, n))
    return [[F[a][b] - sum(F[a][c1] * inv[c1][c2] * F[c2][b] for c1 in range(m) for c2 in range(m)) for b in rest] for a in rest]

print('--- table 6.1 (residue of B after U, W, A) ---')
for sB in [{'p', 'q'}, {'j', 'p'}, {'l', 'p'}, {'i', 'p'}, {'k', 'p'}, {'j', 'l'}, {'i', 'l'}, {'j', 'k'}, {'i', 'k'}]:
    Q = dict(Q0); Q[frozenset(('U', 'B'))] = 1; Q[frozenset(('W', 'B'))] = 1
    M = pmat(['U', 'W', 'A', 'B'], dict(base, B=sB), Q)
    print(sorted(sB), schur(M)[0][0])
print('--- table 6.3 (S on B, V: y = S_VV, x = S_BV, 3y - x^2) ---')
for sV in itertools.combinations(['j', 'l', 'p', 'q', 'r'], 2):
    out = []
    for qbv in (0, 1):
        Qc = dict(Q); Qc[frozenset(('B', 'V'))] = qbv
        Qc.update({frozenset(('U', 'B')): 1, frozenset(('W', 'B')): 1, frozenset(('A', 'V')): 2})
        S = schur(pmat(['U', 'W', 'A', 'B', 'V'], dict(supp, V=set(sV)), Qc))
        out.append(f'Q_BV={qbv}: S_BB={S[0][0]}, x={S[0][1]}, y={S[1][1]}, 3y-x^2={S[0][0]*S[1][1]-S[0][1]**2}')
    print(''.join(sV), ' | ', ' ; '.join(out))
