"""Block line of an A-block (all 6 pairs at c anchored at c): tau_c = iota_c (every orbit through c is in a digon of R_c).
Rank_7 of the 19x19 odd block-line matrix, and when it is 18, the quadratic character of the determinant of the 18x18 block
obtained by deleting d_c (the radical coordinate z_c), i.e. the local sign of NOTE_S27 par. 7."""
from collections import Counter
from blockline_exhaustive import labels, mate, rank_mod, matrix, fpf_involutions
P = 7
def det_mod(M, p=P):
    M = [[x % p for x in r] for r in M]; n = len(M); d = 1
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c]), None)
        if piv is None: return 0
        if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
        d = d * M[c][c] % p; inv = pow(M[c][c], -1, p)
        for i in range(c + 1, n):
            if M[i][c]:
                f = M[i][c] * inv % p
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[c])]
    return d % p
def chi(a): return 0 if a % P == 0 else (1 if pow(a, (P - 1) // 2, P) == 1 else -1)
h = Counter()
for tau in fpf_involutions(labels):
    if any(tau[v] == mate[v] for v in labels): continue      # no sibling pair at c (Sigma = empty)
    M = matrix(tau, tau)
    r = rank_mod(M)
    s = chi(det_mod([row[1:] for row in M[1:]])) if r == 18 else None
    h[(r, s)] += 1
print('A-block (tau = iota, no mate digon): (rank, local sign) ->', dict(h))
