"""S3 (D11): C1's lemmas L21 and L23 rerun on REAL graphs, in their general form, with negative controls.

Part A — L21 mechanism (ball of a line), for srg with lambda = 1, mu = 2 (srg(9,4,1,2), BvLS(243,22,1,2)), p | r - s.
  L = line (triangle), N(L) = points adjacent to exactly one point of L (|N(L)| = 3(k-2)), G = G(L) induced (cubic,
  one neighbour in each part).  With c = r + 1, B = (A - rI)[L] = J - cI, B^{-1} = -(1/c)(I + J/(c-3)), C = incidence:
      S = A_G - rI + (1/c) (D + J/(c-3))        (D = "same part" matrix),
  rank_p (A - rI)[L ∪ N(L)] = 3 + rank_p S,  U = span of the part indicators is A_G-invariant, S = A_G - rI on U^⊥,
  so rank = 3 + [rank(A_G - rI) - rank((A_G - rI)|_U)] + rank(S|_U).
  For 99 (p=7): rank(S|_U) = 1, rank((A_G - 3I)|_U) = 2  -> C1's "2 + rank(A_G - 3I)".
  For BvLS and srg9 (p=3): S|_U = 0, (A_G - rI)|_U invertible -> prediction rank = rank_p(A_G - rI).
  Also checks van Eijl: rank_p(A - rI) >= local rank for every line.
  Negative control: drop the J/(c-3) term (wrong Schur complement) -> identity must fail.
Part B — L23 principle, for srg with p || d, p ∤ v (triangular graphs T(n), n odd):
  Q = A - rI - ((k - r)/v) J (mod p) satisfies Q^2 ≡ 0, rank_p Q = alpha; for any S with rank_p Q_SS = alpha, the
  non-degenerate part of Q_SS has the discriminant class of V (p-modular part of Lambda_s), since alpha is even.
  Negative control: another multiple of J (Q^2 ≢ 0) gives no constant class."""
import itertools
import random
import sys

import numpy as np

sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c99_structures_lab\scripts')
import graphs
import genus_tools as gt
from eigenlattice_check import int_kernel
import genus_scan as gs


def rank_p(M, p):
    M = np.array(M, dtype=np.int64) % p
    m, n = M.shape
    r = 0
    for c in range(n):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0:
            continue
        i = r + piv[0]
        M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, p)) % p
        col = M[:, c].copy(); col[r] = 0
        M = (M - np.outer(col, M[r])) % p
        r += 1
        if r == m:
            break
    return r


def disc_class(M, p):
    """Legendre class of the discriminant of the non-degenerate part of a symmetric matrix over F_p (odd p)."""
    M = [[x % p for x in row] for row in M]
    n = len(M)
    idx = list(range(n))
    prod = 1
    while True:
        # find a non-isotropic basis vector (diagonal) or a hyperbolic pair
        piv = next((i for i in idx if M[i][i] % p), None)
        if piv is None:
            pair = next(((i, j) for i in idx for j in idx if i < j and M[i][j] % p), None)
            if pair is None:
                break
            i, j = pair
            for t in range(n):           # e_i <- e_i + e_j  (M_ii = M_jj = 0, M_ij != 0 gives 2 M_ij != 0)
                M[i][t] = (M[i][t] + M[j][t]) % p
            for t in range(n):
                M[t][i] = (M[t][i] + M[t][j]) % p
            piv = i
        a = M[piv][piv]
        prod = prod * a % p
        inv = pow(a, -1, p)
        rest = [t for t in idx if t != piv]
        for s_ in rest:
            f = M[s_][piv] * inv % p
            if f:
                for t in range(n):
                    M[s_][t] = (M[s_][t] - f * M[piv][t]) % p
        for s_ in rest:
            M[s_][piv] = M[piv][s_] = 0
        idx = rest
    return 1 if pow(prod, (p - 1) // 2, p) == 1 else -1


def part_A(name, A, p):
    A = np.array(A, dtype=np.int64)
    v = len(A); k = int(A[0].sum())
    r, s, f, g = gs.feasible(v, k, 1, 2)
    c = r + 1
    inv = lambda x: pow(x % p, -1, p)
    glob = rank_p(A - r * np.eye(v, dtype=np.int64), p)
    lines = {tuple(sorted((x, y, z))) for x in range(v) for y in np.nonzero(A[x])[0] for z in np.nonzero(A[x] * A[y])[0]}
    ok = neg_fail = 0
    locmax = 0
    for L in lines:
        L = list(L)
        Nset = [y for y in range(v) if y not in L and int(A[y, L].sum()) == 1]
        part = {y: next(a for a in L if A[y, a]) for y in Nset}
        G = A[np.ix_(Nset, Nset)]
        assert all(int(G[i].sum()) == 3 for i in range(len(Nset))), 'G(L) not cubic'
        D = np.array([[1 if part[y] == part[z] else 0 for z in Nset] for y in Nset], dtype=np.int64)
        J = np.ones_like(D)
        S = (G - r * np.eye(len(Nset), dtype=np.int64) + inv(c) * (D + inv(c - 3) * J)) % p
        S_bad = (G - r * np.eye(len(Nset), dtype=np.int64) + inv(c) * D) % p
        idx = L + Nset
        loc = rank_p((A - r * np.eye(v, dtype=np.int64))[np.ix_(idx, idx)], p)
        pred = rank_p(G - r * np.eye(len(Nset), dtype=np.int64), p)      # prediction for S|_U = 0 and (A_G - r)|_U invertible
        ok += (loc == 3 + rank_p(S, p)) and (loc == pred)
        neg_fail += loc != 3 + rank_p(S_bad, p)
        locmax = max(locmax, loc)
    print(f'  {name}: p={p}, {len(lines)} lines, global rank_p(A - {r}I) = {glob}; identity rank = 3 + rank S and '
          f'= rank(A_G - {r}I) on {ok}/{len(lines)} lines; negative control (no J term) fails on {neg_fail}/{len(lines)}; '
          f'max local rank {locmax} <= global {glob}: {locmax <= glob}')
    return ok == len(lines) and locmax <= glob


def triangular(n):
    V = list(itertools.combinations(range(n), 2))
    return [[1 if (a != b and set(a) & set(b)) else 0 for b in V] for a in V]


def part_B(n, samples=60, seed=1):
    A = triangular(n); v = len(A); k = 2 * (n - 2)
    r, s, f, g = gs.feasible(v, k, n - 2, 4)
    d = r - s
    rng = random.Random(seed)
    out_ok = True
    for p in gs.primes(d):
        if d % (p * p) == 0 or v % p == 0 or p == 2:
            continue
        # genus sign of V from the exact lattice Lambda_s
        B = int_kernel([[A[i][j] - (s if i == j else 0) for j in range(v)] for i in range(v)])
        Gm = [[sum(a * b for a, b in zip(u, w)) for w in B] for u in B]
        comp = {j: (nn, eps) for j, nn, eps in gt.jordan_symbol(Gm, p)}
        alpha, epsV = comp.get(1, (0, 1))
        cJ = (-(k - r) * pow(v, -1, p)) % p
        An = np.array(A, dtype=np.int64)
        for cc, label in ((cJ, 'Q'), ((cJ + 1) % p, 'wrong J')):
            Q = (An - r * np.eye(v, dtype=np.int64) + cc * np.ones((v, v), dtype=np.int64)) % p
            sq0 = not ((Q @ Q) % p).any()
            rq = rank_p(Q, p)
            classes = []
            tries = 0
            while len(classes) < samples and tries < 20 * samples:
                tries += 1
                S = sorted(rng.sample(range(v), rng.randint(max(alpha + 2, v // 3), v)))
                QS = Q[np.ix_(S, S)]
                if rank_p(QS, p) == rq:
                    classes.append(disc_class(QS.tolist(), p))
            const = len(set(classes)) == 1
            print(f'  T({n}) p={p} {label}: Q^2 ≡ 0: {sq0}, rank Q = {rq} (alpha = {alpha}), full-rank samples {len(classes)}, '
                  f'disc classes {sorted(set(classes))}, genus eps(V) = {epsV}' +
                  (f' -> {"MATCH" if const and classes and classes[0] == epsV else "MISMATCH"}' if label == 'Q' else ''))
            if label == 'Q':
                out_ok &= sq0 and rq == alpha and const and bool(classes) and classes[0] == epsV
    return out_ok


if __name__ == '__main__':
    allok = True
    print('Part A (L21 mechanism, ball of a line):')
    allok &= part_A('srg(9,4,1,2)', graphs.srg9(), 3)
    allok &= part_A('BvLS(243,22,1,2)', graphs.bvls(), 3)
    print('Part B (L23 principle, sign of full-rank principal blocks):')
    for n in (7, 9, 13, 15):
        allok &= part_B(n)
    print('ALL CONSISTENT' if allok else '*** SOME INCONSISTENCY ***')
