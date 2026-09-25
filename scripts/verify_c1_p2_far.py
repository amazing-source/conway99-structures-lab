"""Independent control (numeric, exact integers where possible) of C1's notes/116 tools, on REAL graphs of the family
srg(v,k,1,2): srg(9,4,1,2) and BvLS srg(243,22,1,2).  The derivations were redone by hand in the lab (NOTE_S12);
this script checks the general-k forms:
  (1) at a vertex x: A_D^2 + A_D - (k-2) I = 2J - d^T d (D = far set, d = label incidence);
      (P2) sum_z cnt(z,R) cnt(z,R') = (k-2)|R n R'| - N(R,R') + sum_{y in R, w in R'} (2 - |y n w|), random R, R';
      diagonal on class-pair blocks: sum cnt = 4(k-2), sum cnt(cnt-1) = 16 - 2 e(R);
  (2) at a line L: M M^T = (k-2) I + 2J - A_G - A_G^2 - K ; A_F M^T = M^T (J/3 - A_G - I);
      exact forms of (U) and (L):  (r-s) |Pi_s 1_S|^2 = r|S| + n^T Y n - 2e(S),
                                   (r-s) |Pi_r 1_S|^2 = 2e(S) - s|S| - n^T Y' n,
      with Y = (MM^T)^+ (J/3 - A_G - (1+r)I), Y' = (MM^T)^+ (J/3 - A_G - (1+s)I); equality at S = F.
For 99 (k = 14, r = 3, s = -4) these give exactly C1's closed forms: 1/12 and 1/5 on 1, -(3-A_G)^-1 and (A_G+4)^-1
on the complement of 1 and U."""
import itertools
import random
import sys

import numpy as np

sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
import graphs


def check_vertex(name, A, trials=200, seed=1):
    v = A.shape[0]
    k = int(A[0].sum())
    rng = random.Random(seed)
    x = 0
    Nx = [u for u in range(v) if A[x, u]]
    far = [u for u in range(v) if u != x and not A[x, u]]
    lab = {w: tuple(sorted(u for u in Nx if A[w, u])) for w in far}
    assert all(len(l) == 2 for l in lab.values())
    AD = A[np.ix_(far, far)].astype(np.int64)
    d = np.array([[1 if u in lab[w] else 0 for w in far] for u in Nx], dtype=np.int64)
    lhs = AD @ AD + AD - (k - 2) * np.eye(len(far), dtype=np.int64)
    rhs = 2 * np.ones((len(far), len(far)), dtype=np.int64) - d.T @ d
    assert (lhs == rhs).all(), name
    idx = {w: t for t, w in enumerate(far)}
    for _ in range(trials):
        R = rng.sample(far, rng.randint(1, min(8, len(far))))
        R2 = rng.sample(far, rng.randint(1, min(8, len(far))))
        cR = AD[:, [idx[w] for w in R]].sum(axis=1)
        cR2 = AD[:, [idx[w] for w in R2]].sum(axis=1)
        left = int(cR @ cR2)
        NRR = int(sum(A[y, w] for y in R for w in R2))
        right = (k - 2) * len(set(R) & set(R2)) - NRR + sum(2 - len(set(lab[y]) & set(lab[w])) for y in R for w in R2)
        assert left == right, (name, left, right)
    # classes of N(x): partners u ~ u'
    classes = []
    for u in Nx:
        if not any(u in c for c in classes):
            (up,) = [w for w in Nx if w != u and A[u, w]]
            classes.append((u, up))
    diag = []
    for c1, c2 in itertools.combinations(classes, 2):
        R = [w for w in far if len(set(lab[w]) & set(c1)) == 1 and len(set(lab[w]) & set(c2)) == 1]
        assert len(R) == 4
        cR = AD[:, [idx[w] for w in R]].sum(axis=1)
        eR = sum(A[y, w] for y, w in itertools.combinations(R, 2))
        assert cR.sum() == 4 * (k - 2) and int((cR * (cR - 1)).sum()) == 16 - 2 * eR, name
        diag.append(eR)
    print(f'(1) {name}: A_D^2 + A_D - (k-2)I = 2J - d^T d exact; (P2) on {trials} random (R,R\') exact; '
          f'diagonal on {len(diag)} class-pair blocks, sum cnt(cnt-1) = 16 - 2e(R) with e(R) in {sorted(set(diag))}')


def check_line(name, A, nlines=4, trials=60, seed=2):
    v = A.shape[0]
    k = int(A[0].sum())
    ev = np.round(np.linalg.eigvalsh(A.astype(float)), 6)
    r, s = sorted(set(ev) - {k})[::-1]
    rng = random.Random(seed)
    lines = [t for t in itertools.combinations(range(v), 3)
             if A[t[0], t[1]] and A[t[0], t[2]] and A[t[1], t[2]]] if v < 30 else None
    done = 0
    for _ in range(nlines):
        if lines is not None:
            L = rng.choice(lines)
        else:
            x = rng.randrange(v)
            a = rng.choice([u for u in range(v) if A[x, u]])
            (ab,) = [u for u in range(v) if A[x, u] and A[a, u]]
            L = (x, a, ab)
        parts = [[u for u in range(v) if A[p, u] and u not in L] for p in L]
        N = parts[0] + parts[1] + parts[2]
        F = [u for u in range(v) if u not in L and u not in N]
        if not F:
            continue
        AG = A[np.ix_(N, N)].astype(float)
        M = A[np.ix_(N, F)].astype(float)
        AF = A[np.ix_(F, F)].astype(float)
        K = np.zeros((len(N), len(N)))
        off = 0
        for P in parts:
            K[off:off + len(P), off:off + len(P)] = 1
            off += len(P)
        J = np.ones((len(N), len(N)))
        I = np.eye(len(N))
        assert np.allclose(M @ M.T, (k - 2) * I + 2 * J - AG - AG @ AG - K), name
        assert np.allclose(AF @ M.T, M.T @ (J / 3 - AG - I)), name
        MMp = np.linalg.pinv(M @ M.T)
        Y = MMp @ (J / 3 - AG - (1 + r) * I)
        Y2 = MMp @ (J / 3 - AG - (1 + s) * I)
        # eigenvectors of A supported on F: M f = 0 and A_F f = theta f
        w_, U_ = np.linalg.eigh(AF)
        def proj(theta):
            cols = U_[:, np.abs(w_ - theta) < 1e-8]
            if cols.shape[1] == 0:
                return np.zeros((len(F), len(F)))
            # keep the part inside ker M
            Z = cols @ np.linalg.svd(M @ cols)[2][np.linalg.matrix_rank(M @ cols, tol=1e-8):].T
            return Z @ np.linalg.pinv(Z)
        Ps, Pr = proj(s), proj(r)
        for t in range(trials + 1):
            S = F if t == trials else rng.sample(F, rng.randint(1, len(F) - 1))
            f = np.array([1.0 if u in set(S) else 0.0 for u in F])
            n = M @ f
            e2 = f @ AF @ f
            assert abs((r - s) * f @ Ps @ f - (r * len(S) + n @ Y @ n - e2)) < 1e-6, (name, 'U exact')
            assert abs((r - s) * f @ Pr @ f - (e2 - s * len(S) - n @ Y2 @ n)) < 1e-6, (name, 'L exact')
        assert abs(r * len(F) + (M @ np.ones(len(F))) @ Y @ (M @ np.ones(len(F))) - np.ones(len(F)) @ AF @ np.ones(len(F))) < 1e-6
        done += 1
    print(f'(2) {name} (r={r:g}, s={s:g}): on {done} lines, MM^T and A_F M^T identities exact; exact forms of (U) and (L) '
          f'on {trials} random S each, equality at S = F')


if __name__ == '__main__':
    A9 = graphs.srg9()
    AB = graphs.bvls()
    check_vertex('srg(9,4,1,2)', np.array(A9, dtype=np.int64))
    check_vertex('BvLS srg(243,22,1,2)', np.array(AB, dtype=np.int64))
    check_line('BvLS srg(243,22,1,2)', np.array(AB, dtype=np.int64))
    # closed forms for 99: k = 14, r = 3, s = -4
    from fractions import Fraction as Fr
    k, r, s = 14, 3, -4
    on1_U = Fr(2 * (k - 6 - r), (k - 2) * (k - 4))
    on1_L = Fr(2 * (k - 6 - s), (k - 2) * (k - 4))
    assert on1_U == Fr(1, 12) and on1_L == Fr(1, 5) and (k - 2) * 1 == -(r * s) and r + s == -1
    print('(3) 99 closed forms: coefficients of |S|^2 are 1/12 in (U) and 1/5 in (L); on the complement '
          '(MM^T)^+ = [(3-A_G)(A_G+4)]^-1, Y = -(3-A_G)^-1, Y\' = (A_G+4)^-1 = Y + 7P')
    print('ALL CHECKS PASSED')
