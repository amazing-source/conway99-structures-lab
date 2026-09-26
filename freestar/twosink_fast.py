"""Vectorised exhaustive search for C1's two-full-sink cross block Sigma (see twosink_search.py for the conditions).
At each node (k chosen columns, in sister order), all 1260 candidate columns are tested at once:
  - row counts (two +1 and two -1 per row at the end), dichotomy (d) for sister columns;
  - mod 7: M = M_D on the chosen columns (symmetric, rank r <= 2).  For a candidate, b = new off-diagonal column, c = new
    diagonal.  Extended rank = r + 2 if b not in col(M); else r + [c != x^T M x] where M x = b.  Need <= 2.  Plus the anisotropy
    rule: a zero diagonal entry forces a zero row (M_D = -K r K^T, r anisotropic);
  - real: G = G_D on the chosen columns (PSD, rank <= 4).  For a candidate, g = new column, gamma = new diagonal:
    need g in col(G), gamma - g^T G^+ g >= 0, and rank(G) + [gamma - g^T G^+ g > 0] <= 4.
Conditions on principal blocks are necessary, so pruning is sound.  Column 0 is taken up to the stabiliser of label 0 in the
signed permutations of the 5 blocks (acting on rows and columns simultaneously), as in twosink_search.py.
At a leaf (10 columns): the same checks for Sigma^T (condition (c)) and (d) for sister rows, done by the slow exact routine."""
import itertools, sys, time, json
import numpy as np
from twosink_search import (COLS, A0, Wm, Mcp, BASE, n, N12, rm_cols, w_cols, rank_mod, canon_col0, check_block)

P7 = 7
INV7 = {1: 1, 2: 4, 3: 5, 4: 2, 5: 3, 6: 6}


def mod7_basis(M):
    """row-reduce symmetric M (k x k) mod 7: return (r, U, left) with U = k x r basis of col(M), left = (k-r) x k with
    left @ M = 0 (mod 7), and a solver matrix Sx (r x k) with M @ (U-coordinates)... we simply return the RREF data."""
    k = M.shape[0]
    A = [list(map(int, row)) for row in (M % P7)]
    # compute basis of column space and left null space via elimination on [M | I]
    aug = [A[i] + [1 if j == i else 0 for j in range(k)] for i in range(k)]
    r = 0
    pivcols = []
    for c in range(k):
        piv = next((i for i in range(r, k) if aug[i][c] % P7), None)
        if piv is None:
            continue
        aug[r], aug[piv] = aug[piv], aug[r]
        inv = INV7[aug[r][c] % P7]
        aug[r] = [(x * inv) % P7 for x in aug[r]]
        for i in range(k):
            if i != r and aug[i][c] % P7:
                f = aug[i][c]
                aug[i] = [(a - f * b) % P7 for a, b in zip(aug[i], aug[r])]
        pivcols.append(c)
        r += 1
    left = np.array([row[k:] for row in aug[r:]], dtype=np.int64).reshape(k - r, k)    # left @ M = 0
    return r, pivcols, left


def solve_in_col(M, r, pivcols, B):
    """for each column b of B (k x m) known to lie in col(M) (M symmetric), return x (k x m) with M x = b (mod 7), and
    the quadratic form value x^T M x (mod 7), computed via the r x r pivot block."""
    k = M.shape[0]
    m = B.shape[1]
    if r == 0:
        return np.zeros((k, m), dtype=np.int64)
    Msub = M[np.ix_(pivcols, pivcols)] % P7          # invertible r x r (symmetric rank r, pivots of rows = cols)
    # inverse mod 7 of Msub
    A = [list(map(int, row)) + [1 if j == i else 0 for j in range(r)] for i, row in enumerate(Msub)]
    for c in range(r):
        piv = next(i for i in range(c, r) if A[i][c] % P7)
        A[c], A[piv] = A[piv], A[c]
        inv = INV7[A[c][c] % P7]
        A[c] = [(x * inv) % P7 for x in A[c]]
        for i in range(r):
            if i != c and A[i][c] % P7:
                f = A[i][c]
                A[i] = [(a - f * b) % P7 for a, b in zip(A[i], A[c])]
    Minv = np.array([row[r:] for row in A], dtype=np.int64)
    X = np.zeros((k, m), dtype=np.int64)
    X[pivcols, :] = (Minv @ (B[pivcols, :] % P7)) % P7
    return X


def run(tlim=3600, report=60):
    t0 = time.time()
    last = [t0]
    col0 = canon_col0()
    stats = {}
    found = []
    nodes = [0]
    timed = [False]
    Rall = {y: rm_cols(COLS, y) for y in range(n)}                 # 1260 x 10 per position y
    Wall = {y: w_cols(COLS, y) for y in range(n)}                  # 1260 x 12
    diagM = {y: (A0[y, y] + np.einsum('ci,ij,cj->c', Rall[y], Wm, Rall[y])) % P7 for y in range(n)}
    diagG = {y: (BASE[y, y] - np.einsum('ci,ij,cj->c', Wall[y], Mcp, Wall[y])) / 4 for y in range(n)}
    colidx = {tuple(v): i for i, v in enumerate(COLS)}

    def bump(key, c):
        stats[key] = stats.get(key, 0) + int(c)

    def dfs(chosen_idx, M, G, rowp, rowm):
        if time.time() - t0 > tlim:
            timed[0] = True
            return
        nodes[0] += 1
        if time.time() - last[0] > report:
            last[0] = time.time()
            print('  .. %.0fs nodes %d depth %d stats %s' % (time.time() - t0, nodes[0], len(chosen_idx), stats), flush=True)
        k = len(chosen_idx)
        y = k
        if k == n:
            S = COLS[chosen_idx].T                                   # rows w, columns y
            ST = S.T
            Rt = np.array([rm_cols(ST[:, yy][None, :], yy)[0] for yy in range(n)]).T
            Wt = np.array([w_cols(ST[:, yy][None, :], yy)[0] for yy in range(n)]).T
            if check_block(Rt, Wt) is None and all(sum(1 for w in range(n) if S[yy, w] and S[yy ^ 1, w]) <= 2
                                                   for yy in range(0, n, 2)):
                found.append(S.copy())
            else:
                bump('leaf_transpose', 1)
            return
        cand = np.array([colidx[tuple(v)] for v in col0]) if k == 0 else np.arange(len(COLS))
        V = COLS[cand]
        p = rowp[None, :] + (V == 1)
        m = rowm[None, :] + (V == -1)
        ok = ~((p > 2).any(1) | (m > 2).any(1) | (((2 - p) + (2 - m)) > n - k - 1).any(1))
        bump('rows', (~ok).sum())
        if y % 2 == 1:
            prev = COLS[chosen_idx[-1]]
            ov = ((V != 0) & (prev[None, :] != 0)).sum(1)
            bad = ok & (ov > 2)
            bump('dich', bad.sum())
            ok &= ~bad
        if not ok.any():
            return
        cand, V = cand[ok], V[ok]
        Ry = Rall[y][cand]                                            # c x 10
        Wy = Wall[y][cand]                                            # c x 12
        cM = diagM[y][cand]
        cG = diagG[y][cand]
        if k:
            Rch = np.array([Rall[yy][ci] for yy, ci in enumerate(chosen_idx)])     # k x 10
            Wch = np.array([Wall[yy][ci] for yy, ci in enumerate(chosen_idx)])     # k x 12
            B = (A0[:k, y][:, None] + Rch @ Wm @ Ry.T) % P7           # k x c
            Gb = (BASE[:k, y][:, None] - Wch @ Mcp @ Wy.T) / 4         # k x c
            # mod 7 rank test
            r, pivcols, left = mod7_basis(M)
            inspan = np.ones(len(cand), bool) if left.shape[0] == 0 else ((left @ B) % P7 == 0).all(0)
            newrank = np.where(inspan, r, r + 2)
            X = solve_in_col(M, r, pivcols, B)
            q = np.einsum('ic,ij,jc->c', X, M % P7, X) % P7
            newrank = np.where(inspan, r + (q != cM % P7), newrank)
            bad = newrank > 2
            # anisotropy: zero diagonal => zero row
            zd = np.array([M[i, i] % P7 == 0 for i in range(k)])
            if zd.any():
                bad |= (B[zd, :] % P7 != 0).any(0)
            bad |= (cM % P7 == 0) & (B % P7 != 0).any(0)
            bump('rank7', bad.sum())
            keep = ~bad
            if not keep.any():
                return
            cand, V, B, Gb, cM, cG = cand[keep], V[keep], B[:, keep], Gb[:, keep], cM[keep], cG[keep]
            # real PSD / rank test
            Gp = np.linalg.pinv(G)
            ev = np.linalg.eigvalsh(G)
            rG = int((ev > 1e-7).sum())
            resid = Gb - G @ (Gp @ Gb)
            inrange = np.abs(resid).max(0) < 1e-6
            schur = cG - np.einsum('ic,ij,jc->c', Gb, Gp, Gb)
            badr = (~inrange) | (schur < -1e-7) | ((rG + (schur > 1e-7)) > 4)
            bump('real', badr.sum())
            keep = ~badr
            cand, V, B, Gb, cM, cG = cand[keep], V[keep], B[:, keep], Gb[:, keep], cM[keep], cG[keep]
        else:
            bad = (cG < -1e-7)
            bump('real', bad.sum())
            keep = ~bad
            cand, V, cM, cG = cand[keep], V[keep], cM[keep], cG[keep]
            B = np.zeros((0, len(cand)), dtype=np.int64)
            Gb = np.zeros((0, len(cand)))
        for t in range(len(cand)):
            M2 = np.zeros((k + 1, k + 1), dtype=np.int64)
            G2 = np.zeros((k + 1, k + 1))
            if k:
                M2[:k, :k] = M
                G2[:k, :k] = G
                M2[:k, k] = M2[k, :k] = B[:, t]
                G2[:k, k] = G2[k, :k] = Gb[:, t]
            M2[k, k] = cM[t]
            G2[k, k] = cG[t]
            v = COLS[cand[t]]
            dfs(chosen_idx + [cand[t]], M2 % P7, G2, rowp + (v == 1), rowm + (v == -1))
            if timed[0]:
                return
    dfs([], np.zeros((0, 0), dtype=np.int64), np.zeros((0, 0)), np.zeros(n, int), np.zeros(n, int))
    return {'complete': not timed[0], 'nodes': nodes[0], 'stats': stats, 'found': found, 'seconds': time.time() - t0}


if __name__ == '__main__':
    tl = float(sys.argv[1]) if len(sys.argv) > 1 else 280
    res = run(tl)
    print('COMPLETE' if res['complete'] else 'TIME LIMIT (incomplete)', '| nodes', res['nodes'], '| stats', res['stats'],
          '| survivors', len(res['found']), '| %.0fs' % res['seconds'], flush=True)
    if res['found']:
        json.dump([S.tolist() for S in res['found']], open('twosink_survivors.json', 'w'))
