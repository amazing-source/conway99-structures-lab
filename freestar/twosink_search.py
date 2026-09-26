"""Exhaustive search for C1's cross block Sigma (c1_recherche notes/131 sec. 6): two full-pivot blocks c, d forming a sib pair.
Sigma: 10 x 10 over the labels R (5 blocks x 2 signs; label i: block i//2, sign parity), entries 0/+-1, exactly two +1 and two
-1 in every row and every column.  Necessary conditions (C1's reduction; not re-derived here):
  (a) M_D = I - J + 2P + Rm^T (5I + 4P) Rm  (mod 7), Rm = Sigma + 2(I - P), has rank <= 2, and a zero diagonal entry forces a
      zero row (M_D = -K r K^T with r anisotropic);
  (b) G_D = (11I - J + P - W^T M_c^+ W)/4 is PSD of real rank <= 4;
  (c) the same for Sigma transposed;
  (d) for every block Y = {y, y^1}: at most 2 rows w with Sigma[w,y] != 0 and Sigma[w,y^1] != 0; same for sister rows.
Every condition restricted to a set of chosen columns is necessary for that principal block, so a column-by-column DFS
(columns in sister order 0,1,2,...) prunes soundly.  Symmetry: the conditions are invariant under simultaneous signed
permutations of the 5 blocks acting on rows and columns; column 0 is taken up to the stabilizer of label 0 (sound reduction:
every Sigma is equivalent to one whose column 0 is the lexicographically smallest in its orbit).
Output: COMPLETE with the list of survivors, or TIME LIMIT."""
import itertools, sys, time
import numpy as np

P7 = 7
n = 10
Ppi = np.zeros((n, n), int)
for i in range(n):
    Ppi[i, i ^ 1] = 1
I = np.eye(n, dtype=int)
A0 = (I - np.ones((n, n), int) + 2 * Ppi) % P7
Wm = (5 * I + 4 * Ppi) % P7
N12 = 12
P12 = np.zeros((N12, N12))
for i in range(N12):
    P12[i, i ^ 1] = 1
Mcp = np.linalg.pinv(11 * np.eye(N12) + P12 - np.ones((N12, N12)))
BASE = 11 * np.eye(n) - np.ones((n, n)) + Ppi.astype(float)

COLS = []
for plus in itertools.combinations(range(n), 2):
    rest = [i for i in range(n) if i not in plus]
    for minus in itertools.combinations(rest, 2):
        v = np.zeros(n, int)
        v[list(plus)] = 1
        v[list(minus)] = -1
        COLS.append(v)
COLS = np.array(COLS)                     # 1260 x 10


def rm_cols(V, y):
    R = V.copy()
    R[:, y] += 2
    R[:, y ^ 1] -= 2
    return R


def w_cols(V, y):
    W = np.zeros((V.shape[0], N12))
    W[:, 0], W[:, 1] = 1, -1
    W[:, 2:] = 4 * V
    W[:, 2 + y] += 1
    W[:, 2 + (y ^ 1)] -= 1
    return W


def rank_mod(M):
    M = [[int(x) % P7 for x in r] for r in M]
    rows = len(M)
    ncol = len(M[0]) if rows else 0
    r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], P7 - 2, P7)
        M[r] = [x * inv % P7 for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % P7 for a, b in zip(M[i], M[r])]
        r += 1
    return r


def signed_perms5():
    for perm in itertools.permutations(range(5)):
        for signs in itertools.product((0, 1), repeat=5):
            yield [2 * perm[i // 2] + ((i % 2) ^ signs[i // 2]) for i in range(n)]


G5 = list(signed_perms5())
STAB0 = [g for g in G5 if g[0] == 0]      # fixes label 0 (hence label 1)


def canon_col0():
    """one representative column vector per orbit of STAB0 acting on column 0's entries (rows permuted by g)."""
    reps = []
    seen = set()
    for v in COLS:
        best = None
        for g in STAB0:
            w = np.zeros(n, int)
            for i in range(n):
                w[g[i]] = v[i]
            t = tuple(w)
            if best is None or t < best:
                best = t
        if best not in seen:
            seen.add(best)
            reps.append(np.array(best))
    return reps


def check_block(Rm, W):
    """Rm: 10 x k (Rm columns), W: 12 x k.  Returns None if the leading principal conditions hold, else a reason."""
    k = Rm.shape[1]
    MD = (A0[:k, :k] + Rm.T @ Wm @ Rm) % P7
    if rank_mod(MD.tolist()) > 2:
        return 'rank7'
    for y in range(k):
        if MD[y, y] == 0 and MD[y, :].any():
            return 'aniso'
    G = (BASE[:k, :k] - W.T @ Mcp @ W) / 4
    ev = np.linalg.eigvalsh(G)
    if ev.min() < -1e-7:
        return 'psd'
    if (ev > 1e-7).sum() > 4:
        return 'rank4'
    return None


def search(tlim=280):
    t0 = time.time()
    stats = {}
    found = []
    col0 = canon_col0()
    nodes = [0]
    timed_out = [False]

    def dfs(chosen, Rm, W, rowp, rowm):
        if time.time() - t0 > tlim:
            timed_out[0] = True
            return
        nodes[0] += 1
        k = len(chosen)
        if k == n:
            S = np.array(chosen).T
            # (c): same conditions for Sigma^T, (d) for sister rows
            ST = S.T
            Rt = np.array([rm_cols(ST[:, y][None, :], y)[0] for y in range(n)]).T
            Wt = np.array([w_cols(ST[:, y][None, :], y)[0] for y in range(n)]).T
            if check_block(Rt, Wt) is None and all(sum(1 for w in range(n) if S[y, w] and S[y ^ 1, w]) <= 2
                                                   for y in range(0, n, 2)):
                found.append(S.copy())
            return
        y = k
        cand = col0 if k == 0 else COLS
        for v in cand:
            p = rowp + (v == 1)
            m = rowm + (v == -1)
            if (p > 2).any() or (m > 2).any():
                continue
            if ((2 - p) + (2 - m) > n - k - 1).any():
                continue
            if y % 2 == 1:                            # (d) for the block {y-1, y}
                prev = chosen[y - 1]
                if int(((prev != 0) & (v != 0)).sum()) > 2:
                    stats['dich'] = stats.get('dich', 0) + 1
                    continue
            rcol = rm_cols(v[None, :], y)[0]
            wcol = w_cols(v[None, :], y)[0]
            Rm2 = np.column_stack([Rm, rcol]) if k else rcol[:, None]
            W2 = np.column_stack([W, wcol]) if k else wcol[:, None]
            why = check_block(Rm2, W2)
            stats[why] = stats.get(why, 0) + 1
            if why is None:
                dfs(chosen + [v], Rm2, W2, p, m)
            if timed_out[0]:
                return
    dfs([], np.zeros((n, 0), int), np.zeros((N12, 0)), np.zeros(n, int), np.zeros(n, int))
    return {'complete': not timed_out[0], 'nodes': nodes[0], 'stats': stats, 'found': found, 'seconds': time.time() - t0,
            'col0_reps': len(col0)}


if __name__ == '__main__':
    tl = float(sys.argv[1]) if len(sys.argv) > 1 else 280
    res = search(tl)
    print('COMPLETE' if res['complete'] else 'TIME LIMIT (incomplete)', '| column-0 representatives', res['col0_reps'],
          '| nodes', res['nodes'], '| prune stats', res['stats'], '| survivors', len(res['found']),
          '| %.0fs' % res['seconds'])
    for S in res['found'][:3]:
        print(S)
