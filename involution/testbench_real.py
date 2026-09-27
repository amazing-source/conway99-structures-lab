"""Test bench for external arguments on the involution case (structures 27/09; C1's control, notes/138 par. 3):
every step that does not use m = 7 must hold on the two EXISTING graphs srg(1 + 2m^2, 2m, 1, 2) with an involution
fixing one vertex: Paley 9 = H(2,3) (m = 2) and BvLS (m = 11), both Cayley on F_3^n with sigma: v -> -v.
Eigenvalues (-1 +- t)/2, t = sqrt(8m - 7): t = 3 (m = 2), 7 (m = 7), 9 (m = 11); the control prime is p = 3 for both
(analogue of N = S - 3I mod 7: N = S - I mod 3 at m = 2, N = S - 4I mod 3 at m = 11).
Objects (C2 notation, gauge = orbit representative containing the + vertex of its lower block):
  orbits: m interior (blocks) then m(m-1) exterior; S[U,W] = A[u,w] - A[u,sigma w] (odd part), Q[U,W] = A[u,w] + A[u,sigma w];
  D = D-matching of exterior orbits (Lemma 1); per block c: star(c) (2(m-1) orbits), far ((m-1)(m-2)), local type
  (mu_c, h_c, tau_c = pi h_c, pivots), ball H_c = -(S - rI) on interior + star(c) (= type_table.ball at m = 7),
  X_F = -(S - rI)[ball, far], T_FF = S[far, far].  Writes testbench_<name>.json.
Usage: python testbench_real.py"""
import itertools, json, sys, math
import numpy as np
sys.path.insert(0, '../famille')
from family_graphs import hamming_9, bvls_243, check_srg


def rank_mod(M, p):
    M = [[int(x) % p for x in row] for row in M]
    r = 0; n = len(M); cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(r, n) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(n):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r


def det_int(M):
    """exact integer determinant (Bareiss)"""
    M = [list(map(int, row)) for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k]), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def bench(name, V, adj):
    n = len(V)
    idx = {v: i for i, v in enumerate(V)}
    neg = [idx[tuple((-x) % 3 for x in v)] for v in V]
    A = np.zeros((n, n), dtype=int)
    for i, nb in enumerate(adj):
        for j in nb:
            A[i, j] = 1
    x0 = idx[tuple([0] * len(V[0]))]
    N0 = [u for u in range(n) if A[x0, u]]
    blocks, seen = [], set()
    for u in sorted(N0):
        if u not in seen:
            assert A[u, neg[u]] == 1 and neg[u] in N0
            blocks.append((u, neg[u])); seen |= {u, neg[u]}
    m = len(blocks); k = 2 * m
    t = int(round(math.sqrt(8 * m - 7))); assert t * t == 8 * m - 7
    r, s = (t - 1) // 2, (-1 - t) // 2
    assert check_srg([set(a) for a in adj], k, 1, 2)
    wv = {}                                        # world vertex (block, sign) -> graph vertex
    for i, (u, v) in enumerate(blocks):
        wv[(i, 0)], wv[(i, 1)] = u, v
    gw = {g: w for w, g in wv.items()}
    roots = [x for x in range(n) if x != x0 and not A[x0, x]]
    pair = {}
    for x in roots:
        nb = [u for u in N0 if A[x, u]]
        assert len(nb) == 2
        a, b = sorted((gw[nb[0]], gw[nb[1]]))
        assert a[0] != b[0]
        pair[x] = (a, b)
    # orbits: interior i -> rep (i,0); exterior {x, -x} -> rep containing (lower block, +)
    reps = [wv[(i, 0)] for i in range(m)]
    ext = sorted({min(x, neg[x], key=lambda y: (pair[y][0][1], pair[y])) for x in roots}, key=lambda y: pair[y])
    for x in ext:
        assert pair[x][0][1] == 0
    reps += ext
    no = len(reps); assert no == m * m
    S = np.array([[A[u, w] - A[u, neg[w]] for w in reps] for u in reps], dtype=int)
    Q = np.array([[A[u, w] + A[u, neg[w]] for w in reps] for u in reps], dtype=int)
    ok_S = (S @ S + S == (2 * m - 2) * np.eye(no, dtype=int)).all()
    oidx = {}
    for i, u in enumerate(reps):
        oidx[u] = i; oidx[neg[u]] = i
    D = {}
    for i, u in enumerate(reps[m:], start=m):
        com = [w for w in range(n) if A[u, w] and A[neg[u], w]]
        assert len(com) == 2 and com[1] == neg[com[0]] or com[0] == neg[com[1]]
        D[i] = oidx[com[0]]
    assert all(D[D[i]] == i for i in D)
    Np = (S - r * np.eye(no, dtype=int)) % 3
    rankN = rank_mod(Np.tolist(), 3)
    N2zero = ((Np @ Np) % 3 == 0).all()
    ones = np.ones(n, dtype=int)
    rankA = rank_mod(((A - r * np.eye(n, dtype=int)) % 3).tolist(), 3)
    out = {'name': name, 'm': m, 'k': k, 't': t, 'eig': [r, s], 'p': 3, 'n_orbits': no,
           'gauge': 'interior orbit i -> (i,+); exterior orbit -> representative containing the + vertex of its lower block',
           'orbits': [['int', i] if j < m else ['ext', [list(pair[u][0]), list(pair[u][1])]] for j, (i, u) in enumerate(zip(list(range(m)) + [None] * (no - m), reps))],
           'S': S.tolist(), 'Q': Q.tolist(), 'D': {str(i): D[i] for i in D},
           'checks': {'S^2 + S = (2m-2) I': bool(ok_S), 'N = S - rI mod 3: rank': rankN, 'N^2 = 0 mod 3': bool(N2zero),
                      'rank_3(A - rI)': rankA}}
    # D-square types (sib / one-block / dis) of the exterior orbits
    def blocks_of(i):
        (a, b) = pair[reps[i]]
        return {a[0], b[0]}
    dtype = {}
    for i, j in D.items():
        inter = len(blocks_of(i) & blocks_of(j))
        dtype[str(i)] = {2: 'sib', 1: 'one-block', 0: 'dis'}[inter]
    out['D_square_type'] = dtype
    out['D_square_counts'] = {x: sum(1 for v in dtype.values() if v == x) // 2 for x in ('sib', 'one-block', 'dis')}
    per = []
    for c in range(m):
        star = [i for i in range(m, no) if c in blocks_of(i)]
        far = [i for i in range(m, no) if c not in blocks_of(i)]
        ball = list(range(m)) + star
        Mr = -(S - r * np.eye(no, dtype=int))
        H = Mr[np.ix_(ball, ball)]
        labels = [(b, e) for b in range(m) if b != c for e in (0, 1)]
        L = len(labels)
        root = lambda sgn, lab: next(x for x in (reps[i] for i in range(m, no)) for y in (x, neg[x]) if False) if False else None
        cp, cm = wv[(c, 0)], wv[(c, 1)]
        def croot(base, lab):
            v = wv[lab]
            return next(x for x in roots if A[x, base] and A[x, v])
        mu = [next(j for j in range(L) if j != l and A[croot(cp, labels[l]), croot(cp, labels[j])]) for l in range(L)]
        h = [next(j for j in range(L) if A[croot(cp, labels[l]), croot(cm, labels[j])]) for l in range(L)]
        pi = [l ^ 1 for l in range(L)]
        tau = [pi[h[l]] for l in range(L)]
        piv = [l for l in range(L) if mu[l] == tau[l]]
        per.append({'c': c, 'labels': [list(x) for x in labels], 'mu': mu, 'tau': tau, 'pivots': piv,
                    'ball_index': ball, 'far_index': far, 'H_c': H.tolist(),
                    'rank_Q(H_c)': int(np.linalg.matrix_rank(H.astype(float))), 'rank_3(H_c)': rank_mod(H.tolist(), 3),
                    'det(H_c)': str(det_int(H.tolist())),
                    'X_F': Mr[np.ix_(ball, far)].tolist(), 'T_FF': S[np.ix_(far, far)].tolist()})
    out['blocks'] = per
    return out


if __name__ == '__main__':
    for name, (V, adj) in (('paley9', hamming_9()), ('bvls', bvls_243())):
        B = bench(name, V, adj)
        json.dump(B, open('testbench_%s.json' % name, 'w'))
        b0 = B['blocks'][0]
        print('%s: m = %d, t = %d, eigenvalues %s, orbits %d | %s | D-squares %s | block 0: ball %d x %d, rank_Q %d, rank_3 %d, '
              'det %s, pivots %d, far %d' % (name, B['m'], B['t'], B['eig'], B['n_orbits'], B['checks'], B['D_square_counts'],
                                            len(b0['ball_index']), len(b0['ball_index']), b0['rank_Q(H_c)'], b0['rank_3(H_c)'],
                                            b0['det(H_c)'], len(b0['pivots']), len(b0['far_index'])))
