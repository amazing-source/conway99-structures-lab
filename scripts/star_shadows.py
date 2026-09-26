"""Reusable star-of-shadows test at k = 22 (own code of 'structures', exact integers, no solver library).
A line ball G: parts A, B, C of 12 vertices (A = 0..11, B = 12..23, C = 24..35), internal involutions a, b, c
(fixed-point-free), cross edges A_j - B_j, A_j - C_j, B_j - C_{hol[j]}.
For u = A_0 and a list of non-partner labels v, the far vertices y_v with A-shadow {u, v} must get shadows n_v
(two non-partner points per part) with:
  lambda/mu against the ball (Bn <= 2 - n), K = MM^T > 0 on the support, g(n,n) = 6 mod 7 (R = 0 on the diagonal),
  pairwise g in {2,3}, a = 3 - g, overlap <= 2 - a, and the star capacity sum_w a(n_v, n_w) <= 1 - (B n_v)_u
  (with equality when full=True and all 10 non-partners are used)."""
import itertools

P7 = 7
N = 36
part = lambda p: p // 12

def build(a_inv, b_inv, c_inv, hol):
    adj = [set() for _ in range(N)]
    for t, f in enumerate([a_inv, b_inv, c_inv]):
        assert all(f[f[i]] == i and f[i] != i for i in range(12))
        for i in range(12):
            adj[12 * t + i].add(12 * t + f[i])
    assert sorted(hol) == list(range(12))
    for i in range(12):
        for p, q in [(i, 12 + i), (i, 24 + i), (12 + i, 24 + hol[i])]:
            adj[p].add(q); adj[q].add(p)
    assert all(len(s) == 3 for s in adj)
    return adj

def rank_mod(M, p=P7):
    M = [[x % p for x in row] for row in M]
    r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r

def bareiss_det(M):
    M = [row[:] for row in M]
    n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def ball_data(adj):
    Bm = [[int(q in adj[p]) for q in range(N)] for p in range(N)]
    H = [[3 * (p == q) - Bm[p][q] for q in range(N)] for p in range(N)]
    tau = bareiss_det([row[1:] for row in H[1:]])
    rk = rank_mod(H)
    return Bm, H, tau, rk

def green_mod7(H):
    """Unique symmetric X with H X = I - J and X 1 = 0 (mod 7): Gauss-Jordan on H X = I - J, then 1^T X = 0."""
    n = N
    Rhs = [[(p == q) - 1 for q in range(n)] for p in range(n)]
    Aug = [[x % P7 for x in H[i]] + [y % P7 for y in Rhs[i]] for i in range(n)]
    pivcols = []; r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if Aug[i][c]), None)
        if piv is None:
            continue
        Aug[r], Aug[piv] = Aug[piv], Aug[r]
        inv = pow(Aug[r][c], -1, P7)
        Aug[r] = [x * inv % P7 for x in Aug[r]]
        for i in range(n):
            if i != r and Aug[i][c]:
                f = Aug[i][c]
                Aug[i] = [(x - f * y) % P7 for x, y in zip(Aug[i], Aug[r])]
        pivcols.append(c); r += 1
    assert all(not any(Aug[i][n:]) for i in range(r, n))
    X = [[0] * n for _ in range(n)]
    for i, c in enumerate(pivcols):
        X[c] = Aug[i][n:]
    inv36 = pow(36, -1, P7)
    colsum = [sum(X[p][q] for p in range(n)) % P7 for q in range(n)]
    L = [[(X[p][q] - colsum[q] * inv36) % P7 for q in range(n)] for p in range(n)]
    assert all(L[p][q] == L[q][p] for p in range(n) for q in range(n))
    assert all(sum(L[p]) % P7 == 0 for p in range(n))
    return L

def star_test(adj, groups, u=0, full=False):
    Bm, H, tau, rk = ball_data(adj)
    assert rk == 35, 'not a generic line (rank_7 H != 35)'
    Lg = green_mod7(H)
    B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
    K = [[12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) for q in range(N)] for p in range(N)]
    assert all(K[p][q] >= 0 for p in range(N) for q in range(N)) and all(K[p][p] == 10 for p in range(N))
    partner = {p: next(q for q in adj[p] if part(q) == part(p)) for p in range(N)}
    def pairs(t):
        return [(12 * t + i, 12 * t + j) for i in range(12) for j in range(i + 1, 12)
                if partner[12 * t + i] != 12 * t + j]
    PB, PC = pairs(1), pairs(2)
    cand = {}
    for v in groups:
        assert part(v) == 0 and v != u and v != partner[u]
        out = []
        for eb in PB:
            for ec in PC:
                S = (u, v) + eb + ec
                Ss = set(S)
                if any(len(adj[p] & Ss) > 2 - (p in Ss) for p in range(N)):
                    continue
                if any(K[p][q] == 0 for p, q in itertools.combinations(S, 2)):
                    continue
                if sum(Lg[p][q] for p in S for q in S) % P7 != 6:
                    continue
                out.append(S)
        cand[v] = out
    counts = [len(cand[v]) for v in groups]
    cap = {v: [1 - len(adj[u] & set(S)) for S in cand[v]] for v in groups}
    hvec = {v: [[sum(Lg[p][q] for p in S) % P7 for q in range(N)] for S in cand[v]] for v in groups}
    sets = {v: [set(S) for S in cand[v]] for v in groups}
    ok, isadj = {}, {}
    for v, w in itertools.permutations(groups, 2):
        for x in range(len(cand[v])):
            hS = hvec[v][x]; okb = adb = 0
            for y, T in enumerate(cand[w]):
                gg = sum(hS[q] for q in T) % P7
                if gg not in (2, 3):
                    continue
                aa = 3 - gg
                if len(sets[v][x] & sets[w][y]) > 2 - aa:
                    continue
                okb |= 1 << y
                if aa:
                    adb |= 1 << y
            ok[v, w, x], isadj[v, w, x] = okb, adb
    nodes = [0]
    def deg(w, y, chosen):
        return sum(1 for w2, y2 in chosen.items() if w2 != w and (isadj[w, w2, y] >> y2) & 1)
    def dfs(chosen, doms):
        nodes[0] += 1
        if not doms:
            if full and any(deg(w, y, chosen) != cap[w][y] for w, y in chosen.items()):
                return None
            return dict(chosen)
        v = min(doms, key=lambda t: bin(doms[t]).count('1'))
        d = doms[v]
        while d:
            x = (d & -d).bit_length() - 1
            d &= d - 1
            new = dict(chosen); new[v] = x
            if any(deg(w, y, new) > cap[w][y] for w, y in new.items()):
                continue
            nd = {}; dead = False
            for w, dw in doms.items():
                if w == v:
                    continue
                for w2, y2 in new.items():
                    dw &= ok[w2, w, y2]
                    if deg(w2, y2, new) == cap[w2][y2]:
                        dw &= ~isadj[w2, w, y2]
                keep = 0; t = dw
                while t:
                    z = (t & -t).bit_length() - 1
                    t &= t - 1
                    if sum(1 for w2, y2 in new.items() if (isadj[w, w2, z] >> y2) & 1) <= cap[w][z]:
                        keep |= 1 << z
                if not keep:
                    dead = True; break
                nd[w] = keep
            if dead:
                continue
            res = dfs(new, nd)
            if res is not None:
                return res
        return None
    sol = dfs({}, {v: (1 << len(cand[v])) - 1 for v in groups})
    if sol is not None:
        sol = {v: cand[v][x] for v, x in sol.items()}
    return {'tau': tau, 'counts': counts, 'solution': sol, 'nodes': nodes[0]}
