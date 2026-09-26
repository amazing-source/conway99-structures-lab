"""Independent check (own code, exact integers, no solver library) of Astra's dossier-5 report, section 1:
an explicit line ball G (36 vertices, cubic) with d(L) = 1 and tau(G) = 6 mod 7 (so L23 does not apply) admits no
k = 22 completion, because seven far vertices y_v with A-shadow {A_0, A_v}, v in {1..6, 8}, cannot be chosen jointly.
Also checks the positive control (v in {1..7}) and the algebraic counterexample of section 3.

Own implementation choices (different from Astra's annex):
  - Green function mod 7 obtained by solving H X = I - J with Gauss-Jordan mod 7 and normalising 1^T X = 0
    (not by inverting H + J);
  - tau(G) computed exactly over the integers (Bareiss) from the reduced Laplacian;
  - search = DFS over groups ordered by current domain size, with bitset (Python int) domains."""
import itertools

P7 = 7
a_inv = [9, 4, 7, 8, 1, 6, 5, 2, 3, 0, 11, 10]
b_inv = [7, 2, 1, 9, 11, 8, 10, 0, 5, 3, 6, 4]
c_inv = [7, 3, 10, 1, 11, 9, 8, 0, 6, 5, 2, 4]
N = 36
part = lambda p: p // 12
adj = [set() for _ in range(N)]
for t, f in enumerate([a_inv, b_inv, c_inv]):
    assert all(f[f[i]] == i and f[i] != i for i in range(12)), 'not a fixed-point-free involution'
    for i in range(12):
        adj[12 * t + i].add(12 * t + f[i])
for i in range(12):
    for p, q in [(i, 12 + i), (i, 24 + i), (12 + i, 24 + (i + 1) % 12)]:
        adj[p].add(q); adj[q].add(p)
partner = {p: next(q for q in adj[p] if part(q) == part(p)) for p in range(N)}

# --- structural checks
cubic = all(len(adj[p]) == 3 for p in range(N))
triangles = sum(1 for p in range(N) for q in adj[p] for r in adj[q] if r in adj[p] and p < q < r)
def is_single_cycle(vs):
    vs = set(vs)
    if any(len(adj[p] & vs) != 2 for p in vs):
        return False
    start = min(vs); seen = {start}; stack = [start]
    while stack:
        p = stack.pop()
        for q in adj[p] & vs:
            if q not in seen:
                seen.add(q); stack.append(q)
    return seen == vs
c24 = [is_single_cycle([p for p in range(N) if part(p) in pr]) for pr in [(0, 1), (0, 2), (1, 2)]]
Bm = [[int(q in adj[p]) for q in range(N)] for p in range(N)]

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

H = [[3 * (p == q) - Bm[p][q] for q in range(N)] for p in range(N)]
rank_H = rank_mod(H)

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
tau = bareiss_det([row[1:] for row in H[1:]])                   # Laplacian of a cubic graph = 3I - B
print(f'cubic {cubic}; triangles {triangles}; two-part unions are single C24: {c24}; rank_7(3I - B) = {rank_H}')
print(f'tau(G) = {tau}; tau mod 7 = {tau % 7}; squares mod 7 are {{1, 2, 4}}')

# --- Green function mod 7: solve H X = I - J, then force 1^T X = 0
def solve_mod(M, Rhs, p=P7):
    n = len(M); m = len(Rhs[0])
    Aug = [[x % p for x in M[i]] + [y % p for y in Rhs[i]] for i in range(n)]
    pivcols = []; r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if Aug[i][c]), None)
        if piv is None:
            continue
        Aug[r], Aug[piv] = Aug[piv], Aug[r]
        inv = pow(Aug[r][c], -1, p)
        Aug[r] = [x * inv % p for x in Aug[r]]
        for i in range(n):
            if i != r and Aug[i][c]:
                f = Aug[i][c]
                Aug[i] = [(x - f * y) % p for x, y in zip(Aug[i], Aug[r])]
        pivcols.append(c); r += 1
    assert all(not any(Aug[i][n:]) for i in range(r, n)), 'inconsistent system'
    X = [[0] * m for _ in range(n)]
    for i, c in enumerate(pivcols):
        X[c] = Aug[i][n:]
    return X
IJ = [[(p == q) - 1 for q in range(N)] for p in range(N)]
X0 = solve_mod(H, IJ)
inv36 = pow(36, -1, P7)
colsum = [sum(X0[p][q] for p in range(N)) % P7 for q in range(N)]
Lg = [[(X0[p][q] - colsum[q] * inv36) % P7 for q in range(N)] for p in range(N)]
assert all(Lg[p][q] == Lg[q][p] for p in range(N) for q in range(N)), 'Green not symmetric'
assert all(sum(Lg[p]) % P7 == 0 for p in range(N)), 'Green 1 != 0'
HL = [[sum(H[p][r] * Lg[r][q] for r in range(N)) % P7 for q in range(N)] for p in range(N)]
assert all(HL[p][q] == ((p == q) - 1) % P7 for p in range(N) for q in range(N)), 'H L != I - J'

# --- K = M M^T predicted by the ball: 12I - B + 2J - B^2 - D (D = same-part indicator)
B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
K = [[12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) for q in range(N)] for p in range(N)]
assert all(K[p][q] >= 0 for p in range(N) for q in range(N)) and all(K[p][p] == 10 for p in range(N))
g = lambda S, T: sum(Lg[p][q] for p in S for q in T) % P7

u = 0
def nonpartner_pairs(t):
    return [(12 * t + i, 12 * t + j) for i in range(12) for j in range(i + 1, 12) if partner[12 * t + i] != 12 * t + j]
PB, PC = nonpartner_pairs(1), nonpartner_pairs(2)

def candidates(v):
    out = []
    for eb in PB:
        for ec in PC:
            S = (u, v) + eb + ec
            Sset = set(S)
            if any(len(adj[p] & Sset) > 2 - (p in Sset) for p in range(N)):
                continue                                      # lambda / mu against the ball
            if any(K[p][q] == 0 for p, q in itertools.combinations(S, 2)):
                continue                                      # no far common neighbour possible
            if g(S, S) != 6:
                continue                                      # R = 0 on the diagonal
            out.append(S)
    return out

def run(groups):
    cand = {v: candidates(v) for v in groups}
    counts = [len(cand[v]) for v in groups]
    cap = {v: [1 - len(adj[u] & set(S)) for S in cand[v]] for v in groups}   # 1 - (B n)_u
    # pairwise data: allowed (g in {2,3}, overlap <= 2 - a) and adjacency a = 3 - g
    ok, isadj = {}, {}
    hvec = {v: [[sum(Lg[p][q] for p in S) % P7 for q in range(N)] for S in cand[v]] for v in groups}
    sets = {v: [set(S) for S in cand[v]] for v in groups}
    for v, w in itertools.permutations(groups, 2):
        for x, S in enumerate(cand[v]):
            okb, adb = 0, 0
            hS = hvec[v][x]
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
    def dfs(chosen, doms):
        nodes[0] += 1
        if not doms:
            return dict(chosen)
        v = min(doms, key=lambda t: bin(doms[t]).count('1'))
        d = doms[v]
        while d:
            x = (d & -d).bit_length() - 1
            d &= d - 1
            # star capacity of the already chosen ones and of x
            deg_x = sum(1 for w, y in chosen.items() if (isadj[v, w, x] >> y) & 1)
            if deg_x > cap[v][x]:
                continue
            bad = False
            degs = {}
            for w, y in chosen.items():
                dy = sum(1 for w2, y2 in chosen.items() if w2 != w and (isadj[w, w2, y] >> y2) & 1)
                dy += (isadj[w, v, y] >> x) & 1
                if dy > cap[w][y]:
                    bad = True; break
                degs[w] = dy
            if bad:
                continue
            new = dict(chosen); new[v] = x
            nd = {}
            for w, dw in doms.items():
                if w == v:
                    continue
                dw &= ok[v, w, x]
                for w2, y2 in new.items():
                    if w2 != v:
                        dw &= ok[w2, w, y2]
                # a chosen vertex at full capacity forbids further adjacency to it
                for w2, y2 in new.items():
                    full = (degs.get(w2) if w2 != v else deg_x) == cap[w2][y2]
                    if full:
                        dw &= ~isadj[w2, w, y2]
                # candidate's own capacity against the chosen ones
                keep = 0; t = dw
                while t:
                    z = (t & -t).bit_length() - 1
                    t &= t - 1
                    if sum(1 for w2, y2 in new.items() if (isadj[w, w2, z] >> y2) & 1) <= cap[w][z]:
                        keep |= 1 << z
                if not keep:
                    bad = True; break
                nd[w] = keep
            if bad:
                continue
            res = dfs(new, nd)
            if res is not None:
                return res
        return None
    sol = dfs({}, {v: (1 << len(cand[v])) - 1 for v in groups})
    if sol is not None:
        sol = {v: cand[v][x] for v, x in sol.items()}
    return counts, sol, nodes[0]

counts, sol, nodes = run([1, 2, 3, 4, 5, 6, 8])
print(f'groups 1..6,8: candidate counts {counts}; joint choice: {sol}; DFS nodes {nodes}')
counts7, sol7, nodes7 = run([1, 2, 3, 4, 5, 6, 7])
print(f'positive control, groups 1..7: candidate counts {counts7}; joint choice found: {sol7 is not None} ({nodes7} nodes)')
if sol7:
    for v in sorted(sol7):
        print('   ', v, ':', *sol7[v])
wit = {1: (0, 1, 15, 19, 31, 34), 2: (0, 2, 17, 18, 32, 35), 3: (0, 3, 15, 18, 26, 29), 4: (0, 4, 18, 20, 27, 33),
       5: (0, 5, 15, 20, 24, 35), 6: (0, 6, 20, 21, 25, 26), 7: (0, 7, 16, 22, 25, 33)}
cand7 = {v: set(candidates(v)) for v in wit}
in_pool = all(tuple(wit[v]) in cand7[v] for v in wit)
pair_ok = all(g(wit[v], wit[w]) in (2, 3) and len(set(wit[v]) & set(wit[w])) <= 2 - (3 - g(wit[v], wit[w]))
              for v, w in itertools.combinations(wit, 2))
star_ok = all(sum(3 - g(wit[v], wit[w]) for w in wit if w != v) <= 1 - len(adj[u] & set(wit[v])) for v in wit)
print(f"Astra's witness for 1..7: all shadows are candidates {in_pool}; pairwise OK {pair_ok}; star OK {star_ok}")

# --- section 3 counterexample over F_7
R4 = [[0, 1, 4, 2], [1, 0, 2, 4], [4, 2, 0, 1], [2, 4, 1, 0]]
uu = [1, 2, 4]
R = [[0] * 7 for _ in range(7)]
for i in range(4):
    for j in range(4):
        R[i][j] = R4[i][j]
for i in range(3):
    for j in range(3):
        R[4 + i][4 + j] = uu[i] * uu[j] % P7
R2 = [[sum(R[i][k] * R[k][j] for k in range(7)) % P7 for j in range(7)] for i in range(7)]
print('section 3 counterexample: rank', rank_mod(R), '; R1 = 0:', all(sum(r) % P7 == 0 for r in R),
      '; tr R = 0:', sum(R[i][i] for i in range(7)) % P7 == 0, '; diag(R^2) = 0:', all(R2[i][i] == 0 for i in range(7)),
      '; R^2 != 0:', any(any(r) for r in R2))
