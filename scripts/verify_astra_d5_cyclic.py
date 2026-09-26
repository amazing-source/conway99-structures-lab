"""Independent check of Astra's cyclic family (dossier 5, follow-up 3, section 1).
Vertices of each part: (i,t) in Z6 x {0,1}, local index i + 6t. Internal matchings (i,0) ~ (i+delta,1) with
(delta_A, delta_B, delta_C) = (0, r, s); cross edges A_{i,t}B_{i,t}, A_{i,t}C_{i,t}, B_{i,t}C_{i+q,t}."""
import itertools, sys, time
sys.path.insert(0, 'scripts')
from star_shadows import build, ball_data, star_test

def inv_from_delta(d):
    f = [None] * 12
    for i in range(6):
        j = (i + d) % 6 + 6
        f[i], f[j] = j, i
    return f

def ball(r, s, q):
    hol = [((i % 6 + q) % 6) + 6 * (i // 6) for i in range(12)]
    return build(inv_from_delta(0), inv_from_delta(r), inv_from_delta(s), hol)

SQ = {1, 2, 4}
survivors = []
tri_free_all = True; c4_free_all = True
for r, s, q in itertools.product(range(1, 6), range(1, 6), range(1, 6)):
    if r == s:
        continue
    adj = ball(r, s, q)
    tri = any(w in adj[p] for p in range(36) for v in adj[p] for w in adj[v] if w != p)
    tri_free_all &= not tri
    # 4-cycles inside the union of two parts
    for P1, P2 in [(0, 1), (0, 2), (1, 2)]:
        vs = {p for p in range(36) if p // 12 in (P1, P2)}
        for p in vs:
            for a_, b_ in itertools.combinations(adj[p] & vs, 2):
                if len((adj[a_] & adj[b_] & vs) - {p}) > 0:
                    c4_free_all = False
    Bm, H, tau, rk = ball_data(adj)
    t7 = tau % 7
    if rk == 35 and t7 != 0 and t7 not in SQ:
        survivors.append((r, s, q, tau, t7))
print(f'100 labelled triples: triangle-free {tri_free_all}; two-part unions without 4-cycles {c4_free_all}')
print('generic with tau a non-square mod 7 (needed at k = 22 by L23):', [(r, s, q, t7) for r, s, q, _, t7 in survivors])
for r, s, q, tau, t7 in survivors:
    if (r, s, q) in [(2, 4, 1), (2, 4, 3)]:
        print(f'   tau(2,4,{q}) = {tau} = {t7} mod 7')
# Astra's 7-group tests on the two representatives
for (r, s, q), groups in [((2, 4, 1), [1, 2, 3, 4, 5, 7, 11]), ((2, 4, 3), [1, 2, 3, 4, 5, 7, 9])]:
    t0 = time.time()
    res = star_test(ball(r, s, q), groups)
    print(f'({r},{s},{q}) groups {groups}: candidates {res["counts"]}; solution {res["solution"]}; nodes {res["nodes"]} ({time.time() - t0:.0f} s)')
# full stars (10 groups, equality) on ALL surviving triples, without relying on the coordinate reduction
for r, s, q, tau, t7 in survivors:
    t0 = time.time()
    groups = [v for v in range(1, 12) if v != 6]
    res = star_test(ball(r, s, q), groups, full=True)
    print(f'full star ({r},{s},{q}): candidates {res["counts"]}; solution {res["solution"] is not None}; nodes {res["nodes"]} ({time.time() - t0:.0f} s)')
