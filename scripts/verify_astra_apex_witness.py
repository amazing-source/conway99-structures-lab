"""Independent check (own code) of Astra's apex witness (reply to follow-up 4, message 1 sec. 1), cyclic ball (2,4,1),
u = A_0: the 10 shadows pass the apex conditions, norms, exact star sum, overlaps and M M^T <= K, but the full Gram fails."""
import sys, itertools
sys.path.insert(0, 'scripts')
from star_shadows import ball_data, green_mod7, N, part, P7
from verify_astra_d5_cyclic import ball
adj = ball(2, 4, 1)
Bm, H, tau, rk = ball_data(adj)
Lg = green_mod7(H)
B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
K = [[12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) for q in range(N)] for p in range(N)]
partner = {p: next(q for q in adj[p] if part(q) == part(p)) for p in range(N)}
W = {1: [0,1,21,23,30,33], 2: [0,2,13,19,32,33], 3: [0,3,13,16,29,31], 4: [0,4,12,23,27,35], 5: [0,5,15,18,27,28],
     7: [0,7,15,22,25,26], 8: [0,8,14,19,31,34], 9: [0,9,14,17,32,35], 10: [0,10,20,22,28,29], 11: [0,11,16,21,24,26]}
u = 0
uB = next(q for q in adj[u] if part(q) == 1); uC = next(q for q in adj[u] if part(q) == 2)
groups = sorted(W)
shape = all(sum(1 for p in S if part(p) == t) == 2 and all(partner[p] not in S for p in S) for S in W.values() for t in range(3))
lm = all(len(adj[p] & set(S)) <= 2 - (p in S) for S in W.values() for p in range(N))
g = lambda S, T: sum(Lg[p][q] for p in S for q in T) % P7
norms = all(g(S, S) == 6 for S in W.values())
apex = [v for v in groups if (uB in W[v]) or (uC in W[v])]
caps = [1 - len(adj[u] & set(W[v])) for v in groups]
apex_rows = all(g(W[a], W[v]) == 3 for a in apex for v in groups if v != a)
star_sum = [0] * N
for S in W.values():
    for p in S:
        star_sum[p] += 1
exact_sum = star_sum == [K[p][u] for p in range(N)]
overl = all(len(set(W[v]) & set(W[w])) <= 2 for v, w in itertools.combinations(groups, 2))
MM = [[sum(1 for S in W.values() if p in S and q in S) for q in range(N)] for p in range(N)]
bound = all(MM[p][q] <= K[p][q] for p in range(N) for q in range(N))
bad = [(v, w, g(W[v], W[w])) for v, w in itertools.combinations(groups, 2) if g(W[v], W[w]) not in (2, 3)]
print(f'shape {shape}; lambda/mu vs ball {lm}; norms 6 {norms}; apexes {apex} (uB={uB}, uC={uC}); capacities {caps}')
print(f'apex rows g = 3 {apex_rows}; exact star sum = K e_u {exact_sum}; overlaps <= 2 {overl}; M M^T <= K {bound}')
print(f'off-diagonal Gram values outside {{2,3}}: {len(bad)} -> {bad}')
