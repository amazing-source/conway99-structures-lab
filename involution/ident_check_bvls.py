"""Check C1's general single-block identity on the real graph BvLS (m = 11), for random label sets S and blocks p:
  a3 + 3 a4 + sum_{Q in S_p} (k_Q - 1)+ + #{Q not in S_p : k_Q = 0}
      = N_Q - |S_p| - (2m-4)|S| + Lambda(S) + C_p(S)        (m = 7: 20 - 10|S| + Lambda(S) + C_p(S))
with k_Q = number of places of the quadruple Phi(z) u Phi(sigma z) in S (with multiplicity),
Lambda(S) = sum_{b<d in S} ([b,d non-mu] + K3[b,d]) + sum_{b in S} K3[b,b]/2, C_p(S) = sum_{b in S} (K2[p+,b] + K2[p-,b])."""
import itertools, random, sys
import numpy as np
sys.path.insert(0, '../famille')
from family_graphs import bvls_243
V, adjl = bvls_243()
V = [tuple(v) for v in V]
index = {v: i for i, v in enumerate(V)}
n = len(V)
A = np.zeros((n, n), dtype=int)
for i, nb in enumerate(adjl):
    for j in nb:
        A[i, j] = 1
neg = [index[tuple((-x) % 3 for x in v)] for v in V]
x0 = index[tuple([0] * len(V[0]))]
N0 = [u for u in range(n) if A[x0, u]]
blocks, seen = [], set()
for u in N0:
    if u not in seen:
        blocks.append((u, neg[u])); seen |= {u, neg[u]}
roots = [r for r in range(n) if r != x0 and not A[x0, r]]
rnb = {r: [u for u in N0 if A[r, u]] for r in roots}
m = 11; NL = 2 * (m - 1); Lr = list(range(NL)); PI = [v ^ 1 for v in Lr]
rng = random.Random(5)
bad = tot = 0
for ci, (cp, cm) in enumerate(blocks[:4]):
    others = [b for j, b in enumerate(blocks) if j != ci]
    lab_of = {}
    for i, (up, um) in enumerate(others):
        lab_of[up] = 2 * i; lab_of[um] = 2 * i + 1
    croot = {lab_of[b if a == cp else a]: r for r in roots for a, b in [rnb[r]] if cp in (a, b)}
    cmroot = {lab_of[b if a == cm else a]: r for r in roots for a, b in [rnb[r]] if cm in (a, b)}
    mu = [next(w for w in Lr if w != v and A[croot[v], croot[w]]) for v in Lr]
    h = [next(u for u in Lr if A[croot[v], cmroot[u]]) for v in Lr]
    tau = [PI[h[v]] for v in Lr]
    far = [r for r in roots if cp not in rnb[r] and cm not in rnb[r]]
    lab = {r: tuple(sorted(lab_of[u] for u in rnb[r])) for r in far}
    rof = {v: r for r, v in lab.items()}
    Phi = {lab[r]: [v for v in Lr if A[r, croot[v]]] for r in far}
    sig = {z: tuple(sorted((PI[z[0]], PI[z[1]]))) for z in Phi}
    K2 = np.full((NL, NL), 2); K3 = np.full((NL, NL), 2)
    for a in Lr:
        K2[a, a] -= 1; K2[a, PI[a]] -= 1; K2[a, mu[a]] -= 1; K2[a, tau[PI[a]]] -= 1
    for b in Lr:
        K3[b, PI[b]] -= 1; K3[b, tau[b]] -= 1; K3[b, tau[mu[b]]] -= 1; K3[b, mu[tau[b]]] -= 1
    Qs = sorted({tuple(sorted((z, sig[z]))) for z in Phi})
    NQ = len(Qs); assert NQ == 90
    for trial in range(300):
        size = rng.randint(1, 8)
        S = set(rng.sample(Lr, size))
        p = rng.randrange(m - 1)
        lhs = 0
        for (z1, z2) in Qs:
            k = sum(1 for b in Phi[z1] if b in S) + sum(1 for b in Phi[z2] if b in S)
            through = p in {a // 2 for a in z1}
            lhs += (k == 3) + 3 * (k == 4)
            lhs += max(k - 1, 0) if through else (k == 0)
        Lam = sum((mu[b] != d) + K3[b, d] for b, d in itertools.combinations(sorted(S), 2)) + sum(K3[b, b] // 2 for b in S)
        Cp = sum(K2[2 * p, b] + K2[2 * p + 1, b] for b in S)
        Sp = 2 * (m - 2)
        rhs = NQ - Sp - (2 * m - 4) * len(S) + Lam + Cp
        tot += 1
        if lhs != rhs:
            bad += 1
            if bad <= 3:
                print('MISMATCH block', ci, 'S', sorted(S), 'p', p, 'lhs', lhs, 'rhs', rhs)
print('BvLS general identity: %d checks (4 blocks x 300 random (S, p)), %d mismatches' % (tot, bad))
