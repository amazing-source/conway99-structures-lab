"""POSITIVE CONTROL of the local system on a real graph: BvLS = srg(243,22,1,2) (same lambda = 1, mu = 2; m = 11 blocks),
Cayley graph on F_3^5, involution sigma: v -> -v (fixes only 0).  At every block c of N(0), compute the true Phi, the mates
m (C2's mu), h and eta = pi o h (C2's tau), and check exactly:
  (i)   Phi: far roots -> non-m pairs is a bijection, |Phi(z)| = 2;
  (ii)  #{z containing a : b in Phi z} = K2[a,b] = 2 - [a=b] - [b=pi a] - [b=mu a] - [b=tau pi a];
  (iii) #{z : b in Phi z, d in Phi sigma z} = K3[b,d] = 2 - [d=pi b] - [d=tau b] - [d=tau mu b] - [d=mu tau b];
  (L1) (S1) (S2) pointwise (C1).
This validates the formulas used in every local model (fullblock_*, iota_design, sat_local, C1's g16/g18 code)."""
import itertools, sys
import numpy as np
sys.path.insert(0, '../famille')
from family_graphs import bvls_243

V, adj = bvls_243()
V = [tuple(v) for v in V]
n = len(V)
index = {v: i for i, v in enumerate(V)}
neg = [index[tuple((-x) % 3 for x in v)] for v in V]
A = np.zeros((n, n), dtype=int)
for _i, _nb in enumerate(adj):
    for _j in _nb:
        A[_i, _j] = 1
x0 = index[tuple([0] * len(V[0]))]
N0 = [u for u in range(n) if A[x0, u]]
assert len(N0) == 22
# blocks: {u, -u}
blocks, seen = [], set()
for u in N0:
    if u in seen:
        continue
    assert A[u, neg[u]] == 1 and neg[u] in N0
    blocks.append((u, neg[u])); seen |= {u, neg[u]}
assert len(blocks) == 11
roots = [r for r in range(n) if r != x0 and not A[x0, r]]
assert len(roots) == 220
rnb = {r: [u for u in N0 if A[r, u]] for r in roots}
assert all(len(v) == 2 for v in rnb.values())
bad_total = 0
for ci, (cp, cm) in enumerate(blocks):
    others = [b for j, b in enumerate(blocks) if j != ci]
    lab_of = {}
    for i, (up, um) in enumerate(others):
        lab_of[up] = 2 * i; lab_of[um] = 2 * i + 1
    Lr = list(range(20))
    PI = [v ^ 1 for v in Lr]
    vert = {v: u for u, v in lab_of.items()}
    # c+-roots (c+, v): the root adjacent to c+ and to vertex(v)
    croot = {}
    for r in roots:
        a, b = rnb[r]
        if cp in (a, b):
            other = b if a == cp else a
            croot[lab_of[other]] = r
    cmroot = {}
    for r in roots:
        a, b = rnb[r]
        if cm in (a, b):
            other = b if a == cm else a
            cmroot[lab_of[other]] = r
    assert len(croot) == 20 and len(cmroot) == 20
    mu = [next(w for w in Lr if w != v and A[croot[v], croot[w]]) for v in Lr]
    h = [next(u for u in Lr if A[croot[v], cmroot[u]]) for v in Lr]
    tau = [PI[h[v]] for v in Lr]                                     # C2's tau = pi o h (C1's eta)
    assert all(mu[mu[v]] == v for v in Lr) and all(tau[tau[v]] == v for v in Lr)
    far = [r for r in roots if cp not in rnb[r] and cm not in rnb[r]]
    assert len(far) == 180
    lab_root = {r: tuple(sorted(lab_of[u] for u in rnb[r])) for r in far}
    root_of = {v: r for r, v in lab_root.items()}
    Phi = {lab_root[r]: tuple(sorted(v for v in Lr if A[r, croot[v]])) for r in far}
    ok_i = all(len(e) == 2 for e in Phi.values()) and len(set(Phi.values())) == 180 and \
        all(mu[e[0]] != e[1] for e in Phi.values())
    sig = {z: tuple(sorted((PI[z[0]], PI[z[1]]))) for z in Phi}
    # check sigma really maps root z to the root with labels pi z
    ok_sig = all(neg[root_of[z]] == root_of[sig[z]] for z in Phi)
    K2 = np.full((20, 20), 2); K3 = np.full((20, 20), 2)
    for a in Lr:
        K2[a, a] -= 1; K2[a, PI[a]] -= 1; K2[a, mu[a]] -= 1; K2[a, tau[PI[a]]] -= 1
    for b in Lr:
        K3[b, PI[b]] -= 1; K3[b, tau[b]] -= 1; K3[b, tau[mu[b]]] -= 1; K3[b, mu[tau[b]]] -= 1
    M2 = np.zeros((20, 20), int); M3 = np.zeros((20, 20), int)
    for z, e in Phi.items():
        for a in z:
            for b in e:
                M2[a, b] += 1
        for b in e:
            for d in Phi[sig[z]]:
                M3[b, d] += 1
    ok_ii = (M2 == K2).all(); ok_iii = (M3 == K3).all()
    pw_bad = 0
    for z, e in Phi.items():
        f = Phi[sig[z]]
        for a in z:
            for aa in (a, PI[a]):
                if aa in e and PI[aa] in f:
                    pw_bad += 1
            l = a
            if l in e and tau[l] in f:
                pw_bad += 1
            if l not in e and mu[l] in e and tau[l] in f:
                pw_bad += 1
    full = all(mu[v] == tau[v] for v in Lr)
    bad = (not ok_i) + (not ok_sig) + (not ok_ii) + (not ok_iii) + (pw_bad > 0)
    bad_total += bad
    print('block %2d: (i) %s  sigma %s  (ii) %s  (iii) %s  pointwise violations %d  | full block: %s  | K2 min %d, K3 diag %s'
          % (ci, ok_i, ok_sig, ok_ii, ok_iii, pw_bad, full, K2.min(), sorted(set(np.diag(K3).tolist()))), flush=True)
print('ALL CHECKS PASSED' if bad_total == 0 else 'FAILURES: %d' % bad_total)
