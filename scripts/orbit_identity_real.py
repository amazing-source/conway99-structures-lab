"""Exact check of the ORBIT-LEVEL second-order identity (S19) on a real graph with the same frame as C2's:
BvLS = srg(243,22,1,2), a Cayley graph on F_3^5, vertex x = 0, involution sigma(v) = -v (fixes only 0, swaps the two
labels s, -s of every block).  With q_U(W) = number of neighbours of a vertex of U in W (U, W far orbits):
    Q^2 + Q = (k - 2) I + 4 J - Sh,   Sh(U,W) = number of common blocks of the block pairs of U and W.
Also checks: q_U(U) = 0; exactly one entry 2 per row (the D-partner); every other entry <= 1 (a far double hits both
orbits of its block pair)."""
import sys
import numpy as np
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
import graphs

A, V, S = graphs.bvls(return_coords=True)
A = np.array(A, dtype=np.int64)
idx = {v: t for t, v in enumerate(V)}
neg = lambda v: tuple((-np.array(v)) % 3)
x = idx[(0, 0, 0, 0, 0)]
k = int(A[x].sum())
Nx = [u for u in range(len(V)) if A[x, u]]
far = [u for u in range(len(V)) if u != x and not A[x, u]]
block = {}
for u in Nx:
    (m,) = [w for w in Nx if w != u and A[u, w]]
    block[u] = frozenset((u, m))
    assert idx[neg(V[u])] == m                      # sigma swaps the two labels of every block
blocks_of = {w: frozenset(block[u] for u in Nx if A[w, u]) for w in far}
assert all(len(b) == 2 for b in blocks_of.values())
orbits, seen = [], set()
for w in far:
    if w in seen:
        continue
    sw = idx[neg(V[w])]
    assert sw != w and sw in blocks_of
    orbits.append((w, sw)); seen |= {w, sw}
n = len(orbits)
Q = np.zeros((n, n), dtype=np.int64)
for a, (u, su) in enumerate(orbits):
    for b, (w, sw) in enumerate(orbits):
        Q[a, b] = A[u, w] + A[u, sw]
        assert Q[a, b] == A[su, w] + A[su, sw]      # the count does not depend on the vertex of U
Sh = np.array([[len(blocks_of[orbits[a][0]] & blocks_of[orbits[b][0]]) for b in range(n)] for a in range(n)])
lhs = Q @ Q + Q
rhs = (k - 2) * np.eye(n, dtype=np.int64) + 4 * np.ones((n, n), dtype=np.int64) - Sh
print(f'BvLS: {len(far)} far vertices, {n} orbits, far degree {k - 2}; Q symmetric: {(Q == Q.T).all()}; '
      f'Q^2 + Q = (k-2)I + 4J - Sh exactly: {(lhs == rhs).all()}')
print(f'q_U(U) = 0 for all U: {(np.diag(Q) == 0).all()}; exactly one 2 per row: {((Q == 2).sum(axis=1) == 1).all()}; '
      f'max entry {Q.max()}; row sums {set(Q.sum(axis=1).tolist())}')
ev = np.round(np.linalg.eigvalsh(Q.astype(float)), 6)
vals, cnts = np.unique(ev, return_counts=True)
print('spectrum of Q:', dict(zip(vals.tolist(), cnts.tolist())))
