"""Quick test for a task proposal to Astra: on real srg(v,k,1,2) graphs, compare
   dim W, W = {edge weights w : [A, L_w] = 0}   (L_w = weighted Laplacian of G)
with the number of components of K (edges of G, adjacent when opposite in a square).
Astra's component indicators lie in W; if dim W is larger, the trace-integrality argument has more room."""
import itertools
import sys
import numpy as np
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
import graphs


def analyse(name, A):
    A = np.array(A, dtype=np.int64)
    n = A.shape[0]
    E = [(u, v) for u in range(n) for v in range(u + 1, n) if A[u, v]]
    m = len(E)
    eid = {e: t for t, e in enumerate(E)}
    # K components (union-find over square-opposite edge pairs)
    par = list(range(m))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    nbr = [set(np.flatnonzero(A[u])) for u in range(n)]
    for (u, v) in E:
        for w in nbr[v] - {u}:
            if w in nbr[u]:
                continue
            for x in (nbr[u] & nbr[w]) - {v}:
                a, b = eid[(min(u, v), max(u, v))], eid[(min(w, x), max(w, x))]
                par[find(a)] = find(b)
    comps = len({find(t) for t in range(m)})
    # W: columns vec([A, L_e]); Gram matrix M^T M, rank by eigenvalues
    cols = []
    for (u, v) in E:
        d = np.zeros(n); d[u], d[v] = 1, -1
        Le = np.outer(d, d)
        C = A @ Le - Le @ A
        cols.append(C.ravel())
    M = np.array(cols).T
    G = M.T @ M
    ev = np.linalg.eigvalsh(G.astype(float))
    rank = int((ev > 1e-6 * ev.max()).sum())
    print(f'{name}: {n} vertices, {m} edges; K components: {comps}; dim W = {m - rank}')


analyse('Paley(9) = srg(9,4,1,2)', graphs.srg9())
analyse('BvLS = srg(243,22,1,2)', graphs.bvls())
