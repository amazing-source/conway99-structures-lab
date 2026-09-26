"""Independent check of Astra's section 3: components of G(L) meeting each part in 4 vertices, no triangle (holonomy
without fixed point) and no prism through the line (no 4-cycle with two opposite internal edges)."""
import itertools
from collections import Counter
M = {'a': [(0, 1), (2, 3)], 'b': [(0, 2), (1, 3)], 'c': [(0, 3), (1, 2)]}
def matching_map(name):
    f = {}
    for x, y in M[name]:
        f[x], f[y] = y, x
    return f
def conj(pi, f):                       # pi f pi^-1 as a matching
    return frozenset(frozenset((pi[x], pi[f[x]])) for x in range(4))
def mset(f):
    return frozenset(frozenset((x, f[x])) for x in range(4))
derangements = [p for p in itertools.permutations(range(4)) if all(p[i] != i for i in range(4))]
assert len(derangements) == 9
def laplacian_tau(edges, n):
    Lm = [[0] * n for _ in range(n)]
    for x, y in edges:
        Lm[x][y] -= 1; Lm[y][x] -= 1; Lm[x][x] += 1; Lm[y][y] += 1
    M_ = [row[1:] for row in Lm[1:]]
    # Bareiss
    n1 = len(M_); sign = 1; prev = 1
    for k in range(n1 - 1):
        if M_[k][k] == 0:
            sw = next((i for i in range(k + 1, n1) if M_[i][k] != 0), None)
            if sw is None:
                return 0
            M_[k], M_[sw] = M_[sw], M_[k]; sign = -sign
        for i in range(k + 1, n1):
            for j in range(k + 1, n1):
                M_[i][j] = (M_[i][j] * M_[k][k] - M_[i][k] * M_[k][j]) // prev
        prev = M_[k][k]
    return sign * M_[n1 - 1][n1 - 1]
adm = []
for mA, mB, mC in itertools.product('abc', repeat=3):
    for pi in derangements:
        fA, fB, fC = matching_map(mA), matching_map(mB), matching_map(mC)
        if mset(fA) == mset(fB) or mset(fA) == mset(fC):
            continue                                   # prism with rung on A-B or A-C (cross maps = identity)
        if conj(pi, fB) == mset(fC):
            continue                                   # prism with rung on B-C
        edges = set()
        for x in range(4):
            edges.add(frozenset((x, fA[x]))); edges.add(frozenset((4 + x, 4 + fB[x]))); edges.add(frozenset((8 + x, 8 + fC[x])))
            edges.add(frozenset((x, 4 + x))); edges.add(frozenset((x, 8 + x))); edges.add(frozenset((4 + x, 8 + pi[x])))
        edges = [tuple(e) for e in edges]
        assert len(edges) == 18
        adm.append(((mA, mB, mC), pi, laplacian_tau(edges, 12)))
print('admissible normalised data:', len(adm), 'out of', 27 * 9)
print('tau values:', Counter(t for _, _, t in adm), '; mod 7:', Counter(t % 7 for _, _, t in adm))
try:
    import networkx as nx
    graphs = []
    for (mA, mB, mC), pi, t in adm:
        fA, fB, fC = matching_map(mA), matching_map(mB), matching_map(mC)
        G = nx.Graph()
        for x in range(4):
            G.add_edge(x, fA[x]); G.add_edge(4 + x, 4 + fB[x]); G.add_edge(8 + x, 8 + fC[x])
            G.add_edge(x, 4 + x); G.add_edge(x, 8 + x); G.add_edge(4 + x, 8 + pi[x])
        graphs.append((G, t))
    classes = []
    for G, t in graphs:
        for c in classes:
            if nx.is_isomorphic(G, c[0]):
                c[2] += 1; break
        else:
            classes.append([G, t, 1])
    print('isomorphism classes of the 12-vertex cubic graphs (ignoring the layers):', [(c[1], c[2]) for c in classes])
except ImportError:
    print('networkx not available')
