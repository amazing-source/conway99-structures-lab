"""Positive control of C1's notes/125 identities on the two real members of the family:
sum_L c3(G(L)) = 2P and sum_L c4(G(L)) = S + 3P, with S = #{(p, E): E at distance 2 from p, E ~ Ebar}, where
Ebar = the vertex whose two neighbours in N(p) are the mates (line partners at p) of E's two neighbours in N(p)."""
import itertools
from family_graphs import hamming_9, bvls_243
from family_invariants import lines_of, holonomy

def run(name, adj, m):
    n = len(adj)
    L = lines_of(adj)
    line_of_edge = {}
    for T in L:
        for a, b in itertools.combinations(T, 2):
            line_of_edge[frozenset((a, b))] = T
    mate = lambda p, a: next(c for c in line_of_edge[frozenset((p, a))] if c not in (p, a))
    P = sum(sum(1 for q, s in holonomy(adj, T).items() if s == q) for T in L) // 2
    # S events
    S = 0
    for p in range(n):
        Np = adj[p]
        for E in range(n):
            if E == p or E in Np:
                continue
            a, b = sorted(adj[E] & Np)
            ab = adj[mate(p, a)] & adj[mate(p, b)] - {p}
            (Eb,) = ab
            S += Eb in adj[E]
    c3 = c4 = 0
    for T in L:
        NL = set().union(*(adj[t] for t in T)) - set(T)
        sub = {x: adj[x] & NL for x in NL}
        assert all(len(sub[x]) == 3 for x in NL)
        # triangles and 4-cycles of the cubic graph G(L)
        for x in NL:
            for y, z in itertools.combinations(sorted(sub[x]), 2):
                if z in sub[y] and x < y < z:
                    c3 += 1
        cyc = set()
        for x in NL:
            for y in sub[x]:
                for z in sub[y]:
                    if z == x:
                        continue
                    for w in sub[z]:
                        if w in (x, y):
                            continue
                        if x in sub[w]:
                            cyc.add(frozenset((x, y, z, w)) if len({x, y, z, w}) == 4 else None)
        cyc.discard(None)
        c4 += len(cyc)
    print(f'{name}: P = {P}, S = {S}; sum_L c3 = {c3} (2P = {2 * P}); sum_L c4 = {c4} (S + 3P = {S + 3 * P})')

run('H(2,3)', hamming_9()[1], 2)
run('BvLS', bvls_243()[1], 11)
