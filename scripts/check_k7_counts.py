"""Independent recount of two entries of k7_cycle_types.log by a different method.
(1) Fano (3^7): labelled Steiner triple systems on 7 points: 7!/168 = 30.
(2) (7,7,7): unordered triples of pairwise edge-disjoint Hamiltonian cycles of K7 (count ordered triples, divide by 6).
(3) (3,3,3,3,3,6): pick 5 edge-disjoint triangles whose complement is a single 6-cycle."""
import itertools
V = range(7)
def ham_cycles():
    out = set()
    for perm in itertools.permutations(range(1, 7)):
        cyc = (0,) + perm
        edges = frozenset(frozenset((cyc[t], cyc[(t + 1) % 7])) for t in range(7))
        out.add(edges)
    return list(out)
H = ham_cycles()
print('Hamiltonian cycles of K7:', len(H), '(expected 360)')
ordered = 0
for a in H:
    for b in H:
        if a & b: continue
        rest = frozenset(frozenset(e) for e in itertools.combinations(V, 2)) - a - b
        # rest has 7 edges, 2-regular; is it a single 7-cycle?
        adj = {v: [] for v in V}
        for e in rest:
            x, y = tuple(e); adj[x].append(y); adj[y].append(x)
        seen, stack = {0}, [0]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in seen: seen.add(w); stack.append(w)
        if len(seen) == 7: ordered += 1
print('(7,7,7) unordered decompositions:', ordered // 6, '(k7_cycle_types.log says 960; ordered pairs (a,b) with Hamiltonian complement:', ordered, ')')
tris = [frozenset(frozenset(e) for e in itertools.combinations(t, 2)) for t in itertools.combinations(V, 3)]
all_edges = frozenset(frozenset(e) for e in itertools.combinations(V, 2))
cnt5 = 0
for five in itertools.combinations(range(len(tris)), 5):
    es = [tris[i] for i in five]
    u = frozenset().union(*es)
    if len(u) != 15: continue
    rest = all_edges - u
    deg = {v: 0 for v in V}
    for e in rest:
        for x in e: deg[x] += 1
    verts = [v for v in V if deg[v]]
    if len(verts) != 6 or any(deg[v] != 2 for v in verts): continue
    adj = {v: [] for v in verts}
    for e in rest:
        x, y = tuple(e); adj[x].append(y); adj[y].append(x)
    seen, stack = {verts[0]}, [verts[0]]
    while stack:
        u2 = stack.pop()
        for w in adj[u2]:
            if w not in seen: seen.add(w); stack.append(w)
    if len(seen) == 6: cnt5 += 1
print('(3,3,3,3,3,6) decompositions:', cnt5, '(k7_cycle_types.log says 420)')
