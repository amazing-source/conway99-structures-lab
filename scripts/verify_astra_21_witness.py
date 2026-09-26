"""Independent check (own code, exact, no solver) of Astra's local (2,1) witness for a dis square (dossier 4, §2).
Coding: interior (c,+) = 2c, (c,-) = 2c+1; sigma = XOR 1 on interiors; a root is a pair of interiors of different blocks.
Five representative roots with their 12 root neighbours are given; their sigma-images follow.  Checks:
  1. degree 14 of the 10 central roots, symmetry of the listed adjacencies;
  2. lambda/mu between the 10 central roots, and between them and x and the 14 interiors;
  3. no overflow of the common-neighbour budget on any pair of the partial graph (99 vertices);
  4. orbit rows: exactly one 2 per row, row sum 12; the (2,1) profile u -> {w, sigma w, w'} but not sigma w';
  5. orbit identity Q^2 + Q = 12I + 4J - Sh on the 5 x 5 block of the central orbits (computable from their rows)."""
import itertools

blk = lambda i: i // 2
sig = lambda r: frozenset(i ^ 1 for i in r)
R = lambda a, b: frozenset((a, b))
data = {
    R(0, 2): [(4, 6), (5, 7), (4, 7), (3, 10), (0, 11), (8, 13), (9, 10), (2, 12), (8, 12), (9, 11), (1, 5), (6, 13)],
    R(4, 6): [(0, 2), (1, 3), (9, 11), (2, 11), (1, 13), (4, 12), (7, 8), (10, 13), (5, 9), (10, 12), (6, 8), (0, 3)],
    R(4, 7): [(0, 2), (9, 10), (8, 11), (7, 12), (4, 13), (2, 13), (6, 12), (1, 8), (3, 11), (1, 9), (3, 5), (0, 10)],
    R(9, 11): [(0, 2), (0, 8), (1, 4), (1, 12), (2, 7), (3, 5), (3, 12), (4, 6), (5, 13), (6, 10), (7, 11), (9, 13)],
    R(2, 11): [(0, 7), (0, 13), (1, 3), (1, 8), (2, 9), (4, 6), (4, 9), (5, 8), (5, 12), (6, 12), (7, 10), (11, 13)],
}
ROOTS = [frozenset((a, b)) for a in range(14) for b in range(a + 1, 14) if blk(a) != blk(b)]
assert len(ROOTS) == 84
N = {}                                                   # root neighbours of the central roots
for r, lst in data.items():
    nb = {frozenset(p) for p in lst}
    assert len(nb) == 12 and all(len(p) == 2 and blk(min(p)) != blk(max(p)) for p in nb), r
    N[r] = nb
    N[sig(r)] = {sig(p) for p in nb}
central = list(N)
assert len(central) == 10
# partial graph on x, interiors, roots
V = ['x'] + list(range(14)) + ROOTS
adj = {v: set() for v in V}
def add(a, b):
    adj[a].add(b); adj[b].add(a)
for i in range(14):
    add('x', i)
    add(i, i ^ 1)                                         # block edge
for r in ROOTS:
    for i in r:
        add(r, i)
known_pairs = set()
for r in central:
    for s in N[r]:
        add(r, s)
problems = []
# 1. symmetry: if s is central, r must be in N[s]; degrees
for r in central:
    for s in N[r]:
        if s in N and r not in N[s]:
            problems.append(('asym', r, s))
    if len(adj[r]) != 14:
        problems.append(('deg', sorted(r), len(adj[r])))
# 2. lambda/mu for central roots vs every vertex whose neighbourhood is known (x, interiors, central roots)
full = ['x'] + list(range(14)) + central                 # vertices with fully known neighbourhoods
for a, b in itertools.combinations(full, 2):
    if not (a in N or b in N):
        continue                                         # interiors/x among themselves: fixed structure
    # the common neighbourhood is exact only if at least one side is fully known
    if isinstance(a, int) or a == 'x':
        a, b = b, a
    common = len(adj[a] & adj[b])
    target = 1 if b in adj[a] else 2
    if common != target:
        problems.append(('lm', a if isinstance(a, str) else sorted(a), b if not isinstance(b, frozenset) else sorted(b), common, target))
# 3. budget on all pairs
over = 0
for a, b in itertools.combinations(V, 2):
    common = len(adj[a] & adj[b])
    if b in adj[a]:
        cap = 1
    else:
        cap = 2
    if common > cap:
        over += 1
# 4. orbit rows
orbits = []
seen = set()
for r in ROOTS:
    if r not in seen:
        orbits.append((r, sig(r))); seen |= {r, sig(r)}
OI = {r: t for t, o in enumerate(orbits) for r in o}
rows = {}
for rep in data:
    q = [0] * 42
    for s in N[rep]:
        q[OI[s]] += 1
    rows[OI[rep]] = q
row_ok = all(sum(q) == 12 and q.count(2) == 1 and max(q) <= 2 for q in rows.values())
u, w, wp = R(0, 2), R(4, 6), R(4, 7)
prof = (rows[OI[u]][OI[w]], rows[OI[u]][OI[wp]], wp in N[u], sig(wp) in N[u])
# 5. orbit identity on the 5 x 5 block
blocks_of = lambda o: frozenset(blk(i) for i in orbits[o][0])
ident_bad = 0
for U, W in itertools.combinations_with_replacement(list(rows), 2):
    lhs = sum(rows[U][Z] * rows[W][Z] for Z in range(42)) + rows[U][W]
    sh = len(blocks_of(U) & blocks_of(W))
    rhs = 12 * (U == W) + 4 - sh
    ident_bad += lhs != rhs
print(f'central roots: {len(central)}; problems in degrees / symmetry / lambda-mu: {len(problems)} {problems[:5]}')
print(f'pairs of the partial graph over their common-neighbour budget: {over}')
print(f'orbit rows (sum 12, a single 2): {row_ok}; profile (q_U(W), q_U(W\'), u~w\', u~sigma w\') = {prof}')
print(f'orbit identity Q^2 + Q = 12I + 4J - Sh on the 5 x 5 central block: {ident_bad} violations')
