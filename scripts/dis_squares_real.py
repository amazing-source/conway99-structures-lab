"""Exploration on a real graph (no solver): D-square types in BvLS = srg(243,22,1,2) at x = 0 with sigma(v) = -v.
For every far orbit U: its D-partner W (the unique orbit with q_U(W) = 2); type sib / one-block / dis from the block
pairs; for dis: q_U(W') with W' the sibling of W ((2,0) or (2,1)); also n_s of every pair."""
import sys
from collections import Counter
import numpy as np
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
import graphs

A, V, S = graphs.bvls(return_coords=True)
A = np.array(A, dtype=np.int64)
idx = {v: t for t, v in enumerate(V)}
neg = lambda v: tuple((-np.array(v)) % 3)
x = idx[(0,) * 5]
Nx = [u for u in range(len(V)) if A[x, u]]
far = [u for u in range(len(V)) if u != x and not A[x, u]]
block = {}
for u in Nx:
    (m,) = [w for w in Nx if w != u and A[u, w]]
    block[u] = frozenset((u, m))
bp = {w: frozenset(block[u] for u in Nx if A[w, u]) for w in far}         # block pair of a far vertex
orb, seen = [], set()
for w in far:
    if w not in seen:
        sw = idx[neg(V[w])]
        orb.append((w, sw)); seen |= {w, sw}
OI = {w: t for t, o in enumerate(orb) for w in o}
q = lambda a, b: int(A[orb[a][0], orb[b][0]] + A[orb[a][0], orb[b][1]])
sib = {}
for a in range(len(orb)):
    for b in range(len(orb)):
        if a != b and bp[orb[a][0]] == bp[orb[b][0]]:
            sib[a] = b
types, dis_profile, ns = Counter(), Counter(), Counter()
for a in range(len(orb)):
    (W,) = [b for b in range(len(orb)) if q(a, b) == 2]
    share = len(bp[orb[a][0]] & bp[orb[W][0]])
    t = 'sib' if share == 2 else ('one-block' if share == 1 else 'dis')
    types[t] += 1
    if t == 'dis':
        dis_profile[(2, q(a, sib[W]))] += 1
for a in range(len(orb)):
    if a < sib[a]:
        ns[q(a, sib[a])] += 1
print(f'BvLS at x=0, sigma = -v: {len(orb)} far orbits; D-square types (per orbit): {dict(types)}')
print(f'dis orbits: profile on the D-partner pair (q_U(W), q_U(W\')): {dict(dis_profile)}')
print(f'n_s over the {sum(ns.values())} pairs: {dict(ns)}')
