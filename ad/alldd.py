"""All-dd configuration (every pair doubly disjoint): pair-level model with a FIXED dd D-graph, to see what survives.
D-graph on pairs = 2-regular multigraph on the 21 pairs with edges between disjoint pairs; here: try a few natural ones."""
import itertools, time, random, sys
from general_model import feasible, PAIRS, pr
def skeleton_from_pair_cycles(cycles):
    """cycles: list of tuples of pairs (length >= 2), consecutive pairs disjoint (cyclically).  Orbit naming as cycle_matching."""
    idx = {(P, s): 2 * i + s for i, P in enumerate(PAIRS) for s in (0, 1)}
    orbs = [{'pair': P, 'pivot': None, 'partner': None} for P in PAIRS for s in (0, 1)]
    for c in cycles:
        L = len(c)
        if L == 2:
            a, b = c
            for s in (0, 1):
                orbs[idx[(a, s)]]['partner'] = idx[(b, s)]; orbs[idx[(b, s)]]['partner'] = idx[(a, s)]
        else:
            for i in range(L):
                a, b = c[i], c[(i + 1) % L]
                orbs[idx[(a, 1)]]['partner'] = idx[(b, 0)]; orbs[idx[(b, 0)]]['partner'] = idx[(a, 1)]
    return orbs
def random_cycles(rng):
    # random 2-factor of the Kneser graph K(7,2) (edges = disjoint pairs), as a permutation without fixed points, P -> disjoint
    while True:
        perm = {}; avail = list(PAIRS); rng.shuffle(avail); ok = True
        for P in PAIRS:
            c = [Q for Q in avail if not set(P) & set(Q)]
            if not c: ok = False; break
            Q = rng.choice(c); perm[P] = Q; avail.remove(Q)
        if ok: break
    # turn the permutation into cycles; a 2-cycle P<->Q is a digon, longer cycles directed
    seen = set(); cycles = []
    for P in PAIRS:
        if P in seen: continue
        c = [P]; seen.add(P); Q = perm[P]
        while Q != P: c.append(Q); seen.add(Q); Q = perm[Q]
        if len(c) == 1: raise ValueError
        cycles.append(tuple(c))
    return cycles
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 5):
    cyc = random_cycles(rng)
    orbs = skeleton_from_pair_cycles(cyc)
    t0 = time.time(); r = feasible(orbs, time_limit=120)
    st = 'dead' if r is None else ('unknown' if isinstance(r, str) else 'ALIVE')
    print('cycle lengths', sorted(len(c) for c in cyc), st, '%.0fs' % (time.time() - t0), flush=True)
    if st == 'ALIVE':
        x, phi, ns = r
        print('   n_s =', {P: int(round(x[j])) for P, j in ns.items() if round(x[j])})
