"""Own solver of the local system C0-C8 of a pure twisted orbit, in the ORIGINAL block names (independent of the
canonical keys of doubles_budget3).  Orbit i free at block j, paired via i to {i,k}, k = other neighbour of i on the
cycle through ij, R_x = V - {i,j,k}; unknowns phi(P) (P pair in R_x) and y_r (r in R_x), all in {0,1,2}."""
import itertools

from zeros_budget import V


def P(a, b):
    return (min(a, b), max(a, b))


def frame(dec, j, i):
    for c in dec:
        m = len(c)
        for t in range(m):
            if {c[t], c[(t + 1) % m]} == {i, j}:
                cyc = c
    m = len(cyc)
    t = cyc.index(i)
    a, b = cyc[(t - 1) % m], cyc[(t + 1) % m]
    k = b if a == j else a
    return cyc, k, [r for r in V if r not in (i, j, k)]


def solve_core(tris, own_tri, i, j, k, R, pair):
    """All solutions (phi, y) of C0-C8 for the orbit i free at j, paired via i to {i,k}; tris excludes the own cycle."""
    n = {r: int(r in pair) for r in V}
    if own_tri and n[k]:
        return []                                                    # C0
    pairs = [P(a, b) for a, b in itertools.combinations(R, 2)]
    phi_fixed, y_fixed, c4, c5 = [], [], [], []
    for t in tris:
        if i in t:
            a, b = sorted(t - {i})
            phi_fixed.append((P(a, b), 1))                           # C1
        elif j in t and k in t:
            (b,) = t - {j, k}
            y_fixed.append((b, 1 - n[k] - n[b]))                     # C3
        elif j in t:
            a, b = sorted(t - {j})
            phi_fixed.append((P(a, b), 1 - n[a] - n[b]))             # C2
        elif k in t:
            c4.append(tuple(sorted(t - {k})))                        # C4
        else:
            c5.append(tuple(sorted(t)))                              # C5
    sols = []
    for yv in itertools.product(range(3), repeat=4):
        if sum(yv) != 2 - n[k]:                                      # C6
            continue
        y = dict(zip(R, yv))
        if any(y[b] != v for b, v in y_fixed):
            continue
        for pv in itertools.product(range(3), repeat=6):
            phi = dict(zip(pairs, pv))
            if any(phi[p] != v for p, v in phi_fixed):
                continue
            if any(y[a] + y[b] + phi[P(a, b)] != 3 for a, b in c4):
                continue
            if any(phi[P(a, b)] + phi[P(a, c)] + phi[P(b, c)] != 3 for a, b, c in c5):
                continue
            if any(n[t] + y[t] + sum(phi[P(t, s)] for s in R if s != t) != 4 for t in R):   # C8
                continue
            sols.append((phi, y))
    return sols


def doubles(sol):
    phi, y = sol
    return sum(v == 2 for v in phi.values()) + sum(v == 2 for v in y.values())
