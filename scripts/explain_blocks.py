"""Human-readable explanation of the 'impossible outright' verdicts of doubles_budget3 (types with >= 3 triangles).
For one labelled decomposition per isomorphism class (canonical form under S7), and each block j, list for every free
orbit (indexed by its other block i) the admissible pairs {l, r} of R_j-neighbours (local system feasible), then count
the cyclic orders compatible with all of them.  A block with 0 compatible orders kills the decomposition."""
import itertools
import sys

sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import decompositions, CYC6, V


def canon_dec(dec):
    best = None
    edges_by_cycle = [frozenset(frozenset((c[t], c[(t + 1) % len(c)])) for t in range(len(c))) for c in dec]
    for perm in itertools.permutations(V):
        form = tuple(sorted(tuple(sorted(tuple(sorted((perm[a], perm[b]))) for a, b in (tuple(e) for e in cyc)))
                            for cyc in edges_by_cycle))
        if best is None or form < best:
            best = form
    return best


def block_report(dec, j):
    cyc_of_edge = {}
    for c in dec:
        m = len(c)
        for t in range(m):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % m]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]

    def other_nb(i, j_):
        c = cyc_of_edge[frozenset((i, j_))]
        m = len(c)
        t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]
        return b if a == j_ else a

    others = [r for r in V if r != j]
    adm = {}
    for i in others:
        k = other_nb(i, j)
        Rx = [r for r in V if r not in (i, j, k)]
        ok = []
        for l, r in itertools.combinations([o for o in others if o != i], 2):
            if d3.local_min(d3.canon(tris, i, j, k, Rx, frozenset((l, r)))) is not None:
                ok.append((l, r))
        adm[i] = (k, ok)
    orders = 0
    for order in CYC6:
        lab = [others[p] for p in order]
        if all(tuple(sorted((lab[t - 1], lab[(t + 1) % 6]))) in adm[lab[t]][1] for t in range(6)):
            orders += 1
    return adm, orders


if __name__ == '__main__':
    seen = {}
    for dec in decompositions((1 << 21) - 1, []):
        if sum(1 for c in dec if len(c) == 3) < 3:
            continue
        typ = tuple(sorted(len(c) for c in dec))
        if typ in seen:          # one representative per TYPE (illustration only; the proof is the exhaustive run)
            continue
        seen[typ] = dec
    print(f'types with >= 3 triangles: {len(seen)} (one representative each)')
    for key, dec in seen.items():
        typ = tuple(sorted(len(c) for c in dec))
        print(f'\n=== type {typ}: cycles {dec}')
        for j in V:
            adm, orders = block_report(dec, j)
            ntri = sum(1 for c in dec if len(c) == 3 and j in c)
            line = f'  block {j} (on {ntri} triangles): compatible cyclic orders = {orders}'
            if orders == 0:
                line += '   <-- KILLS'
                for i, (k, ok) in adm.items():
                    line += f'\n      orbit {i} (k={k}): admissible neighbour pairs {ok}'
            print(line)
