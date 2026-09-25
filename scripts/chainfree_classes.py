"""Isomorphism classes of the chain-free pure twisted decompositions with >= 3 triangles (other than Fano), found by
structural normalisation (small symmetry groups instead of S7), and their critical blocks for the hand proofs.

Families:
  STAR  : exactly 3 triangles, all through one block v; the other 12 edges (octahedron on the 6 other blocks,
          antipodes = the triangle partners) are split into (6,6) or (4,4,4).
  QUAD  : 4 triangles = star of 3 at v + one transversal D (meets each star triangle outside v); a 4-cycle + a 5-cycle.
  PASCH : 3 triangles pairwise meeting in 3 distinct blocks, no common block; the 7th block w lies on no triangle;
          three 4-cycles.
For each class: size (labelled count), one representative, and for every block the number of compatible cyclic
orders (0 = critical, kills outright), with the admissible neighbour pairs, or, when no block is critical, the minimal
doubles per block (budget argument)."""
import itertools
import sys
import time

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions, V
from chain_free import has_chain
from explain_blocks import block_report
import doubles_budget3 as d3


def cyc_edges(c):
    return frozenset(frozenset((c[t], c[(t + 1) % len(c)])) for t in range(len(c)))


def form(dec, perm):
    return tuple(sorted(tuple(sorted(tuple(sorted((perm[a], perm[b]))) for a, b in (tuple(e) for e in cyc_edges(c))))
                        for c in dec))


def family_and_canon(dec):
    tris = [frozenset(c) for c in dec if len(c) == 3]
    deg = {v: sum(1 for t in tris if v in t) for v in V}
    if len(tris) == 3 and any(deg[v] == 3 for v in V):
        v = next(u for u in V if deg[u] == 3)
        pairs = [tuple(sorted(t - {v})) for t in tris]
        perms = []
        for order in itertools.permutations(range(3)):
            for flips in itertools.product((0, 1), repeat=3):
                img = {v: 0}
                for slot, idx in enumerate(order):
                    a, b = pairs[idx]
                    if flips[slot]:
                        a, b = b, a
                    img[a], img[b] = 1 + 2 * slot, 2 + 2 * slot
                perms.append(img)
        return 'STAR', min(form(dec, p) for p in perms)
    if len(tris) == 4 and any(deg[v] == 3 for v in V):
        v = next(u for u in V if deg[u] == 3)
        D = next(t for t in tris if v not in t)
        star = [t for t in tris if v in t]
        perms = []
        for order in itertools.permutations(range(3)):
            img = {v: 0}
            for slot, idx in enumerate(order):
                t = star[idx]
                d = next(iter(t & D))
                o = next(iter(t - {v, d}))
                img[d], img[o] = 1 + 2 * slot, 2 + 2 * slot
            perms.append(img)
        return 'QUAD', min(form(dec, p) for p in perms)
    if len(tris) == 3 and all(deg[v] <= 2 for v in V):
        w = next(u for u in V if deg[u] == 0)
        perms = []
        for order in itertools.permutations(range(3)):
            T = [tris[i] for i in order]
            img = {w: 0}
            # intersections: p = T0∩T1, q = T1∩T2, r = T2∩T0 ; privates a in T0, b in T1, c in T2
            p = next(iter(T[0] & T[1])); q = next(iter(T[1] & T[2])); r = next(iter(T[2] & T[0]))
            a = next(iter(T[0] - {p, r})); b = next(iter(T[1] - {p, q})); c = next(iter(T[2] - {q, r}))
            img.update({p: 1, q: 2, r: 3, a: 4, b: 5, c: 6})
            perms.append(img)
        return 'PASCH', min(form(dec, p) for p in perms)
    return None, None


if __name__ == '__main__':
    t0 = time.time()
    classes = {}
    for dec in decompositions((1 << 21) - 1, []):
        tris = [frozenset(c) for c in dec if len(c) == 3]
        if len(tris) < 3 or len(tris) == 7 or has_chain(tris):
            continue
        fam, key = family_and_canon(dec)
        if fam is None:
            print('UNCLASSIFIED', dec)
            continue
        classes.setdefault((fam, key), []).append(dec)
    print(f'classes: {len(classes)}  [{time.time() - t0:.0f} s]')
    for (fam, key), decs in sorted(classes.items(), key=lambda kv: (kv[0][0], len(kv[1]))):
        dec = decs[0]
        typ = tuple(sorted(len(c) for c in dec))
        total = d3.analyse(dec)[0]
        print(f'\n=== {fam} {typ}: {len(decs)} labelled; representative {dec}; exhaustive verdict: '
              f'{"impossible outright" if total is None else f"lower bound {total}"}')
        for j in V:
            adm, orders = block_report(dec, j)
            ntri = sum(1 for c in dec if len(c) == 3 and j in c)
            line = f'  block {j} (on {ntri} triangles): compatible orders = {orders}'
            if orders == 0:
                line += '  <-- critical'
                for i, (k, ok) in adm.items():
                    line += f'\n      orbit {i} (k={k}): admissible {ok}'
            print(line)
    print(f'\n[{time.time() - t0:.0f} s]')
