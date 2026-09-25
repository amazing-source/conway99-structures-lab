"""Anti-invariant side (C2 §0): S (49x49) has spectrum exactly {3, -4}, so every principal submatrix S_W has spectrum in
[-4, 3] (Cauchy interlacing), and a vector reaching 3 or -4 on S_W is a global S-eigenvector.
In a pure twisted configuration, at a block k the structure of the 12 orbits through k is known: 3 digons (S = 0 inside a
digon) and one free 6-cycle alternating m (+1) / pi (-1) (reps chosen containing (k,0)).  With the 7 interior orbits
(S = -1 on the diagonal, S_{U,i} = C_Ui), the whole block W = 7 interiors + 12 orbits through k is known up to labels.
Compute its spectrum for every labelling: which orbit of each pair kr is free (2^6), the cyclic order of the 6 free orbits
(5!/2), which pairs are digon-matched (15 matchings), and where the m/pi alternation starts (2)."""
import itertools
import numpy as np

blocks = [r for r in range(7) if r != 0]            # k = 0
def build(pfree, order, matching, start):
    # index: interiors 0..6, then orbits (r,p) for r in blocks, p in {0,1}
    idx = {('I', i): i for i in range(7)}
    for r in blocks:
        for p in (0, 1):
            idx[(r, p)] = len(idx)
    n = len(idx)
    S = np.zeros((n, n))
    for i in range(7):
        S[i, i] = -1
    for r in blocks:
        for p in (0, 1):
            u = idx[(r, p)]
            S[u, 0] = S[0, u] = 1                     # rep contains (k,0)
            S[u, r] = S[r, u] = 1 if p == 0 else -1   # C_{U,r}
    # free 6-cycle with alternating signs
    cyc = [(r, pfree[r]) for r in order]
    for t in range(6):
        a, b = idx[cyc[t]], idx[cyc[(t + 1) % 6]]
        s = 1 if (t + start) % 2 == 0 else -1
        S[a, b] = S[b, a] = s
    # digons: S = 0 between D-paired orbits (nothing to add)
    return S

def matchings(xs):
    if not xs:
        yield []
        return
    a = xs[0]
    for b in xs[1:]:
        rest = [x for x in xs if x not in (a, b)]
        for m in matchings(rest):
            yield [(a, b)] + m

spec_min, spec_max, touch3, touchm4 = 9, -9, 0, 0
seen = set()
tot = 0
for pf in itertools.product((0, 1), repeat=6):
    pfree = dict(zip(blocks, pf))
    for perm in itertools.permutations(blocks[1:]):
        order = [blocks[0]] + list(perm)
        if order[1] > order[-1]:
            continue
        for start in (0, 1):
            S = build(pfree, order, None, start)
            ev = np.linalg.eigvalsh(S)
            tot += 1
            spec_min, spec_max = min(spec_min, ev[0]), max(spec_max, ev[-1])
            touch3 += np.isclose(ev[-1], 3, atol=1e-9)
            touchm4 += np.isclose(ev[0], -4, atol=1e-9)
            seen.add(tuple(np.round(ev, 6)))
print(f'labellings {tot}; spectrum range [{spec_min:.6f}, {spec_max:.6f}]; reaches 3: {touch3}; reaches -4: {touchm4}; '
      f'distinct spectra {len(seen)}')
for sp in sorted(seen)[:6]:
    print('  ', sp)
