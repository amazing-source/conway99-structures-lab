"""Independent verification (lab 'structures', own code) of C1's PAIR x BLOCK certificates (X19n/X19q, exported in
the x19k format: <base>.json with tau, sigma, specials, rows, cols, farkas_scaled, D; <base>.npz with the CSR matrix
A, b, y_int, D).  Exact integers only; nothing imported from C1's scripts.
Model (C1's X19n): x[z,t] = type t of the far vertex (pair) z; s[z,t,R,S] = the exact subset S of the 4 pairs of a
block R (disjoint from the classes of z) that z is adjacent to, |S| = t[R].
  1. specials recomputed from tau, sigma and compared; types recomputed by enumerating the label graphs G_z;
     completeness: for every x-column and every block R of its type, the s-columns are exactly ALL subsets of size t[R]
     of the 4 pairs of R (a real configuration's subset is always a column);
  2. every row rebuilt FROM ITS NAME with own semantics:
       conv|z            sum_t x = 1
       cons|i|R|z        sum_S s[z,i,R,S] = x[z,i]
       labdeg|i|l|z      sum_R sum_S s[z,i,R,S] #{y in S : l in y} = (2 - h_z(l)) x[z,i]      (l outside z's classes)
       bal|R|R', lab|u|R, P2|R|R'   block level, as in X19 (constants to the right)
       pairblk|y|R'      sum_{z != y} [z ~ y] cnt(z,R') + cnt(y,R') = sum_{w in R'} (12 [w = y] + (2 - |y cap w|) [w != y])
                         ((B^2)_{yw} = 12 [y = w] + 2 - |y cap w| - [y ~ w]; the [y ~ w] terms moved to the left)
       pairpair|y|y'     sum_{z not in {y,y'}} [z ~ y][z ~ y'] = 2 - |y cap y'| - [y' in spec(y)]  (same block)
     [z ~ y] is the constant [y in spec(z)] when the blocks of y and z meet, else sum of s[z,i,R_y,S] over S containing y;
     a product with cnt(z,R') is carried by the s-columns of R_y (times t[R'] or the special count), or by the
     x-columns when [z ~ y] is constant.  Rebuilt coefficients and right-hand sides are compared with the file;
  3. Farkas on the REBUILT A, b: b.y - sum_j max(0, (A^T y)_j) > 0, all variables in [0,1]."""
import itertools
import json
import sys
from collections import Counter

import numpy as np

base = sys.argv[1]
D = json.load(open(base + '.json'))
Z = np.load(base + '.npz')

lab = lambda s: tuple(int(t) for t in s.split('.'))
mate = lambda u: (u[0], 1 - u[1])
PAIRS = [frozenset({(c, s), (d, t)}) for c, d in itertools.combinations(range(7), 2) for s in (0, 1) for t in (0, 1)]
BLOCKS = [frozenset(b) for b in itertools.combinations(range(7), 2)]
blk = lambda z: frozenset(l[0] for l in z)
pn = lambda z: '-'.join(f'{l[0]}.{l[1]}' for l in sorted(z))
PN = {pn(z): z for z in PAIRS}
Bn = lambda s: frozenset(int(ch) for ch in s)
inblock = {R: [w for w in PAIRS if blk(w) == R] for R in BLOCKS}

# ---- 1. specials, types, completeness
tau = {}
for u, lst in D['tau'].items():
    m = {}
    for v, w in lst:
        m[lab(v)] = lab(w); m[lab(w)] = lab(v)
    tau[lab(u)] = m
sigma = {int(c): {lab(v): lab(w) for v, w in lst} for c, lst in D['sigma'].items()}
sinv = {c: {w: v for v, w in m.items()} for c, m in sigma.items()}
spec = {}
for z in PAIRS:
    e, f = sorted(z)
    S = [frozenset({e, tau[e][f]}), frozenset({f, tau[f][e]})]
    for u, v in ((e, f), (f, e)):
        w = sigma[u[0]][v] if u[1] == 0 else sinv[u[0]][v]
        S.append(frozenset({mate(u), w}))
    spec[z] = S
spec_same = all(Counter(pn(w) for w in spec[z]) == Counter(D['specials'][pn(z)]) for z in PAIRS)
spec_sym = all(z in spec[w] for z in PAIRS for w in spec[z])
spec_distinct = all(len(set(spec[z])) == 4 and z not in spec[z] for z in PAIRS)
need = {}
for z in PAIRS:
    P = blk(z)
    h = Counter(l for w in spec[z] for l in w if l[0] not in P)
    need[z] = {(c, s): 2 - h[(c, s)] for c in range(7) if c not in P for s in (0, 1)}


def types_of(z):
    labs = sorted(need[z])
    nd = dict(need[z])
    if min(nd.values()) < 0:
        return set()
    cand = [frozenset({a, b}) for a, b in itertools.combinations(labs, 2) if a[0] != b[0]]
    found = set()

    def rec(i, chosen):
        if all(v == 0 for v in nd.values()):
            found.add(frozenset(Counter(blk(e) for e in chosen).items()))
            return
        if i == len(cand):
            return
        a, b = tuple(cand[i])
        if nd[a] > 0 and nd[b] > 0:
            nd[a] -= 1; nd[b] -= 1
            rec(i + 1, chosen + [cand[i]])
            nd[a] += 1; nd[b] += 1
        rec(i + 1, chosen)
    rec(0, [])
    return found


cols = D['cols']
CI, xtype = {}, {}
for j, name in enumerate(cols):
    p = name.split('|')
    z = PN[p[1]]; i = int(p[2])
    if p[0] == 'x':
        t = {}
        for item in p[3].split(','):
            R, m = item.split(':')
            t[Bn(R)] = int(m)
        xtype[(z, i)] = t
        CI[('x', z, i)] = j
    elif p[0] == 's':
        R = Bn(p[3]); S = frozenset(PN[q] for q in p[4].split('+'))
        assert all(blk(y) == R for y in S), name
        CI[('s', z, i, R, S)] = j
    else:
        raise ValueError(name)
type_mismatch = dup = 0
for z in PAIRS:
    theirs = [frozenset(t.items()) for (zz, i), t in xtype.items() if zz == z]
    dup += len(theirs) - len(set(theirs))
    type_mismatch += set(theirs) != types_of(z)
xs = {z: [(i, t) for (zz, i), t in xtype.items() if zz == z] for z in PAIRS}
subsets = {(z, i, R): [] for (z, i), t in xtype.items() for R in t}
for k in CI:
    if k[0] == 's':
        subsets.setdefault((k[1], k[2], k[3]), []).append(k[4])
split_bad = 0
for (z, i), t in xtype.items():
    for R in BLOCKS:
        have = set(subsets.get((z, i, R), []))
        want = {frozenset(c) for c in itertools.combinations(inblock[R], t[R])} if R in t else set()
        split_bad += have != want
n_s = sum(1 for k in CI if k[0] == 's')
print(f'specials equal: {spec_same}, symmetric: {spec_sym}, 4 distinct: {spec_distinct}; {len(xtype)} x-columns, type '
      f'sets differing from own enumeration: {type_mismatch}, duplicates: {dup}; {n_s} s-columns, (x-column, block) '
      f'whose subsets are not ALL subsets of size t[R]: {split_bad}', flush=True)


# ---- 2. rows from names
def cb(z, R):
    if R & blk(z):
        return sum(1 for w in spec[z] if blk(w) == R), {}
    return 0, {CI[('x', z, i)]: t[R] for i, t in xs[z] if t.get(R, 0)}


def cnt_val(z, i, t, R):
    """cnt(z,R) for the type t of z (a number)."""
    return sum(1 for w in spec[z] if blk(w) == R) if R & blk(z) else t.get(R, 0)


def add(F, G, s=1):
    for k, v in G.items():
        F[k] = F.get(k, 0) + s * v


def row(name):
    p = name.split('|')
    kind = p[0]
    F, c0 = {}, 0
    if kind == 'conv':
        z = PN[p[1]]
        return {CI[('x', z, i)]: 1 for i, t in xs[z]}, 1
    if kind == 'cons':
        i, R, z = int(p[1]), Bn(p[2]), PN[p[3]]
        F = {CI[('s', z, i, R, S)]: 1 for S in subsets[(z, i, R)]}
        F[CI[('x', z, i)]] = -1
        return F, 0
    if kind == 'labdeg':
        i, l, z = int(p[1]), lab(p[2]), PN[p[3]]
        assert l[0] not in blk(z)
        t = xtype[(z, i)]
        for R in t:
            if l[0] in R:
                for S in subsets[(z, i, R)]:
                    k = sum(1 for y in S if l in y)
                    if k:
                        F[CI[('s', z, i, R, S)]] = F.get(CI[('s', z, i, R, S)], 0) + k
        F[CI[('x', z, i)]] = F.get(CI[('x', z, i)], 0) - need[z][l]
        return F, 0
    if kind == 'bal':
        R, Rp = Bn(p[1]), Bn(p[2])
        for w in inblock[R]:
            c, G = cb(w, Rp); c0 += c; add(F, G)
        for w in inblock[Rp]:
            c, G = cb(w, R); c0 -= c; add(F, G, -1)
        return F, -c0
    if kind == 'lab':
        u, R = lab(p[1]), Bn(p[2])
        for z in PAIRS:
            if u in z:
                c, G = cb(z, R); c0 += c; add(F, G)
        return F, 8 - c0
    if kind == 'P2':
        R, Rp = Bn(p[1]), Bn(p[2])
        for z in PAIRS:
            if (R & blk(z)) and (Rp & blk(z)):
                c0 += cb(z, R)[0] * cb(z, Rp)[0]; continue
            for i, t in xs[z]:
                v = cnt_val(z, i, t, R) * cnt_val(z, i, t, Rp)
                if v:
                    F[CI[('x', z, i)]] = F.get(CI[('x', z, i)], 0) + v
        for w in inblock[R]:
            c, G = cb(w, Rp); c0 += c; add(F, G)
        r = 12 * len(set(inblock[R]) & set(inblock[Rp])) + 32 - sum(len(a & b) for a in inblock[R] for b in inblock[Rp])
        return F, r - c0
    if kind == 'pairblk':
        y, Rp = PN[p[1]], Bn(p[2])
        Ry = blk(y)
        for z in PAIRS:
            if z == y:
                continue
            if Ry & blk(z):                          # [z ~ y] constant
                a = int(y in spec[z])
                if not a:
                    continue
                if Rp & blk(z):
                    c0 += a * cb(z, Rp)[0]
                else:
                    for i, t in xs[z]:
                        if t.get(Rp, 0):
                            F[CI[('x', z, i)]] = F.get(CI[('x', z, i)], 0) + a * t[Rp]
            else:                                    # [z ~ y] on the s-columns of block R_y
                for i, t in xs[z]:
                    cR = cnt_val(z, i, t, Rp)
                    if not cR or Ry not in t:
                        continue
                    for S in subsets[(z, i, Ry)]:
                        if y in S:
                            F[CI[('s', z, i, Ry, S)]] = F.get(CI[('s', z, i, Ry, S)], 0) + cR
        c, G = cb(y, Rp); c0 += c; add(F, G)       # the -[y ~ w] terms, moved to the left
        r = sum(12 if w == y else 2 - len(y & w) for w in inblock[Rp])
        return F, r - c0
    if kind == 'pairpair':
        y, y2 = PN[p[1]], PN[p[2]]
        R = blk(y)
        assert blk(y2) == R and y != y2
        for z in PAIRS:
            if z in (y, y2):
                continue
            if R & blk(z):
                c0 += int(y in spec[z] and y2 in spec[z]); continue
            for i, t in xs[z]:
                if R not in t:
                    continue
                for S in subsets[(z, i, R)]:
                    if y in S and y2 in S:
                        F[CI[('s', z, i, R, S)]] = F.get(CI[('s', z, i, R, S)], 0) + 1
        return F, 2 - len(y & y2) - int(y2 in spec[y]) - c0
    raise ValueError('unknown row ' + name)


data, indices, indptr = Z['data'], Z['indices'], Z['indptr']
b_file = Z['b']
yv = [int(v) for v in Z['y_int']]
coef_mismatch = rhs_mismatch = 0
bad_kinds = Counter()
w = [0] * len(cols)
by = 0
for r_, name in enumerate(D['rows']):
    F, rhs = row(name)
    F = {k: v for k, v in F.items() if v}
    theirs = {int(indices[t]): int(data[t]) for t in range(indptr[r_], indptr[r_ + 1]) if int(data[t])}
    if F != theirs:
        coef_mismatch += 1; bad_kinds[name.split('|')[0]] += 1
    if rhs != int(b_file[r_]):
        rhs_mismatch += 1; bad_kinds[name.split('|')[0] + '(rhs)'] += 1
    if yv[r_]:
        by += yv[r_] * rhs
        for k, v in F.items():
            w[k] += yv[r_] * v
pos = sum(v for v in w if v > 0)
kinds = Counter(n.split('|')[0] for n in D['rows'])
print(f'rows {len(D["rows"])} {dict(kinds)}: rebuilt coefficients differing: {coef_mismatch}; right-hand sides '
      f'differing: {rhs_mismatch} {dict(bad_kinds) if bad_kinds else ""}', flush=True)
print(f'Farkas on the REBUILT system: b.y - sum_j max(0, (A^T y)_j) = {by - pos}  (D = {int(Z["D"][0])}; file claims '
      f'{D["farkas_scaled"]})', flush=True)
ok = (by - pos > 0 and spec_same and spec_sym and spec_distinct and type_mismatch == 0 and dup == 0 and split_bad == 0
      and coef_mismatch == 0 and rhs_mismatch == 0)
print('CERTIFICATE VALID' if ok else 'NOT VALID (see above)')
