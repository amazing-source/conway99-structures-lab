"""Independent verification (lab 'structures', own code) of C1's CELL-level certificate data/x19k_witness_cells.{json,npz}
(X19k, notes/117 draft): the Fano-tight witness (tau, sigma) is impossible once the block-level system is refined by
the antipodal cells.  Exact integers only; nothing is imported from C1's scripts.
  1. specials of every pair recomputed from tau, sigma (own reading) and compared with the file; symmetric;
  2. types of every pair recomputed by enumerating the actual label graphs G_z (8 edges on the 10 labels outside the
     classes of z, no mate edge, simple, label degrees 2 - h_z) and compared with the x-columns; for every x-column and
     every block R of its type, the s-columns must be exactly ALL splits (k_A, k_B), k_A + k_B = t[R], 0 <= k <= 2
     (completeness: the true split of a real configuration is always a column);
  3. every row rebuilt FROM ITS NAME with own semantics, from
         sum_z c_z(C) c_z(C') + e(C,C') = 12 |C cap C'| + 2 |C| |C'| - L(C,C'),
         e(C,C') = sum_{w in C} c_w(C'),   L(C,C') = sum_{w in C, w' in C'} |w cap w'|,
     where c_z(C) is a constant (count of specials) when the block of C meets the classes of z, a type count on the
     x-columns for a whole block, a split count on the s-columns for a cell; constants go to the right-hand side.
     Rebuilt coefficients and right-hand sides are compared with the file's A and b entry by entry;
  4. Farkas on the REBUILT A, b:  b.y - sum_j max(0, (A^T y)_j) > 0  with every variable in [0,1]
     (x <= 1 by conv, s <= x by cons; a real configuration is a 0/1 point)."""
import itertools
import json
import sys
from collections import Counter

import numpy as np

base = sys.argv[1] if len(sys.argv) > 1 else r'C:\Users\RESEARCHER\Desktop\c1_experimental\data\x19k_witness_cells'
D = json.load(open(base + '.json'))
Z = np.load(base + '.npz')

lab = lambda s: tuple(int(t) for t in s.split('.'))
mate = lambda u: (u[0], 1 - u[1])
LABELS = [(c, s) for c in range(7) for s in (0, 1)]
PAIRS = [frozenset({(c, s), (d, t)}) for c, d in itertools.combinations(range(7), 2) for s in (0, 1) for t in (0, 1)]
BLOCKS = [frozenset(b) for b in itertools.combinations(range(7), 2)]
blk = lambda z: frozenset(l[0] for l in z)
cell = lambda w: (blk(w), 'A' if len({l[1] for l in w}) == 1 else 'B')
pn = lambda z: '-'.join(f'{l[0]}.{l[1]}' for l in sorted(z))
PN = {pn(z): z for z in PAIRS}
Bn = lambda s: frozenset(int(ch) for ch in s)
members = {}
for w in PAIRS:
    members.setdefault(cell(w), []).append(w)
inblock = {R: [w for w in PAIRS if blk(w) == R] for R in BLOCKS}

# ---- 1. specials
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


# ---- 2. types by enumerating label graphs
def types_of(z):
    P = blk(z)
    labs = [(c, s) for c in range(7) if c not in P for s in (0, 1)]
    h = {u: sum(1 for w in spec[z] if u in w) for u in labs}
    need = {u: 2 - h[u] for u in labs}
    if min(need.values()) < 0:
        return set()
    cand = [frozenset({a, b}) for a, b in itertools.combinations(labs, 2) if a[0] != b[0]]
    found = set()

    def rec(i, chosen):
        if all(v == 0 for v in need.values()):
            found.add(frozenset(Counter(blk(e) for e in chosen).items()))
            return
        if i == len(cand):
            return
        a, b = tuple(cand[i])
        if need[a] > 0 and need[b] > 0:
            need[a] -= 1; need[b] -= 1
            rec(i + 1, chosen + [cand[i]])
            need[a] += 1; need[b] += 1
        rec(i + 1, chosen)
    rec(0, [])
    return found


cols = D['cols']
CI = {}
xtype = {}                                      # (z, i) -> {R: m}
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
    else:
        R = Bn(p[3]); kA, kB = int(p[4][1]), int(p[4][3])
        CI[('s', z, i, R, (kA, kB))] = j
type_mismatch = 0
dup = 0
for z in PAIRS:
    mine = types_of(z)
    theirs = [frozenset(t.items()) for (zz, i), t in xtype.items() if zz == z]
    dup += len(theirs) - len(set(theirs))
    type_mismatch += set(theirs) != mine
split_bad = 0
for (z, i), t in xtype.items():
    for R in BLOCKS:
        have = {k for k in [(a, b) for a in range(3) for b in range(3)] if ('s', z, i, R, k) in CI}
        want = {(a, t[R] - a) for a in range(3) if 0 <= t[R] - a <= 2} if R in t else set()
        split_bad += have != want
n_s = sum(1 for k in CI if k[0] == 's')
print(f'specials equal to the file: {spec_same}; symmetric: {spec_sym}; types: {len(xtype)} x-columns, pairs whose type '
      f'set differs from own enumeration: {type_mismatch}, duplicates: {dup}; s-columns {n_s}, (x-column, block) with a '
      f'split set different from ALL splits: {split_bad}', flush=True)

# ---- 3. rows from names
xs = {z: [(i, t) for (zz, i), t in xtype.items() if zz == z] for z in PAIRS}


def cb(z, R):
    """c_z(block R): (constant, {col: coef})."""
    if R & blk(z):
        return sum(1 for w in spec[z] if blk(w) == R), {}
    return 0, {CI[('x', z, i)]: t[R] for i, t in xs[z] if t.get(R, 0)}


def cc(z, C):
    """c_z(cell C)."""
    R, h = C
    if R & blk(z):
        return sum(1 for w in spec[z] if cell(w) == C), {}
    F = {}
    for i, t in xs[z]:
        if R in t:
            for a in range(3):
                k = (a, t[R] - a)
                if ('s', z, i, R, k) in CI:
                    v = k[0] if h == 'A' else k[1]
                    if v:
                        F[CI[('s', z, i, R, k)]] = v
    return 0, F


def prod(z, X, Y):
    """c_z(X) * c_z(Y) for X, Y in {('b', R), ('c', (R, h))}; linear because the columns of z are exclusive.
    A cell factor with a variable block is carried by the s-columns, everything else by the x-columns."""
    def const_or_none(W):
        R = W[1] if W[0] == 'b' else W[1][0]
        if R & blk(z):
            return (cb(z, R) if W[0] == 'b' else cc(z, W[1]))[0]
        return None
    cX, cY = const_or_none(X), const_or_none(Y)
    if cX is not None and cY is not None:
        return cX * cY, {}
    F = {}
    if X[0] == 'c' and cX is None or Y[0] == 'c' and cY is None:
        # a variable cell factor: put the product on the s-columns of that cell's block
        if not (X[0] == 'c' and cX is None):
            X, Y, cX, cY = Y, X, cY, cX
        R, h = X[1]
        for i, t in xs[z]:
            if R not in t:
                continue
            for a in range(3):
                k = (a, t[R] - a)
                if ('s', z, i, R, k) not in CI:
                    continue
                vx = k[0] if h == 'A' else k[1]
                if Y[0] == 'b':
                    vy = cY if cY is not None else t.get(Y[1], 0)
                else:
                    if cY is not None:
                        vy = cY
                    else:
                        assert Y[1][0] == R, 'cell x cell across two variable blocks is not linear'
                        vy = k[0] if Y[1][1] == 'A' else k[1]
                if vx * vy:
                    F[CI[('s', z, i, R, k)]] = F.get(CI[('s', z, i, R, k)], 0) + vx * vy
        return 0, F
    for i, t in xs[z]:                     # only block factors (or constant cells) are variable: x-columns
        vx = cX if cX is not None else t.get(X[1], 0)
        vy = cY if cY is not None else t.get(Y[1], 0)
        if vx * vy:
            F[CI[('x', z, i)]] = vx * vy
    return 0, F


def add(F, G, s=1):
    for k, v in G.items():
        F[k] = F.get(k, 0) + s * v


def L(C1, C2):
    return sum(len(a & b) for a in C1 for b in C2)


def second_order(X, Y, setX, setY, eX):
    """sum_z c_z(X) c_z(Y) + e(X,Y) = 12|X cap Y| + 2|X||Y| - L(X,Y); returns (F, const on the left, rhs)."""
    F, c0 = {}, 0
    for z in PAIRS:
        c, G = prod(z, X, Y); c0 += c; add(F, G)
    for w in setX:
        c, G = eX(w); c0 += c; add(F, G)
    r = 12 * len(set(setX) & set(setY)) + 2 * len(setX) * len(setY) - L(setX, setY)
    return F, c0, r


def row(name):
    p = name.split('|')
    kind = p[0]
    F, c0 = {}, 0
    if kind == 'conv':
        z = PN[p[1]]
        return {CI[('x', z, i)]: 1 for i, t in xs[z]}, 1
    if kind == 'cons':
        i, R, z = int(p[1]), Bn(p[2]), PN[p[3]]
        F = {CI[('s', z, i, R, k)]: 1 for k in [(a, xtype[(z, i)][R] - a) for a in range(3)] if ('s', z, i, R, k) in CI}
        F[CI[('x', z, i)]] = -1
        return F, 0
    if kind == 'bal':
        R, Rp = Bn(p[1]), Bn(p[2])
        for w in inblock[R]:
            c, G = cb(w, Rp); c0 += c; add(F, G)
        for w in inblock[Rp]:
            c, G = cb(w, R); c0 -= c; add(F, G, -1)
        return F, 0 - c0
    if kind == 'lab':
        u, R = lab(p[1]), Bn(p[2])
        assert u[0] not in R
        for z in PAIRS:
            if u in z:
                c, G = cb(z, R); c0 += c; add(F, G)
        return F, 8 - c0
    if kind == 'P2':
        R, Rp = Bn(p[1]), Bn(p[2])
        F, c0, r = second_order(('b', R), ('b', Rp), inblock[R], inblock[Rp], lambda w: cb(w, Rp))
        return F, r - c0
    if kind == 'cellcross':
        R, Rp = Bn(p[1]), Bn(p[2])
        FA, cA, rA = second_order(('c', (R, 'A')), ('b', Rp), members[(R, 'A')], inblock[Rp], lambda w: cb(w, Rp))
        FB, cB, rB = second_order(('c', (R, 'B')), ('b', Rp), members[(R, 'B')], inblock[Rp], lambda w: cb(w, Rp))
        add(FA, FB, -1)
        return FA, (rA - rB) - (cA - cB)
    if kind == 'cellsame':
        R, h1, h2 = Bn(p[1]), p[2][0], p[2][1]
        C1, C2 = (R, h1), (R, h2)
        F, c0, r = second_order(('c', C1), ('c', C2), members[C1], members[C2], lambda w: cc(w, C2))
        return F, r - c0
    if kind == 'cellbal':
        C1 = (Bn(p[1][:-1]), p[1][-1]); C2 = (Bn(p[2][:-1]), p[2][-1])
        for w in members[C1]:
            c, G = cc(w, C2); c0 += c; add(F, G)
        for w in members[C2]:
            c, G = cc(w, C1); c0 -= c; add(F, G, -1)
        return F, 0 - c0
    if kind == 'celllab':
        u, C = lab(p[1]), (Bn(p[2][:-1]), p[2][-1])
        assert u[0] not in C[0]
        for z in PAIRS:
            if u in z:
                c, G = cc(z, C); c0 += c; add(F, G)
        return F, 2 * len(members[C]) - c0
    raise ValueError(name)


data, indices, indptr = Z['data'], Z['indices'], Z['indptr']
b_file = Z['b']
y = [int(v) for v in Z['y_int']]
coef_mismatch = rhs_mismatch = 0
w = [0] * len(cols)
by = 0
kinds_bad = Counter()
for r_, name in enumerate(D['rows']):
    F, rhs = row(name)
    F = {k: v for k, v in F.items() if v}
    theirs = {int(indices[t]): int(data[t]) for t in range(indptr[r_], indptr[r_ + 1]) if int(data[t])}
    if F != theirs:
        coef_mismatch += 1; kinds_bad[name.split('|')[0]] += 1
    if rhs != int(b_file[r_]):
        rhs_mismatch += 1; kinds_bad[name.split('|')[0] + '(rhs)'] += 1
    if y[r_]:
        by += y[r_] * rhs
        for k, v in F.items():
            w[k] += y[r_] * v
pos = sum(v for v in w if v > 0)
print(f'rows {len(D["rows"])}: rebuilt coefficients differing from the file: {coef_mismatch}; right-hand sides '
      f'differing: {rhs_mismatch} {dict(kinds_bad) if kinds_bad else ""}', flush=True)
print(f'Farkas on the REBUILT system: b.y - sum_j max(0, (A^T y)_j) = {by - pos}  (scale D = {int(Z["D"][0])}; '
      f'file claims {D["farkas_scaled"]})', flush=True)
ok = by - pos > 0 and spec_same and spec_sym and type_mismatch == 0 and dup == 0 and split_bad == 0 \
    and coef_mismatch == 0 and rhs_mismatch == 0
print('CERTIFICATE VALID' if ok else 'NOT VALID (see above)')
