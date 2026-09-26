"""Independent verification (lab 'structures', own code) of C1's certificate data/x19_cert_s1_c1.json (notes/117 draft):
block-level (P2) infeasibility of one maximal-torsion configuration (tau, sigma) at a vertex x of srg(99,14,1,2).
Checks, in exact rationals:
  1. specials of every pair z (2 tau, 2 sigma) recomputed from tau, sigma;
  2. the set of TYPES of every pair recomputed by enumerating ACTUAL label graphs G_z (8 edges on the 10 labels
     outside the classes of z, no mate edge, simple, label degrees 2 - h_z), and compared with C1's list
     (completeness of the columns is what the certificate needs);
  3. every row rebuilt from its name with own formulas: conv, bal (block balance), lab (label balance = 8),
     P2 (sum_z cnt cnt + N(R,R') = 12|R cap R'| + sum_{y in R, w in R'} (2 - |y cap w|));
  4. b.y - sum over all columns of max(0, col.y) > 0."""
import itertools
import json
import sys
from fractions import Fraction

path = sys.argv[1] if len(sys.argv) > 1 else r'C:\Users\RESEARCHER\Desktop\c1_experimental\data\x19_cert_s1_c1.json'
D = json.load(open(path))
lab = lambda s: tuple(int(t) for t in s.split('.'))
LABELS = [(c, s) for c in range(7) for s in (0, 1)]
mate = lambda u: (u[0], 1 - u[1])

tau = {}
for u, lst in D['tau'].items():
    m = {}
    for v, w in lst:
        m[lab(v)] = lab(w); m[lab(w)] = lab(v)
    tau[lab(u)] = m
sigma = {int(c): {lab(v): lab(w) for v, w in lst} for c, lst in D['sigma'].items()}
sinv = {c: {w: v for v, w in m.items()} for c, m in sigma.items()}
PAIRS = [frozenset({(c, s), (d, t)}) for c, d in itertools.combinations(range(7), 2) for s in (0, 1) for t in (0, 1)]
blk = lambda z: frozenset(l[0] for l in z)
BLOCKS = [frozenset(b) for b in itertools.combinations(range(7), 2)]

# 1. specials
spec = {}
for z in PAIRS:
    e, f = sorted(z)
    S = [frozenset({e, tau[e][f]}), frozenset({f, tau[f][e]})]
    for u, v in ((e, f), (f, e)):
        c = u[0]
        w = sigma[c][v] if u[1] == 0 else sinv[c][v]
        S.append(frozenset({mate(u), w}))
    spec[z] = S
sym_ok = all(z in spec[w] for z in PAIRS for w in spec[z])


# 2. types by enumerating label graphs
def types_of(z):
    P = blk(z)
    others = [c for c in range(7) if c not in P]
    labs = [(c, s) for c in others for s in (0, 1)]
    h = {u: sum(1 for w in spec[z] if u in w) for u in labs}
    need = {u: 2 - h[u] for u in labs}
    if min(need.values()) < 0:
        return set()
    cand = [frozenset({a, b}) for a, b in itertools.combinations(labs, 2) if a[0] != b[0]]
    found = set()

    def rec(i, chosen, need):
        if sum(need.values()) == 0:
            cnts = {}
            for e in chosen:
                R = frozenset(l[0] for l in e)
                cnts[R] = cnts.get(R, 0) + 1
            found.add(frozenset(cnts.items()))
            return
        if i == len(cand):
            return
        # prune: remaining degree must be coverable
        a, b = tuple(cand[i])
        if need[a] > 0 and need[b] > 0:
            need[a] -= 1; need[b] -= 1
            rec(i + 1, chosen + [cand[i]], need)
            need[a] += 1; need[b] += 1
        # skip edge i only if the labels can still be satisfied later
        rec(i + 1, chosen, need)
    rec(0, [], dict(need))
    return found


def types_c1(zname):
    return {frozenset((frozenset(int(ch) for ch in R), m) for R, m in t.items()) for t in D['types_per_pair'][zname]}


pn = lambda z: '-'.join(f'{l[0]}.{l[1]}' for l in sorted(z))
T = {}
mismatch = 0
for z in PAIRS:
    mine = types_of(z)
    theirs = types_c1(pn(z))
    if mine != theirs:
        mismatch += 1
    T[z] = sorted(mine, key=lambda t: sorted((tuple(sorted(R)), m) for R, m in t))
print(f'specials symmetric: {sym_ok}; types: {sum(len(v) for v in T.values())} columns; pairs whose type set differs '
      f'from C1\'s list: {mismatch}')


# 3. rows from names, with C1's convention: a pair whose contribution to a row is CONSTANT (all its counts in the
#    row are near counts, fixed by the specials) goes to the right-hand side; otherwise its FULL contribution goes on
#    each of its type columns.
def cnt_const(z, R):
    """(is_constant, value or None): near block -> constant count of specials in R; disjoint block -> variable."""
    if R & blk(z):
        return True, sum(1 for w in spec[z] if blk(w) == R)
    return False, None


cols = [(z, t) for z in PAIRS for t in T[z]]
colidx = {}
for c, (z, t) in enumerate(cols):
    colidx.setdefault(z, []).append(c)
Bn = lambda s: frozenset(int(ch) for ch in s)
by = Fraction(0)
colval = [Fraction(0)] * len(cols)
rhs_mismatch = 0
for rowd, ys in zip(D['rows'], D['y']):
    y = Fraction(ys)
    name = rowd['name']
    parts = name.split('|')
    kind = parts[0]
    c0 = 0
    coef = {}

    def add_linear(z, R, sign):
        global c0
        isc, v = cnt_const(z, R)
        if isc:
            return sign * v
        for c in colidx[z]:
            a = dict(cols[c][1]).get(R, 0)
            if a:
                coef[c] = coef.get(c, 0) + sign * a
        return 0
    if kind == 'conv':
        z = next(z for z in PAIRS if pn(z) == parts[1])
        for c in colidx[z]:
            coef[c] = 1
        r = 1
    elif kind == 'bal':
        R, Rp = Bn(parts[1]), Bn(parts[2])
        for z in PAIRS:
            if blk(z) == R:
                c0 += add_linear(z, Rp, 1)
            elif blk(z) == Rp:
                c0 += add_linear(z, R, -1)
        r = 0
    elif kind == 'lab':
        u, R = lab(parts[1]), Bn(parts[2])
        for z in PAIRS:
            if u in z:
                c0 += add_linear(z, R, 1)
        r = 8
    elif kind == 'P2':
        R, Rp = Bn(parts[1]), Bn(parts[2])
        for z in PAIRS:
            ia, va = cnt_const(z, R)
            ib, vb = cnt_const(z, Rp)
            if ia and ib:
                c0 += va * vb
            else:
                for c in colidx[z]:
                    t = dict(cols[c][1])
                    av = va if ia else t.get(R, 0)
                    bv = vb if ib else t.get(Rp, 0)
                    if av * bv:
                        coef[c] = coef.get(c, 0) + av * bv
        for yy in PAIRS:
            if blk(yy) == R:
                c0 += add_linear(yy, Rp, 1)
        r = 12 * (4 if R == Rp else 0) + sum(2 - len(yy & ww) for yy in PAIRS if blk(yy) == R
                                             for ww in PAIRS if blk(ww) == Rp and ww != yy)
    else:
        raise ValueError(name)
    b = r - c0
    if b != rowd['rhs']:
        rhs_mismatch += 1
    by += y * b
    for c, a in coef.items():
        colval[c] += y * a
pos = sum(v for v in colval if v > 0)
per_pair = sum(max(colval[c] for c in colidx[z]) for z in PAIRS)
print(f'rows: {len(D["rows"])}; right-hand sides recomputed independently that differ from C1\'s: {rhs_mismatch}')
print(f'b.y = {float(by):.4f};  sum_c max(0, col.y) = {float(pos):.4f};  b.y - that = {float(by - pos):.4f}  '
      f'(C1 claims {float(Fraction(D["farkas_value"])):.4f}; equal: {by - pos == Fraction(D["farkas_value"])})')
print(f'stronger bound with sum_t x = 1: b.y - sum_z max_t col.y = {float(by - per_pair):.4f}')
print('CERTIFICATE VALID' if (by - pos > 0 and mismatch == 0 and sym_ok and rhs_mismatch == 0) else 'NOT VALID (see above)')
