"""Route (C), structures 27/09: explicit rational generators in Q^24 for Borcherds' 156 odd unimodular 24-dim lattices with
min >= 2 (the envelopes U), for C1's trap tool (conway99-c1 code/envelope_trap.py).

Method (Kneser 2-neighbours, everything stays in Q^24 with the standard dot product):
- every such U has two even unimodular neighbours N0, N0' (Niemeier) and U = N0_v + (odd coset), where
  N0_v = {x in N0 : (x,v) even}, v in N0, and the two cosets of N0_v other than N0 make up v/2 + N0;
- min U >= 2  <=>  v/2 + N0 has no norm-1 vector  <=>  the class v + 2N0 has no norm-4 vector;
- WLOG v is a minimal vector of its class; then |(v, alpha)| <= 2 for every root (else v -+ 2 alpha is shorter), so up to
  the Weyl group v = sum over the root components of dominant weights lambda_c with (lambda_c, theta_c) <= 2;
- roots(U) = {alpha in roots(N0) : (alpha, v) even}; roots(N0') = those plus the norm-2 vectors of v/2 + N0 (one CVP).
Breadth-first from D24+ over the Niemeier lattices.  Each U is identified by (roots(U), N0, N0'), the key of
data/borcherds_odd24.json (156 distinct keys).  Every U built is re-verified from its own basis: integral Gram, det 1,
odd, no norm-1 vector, root system = the table's, number of roots = the table's.
Usage: python build156.py [max_seconds]      (resumable; state in state/, one JSON per envelope in U/)"""
import collections, json, math, os, random, re, sys, time
from fractions import Fraction as F
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = 'C:/home/researcher/Desktop/c1_recherche/data/borcherds_odd24.json'
STATE = os.path.join(HERE, 'state')
UDIR = os.path.join(HERE, 'U')
T0 = time.time()
BUDGET = 270.0                                   # set from argv in main() (importable module)


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


# ------------------------------------------------------------------ root-system names
TOK = re.compile(r'([ADEade])_(?:\{(\d+)\}|(\d))(?:\^(?:\{(\d+)\}|(\d)))?')


def parse_type(s):
    c = collections.Counter()
    for m in TOK.finditer(s):
        t, n, e = m.group(1).lower(), int(m.group(2) or m.group(3)), int(m.group(4) or m.group(5) or 1)
        if t == 'd' and n == 3:
            t = 'a'
        c[(t, n)] += e
    return c


def tkey(c):
    return tuple(sorted((k, v) for k, v in c.items() if v))


def tname(key):
    f = lambda n: str(n) if n < 10 else '{%d}' % n
    return ''.join('%s_%s%s' % (t, f(n), '' if e == 1 else '^' + f(e)) for (t, n), e in key) or '(none)'


def type_of(count, rank):
    if count == rank * (rank + 1):
        return ('a', rank)
    if count == 2 * rank * (rank - 1):
        return ('d', rank)
    if (rank, count) in ((6, 72), (7, 126), (8, 240)):
        return ('e', rank)
    raise ValueError((count, rank))


# ------------------------------------------------------------------ exact integer linear algebra
def lcm(a, b):
    return a * b // math.gcd(a, b)


def int_rows(rows):
    d = 1
    for r in rows:
        for x in r:
            d = lcm(d, F(x).denominator)
    return [[int(F(x) * d) for x in r] for r in rows], d


def hnf(rows):
    """a Z-basis (echelon rows) of the span of integer rows"""
    A = [list(r) for r in rows if any(r)]
    n = len(A[0])
    out, col = [], 0
    while A and col < n:
        nz = [r for r in A if r[col]]
        if not nz:
            col += 1
            continue
        while len(nz) > 1:
            nz.sort(key=lambda r: abs(r[col]))
            p = nz[0]
            for r in nz[1:]:
                q = r[col] // p[col]
                for k in range(col, n):
                    r[k] -= q * p[k]
            nz = [r for r in nz if r[col]]
        p = nz[0]
        if p[col] < 0:
            p[:] = [-x for x in p]
        out.append(p)
        A = [r for r in A if r is not p and any(r)]
        col += 1
    return out


def gram_of(M, d):
    n, d2 = len(M), d * d
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            s = sum(a * b for a, b in zip(M[i], M[j]))
            assert s % d2 == 0, 'Gram not integral'
            G[i][j] = G[j][i] = s // d2
    return G


def det_exact(G):
    """Bareiss"""
    M = [list(map(int, r)) for r in G]
    n, sign, prev = len(M), 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k]), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def inv_exact(M):
    n = len(M)
    A = [[F(x) for x in row] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


def lll(M, d, delta=0.99):
    """LLL of the basis rows M/d (exact integer Gram updates, float Gram-Schmidt); returns (M', G')"""
    n = len(M)
    G = gram_of(M, d)
    T = [[int(i == j) for j in range(n)] for i in range(n)]
    chol = lambda: np.linalg.cholesky(np.array(G, dtype=float))
    L = chol()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = int(round(L[k, j] / L[j, j]))
            if q:
                T[k] = [a - q * b for a, b in zip(T[k], T[j])]
                G[k] = [a - q * b for a, b in zip(G[k], G[j])]
                for l in range(n):
                    G[l][k] -= q * G[l][j]
                L = chol()
        if L[k, k] ** 2 >= (delta - (L[k, k - 1] / L[k - 1, k - 1]) ** 2) * L[k - 1, k - 1] ** 2:
            k += 1
        else:
            T[k], T[k - 1] = T[k - 1], T[k]
            G[k], G[k - 1] = G[k - 1], G[k]
            for row in G:
                row[k], row[k - 1] = row[k - 1], row[k]
            L = chol()
            k = max(k - 1, 1)
    M2 = [[sum(T[i][a] * M[a][c] for a in range(n)) for c in range(len(M[0]))] for i in range(n)]
    return M2, G


class Lattice:
    def __init__(self, gens):
        """gens: rational rows in Q^24; builds an LLL-reduced basis"""
        Mi, d = int_rows(gens)
        B = hnf(Mi)
        assert len(B) == 24, len(B)
        self.M, self.G = lll(B, d)
        self.d = d
        self.Gm = np.array(self.G, dtype=np.int64)
        self.Ginv = None

    def rows(self):
        return [[F(x, self.d) for x in r] for r in self.M]

    def coords(self, coef):
        """coefficient vector (rational) -> coordinates in Q^24"""
        return [sum(F(coef[i]) * self.M[i][c] for i in range(24)) / self.d for c in range(24)]


def qform(Gm, w):
    w = np.array(w, dtype=np.int64)
    return int(w @ Gm @ w)


def enum(G, bound, cnum=None, cden=1):
    """Fincke-Pohst: integer z with (z - c)^T G (z - c) <= bound, c = cnum/cden.  c = 0: nonzero z, one per +-pair.
    Returns [(z, exact value)]."""
    n = len(G)
    Gm = np.array(G, dtype=np.int64)
    R = np.linalg.cholesky(np.array(G, dtype=float)).T
    d2 = [R[i, i] ** 2 for i in range(n)]
    off = [[R[i, j] / R[i, i] for j in range(n)] for i in range(n)]
    c = [0.0] * n if cnum is None else [x / cden for x in cnum]
    z = [0] * n
    raw = []
    eps = 1e-6

    def rec(i, rem):
        ctr = c[i] - sum(off[i][j] * (z[j] - c[j]) for j in range(i + 1, n))
        r = math.sqrt(max(rem, 0.0) / d2[i])
        for zi in range(math.ceil(ctr - r - eps), math.floor(ctr + r + eps) + 1):
            val = d2[i] * (zi - ctr) ** 2
            if val > rem + eps:
                continue
            z[i] = zi
            if i == 0:
                raw.append(tuple(z))
            else:
                rec(i - 1, rem - val)
        z[i] = 0
    rec(n - 1, float(bound) + eps)
    out = []
    for zz in raw:
        if cnum is None:
            if not any(zz) or next(x for x in reversed(zz) if x) < 0:
                continue
            val = F(qform(Gm, zz))
        else:
            val = F(qform(Gm, [cden * a - b for a, b in zip(zz, cnum)]), cden * cden)
        if val <= bound:
            out.append((zz, val))
    return out


# ------------------------------------------------------------------ root systems
def comps_of(V, Gm):
    """connected components of a root set (rows of V, one per +-pair) under nonzero inner products"""
    if len(V) == 0:
        return []
    A = (V @ Gm @ V.T) != 0
    seen, out = [False] * len(V), []
    for s in range(len(V)):
        if seen[s]:
            continue
        stack, comp = [s], []
        seen[s] = True
        while stack:
            i = stack.pop()
            comp.append(i)
            for j in np.nonzero(A[i])[0]:
                if not seen[j]:
                    seen[j] = True
                    stack.append(j)
        out.append(comp)
    return out


def rtype(V, Gm):
    c = collections.Counter()
    for comp in comps_of(V, Gm):
        rank = int(np.linalg.matrix_rank(V[comp].astype(float)))
        c[type_of(2 * len(comp), rank)] += 1
    return c


def root_structure(roots, Gm):
    """roots: int array (one per +-pair).  Per component: type, simple roots S, Cartan C and inverse, marks, positive roots
    with their simple-root coordinates."""
    V = np.array(roots, dtype=np.int64)
    rng = random.Random(7)
    while True:
        f = np.array([rng.randint(-10 ** 6, 10 ** 6) for _ in range(V.shape[1])], dtype=np.int64)
        h = V @ f
        if (h != 0).all():
            break
    out = []
    for comp in comps_of(V, Gm):
        pos = [V[i] if h[i] > 0 else -V[i] for i in comp]
        hp = [abs(int(h[i])) for i in comp]
        order = sorted(range(len(pos)), key=lambda t: hp[t])
        pset = {tuple(p) for p in pos}
        S = []
        for t in order:
            a = pos[t]
            if not any(tuple(a - s) in pset for s in S):
                S.append(a)
        S = np.array(S, dtype=np.int64)
        C = (S @ Gm @ S.T).tolist()
        Ci = inv_exact(C)
        typ = type_of(2 * len(pos), len(S))
        SG = S @ Gm
        ncoord = []
        for a in pos:
            w = (SG @ a).tolist()
            nc = [sum(Ci[i][j] * w[j] for j in range(len(S))) for i in range(len(S))]
            assert all(x.denominator == 1 for x in nc)
            ncoord.append([int(x) for x in nc])
        theta = ncoord[order[-1]]
        assert all(x >= 1 for x in theta), theta
        out.append(dict(typ=typ, S=S, C=C, Ci=Ci, marks=theta, pos=np.array(pos), ncoord=np.array(ncoord)))
    return out


def component_options(comp, Gm, maxnorm):
    """dominant weights lambda with (lambda, theta) <= 2 and norm <= maxnorm: coefficient vector, norm, pairing with the
    basis (mod 1), kernel root type, label"""
    r = len(comp['S'])
    mk = comp['marks']
    avecs = [[0] * r]
    for i in range(r):
        if mk[i] in (1, 2):
            avecs.append([int(j == i) for j in range(r)])
    ones = [i for i in range(r) if mk[i] == 1]
    for x in range(len(ones)):
        for y in range(x, len(ones)):
            a = [0] * r
            a[ones[x]] += 1
            a[ones[y]] += 1
            avecs.append(a)
    Ci, S = comp['Ci'], comp['S']
    opts = []
    for a in avecs:
        norm = sum(a[i] * Ci[i][j] * a[j] for i in range(r) for j in range(r))
        if norm > maxnorm:
            continue
        w = [sum(a[i] * Ci[i][j] for i in range(r)) for j in range(r)]          # lambda = sum_j w_j s_j
        lam = [sum(w[j] * int(S[j][k]) for j in range(r)) for k in range(24)]  # coefficient vector (rational)
        pair = [sum(lam[i] * int(Gm[i][k]) for i in range(24)) % 1 for k in range(24)]
        par = (comp['ncoord'] @ np.array(a, dtype=np.int64)) % 2
        ker = comp['pos'][par == 0]
        kt = rtype(ker, Gm) if len(ker) else collections.Counter()
        opts.append(dict(a=a, norm=norm, lam=lam, pair=pair, ktype=kt, label=(comp['typ'], norm, tkey(kt))))
    return opts


# ------------------------------------------------------------------ the neighbour sampler for one Niemeier lattice
class Sampler:
    TARGET_NORMS = (8, 12, 16, 20, 24)

    def __init__(self, comps_opts):
        self.opts = comps_opts
        allopts = [o for co in comps_opts for o in co]
        self.Dn = 1
        for o in allopts:
            self.Dn = lcm(self.Dn, o['norm'].denominator)
        Dp = 1
        for o in allopts:
            for x in o['pair']:
                Dp = lcm(Dp, x.denominator)
        cols = [k for k in range(24) if any(o['pair'][k] for o in allopts)]
        seen, keep = set(), []
        for k in cols:
            sig = tuple(o['pair'][k] for o in allopts)
            if sig not in seen:
                seen.add(sig)
                keep.append(k)
        for o in allopts:
            o['p'] = tuple(int(o['pair'][k] * Dp) % Dp for k in keep)
            o['n'] = int(o['norm'] * self.Dn)
        zero = tuple([0] * len(keep))
        # the finite group of reachable pairings
        idx, elems, frontier = {zero: 0}, [zero], [zero]
        gens = {o['p'] for o in allopts}
        while frontier:
            nf = []
            for g in frontier:
                for p in gens:
                    h = tuple((a + b) % Dp for a, b in zip(g, p))
                    if h not in idx:
                        idx[h] = len(elems)
                        elems.append(h)
                        nf.append(h)
            frontier = nf
        self.ng = len(elems)
        for o in allopts:
            o['add'] = np.array([idx[tuple((a + b) % Dp for a, b in zip(g, o['p']))] for g in elems], dtype=np.int64)
            inv = np.empty(self.ng, dtype=np.int64)
            inv[o['add']] = np.arange(self.ng)
            o['sub'] = inv
        self.maxn = max(self.TARGET_NORMS) * self.Dn
        lay = np.zeros((self.maxn + 1, self.ng))
        lay[0, 0] = 1.0
        self.layers = [lay]
        for co in comps_opts:
            new = np.zeros_like(lay)
            for o in co:
                n = o['n']
                if n > self.maxn:
                    continue
                new[n:, :][:, o['add']] += lay[:self.maxn + 1 - n, :]
            self.layers.append(new)
            lay = new
        self.total = {T: lay[T * self.Dn, 0] for T in self.TARGET_NORMS}

    def sample(self, T, rng):
        nu, g = T * self.Dn, 0
        pick = [None] * len(self.opts)
        for c in range(len(self.opts) - 1, -1, -1):
            prev = self.layers[c]
            ws = []
            for o in self.opts[c]:
                n = o['n']
                ws.append(prev[nu - n, o['sub'][g]] if n <= nu else 0.0)
            tot = sum(ws)
            x = rng.random() * tot
            k = 0
            while k < len(ws) - 1 and x >= ws[k]:
                x -= ws[k]
                k += 1
            o = self.opts[c][k]
            pick[c] = k
            nu -= o['n']
            g = int(o['sub'][g])
        assert nu == 0 and g == 0
        return pick


# ------------------------------------------------------------------ state
def load_json(p, default):
    return json.load(open(p)) if os.path.exists(p) else default


def save_json(p, obj):
    tmp = p + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f)
    for k in range(20):                     # Windows: a scanner may hold the file for a moment
        try:
            os.replace(tmp, p)
            return
        except PermissionError:
            time.sleep(0.25)
    with open(p, 'w') as f:
        json.dump(obj, f)


def rows_str(rows):
    return [[str(x) for x in r] for r in rows]


def rows_frac(rows):
    return [[F(x) for x in r] for r in rows]


def d24plus():
    gens = []
    for i in range(23):
        v = [F(0)] * 24
        v[i], v[i + 1] = F(1), F(-1)
        gens.append(v)
    v = [F(0)] * 24
    v[22], v[23] = F(1), F(1)
    gens.append(v)
    gens.append([F(1, 2)] * 24)
    return gens


def verify_U(rows, target):
    """independent check of an envelope from its basis rows (rational, Q^24)"""
    Mi, d = int_rows(rows)
    G = gram_of(Mi, d)
    det = det_exact(G)
    odd = any(G[i][i] % 2 for i in range(24))
    sv = enum(G, 2)
    n1 = sum(1 for z, v in sv if v == 1)
    R = np.array([z for z, v in sv if v == 2], dtype=np.int64)
    rt = tkey(rtype(R, np.array(G, dtype=np.int64))) if len(R) else ()
    ok = det == 1 and odd and n1 == 0 and 2 * len(R) == target['nroots'] and rt == tkey(parse_type(target['roots']))
    return ok, dict(det=det, odd=odd, norm1=n1, roots=2 * len(R), rtype=tname(rt))


# Transcription fixes to data/borcherds_odd24.json, each with a proof.
# 150: listed as (D10E7^2, D16E8).  Impossible: roots(U) lie in both even neighbours and e_7^2 is not a root subsystem of
# D16 + E8.  E8^3 (same Coxeter number 30, hence the same root count 8(18 + 30 - 2) = 368) has exactly three admissible v
# (w1 w1 0 -> 152, w1 w1 w1 -> 146, w1 w8 w8 -> d8 e7^2 a1^2), so the second neighbour of 150 is E8^3.
TABLE_FIXES = {150: {'nb2': 'E_8^3'}}


def sub_multiset(a, b):
    return all(b[k] >= v for k, v in a.items())


def try_v(N0, name, roots, ch, T, need, st, niem_by_type, tried):
    """one candidate v = sum of the chosen dominant weights: CVP, identification, construction, verification"""
    ut = collections.Counter()
    for o in ch:
        ut.update(o['ktype'])
    ut = tkey(ut)
    if ut not in need:
        return
    lab = (T, tuple(sorted(o['label'] for o in ch if o['a'] != [0] * len(o['a']))))
    if tried[lab] >= 3:
        return
    tried[lab] += 1
    cv = [sum(o['lam'][i] for o in ch) for i in range(24)]
    assert all(x.denominator == 1 for x in cv)
    cv = [int(x) for x in cv]
    if all(x % 2 == 0 for x in cv):
        return
    near = enum(N0.G, 2, cnum=[-x for x in cv], cden=2)
    if any(v < 2 for z, v in near):
        return False                              # class has a norm-4 (or 0) vector: U would have norm 1
    new = [[2 * a + b for a, b in zip(z, cv)] for z, v in near if v == 2]
    keep = []
    seenpm = set()
    for w in new:
        t = tuple(w)
        if tuple(-x for x in t) not in seenpm:
            seenpm.add(t)
            keep.append(w)
    par = (np.array(roots, dtype=np.int64) @ N0.Gm @ np.array(cv, dtype=np.int64)) % 2
    kr = [[2 * x for x in roots[i]] for i in range(len(roots)) if par[i] == 0]
    allr = np.array(kr + keep, dtype=np.int64)
    nt = tkey(rtype(allr, N0.Gm))
    other = niem_by_type.get(nt)
    hit = [(o, e) for o, e in need[ut] if o == other]
    log('  T=%d U=%s N0\'=%s %s' % (T, tname(ut), other or tname(nt), '-> id %d' % hit[0][1]['id'] if hit else ''))
    # the two cosets: A = v/2 + N0_v (odd iff T = 4 mod 8), B = A + b_k0
    a = (N0.Gm @ np.array(cv, dtype=np.int64)) % 2
    k0 = int(np.nonzero(a)[0][0])
    base = []
    for j in range(24):
        e_ = [0] * 24
        if j == k0:
            e_[j] = 2
        else:
            e_[j] = 1
            if a[j]:
                e_[k0] = 1
        base.append(e_)
    half = [F(x, 2) for x in cv]
    halfB = list(half)
    halfB[k0] += 1
    glueU, glueN = (half, halfB) if T % 8 == 4 else (halfB, half)
    if other and other not in st['niemeier']:
        Ln = Lattice([N0.coords(b) for b in base] + [N0.coords(glueN)])
        assert all(Ln.G[i][i] % 2 == 0 for i in range(24)) and det_exact(Ln.G) == 1
        st['niemeier'][other] = rows_str(Ln.rows())
        if other not in st['queue'] and other not in st['processed']:
            st['queue'].append(other)
        log('  new Niemeier lattice %s stored' % other)
    if hit:
        o_, e = hit[0]
        LU = Lattice([N0.coords(b) for b in base] + [N0.coords(glueU)])
        ok, info = verify_U(LU.rows(), e)
        log('  id %d %s: %s' % (e['id'], 'VERIFIED' if ok else 'FAILED', info))
        if ok:
            save_json(os.path.join(UDIR, '%03d.json' % e['id']),
                      {'id': e['id'], 'roots': e['roots'], 'nb1': e['nb1'], 'nb2': e['nb2'], 'nroots': e['nroots'],
                       'built_from': name, 'v_norm': T, 'check': info,
                       'table_fix': {k: v for k, v in e.items() if k.startswith('fix_')}, 'gens': rows_str(LU.rows())})
            st['built'][str(e['id'])] = name
            need[ut] = [(o, x) for o, x in need[ut] if x is not e]
            if not need[ut]:
                del need[ut]
    return True


def main():
    global BUDGET
    if len(sys.argv) > 1:
        BUDGET = float(sys.argv[1])
    table = json.load(open(TABLE))
    for e in table:
        for k, v in TABLE_FIXES.get(e['id'], {}).items():
            e['fix_' + k] = e[k]
            e[k] = v
    bykey = {}
    for e in table:
        bykey[(tkey(parse_type(e['roots'])), e['nb1'], e['nb2'])] = e
        bykey[(tkey(parse_type(e['roots'])), e['nb2'], e['nb1'])] = e
    niem_names = sorted({e['nb1'] for e in table} | {e['nb2'] for e in table})
    niem_by_type = {tkey(parse_type(nm)): nm for nm in niem_names}
    st_path = os.path.join(STATE, 'state.json')
    st = load_json(st_path, {'niemeier': {}, 'processed': [], 'built': {}, 'queue': ['D_{24}']})
    if 'D_{24}' not in st['niemeier']:
        L = Lattice(d24plus())
        st['niemeier']['D_{24}'] = rows_str(L.rows())
        save_json(st_path, st)
    rng = random.Random(2027)
    while st['queue'] and time.time() - T0 < BUDGET:
        name = st['queue'][0]
        targets = [e for e in table if name in (e['nb1'], e['nb2']) and str(e['id']) not in st['built']]
        log('== %s: %d targets left %s' % (name, len(targets), [e['id'] for e in targets]))
        if not targets:
            st['queue'].pop(0)
            st['processed'].append(name)
            save_json(st_path, st)
            continue
        N0 = Lattice(rows_frac(st['niemeier'][name]))
        assert all(N0.G[i][i] % 2 == 0 for i in range(24)) and det_exact(N0.G) == 1
        roots = [z for z, v in enum(N0.G, 2) if v == 2]
        comps = root_structure(roots, N0.Gm)
        got = tkey(collections.Counter(c['typ'] for c in comps))
        assert got == tkey(parse_type(name)), (got, name)
        log('%s: %d roots, components %s' % (name, 2 * len(roots), tname(got)))
        opts = [component_options(c, N0.Gm, 24) for c in comps]
        smp = Sampler(opts)
        log('options per component %s, glue group %d, solutions per norm %s' % (
            [len(o) for o in opts], smp.ng, {T: '%.3g' % smp.total[T] for T in smp.TARGET_NORMS}))
        need = collections.defaultdict(list)
        for e in targets:
            other = e['nb2'] if e['nb1'] == name else e['nb1']
            need[tkey(parse_type(e['roots']))].append((other, e))
        tried = collections.Counter()
        norms = [T for T in smp.TARGET_NORMS if smp.total[T] > 0]
        nsamp = 0
        while need and time.time() - T0 < BUDGET and nsamp < 200000:
            T = rng.choice(norms)
            pick = smp.sample(T, rng)
            nsamp += 1
            if try_v(N0, name, roots, [opts[c][k] for c, k in enumerate(pick)], T, need, st, niem_by_type, tried):
                save_json(st_path, st)
        for ut in list(need):
            if time.time() - T0 > BUDGET:
                break
            utc = collections.Counter(dict(ut))
            fopts = [[dict(o) for o in co if sub_multiset(o['ktype'], utc)] for co in opts]
            if any(not f for f in fopts):
                continue
            s2 = Sampler(fopts)
            n2 = [T for T in s2.TARGET_NORMS if s2.total[T] > 0]
            log('targeted %s: solutions per norm %s' % (tname(ut), {T: '%.3g' % s2.total[T] for T in n2}))
            for _ in range(300):
                if ut not in need or not n2:
                    break
                T = rng.choice(n2)
                pick = s2.sample(T, rng)
                try_v(N0, name, roots, [fopts[c][k] for c, k in enumerate(pick)], T, need, st, niem_by_type, tried)
            save_json(st_path, st)
        if not need:
            st['queue'].pop(0)
            st['processed'].append(name)
        else:
            log('%s: still missing %s after %d samples' % (name, [(tname(k), [(o, e['id']) for o, e in v]) for k, v in need.items()],
                                                         nsamp))
            if nsamp >= 200000:
                st['queue'].append(st['queue'].pop(0))
        save_json(st_path, st)
    log('built %d/156; queue %s; known Niemeier %d' % (len(st['built']), st['queue'], len(st['niemeier'])))


if __name__ == '__main__':
    main()
