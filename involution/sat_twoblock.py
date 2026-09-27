"""Two-block SAT for the T = K7 case of the Sigma = empty 'one-block' theorem (structures, 27/09).

Two local systems (encoding 2, sat_local2.py) at two blocks c and d of N(x), glued along everything they both see.
World frame: N(x) vertices 2k (+) and 2k+1 (-) of world block k; c = world block 0, d = world block 1, sigma = (v -> v^1).
  c-frame  : the labels of the c-system (canonical type rep); c-frame block p is d, the others go to world blocks 2..m-1
             in increasing order, orientation kept (Fc).
  d-frame  : the labels of the d-system (canonical type rep of d); G: d-frame label -> world vertex of N(x) - d,
             G(gamma) = c+ (world 0), G(pi gamma) = c- (world 1); G on the other blocks is either FIXED (--glue) or a
             SAT variable (signed bijection, 50 variables B[beta, k, o] at m = 7).
Every exterior vertex is a pair of world N(x) vertices: c-roots contain c+-, d-roots d+-, the 4 shared ones both, the 40
far-far ones neither. The c-system describes adj(X, Y) whenever X or Y is a c-root (c-root x c-root: the type (mu_c, h_c);
c-root x far-at-c: U_c), the d-system whenever X or Y is a d-root.
  COUPLING : for every pair described by BOTH systems, the two descriptions are equal. That is, mod sigma:
             c-root x d-root 200 (C1/C2's U_c[(d+,w),v] <-> U_d[(c+,v),w] and its c-/d+ twin), c-root x shared 40
             and shared x d-root 40 (the units from the other type), shared x far-far 80 (U_c[Y, d+-] <-> U_d[Y, c+-]),
             shared x shared 3 (type compatibility).
  CROSS    : the two-block (iii). For a non-shared c-root r and a non-shared d-root s, EVERY common neighbour is modelled
             (a common neighbour is adjacent to a c-root and to a d-root), so
               [r and s share their N(x) vertex] + sum_X adj(r,X) adj(s,X) + adj(r,s) = 2       (lambda = 1, mu = 2),
             200 exact equations mod sigma.
Local systems: sat_local2.encode at c (c-frame) and at d (d-frame), with the flags given (--sec --twist for T = K7).
Scope: T = K7 = every pair twisted, every block of type (0,6), Sigma = empty (SEC valid). Pivot consistency on the pair cd:
[d+ pivot at c] = [c- pivot at d] (the orbit that pivots at c does not pivot at d).
Usage:
  python sat_twoblock.py --list                               writes twoblock_instances.json (60 instances, see make_list)
  python sat_twoblock.py twoblock_instances.json K [--dimacs] [--nocross] [--sec] [--twist] [--glue=l0,...,l11] [secs]
  python sat_twoblock.py --selftest                           BvLS (real graph, m = 11) + synthetic m = 7 gluing tests"""
import itertools, json, sys, time, random
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType


def neg(l):
    return (not l) if isinstance(l, bool) else -l


class Builder:
    """CNF with constant folding; literals are nonzero ints or the constants True/False. When every provider returns
    constants (evaluation mode) no variable is created and each violated constraint is recorded in .violations."""

    def __init__(self, start=1):
        self.pool = IDPool(start_from=start)
        self.cnf = CNF()
        self.violations = []
        self.memo = {}

    def new(self, key):
        assert key not in self.pool.obj2id, key
        return self.pool.id(key)

    def clause(self, lits, tag=None):
        out = []
        for l in lits:
            if l is True:
                return
            if l is False:
                continue
            out.append(l)
        if not out:
            self.violations.append(tag)
        self.cnf.append(out)

    def imp(self, conds, lits, tag=None):
        self.clause([neg(c) for c in conds] + list(lits), tag)

    def equiv(self, a, b, conds=(), tag=None):
        self.imp(conds, [neg(a), b], tag)
        self.imp(conds, [a, neg(b)], tag)

    def AND(self, a, b):
        if a is False or b is False:
            return False
        if a is True:
            return b
        if b is True:
            return a
        key = ('and', min(a, b), max(a, b))
        if key not in self.memo:
            v = self.new(key)
            self.cnf.extend([[-v, a], [-v, b], [v, -a, -b]])
            self.memo[key] = v
        return self.memo[key]

    def OR(self, a, b):
        if a is True or b is True:
            return True
        if a is False:
            return b
        if b is False:
            return a
        key = ('or', min(a, b), max(a, b))
        if key not in self.memo:
            v = self.new(key)
            self.cnf.extend([[v, -a], [v, -b], [-v, a, b]])
            self.memo[key] = v
        return self.memo[key]

    def exactly(self, lits, k, tag=None):
        rest, t = [], 0
        for l in lits:
            if l is True:
                t += 1
            elif l is not False:
                rest.append(l)
        k -= t
        if k < 0 or k > len(rest):
            self.violations.append(tag)
            self.cnf.append([])
            return
        if not rest:
            return
        if len(set(rest)) < len(rest):                  # repeated literal: separate copies, so it counts twice
            seen, cp = set(), []
            for l in rest:
                if l in seen:
                    c = self.new(('copy', len(self.memo), l)); self.memo[('copy', len(self.memo), l)] = c
                    self.cnf.extend([[-c, l], [c, -l]]); cp.append(c)
                else:
                    seen.add(l); cp.append(l)
            rest = cp
        if k == 0:
            self.cnf.extend([[-l] for l in rest])
        elif k == len(rest):
            self.cnf.extend([[l] for l in rest])
        else:
            self.cnf.extend(CardEnc.equals(lits=rest, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses)


def sig(Z):
    return tuple(sorted((Z[0] ^ 1, Z[1] ^ 1)))


def root_parts(X, blk):
    """(vertex of world block blk in X, the other vertex), or None"""
    a, b = X
    if a // 2 == blk:
        return a, b
    if b // 2 == blk:
        return b, a
    return None


class TwoBlock:
    """coupling + cross constraints between the c-system and the d-system (generic m), see the module docstring.
    Uc(zc, b) / Ud(zd, lam): literal providers of the two local systems (U[z,b] = [z ~ (base+, b)], z a sorted label pair)."""

    def __init__(self, B, m, muc, tauc, p, mud, taud, gamma, Uc, Ud, glue=None, full=False, FFval=None):
        self.B, self.m = B, m
        self.full, self.FFval, self.FF = full, FFval, {}
        n = 2 * (m - 1); self.n = n
        self.muc, self.hc = muc, [tauc[l] ^ 1 for l in range(n)]          # tau = pi o h  =>  h = pi o tau
        self.mud, self.hd = mud, [taud[l] ^ 1 for l in range(n)]
        self.p, self.gamma = p, gamma
        others = [q for q in range(m - 1) if q != p]
        cblk = {p: 1}
        cblk.update({q: k + 2 for k, q in enumerate(others)})
        self.Fc = [2 * cblk[l // 2] + (l & 1) for l in range(n)]
        self.Fci = {w: l for l, w in enumerate(self.Fc)}
        self.Uc_, self.Ud_ = Uc, Ud
        self.gb = gamma // 2
        self.dfar = [z for z in itertools.combinations(range(n), 2) if z[0] // 2 != z[1] // 2]
        self.ext = [z for z in itertools.combinations(range(2 * m), 2) if z[0] // 2 != z[1] // 2]
        if glue is not None:
            assert glue[gamma] == 0 and glue[gamma ^ 1] == 1, 'glue must send gamma to c+ and pi gamma to c-'
            assert sorted(glue) == [w for w in range(2 * m) if w // 2 != 1]
            assert all(glue[l ^ 1] == glue[l] ^ 1 for l in range(n))
            self.G = list(glue)
        else:
            self.G = None
            self.Bv = {}
            dbl = [b for b in range(m - 1) if b != self.gb]
            for b in dbl:
                for k in range(2, m):
                    for o in (0, 1):
                        self.Bv[b, k, o] = B.new(('B', b, k, o))
            for b in dbl:
                B.exactly([self.Bv[b, k, o] for k in range(2, m) for o in (0, 1)], 1, ('glue-row', b))
            for k in range(2, m):
                B.exactly([self.Bv[b, k, o] for b in dbl for o in (0, 1)], 1, ('glue-col', k))

    # ---- gluing ----
    def g(self, lam, a):
        """[G(lam) = a] for a d-frame label lam and a world vertex a outside block d"""
        if self.G is not None:
            return self.G[lam] == a
        if lam // 2 == self.gb:
            return a == (0 if lam == self.gamma else 1)
        if a // 2 in (0, 1):
            return False
        return self.Bv[lam // 2, a // 2, (lam ^ a) & 1]

    def pre(self, a):
        out = [(lam, self.g(lam, a)) for lam in range(self.n)]
        return [(lam, gl) for lam, gl in out if gl is not False]

    def _transport(self, key, cands, value):
        """literal equal to value(c) for the unique candidate c whose condition literal is true"""
        if key in self.B.memo:
            return self.B.memo[key]
        cands = [(c, cl) for c, cl in cands if cl is not False]
        assert cands, key
        if len(cands) == 1 and cands[0][1] is True:
            r = value(cands[0][0])
        else:
            r = self.B.new(key)
            for c, cl in cands:
                self.B.equiv(r, value(c), conds=[cl], tag=key)
        self.B.memo[key] = r
        return r

    def M(self, x, y):
        """[mu_d(G^-1 x) = G^-1 y], world vertices x, y outside block d"""
        return self._transport(('M', x, y), self.pre(x), lambda lam: self.g(self.mud[lam], y))

    def H(self, x, y):
        """[h_d(G^-1 x) = G^-1 y]"""
        return self._transport(('H', x, y), self.pre(x), lambda lam: self.g(self.hd[lam], y))

    def Y1(self, zd, x):
        """U_d[zd, G^-1 x] (d-frame far root zd, world label x)"""
        return self._transport(('Y1', zd, x), self.pre(x), lambda lam: self.Ud_(zd, lam))

    def R(self, zd, Y):
        """[G maps the d-frame far root zd onto the world vertex Y]"""
        if ('R', zd, Y) in self.B.memo:
            return self.B.memo[('R', zd, Y)]
        l1, l2 = zd
        y1, y2 = Y
        r = self.B.OR(self.B.AND(self.g(l1, y1), self.g(l2, y2)), self.B.AND(self.g(l1, y2), self.g(l2, y1)))
        self.B.memo[('R', zd, Y)] = r
        return r

    def Wd(self, Y, x):
        """[Y ~ (d+, x)] = U_d[G^-1 Y, G^-1 x], for Y far at d (world) and x a world vertex outside block d"""
        if ('W', Y, x) in self.B.memo:
            return self.B.memo[('W', Y, x)]
        return self._transport(('W', Y, x), [(zd, self.R(zd, Y)) for zd in self.dfar], lambda zd: self.Y1(zd, x))

    # ---- adjacency as described by each system ----
    def descC(self, X, Y):
        px, py = root_parts(X, 0), root_parts(Y, 0)
        if px and py:
            (s, x), (t, y) = px, py
            a, b = self.Fci[x], self.Fci[y]
            if s == 0 and t == 0:
                return self.muc[a] == b
            if s == 0 and t == 1:
                return self.hc[a] == b
            if s == 1 and t == 0:
                return self.hc[b] == a
            return self.muc[a ^ 1] == b ^ 1
        if px or py:
            (s, x), Z = (px, Y) if px else (py, X)
            if s == 1:
                Z, x = sig(Z), x ^ 1
            return self.Uc_(tuple(sorted((self.Fci[Z[0]], self.Fci[Z[1]]))), self.Fci[x])
        return None

    def descD(self, X, Y):
        px, py = root_parts(X, 1), root_parts(Y, 1)
        if px and py:
            (s, x), (t, y) = px, py
            if s == 2 and t == 2:
                return self.M(x, y)
            if s == 2 and t == 3:
                return self.H(x, y)
            if s == 3 and t == 2:
                return self.H(y, x)
            return self.M(x ^ 1, y ^ 1)
        if px or py:
            (s, x), Z = (px, Y) if px else (py, X)
            if s == 3:
                Z, x = sig(Z), x ^ 1
            return self.Wd(Z, x)
        return None

    def descF(self, X, Y):
        """full model only: adjacency of two far-far vertices (a variable per sigma-class; {z, sigma z} never adjacent)"""
        if not self.full or self.kind(X) != 'F' or self.kind(Y) != 'F':
            return None
        if Y == sig(X):
            return False
        if self.FFval is not None:
            return self.FFval(X, Y)
        key = min(tuple(sorted((X, Y))), tuple(sorted((sig(X), sig(Y)))))
        if key not in self.FF:
            self.FF[key] = self.B.new(('FF',) + key)
        return self.FF[key]

    def adj(self, X, Y):
        a = self.descC(X, Y)
        if a is None:
            a = self.descD(X, Y)
        if a is None:
            a = self.descF(X, Y)
        return a

    def full_counts(self):
        """full model: exact lambda/mu counts for every pair not already exact (local systems, CROSS), i.e. every pair
        with a far-far member: (N(x) vertex outside c, d ; far-far), (root ; far-far), (far-far ; far-far)."""
        assert self.full
        F = [X for X in self.ext if self.kind(X) == 'F']
        R = [X for X in self.ext if self.kind(X) != 'F']
        seen, n = set(), 0
        for a in range(4, 2 * self.m):
            for z in F:
                key = min((a, z), (a ^ 1, sig(z)))
                if key in seen:
                    continue
                seen.add(key)
                terms = [self.adj(Y, z) for Y in self.ext if a in Y and Y != z] + [(a ^ 1) in z]
                self.B.exactly(terms, 1 if a in z else 2, tag=('full-a', a, z))
                n += 1
        for X in R + F:
            for z in F:
                if X == z:
                    continue
                pair = tuple(sorted((X, z)))
                key = min(pair, tuple(sorted((sig(X), sig(z)))))
                if key in seen:
                    continue
                seen.add(key)
                terms = [self.B.AND(self.adj(X, Y), self.adj(z, Y)) for Y in self.ext if Y != X and Y != z]
                self.B.exactly(terms + [self.adj(X, z)], 2 - len(set(X) & set(z)), tag=('full', X, z))
                n += 1
        return n

    @staticmethod
    def kind(X):
        c = X[0] // 2 == 0 or X[1] // 2 == 0
        d = X[0] // 2 == 1 or X[1] // 2 == 1
        return 'S' if c and d else 'C' if c else 'D' if d else 'F'

    def coupling(self):
        """one equivalence per sigma-class of pairs described by both systems; returns counts per kind of pair"""
        cnt = {}
        seen = set()
        for X, Y in itertools.combinations(self.ext, 2):
            key = min((X, Y), (sig(X), sig(Y)), (Y, X), (sig(Y), sig(X)))
            key = tuple(sorted(key))
            if key in seen:
                continue
            seen.add(key)
            dc, dd = self.descC(X, Y), self.descD(X, Y)
            if dc is None or dd is None:
                continue
            kk = ''.join(sorted(self.kind(X) + self.kind(Y)))
            cnt[kk] = cnt.get(kk, 0) + 1
            self.B.equiv(dc, dd, tag=('couple', X, Y))
        return cnt

    def cross(self):
        """two-block (iii): exact common-neighbour count of every non-shared c-root / non-shared d-root pair, mod sigma"""
        C = [X for X in self.ext if self.kind(X) == 'C']
        D = [X for X in self.ext if self.kind(X) == 'D']
        seen, n = set(), 0
        for r in C:
            for s in D:
                key = min((r, s), (sig(r), sig(s)))
                if key in seen:
                    continue
                seen.add(key)
                terms = [self.B.AND(self.descC(r, X), self.descD(s, X)) for X in self.ext if X != r and X != s]
                self.B.exactly(terms + [self.descC(r, s)], 2 - len(set(r) & set(s)), tag=('cross', r, s))
                n += 1
        return n


# ---------------------------------------------------------------- types, symmetry, instance list (m = 7)
L12 = list(range(12))


def b6():
    for s in itertools.permutations(range(6)):
        for f in itertools.product((0, 1), repeat=6):
            yield tuple(2 * s[l // 2] + ((l & 1) ^ f[l // 2]) for l in L12)


def conj(gp, m):
    r = [0] * len(gp)
    for l in range(len(gp)):
        r[gp[l]] = gp[m[l]]
    return r


def pivots(mu, tau):
    return {l for l in range(len(mu)) if mu[l] == tau[l]}


def load_types():
    return [json.load(open('t06_%d.json' % i)) for i in range(7)]


def make_list(p=0):
    """T = K7 reduction: take c = a block of MINIMAL type index, d = its c-frame block p. Then t_c <= t_d, the gluing of
    the 5 other blocks is a SAT variable, and gamma (the d-frame label of c+) runs over the labels allowed by pivot
    consistency, up to Stab(t_d)."""
    T = load_types()
    G = list(b6())
    stab = [[gp for gp in G if conj(gp, t['mu']) == t['mu'] and conj(gp, t['tau']) == t['tau']] for t in T]
    out = []
    for tc in range(7):
        pc = pivots(T[tc]['mu'], T[tc]['tau'])
        for td in range(tc, 7):
            pd = pivots(T[td]['mu'], T[td]['tau'])
            allowed = [gm for gm in L12 if ((2 * p) in pc) == ((gm ^ 1) in pd)]
            reps = sorted({min(gp[gm] for gp in stab[td]) for gm in allowed})
            for gm in reps:
                orbit = sorted({gp[gm] for gp in stab[td]})
                out.append({'id': len(out), 'tc': tc, 'td': td, 'p': p, 'gamma': gm, 'gamma_orbit': orbit,
                            'stab_c': len(stab[tc]), 'stab_d': len(stab[td]),
                            'mu_c': T[tc]['mu'], 'tau_c': T[tc]['tau'], 'mu_d': T[td]['mu'], 'tau_d': T[td]['tau'],
                            'src_c': T[tc]['note'], 'src_d': T[td]['note']})
    return out


COUPLEC = ('c COUPLING (structures, 27/09): world frame = N(x) with c = block 0, d = block 1; the c-system describes adj(X,Y) '
           'when X or Y is a c-root (type or U_c), the d-system when X or Y is a d-root (type or U_d through the gluing G). '
           'For every pair described by both, the two descriptions are equal: c-root x d-root (C1/C2 coupling), c-root x '
           'shared and shared x d-root (units from the other type), shared x far-far (U_c[Y,d+-] = U_d[Y,c+-]), shared x shared.')
CROSSC = ('c CROSS (structures, 27/09): two-block (iii). For a non-shared c-root r and non-shared d-root s every common '
          'neighbour is modelled (N(x) vertex, c-roots, d-roots, far-far vertices through U_c and U_d), so '
          '[r, s share an N(x) vertex] + sum_X adj(r,X) adj(s,X) + adj(r,s) = 2 (lambda = 1, mu = 2).')
GLUEC = ('c GLUE: G = signed bijection d-frame labels -> N(x) - d, G(gamma) = c+, G(pi gamma) = c-; B[beta,k,o]: d-frame '
         'block beta -> world block k with orientation o; FIXED if the comment below gives glue, else a SAT variable. '
         'Scope T = K7 (all pairs twisted, all blocks (0,6)): c = a block of minimal type, d = c-frame block p.')


FULLC = ('c FULL (structures, 27/09): variables for the adjacency of two far-far vertices (sigma-classes; {z, sigma z} never '
         'adjacent), and exact lambda/mu counts for every pair with a far-far member: (N(x) vertex outside c,d ; far-far), '
         '(root ; far-far), (far-far ; far-far). With the local systems and CROSS this is the complete srg(99,14,1,2) graph '
         'with sigma (every pair of vertices has an exact common-neighbour count).')


def build_instance(inst, flags, cross=True, glue=None, full=False):
    import sat_local2 as S2
    mc, tc_, md, td_ = inst['mu_c'], inst['tau_c'], inst['mu_d'], inst['tau_d']
    p, gm = inst['p'], inst['gamma']
    assert ((2 * p) in pivots(mc, tc_)) == ((gm ^ 1) in pivots(md, td_)), 'pivot consistency fails for this gamma'
    t0 = time.time()
    cnf_c, pool_c, far_c, _ = S2.encode(mc, tc_, **flags)
    cnf_d, pool_d, far_d, _ = S2.encode(md, td_, **flags)
    N1, N2 = pool_c.top, pool_d.top
    assert max(abs(l) for cl in cnf_c.clauses for l in cl) <= N1 and max(abs(l) for cl in cnf_d.clauses for l in cl) <= N2
    fc = {z: i for i, z in enumerate(far_c)}
    fd = {z: i for i, z in enumerate(far_d)}
    B = Builder(start=N1 + N2 + 1)
    Uc = lambda z, b: pool_c.obj2id[('u', fc[z], b)]
    Ud = lambda z, b: pool_d.obj2id[('u', fd[z], b)] + N1
    T = TwoBlock(B, 7, mc, tc_, p, md, td_, gm, Uc, Ud, glue=glue, full=full)
    cnt = T.coupling()
    ncross = T.cross() if cross else 0
    nfull = T.full_counts() if full else 0
    if B.violations:
        print('WARNING: %d constraints are statically violated (trivial UNSAT), e.g. %s' % (len(B.violations), B.violations[:3]))
    cnf = CNF()
    cnf.extend(cnf_c.clauses)
    cnf.extend([[l + N1 if l > 0 else l - N1 for l in cl] for cl in cnf_d.clauses])
    cnf.extend(B.cnf.clauses)
    info = {'vars': B.pool.top, 'clauses': len(cnf.clauses), 'c_system': [N1, len(cnf_c.clauses)],
            'd_system': [N2, len(cnf_d.clauses)], 'coupling_pairs': cnt, 'cross_equations': ncross,
            'full_equations': nfull, 'ff_vars': len(T.FF),
            'glue': 'fixed' if glue is not None else 'variable (B5)', 'static_violations': len(B.violations),
            'encode_s': round(time.time() - t0, 1)}
    build_instance.last = T
    return cnf, info


def glue_assumptions(T, glue):
    """B literals (assumptions) selecting a full gluing (list d-frame label -> world vertex) in a variable-gluing build"""
    out = []
    for b in range(T.m - 1):
        if b == T.gb:
            continue
        k, o = glue[2 * b] // 2, (glue[2 * b] & 1)
        out.append(T.Bv[b, k, o])
    return out


def random_glue(rng, gamma, m=7):
    wb = list(range(2, m)); rng.shuffle(wb)
    G = [None] * (2 * (m - 1))
    G[gamma], G[gamma ^ 1] = 0, 1
    for b, k in zip([b for b in range(m - 1) if b != gamma // 2], wb):
        o = rng.randrange(2)
        for e in (0, 1):
            G[2 * b + e] = 2 * k + (e ^ o)
    return G


# ---------------------------------------------------------------- self-tests
def world_from_graph(A, x0, neg, rng, m):
    """BvLS-type data: world vertex -> graph vertex (random order and orientation of the other blocks)"""
    n = len(A)
    N0 = [u for u in range(n) if A[x0][u]]
    blocks, seen = [], set()
    for u in N0:
        if u not in seen:
            blocks.append((u, neg[u])); seen |= {u, neg[u]}
    assert len(blocks) == m
    order = list(range(m)); rng.shuffle(order)
    phys = {}
    for k, bi in enumerate(order):
        u, v = blocks[bi]
        if rng.random() < 0.5:
            u, v = v, u
        phys[2 * k], phys[2 * k + 1] = u, v
    ext = {}
    for r in range(n):
        if r != x0 and not A[x0][r]:
            nb = [u for u in N0 if A[r][u]]
            assert len(nb) == 2
            inv = {g: w for w, g in phys.items()}
            ext[tuple(sorted((inv[nb[0]], inv[nb[1]])))] = r
    return phys, ext


def local_type(adjw, base, frame):
    """mu, tau of the system based at world vertex base; frame: label -> world vertex (N(x) minus base's block)"""
    n = len(frame)
    root = lambda b, l: tuple(sorted((b, frame[l])))
    mu = [next(k for k in range(n) if k != l and adjw(root(base, l), root(base, k))) for l in range(n)]
    h = [next(k for k in range(n) if adjw(root(base, l), root(base ^ 1, k))) for l in range(n)]
    return mu, [h[l] ^ 1 for l in range(n)]


def selftest_bvls(rng, reps=3):
    sys.path.insert(0, '../famille')
    from family_graphs import bvls_243
    V, adjl = bvls_243()
    V = [tuple(v) for v in V]
    index = {v: i for i, v in enumerate(V)}
    n = len(V)
    A = [[0] * n for _ in range(n)]
    for i, nb in enumerate(adjl):
        for j in nb:
            A[i][j] = 1
    negv = [index[tuple((-x) % 3 for x in v)] for v in V]
    x0 = index[tuple([0] * len(V[0]))]
    m = 11
    ok = True
    for rep in range(reps):
        phys, ext = world_from_graph(A, x0, negv, rng, m)
        adjw = lambda X, Y: A[ext[X]][ext[Y]] == 1
        p = rng.randrange(m - 1)
        others = [q for q in range(m - 1) if q != p]
        cblk = {p: 1}; cblk.update({q: k + 2 for k, q in enumerate(others)})
        Fc = [2 * cblk[l // 2] + (l & 1) for l in range(2 * (m - 1))]
        # random d-frame: G = signed bijection d-frame -> world vertices outside block 1, G(gamma) = c+
        wb = [0] + list(range(2, m)); rng.shuffle(wb)
        flips = [rng.randrange(2) for _ in range(m - 1)]
        G = [2 * wb[l // 2] + ((l & 1) ^ flips[l // 2]) for l in range(2 * (m - 1))]
        gamma = G.index(0)
        muc, tauc = local_type(adjw, 0, Fc)
        mud, taud = local_type(adjw, 2, G)
        Uc = lambda z, b: adjw(tuple(sorted((Fc[z[0]], Fc[z[1]]))), tuple(sorted((0, Fc[b]))))
        Ud = lambda z, b: adjw(tuple(sorted((G[z[0]], G[z[1]]))), tuple(sorted((2, G[b]))))
        Bd = Builder()
        T = TwoBlock(Bd, m, muc, tauc, p, mud, taud, gamma, Uc, Ud, glue=G, full=True, FFval=adjw)
        cnt = T.coupling(); nc = T.cross(); nf = T.full_counts()
        good = not Bd.violations and Bd.pool.top == 0
        print('BvLS rep %d: p=%d gamma=%d  coupling pairs %s  cross equations %d  full equations %d  violations %d  %s'
              % (rep, p, gamma, cnt, nc, nf, len(Bd.violations), 'OK' if good else 'FAIL %s' % Bd.violations[:3]), flush=True)
        ok &= good
        # negative controls: a wrong gluing (swap two other d-frame blocks) and one flipped U_c value
        b1, b2 = [b for b in range(m - 1) if b != gamma // 2][:2]
        G2 = list(G)
        for e in (0, 1):
            G2[2 * b1 + e], G2[2 * b2 + e] = G[2 * b2 + e], G[2 * b1 + e]
        Bn = Builder(); Tn = TwoBlock(Bn, m, muc, tauc, p, mud, taud, gamma, Uc, Ud, glue=G2); Tn.coupling(); Tn.cross()
        zf = (Fc.index(2), Fc.index(4)); bf = Fc.index(5)
        Ucf = lambda z, b: (not Uc(z, b)) if (z, b) == (tuple(sorted(zf)), bf) else Uc(z, b)
        Bf = Builder(); Tf = TwoBlock(Bf, m, muc, tauc, p, mud, taud, gamma, Ucf, Ud, glue=G); Tf.coupling(); Tf.cross()
        print('   negative controls: wrong gluing -> %d violations, one flipped U_c -> %d violations (both must be > 0)'
              % (len(Bn.violations), len(Bf.violations)), flush=True)
        ok &= bool(Bn.violations) and bool(Bf.violations)
    return ok


def selftest_synthetic(rng, reps=2):
    """m = 7, the real (0,6) types, a random consistent world (types exact, everything else random and sigma-invariant):
    (1) with the true fixed gluing every coupling equivalence holds; (2) with the gluing as SAT variables the coupling is
    SAT and its ONLY model is the true gluing; (3) a flipped U value is caught. CROSS is not tested here (random world)."""
    from pysat.solvers import Solver
    T6 = load_types()
    m, n = 7, 12
    ok = True
    for rep in range(reps):
        tc, td = rng.randrange(7), rng.randrange(7)
        muc, tauc, mud, taud = T6[tc]['mu'], T6[tc]['tau'], T6[td]['mu'], T6[td]['tau']
        p = rng.randrange(6)
        pc, pd = pivots(muc, tauc), pivots(mud, taud)
        gamma = rng.choice([gm for gm in L12 if ((2 * p) in pc) == ((gm ^ 1) in pd)])
        wb = list(range(2, m)); rng.shuffle(wb)
        dbl = [b for b in range(6) if b != gamma // 2]
        G = [None] * n
        G[gamma], G[gamma ^ 1] = 0, 1
        for b, k in zip(dbl, wb):
            o = rng.randrange(2)
            for e in (0, 1):
                G[2 * b + e] = 2 * k + (e ^ o)
        others = [q for q in range(6) if q != p]
        cblk = {p: 1}; cblk.update({q: k + 2 for k, q in enumerate(others)})
        Fc = [2 * cblk[l // 2] + (l & 1) for l in range(n)]
        ext = [z for z in itertools.combinations(range(2 * m), 2) if z[0] // 2 != z[1] // 2]
        adj = {}

        def put(X, Y, v):
            X, Y = tuple(sorted(X)), tuple(sorted(Y))
            for P in ((X, Y), (sig(X), sig(Y))):
                k = tuple(sorted(P))
                if k in adj:
                    assert adj[k] == v, ('conflict', k)
                adj[k] = v
        # the two types, by definition: (b+, l) ~ (b+, mu l), (b+, l) ~ (b-, h l); every other root pair of a block: 0
        for base, frame, mu, tau in ((0, Fc, muc, tauc), (2, G, mud, taud)):
            h = [tau[l] ^ 1 for l in range(n)]
            rts = [X for X in ext if base // 2 in (X[0] // 2, X[1] // 2)]
            want = set()
            for l in range(n):
                for (a, b) in (((base, frame[l]), (base, frame[mu[l]])), ((base, frame[l]), (base ^ 1, frame[h[l]]))):
                    X, Y = tuple(sorted(a)), tuple(sorted(b))
                    want.add(tuple(sorted((X, Y)))); want.add(tuple(sorted((sig(X), sig(Y)))))
            for X, Y in itertools.combinations(rts, 2):
                put(X, Y, tuple(sorted((X, Y))) in want)
        for X, Y in itertools.combinations(ext, 2):
            k = tuple(sorted((X, Y)))
            if k not in adj:
                put(X, Y, rng.random() < 0.3)
        adjw = lambda X, Y: adj[tuple(sorted((X, Y)))]
        Uc = lambda z, b: adjw(tuple(sorted((Fc[z[0]], Fc[z[1]]))), tuple(sorted((0, Fc[b]))))
        Ud = lambda z, b: adjw(tuple(sorted((G[z[0]], G[z[1]]))), tuple(sorted((2, G[b]))))
        B1 = Builder()
        T1 = TwoBlock(B1, m, muc, tauc, p, mud, taud, gamma, Uc, Ud, glue=G)
        cnt = T1.coupling()
        r1 = not B1.violations and B1.pool.top == 0
        B2 = Builder()
        T2 = TwoBlock(B2, m, muc, tauc, p, mud, taud, gamma, Uc, Ud, glue=None)
        T2.coupling()
        with Solver(name='glucose4', bootstrap_with=B2.cnf.clauses) as s:
            sat = s.solve()
            model = set(l for l in s.get_model() if l > 0) if sat else set()
            Gm = None
            if sat:
                Gm = [None] * n
                Gm[gamma], Gm[gamma ^ 1] = 0, 1
                for (b, k, o), v in T2.Bv.items():
                    if v in model:
                        for e in (0, 1):
                            Gm[2 * b + e] = 2 * k + (e ^ o)
                s.add_clause([-v for v in T2.Bv.values() if v in model])
                unique = not s.solve()
        r2 = sat and Gm == G and unique
        zf = (Fc.index(2), Fc.index(4)); bf = Fc.index(5)
        Ucf = lambda z, b: (not Uc(z, b)) if (z, b) == (tuple(sorted(zf)), bf) else Uc(z, b)
        B3 = Builder(); T3 = TwoBlock(B3, m, muc, tauc, p, mud, taud, gamma, Ucf, Ud, glue=G); T3.coupling()
        r3 = bool(B3.violations)
        print('synthetic m=7 rep %d: types (%d,%d) p=%d gamma=%d  pairs %s | fixed true glue: %s | variable glue: SAT %s, '
              'model = true glue %s, unique %s | flipped U caught: %s'
              % (rep, tc, td, p, gamma, cnt, 'OK' if r1 else 'FAIL', sat, Gm == G, unique if sat else '-', r3), flush=True)
        ok &= r1 and r2 and r3
    return ok


if __name__ == '__main__':
    if '--list' in sys.argv:
        lst = make_list()
        json.dump(lst, open('twoblock_instances.json', 'w'), indent=1)
        print(len(lst), 'instances written to twoblock_instances.json')
        sys.exit(0)
    if '--selftest' in sys.argv:
        rng = random.Random(20260927)
        a = selftest_bvls(rng)
        b = selftest_synthetic(rng)
        print('SELFTEST', 'PASSED' if a and b else 'FAILED')
        sys.exit(0 if a and b else 1)
    lst = json.load(open(sys.argv[1]))
    inst = lst[int(sys.argv[2])]
    flags = {'sec': '--sec' in sys.argv, 'twist': '--twist' in sys.argv}
    sweepfix = next((a.split('=')[1] for a in sys.argv if a.startswith('--sweepfix=')), None)
    if sweepfix:                                # --sweepfix=N,per_secs,seed: N random gluings, each a FRESH fixed build
        from pysat.solvers import Solver
        import threading
        nsw, per, seed = int(sweepfix.split(',')[0]), float(sweepfix.split(',')[1]), int(sweepfix.split(',')[2])
        rng = random.Random(seed)
        for i in range(nsw):
            G = random_glue(rng, inst['gamma'])
            cnf, info = build_instance(inst, flags, cross='--nocross' not in sys.argv, glue=G, full='--full' in sys.argv)
            with Solver(name='glucose4', bootstrap_with=cnf.clauses) as s:
                timer = threading.Timer(per, lambda: s.interrupt()); timer.start()
                t1 = time.time()
                res = s.solve_limited(expect_interrupt=True); timer.cancel()
            print('fixed glue %s : %s %.1fs' % (G, {True: 'SAT', False: 'UNSAT', None: 'UNKNOWN'}[res], time.time() - t1), flush=True)
        sys.exit(0)
    relax = next((tuple(a.split('=')[1].split(',')) for a in sys.argv if a.startswith('--relax=')), ())
    if relax:                                   # positive controls only: drop (ii) and/or (iii) at both blocks
        flags['relax'] = relax
    glue = next(([int(t) for t in a.split('=')[1].split(',')] for a in sys.argv if a.startswith('--glue=')), None)
    cnf, info = build_instance(inst, flags, cross='--nocross' not in sys.argv, glue=glue, full='--full' in sys.argv)
    tagname = 'twoblock_%d%s%s%s%s' % (inst['id'], '_glue' if glue else '', '_nocross' if '--nocross' in sys.argv else '',
                                       '_full' if '--full' in sys.argv else '',
                                     '_relax' + ''.join(relax) if relax else '')
    print(tagname, json.dumps(info), flush=True)
    if '--dimacs' in sys.argv:
        import sat_local2 as S2
        head = {k: inst[k] for k in ('id', 'tc', 'td', 'p', 'gamma', 'mu_c', 'tau_c', 'mu_d', 'tau_d')}
        com = ['c two-block instance ' + json.dumps(head), 'c flags: ' + json.dumps(flags) + (' glue=' + json.dumps(glue) if glue else ''),
               'c variables 1..%d: c-system (encoding 2, c-frame); next %d: d-system (encoding 2, d-frame); then glue/coupling/cross'
               % (info['c_system'][0], info['d_system'][0]), 'c local systems: sat_local2.encode (i)-(iii) + L1/S1/S2', GLUEC, COUPLEC]
        if '--nocross' not in sys.argv:
            com.append(CROSSC)
        if '--full' in sys.argv:
            com.append(FULLC)
        if flags['sec']:
            com.append(S2.SECC)
        if flags['twist']:
            com.append(S2.TWISTC)
        cnf.to_file(tagname + '.cnf', comments=com)
        print('DIMACS written to', tagname + '.cnf')
        sys.exit(0)
    from pysat.solvers import Solver
    import threading
    secs = float(next((a for a in sys.argv[3:] if a.replace('.', '').isdigit()), 240))
    sweep = next((a.split('=')[1] for a in sys.argv if a.startswith('--sweep=')), None)
    if sweep:                                   # --sweep=N,per_secs: N random gluings under assumptions, one solver
        nsw, per = int(sweep.split(',')[0]), float(sweep.split(',')[1])
        T = build_instance.last
        rng = random.Random(1)
        t0 = time.time()
        with Solver(name='glucose4', bootstrap_with=cnf.clauses) as s:
            for i in range(nsw):
                if time.time() - t0 + per > secs:
                    break
                G = random_glue(rng, inst['gamma'])
                timer = threading.Timer(per, lambda: s.interrupt()); timer.start()
                t1 = time.time()
                res = s.solve_limited(assumptions=glue_assumptions(T, G), expect_interrupt=True); timer.cancel()
                s.clear_interrupt()
                print('glue %s : %s %.1fs' % (G, {True: 'SAT', False: 'UNSAT', None: 'UNKNOWN'}[res], time.time() - t1), flush=True)
        sys.exit(0)
    t0 = time.time()
    with Solver(name='glucose4', bootstrap_with=cnf.clauses) as s:
        timer = threading.Timer(secs, lambda: s.interrupt()); timer.start()
        res = s.solve_limited(expect_interrupt=True); timer.cancel()
    print('status', {True: 'SAT', False: 'UNSAT', None: 'UNKNOWN (interrupted)'}[res], '%.1fs' % (time.time() - t0), flush=True)
