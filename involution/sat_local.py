"""SAT encoding (CaDiCaL via pysat) of the local system at a block c:
  (i) Phi bijection far roots -> non-mu pairs;  (ii) K2 counts;  (iii) K3 counts;  full block extras (iv), (v), no cancellation.
Variables x[z,e]; u[z,b] <-> OR_{e containing b} x[z,e]; w[Z,b,d] <-> u[z,b] & u[sigma z,d].
Cardinality constraints: pysat.card (totalizer / seqcounter).  UNSAT => no block of this type in an involution.
Usage: python sat_local.py cfg.json seconds [solver]"""
import itertools, json, sys, time, threading
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
import iota_design as D

L = list(range(12))
PI = [v ^ 1 for v in L]


def anchored_structures(mu, tau):
    """anchored (non-sib) left blocks at c: p+ and p- pivots (mu = tau on them) and mu p+ != pi p+.
    Returns (digons, paths): digons = list of (p, q) with mu mapping p's labels into block q (anchored digon, q anchored too);
    paths = anchored non-sib blocks whose two mates lie in two different blocks."""
    digons, paths = [], []
    for p in range(6):
        a, b = 2 * p, 2 * p + 1
        if not (mu[a] == tau[a] and mu[b] == tau[b]) or mu[a] == b:
            continue
        if mu[a] // 2 == mu[b] // 2:
            q = mu[a] // 2
            if p < q:
                digons.append((p, q))
        else:
            paths.append(p)
    return digons, paths


def unary_geq(cnf, pool, lits, K, tag):
    """o[j] <-> (sum of lits >= j), j = 1..K, by a sequential counter with full equivalences."""
    prev = {0: None}                                   # None = constant true
    for j in range(1, K + 1):
        prev[j] = False                                 # constant false for i = 0
    for i, x in enumerate(lits):
        cur = {0: None}
        for j in range(1, K + 1):
            v = pool.id((tag, i, j)); cur[j] = v
            a = prev[j]; b = prev[j - 1]
            # v <-> a or (b and x)
            if a is not False and a is not None:
                cnf.append([-a, v])
            if b is None:
                cnf.append([-x, v])
            elif b is not False:
                cnf.append([-b, -x, v])
            rhs = []
            if a is None:
                cnf.append([v]); continue
            if a is not False:
                rhs.append(a)
            if b is None:
                cnf.append([-v] + rhs + [x])
            elif b is False:
                cnf.append([-v] + rhs)
            else:
                cnf.append([-v] + rhs + [b]); cnf.append([-v] + rhs + [x])
        prev = cur
    return [prev[j] for j in range(1, K + 1)]


def path_constraints(cnf, pool, mu, tau, far, sig, U):
    """C1's 'exactly one excess' at every anchored non-sib PATH block p (notes/134 par. 6, NOTE_S29 par. 9):
    with k_Z = |Phi(z) u Phi(sigma z) & Sp|, Sp = {p+, p-, mu p+, mu p-}, S_p = sigma-pairs through block p:
    k <= 2 on S_p, k <= 3 elsewhere, and exactly one of {Z in S_p: k = 2}, {Z not in S_p: k = 3}, {Z not in S_p: k = 0}."""
    _, paths = anchored_structures(mu, tau)
    for p in paths:
        Sp = [2 * p, 2 * p + 1, mu[2 * p], mu[2 * p + 1]]
        events = []
        for z in range(60):
            if z > sig[z]:
                continue
            lits = [U(z, b) for b in Sp] + [U(sig[z], b) for b in Sp]
            o = unary_geq(cnf, pool, lits, 4, ('pathcnt', p, z))
            through = p in {a // 2 for a in far[z]}
            if through:
                cnf.append([-o[2]])                    # k <= 2 on S_p
                events.append(o[1])                    # k = 2  <->  o2 (given not o3)
            else:
                cnf.append([-o[3]])                    # k <= 3 (never 4)
                events.append(o[2])                    # k = 3  <->  o3
                e0 = pool.id(('pathe0', p, z)); cnf.append([-e0, -o[0]]); cnf.append([e0, o[0]])   # e0 <-> k = 0
                events.append(e0)
        cnf.extend(CardEnc.equals(lits=events, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses)
    return paths


PATHC = ('c PATH (C1 notes/134 par. 6, equality form of the single-block lemma): at every anchored non-sib path block p, '
         'k <= 2 on the sigma-pairs through p, never 4, and exactly one of {k = 2 through p}, {k = 3 not through p}, {k = 0 not through p}.')


SECC = ('c SEC (structures, 27/09): in a configuration WITHOUT dis D-squares (e.g. the Sigma = empty all-one-block family), '
        'a far root z with v in Phi(z) and v in Phi(sigma z) has D-partner o_v in star(c); the square is one-block, so z meets '
        'block(v). Clauses: not(U[z,v] and U[sigma z,v]) whenever block(v) is not a block of z. NOT valid if dis squares exist.')


def twist_terms(o, through):
    """literals (with repetition) whose number of true members is the contribution of one sigma-pair to C1's twisted-block
    identity, given o = [k>=1, k>=2, k>=3, k>=4]:  [k>=3] + 2[k>=4]  (+ [k>=2] + [k>=3] + [k>=4] if through p)
    (+ [k = 0] = not o[0] if not through p)."""
    t = [o[2], o[3], o[3]]
    if through:
        t += [o[1], o[2], o[3]]
    else:
        t += [-o[0]]
    return t


def twist_constraints(cnf, pool, mu, tau, far, sig, U):
    """C1's identity at every TWISTED left block p (exactly one pivot label p*, non-pivot p°, q* = mu p* in another block):
    with S = {p*, p°, q*} and k_Z = sum over S of U[z,b] + U[sigma z,b] (multiplicity: p° may be in both Phi's),
    a3 + 3 a4 + sum_{Z through p} (k-1)+ + #{Z not through p : k = 0} = 3."""
    blocks = []
    for p in range(6):
        a, b = 2 * p, 2 * p + 1
        piv = [l for l in (a, b) if mu[l] == tau[l]]
        if len(piv) != 1:
            continue
        ps = piv[0]; po = a + b - ps; qs = mu[ps]
        if qs // 2 == p:
            continue                                   # sib (not a twisted block)
        S = [ps, po, qs]
        lits = []
        for z in range(60):
            if z > sig[z]:
                continue
            o = unary_geq(cnf, pool, [U(z, x) for x in S] + [U(sig[z], x) for x in S], 4, ('twcnt', p, z))
            cnf.append([-o[3], o[3]])                   # (k <= 4 automatically: only 4 slots)
            through = p in {w // 2 for w in far[z]}
            for t in twist_terms(o, through):
                c = pool.id(('twcopy', p, z, len(lits)))  # fresh copy, so that repeated literals are separate inputs
                cnf.append([-c, t]); cnf.append([c, -t])
                lits.append(c)
        cnf.extend(CardEnc.equals(lits=lits, bound=3, vpool=pool, encoding=EncType.seqcounter).clauses)
        blocks.append(p)
    return blocks


TWISTC = ('c TWIST (C1, 27/09; relu par structures): at every twisted left block p (one pivot p*, non-pivot p°, q* = mu p*), '
          'with k counted with multiplicity on S = {p*, p°, q*}: a3 + 3 a4 + sum_{Z through p} (k-1)+ + #{Z not through p: k = 0} = 3. '
          'Uses only (i)-(iii) and block-level coverage (no dis-square hypothesis).')


def encode(mu, tau, cycles=(), digon=False, path=False, sec=False, twist=False, ident=None, beta=None):
    K2, K3 = D.K2K3(mu, tau)
    far = [z for z in itertools.combinations(L, 2) if z[0] // 2 != z[1] // 2]
    fidx = {z: i for i, z in enumerate(far)}
    sig = [fidx[tuple(sorted((PI[z[0]], PI[z[1]])))] for z in far]
    edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    pool = IDPool()
    X = lambda z, e: pool.id(('x', z, e))
    U = lambda z, b: pool.id(('u', z, b))
    Wv = lambda z, b, d: pool.id(('w', z, b, d))
    cnf = CNF()

    def eq(lits, k):
        if k < 0 or k > len(lits):
            cnf.append([])            # impossible
            return
        enc = CardEnc.equals(lits=lits, bound=k, vpool=pool, encoding=EncType.seqcounter)
        cnf.extend(enc.clauses)
    for z in range(60):
        eq([X(z, e) for e in range(60)], 1)
    for e in range(60):
        eq([X(z, e) for z in range(60)], 1)
    for z in range(60):
        for b in L:
            es = [e for e in range(60) if b in edges[e]]
            cnf.append([-U(z, b)] + [X(z, e) for e in es])
            for e in es:
                cnf.append([-X(z, e), U(z, b)])
    for a in L:                                               # (ii)
        for b in L:
            eq([U(z, b) for z in range(60) if a in far[z]], int(K2[a, b]))
    reps = [z for z in range(60) if z < sig[z]]
    for z in reps:                                            # w = u(z,b) & u(sz,d)
        for b in L:
            for d in L:
                w = Wv(z, b, d)
                cnf.append([-w, U(z, b)]); cnf.append([-w, U(sig[z], d)])
                cnf.append([w, -U(z, b), -U(sig[z], d)])
    for b in L:                                               # (iii)
        for d in L:
            eq([Wv(z, b, d) for z in reps] + [Wv(z, d, b) for z in reps], int(K3[b, d])) if b != d else \
                eq([Wv(z, b, b) for z in reps], int(K3[b, b]) // 2)
    # pointwise conditions (C1, 27/09; re-derived from lambda = 1, mu = 2; valid at every block type), for z and for sigma z:
    #  (L1) a in blocks(z): not (a in Phi z and pi a in Phi sigma z)
    #  (S1) l in z: l in Phi z  => eta l not in Phi sigma z          (eta = tau in C2 notation = pi o h)
    #  (S2) l in z: l not in Phi z => not (mu l in Phi z and eta l in Phi sigma z)
    for z in range(60):
        sz = sig[z]
        for a in far[z]:
            for aa in (a, PI[a]):
                cnf.append([-U(z, aa), -U(sz, PI[aa])])
            l = a
            cnf.append([-U(z, l), -U(sz, tau[l])])
            cnf.append([U(z, l), -U(z, mu[l]), -U(sz, tau[l])])
    if digon:                                                 # C1's digon balance (anchored non-sib digons, any type)
        dg, _ = anchored_structures(mu, tau)
        for (p, q) in dg:
            Sp = [2 * p, 2 * p + 1, 2 * q, 2 * q + 1]
            for z in range(60):
                if z > sig[z]:
                    continue
                eq([U(z, b) for b in Sp] + [U(sig[z], b) for b in Sp], 2 - len({a // 2 for a in far[z]} & {p, q}))
        print('digon balance added for anchored digons', dg, flush=True)
    if sec:                                                   # (SEC): no dis D-square between far(c) and star(c)
        nsec = 0
        for z in range(60):
            zb = {a // 2 for a in far[z]}
            for v in L:
                if v // 2 not in zb:
                    cnf.append([-U(z, v), -U(sig[z], v)]); nsec += 1
        print('SEC clauses added:', nsec, flush=True)
    if twist:
        tb = twist_constraints(cnf, pool, mu, tau, far, sig, U)
        print('twist identities added for twisted blocks', tb, flush=True)
    if beta is not None:                                      # C1's finish line (shared code in sat_local2)
        from sat_local2 import beta_constraints
        npv, nnp, need = beta_constraints(cnf, pool, mu, tau, far, sig, U, beta)
        print('beta (%s): %d pivots, %d non-pivot labels, need sum beta >= %d' % (beta, npv, nnp, need), flush=True)
    if ident is not None:                                     # C1's general identity (shared code in sat_local2)
        from sat_local2 import ident_list, ident_constraints
        ids = ident_list(mu, tau, ident)
        ident_constraints(cnf, pool, far, sig, U, ids)
        print('general identities added (slack <= %d): %d' % (ident, len(ids)), flush=True)
    if path:
        pb = path_constraints(cnf, pool, mu, tau, far, sig, U)
        print('path constraints (exactly one excess) added for anchored path blocks', pb, flush=True)
    if cycles:                                                # full block extras
        for z in range(60):
            sz = sig[z]
            for b in L:
                cnf.append([-U(z, b), -U(sz, b)])             # no cancellation
            for a in far[z]:
                cnf.append([-U(z, a), -U(sz, PI[a])])        # (iv)
                cnf.append([-U(z, PI[a]), -U(sz, a)])
            for C in cycles:                                  # (v)
                lc = [2 * p + s for p in C for s in (0, 1)]
                need = len(C) - len({a // 2 for a in far[z]} & set(C))
                eq([U(z, v) for v in lc] + [U(sz, v) for v in lc], need)
    return cnf, pool, far, edges


if __name__ == '__main__':
    cfg = json.load(open(sys.argv[1]))
    secs = float(sys.argv[2]) if len(sys.argv) > 2 else 120
    sname = sys.argv[3] if len(sys.argv) > 3 else 'cadical195'
    if 'cycles' in cfg:
        from fullblock_rigid import mates_from_D
        cycles = [tuple(c) for c in cfg['cycles']]
        m = mates_from_D(cycles); mu = tau = m
    else:
        cycles = []; mu, tau = cfg['mu'], cfg['tau']
    t0 = time.time()
    ident = next((int(a.split('=')[1]) if '=' in a else 0 for a in sys.argv if a.startswith('--ident')), None)
    beta = next(((a.split('=')[1] if '=' in a else 'card') for a in sys.argv if a.startswith('--beta')), None)
    cnf, pool, far, edges = encode(mu, tau, cycles, digon='--digon' in sys.argv, path='--path' in sys.argv, sec='--sec' in sys.argv, twist='--twist' in sys.argv, ident=ident, beta=beta)
    if '--symbreak' in sys.argv and cycles:        # C1's lex-leader symmetry breaking (conway99-c1 code/fullblock_symmetry.py)
        sys.path.insert(0, 'C:/home/researcher/Desktop/c1_recherche/code')
        from fullblock_symmetry import lexleader_clauses
        pref = int(next((a.split('=')[1] for a in sys.argv if a.startswith('--prefix=')), 10))
        cnt = itertools.count()
        lex = lexleader_clauses(mu, lambda z, e: pool.id(('x', z, e)), lambda: pool.id(('lex', next(cnt))), prefix_roots=pref)
        cnf.extend(lex)
        print('lex-leader clauses added:', len(lex), '(prefix %d)' % pref, flush=True)
    print(sys.argv[1], 'vars', pool.top, 'clauses', len(cnf.clauses), 'encode %.1fs' % (time.time() - t0), flush=True)
    if '--dimacs' in sys.argv:                 # for an external solver with a DRAT proof: cadical in.cnf proof.drat; drat-trim
        out = sys.argv[1].replace('.json', '.cnf')
        com = ['c local system (i)-(iii) + full-block extras + C1 pointwise (L1,S1,S2)', 'c cfg ' + json.dumps(cfg),
               'c flags: ' + ' '.join(a for a in sys.argv[3:] if a.startswith('--'))]
        if beta is not None:
            from sat_local2 import BETAC
            com += [BETAC, 'c BETA mode: %s' % beta]
        if ident is not None:
            from sat_local2 import IDENTC, ident_list
            com += [IDENTC, 'c IDENT list (p, S, slack): ' + json.dumps(ident_list(mu, tau, ident))]
        cnf.to_file(out, comments=com)
        print('DIMACS written to', out, flush=True)
        sys.exit(0)
    with Solver(name=sname, bootstrap_with=cnf.clauses) as s:
        timer = threading.Timer(secs, lambda: s.interrupt())
        timer.start()
        res = s.solve_limited(expect_interrupt=True)
        timer.cancel()
        st = {True: 'SAT', False: 'UNSAT', None: 'UNKNOWN (interrupted)'}[res]
        print('status', st, '%.1fs' % (time.time() - t0), flush=True)
        if res:
            model = set(l for l in s.get_model() if l > 0)
            phi = {str(list(far[z])): list(edges[e]) for z in range(60) for e in range(60) if pool.id(('x', z, e)) in model}
            json.dump({'cfg': cfg, 'Phi': phi}, open(sys.argv[1].replace('.json', '_satsol.json'), 'w'), indent=1)
