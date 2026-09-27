"""Second SAT encoding of the local system at a block (any type (mu, tau), C2 notation), built around the iota-pairing:
  x[z,e]  : Phi(z) = e (bijection far roots -> non-mu pairs)                        (i)
  U[z,b]  : b in Phi(z); (ii) as cardinality over U (10 literals each)               (ii)
  q[{e,f}]: iota pairs e with f (Phi sigma z = f when Phi z = e); perfect matching on the non-mu pairs;
            (iii) as cardinality over q (cross counts; diagonal = secant pairs)      (iii)
  linking : q[{e,f}] & x[z,e] -> x[sigma z, f]   and   x[z,e] & x[sigma z, f] -> q[{e,f}]   (ternary clauses)
  pointwise (C1): L1, S1, S2 on (x, U); full block extras (no cancellation, (iv), (v)) if cfg has "cycles".
Same logical content as sat_local.py (UNSAT here <=> UNSAT there); only the encoding of (iii) differs.
Usage: python sat_local2.py cfg.json seconds [--dimacs] [--solver=name]"""
import itertools, json, sys, time, threading
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
import iota_design as D

L = list(range(12)); PI = [v ^ 1 for v in L]


LEMMA = ('c LEMMA (structures NOTE_S29 par. 9, generalising C1 notes/134): at a block c of any type (mu, tau), let p be a '
         'left block anchored at c and NOT sib (p+ and p- pivots: mu = tau on them, and mu p+ != p-), Sp = {p+, p-, mu p+, mu p-}; then among the 30 '
         'iota-quadruples none contains the 4 labels of Sp and at most one contains 3 of them (none if mu p+, mu p- lie in '
         'one block). Proof: block-level (ii) gives coverage 2 on Sp over the 10 quadruples through p, so at least 2 of '
         'them avoid Sp, while inclusion-exclusion over the near-design gives a0 = 3 - a3 - 3 a4 (2 - a3 - 3 a4 for a digon).')


def lemma_cuts(mu, tau, edges, Q):
    """returns (unit clauses, [(literals, bound)]) for the generalised single-block lemma (anchored non-sib left blocks)."""
    units, amos = [], []
    for p in range(6):
        a, b = 2 * p, 2 * p + 1
        if not (mu[a] == tau[a] and mu[b] == tau[b]) or mu[a] == b:
            continue                                   # not anchored, or sib (the lemma does not apply)
        Sp = {a, b, mu[a], mu[b]}
        digon = mu[a] // 2 == mu[b] // 2
        three = []
        for (i, j), v in Q.items():
            k = len((set(edges[i]) | set(edges[j])) & Sp)
            if k == 4:
                units.append([-v])
            elif k == 3:
                three.append(v)
        if digon:
            units.extend([[-v] for v in three])
        elif len(three) > 1:
            amos.append((three, 1))
    return units, amos


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


DIGON = ('c DIGON (C1, 27/09, derived from the equality form of the single-block lemma): for an anchored non-sib digon {p, q} '
         'at c (p, q anchored, mu maps the labels of p onto those of q), every sigma-pair Z = {z, sigma z} has exactly '
         '2 - |blocks(z) & {p, q}| labels of Sp = {p+, p-, q+, q-} in Phi(z) u Phi(sigma z) (the D-cycle balance, in any type).')


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


def encode(mu, tau, cycles=(), cuts=False, digon=False, path=False, sec=False, twist=False):
    K2, K3 = D.K2K3(mu, tau)
    far = [z for z in itertools.combinations(L, 2) if z[0] // 2 != z[1] // 2]
    fidx = {z: i for i, z in enumerate(far)}
    sig = [fidx[tuple(sorted((PI[z[0]], PI[z[1]])))] for z in far]
    edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    E = len(edges)
    pool = IDPool(); cnf = CNF()
    X = lambda z, e: pool.id(('x', z, e))
    U = lambda z, b: pool.id(('u', z, b))

    def eq(lits, k):
        if k < 0 or k > len(lits):
            cnf.append([]); return
        cnf.extend(CardEnc.equals(lits=lits, bound=k, vpool=pool, encoding=EncType.seqcounter).clauses)
    for z in range(60):
        eq([X(z, e) for e in range(E)], 1)
    for e in range(E):
        eq([X(z, e) for z in range(60)], 1)
    for z in range(60):
        for b in L:
            es = [e for e in range(E) if b in edges[e]]
            cnf.append([-U(z, b)] + [X(z, e) for e in es])
            for e in es:
                cnf.append([-X(z, e), U(z, b)])
    for a in L:                                                        # (ii)
        for b in L:
            eq([U(z, b) for z in range(60) if a in far[z]], int(K2[a, b]))
    # iota-pair variables: allowed iff disjoint, or secant at a label b with K3[b,b] = 2
    Q = {}
    for i in range(E):
        for j in range(i + 1, E):
            com = set(edges[i]) & set(edges[j])
            if any(K3[b, b] < 2 for b in com):
                continue
            Q[i, j] = pool.id(('q', i, j))
    qv = lambda i, j: Q.get((min(i, j), max(i, j)))
    for i in range(E):
        eq([qv(i, j) for j in range(E) if j != i and qv(i, j)], 1)
    acc = {}
    for (i, j), v in Q.items():
        for b in edges[i]:
            for d in edges[j]:
                acc.setdefault((min(b, d), max(b, d)), []).append(v)
    for b in L:                                                        # (iii) through iota
        for d in L:
            if d < b:
                continue
            eq(acc.get((b, d), []), int(K3[b, d]) if b != d else int(K3[b, b]) // 2)
    for z in range(60):                                                # linking x <-> q
        sz = sig[z]
        for e in range(E):
            for f in range(E):
                if f == e:
                    continue
                v = qv(e, f)
                if v is None:
                    cnf.append([-X(z, e), -X(sz, f)])                  # forbidden pairing
                else:
                    cnf.append([-v, -X(z, e), X(sz, f)])
                    cnf.append([-X(z, e), -X(sz, f), v])
    for z in range(60):                                                # pointwise (C1)
        sz = sig[z]
        for a in far[z]:
            for aa in (a, PI[a]):
                cnf.append([-U(z, aa), -U(sz, PI[aa])])
            cnf.append([-U(z, a), -U(sz, tau[a])])
            cnf.append([U(z, a), -U(z, mu[a]), -U(sz, tau[a])])
    if cuts:
        units, amos = lemma_cuts(mu, tau, edges, Q)
        cnf.extend(units)
        for lits, bnd in amos:
            cnf.extend(CardEnc.atmost(lits=lits, bound=bnd, vpool=pool, encoding=EncType.seqcounter).clauses)
        print('lemma cuts added: %d unit clauses, %d at-most-one groups' % (len(units), len(amos)), flush=True)
    if digon:
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
    if path:
        pb = path_constraints(cnf, pool, mu, tau, far, sig, U)
        print('path constraints (exactly one excess) added for anchored path blocks', pb, flush=True)
    if cycles:
        for z in range(60):
            sz = sig[z]
            for b in L:
                cnf.append([-U(z, b), -U(sz, b)])
            for C in cycles:
                lc = [2 * p + s for p in C for s in (0, 1)]
                need = len(C) - len({a // 2 for a in far[z]} & set(C))
                eq([U(z, v) for v in lc] + [U(sz, v) for v in lc], need)
    return cnf, pool, far, edges


if __name__ == '__main__':
    cfg = json.load(open(sys.argv[1]))
    secs = float(sys.argv[2]) if len(sys.argv) > 2 else 120
    sname = next((a.split('=')[1] for a in sys.argv if a.startswith('--solver=')), 'cadical195')
    if 'cycles' in cfg:
        from fullblock_rigid import mates_from_D
        cycles = [tuple(c) for c in cfg['cycles']]
        m = mates_from_D(cycles); mu = tau = m
    else:
        cycles = []; mu, tau = cfg['mu'], cfg['tau']
    t0 = time.time()
    cuts = '--cuts' in sys.argv; digon = '--digon' in sys.argv; path = '--path' in sys.argv; sec = '--sec' in sys.argv; twist = '--twist' in sys.argv
    cnf, pool, far, edges = encode(mu, tau, cycles, cuts=cuts, digon=digon, path=path, sec=sec, twist=twist)
    print(sys.argv[1], '[encoding 2: iota-pairs]' + (' + lemma cuts' if cuts else '') + (' + digon balance' if digon else '') + (' + path' if path else '') + (' + SEC' if sec else '') + (' + twist' if twist else ''), 'vars', pool.top, 'clauses', len(cnf.clauses), 'encode %.1fs' % (time.time() - t0), flush=True)
    if '--dimacs' in sys.argv:
        out = sys.argv[1].replace('.json', '_enc2.cnf')
        com = ['c local system, encoding 2 (iota-pairs), ' + json.dumps(cfg)]
        if cuts:
            com.append(LEMMA)
        if digon:
            com.append(DIGON)
        if path:
            com.append(PATHC)
        if sec:
            com.append(SECC)
        if twist:
            com.append(TWISTC)
        cnf.to_file(out, comments=com)
        print('DIMACS written to', out); sys.exit(0)
    with Solver(name=sname, bootstrap_with=cnf.clauses) as s:
        timer = threading.Timer(secs, lambda: s.interrupt()); timer.start()
        res = s.solve_limited(expect_interrupt=True); timer.cancel()
        print('status', {True: 'SAT', False: 'UNSAT', None: 'UNKNOWN (interrupted)'}[res], '%.1fs' % (time.time() - t0), flush=True)
        if res:
            model = set(l for l in s.get_model() if l > 0)
            phi = {str(list(far[z])): list(edges[e]) for z in range(60) for e in range(len(edges)) if pool.id(('x', z, e)) in model}
            json.dump({'cfg': cfg, 'Phi': phi}, open(sys.argv[1].replace('.json', '_satsol2.json'), 'w'), indent=1)
