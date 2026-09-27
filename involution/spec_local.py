"""Declarable symmetry breaking for the local systems (structures 27/09, for C2's (A) lanes; C1's lex-leader).

Chain of the certificate:  (CNF = compile(SPEC) + LEX, DRAT VERIFIED)  +  (SPEC invariant under G, VERIFIED here)  +
(LEX = C1's lexleader_clauses over permutations of G = C1's pivot_types.stabiliser(mu, tau), re-generated here).
If every g in G maps the constraint system to itself, every G-orbit of models contains its lex-minimum, which satisfies LEX;
so CNF + LEX UNSAT implies CNF UNSAT (the auxiliary variables are only existentially quantified).

SPEC = the ordered list of events that sat_local2.encode performs, at the level of PRIMARY variables (x, u, q, beta):
  ('alloc', key) | ('CL', lits) | ('EQ'|'LE'|'GE', lits, k) | ('TWIST', p) | ('PATH', p) | ('IDENT', nid, p, S, slack)
with lits = [(sign, key)].  compile_spec replays it with a fresh pool: it must give EXACTLY sat_local2.encode's CNF
(tested byte for byte), so the spec is complete.  The invariance test works on label-level canonical forms:
x(z, e), u(z, b), q({e, f}), beta(v); IDENT/TWIST as (block, label set, slack), PATH as (block, Sp).
Usage:
  python spec_local.py export types.json i OUT.cnf [--cuts --digon --path --twist --ident=1 --balance --beta[=m] --sec]
                                                  [--lex=all|gens|sample:N:seed] [--prefix=6]
  python spec_local.py verify OUT.cnf types.json i
  python spec_local.py selftest"""
import itertools, json, sys, hashlib, random, os, io, contextlib
from collections import Counter
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType
import iota_design as D
import sat_local2 as S2

L = list(range(12)); PI = [v ^ 1 for v in L]
C1CODE = os.environ.get('C1_CODE', 'C:/home/researcher/Desktop/c1_recherche/code')


def c1_modules():
    sys.path.insert(0, C1CODE)
    import pivot_types as PT, fullblock_symmetry as FS
    shas = {f: hashlib.sha256(open(os.path.join(C1CODE, f), 'rb').read()).hexdigest()
            for f in ('pivot_types.py', 'fullblock_symmetry.py')}
    return PT, FS, shas


def frame(mu):
    far = [z for z in itertools.combinations(L, 2) if z[0] // 2 != z[1] // 2]
    fidx = {z: i for i, z in enumerate(far)}
    sig = [fidx[tuple(sorted((PI[z[0]], PI[z[1]])))] for z in far]
    edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    return far, fidx, sig, edges


# ------------------------------------------------------------------ spec (mirror of sat_local2.encode)
def build_spec(mu, tau, cuts=False, digon=False, path=False, sec=False, twist=False, relax=(), ident=None, fullcyc=False,
               balance=False, beta=None):
    K2, K3 = D.K2K3(mu, tau)
    far, fidx, sig, edges = frame(mu)
    E = len(edges)
    ev, seen = [], set()

    def key(k):
        if k not in seen:
            seen.add(k); ev.append(('alloc', k))
        return k
    X = lambda z, e: key(('x', z, e))
    U = lambda z, b: key(('u', z, b))
    eq = lambda lits, k: ev.append(('EQ', [(1, l) for l in lits], k))
    cl = lambda lits: ev.append(('CL', list(lits)))
    for z in range(60):
        eq([X(z, e) for e in range(E)], 1)
    for e in range(E):
        eq([X(z, e) for z in range(60)], 1)
    for z in range(60):
        for b in L:
            es = [e for e in range(E) if b in edges[e]]
            cl([(-1, U(z, b))] + [(1, X(z, e)) for e in es])
            for e in es:
                cl([(-1, X(z, e)), (1, U(z, b))])
    for a in L:
        for b in L:
            if 'ii' not in relax:
                eq([U(z, b) for z in range(60) if a in far[z]], int(K2[a, b]))
    Q = {}
    for i in range(E):
        for j in range(i + 1, E):
            com = set(edges[i]) & set(edges[j])
            if any(K3[b, b] < 2 for b in com):
                continue
            Q[i, j] = key(('q', i, j))
    qv = lambda i, j: Q.get((min(i, j), max(i, j)))
    for i in range(E):
        eq([qv(i, j) for j in range(E) if j != i and qv(i, j)], 1)
    acc = {}
    for (i, j), v in Q.items():
        for b in edges[i]:
            for d in edges[j]:
                acc.setdefault((min(b, d), max(b, d)), []).append(v)
    for b in L:
        for d in L:
            if d < b or 'iii' in relax:
                continue
            eq(acc.get((b, d), []), int(K3[b, d]) if b != d else int(K3[b, b]) // 2)
    for z in range(60):
        sz = sig[z]
        for e in range(E):
            for f in range(E):
                if f == e:
                    continue
                v = qv(e, f)
                if v is None:
                    cl([(-1, X(z, e)), (-1, X(sz, f))])
                else:
                    cl([(-1, v), (-1, X(z, e)), (1, X(sz, f))])
                    cl([(-1, X(z, e)), (-1, X(sz, f)), (1, v)])
    for z in range(60):
        sz = sig[z]
        for a in far[z]:
            for aa in (a, PI[a]):
                cl([(-1, U(z, aa)), (-1, U(sz, PI[aa]))])
            cl([(-1, U(z, a)), (-1, U(sz, tau[a]))])
            cl([(1, U(z, a)), (-1, U(z, mu[a])), (-1, U(sz, tau[a]))])
    if cuts:                                           # sat_local2.lemma_cuts, same order
        units, amos = [], []
        for p in range(6):
            a, b = 2 * p, 2 * p + 1
            if not (mu[a] == tau[a] and mu[b] == tau[b]) or mu[a] == b:
                continue
            Sp = {a, b, mu[a], mu[b]}
            dgn = mu[a] // 2 == mu[b] // 2
            three = []
            for (i, j), v in Q.items():
                k = len((set(edges[i]) | set(edges[j])) & Sp)
                if k == 4:
                    units.append([(-1, v)])
                elif k == 3:
                    three.append(v)
            if dgn:
                units.extend([[(-1, v)] for v in three])
            elif len(three) > 1:
                amos.append((three, 1))
        for u in units:
            cl(u)
        for lits, bnd in amos:
            ev.append(('LE', [(1, l) for l in lits], bnd))
    if digon:
        dg, _ = S2.anchored_structures(mu, tau)
        for (p, q) in dg:
            Sp = [2 * p, 2 * p + 1, 2 * q, 2 * q + 1]
            for z in range(60):
                if z > sig[z]:
                    continue
                eq([U(z, b) for b in Sp] + [U(sig[z], b) for b in Sp], 2 - len({a // 2 for a in far[z]} & {p, q}))
    if sec:
        for z in range(60):
            zb = {a // 2 for a in far[z]}
            for v in L:
                if v // 2 not in zb:
                    cl([(-1, U(z, v)), (-1, U(sig[z], v))])
    if twist:
        for p in range(6):
            if S2.twist_S3(mu, tau, p) is not None:
                ev.append(('TWIST', p))
    if path:
        _, paths = S2.anchored_structures(mu, tau)
        for p in paths:
            ev.append(('PATH', p))
    if fullcyc and all(mu[l] == tau[l] for l in L):
        for C in [C for C in S2.block_cycles(mu) if 3 <= len(C) <= 4]:
            lc = [2 * p + s for p in C for s in (0, 1)]
            for z in range(60):
                if z > sig[z]:
                    continue
                eq([U(z, v) for v in lc] + [U(sig[z], v) for v in lc], len(C) - len({a // 2 for a in far[z]} & set(C)))
    if beta is not None:                              # sat_local2.beta_constraints, same order
        nonpiv = [v for v in L if mu[v] != tau[v]]
        npiv = 12 - len(nonpiv)
        B = {v: key(('beta', v)) for v in nonpiv}
        for v in nonpiv:
            for z in range(60):
                meets = v // 2 in {a // 2 for a in far[z]}
                cl([((-1 if meets else 1), B[v]), (-1, U(z, v)), (-1, U(sig[z], v))])
        for v in nonpiv:
            if mu[v] == PI[v] or tau[v] == PI[v]:
                cl([(1, B[v])])
                if PI[v] in B:
                    cl([(1, B[PI[v]])])
        need = max(0, 12 - 2 * npiv)
        if beta == 'all':
            for v in nonpiv:
                cl([(1, B[v])])
        elif beta == 'zero':
            for v in nonpiv:
                cl([(-1, B[v])])
        elif need > 0:
            if need > len(nonpiv):
                cl([])
            else:
                ev.append(('GE', [(1, B[v]) for v in nonpiv], need))
    if balance:
        obs = S2.pmt_orbits(mu, tau)
        if len(obs) > 1:
            for S in obs:
                for z in range(60):
                    if z > sig[z]:
                        continue
                    eq([U(z, b) for b in S] + [U(sig[z], b) for b in S], len(S) // 2 - len(set(far[z]) & set(S)))
    if ident is not None:
        for nid, (p, S, sl) in enumerate(S2.ident_list(mu, tau, ident)):
            ev.append(('IDENT', nid, p, tuple(S), sl))
    return ev


def compile_spec(ev, mu, tau):
    far, fidx, sig, edges = frame(mu)
    pool = IDPool(); cnf = CNF()
    U = lambda z, b: pool.id(('u', z, b))
    lit = lambda s, k: s * pool.obj2id[k]
    for e in ev:
        t = e[0]
        if t == 'alloc':
            assert e[1] not in pool.obj2id, e
            pool.id(e[1])
        elif t == 'CL':
            cnf.append([lit(s, k) for s, k in e[1]])
        elif t in ('EQ', 'LE', 'GE'):
            lits = [lit(s, k) for s, k in e[1]]; k = e[2]
            if t == 'EQ':
                if k < 0 or k > len(lits):
                    cnf.append([])
                else:
                    cnf.extend(CardEnc.equals(lits=lits, bound=k, vpool=pool, encoding=EncType.seqcounter).clauses)
            elif t == 'LE':
                cnf.extend(CardEnc.atmost(lits=lits, bound=k, vpool=pool, encoding=EncType.seqcounter).clauses)
            else:
                cnf.extend(CardEnc.atleast(lits=lits, bound=k, vpool=pool, encoding=EncType.seqcounter).clauses)
        elif t == 'TWIST':
            assert S2.twist_block(cnf, pool, mu, tau, far, sig, U, e[1])
        elif t == 'PATH':
            S2.path_block(cnf, pool, mu, far, sig, U, e[1])
        elif t == 'IDENT':
            S2.ident_one(cnf, pool, far, sig, U, e[1], e[2], list(e[3]), e[4])
        else:
            raise ValueError(e)
    return cnf, pool


# ------------------------------------------------------------------ canonical forms and the action of B6
def canon_spec(ev, mu, tau):
    far, fidx, sig, edges = frame(mu)

    def name(k):
        if k[0] == 'x':
            return ('x', far[k[1]], edges[k[2]])
        if k[0] == 'u':
            return ('u', far[k[1]], k[2])
        if k[0] == 'q':
            return ('q', tuple(sorted((edges[k[1]], edges[k[2]]))))
        if k[0] == 'beta':
            return ('beta', k[1])
        raise ValueError(k)
    out = []
    for e in ev:
        t = e[0]
        if t == 'alloc':
            continue
        if t == 'CL':
            out.append(('CL', tuple(sorted({(s, name(k)) for s, k in e[1]}))))
        elif t in ('EQ', 'LE', 'GE'):
            out.append((t, tuple(sorted((s, name(k)) for s, k in e[1])), e[2]))
        elif t == 'TWIST':
            out.append(('IDENT', e[1], tuple(sorted(S2.twist_S3(mu, tau, e[1]))), 3))
        elif t == 'PATH':
            p = e[1]
            out.append(('PATH', p, tuple(sorted([2 * p, 2 * p + 1, mu[2 * p], mu[2 * p + 1]]))))
        elif t == 'IDENT':
            out.append(('IDENT', e[2], tuple(sorted(e[3])), e[4]))
    return Counter(out)


def act(g, c):
    gp = lambda e: tuple(sorted((g[e[0]], g[e[1]])))

    def gn(n):
        if n[0] == 'x':
            return ('x', gp(n[1]), gp(n[2]))
        if n[0] == 'u':
            return ('u', gp(n[1]), g[n[2]])
        if n[0] == 'q':
            return ('q', tuple(sorted((gp(n[1][0]), gp(n[1][1])))))
        return ('beta', g[n[1]])
    t = c[0]
    if t == 'CL':
        return ('CL', tuple(sorted({(s, gn(n)) for s, n in c[1]})))
    if t in ('EQ', 'LE', 'GE'):
        return (t, tuple(sorted((s, gn(n)) for s, n in c[1])), c[2])
    if t == 'IDENT':
        return ('IDENT', g[2 * c[1]] // 2, tuple(sorted(g[b] for b in c[2])), c[3])
    if t == 'PATH':
        return ('PATH', g[2 * c[1]] // 2, tuple(sorted(g[b] for b in c[2])))
    raise ValueError(c)


def invariant(C, g):
    img = Counter()
    for c, m in C.items():
        img[act(g, c)] += m
    return img == C


def spec_sha(C):
    return hashlib.sha256(json.dumps(sorted([[repr(c), m] for c, m in C.items()])).encode()).hexdigest()


def compose(a, b):                                   # (a o b)(l) = a[b[l]]
    return tuple(a[b[l]] for l in L)


def closure(gens):
    idp = tuple(L)
    seen = {idp}; frontier = [idp]
    while frontier:
        nxt = []
        for h in frontier:
            for g in gens:
                k = compose(g, h)
                if k not in seen:
                    seen.add(k); nxt.append(k)
        frontier = nxt
    return seen


def generators(G):
    gens, cl_ = [], {tuple(L)}
    for g in G:
        g = tuple(g)
        if g not in cl_:
            gens.append(g); cl_ = closure(gens)
        if len(cl_) == len(G):
            break
    return gens


def lex_perms(G, mode):
    G = [tuple(g) for g in G]
    if mode == 'all':
        return G
    if mode == 'gens':
        return generators(G)
    if mode.startswith('sample:'):
        _, n, seed = mode.split(':')
        n = min(int(n), len(G))
        return sorted(random.Random(int(seed)).sample(G, n))
    raise ValueError(mode)


FLAGS = ('cuts', 'digon', 'path', 'sec', 'twist', 'balance')


def parse_flags(argv):
    kw = {f: ('--' + f) in argv for f in FLAGS}
    kw['ident'] = next((int(a.split('=')[1]) if '=' in a else 0 for a in argv if a.startswith('--ident')), None)
    kw['beta'] = next(((a.split('=')[1] if '=' in a else 'card') for a in argv if a.startswith('--beta')), None)
    return kw


def build_all(mu, tau, kw, lexmode, prefix):
    """spec, canonical multiset, compiled CNF (+ LEX clauses), G, the permutations used by LEX, C1 file hashes"""
    PT, FS, shas = c1_modules()
    ev = build_spec(mu, tau, **kw)
    C = canon_spec(ev, mu, tau)
    cnf, pool = compile_spec(ev, mu, tau)
    G = PT.stabiliser(mu, tau)
    perms = lex_perms(G, lexmode)
    lex = []
    if len(G) > 1 and perms:
        cnt = itertools.count()
        lex = FS.lexleader_clauses(mu, lambda z, e: pool.obj2id[('x', z, e)], lambda: pool.id(('lex', next(cnt))),
                                   prefix_roots=prefix, G=[list(g) for g in perms])
    return ev, C, cnf, pool, G, perms, lex, shas


def export(types, i, out, argv):
    rec = types[i]; mu, tau = rec['mu'], rec['tau']
    kw = parse_flags(argv)
    lexmode = next((a.split('=', 1)[1] for a in argv if a.startswith('--lex=')), 'all')
    prefix = int(next((a.split('=')[1] for a in argv if a.startswith('--prefix=')), 6))
    with contextlib.redirect_stdout(io.StringIO()):
        ev, C, cnf, pool, G, perms, lex, shas = build_all(mu, tau, kw, lexmode, prefix)
    gens = generators(G)
    for g in gens:
        assert invariant(C, g), 'spec NOT invariant under a generator of Stab: refusing to export'
    full = CNF(); full.extend(cnf.clauses); full.extend(lex)
    head = {'mu': mu, 'tau': tau, 'flags': kw, 'lex': lexmode, 'prefix': prefix, 'G_order': len(G),
            'generators': [list(g) for g in gens], 'spec_sha256': spec_sha(C), 'base_clauses': len(cnf.clauses),
            'lex_clauses': len(lex), 'c1_files_sha256': shas}
    full.to_file(out, comments=['c LEXSPEC v1 (structures, spec_local.py): CNF = compile(SPEC) + LEX; verify with '
                                '"python spec_local.py verify FILE types.json i"', 'c LEXSPEC ' + json.dumps(head)])
    print('%s: |G| = %d, %d generators, spec invariant; base %d clauses + LEX %d (%s, prefix %d); vars %d'
          % (out, len(G), len(gens), len(cnf.clauses), len(lex), lexmode, prefix, pool.top))


def read_dimacs(path):
    head, clauses, cur = None, [], []
    for line in open(path):
        if line.startswith('c LEXSPEC {'):
            head = json.loads(line[len('c LEXSPEC '):])
        if line[0] in 'cp':
            continue
        for tok in line.split():
            v = int(tok)
            if v == 0:
                clauses.append(cur); cur = []
            else:
                cur.append(v)
    return head, clauses


def verify(path, types, i, gens_override=None, spec_mutator=None):
    """independent re-check; returns (ok, message)"""
    head, clauses = read_dimacs(path)
    rec = types[i]; mu, tau = rec['mu'], rec['tau']
    if head is None or head['mu'] != mu or head['tau'] != tau:
        return False, 'header missing or type mismatch'
    PT, FS, shas = c1_modules()
    if shas != head['c1_files_sha256']:
        return False, 'C1 source files differ from those used at export'
    kw = head['flags']
    with contextlib.redirect_stdout(io.StringIO()):
        ev = build_spec(mu, tau, **kw)
    if spec_mutator:
        ev = spec_mutator(ev)
    C = canon_spec(ev, mu, tau)
    if spec_sha(C) != head['spec_sha256']:
        return False, 'spec hash mismatch'
    cnf, pool = compile_spec(ev, mu, tau)
    G = [tuple(g) for g in PT.stabiliser(mu, tau)]
    gens = [tuple(g) for g in (gens_override or head['generators'])]
    Gs = set(G)
    if any(g not in Gs for g in gens):
        return False, 'a generator is not in Stab(mu, tau)'
    if closure(gens) != Gs:
        return False, 'the generators do not generate Stab(mu, tau)'
    for g in gens:
        if not invariant(C, g):
            return False, 'spec NOT invariant under generator %s' % (g,)
    perms = lex_perms(G, head['lex'])
    if any(tuple(g) not in Gs for g in perms):
        return False, 'LEX uses a permutation outside Stab(mu, tau)'
    lex = []
    if len(G) > 1 and perms:
        cnt = itertools.count()
        lex = FS.lexleader_clauses(mu, lambda z, e: pool.obj2id[('x', z, e)], lambda: pool.id(('lex', next(cnt))),
                                   prefix_roots=head['prefix'], G=[list(g) for g in perms])
    if [list(c) for c in cnf.clauses] + [list(c) for c in lex] != clauses:
        return False, 'the DIMACS clauses are not compile(SPEC) + LEX'
    return True, 'VERIFIED: CNF = compile(SPEC) + LEX; SPEC invariant under %d generators of Stab (|G| = %d); LEX over %d perms of G' % (
        len(gens), len(G), len(perms))


def selftest(part='all'):
    T = json.load(open('all_types_3308.json'))
    R = json.load(open('C:/home/researcher/Desktop/c1_recherche/data/allonblock_ge6_rel.json'))
    ok = True
    cases = [(T[i]['mu'], T[i]['tau']) for i in (0, 1, 674, 3307)] + [(R[k]['mu'], R[k]['tau']) for k in (2, 25)]
    flagsets = [dict(), dict(cuts=True, digon=True, path=True, twist=True, ident=1, balance=True, beta='card'),
                dict(sec=True, twist=True, fullcyc=True, beta='all', cuts=True)]
    if part in ('all', 'bytes'):                      # 1. compile(build_spec) == sat_local2.encode, byte for byte
        same = 0
        for mu, tau in cases:
            for fl in flagsets:
                with contextlib.redirect_stdout(io.StringIO()):
                    a, pa, _, _ = S2.encode(mu, tau, **fl)
                    b, pb = compile_spec(build_spec(mu, tau, **fl), mu, tau)
                same += (a.clauses == b.clauses and pa.top == pb.top)
        print('compile(spec) == encode byte for byte: %d / %d' % (same, len(cases) * len(flagsets)), flush=True)
        ok &= same == len(cases) * len(flagsets)
    if part in ('all', 'inv'):                        # 2. invariance under Stab; 3. negative controls
        PT, FS, shas = c1_modules()
        inv = tot = 0; neg_gen = neg_spec = nneg = 0
        rng = random.Random(3)
        allb6 = [tuple(g) for g in FS.b6()]
        for mu, tau in cases:
            G = PT.stabiliser(mu, tau)
            gens = generators(G)
            Gs = {tuple(g) for g in G}
            for fl in flagsets[1:]:
                with contextlib.redirect_stdout(io.StringIO()):
                    ev = build_spec(mu, tau, **fl)
                C = canon_spec(ev, mu, tau)
                tot += 1
                inv += all(invariant(C, g) for g in gens)
                if len(G) == 1:
                    continue                           # no symmetry: the negative controls do not apply
                nneg += 1
                bad = next((g for g in (rng.choice(allb6) for _ in range(1000)) if g not in Gs),
                           (2, 1, 0) + tuple(range(3, 12)))   # G = all of B6: use a permutation outside B6
                neg_gen += not invariant(C, bad)
                ev2 = list(ev)
                k = next(j for j, e in enumerate(ev2) if e[0] == 'EQ' and e[1])
                ev2[k] = ('EQ', ev2[k][1][:-1], ev2[k][2])     # drop one literal of one cardinality constraint
                C2 = canon_spec(ev2, mu, tau)
                neg_spec += not all(invariant(C2, g) for g in gens)
            print('  case |G| = %d, %d generators done' % (len(G), len(gens)), flush=True)
        print('spec invariant under Stab generators: %d / %d' % (inv, tot), flush=True)
        print('negative controls (cases with |G| > 1): generator outside G caught %d / %d; broken constraint caught %d / %d'
              % (neg_gen, nneg, neg_spec, nneg), flush=True)
        ok &= inv == tot and neg_gen == nneg and neg_spec == nneg
    if part in ('all', 'cycle'):                      # 4. export -> verify -> tamper -> verify fails
        import tempfile
        tmp = os.path.join(tempfile.gettempdir(), 'lexspec_selftest.cnf')
        rel = [{'mu': R[k]['mu'], 'tau': R[k]['tau']} for k in range(len(R))]
        k = 25                                         # a (0,6) type, |Stab| = 6
        export(rel, k, tmp, ['--twist', '--ident=1', '--balance', '--beta', '--lex=all', '--prefix=2'])
        ok1, m1 = verify(tmp, rel, k)
        print('  verify genuine :', m1, flush=True)
        lines = open(tmp).read().split(chr(10))
        j = next(n for n, l in enumerate(lines) if l and l[0] not in 'cp')
        tam = list(lines); tam[j] = '-' + tam[j] if not tam[j].startswith('-') else tam[j][1:]
        open(tmp + '.tam', 'w').write(chr(10).join(tam))
        ok2, m2 = verify(tmp + '.tam', rel, k)
        print('  verify tampered clause :', m2, flush=True)
        ok3, m3 = verify(tmp, rel, k, gens_override=[[1, 0] + list(range(2, 12))])
        print('  verify with a generator outside G :', m3, flush=True)
        ok4, m4 = verify(tmp, rel, k, spec_mutator=lambda ev: ev[:-1])
        print('  verify with a truncated spec :', m4, flush=True)
        ok &= ok1 and not ok2 and not ok3 and not ok4
    print('SELFTEST', part, 'PASSED' if ok else 'FAILED')
    return ok


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'selftest':
        sys.exit(0 if selftest(sys.argv[2] if len(sys.argv) > 2 else 'all') else 1)
    if cmd == 'export':
        export(json.load(open(sys.argv[2])), int(sys.argv[3]), sys.argv[4], sys.argv[5:])
    elif cmd == 'verify':
        ok, msg = verify(sys.argv[2], json.load(open(sys.argv[3])), int(sys.argv[4]))
        print(msg)
        sys.exit(0 if ok else 1)
