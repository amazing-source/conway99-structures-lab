"""Two-step search for a LOCAL SOLUTION at any block type: step 1 (CP-SAT) gives a near-design iota satisfying (i)+(iii);
step 2 is then solved EXACTLY as a small SAT problem (CaDiCaL): y[Z, P, o] = sigma-pair Z -> iota-pair P with orientation o,
exactly one per Z and per P, (ii) as cardinality equalities, C1's pointwise conditions as variable filters.
SAT  => a full local solution of (i)+(ii)+(iii)(+pointwise) for this type (saved);  UNSAT => this iota dies.
Usage: python step2_sat.py cfg.json n_iota step1_seconds step2_seconds"""
import itertools, json, sys, time, threading
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
import iota_design as D

cfg = json.load(open(sys.argv[1])); n_iota = int(sys.argv[2]); t1 = float(sys.argv[3]); t2 = float(sys.argv[4])
if 'cycles' in cfg:
    from fullblock_rigid import mates_from_D
    m = mates_from_D([tuple(c) for c in cfg['cycles']]); mu = tau = m
else:
    mu, tau = cfg['mu'], cfg['tau']
K2, K3 = D.K2K3(mu, tau)
L = D.L; PI = D.PI
far = [z for z in itertools.combinations(L, 2) if z[0] // 2 != z[1] // 2]
reps, seen = [], set()
for z in far:
    if z in seen:
        continue
    sz = tuple(sorted((PI[z[0]], PI[z[1]]))); seen |= {z, sz}; reps.append((z, sz))


def pointwise_ok(z, e, f):
    se, sf = set(e), set(f)
    for a in z:
        for aa in (a, PI[a]):
            if aa in se and PI[aa] in sf:
                return False
        if a in se and tau[a] in sf:
            return False
        if a not in se and mu[a] in se and tau[a] in sf:
            return False
    return True
found = []
t0 = time.time()
for k in range(n_iota):
    st, sol, edges = D.step1(mu, tau, 1, t1, forbid=found)
    print('step 1 #%d' % k, st, '%.0fs' % (time.time() - t0), flush=True)
    if sol is None:
        print('step 1', st, '%.0fs' % (time.time() - t0), flush=True)
        break
    found.append(sol)
    pool = IDPool(); cnf = CNF()
    Y = {}
    for Zi, (z, sz) in enumerate(reps):
        for Pj, (i, j) in enumerate(sol):
            for o, (e, f) in enumerate(((edges[i], edges[j]), (edges[j], edges[i]))):
                if pointwise_ok(z, e, f) and pointwise_ok(sz, f, e):
                    Y[Zi, Pj, o] = (pool.id(('y', Zi, Pj, o)), e, f)

    def eq(lits, kk):
        if kk < 0 or kk > len(lits):
            cnf.append([]); return
        cnf.extend(CardEnc.equals(lits=lits, bound=kk, vpool=pool, encoding=EncType.seqcounter).clauses)
    for Zi in range(30):
        eq([v for (a, b, c), (v, e, f) in Y.items() if a == Zi], 1)
    for Pj in range(30):
        eq([v for (a, b, c), (v, e, f) in Y.items() if b == Pj], 1)
    acc = {}
    for (Zi, Pj, o), (v, e, f) in Y.items():
        z, sz = reps[Zi]
        for a in z:
            for b in e: acc.setdefault((a, b), []).append(v)
        for a in sz:
            for b in f: acc.setdefault((a, b), []).append(v)
    for a in L:
        for b in L:
            eq(acc.get((a, b), []), int(K2[a, b]))
    with Solver(name='glucose4', bootstrap_with=cnf.clauses) as s:
        timer = threading.Timer(t2, lambda: s.interrupt()); timer.start()
        res = s.solve_limited(expect_interrupt=True); timer.cancel()
        tag = {True: 'SAT', False: 'UNSAT', None: 'UNKNOWN'}[res]
        print('iota #%d: step 2 SAT solver: %s' % (k, tag), '| y vars', len(Y), '%.0fs' % (time.time() - t0), flush=True)
        if res:
            model = set(l for l in s.get_model() if l > 0)
            asg = [[list(reps[Zi][0]), list(e), list(f)] for (Zi, Pj, o), (v, e, f) in Y.items() if v in model]
            json.dump({'cfg': cfg, 'Phi (z -> Phi z, Phi sigma z)': asg}, open(sys.argv[1].replace('.json', '_LOCALSOL.json'), 'w'), indent=1)
            print('   !!! LOCAL SOLUTION FOUND and saved', flush=True)
            break
