"""Quick triage of the 79 relevant types (C1's allonblock_ge6_rel.json) with the single-block encoding 2 and every flag
that applies (--cuts --digon --path --sec, + --twist for twisted blocks), glucose4 with a short time limit per type.
Laptop use: a slice of indices per process (<= 5 min). Output: one line per type, appended to sweep79.log.
Usage: python sweep79.py first last seconds [--ident=s]"""
import json, sys, time, threading
from pysat.solvers import Solver
import sat_local2 as S2

rel = json.load(open('C:/home/researcher/Desktop/c1_recherche/data/allonblock_ge6_rel.json'))
a, b, secs = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
ident = next((int(x.split('=')[1]) for x in sys.argv if x.startswith('--ident=')), None)
for i in range(a, min(b, len(rel))):
    d = rel[i]
    mu, tau = d['mu'], d['tau']
    twisted = any(sum(mu[l] == tau[l] for l in (2 * p, 2 * p + 1)) == 1 for p in range(6))
    t0 = time.time()
    cnf, pool, far, edges = S2.encode(mu, tau, cuts=True, digon=True, path=True, sec=True, twist=twisted, ident=ident)
    with Solver(name='glucose4', bootstrap_with=cnf.clauses) as s:
        timer = threading.Timer(secs, lambda: s.interrupt()); timer.start()
        t1 = time.time()
        res = s.solve_limited(expect_interrupt=True); timer.cancel()
    line = 'rel %2d profile %s  %s  %.1fs  (vars %d, clauses %d, flags cuts digon path sec%s%s)' % (
        i, d.get('anch_tw'), {True: 'SAT', False: 'UNSAT', None: 'unknown'}[res], time.time() - t1, pool.top,
        len(cnf.clauses), ' twist' if twisted else '', ' ident<=%d' % ident if ident is not None else '')
    print(line, flush=True)
    with open('sweep79.log', 'a') as f:
        f.write(line + '\n')
