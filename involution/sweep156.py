"""Quick triage for rung 1 of the ladder (structures 27/09): the 156 classes with >= 6 pivots (C2's all_ge6.json), encoding 2
with ONLY the generally valid flags of lemma_scope_156.json (no SEC), glucose4 with a short limit per class.
Usage: python sweep156.py first last seconds      (positions in lemma_scope_156.json 'order'); appends to sweep156.log"""
import json, sys, time, threading
from pysat.solvers import Solver
import sat_local2 as S2

T = json.load(open('C:/home/researcher/Desktop/c2_experimental_ro/scripts/c1_tiers/data/all_ge6.json'))
D = json.load(open('lemma_scope_156.json'))
byi = {o['i']: o for o in D['classes']}
a, b, secs = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
for pos in range(a, min(b, len(D['order']))):
    i = D['order'][pos]; o = byi[i]; t = T[i]
    fl = o['flags_general']
    kw = dict(cuts='--cuts' in fl, digon='--digon' in fl, path='--path' in fl, twist='--twist' in fl,
              ident=1 if '--ident=1' in fl else None)
    cnf, pool, far, edges = S2.encode(t['mu'], t['tau'], **kw)
    with Solver(name='glucose4', bootstrap_with=cnf.clauses) as s:
        timer = threading.Timer(secs, lambda: s.interrupt()); timer.start()
        t1 = time.time(); res = s.solve_limited(expect_interrupt=True); timer.cancel()
    line = 'pos %3d  i %3d  %s  rel %s  %s  %.1fs  flags %s' % (pos, i, tuple(o['anch_tw']), o['relevant_index'],
                                                            {True: 'SAT', False: 'UNSAT', None: 'unknown'}[res],
                                                            time.time() - t1, ' '.join(fl))
    print(line, flush=True)
    with open('sweep156.log', 'a') as f:
        f.write(line + '\n')
