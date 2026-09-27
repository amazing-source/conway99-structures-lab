"""Is a sib class's local system (i)-(iii) + L1/S1/S2 SATISFIABLE at m = 7? (structures 27/09, rung 1 check)
Two-step: near-design iota (step 1, CP-SAT), then step 2 (MILP: sigma-pairs -> oriented iota-pairs, (ii) exact) with the
pointwise conditions filtered in (both roots of each sigma-pair).  Usage: python sib_twostep.py i[,i...] seconds n_iota"""
import json, sys, time
import iota_design as D
T = json.load(open('C:/home/researcher/Desktop/c2_experimental_ro/scripts/c1_tiers/data/all_ge6.json'))
PI = [v ^ 1 for v in range(12)]
def pointwise_ok(mu, tau, z, e, f):
    for a in z:
        for aa in (a, PI[a]):
            if aa in e and PI[aa] in f:
                return False
        if a in e and tau[a] in f:
            return False
        if a not in e and mu[a] in e and tau[a] in f:
            return False
    return True
tl = float(sys.argv[2]); nio = int(sys.argv[3])
for i in map(int, sys.argv[1].split(',')):
    mu, tau = T[i]['mu'], T[i]['tau']
    ex = lambda z, sz, e, f: pointwise_ok(mu, tau, z, e, f) and pointwise_ok(mu, tau, sz, f, e)
    found = []; t0 = time.time()
    for k in range(nio):
        st, sol, edges = D.step1(mu, tau, 1, tl, forbid=found)
        if sol is None:
            print('i %d step1 #%d %s (%.0fs)' % (i, k, st, time.time() - t0), flush=True); break
        found.append(sol)
        st2, nv = D.step2(mu, tau, sol, edges, tl=tl, extras=ex)
        print('i %d %s: iota #%d -> step 2 with pointwise: %s (vars %d, %.0fs)' % (i, tuple(T[i]['anch_tw']), k, st2, nv, time.time() - t0), flush=True)
