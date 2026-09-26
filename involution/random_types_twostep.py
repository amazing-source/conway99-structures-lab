"""Two-step test (step 1 near-design / iota with K3, step 2 = (ii) as an assignment, LP then MILP) on RANDOM block types
(mu, tau) with K2, K3 >= 0 (C2 notation).  Evidence for 'the local system is infeasible for every type'."""
import random, sys, time, json
import numpy as np
import iota_design as D
L = list(range(12)); PI = [v ^ 1 for v in L]
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
ntypes = int(sys.argv[2]) if len(sys.argv) > 2 else 4


def rand_fpf():
    pts = L[:]; rng.shuffle(pts); m = [0] * 12
    for i in range(0, 12, 2):
        m[pts[i]] = pts[i + 1]; m[pts[i + 1]] = pts[i]
    return m
done = 0
t0 = time.time()
while done < ntypes:
    mu, tau = rand_fpf(), rand_fpf()
    K2, K3 = D.K2K3(mu, tau)
    if (K2 < 0).any() or (K3 < 0).any():
        continue
    done += 1
    kind = 'full' if mu == tau else ('free' if all(mu[v] != tau[v] for v in L) else 'mixed')
    st, sol, edges = D.step1(mu, tau, 1, 60)
    if sol is None:
        print('type %d (%s): step 1 %s' % (done, kind, st), '%.0fs' % (time.time() - t0), flush=True)
        continue
    st2, nv = D.step2(mu, tau, sol, edges, tl=60)
    print('type %d (%s): step 1 FEASIBLE, step 2 %s' % (done, kind, st2), '%.0fs' % (time.time() - t0), flush=True)
    json.dump({'mu': mu, 'tau': tau}, open('rand_type_%s_%d.json' % (sys.argv[1] if len(sys.argv) > 1 else '0', done), 'w'))
