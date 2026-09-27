"""Step 2 (assignment sigma-pairs -> oriented iota-pairs under (ii)) in CP-SAT for one iota of a saved type; several iota."""
import json, sys, time, itertools
from ortools.sat.python import cp_model
import iota_design as D
cfg = json.load(open(sys.argv[1])); mu, tau = cfg['mu'], cfg['tau']; n_iota = int(sys.argv[2]); tl = float(sys.argv[3])
K2, K3 = D.K2K3(mu, tau)
far = [z for z in itertools.combinations(D.L, 2) if z[0] // 2 != z[1] // 2]
reps, seen = [], set()
for z in far:
    if z in seen: continue
    sz = tuple(sorted((D.PI[z[0]], D.PI[z[1]]))); seen |= {z, sz}; reps.append((z, sz))
found = []; t0 = time.time()
for k in range(n_iota):
    st, sol, edges = D.step1(mu, tau, 1, 60, forbid=found)
    if sol is None:
        print('step1', st); break
    found.append(sol)
    md = cp_model.CpModel()
    y = {}
    for Zi in range(30):
        for Pj, (i, j) in enumerate(sol):
            for o, (e, f) in enumerate(((edges[i], edges[j]), (edges[j], edges[i]))):
                y[Zi, Pj, o] = (md.NewBoolVar(''), e, f)
    for Zi in range(30):
        md.AddExactlyOne(v for (a, b, c), (v, e, f) in y.items() if a == Zi)
    for Pj in range(30):
        md.AddExactlyOne(v for (a, b, c), (v, e, f) in y.items() if b == Pj)
    acc = {}
    for (Zi, Pj, o), (v, e, f) in y.items():
        z, sz = reps[Zi]
        for a in z:
            for b in e: acc.setdefault((a, b), []).append(v)
        for a in sz:
            for b in f: acc.setdefault((a, b), []).append(v)
    for a in D.L:
        for b in D.L:
            md.Add(sum(acc.get((a, b), [])) == int(K2[a, b]))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 2
    r = s.Solve(md)
    print('iota #%d: step 2 CP-SAT %s' % (k, s.StatusName(r)), '%.0fs' % (time.time() - t0), flush=True)
    if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        print('   !!! a full LOCAL SOLUTION of (i)+(ii)+(iii) exists for this type', flush=True)
        json.dump({'cfg': cfg, 'iota': [[list(edges[i]), list(edges[j])] for i, j in sol],
                   'assignment': [[list(reps[Zi][0]), list(e), list(f)] for (Zi, Pj, o), (v, e, f) in y.items() if s.Value(v)]},
                  open(sys.argv[1].replace('.json', '_LOCALSOL_%d.json' % k), 'w'), indent=1)
