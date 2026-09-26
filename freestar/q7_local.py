"""C1's q:7 + RSF local test (26/09): at a block with t_c >= 3, q:7 forces a nonzero BALANCED torsion vector phi in ker Q_c:
    s_f . phi = 0 (mod 7) for the 30 far orbits f = {z, sigma z},  s_f = 1_{Phi(z)} - 1_{Phi(sigma z)} + 2(1_z - 1_{sigma z}).
Up to the stabilizer of (mu, tau), phi is one of C1's orbit representatives.  For each representative, solve the full
compact formulation (i)+(ii)+(iii) (full_local.py) plus these 30 congruences (linear with a free integer y_f per sigma-pair).
All infeasible  =>  no star of this type (q:7 + RSF).  Each run is capped at 5 minutes (solver rule)."""
import json, sys, time, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp
from full_local import K2, K3, sig, L
REPS = r"C:\Users\RESEARCHER\Desktop\c1_experimental\q7\A33_torsion_reps.json"

def solve_phi(mu, tau, phi, tl=300):
    k2, k3 = K2(mu, tau), K3(mu, tau)
    far = [e for e in itertools.combinations(L, 2) if sig[e[0]] != e[1]]
    reps = []; seen = set()
    for z in far:
        if z in seen: continue
        sz = tuple(sorted((sig[z[0]], sig[z[1]]))); seen |= {z, sz}; reps.append((z, sz))
    nonmu = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    def allowed(z, e): return all(k2[a, b] > 0 for a in z for b in e)
    var = []
    for Zi, (z, sz) in enumerate(reps):
        for e in nonmu:
            if not allowed(z, e): continue
            for e2 in nonmu:
                if e2 != e and allowed(sz, e2): var.append((Zi, e, e2))
    nv0 = len(var); nZ = len(reps); nv = nv0 + nZ          # + y_Z integers
    ri, ci, vv, lo = [], [], [], []; r = 0
    def row(entries, val):
        nonlocal r
        for j, c in entries: ri.append(r); ci.append(j); vv.append(c)
        lo.append(val); r += 1
    byZ, byE, acc, acc3 = {}, {}, {}, {}
    for j, (Zi, e, e2) in enumerate(var):
        z, sz = reps[Zi]
        byZ.setdefault(Zi, []).append(j); byE.setdefault(e, []).append(j); byE.setdefault(e2, []).append(j)
        for a in z:
            for b in e: acc.setdefault((a, b), []).append(j)
        for a in sz:
            for b in e2: acc.setdefault((a, b), []).append(j)
        for b in e:
            for d in e2: acc3.setdefault((b, d), []).append(j)
        for b in e2:
            for d in e: acc3.setdefault((b, d), []).append(j)
    for Zi in range(nZ): row([(j, 1) for j in byZ.get(Zi, [])], 1)
    for e in nonmu: row([(j, 1) for j in byE.get(e, [])], 1)
    for a in L:
        for b in L: row([(j, 1) for j in acc.get((a, b), [])], int(k2[a, b]))
    for b in L:
        for d in L: row([(j, 1) for j in acc3.get((b, d), [])], int(k3[b, d]))
    # q:7 balance: sum_{vars of Z} c_var v - 7 y_Z = -const_Z
    for Zi, (z, sz) in enumerate(reps):
        const = 2 * (sum(phi[x] for x in z) - sum(phi[x] for x in sz))
        ent = []
        for j in byZ.get(Zi, []):
            _, e, e2 = var[j]
            c = sum(phi[b] for b in e) - sum(phi[b] for b in e2)
            if c: ent.append((j, c))
        ent.append((nv0 + Zi, -7))
        row(ent, -const)
    A = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r, nv))
    v = np.array(lo, float)
    lb = np.concatenate([np.zeros(nv0), -8 * np.ones(nZ)]); ub = np.concatenate([np.ones(nv0), 8 * np.ones(nZ)])
    t0 = time.time()
    res = milp(np.zeros(nv), constraints=[LinearConstraint(A, v, v)], integrality=np.ones(nv), bounds=Bounds(lb, ub),
               options={'time_limit': tl})
    st = {0: 'FEASIBLE', 2: 'INFEASIBLE'}.get(res.status, 'unknown(%d)' % res.status)
    sol = None
    if res.status == 0:
        sol = {}
        for j in np.nonzero(res.x[:nv0] > 0.5)[0]:
            Zi, e, e2 = var[j]; z, sz = reps[Zi]
            sol[str(list(z))] = list(e); sol[str(list(sz))] = list(e2)
    return st, sol, time.time() - t0

if __name__ == '__main__':
    d = json.load(open(REPS))
    mu = {i: d['mu'][i] for i in L}; tau = {i: d['tau'][i] for i in L}
    which = sys.argv[1] if len(sys.argv) > 1 else 'mixed'
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 300
    out = []
    for k, rep in enumerate(d['torsion_reps']):
        if which == 'mixed' and rep['single_component']: continue
        st, sol, dt = solve_phi(mu, tau, rep['phi'], tl)
        print('rep', k, 'phi', rep['phi'], 'single' if rep['single_component'] else 'mixed', '->', st, '%.0fs' % dt, flush=True)
        out.append({'rep': k, 'phi': rep['phi'], 'status': st, 'seconds': round(dt), 'Phi': sol})
        json.dump(out, open('q7_local_A33_results.json', 'w'), indent=1)
