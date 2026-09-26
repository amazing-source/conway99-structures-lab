"""Code control of STEP 2 (iota_design.step2 logic, generalised to 2(m-1) labels) on BvLS (m = 11): with the TRUE iota of
the real graph at a block, the step-2 LP / MILP ((ii) as an assignment sigma-pairs -> oriented iota-pairs) must be feasible.
Also: a RANDOM re-pairing test is not meaningful here; we only check that the true instance is accepted."""
import itertools, sys
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
exec(open('bvls_local_control.py').read().split("bad_total = 0")[0])       # reuse the BvLS set-up (graph, blocks, roots)
ci = 0
cp, cm = blocks[ci]
others = [b for j, b in enumerate(blocks) if j != ci]
lab_of = {}
for i, (up, um) in enumerate(others):
    lab_of[up] = 2 * i; lab_of[um] = 2 * i + 1
Lr = list(range(20)); PI = [v ^ 1 for v in Lr]
croot = {}; cmroot = {}
for r in roots:
    a, b = rnb[r]
    if cp in (a, b): croot[lab_of[b if a == cp else a]] = r
    if cm in (a, b): cmroot[lab_of[b if a == cm else a]] = r
mu = [next(w for w in Lr if w != v and A[croot[v], croot[w]]) for v in Lr]
h = [next(u for u in Lr if A[croot[v], cmroot[u]]) for v in Lr]
tau = [PI[h[v]] for v in Lr]
far = [r for r in roots if cp not in rnb[r] and cm not in rnb[r]]
lab_root = {r: tuple(sorted(lab_of[u] for u in rnb[r])) for r in far}
Phi = {lab_root[r]: tuple(sorted(v for v in Lr if A[r, croot[v]])) for r in far}
sig = {z: tuple(sorted((PI[z[0]], PI[z[1]]))) for z in Phi}
K2 = np.full((20, 20), 2)
for a in Lr:
    K2[a, a] -= 1; K2[a, PI[a]] -= 1; K2[a, mu[a]] -= 1; K2[a, tau[PI[a]]] -= 1
reps, seen = [], set()
for z in sorted(Phi):
    if z in seen: continue
    seen |= {z, sig[z]}; reps.append((z, sig[z]))
iota_pairs = []
for z, sz in reps:
    iota_pairs.append((Phi[z], Phi[sz]))
# shuffle the iota-pairs' order so that the identity assignment is not trivially index-aligned
rng = np.random.default_rng(1); perm = rng.permutation(len(iota_pairs)); iota_pairs = [iota_pairs[i] for i in perm]
var = [(Zi, Pj, e, f) for Zi in range(len(reps)) for Pj, (e0, f0) in enumerate(iota_pairs) for (e, f) in ((e0, f0), (f0, e0))]
nv = len(var)
ri, ci_, vv, lo = [], [], [], []
r_ = [0]
def row(ent, val):
    for j, c in ent:
        ri.append(r_[0]); ci_.append(j); vv.append(c)
    lo.append(val); r_[0] += 1
byZ, byP, acc = {}, {}, {}
for k, (Zi, Pj, e, f) in enumerate(var):
    byZ.setdefault(Zi, []).append(k); byP.setdefault(Pj, []).append(k)
    z, sz = reps[Zi]
    for a in z:
        for b in e: acc.setdefault((a, b), []).append(k)
    for a in sz:
        for b in f: acc.setdefault((a, b), []).append(k)
for Zi in range(len(reps)): row([(k, 1) for k in byZ[Zi]], 1)
for Pj in range(len(iota_pairs)): row([(k, 1) for k in byP[Pj]], 1)
for a in Lr:
    for b in Lr: row([(k, 1) for k in acc.get((a, b), [])], int(K2[a, b]))
Am = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci_))), shape=(r_[0], nv))
lp = linprog(np.zeros(nv), A_eq=Am, b_eq=np.array(lo, float), bounds=(0, 1), method='highs')
ip = milp(np.zeros(nv), constraints=[LinearConstraint(Am, lo, lo)], integrality=np.ones(nv), bounds=Bounds(0, 1),
          options={'time_limit': 120})
print('BvLS block 0, TRUE iota: step-2 LP', {0: 'feasible', 2: 'INFEASIBLE'}.get(lp.status, lp.status),
      '| MILP', {0: 'feasible', 2: 'INFEASIBLE'}.get(ip.status, ip.status), '| vars', nv)
