"""Root LP of the FULL local system (structures 27/09): iota-pair variables p[n] AND sigma-pair assignments x[Z,n,o] together,
with label-level (ii) (step 2) instead of the block-level (ii) of step 1'.  Constraints: each sigma-pair gets one oriented
iota-pair; sum_{Z,o} x[Z,n,o] = p[n]; iota perfect matching of the non-mu pairs; K3 cross counts on p; label-level (ii):
sum x [a in root][b in Phi(root)] = K2[a,b]; pointwise L1/S1/S2 as variable filters.  Infeasible root => the type is dead.
Usage: python full_lp_root.py types.json i   (types.json: list of {mu, tau} records)"""
import json, sys, itertools, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import iota_design as D
T = json.load(open(sys.argv[1])); i0 = int(sys.argv[2]); mu, tau = T[i0]['mu'], T[i0]['tau']
L = range(12); PI = [v ^ 1 for v in L]
K2, K3 = D.K2K3(mu, tau)
edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]; E = len(edges)
pairs = [(i, j) for i in range(E) for j in range(i + 1, E) if not any(K3[b, b] < 2 for b in set(edges[i]) & set(edges[j]))]
far = [z for z in itertools.combinations(L, 2) if z[0] // 2 != z[1] // 2]
reps, seen = [], set()
for z in far:
    if z not in seen:
        sz = tuple(sorted((PI[z[0]], PI[z[1]]))); seen |= {z, sz}; reps.append((z, sz))
def pw(z, e, f):
    for a in z:
        for aa in (a, PI[a]):
            if aa in e and PI[aa] in f: return False
        if a in e and tau[a] in f: return False
        if a not in e and mu[a] in e and tau[a] in f: return False
    return True
X = []
for Zi, (z, sz) in enumerate(reps):
    for n, (i, j) in enumerate(pairs):
        for o, (e, f) in enumerate(((edges[i], edges[j]), (edges[j], edges[i]))):
            if pw(z, e, f) and pw(sz, f, e):
                X.append((Zi, n, e, f))
npair = len(pairs); nx = len(X); nv = npair + nx
ri, ci, vv, rhs = [], [], [], []; r = [0]
def row(ent, b):
    for j, c in ent:
        ri.append(r[0]); ci.append(j); vv.append(c)
    rhs.append(b); r[0] += 1
byZ, byN, acc2 = {}, {}, {}
for k, (Zi, n, e, f) in enumerate(X):
    byZ.setdefault(Zi, []).append(npair + k); byN.setdefault(n, []).append(npair + k)
    z, sz = reps[Zi]
    for a in z:
        for b in e: acc2.setdefault((a, b), []).append(npair + k)
    for a in sz:
        for b in f: acc2.setdefault((a, b), []).append(npair + k)
for Zi in range(30): row([(j, 1) for j in byZ.get(Zi, [])], 1)
for n in range(npair): row([(j, 1) for j in byN.get(n, [])] + [(n, -1)], 0)
for ei in range(E): row([(n, 1) for n, p_ in enumerate(pairs) if ei in p_], 1)
acc3 = {}
for n, (i, j) in enumerate(pairs):
    for b in edges[i]:
        for d in edges[j]:
            acc3.setdefault((min(b, d), max(b, d)), []).append(n)
for b in L:
    for d in L:
        if d >= b: row([(n, 1) for n in acc3.get((b, d), [])], int(K3[b, d]) if b != d else int(K3[b, b]) // 2)
for a in L:
    for b in L: row([(j, 1) for j in acc2.get((a, b), [])], int(K2[a, b]))
A = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r[0], nv))
t0 = time.time()
b_eq = np.array(rhs, float)
res = linprog(np.zeros(nv), A_eq=A, b_eq=b_eq, bounds=(0, 1), method='highs')
print('i %d mu %s: vars %d (iota %d + assign %d), rows %d | root LP: %s (%.1fs)' % (i0, mu, nv, npair, nx, r[0],
      {0: 'FEASIBLE', 2: 'INFEASIBLE'}.get(res.status, res.status), time.time() - t0), flush=True)
if any(a.startswith('--depth') for a in sys.argv):                # LP death depth along true near-designs (step 1), random fixing order
    import random
    pv = {p_: n for n, p_ in enumerate(pairs)}; rng = random.Random(0); found = []
    for it in range(int(next((a.split('=')[1] for a in sys.argv if a.startswith('--depth=')), 3))):
        st, sol, _ = D.step1(mu, tau, 1, 40, forbid=found)
        if sol is None:
            print('step1', st); break
        found.append(sol); order = list(range(30)); rng.shuffle(order); out = []
        for k in range(0, 31):
            lb = np.zeros(nv)
            for q in order[:k]:
                i, j = sol[q]; lb[pv[(min(i, j), max(i, j))]] = 1
            st2 = linprog(np.zeros(nv), A_eq=A, b_eq=b_eq, bounds=list(zip(lb, np.ones(nv))), method='highs').status
            if st2 == 2:
                out.append('dies at %d' % k); break
        else:
            out.append('survives all 30')
        print('  iota #%d: %s' % (it, out[0]), flush=True)
if any(a.startswith('--milp') for a in sys.argv):   # the whole model as a MILP (HiGHS), time-limited
    from scipy.optimize import milp, LinearConstraint, Bounds
    tl = float(next((a.split('=')[1] for a in sys.argv if a.startswith('--milp=')), 200))
    t0 = time.time()
    rm = milp(np.zeros(nv), constraints=[LinearConstraint(A, b_eq, b_eq)], integrality=np.ones(nv), bounds=Bounds(0, 1),
              options={'time_limit': tl, 'disp': False})
    print('  full MILP: %s (%.1fs)' % ({0: 'FEASIBLE', 2: 'INFEASIBLE'}.get(rm.status, 'status %s: %s' % (rm.status, rm.message)),
                                       time.time() - t0), flush=True)
