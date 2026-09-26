"""Step 1' as a MILP (HiGHS): same model as step1prime.py (near-design + block-level (ii)); LP-based pruning.
Usage: python step1prime_milp.py cfg.json time_limit"""
import itertools, json, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp
import iota_design as D

L = list(range(12))
cfg = json.load(open(sys.argv[1])); tl = float(sys.argv[2])
if 'cycles' in cfg:
    from fullblock_rigid import mates_from_D
    m = mates_from_D([tuple(c) for c in cfg['cycles']]); mu = tau = m
else:
    mu, tau = cfg['mu'], cfg['tau']
K2, K3 = D.K2K3(mu, tau)
edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
E = len(edges); BP = list(itertools.combinations(range(6), 2))
pairs = [(i, j) for i in range(E) for j in range(i + 1, E) if not any(K3[b, b] < 2 for b in set(edges[i]) & set(edges[j]))]
npair = len(pairs)
pv = {pr: n for n, pr in enumerate(pairs)}
av = lambda n, k: npair + 15 * n + k
nv = npair * 16
ri, ci, vv, lo, hi = [], [], [], [], []
r = [0]
def row(ent, a, b):
    for j, c in ent:
        ri.append(r[0]); ci.append(j); vv.append(c)
    lo.append(a); hi.append(b); r[0] += 1
for i in range(E):
    row([(pv[(min(i, j), max(i, j))], 1) for j in range(E) if j != i and (min(i, j), max(i, j)) in pv], 1, 1)
acc = {}
for n, (i, j) in enumerate(pairs):
    for b in edges[i]:
        for d in edges[j]:
            acc.setdefault((min(b, d), max(b, d)), []).append(n)
for b in L:
    for d in L:
        if d < b: continue
        t = int(K3[b, d]) if b != d else int(K3[b, b]) // 2
        row([(n, 1) for n in acc.get((b, d), [])], t, t)
for n in range(npair):
    row([(av(n, k), 1) for k in range(15)] + [(n, -1)], 0, 0)
for k in range(15):
    row([(av(n, k), 1) for n in range(npair)], 2, 2)
for rr in range(6):
    ks = [k for k, bp in enumerate(BP) if rr in bp]
    for b in L:
        ent = [(av(n, k), 1) for n, (i, j) in enumerate(pairs) if b in edges[i] or b in edges[j] for k in ks]
        t = int(K2[2 * rr, b] + K2[2 * rr + 1, b])
        row(ent, t, t)
if 'cycles' in cfg and len(sys.argv) > 3 and sys.argv[3] == 'lemma':      # C1's full-block lemma as cuts
    for C in cfg['cycles']:
        for p_ in C:
            Sp = {2 * p_, 2 * p_ + 1, m[2 * p_], m[2 * p_ + 1]}
            big = [n for n, (i, j) in enumerate(pairs) if len((set(edges[i]) | set(edges[j])) & Sp) >= 3]
            four = [n for n, (i, j) in enumerate(pairs) if len((set(edges[i]) | set(edges[j])) & Sp) == 4]
            for n in four:
                row([(n, 1)], 0, 0)
            row([(n, 1) for n in big], 0, 1 if len(C) >= 3 else 0)
    # per block p: #{Q in S_p : k_p(Q) = 0} - sum_{Q in S_p} (k_p(Q) - 1)^+ = 2, i.e. sum_{Q in S_p} (1 - k_p(Q)) = 2
    for p_ in range(6):
        Sp = {2 * p_, 2 * p_ + 1, m[2 * p_], m[2 * p_ + 1]}
        ks = [k for k, bp in enumerate(BP) if p_ in bp]
        row([(av(n, k), 1 - len((set(edges[i]) | set(edges[j])) & Sp)) for n, (i, j) in enumerate(pairs) for k in ks], 2, 2)
A = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r[0], nv))
t0 = time.time()
res = milp(np.zeros(nv), constraints=[LinearConstraint(A, np.array(lo, float), np.array(hi, float))],
           integrality=np.ones(nv), bounds=Bounds(0, 1), options={'time_limit': tl, 'disp': False})
print(sys.argv[1], "step 1' MILP:", {0: 'FEASIBLE', 2: 'INFEASIBLE'}.get(res.status, 'unknown(%d)' % res.status),
      '%.1fs' % (time.time() - t0), '| vars', nv, '| rows', r[0], flush=True)
