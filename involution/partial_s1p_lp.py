"""Same test with the step-1' model (near-design + block-level (ii); 19 680 variables): LP with k iota-pairs fixed."""
import json, sys, time, random, itertools
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
import iota_design as D
from fullblock_rigid import mates_from_D

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
m = mates_from_D([tuple(c) for c in cfg['cycles']]); mu = tau = m
K2, K3 = D.K2K3(mu, tau)
L = list(range(12))
edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
eidx = {e: i for i, e in enumerate(edges)}
E = len(edges); BP = list(itertools.combinations(range(6), 2))
pairs = [(i, j) for i in range(E) for j in range(i + 1, E) if not (set(edges[i]) & set(edges[j]))]
npair = len(pairs); pv = {pr: n for n, pr in enumerate(pairs)}
av = lambda n, k: npair + 15 * n + k
nv = npair * 16
ri, ci, vv, lo = [], [], [], []
r = [0]
def row(ent, a):
    for j, c in ent:
        ri.append(r[0]); ci.append(j); vv.append(c)
    lo.append(a); r[0] += 1
for i in range(E):
    row([(pv[(min(i, j), max(i, j))], 1) for j in range(E) if j != i and (min(i, j), max(i, j)) in pv], 1)
acc = {}
for n, (i, j) in enumerate(pairs):
    for b in edges[i]:
        for d in edges[j]:
            acc.setdefault((min(b, d), max(b, d)), []).append(n)
for b in L:
    for d in L:
        if d <= b: continue
        row([(n, 1) for n in acc.get((b, d), [])], int(K3[b, d]))
for n in range(npair):
    row([(av(n, k), 1) for k in range(15)] + [(n, -1)], 0)
for k in range(15):
    row([(av(n, k), 1) for n in range(npair)], 2)
for rr in range(6):
    ks = [k for k, bp in enumerate(BP) if rr in bp]
    for b in L:
        row([(av(n, k), 1) for n, (i, j) in enumerate(pairs) if b in edges[i] or b in edges[j] for k in ks],
            int(K2[2 * rr, b] + K2[2 * rr + 1, b]))
A = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r[0], nv)); b_ = np.array(lo, float)
t0 = time.time()
res = linprog(np.zeros(nv), A_eq=A, b_eq=b_, bounds=(0, 1), method='highs')
print("step-1' LP root:", {0: 'feasible', 2: 'INFEASIBLE'}.get(res.status, res.status), '%.1fs' % (time.time() - t0), flush=True)
rng = random.Random(0)
for k_iota in range(2):
    ip = [(eidx[tuple(e)], eidx[tuple(f)]) for e, f in iotas[k_iota]]
    order = list(range(30)); rng.shuffle(order)
    for kfix in (3, 6, 8, 10, 12, 15, 20):
        lb = np.zeros(nv)
        for q in order[:kfix]:
            i, j = ip[q]
            lb[pv[(min(i, j), max(i, j))]] = 1
        t0 = time.time()
        res = linprog(np.zeros(nv), A_eq=A, b_eq=b_, bounds=list(zip(lb, np.ones(nv))), method='highs')
        print("iota #%d, %2d fixed: step-1' LP" % (k_iota, kfix), {0: 'feasible', 2: 'INFEASIBLE'}.get(res.status, res.status),
              '%.2fs' % (time.time() - t0), flush=True)
