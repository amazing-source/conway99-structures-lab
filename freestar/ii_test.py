"""Free star, constraints (i)+(ii) only: Phi = bijection E(CP(sigma)) -> E(CP(mu)) with |Phi(star_a) n star_b| = K2[a,b],
K2 = 2J - I - P_s - P_mu - P_{tau s} (entrywise, overlapping indicators subtract separately).  ILP on the 60x60 permutation."""
import itertools, random, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp
L = list(range(12)); sig = {i: i ^ 1 for i in L}
def rand_fpf(rng):
    xs = L[:]; rng.shuffle(xs); f = {}
    for i in range(0, 12, 2): f[xs[i]], f[xs[i + 1]] = xs[i + 1], xs[i]
    return f
def K2(mu, tau):
    K = np.full((12, 12), 2); K -= np.eye(12, dtype=int)
    for a in L:
        K[a, sig[a]] -= 1; K[a, mu[a]] -= 1; K[a, tau[sig[a]]] -= 1
    return K
def solve_ii(mu, tau, tl=120):
    far = [e for e in itertools.combinations(L, 2) if sig[e[0]] != e[1]]
    nonmu = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    n = 60; var = lambda i, j: i * n + j
    rows = []; vals = []
    ri = []; ci = []; dv = []; r = 0
    for i in range(n):
        for j in range(n): ri.append(r); ci.append(var(i, j)); dv.append(1)
        vals.append(1); r += 1
    for j in range(n):
        for i in range(n): ri.append(r); ci.append(var(i, j)); dv.append(1)
        vals.append(1); r += 1
    K = K2(mu, tau)
    for a in L:
        for b in L:
            for i, z in enumerate(far):
                if a not in z: continue
                for j, e in enumerate(nonmu):
                    if b in e: ri.append(r); ci.append(var(i, j)); dv.append(1)
            vals.append(K[a, b]); r += 1
    A = sp.csr_matrix((dv, (ri, ci)), shape=(r, n * n))
    v = np.array(vals, float)
    res = milp(np.zeros(n * n), constraints=[LinearConstraint(A, v, v)], integrality=np.ones(n * n), bounds=Bounds(0, 1),
               options={'time_limit': tl})
    return {0: 'feasible', 2: 'INFEASIBLE'}.get(res.status, 'unknown(%d)' % res.status)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
from collections import Counter
cnt = Counter(); t0 = time.time()
for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 5):
    while True:
        mu, tau = rand_fpf(rng), rand_fpf(rng)
        if any(mu[b] == tau[b] for b in L): continue           # fully free
        break
    r = solve_ii(mu, tau)
    cnt[r] += 1
    print(t, r, '%.0fs' % (time.time() - t0), flush=True)
print(dict(cnt))
