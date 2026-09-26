"""Free star (C2 note 20 par. 3), sub-problem (iii) alone.  Phi is a bijection from the 60 edges of K12 - sigma (far roots) onto
the 60 edges of K12 - mu (non-mu pairs); rho := Phi o sigma o Phi^-1 is a fixed-point-free involution of the edges of K12 - mu,
and (iii) says:  #{edges e containing b : d in rho(e)} = K3[b,d] = 2J - P_sigma - P_{tau mu} - P_{mu tau} - P_tau.
Test: does such a rho exist?  (small ILP: pairs of edges.)"""
import itertools, random, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
L = list(range(12)); sig = {i: i ^ 1 for i in L}          # labels 2e, 2e+1 = (e,+), (e,-)
def rand_fpf(rng):
    xs = L[:]; rng.shuffle(xs); f = {}
    for i in range(0, 12, 2): f[xs[i]], f[xs[i + 1]] = xs[i + 1], xs[i]
    return f
def K3(mu, tau):
    K = np.full((12, 12), 2)
    for b in L:
        for f in (sig, {x: tau[mu[x]] for x in L}, {x: mu[tau[x]] for x in L}, tau):
            K[b, f[b]] -= 1
    return K
def rho_exists(mu, tau, tl=60):
    edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    E = len(edges); assert E == 60
    pairs = list(itertools.combinations(range(E), 2))
    K = K3(mu, tau)
    rows = []; lo = []; hi = []
    # each edge in exactly one pair
    A1 = np.zeros((E, len(pairs)))
    for k, (e, f) in enumerate(pairs): A1[e, k] = 1; A1[f, k] = 1
    # incidence counts
    A2 = np.zeros((144, len(pairs)))
    for k, (e, f) in enumerate(pairs):
        for b in edges[e]:
            for d in edges[f]: A2[b * 12 + d, k] += 1
        for b in edges[f]:
            for d in edges[e]: A2[b * 12 + d, k] += 1
    A = np.vstack([A1, A2]); rhs = np.concatenate([np.ones(E), K.flatten()])
    res = milp(c=np.zeros(len(pairs)), constraints=[LinearConstraint(A, rhs, rhs)], integrality=np.ones(len(pairs)),
               bounds=Bounds(0, 1), options={'time_limit': tl})
    return {0: 'feasible', 2: 'INFEASIBLE'}.get(res.status, 'unknown(%d)' % res.status)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
mode = sys.argv[3] if len(sys.argv) > 3 else 'nodigon'
from collections import Counter
cnt = Counter(); t0 = time.time()
for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 20):
    while True:
        mu, tau = rand_fpf(rng), rand_fpf(rng)
        if any(mu[b] == tau[b] == sig[b] for b in L): continue          # sibling pair: not a free star
        if mode == 'nodigon' and any(mu[b] == tau[b] for b in L): continue
        break
    r = rho_exists(mu, tau)
    cnt[r] += 1
print(mode, dict(cnt), '%.0fs' % (time.time() - t0))
