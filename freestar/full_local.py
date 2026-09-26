"""Free star, FULL compact formulation (i)+(ii)+(iii) (C2 note 20 par. 3) as one ILP, sigma-pair formulation.
Labels 0..11: block i//2, sign + if i even; sigma(i) = i ^ 1.  Far roots = pairs {a,a'} of labels from different blocks (60);
they come in 30 sigma-pairs Z = {z, sigma z}.  Unknown: for each Z (with a fixed representative z) the ORDERED pair of distinct
non-mu pairs (e, e') = (Phi(z), Phi(sigma z)).  Then
  (i)   every non-mu pair is used exactly once;
  (ii)  sum over far roots z containing a of [b in Phi(z)] = K2[a,b];
  (iii) sum over far roots z of [b in Phi(z)][d in Phi(sigma z)] = K3[b,d]  (both z and sigma z counted).
Pre-pruning (exact consequences of (ii) where K2 = 0): Phi(z) avoids every b with K2[a,b] = 0 for a in z."""
import itertools, json, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp

L = list(range(12))
sig = {i: i ^ 1 for i in L}


def K2(mu, tau):
    K = np.full((12, 12), 2) - np.eye(12, dtype=int)
    for a in L:
        K[a, sig[a]] -= 1
        K[a, mu[a]] -= 1
        K[a, tau[sig[a]]] -= 1
    return K


def K3(mu, tau):
    K = np.full((12, 12), 2)
    for b in L:
        for f in (sig, {x: tau[mu[x]] for x in L}, {x: mu[tau[x]] for x in L}, tau):
            K[b, f[b]] -= 1
    return K


def solve(mu, tau, tl=300, verbose=True):
    k2, k3 = K2(mu, tau), K3(mu, tau)
    assert (k2 >= 0).all() and (k3 >= 0).all()
    far = [e for e in itertools.combinations(L, 2) if sig[e[0]] != e[1]]
    reps = []
    seen = set()
    for z in far:
        if z in seen:
            continue
        sz = tuple(sorted((sig[z[0]], sig[z[1]])))
        seen |= {z, sz}
        reps.append((z, sz))
    nonmu = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    eidx = {e: i for i, e in enumerate(nonmu)}

    def allowed(z, e):
        return all(k2[a, b] > 0 for a in z for b in e)
    var = []
    for Zi, (z, sz) in enumerate(reps):
        for e in nonmu:
            if not allowed(z, e):
                continue
            for e2 in nonmu:
                if e2 != e and allowed(sz, e2):
                    var.append((Zi, e, e2))
    nv = len(var)
    if verbose:
        print('variables', nv, flush=True)
    ri, ci, vv, lo, hi = [], [], [], [], []
    r = 0

    def row(entries, val):
        nonlocal r
        for j, c in entries:
            ri.append(r); ci.append(j); vv.append(c)
        lo.append(val); hi.append(val); r += 1
    byZ = {}
    byE = {}
    for j, (Zi, e, e2) in enumerate(var):
        byZ.setdefault(Zi, []).append(j)
        byE.setdefault(e, []).append(j)
        byE.setdefault(e2, []).append(j)
    for Zi in range(len(reps)):
        row([(j, 1) for j in byZ.get(Zi, [])], 1)
    for e in nonmu:
        row([(j, 1) for j in byE.get(e, [])], 1)
    # (ii): for each a, b: sum_{z containing a} [b in Phi(z)]
    acc = {}
    for j, (Zi, e, e2) in enumerate(var):
        z, sz = reps[Zi]
        for a in z:
            for b in e:
                acc.setdefault((a, b), []).append(j)
        for a in sz:
            for b in e2:
                acc.setdefault((a, b), []).append(j)
    for a in L:
        for b in L:
            row([(j, 1) for j in acc.get((a, b), [])], int(k2[a, b]))
    # (iii): sum over far roots z of [b in Phi(z)][d in Phi(sigma z)]; each variable covers z -> (e, e2) and sigma z -> (e2, e)
    acc3 = {}
    for j, (Zi, e, e2) in enumerate(var):
        for b in e:
            for d in e2:
                acc3.setdefault((b, d), []).append(j)
        for b in e2:
            for d in e:
                acc3.setdefault((b, d), []).append(j)
    for b in L:
        for d in L:
            row([(j, 1) for j in acc3.get((b, d), [])], int(k3[b, d]))
    A = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r, nv))
    v_lo = np.array(lo, float)
    t0 = time.time()
    res = milp(np.zeros(nv), constraints=[LinearConstraint(A, v_lo, v_lo)], integrality=np.ones(nv), bounds=Bounds(0, 1),
               options={'time_limit': tl})
    st = {0: 'FEASIBLE', 2: 'INFEASIBLE'}.get(res.status, 'unknown(%d)' % res.status)
    if verbose:
        print('status', st, '%.0fs' % (time.time() - t0), flush=True)
    sol = None
    if res.status == 0:
        sol = {}
        for j in np.nonzero(res.x > 0.5)[0]:
            Zi, e, e2 = var[j]
            z, sz = reps[Zi]
            sol[str(list(z))] = list(e)
            sol[str(list(sz))] = list(e2)
    return st, sol


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'A111111_222'
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 300
    if which == 'A111111_222':
        mu = dict(sig)
        tau = {}
        for a, b in ((0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11)):
            tau[a] = b; tau[b] = a
    else:
        d = json.load(open(which))
        mu = {int(k): v for k, v in d['mu'].items()}
        tau = {int(k): v for k, v in d['tau'].items()}
    st, sol = solve(mu, tau, tl=tl)
    if sol:
        json.dump({'labels': 'i in 0..11: block i//2, sign + if i even; sigma(i) = i^1',
                   'mu': mu, 'tau': tau, 'Phi (far root as label pair -> non-mu pair)': sol,
                   'satisfies': '(i), (ii), (iii) of the compact formulation'}, open('full_local_%s.json' % which.replace('/', '_'), 'w'), indent=1)
        print('saved')
