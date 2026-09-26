"""Human-readable single-root contradictions for AD configurations.
For orbit W (support {a,b}) the pair-level facts are listed as equations in the unknowns phi(P); a deletion filter keeps a
minimal infeasible subset (the variable bounds / forced zeros are kept as 'facts')."""
import itertools, json, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from pair_level import N, PAIRS, PIDX, pr, cycle_matching, two_regular_graphs, parse_survivors

def facts_for_orbit(In, dd, ddm, Ds, o, A_known=True):
    """returns (eqs, lb, ub, names): eqs = list of (coef dict over pairs, value, label)."""
    anchor = {}
    for k in range(N):
        for m in In[k]: anchor[pr(k, m)] = k
    for p in dd: anchor[p] = None
    D = {}
    for k, cycles in Ds.items():
        for cyc in cycles:
            for (m, s), (m2, s2) in cycle_matching(cyc).items(): D[(pr(k, m), s)] = (pr(k, m2), s2)
    D.update(ddm)
    P, s = o; a, b = P
    lb = {Q: 0 for Q in PAIRS}; ub = {Q: 3 for Q in PAIRS}; why = {}
    eqs = []
    if anchor[P] is not None: ub[P] = 0; why[P] = 'own pair anchored => n_s = 0'
    else: ub[P] = 1
    for t in (a, b):
        eqs.append(({pr(t, m): 1 for m in range(N) if m != t}, 2, 'F1 at own block %d' % t))
        if anchor[P] == t:
            P2 = D[o][0]; m2 = P2[0] + P2[1] - t
            for m in range(N):
                if m == t: continue
                Q = pr(t, m)
                if m == m2: lb[Q] = ub[Q] = 2; why[Q] = 'D(W) in %d%d (pivot %d)' % (Q[0], Q[1], t)
                else: ub[Q] = 0; why[Q] = 'through-%d nbrs of w are in D(W)' % t
        else:
            for m in In[t]:
                Q = pr(t, m); ub[Q] = 0; why[Q] = 'W free at %d, %d%d anchored at %d' % (t, Q[0], Q[1], t)
    for k in range(N):
        if k in P: continue
        comps = [set(c) for c in Ds.get(k, [])]
        Fk = set(range(N)) - {k} - set(In[k])
        if Fk: comps.append(Fk)
        for S in comps:
            kind = 'F_%d' % k if S == Fk else 'M-cycle at %d' % k
            eqs.append(({pr(k, m): 1 for m in S}, len(S - {a, b}), 'balance at %d on %s %s (|S|=%d, own interiors %d)' % (k, kind, sorted(S), len(S), len(S & {a, b}))))
    if anchor[P] is None:
        Q = D[o][0]; lb[Q] = 2; why[Q] = 'dis partner D(W) in %d%d' % Q
    return eqs, lb, ub, why

def infeasible(eqs, lb, ub):
    idx = {Q: i for i, Q in enumerate(PAIRS)}
    if not eqs: return any(lb[Q] > ub[Q] for Q in PAIRS)
    A = np.zeros((len(eqs), 21)); v = np.zeros(len(eqs))
    for r, (c, val, _) in enumerate(eqs):
        for Q, x in c.items(): A[r, idx[Q]] = x
        v[r] = val
    res = milp(c=np.zeros(21), constraints=[LinearConstraint(A, v, v)], integrality=np.ones(21),
               bounds=Bounds([lb[Q] for Q in PAIRS], [ub[Q] for Q in PAIRS]))
    return res.status == 2

def minimal(eqs, lb, ub):
    keep = list(eqs)
    i = 0
    while i < len(keep):
        trial = keep[:i] + keep[i + 1:]
        if infeasible(trial, lb, ub): keep = trial
        else: i += 1
    return keep

def show(o, eqs, lb, ub, why):
    used = set()
    for c, _, _ in eqs: used |= set(c)
    print('  root w in orbit %d%d/%d' % (o[0][0], o[0][1], o[1]))
    for Q in sorted(used):
        if lb[Q] == ub[Q]: print('     phi(%d%d) = %d   [%s]' % (Q[0], Q[1], lb[Q], why.get(Q, '')))
        elif lb[Q] > 0: print('     phi(%d%d) >= %d   [%s]' % (Q[0], Q[1], lb[Q], why.get(Q, '')))
        elif ub[Q] < 3: print('     phi(%d%d) <= %d   [%s]' % (Q[0], Q[1], ub[Q], why.get(Q, '')))
    for c, val, lab in eqs:
        print('     %s = %d    <- %s' % (' + '.join('phi(%d%d)' % Q for Q in sorted(c)), val, lab))
