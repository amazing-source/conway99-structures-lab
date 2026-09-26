"""Pair-level model with SIBLING pairs (Sigma != empty), for C2's surviving motifs.  Written by C1 at the request of
"structures"; imports general_model without modifying it.

Skeleton format = general_model's (42 orbits {pair, pivot, partner}), plus pivot 'sib': the two orbits of a sibling pair are
each other's D-partner (D(U) = sister, n_s = 2).  Everything imposed is proved:
  (S1)  at each end k of a sib pair {k,m}, {m} is a closed component of G_k (prism) -> balance gives C2's 8.1;
  (S2)  a sib orbit has phi(own pair) = 2 and zeros on every other pair through either of its blocks (F1 saturated);
  (S3)  an orbit W not in the sib pair {t,m} but through t has phi_W(tm) = 0 (its through-t neighbours would be D-partners
        of the sib roots, i.e. in the sister orbit) -- also implied by (S1);
  (dich) for a sib pair P the dichotomy constant is 2 - n_s(P) = 0: (Q^2)[U,U'] = 4 - Sh[U,U'] - n_s = 2 - 2 = 0, so every
        orbit Z outside P has phi_Z(P) <= 1 (and never 3);
  all general_model facts (F1, pivot zeros, free caps, balance on certain unions of components, dis, symmetry, dichotomy)
  for the non-sib orbits.
"""
import itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from general_model import N, PAIRS, pr, closed_classes


def sib_pairs(orbs):
    return {o['pair'] for o in orbs if o['pivot'] == 'sib'}


def check_skeleton_sib(orbs):
    assert len(orbs) == 42
    cnt = {}
    for i, o in enumerate(orbs):
        cnt[o['pair']] = cnt.get(o['pair'], 0) + 1
        j = o['partner']
        if o['pivot'] == 'sib':
            assert j is not None and orbs[j]['pivot'] == 'sib' and orbs[j]['pair'] == o['pair'] and j != i
            continue
        if j is None:
            assert o['pivot'] is not None
            continue
        assert orbs[j]['partner'] == i and j != i
        P, Q = set(o['pair']), set(orbs[j]['pair'])
        if o['pivot'] is None:
            assert orbs[j]['pivot'] is None and not (P & Q)
        else:
            assert orbs[j]['pivot'] == o['pivot'] and P & Q == {o['pivot']}
    assert all(cnt.get(P) == 2 for P in PAIRS)


def components_sib(orbs, k, sibs):
    """certain unions of components of G_k: closed D-classes, sib prisms {m}, and the complement."""
    closed, rest = closed_classes(orbs, k)
    prisms = [{m} for m in range(N) if m != k and pr(k, m) in sibs]
    pm = set().union(*prisms) if prisms else set()
    closed = [c for c in closed if not (c & pm)]
    rest = set(rest) - pm
    return closed + prisms, rest


def build_sib(orbs, dichotomy=True, balance=True):
    check_skeleton_sib(orbs)
    sibs = sib_pairs(orbs)
    n_orb = len(orbs)
    by_pair = {P: [i for i, o in enumerate(orbs) if o['pair'] == P] for P in PAIRS}
    both_dis = {P: all(orbs[i]['pivot'] is None for i in by_pair[P]) for P in PAIRS}
    cnt = [0]

    def new(n=1):
        a = cnt[0]
        cnt[0] += n
        return a
    phi = {(i, P): new() for i in range(n_orb) for P in PAIRS}
    ns = {P: new() for P in PAIRS if both_dis[P]}
    y = {}
    if dichotomy:
        for i in range(n_orb):
            for P in PAIRS:
                if P != orbs[i]['pair']:
                    y[(i, P)] = new(4)
    xv = {}
    freeorb = [i for i, o in enumerate(orbs) if o['partner'] is None]
    for i, j in itertools.combinations(freeorb, 2):
        if orbs[i]['pivot'] == orbs[j]['pivot'] and orbs[i]['pair'] != orbs[j]['pair']:
            xv[(i, j)] = new()
    nvar = cnt[0]
    lb = np.zeros(nvar)
    ub = np.full(nvar, 3.0)
    for key in xv.values():
        ub[key] = 1
    rows, lo, hi = [], [], []

    def add(coefs, l, h):
        rows.append(coefs)
        lo.append(l)
        hi.append(h)
    comps = {k: components_sib(orbs, k, sibs) for k in range(N)}
    for i, o in enumerate(orbs):
        P = o['pair']
        a, b = P
        if o['pivot'] == 'sib':
            lb[phi[(i, P)]] = ub[phi[(i, P)]] = 2
        elif P in ns:
            add([(phi[(i, P)], 1), (ns[P], -1)], 0, 0)
        else:
            ub[phi[(i, P)]] = 0
        for t in (a, b):
            add([(phi[(i, pr(t, m))], 1) for m in range(N) if m != t], 2, 2)
            if o['pivot'] == 'sib':
                for m in range(N):
                    if m != t and pr(t, m) != P:
                        ub[phi[(i, pr(t, m))]] = 0
                continue
            for m in range(N):                       # (S3)
                if m != t and pr(t, m) in sibs and pr(t, m) != P:
                    ub[phi[(i, pr(t, m))]] = 0
            if o['pivot'] == t and o['partner'] is None:
                xs = {j: xv[pr(i, j)] for j in freeorb if pr(i, j) in xv}
                add([(v, 1) for v in xs.values()], 1, 1)
                for m in range(N):
                    if m == t:
                        continue
                    Q = pr(t, m)
                    js = [v for j, v in xs.items() if orbs[j]['pair'] == Q]
                    if js:
                        add([(phi[(i, Q)], 1)] + [(v, -2) for v in js], 0, 0)
                    else:
                        ub[phi[(i, Q)]] = 0
            elif o['pivot'] == t:
                q = orbs[o['partner']]['pair']
                for m in range(N):
                    if m == t:
                        continue
                    j = phi[(i, pr(t, m))]
                    if pr(t, m) == q:
                        lb[j] = ub[j] = 2
                    else:
                        ub[j] = 0
            else:
                for m in range(N):
                    if m == t or pr(t, m) == P or pr(t, m) in sibs:
                        continue
                    cap = sum(1 for j2 in by_pair[pr(t, m)] if orbs[j2]['pivot'] != t)
                    ub[phi[(i, pr(t, m))]] = min(ub[phi[(i, pr(t, m))]], cap)
        if balance:
            for k in range(N):
                if k in P:
                    continue
                closed, rest = comps[k]
                for S in closed + ([rest] if rest else []):
                    add([(phi[(i, pr(k, m))], 1) for m in S], len(S - {a, b}), len(S - {a, b}))
        if o['pivot'] is None:
            j = phi[(i, orbs[o['partner']]['pair'])]
            lb[j] = max(lb[j], 2)
    for P, P2 in itertools.combinations(PAIRS, 2):
        add([(phi[(i, P2)], 1) for i in by_pair[P]] + [(phi[(i, P)], -1) for i in by_pair[P2]], 0, 0)
    for P in ns:
        ub[ns[P]] = 1
    if dichotomy:
        for (i, P), yy in y.items():
            for f in range(4):
                ub[yy + f] = 1
            add([(yy + f, 1) for f in range(4)], 1, 1)
            add([(phi[(i, P)], 1)] + [(yy + f, -f) for f in range(4)], 0, 0)
            if orbs[i]['partner'] is None:
                ub[yy + 3] = 0
                continue
            is_partner = orbs[orbs[i]['partner']]['pair'] == P
            if not (is_partner and orbs[i]['pivot'] is None):
                ub[yy + 3] = 0
        for P in PAIRS:
            coefs = []
            partners = {orbs[j]['partner'] for j in by_pair[P] if orbs[j]['partner'] is not None}
            for i in range(n_orb):
                if orbs[i]['pair'] == P:
                    continue
                yy = y[(i, P)]
                if orbs[i]['partner'] is None:
                    coefs.append((yy + 2, 1))
                    for j in by_pair[P]:
                        if pr(i, j) in xv:
                            coefs.append((xv[pr(i, j)], -1))
                    continue
                if i not in partners:
                    coefs.append((yy + 2, 1))
                coefs.append((yy + 3, 2))
            if P in sibs:
                add(coefs, 0, 0)
                continue
            if P in ns:
                coefs.append((ns[P], 1))
            add(coefs, 2, 2)
    A = np.zeros((len(rows), nvar))
    for r, cs in enumerate(rows):
        for j, c in cs:
            A[r, j] += c
    return A, np.array(lo, float), np.array(hi, float), lb, ub, nvar, phi, ns


def feasible_sib(orbs, dichotomy=True, balance=True, time_limit=120):
    A, lo, hi, lb, ub, nvar, phi, ns = build_sib(orbs, dichotomy, balance)
    if np.any(lb > ub):
        return None
    res = milp(c=np.zeros(nvar), constraints=[LinearConstraint(A, lo, hi)], integrality=np.ones(nvar),
               bounds=Bounds(lb, ub), options={'time_limit': time_limit})
    if res.status == 0:
        return res.x, phi, ns
    if res.status == 2:
        return None
    return 'unknown(%d)' % res.status


def from_anchored_motif(sig, In, Dg):
    """Anchored motif skeleton (C1, notes/130): sig = sib pairs (0-based), In[k] = blocks m with pair km anchored at k,
    Dg[k] = D-graph at pivot k as a list of edges (m, m') with multiplicity (2-regular loopless on In[k]).
    Orbit s = 0/1 of pair (k,m) takes the s-th D-edge of m at k (WLOG: swapping the two orbits of a pair is a symmetry)."""
    idx = {}
    orbs = []
    for P in PAIRS:
        for s in (0, 1):
            idx[(P, s)] = len(orbs)
            orbs.append({'pair': P, 'pivot': None, 'partner': None})
    for P in sig:
        P = pr(*P)
        i0, i1 = idx[(P, 0)], idx[(P, 1)]
        orbs[i0].update(pivot='sib', partner=i1)
        orbs[i1].update(pivot='sib', partner=i0)
    for k, E in Dg.items():
        used = {}
        for (m, m2) in E:
            s1 = used.get(m, 0); used[m] = s1 + 1
            s2 = used.get(m2, 0); used[m2] = s2 + 1
            i, j = idx[(pr(k, m), s1)], idx[(pr(k, m2), s2)]
            orbs[i].update(pivot=k, partner=j)
            orbs[j].update(pivot=k, partner=i)
    return orbs
