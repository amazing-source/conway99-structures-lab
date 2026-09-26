"""C2's generalized balance (*) added to the pair-level model (general_model.build), as valid linear inequalities.

(*) (C2, 26/09; re-derived here): for a block k, a set S of blocks (k not in S, s = |S|), Lab = orbits of the pairs {k,m},
m in S, and g_o = sum_{m in S} phi_o(km) for o in Lab (own pair included):
    sum_{W avoiding k} delta_W^2 + sum_{W through k, not in Lab} eps_W^2 = sum_{o in Lab} (2 - g_o)(g_o + 5),
    delta_W = |P_W n S| + sum_{m in S} phi_W(km) - s,   eps_W = sum_{m in S} phi_W(km).
Proof (per root, then /2): for C = C_S in N(l_k), with B the incidence N(l_k) x far roots,
    BB^T = 12I + 2J - A - Cls - A^2   (lambda = 1, mu = 2, one l-neighbour per vertex, G(l_k) 3-regular),
so sum_z (n_C(z) - s)^2 = sum_{v in C} L_v (7 - L_v) - sum_{v not in C} d_C(v)^2 with L_v = #G(l)-neighbours of v outside C;
interiors of S have L = 0, a star root of an orbit o in Lab has L = 2 - g_o, a star root outside C has d_C = eps.
Linear use (exact at integer points): q_W >= tangent lines of delta^2 at integers, p_W >= tangents of eps^2,
r_o <= min(10 - 4g, 12 - 6g) (the concave (2-g)(g+5) on {0,1,2}), and sum q + sum p <= sum r."""
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.optimize import milp, LinearConstraint, Bounds
from general_model import build, N, PAIRS, pr


def build_star(orbs, subsets='all', **kw):
    A, lo, hi, lb, ub, nvar, phi, ns = build(orbs, **kw)
    A = A.tocoo()
    ri, ci, vv = list(A.row), list(A.col), list(A.data)
    lo, hi = list(lo), list(hi)
    lb, ub = list(lb), list(ub)
    r = A.shape[0]
    nv = nvar
    n_orb = len(orbs)

    def newvar(l=-np.inf, u=np.inf):
        nonlocal nv
        lb.append(l); ub.append(u); nv += 1
        return nv - 1

    def addrow(coefs, l, h):
        nonlocal r
        for j, c in coefs:
            ri.append(r); ci.append(j); vv.append(c)
        lo.append(l); hi.append(h); r += 1
    for k in range(N):
        others = [m for m in range(N) if m != k]
        for s in range(1, 6):
            for S in itertools.combinations(others, s):
                Sset = set(S)
                lab = [i for i, o in enumerate(orbs) if k in o['pair'] and (o['pair'][0] + o['pair'][1] - k) in Sset]
                terms = []
                for i, o in enumerate(orbs):
                    P = o['pair']
                    if k in P:
                        continue
                    # delta = |P n S| + sum_{m in S} phi(km) - s
                    base = len(set(P) & Sset) - s
                    lin = [(phi[(i, pr(k, m))], 1) for m in S]
                    q = newvar(0, np.inf)
                    for a, b in ((1, 0), (3, 2), (5, 6)):
                        # q >= a*delta - b  and  q >= -a*delta - b
                        addrow([(q, 1)] + [(j, -a * c) for j, c in lin], a * base - b, np.inf)
                        addrow([(q, 1)] + [(j, a * c) for j, c in lin], -a * base - b, np.inf)
                    terms.append(q)
                for i, o in enumerate(orbs):
                    P = o['pair']
                    if k not in P or i in lab:
                        continue
                    lin = [(phi[(i, pr(k, m))], 1) for m in S]
                    p = newvar(0, np.inf)
                    addrow([(p, 1)] + [(j, -c) for j, c in lin], 0, np.inf)
                    addrow([(p, 1)] + [(j, -3 * c) for j, c in lin], -2, np.inf)
                    terms.append(p)
                rvars = []
                for i in lab:
                    lin = [(phi[(i, pr(k, m))], 1) for m in S]
                    rv = newvar(-np.inf, np.inf)
                    addrow([(rv, 1)] + [(j, 4 * c) for j, c in lin], -np.inf, 10)
                    addrow([(rv, 1)] + [(j, 6 * c) for j, c in lin], -np.inf, 12)
                    rvars.append(rv)
                addrow([(t, 1) for t in terms] + [(t, -1) for t in rvars], -np.inf, 0)
    A2 = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r, nv))
    integ = np.concatenate([np.ones(nvar), np.zeros(nv - nvar)])
    return A2, np.array(lo, float), np.array(hi, float), np.array(lb, float), np.array(ub, float), nv, integ, phi, ns


def feasible_star(orbs, time_limit=120, lp=False, **kw):
    A, lo, hi, lb, ub, nv, integ, phi, ns = build_star(orbs, **kw)
    if np.any(lb > ub):
        return None
    res = milp(np.zeros(nv), constraints=[LinearConstraint(A, lo, hi)], integrality=(np.zeros(nv) if lp else integ),
               bounds=Bounds(lb, ub), options={'time_limit': time_limit})
    if res.status == 0:
        return 'feasible'
    if res.status == 2:
        return None
    return 'unknown(%d)' % res.status
