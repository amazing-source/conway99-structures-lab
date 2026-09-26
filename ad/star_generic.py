"""Generic (*) augmentation: add C2's generalized balance (*) (see star_lemma.py) to any pair-level model given as
(A, lo, hi, lb, ub, nvar) with phi[(orbit_key, pair)] -> variable index and orbit_pair[orbit_key] -> its pair."""
import itertools
import numpy as np
import scipy.sparse as sp
N = 7
def pr(a, b): return (a, b) if a < b else (b, a)
def augment(A, lo, hi, lb, ub, nvar, phi, orbit_pair, blocks=range(7)):
    A = sp.coo_matrix(A)
    ri, ci, vv = list(A.row), list(A.col), list(A.data)
    lo, hi, lb, ub = list(lo), list(hi), list(lb), list(ub)
    r = A.shape[0]; nv = nvar
    def newvar(l, u):
        nonlocal nv
        lb.append(l); ub.append(u); nv += 1; return nv - 1
    def addrow(coefs, l, h):
        nonlocal r
        for j, c in coefs: ri.append(r); ci.append(j); vv.append(c)
        lo.append(l); hi.append(h); r += 1
    orbs = list(orbit_pair.keys())
    for k in blocks:
        others = [m for m in range(N) if m != k]
        for s in range(1, 6):
            for S in itertools.combinations(others, s):
                Sset = set(S)
                lab = [o for o in orbs if k in orbit_pair[o] and (sum(orbit_pair[o]) - k) in Sset]
                labset = set(lab)
                terms = []
                for o in orbs:
                    P = orbit_pair[o]
                    if k in P: continue
                    base = len(set(P) & Sset) - s
                    lin = [(phi[(o, pr(k, m))], 1) for m in S]
                    q = newvar(0, np.inf)
                    for a, b in ((1, 0), (3, 2), (5, 6)):
                        addrow([(q, 1)] + [(j, -a * c) for j, c in lin], a * base - b, np.inf)
                        addrow([(q, 1)] + [(j, a * c) for j, c in lin], -a * base - b, np.inf)
                    terms.append(q)
                for o in orbs:
                    P = orbit_pair[o]
                    if k not in P or o in labset: continue
                    lin = [(phi[(o, pr(k, m))], 1) for m in S]
                    p = newvar(0, np.inf)
                    addrow([(p, 1)] + [(j, -c) for j, c in lin], 0, np.inf)
                    addrow([(p, 1)] + [(j, -3 * c) for j, c in lin], -2, np.inf)
                    terms.append(p)
                rv_list = []
                for o in lab:
                    lin = [(phi[(o, pr(k, m))], 1) for m in S]
                    rv = newvar(-np.inf, np.inf)
                    addrow([(rv, 1)] + [(j, 4 * c) for j, c in lin], -np.inf, 10)
                    addrow([(rv, 1)] + [(j, 6 * c) for j, c in lin], -np.inf, 12)
                    rv_list.append(rv)
                addrow([(t, 1) for t in terms] + [(t, -1) for t in rv_list], -np.inf, 0)
    A2 = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r, nv))
    return A2, np.array(lo, float), np.array(hi, float), np.array(lb, float), np.array(ub, float), nv
