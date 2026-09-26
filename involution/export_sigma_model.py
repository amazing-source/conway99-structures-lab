"""Strengthened local model for C2's iota-DFS: sigma-pair variables y[Z, e, f] (Phi(z) = e, Phi(sigma z) = f) with
  (i) bijection, (ii) at LABEL level (orientation-aware), (iii) exact (linear in y), full-block extras (iv) (v),
  C1's pointwise conditions (L1) (S1) (S2) as variable filters (for z and sigma z),
plus explicit iota-pair variables p[n] (n = unordered pair {e, f}) with p[n] = sum_{Z, orientation} y  (branch on p),
plus C1's single-block lemma cuts on p (full block): no quadruple with 4 labels of Sp_p; at most one (cycle >= 3) / none
(digon) with 3 labels of Sp_p = {p+, p-, m(p+), m(p-)}.
Exports A (CSR), lo, hi, lb, ub, integrality (.npz) and an index (.json).  Infeasible => this block type is impossible."""
import itertools, json, sys
import numpy as np
import scipy.sparse as sp
from fullblock_rigid import mates_from_D, far_roots, K2K3, lab, pi, L

cfg = json.load(open(sys.argv[1])); out = sys.argv[2]
cycles = [tuple(c) for c in cfg['cycles']]
m = mates_from_D(cycles); mu = eta = m
K2, K3 = K2K3(m)
far, reps = far_roots()
nonm = [e for e in itertools.combinations(L, 2) if m[e[0]] != e[1]]


def pointwise_ok(z, e, f):
    """(L1), (S1), (S2) for the root z with Phi(z) = e, Phi(sigma z) = f."""
    se, sf = set(e), set(f)
    for a in z:
        for aa in (a, pi[a]):
            if aa in se and pi[aa] in sf:
                return False
        l = a
        if l in se and eta[l] in sf:
            return False
        if l not in se and mu[l] in se and eta[l] in sf:
            return False
    return True


def admissible(z, sz, e, f):
    if set(e) & set(f):
        return False
    if not (pointwise_ok(z, e, f) and pointwise_ok(sz, f, e)):
        return False
    for C in cycles:
        lc = lab(C)
        if len({a // 2 for a in z} & set(C)) + len(set(e) & lc) + len(set(f) & lc) != len(C):
            return False
    return True
var = []
for Zi, (z, sz) in enumerate(reps):
    for e in nonm:
        for f in nonm:
            if e != f and admissible(z, sz, e, f):
                var.append((Zi, e, f))
ny = len(var)
pk = sorted({tuple(sorted((e, f))) for _, e, f in var})
pidx = {q: ny + i for i, q in enumerate(pk)}
nv = ny + len(pk)
ri, ci, vv, lo, hi = [], [], [], [], []
r = [0]
def row(ent, a, b):
    for j, c in ent:
        ri.append(r[0]); ci.append(j); vv.append(c)
    lo.append(a); hi.append(b); r[0] += 1
byZ, byE, acc2, acc3, byP = {}, {}, {}, {}, {}
for j, (Zi, e, f) in enumerate(var):
    z, sz = reps[Zi]
    byZ.setdefault(Zi, []).append(j); byE.setdefault(e, []).append(j); byE.setdefault(f, []).append(j)
    byP.setdefault(tuple(sorted((e, f))), []).append(j)
    for a in z:
        for b in e: acc2.setdefault((a, b), []).append(j)
    for a in sz:
        for b in f: acc2.setdefault((a, b), []).append(j)
    for b in e:
        for d in f: acc3.setdefault((b, d), []).append(j)
    for b in f:
        for d in e: acc3.setdefault((b, d), []).append(j)
for Zi in range(30): row([(j, 1) for j in byZ.get(Zi, [])], 1, 1)
for e in nonm: row([(j, 1) for j in byE.get(e, [])], 1, 1)
for a in L:
    for b in L: row([(j, 1) for j in acc2.get((a, b), [])], int(K2[a, b]), int(K2[a, b]))
for b in L:
    for d in L: row([(j, 1) for j in acc3.get((b, d), [])], int(K3[b, d]), int(K3[b, d]))
for q in pk:                                             # p[n] = sum y
    row([(pidx[q], 1)] + [(j, -1) for j in byP[q]], 0, 0)
for C in cycles:                                         # C1's single-block lemma
    for p_ in C:
        Sp = {2 * p_, 2 * p_ + 1, m[2 * p_], m[2 * p_ + 1]}
        big = [pidx[q] for q in pk if len((set(q[0]) | set(q[1])) & Sp) == 3]
        four = [pidx[q] for q in pk if len((set(q[0]) | set(q[1])) & Sp) == 4]
        for j in four: row([(j, 1)], 0, 0)
        row([(j, 1) for j in big], 0, 1 if len(C) >= 3 else 0)
A = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r[0], nv))
np.savez_compressed(out + '.npz', A_data=A.data, A_indices=A.indices, A_indptr=A.indptr, A_shape=np.array(A.shape),
                    lo=np.array(lo, float), hi=np.array(hi, float), lb=np.zeros(nv), ub=np.ones(nv), integrality=np.ones(nv))
json.dump({'cfg': cfg, 'labels': 'v = 2i+s (block i, sign s)', 'mates m': m,
           'y variables 0..%d' % (ny - 1): 'y[j] = (sigma-pair Zi, Phi(z) = e, Phi(sigma z) = f)', 'y': [[Zi, list(e), list(f)] for Zi, e, f in var],
           'sigma-pairs': [[list(z), list(sz)] for z, sz in reps],
           'p variables %d..%d' % (ny, nv - 1): 'p[n] = 1 iff {e, f} is an iota-pair', 'p': [[list(q[0]), list(q[1])] for q in pk],
           'rows': '30 Z, 60 edges, 144 (ii) label level, 144 (iii), p-links, C1 lemma'}, open(out + '_index.json', 'w'))
print(out, 'y vars', ny, 'p vars', len(pk), 'rows', A.shape[0], 'nnz', A.nnz)
