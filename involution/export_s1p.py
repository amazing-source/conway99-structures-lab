"""Export the step-1' MILP (near-design + block-level (ii)) of a block type as .npz for C2's bb_certify.py / structures'
certify.py: A (CSR), lo, hi, lb, ub, integrality; plus a JSON index of the variables.
Infeasible => (i)+(ii)+(iii) infeasible for this block type => no block of this type in an involution of srg(99,14,1,2)."""
import itertools, json, sys
import numpy as np
import scipy.sparse as sp
import iota_design as D

L = list(range(12))
cfg = json.load(open(sys.argv[1])); out = sys.argv[2]
if 'cycles' in cfg:
    from fullblock_rigid import mates_from_D
    m = mates_from_D([tuple(c) for c in cfg['cycles']]); mu = tau = m
else:
    mu, tau = cfg['mu'], cfg['tau']
K2, K3 = D.K2K3(mu, tau)
edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
E = len(edges); BP = list(itertools.combinations(range(6), 2))
pairs = [(i, j) for i in range(E) for j in range(i + 1, E) if not any(K3[b, b] < 2 for b in set(edges[i]) & set(edges[j]))]
npair = len(pairs); pv = {pr: n for n, pr in enumerate(pairs)}
av = lambda n, k: npair + 15 * n + k
nv = npair * 16
ri, ci, vv, lo = [], [], [], []
CUT3 = []
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
        if d < b: continue
        t = int(K3[b, d]) if b != d else int(K3[b, b]) // 2
        row([(n, 1) for n in acc.get((b, d), [])], t)
for n in range(npair):
    row([(av(n, k), 1) for k in range(15)] + [(n, -1)], 0)
for k in range(15):
    row([(av(n, k), 1) for n in range(npair)], 2)
for rr in range(6):
    ks = [k for k, bp in enumerate(BP) if rr in bp]
    for b in L:
        # multiplicity [b in e_i] + [b in e_j]: a secant iota-pair (common label b) contributes 2 (C1's caution, 27/09)
        row([(av(n, k), (b in edges[i]) + (b in edges[j])) for n, (i, j) in enumerate(pairs)
             if b in edges[i] or b in edges[j] for k in ks], int(K2[2 * rr, b] + K2[2 * rr + 1, b]))
ncut = 0
if '--cuts' in sys.argv:          # generalised single-block lemma (NOTE_S29 par. 9): anchored left blocks p (p+ and p- pivots)
    for p_ in range(6):
        a_, b_ = 2 * p_, 2 * p_ + 1
        if not (mu[a_] == tau[a_] and mu[b_] == tau[b_]):
            continue
        Sp = {a_, b_, mu[a_], mu[b_]}
        digon = mu[a_] // 2 == mu[b_] // 2
        three = []
        for n, (i, j) in enumerate(pairs):
            k_ = len((set(edges[i]) | set(edges[j])) & Sp)
            if k_ == 4:
                row([(n, 1)], 0); ncut += 1
            elif k_ == 3:
                three.append(n)
        # at most one (none for a digon) iota-pair with 3 labels of Sp: row sum <= bound, written as equality with a slack
        # variable is not available here, so we use lo = 0, hi = bound through a separate list
        CUT3.append((three, 0 if digon else 1)); ncut += 1
A = sp.csr_matrix((np.array(vv, float), (np.array(ri), np.array(ci))), shape=(r[0], nv))
lo = np.array(lo, float)
hi = lo.copy()
if CUT3:                            # append the '<= bound' rows (lo = 0, hi = bound)
    extra = sp.csr_matrix((np.ones(sum(len(t) for t, _ in CUT3)),
                           ([k for k, (t, _) in enumerate(CUT3) for _ in t], [n for t, _ in CUT3 for n in t])),
                          shape=(len(CUT3), nv))
    A = sp.vstack([A, extra]).tocsr()
    lo = np.concatenate([lo, np.zeros(len(CUT3))]); hi = np.concatenate([hi, np.array([b for _, b in CUT3], float)])
np.savez_compressed(out + '.npz', A_data=A.data, A_indices=A.indices, A_indptr=A.indptr, A_shape=np.array(A.shape),
                    lo=lo, hi=hi, lb=np.zeros(nv), ub=np.ones(nv), integrality=np.ones(nv))
json.dump({'block_type': cfg, 'mu': list(mu), 'tau': list(tau), 'edges (non-mu pairs)': [list(e) for e in edges],
           'var n < %d' % npair: 'p[i,j] = 1 iff non-mu pairs i, j are iota-partners; pairs listed in iota_pairs',
           'iota_pairs': [list(p) for p in pairs],
           'var npair + 15 n + k': 'a[n,k] = iota-pair n is the image of a sigma-pair of block pair k',
           'block_pairs': [list(bp) for bp in BP],
           'rows': ['60 iota perfect matching', '78 K3 cross counts (diagonal: K3[b,b]/2)', 'npair: sum_k a[n,k] = p[n]',
                    '15: each block pair gets 2 iota-pairs', '72 block-level (ii): (block r, label b)'],
           'branching hint': 'branch on p (partner of the lowest unpaired non-mu pair); LP dies at ~12 fixed pairs along '
                             'a true near-design (0.2 s per LP); symmetry = Aut(pi, mu, tau) on labels'},
          open(out + '_index.json', 'w'), indent=1)
print(out, 'vars', nv, 'rows', A.shape[0], 'nnz', A.nnz, '| lemma cut rows', ncut)
