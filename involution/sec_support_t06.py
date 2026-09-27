"""Which constraints kill the 6-secant configuration? (structures 27/09) Phase-1 LP (min total artificial slack) on the
step-1' model of rel 25 with SEC and the 6 secant iota-pairs of a true iota fixed; the equality duals give a Farkas-type
certificate; we print its support by row family (iota matching / K3 / linking / block-pair / block-level (ii))."""
import json, collections
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import iota_design as D
cfg = json.load(open('t06_0.json')); mu, tau = cfg['mu'], cfg['tau']
d = np.load('s1p_t06_rel25.npz'); idx = json.load(open('s1p_t06_rel25_index.json'))
A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape'])); lo, ub = d['lo'], d['ub']
edges = [tuple(e) for e in idx['edges (non-mu pairs)']]; pairs = [tuple(p) for p in idx['iota_pairs']]
pv = {p: n for n, p in enumerate(pairs)}; npair = len(pairs); E = len(edges)
fam = ['iota-match'] * E + ['K3'] * 78 + ['link'] * npair + ['blockpair'] * 15 + ['(ii)blk'] * 72
assert len(fam) == A.shape[0]
K3rows = [(b, dd) for b in range(12) for dd in range(12) if dd >= b]
st, sol, _ = D.step1(mu, tau, 1, 60)
sec = [q for q in sol if set(edges[q[0]]) & set(edges[q[1]])]
lb = np.zeros(A.shape[1])
for (i, j) in sec:
    lb[pv[(min(i, j), max(i, j))]] = 1
m, nv = A.shape
Aph = sp.hstack([A, sp.identity(m), -sp.identity(m)]).tocsr()
c = np.concatenate([np.zeros(nv), np.ones(2 * m)])
bounds = list(zip(lb, ub)) + [(0, None)] * (2 * m)
r = linprog(c, A_eq=Aph, b_eq=lo, bounds=bounds, method='highs')
print('phase-1 optimum (total infeasibility):', round(r.fun, 6))
y = r.eqlin.marginals
sup = [i for i in range(m) if abs(y[i]) > 1e-9]
cnt = collections.Counter(fam[i] for i in sup)
print('certificate support: %d rows' % len(sup), dict(cnt))
for i in sup:
    if fam[i] == '(ii)blk':
        k = i - (E + 78 + npair + 15); print('   (ii) block %d label %d  y=%.3f  rhs=%g' % (k // 12, k % 12, y[i], lo[i]))
    elif fam[i] == 'K3':
        k = i - E; print('   K3 %s  y=%.3f  rhs=%g' % (K3rows[k], y[i], lo[i]))
    elif fam[i] == 'blockpair':
        k = i - (E + 78 + npair); print('   block pair %s  y=%.3f' % (idx['block_pairs'][k], y[i]))
print('secant pairs fixed:', [(edges[i], edges[j]) for (i, j) in sec])
