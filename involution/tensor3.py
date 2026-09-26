"""3-way consistency test of the local system at a block: T3[a,b,d] = #{far z containing a : b in Phi(z), d in Phi(sigma z)}
must be a nonnegative integer tensor with marginals
   sum_a T3 = 2 K3[b,d],  sum_d T3 = 2 K2[a,b],  sum_b T3 = 2 K2[pi a, d],
and bounds T3 <= min(K2[a,b], K2[pi a, d], K3[b,d]); plus T3[a,b,d] = 0 if b = d (Phi z, Phi sigma z disjoint at a full block
only) -- here we only use the general facts.  LP and ILP over all block types (mu up to B6, all tau), K >= 0."""
import itertools, sys, time
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
sys.path.insert(0, '.')
from local_gram_psd import reps_mu, ALL, P, comp, PI, I, J, Ppi

L = list(range(12))
idx = lambda a, b, d: (a * 12 + b) * 12 + d
rows_tpl = []
for b in L:
    for d in L:
        rows_tpl.append(('3', b, d, [idx(a, b, d) for a in L]))
for a in L:
    for b in L:
        rows_tpl.append(('2', a, b, [idx(a, b, d) for d in L]))
for a in L:
    for d in L:
        rows_tpl.append(('2p', a, d, [idx(a, b, d) for b in L]))
ri = [r for r, t in enumerate(rows_tpl) for _ in t[3]]; ci = [c for t in rows_tpl for c in t[3]]
A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(rows_tpl), 1728))
stats = {'types': 0, 'LP_inf': 0, 'ILP_inf': 0}
t0 = time.time()
limit = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 9
for mu in reps_mu:
    Pmu = P(mu)
    for tau in ALL:
        K2 = 2 * J - I - Ppi - Pmu - P(comp(tau, PI))
        K3 = 2 * J - Ppi - P(tau) - P(comp(tau, mu)) - P(comp(mu, tau))
        if (K2 < 0).any() or (K3 < 0).any():
            continue
        stats['types'] += 1
        rhs = []
        for kind, x, y, _ in rows_tpl:
            rhs.append(2 * K3[x, y] if kind == '3' else (2 * K2[x, y] if kind == '2' else 2 * K2[PI[x], y]))
        rhs = np.array(rhs, float)
        ub = np.array([min(K2[a, b], K2[PI[a], d], K3[b, d]) for a in L for b in L for d in L], float)
        lp = linprog(np.zeros(1728), A_eq=A, b_eq=rhs, bounds=list(zip(np.zeros(1728), ub)), method='highs')
        if lp.status == 2:
            stats['LP_inf'] += 1
            print('LP-infeasible type: mu', mu, 'tau', tau, flush=True)
        if stats['types'] >= limit:
            break
    if stats['types'] >= limit:
        break
print(stats, '%.0fs' % (time.time() - t0))
