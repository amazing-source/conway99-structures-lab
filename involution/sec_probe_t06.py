"""C1's probes on the (0,6) step-1' LP with SEC (structures 27/09), rel 25 = t06_0, along true near-designs from step 1:
  (i)   fix only the 6 secant iota-pairs of the true iota;
  (ii)  (i) + the crossed quadruples (p* and q* = mu p* on the two halves; one per transition, so 3);
  (iii) 8 random NON-secant iota-pairs;   (iv) 8 random iota-pairs that include every secant.
Reports F (LP feasible) or X (infeasible)."""
import json, random
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import iota_design as D
cfg = json.load(open('t06_0.json')); mu, tau = cfg['mu'], cfg['tau']
piv = [l for l in range(12) if mu[l] == tau[l]]
d = np.load('s1p_t06_rel25.npz'); idx = json.load(open('s1p_t06_rel25_index.json'))
A = sp.csr_matrix((d['A_data'], d['A_indices'], d['A_indptr']), shape=tuple(d['A_shape'])); lo, ub = d['lo'], d['ub']
edges = [tuple(e) for e in idx['edges (non-mu pairs)']]
pv = {tuple(p): n for n, p in enumerate(idx['iota_pairs'])}
def lp(fixed):
    lb = np.zeros(A.shape[1])
    for (i, j) in fixed:
        lb[pv[(min(i, j), max(i, j))]] = 1
    if (lb > ub).any():
        return 'X(fix)'
    r = linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=lo, bounds=list(zip(lb, ub)), method='highs')
    return 'F' if r.status == 0 else 'X' if r.status == 2 else str(r.status)
found, rng = [], random.Random(2)
for it in range(2):
    st, sol, E2 = D.step1(mu, tau, 1, 60, forbid=found)
    if sol is None:
        print('step 1', st); break
    found.append(sol)
    sec = [(i, j) for (i, j) in sol if set(edges[i]) & set(edges[j])]
    crossed = []
    for ps in piv:
        qs = mu[ps]
        for (i, j) in sol:
            e, f = set(edges[i]), set(edges[j])
            if ((ps in e and qs in f) or (ps in f and qs in e)) and (i, j) not in crossed:
                crossed.append((i, j))
    nonsec = [q for q in sol if q not in sec]
    print('iota #%d: %d secant pairs, %d crossed pairs (overlap with secants: %d)' %
          (it, len(sec), len(crossed), len(set(sec) & set(crossed))), flush=True)
    print('  (i)   6 secants            :', lp(sec), flush=True)
    print('  (ii)  secants + crossed    :', lp(sec + [q for q in crossed if q not in sec]), flush=True)
    out = []
    for t in range(4):
        out.append(lp(rng.sample(nonsec, 8)))
    print('  (iii) 8 random non-secant  :', ' '.join(out), flush=True)
    out = []
    for t in range(3):
        out.append(lp(sec + rng.sample(nonsec, 2)))
    print('  (iv)  6 secants + 2 random :', ' '.join(out), flush=True)
