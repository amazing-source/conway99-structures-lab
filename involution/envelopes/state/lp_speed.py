import json, sys, time, os
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, '.')
from build156 import rows_frac
import frame_lp as FL
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, identity
rec = json.load(open('U/132.json'))
P = FL.prepare(rows_frac(rec['gens']), verbose=False)
data = json.load(open('state/pairs_132.json'))
G = P['G']
for k in (1210, 2020, 410):
    y1, y2 = FL.coef_of(P, [F(x) for x in data['pairs'][k][0]]), FL.coef_of(P, [F(x) for x in data['pairs'][k][1]])
    g1, g2 = G @ y1, G @ y2
    N3, N4, CH = P['N3'], P['N4'], P['CH']
    Us = N3[(N3 @ g1 == 0) & (N3 @ g2 == 0)]
    H0 = N4[(N4 @ g1 == 0) & (N4 @ g2 == 0)]
    Hs = H0[(np.abs(CH @ G @ H0.T) <= 2).all(axis=0)]
    # star filter: at least 12 u-candidates with |u.h| = 1
    UH = np.abs(Us @ G @ Hs.T)
    star = (UH == 1).sum(axis=0)
    Hf = Hs[star >= 12]
    for label, HH in (('all h', Hs), ('star>=12', Hf)):
        iu = [(i, j) for i in range(24) for j in range(i, 24)]
        Ue, He = Us @ G, HH @ G
        nu, nh = len(Us), len(HH)
        Aeq = np.vstack([np.hstack([np.array([Ue[:, i] * Ue[:, j] for i, j in iu]), np.array([He[:, i] * He[:, j] for i, j in iu])]),
                         np.concatenate([np.zeros(nu), np.ones(nh)])[None, :], np.concatenate([np.ones(nu), np.zeros(nh)])[None, :]]).astype(float)
        beq = np.array([7 * G[i, j] - g1[i] * g1[j] - g2[i] * g2[j] for i, j in iu] + [7, 42], dtype=float)
        n = nu + nh; m = len(beq)
        A = hstack([csr_matrix(Aeq), identity(m, format='csr'), -identity(m, format='csr')]).tocsr()
        for meth in ('highs-ds', 'highs-ipm'):
            t0 = time.time()
            lp = linprog(np.concatenate([np.zeros(n), np.ones(2 * m)]), A_eq=A, b_eq=beq, bounds=[(0, 1)] * n + [(0, None)] * (2 * m), method=meth)
            print('pair %d %-9s %-9s n_u %d n_h %d  phase1 %.4f  %.1fs' % (k, label, meth, nu, nh, lp.fun if lp.status == 0 else -1, time.time() - t0), flush=True)
