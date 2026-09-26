"""Single-block covering filter (C1's idea) on saved near-designs: for block p, is there w in [0,1]^30 with sum w = 10 and
sum_Q w_Q (mult of b in Q) = t_p(b) = K2[p+,b] + K2[p-,b] for every label b?  (w_Q = [quadruple Q assigned to a block pair
through p].)  Also the refined version with the block-pair structure: w_Q = sum_{q != p} a[Q,{p,q}], each {p,q} gets exactly 2."""
import json, sys, itertools
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
from fullblock_rigid import mates_from_D
import iota_design as D

cfg = json.load(open(sys.argv[1])); iotas = json.load(open(sys.argv[1].replace('.json', '_iotas.json')))
m = mates_from_D([tuple(c) for c in cfg['cycles']]); K2, K3 = D.K2K3(m, m)
for k, pairs in enumerate(iotas):
    mult = [[(b in e) + (b in f) for b in range(12)] for e, f in pairs]
    out = []
    for p in range(6):
        A = [[1] * 30] + [[mult[j][b] for j in range(30)] for b in range(12)]
        bb = [10] + [int(K2[2 * p, b] + K2[2 * p + 1, b]) for b in range(12)]
        lp = linprog(np.zeros(30), A_eq=np.array(A, float), b_eq=np.array(bb, float), bounds=(0, 1), method='highs')
        # refined: a[j, q] for q != p, sum_q a = w_j <= 1, sum_j a[j,q] = 2
        nv = 30 * 5
        qs = [q for q in range(6) if q != p]
        rows, rhs = [], []
        for b in range(12):
            rows.append([mult[j][b] for j in range(30) for q in qs]); rhs.append(int(K2[2 * p, b] + K2[2 * p + 1, b]))
        for qi in range(5):
            rows.append([1 if qq == qi else 0 for j in range(30) for qq in range(5)]); rhs.append(2)
        Aub = [[1 if jj == j else 0 for jj in range(30) for qq in range(5)] for j in range(30)]
        lp2 = linprog(np.zeros(nv), A_eq=np.array(rows, float), b_eq=np.array(rhs, float), A_ub=np.array(Aub, float),
                      b_ub=np.ones(30), bounds=(0, 1), method='highs')
        out.append(('%d:%s/%s' % (p, 'ok' if lp.status == 0 else 'DEAD', 'ok' if lp2.status == 0 else 'DEAD')))
    print('iota #%d single-block (plain/refined):' % k, ' '.join(out))
