"""Step 1 with C1's full-block lemma (a3(p) + 3 a4(p) <= 1 on cycles of length >= 3, a3 = a4 = 0 on digons), then the
single-block filter and step 2 on each iota found."""
import json, sys, time
import numpy as np
from scipy.optimize import linprog
import iota_design as D
from fullblock_rigid import mates_from_D

cfg = json.load(open(sys.argv[1])); n = int(sys.argv[2]); tl = float(sys.argv[3])
cycles = [tuple(c) for c in cfg['cycles']]
m = mates_from_D(cycles); mu = tau = m
K2, K3 = D.K2K3(mu, tau)
lemma = []
for C in cycles:
    for p in C:
        Sp = {2 * p, 2 * p + 1, m[2 * p], m[2 * p + 1]}
        lemma.append((Sp, 1 if len(C) >= 3 else 0))
found = []; t0 = time.time()
for k in range(n):
    st, sol, edges = D.step1(mu, tau, 1, tl, forbid=found, c1lemma=lemma)
    print('step1+lemma #%d' % k, st, '%.1fs' % (time.time() - t0), flush=True)
    if sol is None: break
    found.append(sol)
    pairs = [(edges[i], edges[j]) for i, j in sol]
    mult = [[(b in e) + (b in f) for b in range(12)] for e, f in pairs]
    dead = []
    for p in range(6):
        A = [[1] * 30] + [[mult[j][b] for j in range(30)] for b in range(12)]
        bb = [10] + [int(K2[2 * p, b] + K2[2 * p + 1, b]) for b in range(12)]
        if linprog(np.zeros(30), A_eq=np.array(A, float), b_eq=np.array(bb, float), bounds=(0, 1), method='highs').status == 2:
            dead.append(p)
    st2, _ = D.step2(mu, tau, sol, edges, tl=30)
    print('   single-block dead at', dead, '| step 2:', st2, flush=True)
