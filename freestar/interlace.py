"""Interlacing (real PSD) test on the free star: 3I - S is PSD (Gram of the odd lattice N, eigenvalues 0 and 7), so for the
19 block-line coordinates I (d_c, 6 other blocks d_e, 12 orbits o_v through c) and any far orbit f:
    q(f) = n_f^T G^+ n_f <= 3,   G = (3I - S)_{II},  n_f = -S_{I,f}.
S on I is fixed by (mu, tau) (C2's notation: mu = line partners, tau = pi partners);  the far column is fixed by the far root
z = {a, a'} (signs at its two blocks) and by Phi(z), Phi(sigma z) (non-mu pairs of labels): S_{f,o_v} = [v in Phi(z)] - [v in Phi(sz)].
Question: how restrictive is q(f) <= 3 on the possible (Phi(z), Phi(sz)) for a far root z?"""
import itertools, random, sys
import numpy as np
L = [(e, s) for e in range(6) for s in (0, 1)]          # labels: block e (0..5 = the 6 other blocks), sign s
LI = {l: i for i, l in enumerate(L)}
neg = {l: (l[0], 1 - l[1]) for l in L}
def S_ball(mu, tau):
    n = 19; S = np.zeros((n, n))
    S[0, 0] = -1
    for e in range(6): S[1 + e, 1 + e] = -1
    for v in L:
        o = 7 + LI[v]
        S[0, o] = S[o, 0] = 1
        sg = 1 if v[1] == 0 else -1
        S[1 + v[0], o] = S[o, 1 + v[0]] = sg
        S[o, 7 + LI[mu[v]]] += 1
        S[o, 7 + LI[tau[v]]] -= 1
    return S
def far_col(z, phiz, phisz):
    """z = (a, a') labels in different blocks; phiz, phisz = 2-sets of labels"""
    col = np.zeros(19)
    for a in z:
        col[1 + a[0]] = 1 if a[1] == 0 else -1
    for v in phiz: col[7 + LI[v]] += 1
    for v in phisz: col[7 + LI[v]] -= 1
    return col
def rand_fpf(rng):
    xs = L[:]; rng.shuffle(xs); f = {}
    for i in range(0, 12, 2): f[xs[i]], f[xs[i + 1]] = xs[i + 1], xs[i]
    return f
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3):
    while True:
        mu, tau = rand_fpf(rng), rand_fpf(rng)
        if any(mu[b] == tau[b] for b in L): continue          # fully free star
        break
    S = S_ball(mu, tau)
    G = 3 * np.eye(19) - S
    Gp = np.linalg.pinv(G)
    rk = np.linalg.matrix_rank(G)
    nonmu = [p for p in itertools.combinations(L, 2) if mu[p[0]] != p[1]]
    tot = ok = 0; qs = []
    z = (L[0], L[2])          # a far root on blocks 0 and 1
    for p1 in nonmu:
        for p2 in nonmu:
            if p1 == p2: continue
            col = far_col(z, p1, p2)
            n = -col
            q = n @ Gp @ n
            tot += 1; ok += q <= 3 + 1e-9; qs.append(q)
    print('rank G =', rk, '| far root z on blocks 0,1: allowed (Phi(z),Phi(sz)) pairs', ok, 'of', tot,
          '| q quantiles', np.round(np.quantile(qs, [0, 0.1, 0.5, 0.9, 1]), 3))
