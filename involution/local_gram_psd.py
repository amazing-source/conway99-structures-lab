"""Necessary PSD test of the local system (i)-(iii) at a block c, for EVERY block type (mu, tau) (C2 notation: mu = mates,
tau = pi o h; C1's (tau, eta)).  For the 60 far roots z let v_z = (1_z ; 1_{Phi(z)}) in R^24 (N(x)-labels ; star labels) and
R the fixed-point-free involution z -> sigma z.  Then M_pm = sum_z v_z (v_z +- v_{sigma z})^T is a Gram matrix (PSD), and by
(i)-(iii) it is KNOWN:
  M_pm = [[W W^T +- W R W^T, K2 +- P_pi K2], [(.)^T, Y Y^T +- K3]],
  W W^T = 9I + J - P_pi,  W R W^T = J - I + 9 P_pi,  Y Y^T = 9I + J - P_mu,  W Y^T = K2,  W R Y^T = P_pi K2,  Y R Y^T = K3.
If M_+ or M_- is not PSD, (i)-(iii) have no solution for that block type.  Also rank M_pm <= 30.
All (mu, tau) pairs of fixed-point-free involutions on the 12 labels, mu up to the hyperoctahedral group; K2, K3 >= 0 kept."""
import itertools
import numpy as np

L = list(range(12))
PI = [x ^ 1 for x in L]


def fpf(pts):
    if not pts:
        yield {}
        return
    a = pts[0]
    for i in range(1, len(pts)):
        b = pts[i]
        rest = pts[1:i] + pts[i + 1:]
        for m in fpf(rest):
            m = dict(m); m[a] = b; m[b] = a
            yield m


ALL = [tuple(m[x] for x in L) for m in fpf(L)]
G = []
for perm in itertools.permutations(range(6)):
    for flips in range(64):
        G.append(tuple(2 * perm[i // 2] + ((i % 2) ^ ((flips >> (i // 2)) & 1)) for i in range(12)))


def conj(t, g):
    ginv = [0] * 12
    for i, gi in enumerate(g):
        ginv[gi] = i
    return tuple(g[t[ginv[x]]] for x in L)


def P(f):
    M = np.zeros((12, 12))
    for i in L:
        M[i, f[i]] = 1
    return M


def comp(f, g):                    # (f o g)(x) = f(g(x))
    return tuple(f[g[x]] for x in L)


Ppi = P(PI)
I = np.eye(12); J = np.ones((12, 12))
reps_mu = []
seen = set()
for t in ALL:
    if t in seen:
        continue
    seen |= {conj(t, g) for g in G}
    reps_mu.append(t)
print('mu classes', len(reps_mu), flush=True)
WW = 9 * I + J - Ppi
WRW = J - I + 9 * Ppi
stats = {'pairs': 0, 'K>=0': 0, 'Mplus_notPSD': 0, 'Mminus_notPSD': 0, 'rank>30': 0}
worst = []
for mu in reps_mu:
    Pmu = P(mu)
    for tau in ALL:
        stats['pairs'] += 1
        K2 = 2 * J - I - Ppi - Pmu - P(comp(tau, PI))
        K3 = 2 * J - Ppi - P(tau) - P(comp(tau, mu)) - P(comp(mu, tau))
        if (K2 < 0).any() or (K3 < 0).any():
            continue
        stats['K>=0'] += 1
        YY = 9 * I + J - Pmu
        out = []
        for sgn, key in ((1, 'Mplus_notPSD'), (-1, 'Mminus_notPSD')):
            M = np.block([[WW + sgn * WRW, K2 + sgn * Ppi @ K2], [(K2 + sgn * Ppi @ K2).T, YY + sgn * K3]])
            ev = np.linalg.eigvalsh(M)
            if ev.min() < -1e-9:
                stats[key] += 1
                out.append((key, round(float(ev.min()), 4)))
            if (ev > 1e-9).sum() > 30:
                stats['rank>30'] += 1
        if out and len(worst) < 5:
            worst.append((mu, tau, out))
print(stats)
for w in worst:
    print(w)
