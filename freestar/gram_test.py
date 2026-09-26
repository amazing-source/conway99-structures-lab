"""Free-star compact formulation (C2 note 20 par. 3): unknown Phi, i.e. the 60 x 12 matrix Y (row z = Phi(z)), with the fixed
incidence N (row z = the two interiors of the far root z).  Known exactly from (ii), (iii):
  N^T N = 9I + J - P_s,  Y^T Y = 9I + J - P_mu,  N^T Y = K2 = 2J - I - P_s - P_mu - P_{tau s},
  Y^T P_s Y = K3 = 2J - P_s - P_{tau mu} - P_{mu tau} - P_tau,  N^T P_s = P_s N^T  (s = sigma).
Hence M+- = Z^T (I +- P_s) Z with Z = [N | Y] is determined, and must be PSD (I +- P_s >= 0).  Test over random (mu, tau)."""
import numpy as np, random, sys, itertools
L = [(e, s) for e in range(6) for s in (0, 1)]
LI = {l: i for i, l in enumerate(L)}
sig = {l: (l[0], 1 - l[1]) for l in L}
def P(f):
    M = np.zeros((12, 12))
    for l in L: M[LI[l], LI[f[l]]] = 1
    return M
def rand_fpf(rng):
    xs = L[:]; rng.shuffle(xs); f = {}
    for i in range(0, 12, 2): f[xs[i]], f[xs[i + 1]] = xs[i + 1], xs[i]
    return f
def comp(f, g): return {l: f[g[l]] for l in L}          # (f o g)(l) = f(g(l))
I = np.eye(12); J = np.ones((12, 12))
def mats(mu, tau):
    Ps, Pm, Pt = P(sig), P(mu), P(tau)
    NN = 9 * I + J - Ps
    YY = 9 * I + J - Pm
    K2 = 2 * J - I - Ps - Pm - P(comp(tau, sig))
    K3 = 2 * J - Ps - P(comp(tau, mu)) - P(comp(mu, tau)) - Pt
    G = np.block([[NN, K2], [K2.T, YY]])
    Gs = np.block([[NN @ Ps, Ps @ K2], [K2.T @ Ps, K3]])
    return G, Gs
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
worst = []
bad = 0; tried = 0
for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 2000):
    mu, tau = rand_fpf(rng), rand_fpf(rng)
    if any(mu[b] == tau[b] == sig[b] for b in L): continue       # sibling pair at the block: not a free star
    tried += 1
    G, Gs = mats(mu, tau)
    ep = np.linalg.eigvalsh((G + Gs + (G + Gs).T) / 2).min()
    em = np.linalg.eigvalsh((G - Gs + (G - Gs).T) / 2).min()
    sym_err = max(abs(Gs - Gs.T).max(), 0)
    if min(ep, em) < -1e-9: bad += 1
    worst.append((min(ep, em), ep, em, sym_err))
worst.sort()
print('tried', tried, 'non-PSD', bad, 'smallest eigenvalues (min, +, -, asym):', [tuple(round(v, 4) for v in w) for w in worst[:5]])
