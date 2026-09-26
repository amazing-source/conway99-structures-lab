"""Orbit-level EVEN model with PSD cutting planes (no columns).

Unknown: the 42 x 42 orbit matrix Q (Q_XY = number of neighbours of a representative of X in orbit Y), symmetric, Q_XX = 0,
Q_{X,D(X)} = 2, other entries in {0,1} (a non-partner orbit has at most one neighbour of x).  Exact reformulation of the even
identity Q^2 + Q = 12I + 4J - Sh:
    P := 6I - 2Q - Sh + J  satisfies  P^2 = 14 P   <=>   0 <= P <= 14 I  and  ||P||_F^2 = 14 tr P = 2940,
because tr P = 210 is automatic, and eigenvalues in [0,14] with sum 210 and square-sum 14*210 must all be 0 or 14.
(With row sums 12 and the block sums F1, QJ = 12J and Q Sh = Sh Q = 8J - 2Sh hold, so P^2 = 14P is equivalent to the identity.)
The Frobenius condition is linear in the 0/1 variables; the two LMIs are imposed lazily by eigenvector cuts.

Structure (from a skeleton, general_model format): D-partners fixed (Q = 2), pivots give the through-t zeros, sisters Q = n_s
(allowed to be 1 only for dd pairs), F1 at every block, row sums 12.  Optional near structure R_t fixes near entries."""
import itertools, time, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp

N = 7
PAIRS = [(a, b) for a in range(N) for b in range(a + 1, N)]


def build_orbit(orbs):
    n = len(orbs)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    vid = {p: k for k, p in enumerate(pairs)}
    nv = len(pairs)
    lb = np.zeros(nv)
    ub = np.ones(nv)
    rows, lo, hi = [], [], []

    def var(i, j):
        return vid[(i, j) if i < j else (j, i)]
    # value of Q_ij = x_ij + [partner] (partners fixed to 2: we store x = 1 and add 1 as a constant)
    partner = {i: o['partner'] for i, o in enumerate(orbs)}
    const = np.zeros((n, n))
    for i in range(n):
        j = partner[i]
        if j is not None:
            const[i, j] = 1          # Q = x + const, with x fixed to 1 below
            lb[var(i, j)] = ub[var(i, j)] = 1
    sh = np.array([[len(set(orbs[i]['pair']) & set(orbs[j]['pair'])) for j in range(n)] for i in range(n)])
    both_dis = {P: all(o['pivot'] is None for o in orbs if o['pair'] == P) for P in PAIRS}
    for i in range(n):
        for j in range(i + 1, n):
            if orbs[i]['pair'] == orbs[j]['pair'] and not both_dis[orbs[i]['pair']]:
                ub[var(i, j)] = 0            # n_s = 0 unless dd (39.3)
    # through-t zeros for pivots: if i pivots at t, its only neighbours through t are in D(i)
    for i, o in enumerate(orbs):
        t = o['pivot']
        if t is None:
            continue
        for j in range(n):
            if j != i and j != partner[i] and t in orbs[j]['pair']:
                ub[var(i, j)] = 0
    # F1 and row sums (on Q = x + const)
    for i in range(n):
        for t in range(N):
            coefs = [(var(i, j), 1) for j in range(n) if j != i and t in orbs[j]['pair']]
            c0 = sum(const[i, j] for j in range(n) if j != i and t in orbs[j]['pair'])
            val = 2 if t in orbs[i]['pair'] else 4
            rows.append(coefs); lo.append(val - c0); hi.append(val - c0)
        coefs = [(var(i, j), 1) for j in range(n) if j != i]
        c0 = sum(const[i, j] for j in range(n) if j != i)
        rows.append(coefs); lo.append(12 - c0); hi.append(12 - c0)
    # Frobenius: sum_{X != Y} (1 - 2Q - Sh)^2 = 2940 - 42*25 = 1890, i.e. over X<Y: 945.
    # For a 0/1 variable x (Q = x + c): (1 - 2x - 2c - sh)^2 = (1-2c-sh)^2 + x*((1-2c-sh-2)^2 - (1-2c-sh)^2)
    fro_c = 0.0
    fro = []
    for (i, j), k in vid.items():
        a0 = 1 - 2 * const[i, j] - sh[i, j]
        fro_c += a0 * a0
        fro.append((k, (a0 - 2) ** 2 - a0 * a0))
    rows.append(fro); lo.append(945 - fro_c); hi.append(945 - fro_c)
    A = sp.lil_matrix((len(rows), nv))
    for r, cs in enumerate(rows):
        for k, c in cs:
            A[r, k] += c
    return A.tocsr(), np.array(lo, float), np.array(hi, float), lb, ub, vid, const, sh


def P_of(x, vid, const, sh, n=42):
    Q = const.copy()
    for (i, j), k in vid.items():
        Q[i, j] += x[k]
        Q[j, i] = Q[i, j] if (i, j) != (j, i) else Q[i, j]
    Q = np.triu(Q, 1)
    Q = Q + Q.T
    return 6 * np.eye(n) - 2 * Q - sh + np.ones((n, n)) - np.diag(np.diag(sh)) + np.diag(np.diag(sh)) * 0


def cut(v, vid, const, sh, n=42):
    """v^T P v as (coef over variables, constant)."""
    coef = np.zeros(len(vid))
    Pc = 6 * np.eye(n) - sh + np.ones((n, n))
    Qc = np.triu(const, 1)
    Qc = Qc + Qc.T
    Pc = Pc - 2 * Qc
    c0 = v @ Pc @ v
    for (i, j), k in vid.items():
        coef[k] = -2 * 2 * v[i] * v[j]          # -2 * (x_ij in both (i,j) and (j,i))
    return coef, c0


def run(orbs, rounds=200, tl=120, verbose=True, integral=True):
    A, lo, hi, lb, ub, vid, const, sh = build_orbit(orbs)
    nv = len(vid)
    C, clo, chi = [], [], []
    for rnd in range(rounds):
        cons = [LinearConstraint(A, lo, hi)]
        if C:
            cons.append(LinearConstraint(np.array(C), np.array(clo), np.array(chi)))
        r = milp(np.zeros(nv), constraints=cons, integrality=np.ones(nv) if integral else np.zeros(nv),
                 bounds=Bounds(lb, ub), options={'time_limit': tl})
        if r.status == 2:
            return 'dead', rnd
        if r.status != 0:
            return 'unknown(%d)' % r.status, rnd
        Q = const.copy()
        for (i, j), k in vid.items():
            Q[i, j] += r.x[k]
        Q = np.triu(Q, 1)
        Q = Q + Q.T
        P = 6 * np.eye(42) - 2 * Q - sh + np.ones((42, 42))
        np.fill_diagonal(P, 5)
        ev, V = np.linalg.eigh(P)
        added = 0
        for kk in range(42):
            if ev[kk] < -1e-6:
                c, c0 = cut(V[:, kk], vid, const, sh)
                C.append(c); clo.append(-c0); chi.append(np.inf); added += 1
            if ev[kk] > 14 + 1e-6:
                c, c0 = cut(V[:, kk], vid, const, sh)
                C.append(c); clo.append(-np.inf); chi.append(14 * (V[:, kk] @ V[:, kk]) - c0); added += 1
        if verbose:
            print('  round', rnd, 'eig range [%.3f, %.3f]' % (ev[0], ev[-1]), 'cuts', added, 'total', len(C), flush=True)
        if not added:
            return 'ALIVE (P^2 = 14P holds)', rnd
    return 'unknown(rounds)', rounds


if __name__ == '__main__':
    import random
    from alldd import skeleton_from_pair_cycles, random_cycles
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    integral = not (len(sys.argv) > 3 and sys.argv[3] == 'lp')
    for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 1):
        orbs = skeleton_from_pair_cycles(random_cycles(rng))
        t0 = time.time()
        print('all-dd trial', t, run(orbs, integral=integral), '%.0fs' % (time.time() - t0), flush=True)
