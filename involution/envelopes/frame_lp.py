"""General frame LP (structures 27/09): for an explicit envelope U and a candidate pair (y1, y2), decide whether the frame
identity can hold, with an exact integer certificate when it cannot.  Generalises o24_globallp.py to any U.
Data: an integral basis (rows, rational coordinates) with Gram G = rows rows^T / scale (scale = 1 for U/NNN.json; 8 for O_24
written as a in Z^24 with norm a^2/8).  Vectors are integer coefficient vectors c; (c, c') = c^T G c'.
Variables: x_u in [0,1] on the norm-3 vectors orthogonal to y1, y2 (up to sign); z_h in [0,1] on the norm-4 vectors
orthogonal to y1, y2 with |(chi,h)| <= 2 for every characteristic vector chi of norm 8 (budget at the shadow roots).
Equations (exact, integers): for all basis vectors e_i, e_j (i <= j)
    sum_u x_u (u.e_i)(u.e_j) + sum_h z_h (h.e_i)(h.e_j) = 7 (e_i.e_j) - (y1.e_i)(y1.e_j) - (y2.e_i)(y2.e_j),
    sum x = 42, sum z = 7.
Any frame with these y gives a 0/1 solution.  Phase-1 LP; if infeasible, integer certificate y (one per equation):
    b.y + sum_k min(0, (-A^T y)_k) > 0  proves it.  Usage (module): P = prepare(rows, scale); frame_lp(P, y1c, y2c)."""
import json, sys, time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, identity
sys.path.insert(0, '.')
from build156 import int_rows, hnf, gram_of, lll, enum, inv_exact


def prepare(rows, scale=1, verbose=True):
    t0 = time.time()
    Mi, d = int_rows(rows)
    B = hnf(Mi)
    assert len(B) == 24
    M, _ = lll(B, d)
    G = [[F(sum(a * b for a, b in zip(M[i], M[j])), d * d * scale) for j in range(24)] for i in range(24)]
    assert all(x.denominator == 1 for r in G for x in r), 'not integral'
    G = [[int(x) for x in r] for r in G]
    sv = enum(G, 4)
    assert not any(v == 1 for z, v in sv)
    Gm = np.array(G, dtype=np.int64)
    Gi = inv_exact(G)
    dg = [G[k][k] % 2 for k in range(24)]
    c = [int(sum(Gi[a][b] * dg[b] for b in range(24))) % 2 for a in range(24)]
    near = enum(G, 2, cnum=[-x for x in c], cden=2)
    CH, seen = [], set()
    for z, v in near:
        w = tuple(2 * a + b for a, b in zip(z, c))
        if tuple(-t for t in w) not in seen:
            seen.add(w)
            CH.append(w)
    P = {'M': M, 'd': d, 'scale': scale, 'G': Gm,
         'R': np.array([z for z, v in sv if v == 2], dtype=np.int64).reshape(-1, 24),
         'N3': np.array([z for z, v in sv if v == 3], dtype=np.int64),
         'N4': np.array([z for z, v in sv if v == 4], dtype=np.int64),
         'CH': np.array(CH, dtype=np.int64)}
    if verbose:
        print('prepared: %d roots, %d norm-3, %d norm-4, %d chi8 (48 + r = %d)  %.1fs' % (
            2 * len(P['R']), 2 * len(P['N3']), 2 * len(P['N4']), 2 * len(P['CH']), 48 + 2 * len(P['R']),
            time.time() - t0), flush=True)
    return P


def coef_of(P, coords):
    """coordinates (rational, same frame as the input rows) -> integer coefficients in the reduced basis"""
    M, d = P['M'], P['d']
    Bq = [[F(x, d) for x in r] for r in M]
    # solve c B = coords  <=>  B^T c^T = coords^T
    A = [[Bq[j][i] for j in range(24)] + [F(coords[i])] for i in range(24)]
    n = 24
    for col in range(n):
        p = next(i for i in range(col, n) if A[i][col] != 0)
        A[col], A[p] = A[p], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for i in range(n):
            if i != col and A[i][col] != 0:
                f = A[i][col]
                A[i] = [x - f * y for x, y in zip(A[i], A[col])]
    c = [A[i][n] for i in range(n)]
    assert all(x.denominator == 1 for x in c), 'not a lattice vector'
    return np.array([int(x) for x in c], dtype=np.int64)


def frame_lp(P, y1, y2, save=None):
    t0 = time.time()
    G = P['G']
    y1, y2 = np.asarray(y1, dtype=np.int64), np.asarray(y2, dtype=np.int64)
    assert y1 @ G @ y1 == 7 and y2 @ G @ y2 == 7 and y1 @ G @ y2 == 0
    g1, g2 = G @ y1, G @ y2                                  # (y . e_i)
    N3, N4, CH = P['N3'], P['N4'], P['CH']
    Us = N3[(N3 @ g1 == 0) & (N3 @ g2 == 0)]
    H0 = N4[(N4 @ g1 == 0) & (N4 @ g2 == 0)]
    Hs = H0[(np.abs(CH @ G @ H0.T) <= 2).all(axis=0)] if len(H0) else H0
    nu, nh = len(Us), len(Hs)
    out = {'n_u': nu, 'n_h': nh}
    if nu < 42 or nh < 7:
        out['verdict'] = 'EXCLUDED (too few candidates)'
        return out
    iu = [(i, j) for i in range(24) for j in range(i, 24)]
    Ue, He = Us @ G, Hs @ G                                   # rows: (v . e_i)
    Aeq = np.vstack([np.hstack([np.array([Ue[:, i] * Ue[:, j] for i, j in iu]), np.array([He[:, i] * He[:, j] for i, j in iu])]),
                     np.concatenate([np.zeros(nu), np.ones(nh)])[None, :].astype(np.int64),
                     np.concatenate([np.ones(nu), np.zeros(nh)])[None, :].astype(np.int64)]).astype(np.int64)
    beq = np.array([7 * G[i, j] - g1[i] * g1[j] - g2[i] * g2[j] for i, j in iu] + [7, 42], dtype=np.int64)
    n = nu + nh
    m = len(beq)
    A = hstack([csr_matrix(Aeq.astype(float)), identity(m, format='csr'), -identity(m, format='csr')]).tocsr()
    lp = linprog(np.concatenate([np.zeros(n), np.ones(2 * m)]), A_eq=A, b_eq=beq.astype(float),
                 bounds=[(0, 1)] * n + [(0, None)] * (2 * m), method='highs')
    out['phase1'] = float(lp.fun) if lp.status == 0 else None
    if lp.status != 0:
        out['verdict'] = 'solver status %d' % lp.status
    elif lp.fun < 1e-7:
        out['verdict'] = 'LP FEASIBLE'
    else:
        yv = lp.eqlin.marginals
        out['verdict'] = 'LP infeasible numerically, no certificate'
        for Dn in (1, 2, 4, 8, 16, 64, 256, 1024, 4096, 65536, 10 ** 6):
            yi = np.rint(yv * Dn).astype(np.int64)
            di = -(Aeq.T @ yi)
            val = int(beq @ yi) + int(np.minimum(di, 0).sum())
            if val > 0:
                out['verdict'] = 'EXCLUDED (frame LP infeasible, exact certificate)'
                out['certificate'] = {'scale': Dn, 'value_times_scale': val}
                if save:
                    np.savez_compressed(save, y=yi, scale=Dn, y1=y1, y2=y2)
                break
    out['seconds'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    # control: O_24 class k of C1's list, written as a in Z^24 with norm a^2/8 (Golay coordinates of C1)
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    D = json.load(open('C:/home/researcher/Desktop/c1_recherche/data/finishC/o24_ypairs.json'))
    basis = D['golay_basis']
    gens = [[1] * 24]
    for w in basis:
        gens.append([2 if w >> i & 1 else 0 for i in range(24)])
    for i in range(24):
        for j in range(i + 1, 24):
            v = [0] * 24; v[i] = 4; v[j] = 4; gens.append(v)
            v = [0] * 24; v[i] = 4; v[j] = -4; gens.append(v)
    P = prepare(gens, scale=8)
    o = D['orbits'][k]
    y1, y2 = coef_of(P, o['a1']), coef_of(P, o['a2'])
    print('class %d:' % k, frame_lp(P, y1, y2))
