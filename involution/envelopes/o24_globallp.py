"""O_24 (id 0), whole-class LP (structures 27/09): for one (y1, y2) class, relax the frame to x_E in [0,1] (sign words
perp y, up to sign) and z_M in [0,1] (signed octads perp y, up to sign), with the EXACT linear frame identity
    sum_E x_E E^T E + 4 sum_M z_M M^T M = 56 I - B1^T B1 - B2^T B2        (300 equations)
plus valid linear cuts:  sum z = 7, sum x = 42;
  (i)   x_E + z_M <= 1 when |E.M| = 8            ((u,h) in {0,+-1});
  (ii)  sum_{E : |E.M| = 4} x_E >= 12 z_M         (the star of h has 12 u);
  (iii) sum_{M : |E.M| = 4} z_M >= 2 x_E          (u lies in exactly two stars);
  (iv)  sum_{M : E.M != 0} z_M <= 2 + 7 (1 - x_E);
  (v)   sum_{E : |E.M| = 4} x_E <= 12 + 30 (1 - z_M).
If the LP is infeasible, an exact certificate is extracted from the phase-1 duals and verified in rationals:
with c = A^T y + G^T lam (lam >= 0 on the <= rows), every feasible point has c.x <= y.b + lam.h, while
c.x >= sum min(0, c_k) on [0,1]^n; so sum_k min(0, c_k) > y.b + lam.h is a proof that the class has no frame.
Usage: python o24_globallp.py class_index [class_index ...]"""
import json, sys, time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix, vstack, hstack, identity
from o24_ustage import C1_YPAIRS, golay_from_basis, sign_vectors


def build(k, D, code, SV):
    OCT = [W for W in code if bin(W).count('1') == 8]
    orb = D['orbits'][k]
    B1, B2 = np.array(orb['a1'], dtype=np.int64), np.array(orb['a2'], dtype=np.int64)
    Es = SV[(SV @ B1 == 0) & (SV @ B2 == 0)]
    Ms = []
    for W in OCT:
        pts = [i for i in range(24) if W >> i & 1]
        for s in range(128):
            signs = [1] + [-1 if s >> j & 1 else 1 for j in range(7)]
            if signs.count(-1) % 2:
                continue
            v = np.zeros(24, dtype=np.int64)
            v[pts] = signs
            if v @ B1 == 0 and v @ B2 == 0:
                Ms.append(v)
    Ms = np.array(Ms)
    return B1, B2, Es, Ms


def solve(k, D, code, SV, verbose=True, cuts=True):
    t0 = time.time()
    B1, B2, Es, Ms = build(k, D, code, SV)
    nE, nM = len(Es), len(Ms)
    n = nE + nM
    iu = [(i, j) for i in range(24) for j in range(i, 24)]
    R0 = 56 * np.eye(24, dtype=np.int64) - np.outer(B1, B1) - np.outer(B2, B2)
    # equality rows
    AE = np.array([Es[:, i] * Es[:, j] for i, j in iu], dtype=np.int64)
    AM = np.array([4 * Ms[:, i] * Ms[:, j] for i, j in iu], dtype=np.int64)
    Aeq = np.vstack([np.hstack([AE, AM]),
                     np.concatenate([np.zeros(nE), np.ones(nM)])[None, :].astype(np.int64),
                     np.concatenate([np.ones(nE), np.zeros(nM)])[None, :].astype(np.int64)])
    beq = np.concatenate([np.array([R0[i, j] for i, j in iu]), [7, 42]]).astype(np.int64)
    # inequality rows G x <= h (sparse)
    PE = Es @ Ms.T                                 # nE x nM
    rows, cols, vals, hs = [], [], [], []
    r = 0
    e8, m8 = np.nonzero(np.abs(PE) == 8)
    for e, m in zip(e8, m8):                       # (i)
        rows += [r, r]; cols += [int(e), nE + int(m)]; vals += [1, 1]; hs.append(1); r += 1
    A4 = (np.abs(PE) == 4)
    NZ = (PE != 0)
    for m in range(nM):                            # (ii): 12 z_M - sum x_E <= 0 ; (v): sum x_E + 30 z_M <= 42
        es = np.nonzero(A4[:, m])[0]
        rows += [r] * (len(es) + 1); cols += [int(e) for e in es] + [nE + m]; vals += [-1] * len(es) + [12]; hs.append(0); r += 1
        rows += [r] * (len(es) + 1); cols += [int(e) for e in es] + [nE + m]; vals += [1] * len(es) + [30]; hs.append(42); r += 1
    for e in range(nE):                            # (iii): 2 x_E - sum z_M <= 0 ; (iv): sum_{nz} z_M + 7 x_E <= 9
        ms = np.nonzero(A4[e])[0]
        rows += [r] * (len(ms) + 1); cols += [nE + int(m) for m in ms] + [e]; vals += [-1] * len(ms) + [2]; hs.append(0); r += 1
        ms2 = np.nonzero(NZ[e])[0]
        rows += [r] * (len(ms2) + 1); cols += [nE + int(m) for m in ms2] + [e]; vals += [1] * len(ms2) + [7]; hs.append(9); r += 1
    if not cuts:
        rows, cols, vals, hs, r = [], [], [], [], 0
    G = coo_matrix((vals, (rows, cols)), shape=(r, n)).tocsr()
    h = np.array(hs, dtype=np.int64)
    m_eq = Aeq.shape[0]
    # phase 1: slacks on the equality rows only
    Aeq_s = hstack([csr_matrix(Aeq.astype(float)), identity(m_eq, format='csr'), -identity(m_eq, format='csr')]).tocsr()
    G_s = hstack([G.astype(float), csr_matrix((r, 2 * m_eq))]).tocsr()
    cobj = np.concatenate([np.zeros(n), np.ones(2 * m_eq)])
    bounds = [(0, 1)] * n + [(0, None)] * (2 * m_eq)
    if r:
        lp = linprog(cobj, A_ub=G_s, b_ub=h.astype(float), A_eq=Aeq_s, b_eq=beq.astype(float), bounds=bounds, method='highs')
    else:
        lp = linprog(cobj, A_eq=Aeq_s, b_eq=beq.astype(float), bounds=bounds, method='highs')
    out = {'class': k, 'cuts': cuts, 'nE': nE, 'nM': nM, 'ineq_rows': r, 'status': lp.status, 'phase1': float(lp.fun) if lp.status == 0 else None}
    if lp.status == 0 and lp.fun < 1e-7:
        out['verdict'] = 'LP FEASIBLE'
        z = lp.x[nE:n]
        out['z_support'] = int((z > 1e-6).sum())
    elif lp.status == 0:
        # duals: y = d f/d b_eq, mu = d f/d h <= 0.  d = -A^T y - G^T mu (reduced costs of the x columns).
        # For every point with A x = b, G x <= h, 0 <= x <= 1:  sum min(0, d_k) <= d.x <= -b.y - h.mu.
        # Certificate: b.y + h.mu + sum min(0, d_k) > 0 (at the phase-1 optimum it equals the optimum, here > 0).
        y = lp.eqlin.marginals
        mu = np.minimum(lp.ineqlin.marginals, 0.0) if r else np.zeros(0)
        Gt = G.T.tocsr()
        dfl = -(Aeq.T.astype(float) @ y) - Gt @ mu
        valf = float(beq @ y + h @ mu + np.minimum(dfl, 0).sum())
        out['float_certificate'] = valf
        ok = False
        GI = G.astype(np.int64)
        for Dn in (1, 2, 4, 8, 16, 64, 256, 1024, 4096, 65536, 10 ** 6, 10 ** 9):
            yi = np.rint(y * Dn).astype(np.int64)
            mi = np.minimum(np.rint(mu * Dn).astype(np.int64), 0)
            di = -(Aeq.T @ yi) - (GI.T @ mi)
            val = int(beq @ yi) + int(h @ mi) + int(np.minimum(di, 0).sum())      # exact integers (scaled by Dn)
            if val > 0:
                ok = True
                out['verdict'] = 'EXCLUDED (whole-class LP infeasible, exact certificate)'
                out['certificate'] = {'scale': Dn, 'value_times_scale': val, 'n_eq_used': int((yi != 0).sum()),
                                      'n_ineq_used': int((mi != 0).sum())}
                np.savez_compressed('state/o24_cert_%03d.npz' % k, y=yi, mu=mi, scale=Dn)
                break
        if not ok:
            out['verdict'] = 'LP infeasible numerically, certificate not found'
    else:
        out['verdict'] = 'solver status %d: %s' % (lp.status, lp.message)
    out['seconds'] = round(time.time() - t0, 1)
    if verbose:
        print(out, flush=True)
    return out


if __name__ == '__main__':
    D = json.load(open(C1_YPAIRS))
    code = golay_from_basis(D['golay_basis'])
    SV = sign_vectors(code)
    res_path = 'o24_globallp_results.json'
    import os
    allres = json.load(open(res_path)) if os.path.exists(res_path) else {}
    use_cuts = '--cuts' in sys.argv
    for a in [x for x in sys.argv[1:] if not x.startswith('--')]:
        k = int(a)
        if str(k) in allres and (allres[str(k)]['verdict'].startswith('EXCLUDED') or not use_cuts):
            continue
        r = solve(k, D, code, SV, cuts=False)
        if use_cuts and not r['verdict'].startswith('EXCLUDED'):
            r = solve(k, D, code, SV, cuts=True)
        allres[str(k)] = r
        json.dump(allres, open(res_path, 'w'), indent=0)
    print('excluded so far: %d / 141' % sum(1 for v in allres.values() if v['verdict'].startswith('EXCLUDED')))
