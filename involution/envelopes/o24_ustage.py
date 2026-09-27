"""O_24 (id 0), u-stage (structures 27/09).  Given the two y and the 7 signed octads of a frame candidate in O_24, decide
whether 42 norm-3 vectors u can complete the frame.  Coordinates scaled by sqrt8: u = E (sign vector, negative set in the
Golay code), h_c = 2 M_c, y_j = B_j.
Exact necessary conditions used:
  E.B_j = 0;  E.M_c in {0, +-4} with exactly two +-4 (u lies in exactly the stars of its two blocks, (u,h_c) = -C_cu);
  for every pair {c,d}: exactly one u with (E.M_c)(E.M_d) = +16 and one with -16 (the two orbits over {c,d} have
  c-vectors e_c + e_d and e_c - e_d);
  the frame identity  sum_u E^T E = R := 56 I - 4 M^T M - B^T B  (all 300 entries; the star relations follow from it).
LP relaxation over x_E in [0,1] (one variable per candidate up to sign).  If the LP is infeasible, the dual of the phase-1 LP
is rationalised and verified EXACTLY: b.y - sum_E max(0, (A^T y)_E) > 0 proves that no 0/1 (indeed no [0,1]) solution exists.
Golay code = C1's (data/finishC/o24_ypairs.json golay_basis), the coordinates of C1's and C2's files.
Usage (module): ustage(B1, B2, M) -> dict;  python o24_ustage.py --selftest"""
import json, sys, time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog

C1_YPAIRS = 'C:/home/researcher/Desktop/c1_recherche/data/finishC/o24_ypairs.json'


def golay_from_basis(basis):
    code = {0}
    for b in basis:
        code |= {c ^ b for c in code}
    assert len(code) == 4096
    return sorted(code)


def sign_vectors(code):
    """one sign vector per +- pair: negative set W with bit 23 clear (W and its complement give E and -E)"""
    reps = [W for W in code if not W >> 23 & 1]
    assert len(reps) == 2048
    return np.array([[-1 if W >> i & 1 else 1 for i in range(24)] for W in reps], dtype=np.int64)


def candidates(SV, B1, B2, M):
    """rows of SV admissible as u: returns list of (row index, c, d, type) with type +1 if E.M_c, E.M_d have the same sign"""
    p1, p2 = SV @ B1, SV @ B2
    PM = SV @ M.T                                     # 2048 x 7
    out = []
    ok = (p1 == 0) & (p2 == 0) & np.isin(PM, (0, 4, -4)).all(axis=1) & ((PM != 0).sum(axis=1) == 2)
    for r in np.nonzero(ok)[0]:
        c, d = np.nonzero(PM[r])[0]
        out.append((int(r), int(c), int(d), 1 if PM[r, c] * PM[r, d] > 0 else -1))
    return out


def ustage(B1, B2, M, SV, verbose=True):
    B1, B2, M = np.array(B1, dtype=np.int64), np.array(B2, dtype=np.int64), np.array(M, dtype=np.int64)
    R = 56 * np.eye(24, dtype=np.int64) - 4 * M.T @ M - np.outer(B1, B1) - np.outer(B2, B2)
    cand = candidates(SV, B1, B2, M)
    res = {'n_candidates': len(cand)}
    groups = {}
    for k, (r, c, d, t) in enumerate(cand):
        groups.setdefault((c, d, t), []).append(k)
    missing = [(c, d, t) for c in range(7) for d in range(c + 1, 7) for t in (1, -1) if (c, d, t) not in groups]
    res['empty_pair_types'] = len(missing)
    if missing:
        res['verdict'] = 'EXCLUDED (a pair {c,d} has no candidate of one type)'
        res['proof'] = {'missing': missing[:5]}
        return res
    # equality system A x = b
    iu = [(i, j) for i in range(24) for j in range(i, 24)]
    E = SV[[r for r, c, d, t in cand]]                # n x 24
    A1 = np.array([E[:, i] * E[:, j] for i, j in iu], dtype=np.int64)       # 300 x n
    b1 = np.array([R[i, j] for i, j in iu], dtype=np.int64)
    keys = sorted(groups)
    A2 = np.zeros((len(keys), len(cand)), dtype=np.int64)
    for q, key in enumerate(keys):
        A2[q, groups[key]] = 1
    A = np.vstack([A1, A2])
    b = np.concatenate([b1, np.ones(len(keys), dtype=np.int64)])
    n = len(cand)
    m = A.shape[0]
    # phase 1: min sum(s+ + s-)  s.t. A x + s+ - s- = b, 0 <= x <= 1, s >= 0
    Aeq = np.hstack([A, np.eye(m), -np.eye(m)]).astype(float)
    cobj = np.concatenate([np.zeros(n), np.ones(2 * m)])
    bounds = [(0, 1)] * n + [(0, None)] * (2 * m)
    lp = linprog(cobj, A_eq=Aeq, b_eq=b.astype(float), bounds=bounds, method='highs')
    res['phase1'] = float(lp.fun)
    if lp.fun < 1e-7:
        res['verdict'] = 'LP FEASIBLE (integer stage needed)'
        res['lp_x'] = lp.x[:n].tolist()
        return res
    # exact Farkas-type certificate from the duals of the equality rows
    y = lp.eqlin.marginals
    for den in (1, 2, 4, 8, 16, 32, 64, 128, 256, 1024, 4096, 10 ** 6):
        yq = [F(float(v)).limit_denominator(den) for v in y]
        yq = [max(F(-1), min(F(1), v)) for v in yq]
        Aty = [sum(yq[i] * int(A[i, k]) for i in range(m) if A[i, k]) for k in range(n)]
        val = sum(yq[i] * int(b[i]) for i in range(m)) - sum(max(F(0), v) for v in Aty)
        if val > 0:
            res['verdict'] = 'EXCLUDED (LP infeasible, exact certificate)'
            res['proof'] = {'value': str(val), 'den': den, 'y': [str(v) for v in yq]}
            return res
    res['verdict'] = 'LP infeasible numerically but no exact certificate found (inspect)'
    return res


if __name__ == '__main__':
    D = json.load(open(C1_YPAIRS))
    code = golay_from_basis(D['golay_basis'])
    SV = sign_vectors(code)
    print('2048 sign vectors ready')
