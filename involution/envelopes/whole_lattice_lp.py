"""Whole-lattice relaxation of the frame (structures 27/09), FEASIBILITY STUDY (floats; no certificate yet).
For an explicit envelope U (U/NNN.json): variables x_u in [0,1] on the norm-3 vectors (up to sign), z_h in [0,1] on the
admissible norm-4 vectors (|(beta,h)| <= 1 for every shadow root beta = chi_8/2; +-2 is excluded by the budget).
Y := 7I - sum x uu^T - sum z hh^T stands for y1 y1^T + y2 y2^T (unknown).  Exact facts used (all proven):
  sum x = 42, sum z = 7;
  per shadow root beta: Y(beta) = beta^T Y beta = 3.5 - 2 U3(beta) - H(beta) in [1/2, 5/2] (budget 14 at beta:
     u's give 10.5 + 2 U3, U3 = number of u at |.| = 3/2; h's give H; y's give 1/2 or 5/2), hence U3 <= 1 (packing);
  per root alpha: alpha^T Y alpha >= 0;
  Y >= 0 (PSD) and sum z hh^T <= 4 I (7 orthogonal norm-4 vectors), added as eigenvector cuts.
Usage: python whole_lattice_lp.py id [rounds]"""
import json, sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack as svstack
sys.path.insert(0, '.')
from build156 import Lattice, enum, rows_frac, inv_exact

T0 = time.time()


def log(*a):
    print('[%5.1fs]' % (time.time() - T0), *a, flush=True)


def setup(i):
    rec = json.load(open('U/%03d.json' % i))
    L = Lattice(rows_frac(rec['gens']))
    G = np.array(L.G, dtype=np.int64)
    sv = enum(L.G, 4)
    R = np.array([z for z, v in sv if v == 2], dtype=np.int64).reshape(-1, 24)
    N3 = np.array([z for z, v in sv if v == 3], dtype=np.int64)
    N4 = np.array([z for z, v in sv if v == 4], dtype=np.int64)
    # characteristic vector chi0 (coefficients mod 2): G c = diag(G) mod 2
    Gi = inv_exact(L.G)
    dg = [L.G[k][k] % 2 for k in range(24)]
    c = [int(sum(Gi[a][b] * dg[b] for b in range(24))) % 2 for a in range(24)]
    near = enum(L.G, 2, cnum=[-x for x in c], cden=2)           # z with |z + c/2|^2 <= 2
    S2 = []
    seen = set()
    for z, v in near:
        assert v == 2, v                                          # shadow norms are even and nonzero
        w = tuple(2 * a + b for a, b in zip(z, c))                # doubled shadow vector (= a characteristic vector chi)
        if tuple(-t for t in w) not in seen:
            seen.add(w)
            S2.append(w)
    CH = np.array(S2, dtype=np.int64)                              # chi = 2 beta, norm 8
    B = np.array([[float(x) for x in r] for r in L.rows()])       # coordinates of the basis (Q^24)
    return rec, G, R, N3, N4, CH, B


def main():
    i = int(sys.argv[1])
    rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    rec, G, R, N3, N4, CH, B = setup(i)
    r = 2 * len(R)
    log('id %d %s: r = %d, norm3 %d, norm4 %d, chi8 %d (48 + r = %d)' % (i, rec['roots'], r, 2 * len(N3), 2 * len(N4),
                                                                      2 * len(CH), 48 + r))
    PCh3 = CH @ G @ N3.T                                           # (chi, u) in odd integers; |.| = 3 <=> |beta.u| = 3/2
    PCh4 = CH @ G @ N4.T                                           # (chi, h) even
    adm = (np.abs(PCh4) <= 2).all(axis=0)
    H = N4[adm]
    PCh4 = PCh4[:, adm]
    log('admissible h: %d of %d (up to sign)' % (len(H), len(N4)))
    nu, nh = len(N3), len(H)
    n = nu + nh
    rows, lo, hi = [], [], []
    # sum x = 42, sum z = 7
    rows.append(np.concatenate([np.ones(nu), np.zeros(nh)])); lo.append(42); hi.append(42)
    rows.append(np.concatenate([np.zeros(nu), np.ones(nh)])); lo.append(7); hi.append(7)
    # shadow roots: 1 <= H + 2 U3 <= 3  (H = sum z (beta.h)^2 = sum z (chi.h)^2 / 4)
    for k in range(len(CH)):
        row = np.concatenate([2.0 * (np.abs(PCh3[k]) == 3), (PCh4[k] ** 2) / 4.0])
        rows.append(row); lo.append(1); hi.append(3)
        rows.append(np.concatenate([np.zeros(nu), (PCh4[k] ** 2) / 4.0])); lo.append(1); hi.append(3)     # H in [1,3]
    # roots: sum x (a.u)^2 + sum z (a.h)^2 <= 14
    if len(R):
        PR3 = (R @ G @ N3.T) ** 2
        PR4 = (R @ G @ H.T) ** 2
        for k in range(len(R)):
            rows.append(np.concatenate([PR3[k], PR4[k]]).astype(float)); lo.append(-np.inf); hi.append(14)
    # Coxeter numbers of the even neighbours
    import re
    def coxeter(name):
        m = re.match(r'([ADE])_(?:\{(\d+)\}|(\d))', name)
        if not m:
            return 0                                               # Leech
        t, k = m.group(1), int(m.group(2) or m.group(3))
        return k + 1 if t == 'A' else 2 * k - 2 if t == 'D' else {6: 12, 7: 18, 8: 30}[k]
    hj = [coxeter(rec['nb1']), coxeter(rec['nb2'])]
    assert r == 8 * (hj[0] + hj[1] - 2), (r, hj)
    if len(R):
        Eu = PR3.sum(axis=0)                                        # E(u) = sum_{alpha>0} (alpha,u)^2 (R holds one per +-)
        Eh = PR4.sum(axis=0)
        mu_ = r / 8.0 - Eu
        assert (mu_ >= 0).all()                                     # lattice identity m(u) = r/8 - E(u) >= 0
        # sum_{alpha>0} alpha^T Y alpha = 7 r - sum x E(u) - sum z E(h) =: E(y1) + E(y2)
        for h_ in hj:
            # sum_u x m(u) + 2 (4 h + r/8) - (E1 + E2) <= 24 h - r
            row = np.concatenate([mu_ + Eu, Eh]).astype(float)
            rows.append(row); lo.append(-np.inf); hi.append(24 * h_ - r - 8 * h_ - r / 4.0 + 7 * r)
        # E(y) <= 4 h_min + r/8 for each y:  E1 + E2 <= 8 h_min + r/4
        row = np.concatenate([Eu, Eh]).astype(float)
        rows.append(row); lo.append(7 * r - 8 * min(hj) - r / 4.0); hi.append(np.inf)
    U3c = N3 @ B                                                   # coordinates
    H4c = H @ B
    cuts = 0
    for rd in range(rounds):
        Rw = np.array(rows)
        lo_a, hi_a = np.array(lo, dtype=float), np.array(hi, dtype=float)
        eq = np.isfinite(lo_a) & np.isfinite(hi_a) & (lo_a == hi_a)
        up = np.isfinite(hi_a) & ~eq
        dn = np.isfinite(lo_a) & ~eq
        Aub = svstack([csr_matrix(Rw[up]), csr_matrix(-Rw[dn])])
        bub = np.concatenate([hi_a[up], -lo_a[dn]])
        lp = linprog(np.zeros(n), A_ub=Aub, b_ub=bub, A_eq=csr_matrix(Rw[eq]), b_eq=lo_a[eq], bounds=[(0, 1)] * n,
                     method='highs')
        if lp.status != 0:
            log('round %d: LP INFEASIBLE (%s) after %d eigen-cuts -> id %d would be EXCLUDED (certificate to extract)'
                % (rd, lp.message[:60], cuts, i))
            return
        x, z = lp.x[:nu], lp.x[nu:]
        Xu = (U3c.T * x) @ U3c
        Xh = (H4c.T * z) @ H4c
        Y = 7 * np.eye(24) - Xu - Xh
        ev, V = np.linalg.eigh(Y)
        evh, Vh = np.linalg.eigh(Xh)
        added = 0
        for k in range(24):
            if ev[k] < -1e-6:                                      # PSD cut: sum x (u.w)^2 + sum z (h.w)^2 <= 7
                w = V[:, k]
                rows.append(np.concatenate([(U3c @ w) ** 2, (H4c @ w) ** 2])); lo.append(-np.inf); hi.append(7.0)
                added += 1
            if evh[k] > 4 + 1e-6:                                  # Bessel cut on the h: sum z (h.w)^2 <= 4
                w = Vh[:, k]
                rows.append(np.concatenate([np.zeros(nu), (H4c @ w) ** 2])); lo.append(-np.inf); hi.append(4.0)
                added += 1
        cuts += added
        log('round %d: min eig Y %.4f, max eig Xh %.4f, cuts added %d, z support %d' % (rd, ev[0], evh[-1], added,
                                                                                     int((z > 1e-7).sum())))
        if not added:
            log('LP FEASIBLE with all eigen-cuts satisfied: relaxation too weak for id %d' % i)
            return
    log('stopped after %d rounds (still feasible)' % rounds)


if __name__ == '__main__':
    main()
