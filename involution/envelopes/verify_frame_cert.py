"""Standalone verifier for frame_lp certificates (structures 27/09).  For a certificate file (npz: y = integer multipliers,
scale, y1, y2 = coefficient vectors, basis = the reduced basis rows used (strings), den):
  1. the stored basis spans the same lattice as U/NNN.json (HNF equality), Gram integral, det 1, odd;
  2. y1, y2 have norm 7 and are orthogonal;
  3. candidates rebuilt from scratch: all norm-3 vectors orthogonal to y1, y2 (up to sign) and all norm-4 vectors
     orthogonal to y1, y2 with |(chi,h)| <= 2 for every characteristic vector chi of norm 8 (sanity: 4096 norm-3 vectors,
     48 + r characteristic vectors of norm 8, as the theta series forces);
  4. the equations  sum_u x (u.e_i)(u.e_j) + sum_h z (h.e_i)(h.e_j) = 7 (e_i.e_j) - sum_y (y.e_i)(y.e_j)  (i <= j),
     sum z = 7, sum x = 42, rebuilt in Python integers; the certificate is valid iff b.y + sum_k min(0, (-A^T y)_k) > 0.
Usage: python verify_frame_cert.py certfile.npz [...]"""
import json, os, sys
from fractions import Fraction as F
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build156 import int_rows, hnf, gram_of, det_exact, enum, inv_exact, rows_frac


def check(path):
    z = np.load(path, allow_pickle=False)
    i = int(os.path.basename(path)[:3])
    zb = z if 'basis' in z.files else np.load(os.path.join(os.path.dirname(path), '%03d_basis.npz' % i))
    den0 = int(zb['den'])
    rows = [[F(int(x), den0) for x in r] for r in zb['basis']]          # stored as integer numerators over den
    rec = json.load(open(os.path.join(HERE, 'U', '%03d.json' % i)))
    a, da = int_rows(rows)
    b, db = int_rows(rows_frac(rec['gens']))
    den = da * db // np.gcd(da, db)
    Ha = hnf([[x * (den // da) for x in r] for r in a])
    Hb = hnf([[x * (den // db) for x in r] for r in b])
    # reduced echelon forms of the same lattice agree after making pivots/rows canonical: compare via mutual membership
    def span_contains(H, v):
        v = list(v)
        for r in H:
            c = next(k for k, x in enumerate(r) if x)
            if v[c] % r[c]:
                return False
            q = v[c] // r[c]
            v = [p - q * s for p, s in zip(v, r)]
        return not any(v)
    same = all(span_contains(Ha, r) for r in Hb) and all(span_contains(Hb, r) for r in Ha)
    G = gram_of(a, da)
    ok_lat = same and det_exact(G) == 1 and any(G[k][k] % 2 for k in range(24))
    y1 = [int(v) for v in z['y1']]; y2 = [int(v) for v in z['y2']]
    ip = lambda u, v: sum(u[p] * G[p][q] * v[q] for p in range(24) for q in range(24) if u[p] and v[q])
    ok_y = ip(y1, y1) == 7 and ip(y2, y2) == 7 and ip(y1, y2) == 0
    sv = enum(G, 4)
    N3 = [list(v) for v, n in sv if n == 3]
    N4 = [list(v) for v, n in sv if n == 4]
    r = 2 * sum(1 for v, n in sv if n == 2)
    Gi = inv_exact(G)
    dg = [G[k][k] % 2 for k in range(24)]
    c = [int(sum(Gi[p][q] * dg[q] for q in range(24))) % 2 for p in range(24)]
    near = enum(G, 2, cnum=[-x for x in c], cden=2)
    CH = []
    seen = set()
    for w, n in near:
        t = tuple(2 * p + q for p, q in zip(w, c))
        if tuple(-s for s in t) not in seen:
            seen.add(t); CH.append(list(t))
    sane = (2 * len(N3) == 4096) and (2 * len(CH) == 48 + r)
    Ge = lambda v: [sum(G[p][q] * v[q] for q in range(24) if v[q]) for p in range(24)]      # (v . e_p)
    g1, g2 = Ge(y1), Ge(y2)
    Us = [u for u in N3 if sum(p * q for p, q in zip(u, g1)) == 0 and sum(p * q for p, q in zip(u, g2)) == 0]
    Hs = []
    CHe = [Ge(ch) for ch in CH]
    for h in N4:
        if sum(p * q for p, q in zip(h, g1)) or sum(p * q for p, q in zip(h, g2)):
            continue
        if all(abs(sum(p * q for p, q in zip(h, ce))) <= 2 for ce in CHe):
            Hs.append(h)
    iu = [(p, q) for p in range(24) for q in range(p, 24)]
    b = [7 * G[p][q] - g1[p] * g1[q] - g2[p] * g2[q] for p, q in iu] + [7, 42]
    yv = [int(v) for v in z['y']]
    assert len(yv) == len(b)
    dmin = 0
    for u in Us:
        e = Ge(u)
        col = sum(yv[k] * e[p] * e[q] for k, (p, q) in enumerate(iu)) + yv[-1]
        dmin += min(0, -col)
    for h in Hs:
        e = Ge(h)
        col = sum(yv[k] * e[p] * e[q] for k, (p, q) in enumerate(iu)) + yv[-2]
        dmin += min(0, -col)
    val = sum(p * q for p, q in zip(b, yv)) + dmin
    verdict = 'VALID' if (ok_lat and ok_y and sane and val > 0) else 'INVALID'
    return {'file': os.path.basename(path), 'lattice_ok': ok_lat, 'y_ok': ok_y, 'sanity(4096, 48+r)': sane,
            'n_u': len(Us), 'n_h': len(Hs), 'value_x_scale': val, 'verdict': verdict}


if __name__ == '__main__':
    for p in sys.argv[1:]:
        print(check(p), flush=True)
