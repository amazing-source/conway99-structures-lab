"""Standalone verifier for the O_24 whole-class LP certificates (structures 27/09).
For class k (C1's list, data/finishC/o24_ypairs.json), rebuild from scratch:
  Es = the sign vectors of Golay words (one per +-pair) orthogonal to y1, y2   (= all norm-3 vectors of O_24 perp y, /+-);
  Ms = the signed octads (+-1 on an octad, even number of -1, first point +1) orthogonal to y1, y2 (= all norm-4 vectors of
       O_24 perp y usable as h: the (+-4,+-4,0^22) type would give 42 + 16 > 56 on the diagonal of the frame identity);
  the 302 equations  sum_E x_E E_i E_j + 4 sum_M z_M M_i M_j = (56 I - B1 B1^T - B2 B2^T)_ij  (i <= j),  sum z = 7,
  sum x = 42.  Any frame gives a 0/1 solution.
The certificate is an integer vector y (one entry per equation).  With d = -A^T y, every x, z in [0,1] satisfying the
equations has  sum_k min(0, d_k) <= d.(x,z) = -b.y ;  the certificate is valid iff  b.y + sum_k min(0, d_k) > 0.
Integer arithmetic throughout (Python ints).  Usage: python verify_o24_cert.py k [k ...]   (reads state/o24_cert_KKK.npz)"""
import json, sys
import numpy as np

YP = 'C:/home/researcher/Desktop/c1_recherche/data/finishC/o24_ypairs.json'


def main():
    D = json.load(open(YP))
    code = {0}
    for b in D['golay_basis']:
        code |= {c ^ b for c in code}
    assert len(code) == 4096
    words = sorted(code)
    for a in sys.argv[1:]:
        k = int(a)
        o = D['orbits'][k]
        B1, B2 = [int(v) for v in o['a1']], [int(v) for v in o['a2']]
        dot = lambda u, v: sum(p * q for p, q in zip(u, v))
        assert dot(B1, B1) == 56 == dot(B2, B2) and dot(B1, B2) == 0
        Es = []
        for W in words:
            if W >> 23 & 1:
                continue
            E = [-1 if W >> i & 1 else 1 for i in range(24)]
            if dot(E, B1) == 0 and dot(E, B2) == 0:
                Es.append(E)
        Ms = []
        for W in words:
            if bin(W).count('1') != 8:
                continue
            pts = [i for i in range(24) if W >> i & 1]
            for s in range(128):
                sg = [1] + [-1 if s >> j & 1 else 1 for j in range(7)]
                if sg.count(-1) % 2:
                    continue
                M = [0] * 24
                for p, v in zip(pts, sg):
                    M[p] = v
                if dot(M, B1) == 0 and dot(M, B2) == 0:
                    Ms.append(M)
        iu = [(i, j) for i in range(24) for j in range(i, 24)]
        R0 = [[(56 if i == j else 0) - B1[i] * B1[j] - B2[i] * B2[j] for j in range(24)] for i in range(24)]
        b = [R0[i][j] for i, j in iu] + [7, 42]
        cert = np.load('state/o24_cert_%03d.npz' % k)
        y = [int(v) for v in cert['y']]
        assert len(y) == len(b), (len(y), len(b))
        # d = -A^T y, column by column
        dmin = 0
        for E in Es:
            col = sum(y[r] * E[i] * E[j] for r, (i, j) in enumerate(iu)) + y[-1]          # sum x = 42 row
            dmin += min(0, -col)
        for M in Ms:
            col = sum(y[r] * 4 * M[i] * M[j] for r, (i, j) in enumerate(iu) if M[i] and M[j]) + y[-2]   # sum z = 7 row
            dmin += min(0, -col)
        val = sum(p * q for p, q in zip(b, y)) + dmin
        print('class %3d: |Es| = %d, |Ms| = %d, certificate value (x scale %d) = %d -> %s'
              % (k, len(Es), len(Ms), int(cert['scale']), val, 'VALID' if val > 0 else 'INVALID'), flush=True)


if __name__ == '__main__':
    main()
