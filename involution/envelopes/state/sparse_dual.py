"""sparsest dual certificate for an O24 class: find symmetric Y with E Y E^T >= 1 (u-candidates), M Y M^T >= 0 (h-candidates),
<Y, R0> <= 42 - gap, minimising sum |Y_ij| (L1).  Prints the support structure relative to T1, T2, D."""
import json, sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, vstack, identity
sys.path.insert(0, '.')
from o24_ustage import C1_YPAIRS, golay_from_basis, sign_vectors
from o24_globallp import build
D = json.load(open(C1_YPAIRS)); code = golay_from_basis(D['golay_basis']); SV = sign_vectors(code)
k = int(sys.argv[1]); gap = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
o = D['orbits'][k]
B1, B2, Es, Ms = build(k, D, code, SV)
iu = [(i, j) for i in range(24) for j in range(i, 24)]
m = len(iu)
R0 = 56 * np.eye(24) - np.outer(B1, B1) - np.outer(B2, B2)
# Y_ij for i<=j; E Y E^T = sum_i Y_ii + 2 sum_{i<j} Y_ij E_i E_j
def row(v, w=1.0):
    return np.array([w * (v[i] * v[j] if i == j else 2 * v[i] * v[j]) for i, j in iu])
AE = np.array([row(E) for E in Es]); AM = np.array([row(M) for M in Ms])
rR = np.array([R0[i, j] if i == j else 2 * R0[i, j] for i, j in iu])
# variables: Y (m, free) = Yp - Yn, Yp, Yn >= 0 ; minimise sum(Yp + Yn)
A_ub = vstack([csr_matrix(np.hstack([-AE, AE])), csr_matrix(np.hstack([-AM, AM])), csr_matrix(np.concatenate([rR, -rR])[None, :])])
b_ub = np.concatenate([-np.ones(len(Es)), np.zeros(len(Ms)), [42 - gap]])
lp = linprog(np.ones(2 * m), A_ub=A_ub, b_ub=b_ub, bounds=[(0, None)] * (2 * m), method='highs')
print('class', k, 'status', lp.status, lp.message[:50])
if lp.status == 0:
    Y = lp.x[:m] - lp.x[m:]
    nz = [(iu[t], round(float(Y[t]), 4)) for t in range(m) if abs(Y[t]) > 1e-7]
    T1, T2 = set(o['T1']), set(o['T2'])
    lab = lambda i: 'T1' if i in T1 else 'T2' if i in T2 else '..'
    print('L1 = %.4f, nonzeros %d' % (lp.fun, len(nz)))
    for (i, j), v in sorted(nz, key=lambda t: -abs(t[1]))[:40]:
        print('  Y[%2d,%2d] %s-%s  B1=(%+d,%+d) B2=(%+d,%+d)  %+.4f' % (i, j, lab(i), lab(j), B1[i], B1[j], B2[i], B2[j], v))
