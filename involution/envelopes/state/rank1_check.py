"""for several O24 classes: sparsest dual (L1) and its rank; report whether a single-vector certificate Y = g g^T / c exists"""
import json, sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack
sys.path.insert(0, '.')
from o24_ustage import C1_YPAIRS, golay_from_basis, sign_vectors
from o24_globallp import build
D = json.load(open(C1_YPAIRS)); code = golay_from_basis(D['golay_basis']); SV = sign_vectors(code)
for k in [int(a) for a in sys.argv[1:]]:
    B1, B2, Es, Ms = build(k, D, code, SV)
    iu = [(i, j) for i in range(24) for j in range(i, 24)]; m = len(iu)
    R0 = 56 * np.eye(24) - np.outer(B1, B1) - np.outer(B2, B2)
    row = lambda v: np.array([v[i] * v[j] if i == j else 2 * v[i] * v[j] for i, j in iu])
    AE = np.array([row(E) for E in Es]); AM = np.array([row(M) for M in Ms])
    rR = np.array([R0[i, j] if i == j else 2 * R0[i, j] for i, j in iu])
    A_ub = vstack([csr_matrix(np.hstack([-AE, AE])), csr_matrix(np.hstack([-AM, AM])), csr_matrix(np.concatenate([rR, -rR])[None, :])])
    b_ub = np.concatenate([-np.ones(len(Es)), np.zeros(len(Ms)), [42 - 1e-3]])
    lp = linprog(np.ones(2 * m), A_ub=A_ub, b_ub=b_ub, bounds=[(0, None)] * (2 * m), method='highs')
    if lp.status != 0:
        print('class %d: no certificate (LP feasible class?)' % k, flush=True); continue
    y = lp.x[:m] - lp.x[m:]
    Y = np.zeros((24, 24))
    for (i, j), v in zip(iu, y):
        Y[i, j] = Y[j, i] = v
    ev, V = np.linalg.eigh(Y)
    big = [e for e in ev if abs(e) > 1e-6]
    msg = 'class %d: L1 %.3f, nonzero eigenvalues %s' % (k, lp.fun, np.round(big, 4).tolist())
    if len(big) == 1:
        g = V[:, np.argmax(np.abs(ev))] * np.sqrt(abs(big[0]))
        g = g / np.min(np.abs(g[np.abs(g) > 1e-6]))
        msg += ' -> RANK 1, g = %s ; (B1.g, B2.g) = (%.1f, %.1f) ; min |E.g| over u-cand = %.2f ; |g|^2 = %.2f' % (
            np.round(g, 2).tolist(), B1 @ g, B2 @ g, np.min(np.abs(Es @ g)), g @ g)
    print(msg, flush=True)
