"""export optimal dual matrices Y (max sum x under the frame identity, x,z >= 0) for O24 classes, and the L1-sparsest certificate"""
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
    if lp.status:
        print('class', k, 'no certificate'); continue
    y = lp.x[:m] - lp.x[m:]
    Y = np.zeros((24, 24))
    for (i, j), v in zip(iu, y):
        Y[i, j] = Y[j, i] = v
    o = D['orbits'][k]
    np.savez_compressed('state/o24_Y_%03d.npz' % k, Y=Y, B1=B1, B2=B2, T1=o['T1'], T2=o['T2'], E_candidates=Es, value=float(np.sum(Y * R0)))
    print('class %d: saved, <Y,R0> = %.4f, rank %d' % (k, np.sum(Y * R0), int(np.sum(np.abs(np.linalg.eigvalsh(Y)) > 1e-6))))
