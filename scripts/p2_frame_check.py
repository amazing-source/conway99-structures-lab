"""Check of the reformulation of (P2): with U = V + B/2 (V: 42x21 orbit values, B: orbit-pair incidence), the 210
identities (P2) and the diagonal (PPE) say  U^T U = M := 16.5 I + 16 J - 4 A_T7  (A_T7: pairs sharing a block).
On a solution of the FIRST-ORDER model (no P2), the components of U^T U on 1 and on the 6-dim block-sum space must
already equal those of M (automatic by mu); the 14-dim part carries the content of (P2)."""
import itertools, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, 'scripts')
import config_model as cm
from ppe_typed_lp import SAMPLES

dec = SAMPLES[2]
cols = cm.columns(dec)
A, lo, hi, _ = cm.build(dec, p2=False, cols=cols)
res = milp(c=np.zeros(A.shape[1]), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(A.shape[1]),
           bounds=Bounds(0, 1), options={'time_limit': 150})
assert res.status == 0
chosen = [cols[c] for c in range(len(cols)) if round(res.x[c])]
PAIRS = cm.PAIRS
orbits = [c[0] for c in chosen]
V = np.array([[c[2][Q] for Q in PAIRS] for c in chosen], float)
Bm = np.array([[1.0 if frozenset(Z) == Q else 0.0 for Q in PAIRS] for Z in orbits])
U = V + Bm / 2
AT = np.array([[1.0 if len(P & Q) == 1 else 0.0 for Q in PAIRS] for P in PAIRS])
M = 16.5 * np.eye(21) + 16 * np.ones((21, 21)) - 4 * AT
G = U.T @ U
w, E = np.linalg.eigh(AT)
proj = {lam: E[:, np.abs(w - lam) < 1e-9] for lam in (10.0, 3.0, -2.0)}
for lam, P in proj.items():
    D = P.T @ (G - M) @ P
    print(f'A_T7 eigenspace {lam:+.0f} (dim {P.shape[1]}): max |(U^TU - M)| there = {np.abs(D).max():.3e}')
X = proj[10.0].T @ (G - M) @ proj[-2.0]
Y = proj[3.0].T @ (G - M) @ proj[-2.0]
print(f'cross blocks: 1 x W14 {np.abs(X).max():.3e};  W6 x W14 {np.abs(Y).max():.3e}')
d = np.array([sum(1 for Q in PAIRS if V[z, PAIRS.index(Q)] == 2 and not (Q & set(orbits[z]))) for z in range(42)])
W14 = U @ proj[-2.0]
print('check |w_Z|^2 = 37/6 + 2 d_Z :', np.allclose((W14 ** 2).sum(axis=1), 37 / 6 + 2 * d))
