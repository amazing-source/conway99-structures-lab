"""Independent exact verification (own code, Fractions) of C1's triangle lemma (message of 25/09).

Rooted model at x: 14 labels, classes {2c, 2c+1} (c = 0..6); pairs {u,v} of labels in different classes.
Q = the 12 pairs whose labels lie in two distinct classes among alpha = 0, beta = 1, gamma = 2 (labels 0..5).
R = {x} ∪ Γ(x) (15 vertices), W = R ∪ Q.  S = 3I + (11/99)J - A  ( = 7 E_{-4} >= 0 ).
Claims:
 (1) S_RR >= 0, S_RQ ⊂ range(S_RR), and 3I + (11/99)J_Q - S_QR S_RR^+ S_RQ = 3I - P_d, where P_d is the orthogonal
     projection on {pairs -> f(u)+f(v)} (d = 6x12 label/pair incidence); entries of P_d in {-1/8,0,1/8,1/4,1/2};
     spectrum of 3I - P_d = 2^6, 3^6.
 (2) q = (-6 at x; 1 on the 6 labels of alpha,beta,gamma; -3 on the other 8 labels; 2 on Q):
     q^T S_WW q = 0  iff  1^T T 1 = 24 (e = 12), and then A q = 3 q with the stated consequences."""
from fractions import Fraction as Fr
import itertools
import sympy as sp

labels = list(range(14))
cls = lambda u: u // 2
Q = [(u, v) for u, v in itertools.combinations(range(6), 2) if cls(u) != cls(v)]
assert len(Q) == 12
R = ['x'] + labels
W = R + Q
idx = {w: i for i, w in enumerate(W)}
n = len(W)

def A_known(a, b):
    """adjacency among W except Q-Q (which is T, unknown)."""
    if a == b: return 0
    if a == 'x' or b == 'x':
        o = b if a == 'x' else a
        return 1 if isinstance(o, int) else 0
    if isinstance(a, int) and isinstance(b, int):
        return 1 if cls(a) == cls(b) else 0
    if isinstance(a, int) or isinstance(b, int):
        lab, pr = (a, b) if isinstance(a, int) else (b, a)
        return 1 if lab in pr else 0
    return None   # Q-Q

c = Fr(11, 99)
def S_entry(a, b, T=None):
    base = (3 if a == b else 0) + c
    ad = A_known(a, b)
    if ad is None:
        ad = T[(a, b)] if T else 0
    return base - ad

SRR = sp.Matrix([[S_entry(a, b) for b in R] for a in R])
SRQ = sp.Matrix([[S_entry(a, b) for b in Q] for a in R])
# (1) PSD, range, Schur complement
ev = SRR.eigenvals()
print('(1) S_RR eigenvalues:', {sp.nsimplify(k): v for k, v in ev.items()})
pinv = SRR.pinv()
resid = SRR * pinv * SRQ - SRQ
print('    S_RQ in range(S_RR):', resid == sp.zeros(*resid.shape))
schur_part = sp.Matrix(12, 12, lambda i, j: (3 if i == j else 0) + c) - SRQ.T * pinv * SRQ
d = sp.Matrix([[1 if u in pr else 0 for pr in Q] for u in range(6)])
Pd = d.T * (d * d.T).inv() * d
diff = sp.simplify(schur_part - (3 * sp.eye(12) - Pd))
print('    3I + (11/99)J - S_QR S_RR^+ S_RQ == 3I - P_d :', diff == sp.zeros(12, 12))
print('    entries of P_d:', sorted(set(Pd)))
print('    spectrum of 3I - P_d:', (3 * sp.eye(12) - Pd).eigenvals())
# (2) tight case
q = {}
q['x'] = -6
for u in labels:
    q[u] = 1 if u < 6 else -3
for pr in Q:
    q[pr] = 2
qv = sp.Matrix([q[w] for w in W])
# q^T S_WW q as a function of the Q-Q adjacency T: S_QQ = 3I + cJ - T, so q^T S q = (q^T S q)|_{T=0} - q_Q^T T q_Q
S0 = sp.Matrix([[S_entry(a, b) for b in W] for a in W])
val0 = (qv.T * S0 * qv)[0]
print('(2) q^T S_WW q at T = 0:', val0, ' and q_Q^T T q_Q = 4 * 1^T T 1, so q^T S q = 0 iff 1^T T 1 =', val0 / 4)
print('    1^T q =', sum(qv), ' (so S q = 0 gives A q = 3 q)')
