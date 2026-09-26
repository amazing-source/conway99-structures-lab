"""Independent check of C1's formula (notes/131 par. 7) at a FULL block c (mu_c = tau_c, so S vanishes on orbit-orbit entries of
the block line):  X^T P^+ X == 5 B^T B + st^T (5I + 4P_pi) st  (mod 7),  st = s + 2(1_z - 1_{sigma z}),
where P = S_II - 3I on the 19 block-line coordinates (d_c, 6 other blocks d_e, 12 orbits o_v through c), X = far columns
(far orbit f with rep z: X[d_e,f] = sign of z at block e, X[o_v,f] = [v in Phi(z)] - [v in Phi(sigma z)]), and B^T B is the
Gram of the far block-sign vectors.  X^T P^+ X is computed by exact elimination mod 7 (P has rank 18, kernel z_c, and
z_c^T X = 0).  Random admissible shadows; exact comparison."""
import random, itertools
import numpy as np
P7 = 7
L = list(range(12))                     # labels v = 2e + s, e = 0..5 (the six other blocks), s = 0 (+), 1 (-)
def build_P():
    S = np.zeros((19, 19), int)
    S[0, 0] = -1
    for e in range(6): S[1 + e, 1 + e] = -1
    for v in L:
        o = 7 + v
        S[0, o] = S[o, 0] = 1
        sg = 1 if v % 2 == 0 else -1
        S[1 + v // 2, o] = S[o, 1 + v // 2] = sg
    return (S - 3 * np.eye(19, dtype=int)) % P7
def solve_mod(Pm, b):
    """one solution u of P u = b (mod 7), P symmetric of rank 18; raises if none"""
    n = Pm.shape[0]
    A = [list(map(int, Pm[i])) + [int(b[i]) % P7] for i in range(n)]
    r = 0; piv = []
    for c in range(n):
        p = next((i for i in range(r, n) if A[i][c] % P7), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = pow(A[r][c], P7 - 2, P7)
        A[r] = [x * inv % P7 for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] % P7:
                f = A[i][c]; A[i] = [(a - f * bb) % P7 for a, bb in zip(A[i], A[r])]
        piv.append(c); r += 1
    for i in range(r, n):
        assert A[i][n] % P7 == 0, 'inconsistent'
    u = [0] * n
    for i, c in enumerate(piv): u[c] = A[i][n]
    return np.array(u)
def rand_far(rng):
    # far root z = {a, a'} from two different blocks; shadow: Phi(z), Phi(sigma z) two distinct 2-sets (arbitrary here)
    e1, e2 = rng.sample(range(6), 2)
    a, a2 = 2 * e1 + rng.randint(0, 1), 2 * e2 + rng.randint(0, 1)
    phi1 = rng.sample(L, 2)
    while True:
        phi2 = rng.sample(L, 2)
        if set(phi2) != set(phi1): break
    col = np.zeros(19, int)
    for x in (a, a2):
        col[1 + x // 2] = 1 if x % 2 == 0 else -1
    for v in phi1: col[7 + v] += 1
    for v in phi2: col[7 + v] -= 1
    z1 = np.zeros(12, int); z1[a] += 1; z1[a2] += 1; z1[a ^ 1] -= 1; z1[a2 ^ 1] -= 1
    return col, z1
Pm = build_P()
Ppi = np.zeros((12, 12), int)
for v in L: Ppi[v, v ^ 1] = 1
rng = random.Random(1)
bad = 0
for trial in range(300):
    cols = [rand_far(rng) for _ in range(6)]
    X = np.array([c for c, _ in cols]).T                      # 19 x 6
    U = np.array([solve_mod(Pm, X[:, j]) for j in range(X.shape[1])]).T
    lhs = (X.T @ U) % P7
    Bf = X[1:7, :]
    s = X[7:, :]
    st = (s + 2 * np.array([z for _, z in cols]).T) % P7
    rhs = (5 * Bf.T @ Bf + st.T @ (5 * np.eye(12, dtype=int) + 4 * Ppi) @ st) % P7
    if not (lhs == rhs).all(): bad += 1
print('trials 300, mismatches:', bad)
