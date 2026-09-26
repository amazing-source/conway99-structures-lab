"""Universal kernel of the 19x19 block-line odd submatrix (S - 3I mod 7): is there a vector in the kernel for every
(tau, iota)?  Print kernel bases (mod 7) for a few random configurations with rank 18."""
import random, sys
sys.path.insert(0, '.')
from blockline_rank import matrix, rand_involution, labels, mate
import sympy as sp
P = 7
def nullspace_mod(M, p=P):
    m, n = len(M), len(M[0])
    A = [[x % p for x in r] for r in M]; pivs = []; r = 0
    for c in range(n):
        piv = next((i for i in range(r, m) if A[i][c]), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]; inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        pivs.append(c); r += 1
    out = []
    for f in [c for c in range(n) if c not in pivs]:
        x = [0] * n; x[f] = 1
        for i, c in enumerate(pivs): x[c] = (-A[i][f]) % p
        out.append(x)
    return out
rng = random.Random(3)
names = ['c'] + [f'e{e}' for e in range(6)] + [f'o{l[0]}{"+-"[l[1]]}' for l in labels]
for t in range(4):
    tau = rand_involution(rng, True); iota = rand_involution(rng, True)
    ns = nullspace_mod(matrix(tau, iota))
    print(f'config {t}: kernel dim {len(ns)}')
    for v in ns:
        print('   ', ' '.join(f'{n}:{x}' for n, x in zip(names, v) if x))
