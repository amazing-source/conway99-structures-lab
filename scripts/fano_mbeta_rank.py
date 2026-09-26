"""Exact rank of C1's M_beta at the Fano-tight witness (verification, no solver).
M_beta = (C - 3 Ind) restricted to E_{-2} = {beta in Q^21 : sum_{R containing c} beta_R = 0 for every class c},
C[z][R] = c_R(z) (block counts, solved exactly in fano_cellcross_rank), Ind[z][R] = [R = blk(z)].
rank M_beta = 14 - dim ker [Inc ; C - 3 Ind] (Inc = 7 x 21 class/block incidence, rank 7 over Q and F_7)."""
import io
import contextlib
import sys
from fractions import Fraction as Fr
sys.path.insert(0, 'scripts')
with contextlib.redirect_stdout(io.StringIO()):
    import fano_cellcross_rank as F


def rank_mod(rows, p=None):
    A = [[(Fr(v) if p is None else v % p) for v in row] for row in rows]
    r = 0
    ncol = len(A[0])
    for col in range(ncol):
        piv = next((i for i in range(r, len(A)) if A[i][col] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = (1 / A[r][col]) if p is None else pow(A[r][col], p - 2, p)
        A[r] = [(v * inv) if p is None else (v * inv) % p for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][col] != 0:
                f = A[i][col]
                A[i] = [(a - f * b) if p is None else (a - f * b) % p for a, b in zip(A[i], A[r])]
        r += 1
    return r


B = F.BLOCKS
Inc = [[1 if c in R else 0 for R in B] for c in range(7)]
Cm3 = [[F.c[z][R] - 3 * (R == F.blk(z)) for R in B] for z in F.PAIRS]
for p in (None, 7):
    rI = rank_mod(Inc, p)
    rS = rank_mod(Inc + Cm3, p)
    ker = 21 - rS
    print(f'{"Q" if p is None else "F_7"}: rank Inc = {rI}; dim ker [Inc; C - 3 Ind] = {ker}; '
          f'rank M_beta = {(21 - rI) - ker} (out of {21 - rI})')
