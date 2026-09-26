"""Smallest integer certificates x (x^T P x < 0) for the non-PSD compressions of dossier 4, section 6 (brute force, numpy)."""
import itertools
import numpy as np

def pm(Qbv, suppV):
    # order U, W, A, B, V ; supports U={i,j}, W={k,l}, A={i,k}, B={p,q}
    S = {'U': {'i', 'j'}, 'W': {'k', 'l'}, 'A': {'i', 'k'}, 'B': {'p', 'q'}, 'V': set(suppV)}
    Q = {('U', 'W'): 2, ('U', 'A'): 1, ('W', 'A'): 1, ('U', 'B'): 1, ('W', 'B'): 1, ('A', 'B'): 0,
         ('A', 'V'): 2, ('U', 'V'): 0, ('W', 'V'): 0, ('B', 'V'): Qbv}
    n = 'UWABV'
    M = np.zeros((5, 5), dtype=np.int64)
    for a in range(5):
        for b in range(5):
            if a == b:
                M[a, b] = 5
            else:
                q = Q.get((n[a], n[b]), Q.get((n[b], n[a])))
                M[a, b] = -2 * q - len(S[n[a]] & S[n[b]]) + 1
    return M

R = 8
g = np.array(list(itertools.product(range(-R, R + 1), repeat=5)), dtype=np.int64)
g = g[np.abs(g).sum(1) > 0]
l1 = np.abs(g).sum(1)
for sV, qbv in [('jl', 0), ('jp', 0), ('jp', 1), ('jr', 1), ('lp', 0), ('lp', 1), ('lr', 1)]:
    M = pm(qbv, sV)
    v = np.einsum('ni,ij,nj->n', g, M, g)
    neg = np.where(v < 0)[0]
    best = neg[np.lexsort((v[neg], l1[neg]))[0]]
    print(f'supp V = {sV}, Q_BV = {qbv}: x (U,W,A,B,V) = {tuple(int(t) for t in g[best])}, x^T P x = {int(v[best])}')
    if (sV, qbv) == ('jp', 0):
        print(M)
