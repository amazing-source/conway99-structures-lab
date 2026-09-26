"""Per-line P-Psi status of C1's 9 stored Fano-tight sigmas (finite exact count, no solver), own code.
For a line L = (X,Y,Z) (C1's FANO order) and its 6 labels, each with one flip bit f(x,L): tau_x on L sends class j
bit s to class k bit s xor f (j, k = the other two classes of L in tuple order), and back.  For every pair y = {x, y'}
of a block of L with third point W: z_X(y) = {x, tau_x(y')} and z_Y(y) = {y', tau_y'(x)} both contain a W-label; their
landings via W (class of sigma_W^{(bit of the W-label)} of the other label) must lie on DIFFERENT lines through W
(P-Psi).  Counts, for each line, the flip settings (out of 64) passing P-Psi on its 12 pairs, and among them those whose
tau-graph T is triangle-free (admissible)."""
import itertools
import os
import pickle
import sys

FANO = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2)]
LINE = {frozenset(p): L for L in FANO for p in itertools.combinations(L, 2)}
third = lambda a, b: next(iter(set(LINE[frozenset((a, b))]) - {a, b}))
base = r'C:\Users\RESEARCHER\Desktop\c1_experimental\data'
names = ['x15_witness'] + [f'x26_fano_s41_{k}' for k in range(8)]


def status(sigma):
    sinv = {c: {w: v for v, w in m.items()} for c, m in sigma.items()}
    out = []
    for L in FANO:
        labs = [(c, s) for c in L for s in (0, 1)]
        ok_psi = ok_adm = 0
        for bits in itertools.product((0, 1), repeat=6):
            f = dict(zip(labs, bits))

            def tau(x, v):                                   # tau_x(v) for v on L, class(v) != class(x)
                j, k = [c for c in L if c != x[0]]
                if v[0] == j:
                    return (k, v[1] ^ f[x])
                return (j, v[1] ^ f[x])

            good = True
            for X, Y in itertools.combinations(L, 2):
                W = third(X, Y)
                for s, t in itertools.product((0, 1), repeat=2):
                    x, yy = (X, s), (Y, t)
                    lands = []
                    for a, b in ((x, yy), (yy, x)):
                        wl = tau(a, b)                       # the W-label of z = {a, tau_a(b)}
                        img = sigma[W][a] if wl[1] == 0 else sinv[W][a]
                        lands.append(img[0])
                    if LINE[frozenset((W, lands[0]))] == LINE[frozenset((W, lands[1]))]:
                        good = False
                        break
                if not good:
                    break
            if not good:
                continue
            ok_psi += 1
            # admissibility: the tau-graph on the 12 pairs of L has no triangle
            pairs = [frozenset(p) for X, Y in itertools.combinations(L, 2)
                     for p in itertools.product([(X, 0), (X, 1)], [(Y, 0), (Y, 1)])]
            adj = {p: set() for p in pairs}
            for p in pairs:
                for a in p:
                    (b,) = p - {a}
                    q = frozenset((a, tau(a, b)))
                    adj[p].add(q)
            tri = any(q2 in adj[q1] for p in pairs for q1 in adj[p] for q2 in adj[p] if q1 != q2)
            ok_adm += not tri
        out.append((ok_psi, ok_adm))
    return out


for n in names:
    W = pickle.load(open(os.path.join(base, n + '.pkl'), 'rb'))
    st = status(W['sigma'])
    alive = sum(1 for a, b in st if b > 0)
    print(f'{n}: per line (P-Psi settings / admissible) = ' + ' '.join(f'{a}/{b}' for a, b in st)
          + f'  -> lines alive: {alive}/7')
