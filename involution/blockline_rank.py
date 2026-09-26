"""Rank mod 7 of S - 3I restricted to the odd part (19 coordinates) of the ball of a block line L_c = {x, (c,+), (c,-)},
in C2's involution frame.  Coordinates: d_c ; d_e for the 6 other blocks ; 12 root orbits o_v = {{u,v}, {ubar, vbar}} (v runs
over the 12 labels outside block c, rep = the u-pair).  Entries (S_ab = A(rep a, rep b) - A(rep a, sigma rep b)):
  S_cc = S_ee = -1 ; S_c,o = 1 ; S_e,o_v = +1 if v = (e,+), -1 if v = (e,-), else 0 ;
  S_{o_v,o_w} = [tau_u(v) = w] - [iota(v) = w], iota = mate o sigma_U, an involution (sigma-symmetry of the ball).
tau_u: fixed-point-free involution on the 12 labels; iota: involution (fixed points of iota = S-events sigma_U(v) = vbar).
The rank of this principal submatrix is a lower bound for g- = rank_7(S - 3I)."""
import random, sys
from collections import Counter
P = 7
labels = [(e, s) for e in range(6) for s in (0, 1)]          # 6 other blocks, signs 0 -> +, 1 -> -
mate = {l: (l[0], 1 - l[1]) for l in labels}

def rank_mod(M, p=P):
    M = [[x % p for x in r] for r in M]; r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], -1, p); M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r

def matrix(tau, iota):
    n = 19; M = [[0] * n for _ in range(n)]
    M[0][0] = -1
    for e in range(6): M[1 + e][1 + e] = -1
    for a, v in enumerate(labels):
        o = 7 + a
        M[0][o] = M[o][0] = 1
        e, s = v; sg = 1 if s == 0 else -1
        M[1 + e][o] = M[o][1 + e] = sg
        M[o][7 + labels.index(tau[v])] += 1
        M[o][7 + labels.index(iota[v])] -= 1
    return [[M[i][j] - 3 * (i == j) for j in range(n)] for i in range(n)]

def rand_involution(rng, fpf=True, nfix=None):
    xs = labels[:]; rng.shuffle(xs); f = {}
    k = 0 if fpf else (nfix if nfix is not None else rng.choice([0, 2, 4]))
    fixed, rest = xs[:k], xs[k:]
    for l in fixed: f[l] = l
    for i in range(0, len(rest), 2): f[rest[i]], f[rest[i + 1]] = rest[i + 1], rest[i]
    return f

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
for label, fpf_iota in [('iota fixed-point-free (S = 0 at x on this line)', True), ('iota with 0/2/4 fixed points', False)]:
    cnt = Counter()
    for _ in range(N):
        tau = rand_involution(rng, True)
        iota = rand_involution(rng, fpf_iota)
        cnt[rank_mod(matrix(tau, iota))] += 1
    print(label, '-> rank_7 distribution of the 19x19 block-line submatrix:', dict(sorted(cnt.items())), flush=True)
# the Hamming-like (untwisted) point: tau = iota = mate
print('tau = mate, iota = mate:', rank_mod(matrix(mate, mate)))
print('tau = mate, iota = id  :', rank_mod(matrix(mate, {l: l for l in labels})))
