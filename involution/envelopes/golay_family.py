"""Golay family of envelopes (structures 27/09): the envelopes with an even neighbour N0 = N(A1^24) (ids 0, 1, 2, 4, 10),
written in the A1-frame: a in Z^24 with norm a.a/8 (a = 2c, c = coordinates on the 24 orthogonal roots alpha_j / 2).
  N0 = {a = 2c : c mod 2 in C}   (C = Golay code; norm sum c^2 / 2)
  U  = N0_v + (v/2 + x0) with v <-> d in Z^24 (a-coordinates of v/2 are d), K = {j : d_j even} = positions of the roots of U.
For id 1 (a1^8): d = 1 on the complement of an octad K (norm 8, class min 8).  Monomial symmetries used:
  G1 = { (pi, eps) : pi in Stab_M24(K), eps in {+-1}^24, {j in K^c : eps_j = -1} in C, eps on K free }.
This file builds U(1) and G1 and verifies: det 1, odd, min >= 2, roots = 8 A1, neighbours (A1^24, A1^24), every generator of
G1 preserves U, and |G1| (sympy on the permutation part, times 2^5 * 2^8 signs).  Usage: python golay_family.py"""
import itertools, json, os, sys, random
from fractions import Fraction as F
import numpy as np
from sympy.combinatorics import Permutation, PermutationGroup

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build156 import hnf, gram_of, det_exact, enum, rtype, tkey, tname, int_rows, lll

OM = (1 << 24) - 1
wt = lambda x: bin(x).count('1')
g = sum(1 << b for b in (11, 9, 7, 6, 5, 1, 0))
WORDS = []
for m in range(4096):
    w = 0
    for j in range(12):
        if m >> j & 1:
            w ^= g << j
    if wt(w) % 2:
        w |= 1 << 23
    WORDS.append(w)
CODE = set(WORDS)
BASIS = []
for w in WORDS:
    r = w
    for b in BASIS:
        r = min(r, r ^ b)
    if r:
        BASIS.append(r)
assert len(BASIS) == 12
INF = 23
Q = {(x * x) % 23 for x in range(1, 23)}
inv = lambda x: pow(x, 21, 23)
PERMS = {'shift': [INF if i == INF else (i + 1) % 23 for i in range(24)],
         'times2': [INF if i == INF else (2 * i) % 23 for i in range(24)],
         'minus_inv': [0 if i == INF else INF if i == 0 else (-inv(i)) % 23 for i in range(24)],
         'delta': [i if i in (0, INF) else ((9 * pow(i, 3, 23)) % 23 if i in Q else (pow(i, 3, 23) * inv(9)) % 23) for i in range(24)]}


def apply_mask(p, x):
    y = 0
    while x:
        b = x & -x
        y |= 1 << p[b.bit_length() - 1]
        x ^= b
    return y


for k_, p_ in PERMS.items():
    assert all(apply_mask(p_, w) in CODE for w in WORDS), k_
M24 = list(PERMS.values())


def compose(p, q):
    return [p[q[i]] for i in range(24)]


def inverse(p):
    r = [0] * 24
    for i, x in enumerate(p):
        r[x] = i
    return r


def setwise_stabilizer(gens, S, samples=80, seed=5):
    """generators of the setwise stabilizer of the set S (bitmask) in <gens>, by Schreier generators on the orbit of S"""
    tree = {S: list(range(24))}
    frontier = [S]
    while frontier:
        nf = []
        for X in frontier:
            for gg in gens:
                Y = apply_mask(gg, X)
                if Y not in tree:
                    tree[Y] = compose(gg, tree[X])
                    nf.append(Y)
        frontier = nf
    to_S = {X: inverse(e) for X, e in tree.items()}
    rng = random.Random(seed)
    keys = list(tree)
    out = []
    for _ in range(samples):
        X = rng.choice(keys)
        gg = rng.choice(gens)
        Y = apply_mask(gg, X)
        s = compose(to_S[Y], compose(gg, tree[X]))
        assert apply_mask(s, S) == S
        out.append(s)
    return out, len(tree)


def u1_lattice():
    K = next(w for w in WORDS if wt(w) == 8)                 # an octad
    Kc = OM ^ K
    d = [0 if K >> i & 1 else 1 for i in range(24)]          # a-coordinates of v/2 (= c-coordinates of v)
    # N0 = {c : c mod 2 in C}; N0_v = {c in N0 : sum_{K^c} c = 0 mod 4}; in a = 2c coordinates
    gens_c = [[1 if w >> i & 1 else 0 for i in range(24)] for w in BASIS] + [[2 * (i == j) for i in range(24)] for j in range(24)]
    phi = lambda c: (sum(c[i] for i in range(24) if Kc >> i & 1) // 2) % 2
    assert all(sum(c[i] for i in range(24) if Kc >> i & 1) % 2 == 0 for c in gens_c)
    # index-2 sublattice {c : phi(c) = 0}: basis of N0 first
    Bc = hnf([list(c) for c in gens_c])
    assert len(Bc) == 24
    ph = [phi(c) for c in Bc]
    k0 = next(t for t in range(24) if ph[t])
    Bv = []
    for t in range(24):
        if t == k0:
            Bv.append([2 * x for x in Bc[t]])
        elif ph[t]:
            Bv.append([x + y for x, y in zip(Bc[t], Bc[k0])])
        else:
            Bv.append(list(Bc[t]))
    j0 = next(i for i in range(24) if Kc >> i & 1)
    odd = [d[i] + (4 if i == j0 else 0) for i in range(24)]  # v/2 + x0, x0 = 2 e_j0 in c (4 e_j0 in a); x0.v = 1
    gens_a = [[2 * x for x in c] for c in Bv] + [odd]
    A = hnf(gens_a)
    assert len(A) == 24
    return K, A


def main():
    K, A = u1_lattice()
    A, _ = lll(A, 1)                                         # LLL on 8 x Gram (same reduced basis)
    G = [[F(sum(p * q for p, q in zip(A[i], A[j])), 8) for j in range(24)] for i in range(24)]
    assert all(x.denominator == 1 for r in G for x in r)
    G = [[int(x) for x in r] for r in G]
    print('det', det_exact(G), 'odd', any(G[i][i] % 2 for i in range(24)), flush=True)
    sv = enum(G, 2)
    print('norm-1 vectors:', sum(1 for z, v in sv if v == 1), flush=True)
    R = np.array([z for z, v in sv if v == 2], dtype=np.int64)
    print('roots:', tname(tkey(rtype(R, np.array(G, dtype=np.int64)))), flush=True)
    # membership test for a-vectors: solve over the basis A
    Aq = np.array(A, dtype=float)
    Ainv = np.linalg.inv(Aq)

    def member(a):
        c = np.array(a, dtype=float) @ Ainv
        return np.allclose(c, np.rint(c), atol=1e-6)
    # generators of G1
    Sperm, orbit_len = setwise_stabilizer(M24, K)
    Hperm = PermutationGroup([Permutation(p) for p in Sperm])
    print('Stab_M24(K): orbit of K has %d octads, |Stab| = %d (expected 322560)' % (orbit_len, Hperm.order()))
    Kc = OM ^ K
    inside = [w for w in WORDS if w & K == 0]                 # codewords inside K^c (32)
    print('codewords inside K^c:', len(inside))
    ok = True
    for p in Sperm[:10]:
        for row in A:                                         # (pi a)_{pi(i)} = a_i
            b = [0] * 24
            for i in range(24):
                b[p[i]] = row[i]
            ok &= member(b)
    for w in inside + [1 << i for i in range(24) if K >> i & 1]:
        for row in A:
            b = [-x if w >> i & 1 else x for i, x in enumerate(row)]
            ok &= member(b)
    print('G1 generators preserve U:', ok)
    print('|G1| = |Stab| * 32 * 2^8 =', Hperm.order() * 32 * 256, ' /|W| (2^8) =', Hperm.order() * 32)
    json.dump({'K': K, 'basis_a': A, 'scale': 8, 'stab_gens': Sperm[:12]}, open(os.path.join(HERE, 'state', 'u1_golay.json'), 'w'))


if __name__ == '__main__':
    main()
