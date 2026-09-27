"""id 1 (a1^8) in the Golay frame (golay_family.py): admissible y (norm 7) up to G1 = Stab_M24(K) x signs (structures 27/09).
U(1)_odd = {a in Z^24 : a odd on K^c, even on K, W(a) = {j in K^c : a_j = 3 mod 4} u {j in K : a_j = 2 mod 4} in C,
            sum_{K^c} a = 4 mod 8}  (checked: gives exactly the 4096 norm-3 vectors).
A norm-7 y has sum a^2 = 56; |a_j| <= 4 on K (root filter |alpha.y| <= 2 with alpha = 4 e_j).  Shapes (|P2|,|P4| on K;
|P3|,|P5| on K^c) with |P2| + 4|P4| + 2|P3| + 6|P5| = 10, |P2| even.  Up to the sign part of G1 (free on K, inner codewords on
K^c) a y is its absolute pattern (the admissible sign patterns form one orbit).  Single-y filters (exact, from the trap and
the budget): |chi.y| <= 3 for the characteristic vectors chi of norm 8, |u.y| <= 3 for the norm-3 vectors u.
Step 1 (this file): list the sign-normalised admissible y and their orbits under Stab_M24(K).  Usage: python golay_y.py"""
import itertools, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from golay_family import WORDS, CODE, M24, u1_lattice, setwise_stabilizer, apply_mask, OM, wt
from build156 import lll, enum, inv_exact

T0 = time.time()


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def lattice_data():
    K, A = u1_lattice()
    A, _ = lll(A, 1)
    G = [[sum(p * q for p, q in zip(A[i], A[j])) // 8 for j in range(24)] for i in range(24)]
    Aa = np.array(A, dtype=np.int64)
    sv = enum(G, 3)
    N3 = np.array([z for z, v in sv if v == 3], dtype=np.int64) @ Aa          # a-coordinates, one per +-
    Gi = inv_exact(G)
    dg = [G[k][k] % 2 for k in range(24)]
    c = [int(sum(Gi[p][q] * dg[q] for q in range(24))) % 2 for p in range(24)]
    near = enum(G, 2, cnum=[-x for x in c], cden=2)
    CH = []
    for z, v in near:
        w = [2 * p + q for p, q in zip(z, c)]                                  # coefficients of chi (norm 8)
        CH.append(np.array(w, dtype=np.int64) @ Aa)
    CH = np.array(CH) // 2                                                      # a-coords of chi = (2z+c) A / ... (see check)
    return K, A, N3, CH


def main():
    K, A, N3, CH = lattice_data()
    Kc = OM ^ K
    Kl = [i for i in range(24) if K >> i & 1]
    Kcl = [i for i in range(24) if Kc >> i & 1]
    log('norm-3 vectors (a-coords, one per +-): %d; chi8: %d' % (len(N3), len(CH)))
    # sanity: norms
    assert all(int(v @ v) == 24 for v in N3)
    chn = sorted({int(v @ v) for v in CH})
    log('shadow-root a-norms:', chn, '(expect 16 = norm 2 * 8)')
    inner = [w & Kc for w in WORDS if w & K == 0]            # inner codewords (subsets of K^c)
    assert len(inner) == 32
    coset_of = {}                                             # P2 mask -> sorted list of X masks (subsets of K^c)
    for w in WORDS:
        coset_of.setdefault(w & K, []).append(w & Kc)
    shapes = []
    for p5 in (0, 1):
        for p4 in range(3):
            for p2 in range(0, 9, 2):
                rest = 10 - p2 - 4 * p4 - 6 * p5
                if rest >= 0 and rest % 2 == 0:
                    shapes.append((p2, p4, rest // 2, p5))
    log('shapes (|P2|,|P4|,|P3|,|P5|):', shapes)
    ys = []
    for p2, p4, p3, p5 in shapes:
        for P2 in itertools.combinations(Kl, p2):
            P2m = sum(1 << i for i in P2)
            if P2m not in coset_of:
                continue
            Xs = coset_of[P2m]
            rK = [i for i in Kl if i not in P2]
            for P4 in itertools.combinations(rK, p4):
                for P3 in itertools.combinations(Kcl, p3):
                    rc = [i for i in Kcl if i not in P3]
                    for P5 in itertools.combinations(rc, p5):
                        absv = np.zeros(24, dtype=np.int64)
                        absv[list(P2)] = 2; absv[list(P4)] = 4
                        absv[Kcl] = 1; absv[list(P3)] = 3; absv[list(P5)] = 5
                        X = min(Xs)                                    # canonical sign coset representative
                        a = absv.copy()
                        for j in Kcl:
                            inX = X >> j & 1
                            m = absv[j]
                            # a_j = 3 mod 4 iff j in X:  m=1 -> -1 ; m=3 -> +3 ; m=5 -> -5
                            if m == 1:
                                a[j] = -1 if inX else 1
                            elif m == 3:
                                a[j] = 3 if inX else -3
                            else:
                                a[j] = -5 if inX else 5
                        if int(a[Kcl].sum()) % 8 != 4:
                            continue
                        ys.append(a)
    Y = np.array(ys)
    log('sign-normalised y (membership conditions): %d' % len(Y))
    # CH holds the shadow roots beta = chi/2 in a-coordinates (a-norm 16): |beta.y| <= 3/2  <=>  |a_beta . a_y| <= 12;
    # norm-3 u: |u.y| <= 3  <=>  |a_u . a_y| <= 24.  Chunked int32 products (memory).
    ok = np.ones(len(Y), dtype=bool)
    B32, N32 = CH.T.astype(np.int32), N3.T.astype(np.int32)
    for s0 in range(0, len(Y), 8192):
        Yc = Y[s0:s0 + 8192].astype(np.int32)
        ok[s0:s0 + 8192] = (np.abs(Yc @ B32) <= 12).all(axis=1) & (np.abs(Yc @ N32) <= 24).all(axis=1)
    Y = Y[ok]
    log('after single-y filters (chi8, norm 3): %d' % len(Y))
    np.save(os.path.join(HERE, 'state', 'u1_ynorm.npy'), Y)
    # orbits under Stab_M24(K) (permutations), re-normalising signs
    Sperm, _ = setwise_stabilizer(M24, K, samples=40)
    Sperm = Sperm[:6]
    key_index = {tuple(v): i for i, v in enumerate(Y)}
    inner_set = inner

    def normalise(a):
        a = a.copy()
        a[Kl] = np.abs(a[Kl])
        X = 0
        for j in Kcl:
            if (a[j] % 4) == 3:
                X |= 1 << j
        best = None
        for w in inner_set:
            Xp = X ^ w
            if best is None or Xp < best:
                best, bw = Xp, w
        for j in Kcl:
            if bw >> j & 1:
                a[j] = -a[j]
        return tuple(a)
    parent = list(range(len(Y)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    missing = 0
    for p in Sperm:
        for i, v in enumerate(Y):
            b = np.zeros(24, dtype=np.int64)
            b[p] = v                                              # (pi a)_{pi(i)} = a_i
            j = key_index.get(normalise(b))
            if j is None:
                missing += 1
                continue
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
    reps = sorted({find(i) for i in range(len(Y))})
    log('y1 orbits under G1: %d (maps leaving the set: %d, must be 0)' % (len(reps), missing))
    sizes = {}
    for i in range(len(Y)):
        sizes[find(i)] = sizes.get(find(i), 0) + 1
    shape_of = lambda v: (int((np.abs(v[Kl]) == 2).sum()), int((np.abs(v[Kl]) == 4).sum()), int((np.abs(v[Kcl]) == 3).sum()), int((np.abs(v[Kcl]) == 5).sum()))
    json.dump({'reps': [Y[r].tolist() for r in reps], 'sizes': [sizes[r] for r in reps], 'shapes': [shape_of(Y[r]) for r in reps]},
              open(os.path.join(HERE, 'state', 'u1_y1reps.json'), 'w'))
    cls_of = {r: t for t, r in enumerate(reps)}
    np.save(os.path.join(HERE, 'state', 'u1_ylabels.npy'), np.array([cls_of[find(i)] for i in range(len(Y))], dtype=np.int8))
    from collections import Counter
    log('orbit shapes:', Counter(shape_of(Y[r]) for r in reps))


if __name__ == '__main__':
    main()
