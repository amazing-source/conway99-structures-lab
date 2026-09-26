"""Independent re-implementation (from the statement of c1_recherche notes/131 par. 7 bis only) of the q:7 kill of the 24
'resistant' step-1 cases of 3K2_or14: a sib orbit whose shadow takes BOTH labels of one block of the NEXT group.

Setting: sink c = 2 full, D-structure three digons; sib pairs (groups) A = {0,4}, B = {1,3}, C = {5,6}; cyclic order
A -> B -> C -> A.  Star(2) labels (b, e), b in K = {0,1,3,4,5,6}, e = +-1; pi flips e.
Facts used (checked separately, see review/q7_3k2_independent.py and NOTE_S29):
 * full block => g- = 20, R_c = K r K^T with r anisotropic of square discriminant, so r = I_2 over F_7 after a change of basis;
 * T_FF - Z = K K^T (mod 7), Z = 3 delta + 5 B^T B + st^T (5I + 4P_pi) st, st = s + 2(1_z - 1_{sigma z});
   so T_ff = 0 gives |k_f|^2 = -Z_ff, and T_fg = Z_fg + k_f . k_g;
 * sink balance: a sib orbit's shadow has 0 labels in its own group and 2 in each other group; two +1, two -1;
 * C2's exact certificates (24 cases): on the PREVIOUS group the shadow has one label per block;
 * on the NEXT group: one label per block (perm type) or both labels of one block m (type (2,0));
 * 8.1 (sib pairs are perfect grids): the sister's shadow occupies the complementary 4 labels outside the pair; T = 0 between
   the two orbits of a sib pair (D-partners); between two disjoint sib pairs the 2 x 2 block of T is a signed permutation.
Claim to re-check: no option of a pair with (2,0) type is compatible (signed permutation for Z + k.k, for SOME k on the
circles) with any option of the NEXT pair; perm-perm regression: 384 options per pair, 2048 compatible couples."""
import itertools
import numpy as np

P7 = 7
K = [0, 1, 3, 4, 5, 6]
GROUPS = {'A': (0, 4), 'B': (1, 3), 'C': (5, 6)}
NEXT = {'A': 'B', 'B': 'C', 'C': 'A'}
PREV = {v: k for k, v in NEXT.items()}


def li(b, e):
    return 2 * K.index(b) + (0 if e == 1 else 1)


W = 5 * np.eye(12, dtype=int)
for v in range(12):
    W[v, v ^ 1] = 4
CIRCLE = {n: [(a, b) for a in range(7) for b in range(7) if (a * a + b * b) % 7 == n] for n in range(1, 7)}
assert all(len(c) == 8 for c in CIRCLE.values())


def reps(g):
    b1, b2 = GROUPS[g]
    return [((b1, 1), (b2, 1)), ((b1, 1), (b2, -1))]


def dvec(rep):
    d = np.zeros(12, int)
    for (b, e) in rep:
        d[li(b, e)] += 1
        d[li(b, -e)] -= 1
    return d


def bsign(r1, r2):
    return sum(e1 * e2 for (b1, e1) in r1 for (b2, e2) in r2 if b1 == b2)


def shadows(g):
    """list of (support list, type) for orbit 0 of pair g; type 'perm' or ('20', m)."""
    pb = GROUPS[PREV[g]]
    nb = GROUPS[NEXT[g]]
    out = []
    for ep in itertools.product((1, -1), repeat=2):
        prevlabs = [li(pb[0], ep[0]), li(pb[1], ep[1])]
        for en in itertools.product((1, -1), repeat=2):
            out.append((prevlabs + [li(nb[0], en[0]), li(nb[1], en[1])], 'perm'))
        for m in nb:
            out.append((prevlabs + [li(m, 1), li(m, -1)], ('20', m)))
    return out


def signed(supp):
    for plus in itertools.combinations(range(4), 2):
        s = np.zeros(12, int)
        for i, v in enumerate(supp):
            s[v] = 1 if i in plus else -1
        yield s


def outside(g):
    return sorted(li(b, e) for b in K if b not in GROUPS[g] for e in (1, -1))


def options(g, with20=True):
    """all admissible (orbit data) for sib pair g: list of dicts with st (2 x 12), Z-diagonal, k-choices."""
    R = reps(g)
    D = [dvec(r) for r in R]
    out = []
    for supp, typ in shadows(g):
        if typ != 'perm' and not with20:
            continue
        comp = [v for v in outside(g) if v not in supp]
        for s0 in signed(supp):
            for s1 in signed(comp):
                st = [s0 + 2 * D[0], s1 + 2 * D[1]]
                z00 = (3 + 5 * bsign(R[0], R[0]) + st[0] @ W @ st[0]) % P7
                z11 = (3 + 5 * bsign(R[1], R[1]) + st[1] @ W @ st[1]) % P7
                z01 = (5 * bsign(R[0], R[1]) + st[0] @ W @ st[1]) % P7
                ks0 = [(0, 0)] if z00 == 0 else CIRCLE[(-z00) % 7]
                ks1 = [(0, 0)] if z11 == 0 else CIRCLE[(-z11) % 7]
                kp = [(k0, k1) for k0 in ks0 for k1 in ks1 if (z01 + k0[0] * k1[0] + k0[1] * k1[1]) % 7 == 0]
                if kp:
                    out.append({'st': st, 'type': typ, 'kpairs': kp, 'rigid': (z00 == 0, z11 == 0)})
    return out


def is_signed_perm(M):
    for p in ((0, 1), (1, 0)):
        if all((M[i][j] in (1, 6)) if j == p[i] else (M[i][j] == 0) for i in range(2) for j in range(2)):
            return True
    return False


def canon_k(o):
    """O(2, F_7) is transitive on each circle and acts on all k's at once: fix k of the first non-rigid orbit of the
    (2,0) option to the first point of its circle (the compatibility test is invariant under this global action)."""
    if o['type'] == 'perm':
        return o['kpairs']
    k0s = [kp for kp in o['kpairs'] if kp[0] == (0, 0) or kp[0] == CIRCLE[(kp[0][0] ** 2 + kp[0][1] ** 2) % 7][0]]
    return k0s


def compatible(g, og, h, oh, canon=True):
    Rg, Rh = reps(g), reps(h)
    Z = [[(5 * bsign(Rg[i], Rh[j]) + og['st'][i] @ W @ oh['st'][j]) % P7 for j in range(2)] for i in range(2)]
    for kg in (canon_k(og) if canon else og['kpairs']):
        for kh in oh['kpairs']:
            M = [[(Z[i][j] + kg[i][0] * kh[j][0] + kg[i][1] * kh[j][1]) % P7 for j in range(2)] for i in range(2)]
            if is_signed_perm(M):
                return True
    return False


def main():
    opts = {g: options(g) for g in GROUPS}
    for g in GROUPS:
        n20 = sum(1 for o in opts[g] if o['type'] != 'perm')
        nperm = len(opts[g]) - n20
        nk = sum(len(o['kpairs']) for o in opts[g] if o['type'] != 'perm')
        print('pair', g, GROUPS[g], ': perm options', nperm, '| (2,0) shadow options', n20, '| (2,0) options with k-pairs', nk,
              flush=True)
    for g in GROUPS:
        h = NEXT[g]
        c20 = 0
        cpp = 0
        for og in opts[g]:
            for oh in opts[h]:
                ok = compatible(g, og, h, oh)
                if og['type'] != 'perm' and ok:
                    c20 += 1
                if og['type'] == 'perm' and oh['type'] == 'perm' and ok:
                    cpp += 1
        print('pair', g, '-> next', h, ': compatible couples with a (2,0) option of', g, ':', c20,
              '| perm-perm compatible couples:', cpp, flush=True)


if __name__ == '__main__':
    main()
