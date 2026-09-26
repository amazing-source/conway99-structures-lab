"""Independent re-implementation of C1's q:7 kill of 3K2_or14 (c1_recherche notes/131 par. 7), written from the statement only.

Setting (0-based blocks 0..6): the sink c = 2 is FULL (every orbit through 2 pivots at 2), with D-structure three digons
{0,4}, {1,3}, {5,6} at 2; sib pairs A = {0,4}, B = {1,3}, C = {5,6}.  Far orbits = the 30 orbits of the 15 pairs inside
K = {0,1,3,4,5,6}.  Labels of star(2): x = (block b in K, sign e in {+1,-1}).

Facts used (all checked separately):
 * full block => g- = 20 and the far residue R_c = K r K^T with r anisotropic (disc square at g = 20, local sign +1);
 * T_FF - Z = R_c (mod 7) with Z = 3 delta + 5 B^T B + st^T (5I + 4P) st, st = s + 2 delta_f   [review/z_formula_check.py];
 * Z_ff = 0 => k_f = 0 => T_{f,f'} = Z_{f,f'} for every far f' ('rigid' row);
 * step 1 of C1 (pair-level MILP, to be certified): each sib orbit's shadow has exactly one label in each of its 4 blocks;
 * 8.1 (sib pairs are perfect grids): a root disjoint from a sib pair P has exactly one neighbour among the 4 roots of P, a root
   meeting P has none.  Consequences: shadows of the two sib orbits of P partition the 8 labels outside P; T = 0 between the
   sister orbits of P (D-partners); between two disjoint sib pairs T is a signed permutation; a cross orbit meeting P has T = 0
   with both orbits of P, and exactly one +-1 with the orbits of the sib pair it avoids;
 * balance at the sink (components {0,4}, {1,3}, {5,6}): a cross orbit's shadow has 1 label in each group it meets and 2 in
   the third, two +1 and two -1.
Everything below is exact modular arithmetic (integers mod 7)."""
import itertools

P7 = 7
K = [0, 1, 3, 4, 5, 6]
GROUPS = {'A': (0, 4), 'B': (1, 3), 'C': (5, 6)}
LABELS = [(b, e) for b in K for e in (1, -1)]


def W(x, y):                      # (5I + 4P_pi)
    if x[0] != y[0]:
        return 0
    return 5 if x[1] == y[1] else 4


def quad(u, v):
    return sum(cu * cv * W(x, y) for x, cu in u.items() for y, cv in v.items())


def delta(rep):
    d = {}
    for (b, e) in rep:
        d[(b, e)] = d.get((b, e), 0) + 1
        d[(b, -e)] = d.get((b, -e), 0) - 1
    return d


def lin(u, v, k):
    w = dict(u)
    for x, c in v.items():
        w[x] = w.get(x, 0) + k * c
    return {x: c for x, c in w.items() if c}


def bsign(rep1, rep2):
    return sum(e1 * e2 for (b1, e1) in rep1 for (b2, e2) in rep2 if b1 == b2)


def Zval(rep1, sh1, rep2, sh2, same):
    v = 5 * bsign(rep1, rep2) + quad(lin(sh1, delta(rep1), 2), lin(sh2, delta(rep2), 2)) + (3 if same else 0)
    return v % P7


def sib_reps(g):
    b1, b2 = GROUPS[g]
    return [((b1, 1), (b2, 1)), ((b1, 1), (b2, -1))]   # the two orbits of the pair (representatives)


def sib_shadows(g):
    """shadows of a sib orbit: one label per block outside the pair (and outside the sink), two +1 and two -1."""
    others = [b for b in K if b not in GROUPS[g]]
    out = []
    for es in itertools.product((1, -1), repeat=4):
        labs = [(b, e) for b, e in zip(others, es)]
        for plus in itertools.combinations(range(4), 2):
            out.append({labs[i]: (1 if i in plus else -1) for i in range(4)})
    return out


def sister_ok(g, s0, s1):
    reps = sib_reps(g)
    if {(b, -e) for (b, e) in s0} != set(s1):                      # complementary supports (8.1)
        return False
    if Zval(reps[0], s0, reps[0], s0, True) or Zval(reps[1], s1, reps[1], s1, True):
        return False                                              # rigid (automatic, kept as a check)
    return Zval(reps[0], s0, reps[1], s1, False) == 0             # T = 0 between sisters


def perm_ok(gx, sx, gy, sy):
    rx, ry = sib_reps(gx), sib_reps(gy)
    Zm = [[Zval(rx[i], sx[i], ry[j], sy[j], False) for j in range(2)] for i in range(2)]
    for p in ((0, 1), (1, 0)):
        if all((Zm[i][j] in (1, 6)) if j == p[i] else (Zm[i][j] == 0) for i in range(2) for j in range(2)):
            return True
    return False


def group_of(b):
    return next(g for g, bs in GROUPS.items() if b in bs)


def cross_shadows(p):
    g1, g2 = group_of(p[0]), group_of(p[1])
    g3 = next(g for g in GROUPS if g not in (g1, g2))
    L = lambda g: [(b, e) for b in GROUPS[g] for e in (1, -1)]
    out = []
    for a in L(g1):
        for b in L(g2):
            for cc in itertools.combinations(L(g3), 2):
                labs = [a, b, cc[0], cc[1]]
                for plus in itertools.combinations(range(4), 2):
                    out.append({labs[i]: (1 if i in plus else -1) for i in range(4)})
    return out


def main():
    per = {}
    for g in GROUPS:
        opts = sib_shadows(g)
        per[g] = [(s0, s1) for s0 in opts for s1 in opts if sister_ok(g, s0, s1)]
        print('sib pair', g, GROUPS[g], ': admissible sister-shadow pairs', len(per[g]))
    assigns = []
    for a in per['A']:
        for b in per['B']:
            if not perm_ok('A', a, 'B', b):
                continue
            for c in per['C']:
                if perm_ok('B', b, 'C', c) and perm_ok('C', c, 'A', a):
                    assigns.append({'A': a, 'B': b, 'C': c})
    print('complete sib assignments:', len(assigns))
    cross_pairs = [(x, y) for x in K for y in K if x < y and group_of(x) != group_of(y)]
    assert len(cross_pairs) == 12
    shopts = {p: cross_shadows(p) for p in cross_pairs}
    killed = 0
    min_seen = 10 ** 9
    for asg in assigns:
        dead = False
        for p in cross_pairs:
            for rep in (((p[0], 1), (p[1], 1)), ((p[0], 1), (p[1], -1))):
                n_ok = 0
                for sh in shopts[p]:
                    good = True
                    for g in GROUPS:
                        zs = [Zval(sib_reps(g)[t], asg[g][t], rep, sh, False) for t in range(2)]
                        if set(p) & set(GROUPS[g]):
                            good = zs == [0, 0]
                        else:
                            good = sorted(0 if z == 0 else (1 if z in (1, 6) else 9) for z in zs) == [0, 1]
                        if not good:
                            break
                    n_ok += good
                min_seen = min(min_seen, n_ok)
                if n_ok == 0:
                    dead = True
                    break
            if dead:
                break
        killed += dead
    print('assignments killed by a cross orbit with no admissible shadow:', killed, 'of', len(assigns),
          '| min admissible options seen:', min_seen)


if __name__ == '__main__':
    main()
