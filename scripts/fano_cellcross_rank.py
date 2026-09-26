"""Independent recomputation (own code, exact rationals) of C1's rank obstruction for the cellcross identity at the
Fano-tight witness (tau, sigma read from data/x19k_witness_cells.json, the same witness as data/x15_witness.pkl).
  1. Fano lines read from tau (tau_u pairs the other classes along the lines through class(u)); checked to be a Fano
     plane, identical for both labels of every class.
  2. Specials of every far pair z from tau, sigma.
  3. Block counts c_R(z) at a Fano-tight vertex: specials in blocks meeting blk(z); generic counts g_z(R) on the 10
     blocks disjoint from blk(z) solved from
       tight family counts   sum_{R in L} c_R(z) = 2 (blk(z) in L) / 1 (one common class) / 3 (disjoint), 7 lines L,
       slot equations        sum_{R containing c} c_R(z) = 4 for every class c outside blk(z);
     the solution must be unique, integral, >= 0.  Controls: sum_R c_R(z) = 12, block balance N(R,R') = N(R',R),
     label balance, no internal edge (c_{blk(z)}(z) = 0).
  4. For every block R: cellcross  sum_z delta_R(z) c_R'(z) = - sum_{y in R} eps(y) c_R'(y)  (one equation per R'),
     eps = +1 on the cell {s = t}, -1 on {s != t}; unknowns delta_R(z) for blk(z) cap R empty and c_R(z) > 0; the
     z with |blk(z) cap R| = 1 have delta fixed by their specials (moved to the right).  Exact rank of A and [A | b]."""
import itertools
import json
import sys
from collections import Counter
from fractions import Fraction as Fr

base = sys.argv[1] if len(sys.argv) > 1 else r'C:\Users\RESEARCHER\Desktop\c1_experimental\data\x19k_witness_cells'
D = json.load(open(base + '.json'))
lab = lambda s: tuple(int(t) for t in s.split('.'))
mate = lambda u: (u[0], 1 - u[1])
PAIRS = [frozenset({(c, s), (d, t)}) for c, d in itertools.combinations(range(7), 2) for s in (0, 1) for t in (0, 1)]
BLOCKS = [frozenset(b) for b in itertools.combinations(range(7), 2)]
blk = lambda z: frozenset(l[0] for l in z)
eps = lambda w: 1 if len({l[1] for l in w}) == 1 else -1
bn = lambda R: ''.join(map(str, sorted(R)))

tau = {}
for u, lst in D['tau'].items():
    m = {}
    for v, w in lst:
        m[lab(v)] = lab(w); m[lab(w)] = lab(v)
    tau[lab(u)] = m
sigma = {int(c): {lab(v): lab(w) for v, w in lst} for c, lst in D['sigma'].items()}
sinv = {c: {w: v for v, w in m.items()} for c, m in sigma.items()}

# 1. Fano lines
lines = set()
consistent = True
for c in range(7):
    per_label = []
    for s in (0, 1):
        u = (c, s)
        cls = frozenset(frozenset({c, v[0], w[0]}) for v, w in tau[u].items())
        per_label.append(cls)
    consistent &= per_label[0] == per_label[1]
    lines |= per_label[0]
fano = len(lines) == 7 and all(len(L) == 3 for L in lines) and \
    all(sum(1 for L in lines if {a, b} <= L) == 1 for a, b in itertools.combinations(range(7), 2))
print(f'Fano lines from tau: {sorted(sorted(L) for L in lines)}; consistent for both labels: {consistent}; '
      f'Fano plane: {fano}')

# 2. specials
spec = {}
for z in PAIRS:
    e, f = sorted(z)
    S = [frozenset({e, tau[e][f]}), frozenset({f, tau[f][e]})]
    for u, v in ((e, f), (f, e)):
        w = sigma[u[0]][v] if u[1] == 0 else sinv[u[0]][v]
        S.append(frozenset({mate(u), w}))
    spec[z] = S


# 3. block counts
def solve_exact(M, rhs):
    """Unique solution of M x = rhs over Q (M: list of rows), or None."""
    n = len(M[0])
    A = [[Fr(v) for v in row] + [Fr(r)] for row, r in zip(M, rhs)]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        A[r] = [v / A[r][c] for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                A[i] = [a - A[i][c] * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    if any(all(v == 0 for v in row[:n]) and row[n] != 0 for row in A):
        return 'inconsistent'
    if r < n:
        return 'not unique'
    x = [Fr(0)] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    return x


c = {}
bad_solve = 0
for z in PAIRS:
    P = blk(z)
    cz = {R: sum(1 for w in spec[z] if blk(w) == R) for R in BLOCKS if R & P}
    free = [R for R in BLOCKS if not (R & P)]
    M, rhs = [], []
    for L in lines:
        k = len(L & P)
        target = 2 if k == 2 else (1 if k == 1 else 3)
        inL = [R for R in BLOCKS if R <= L]
        M.append([1 if R in inL else 0 for R in free]); rhs.append(target - sum(cz.get(R, 0) for R in inL if R & P))
    for cl in range(7):
        if cl in P:
            continue
        M.append([1 if cl in R else 0 for R in free])
        rhs.append(4 - sum(v for R, v in cz.items() if cl in R))
    x = solve_exact(M, rhs)
    if not isinstance(x, list) or any(v.denominator != 1 or v < 0 for v in x):
        bad_solve += 1
        continue
    for R, v in zip(free, x):
        cz[R] = int(v)
    c[z] = cz
print(f'block counts solved uniquely, integral, >= 0 for {len(c)}/84 pairs (failures {bad_solve})')
tot = Counter(sum(c[z].values()) for z in c)
internal = sum(c[z][blk(z)] for z in c)
bal = all(sum(c[z][Rp] for z in PAIRS if blk(z) == R) == sum(c[z][R] for z in PAIRS if blk(z) == Rp)
          for R, Rp in itertools.combinations(BLOCKS, 2))
print(f'controls: sum_R c_R(z) distribution {dict(tot)}; internal edges {internal}; block balance {bal}')

# 4. cellcross systems, exact ranks
def rank(rows):
    A = [[Fr(v) for v in row] for row in rows]
    r = 0
    for col in range(len(A[0])):
        p = next((i for i in range(r, len(A)) if A[i][col] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        for i in range(len(A)):
            if i != r and A[i][col] != 0:
                f = A[i][col] / A[r][col]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
    return r


inconsistent = []
for R in BLOCKS:
    free = [z for z in PAIRS if not (blk(z) & R) and c[z][R] > 0]
    one = [z for z in PAIRS if len(blk(z) & R) == 1]
    own = [z for z in PAIRS if blk(z) == R]
    A, b = [], []
    for Rp in BLOCKS:
        A.append([c[z][Rp] for z in free])
        rhs = -sum(eps(y) * c[y][Rp] for y in own)
        for z in one:
            ds = sum(eps(w) for w in spec[z] if blk(w) == R)
            rhs -= ds * c[z][Rp]
        for z in own:                              # own block: specials in R (none if no internal edge)
            ds = sum(eps(w) for w in spec[z] if blk(w) == R)
            rhs -= ds * c[z][Rp]
        b.append(rhs)
    rA = rank(A)
    rAb = rank([row + [bb] for row, bb in zip(A, b)])
    if rAb > rA:
        inconsistent.append(bn(R))
    print(f'block {bn(R)}: {len(free)} unknowns, rank A = {rA}, rank [A|b] = {rAb}')
print(f'blocks whose cellcross system is inconsistent even without bounds: {inconsistent}')


# 5. universal kernel (C1's parametrisation of the permitted landings), recomputed with the same own solver
def third(p, q):
    (L,) = [L for L in lines if {p, q} <= L]
    (m,) = L - {p, q}
    return m


def counts_abstract(P, Q, a, b):
    """Block counts of a pair of block {P,Q} whose sigma-specials land in {P,a} and {Q,b} (tau-specials in PM, QM)."""
    M = third(P, Q)
    PQ = frozenset((P, Q))
    sp = Counter([frozenset((P, M)), frozenset((Q, M)), frozenset((P, a)), frozenset((Q, b))])
    cz = {R: sp[R] for R in BLOCKS if R & PQ}
    free = [R for R in BLOCKS if not (R & PQ)]
    Mx, rhs = [], []
    for L in lines:
        k = len(L & PQ)
        target = 2 if k == 2 else (1 if k == 1 else 3)
        inL = [R for R in BLOCKS if R <= L]
        Mx.append([1 if R in inL else 0 for R in free]); rhs.append(target - sum(cz.get(R, 0) for R in inL if R & PQ))
    for cl in range(7):
        if cl in PQ:
            continue
        Mx.append([1 if cl in R else 0 for R in free]); rhs.append(4 - sum(v for R, v in cz.items() if cl in R))
    x = solve_exact(Mx, rhs)
    if not isinstance(x, list) or any(v.denominator != 1 or v < 0 for v in x):
        return None
    for R, v in zip(free, x):
        cz[R] = int(v)
    return [cz[R] for R in BLOCKS]


def nullspace_left(vectors, n=21):
    """Basis of {phi in Q^n : <v, phi> = 0 for all v}."""
    rows = [[Fr(x) for x in v] for v in vectors]
    piv, r = [], 0
    for col in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][col] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        rows[r] = [x / rows[r][col] for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col] != 0:
                rows[i] = [a - rows[i][col] * bb for a, bb in zip(rows[i], rows[r])]
        piv.append(col); r += 1
    basis = []
    for fc in [cc for cc in range(n) if cc not in piv]:
        phi = [Fr(0)] * n
        phi[fc] = Fr(1)
        for i, pc in enumerate(piv):
            phi[pc] = -rows[i][fc]
        basis.append(phi)
    return basis


bad_abs, witness_in_family = 0, True
univ_dims, witness_orth = {}, {}
fam = {}
for P, Q in itertools.combinations(range(7), 2):
    M = third(P, Q)
    C4 = [x for x in range(7) if x not in (P, Q, M)]
    fam[frozenset((P, Q))] = []
    for a in C4:
        for b in {third(P, a), third(Q, a)}:
            v = counts_abstract(P, Q, a, b)
            if v is None:
                bad_abs += 1
            else:
                fam[frozenset((P, Q))].append(v)
for z in PAIRS:                                    # the witness's own counts must be one of the 8 permitted vectors
    witness_in_family &= [c[z][R] for R in BLOCKS] in fam[blk(z)]
print(f'abstract counts: {sum(len(v) for v in fam.values())} (block pair, landing) vectors, failures {bad_abs}; '
      f'every witness count vector is a permitted one: {witness_in_family}')
for R in BLOCKS:
    vecs = [v for PQ, vs in fam.items() if not (PQ & R) for v in vs]
    K = nullspace_left(vecs)
    univ_dims[bn(R)] = len(K)
    own = [z for z in PAIRS if blk(z) == R]
    one = [z for z in PAIRS if len(blk(z) & R) == 1]
    bR = []
    for Rp in BLOCKS:
        rhs = -sum(eps(y) * c[y][Rp] for y in own)
        for z in one + own:
            rhs -= sum(eps(w) for w in spec[z] if blk(w) == R) * c[z][Rp]
        bR.append(rhs)
    witness_orth[bn(R)] = all(sum(p * q for p, q in zip(phi, bR)) == 0 for phi in K)
print(f'universal kernel dimension per block: {univ_dims}')
print(f'witness right-hand side orthogonal to the universal kernel, per block: '
      f'{all(witness_orth.values())} {[k for k, v in witness_orth.items() if not v]}')


# 6. C1's functional Psi_R (26/09, "PISTE a verifier"), checked with own data
def Psi(R):
    P, Q = sorted(R)
    M = third(P, Q)
    C4 = [x for x in range(7) if x not in (P, Q, M)]
    MM = [frozenset(pr) for pr in itertools.combinations(C4, 2) if third(*pr) == M]
    assert len(MM) == 2 and MM[0] | MM[1] == set(C4)
    s = {x: (1 if x in MM[0] else -1) for x in C4}   # sign choice: the first M_M pair gets +1 (global sign free)
    psi = {B: Fr(0) for B in BLOCKS}
    for x in C4:
        psi[frozenset((M, x))] += 9 * s[x]
        psi[frozenset((P, x))] += s[x]
        psi[frozenset((Q, x))] += s[x]
    psi[MM[0]] += 4
    psi[MM[1]] -= 4
    return [psi[B] for B in BLOCKS], s, M


RI = {B: t for t, B in enumerate(BLOCKS)}
eff_dims, psi_ok, psi_off = {}, True, Counter()
for R in BLOCKS:
    vecs = [v for PQ, vs in fam.items() if not (PQ & R) for v in vs if v[RI[R]] > 0]
    eff_dims[bn(R)] = len(nullspace_left(vecs))
    psi, s, M = Psi(R)
    for PQ, vs in fam.items():
        if PQ & R:
            continue
        for v in vs:
            val = sum(p * q for p, q in zip(v, psi))
            if v[RI[R]] > 0:
                psi_ok &= val == 0
            else:
                psi_off[val] += 1
print(f'effective kernel (only landings with c_R > 0) dimension per block: {set(eff_dims.values())} {eff_dims}')
print(f'<c(z;a,b), Psi_R> = 0 whenever c_R(z;a,b) > 0, all 21 blocks x all landings: {psi_ok}; values when c_R = 0: '
      f'{dict(psi_off)}')


def landing(z, via_class):
    """Class where the sigma-special of z through class via_class lands."""
    (e,) = [l for l in z if l[0] == via_class]
    (f,) = [l for l in z if l[0] != via_class]
    w = sigma[e[0]][f] if e[1] == 0 else sinv[e[0]][f]
    return w[0]


vals, closed_ok, I = {}, True, {}
for R in BLOCKS:
    psi, s, M = Psi(R)
    P, Q = sorted(R)
    own = [z for z in PAIRS if blk(z) == R]
    one = [z for z in PAIRS if len(blk(z) & R) == 1]
    bR = []
    for Rp in BLOCKS:
        rhs = -sum(eps(y) * c[y][Rp] for y in own)
        for z in one + own:
            rhs -= sum(eps(w) for w in spec[z] if blk(w) == R) * c[z][Rp]
        bR.append(rhs)
    val = sum(p * q for p, q in zip(bR, psi))
    vals[bn(R)] = val
    # closed form: -6 sum_{y in R} eps(y) [theta(z_P(y)) + theta(z_Q(y))]
    tot = 0
    for y in own:
        for X in (P, Q):
            (e,) = [l for l in y if l[0] == X]
            (f,) = [l for l in y if l[0] != X]
            zX = frozenset({e, tau[e][f]})                  # tau-special of y through X: block {X, M}
            assert blk(zX) == frozenset((X, M))
            a = landing(zX, X)
            bM = landing(zX, M)
            if bM == third(M, a):
                eta = 1
            elif bM == third(X, a):
                eta = -1
            else:
                raise AssertionError('landing outside the permitted set')
            tot += eps(y) * s[a] * eta
    closed_ok &= val == -6 * tot
    I[bn(R)] = val / 24
print(f'<rhs_R, Psi_R> at the witness: {dict((k, int(v)) for k, v in vals.items() if v)} (zero elsewhere)')
print(f'blocks with <rhs_R, Psi_R> != 0: {sorted(k for k, v in vals.items() if v)}')
print(f'closed form -6 sum eps(y)[theta(z_P(y)) + theta(z_Q(y))] equal to <rhs_R, Psi_R> on all 21 blocks: {closed_ok}')
print(f'I_R = <rhs_R, Psi_R>/24: {dict((k, str(v)) for k, v in I.items() if v)}')
