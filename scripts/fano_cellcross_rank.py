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
