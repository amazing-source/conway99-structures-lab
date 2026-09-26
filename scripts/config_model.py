"""Configuration model of the pure twisted orbit system (C2's frame), where C1's (P2) identities become LINEAR.
Variables: x[Z, p, s] = 1 if orbit Z = (j;i) has free neighbours p at block j and local solution s (exact solutions of
C0-C8, from local_system.solve_core).  Every value val_Z(Q) is then a constant of (Z, p, s):
  local phi / y for Q disjoint from {i,j};  n_r = [r in p] for Q = {j,r};  2[r = k] for Q = {i,r};  0 for Q = {i,j}.
Constraints:
  (O) one choice per orbit;
  (Y) the free graph at each block is a single 6-cycle (edge symmetry, and not two triangles);
  (PPE) each pair doubled by exactly 2 orbits;
  (E) coupling E[P][Q] = E[Q][P] for disjoint pairs;
  (P2) optional: sum_Z val_Z(Q1) val_Z(Q2) = 16 - E[Q1][Q2] (disjoint) or 12 - sum_{Z in Q1} val_Z(Q2) (one
       common block).
TEST tool (runs in seconds): it orients the maths, it proves nothing by itself."""
import itertools
import sys

import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

sys.path.insert(0, 'scripts')
from zeros_budget import V
from local_system import frame, solve_core, P

PAIRS = [frozenset(q) for q in itertools.combinations(V, 2)]
LABELS = []


def columns(dec):
    cols = []            # (Z, pair, vals dict Q -> value)
    for j in V:
        for i in V:
            if i == j:
                continue
            cyc, k, R = frame(dec, j, i)
            tris = [set(c) for c in dec if len(c) == 3 and set(c) != set(cyc)]
            for pair in itertools.combinations([r for r in V if r not in (i, j)], 2):
                for phi, y in solve_core(tris, len(cyc) == 3, i, j, k, R, set(pair)):
                    vals = {}
                    for Q in PAIRS:
                        if Q == frozenset((i, j)):
                            vals[Q] = 0
                        elif not (Q & {i, j}):
                            vals[Q] = y[next(iter(Q - {k}))] if k in Q else phi[P(*Q)]
                        elif j in Q:
                            vals[Q] = int(next(iter(Q - {j})) in pair)
                        else:
                            vals[Q] = 2 if next(iter(Q - {i})) == k else 0
                    cols.append(((j, i), frozenset(pair), vals))
    return cols


def typed_columns(cols):
    """Signed refinement: every double of a column gets a type (the block its two neighbours share).
    (L1) a y-double {k,r} is of type r; (L2) at a block t of R_x, the doubles incident to t are all of type t or all
    of another type.  Returns columns (Z, pair, vals, types) with types = {Q: block}."""
    out = []
    for Z, pair, vals in cols:
        j, i = Z
        cyc, k, R = frame_cache(Z, vals)
        dbl = [Q for Q, v in vals.items() if v == 2 and not (Q & {i, j})]
        for choice in itertools.product(*[sorted(Q) for Q in dbl]):
            types = dict(zip(dbl, choice))
            if any(k in Q and types[Q] == k for Q in dbl):                    # (L1)
                continue
            ok = True
            for t in R:
                inc = [Q for Q in dbl if t in Q]
                a = sum(1 for Q in inc if types[Q] == t)
                if a and a < len(inc):                                       # (L2): mixed at t
                    ok = False
                    break
            if ok:
                out.append((Z, pair, vals, types))
    return out


_FR = {}


def frame_cache(Z, vals):
    return _FR[Z]


def build(dec, p2=True, ppe=True, coupling=True, cols=None, typed=False):
    cols = cols if cols is not None else columns(dec)
    if typed:
        for j in V:
            for i in V:
                if i != j:
                    _FR[(j, i)] = frame(dec, j, i)
        cols = typed_columns([c[:3] for c in cols])
    rows, lo, hi, data = [], [], [], []
    r = 0
    LABELS.clear()
    tag = ['']

    def row(entries, l, h):
        nonlocal r
        LABELS.append(tag[0])
        for c, v in entries:
            if v:
                data.append((r, c, v))
        lo.append(l); hi.append(h); r += 1
    orbits = sorted({c[0] for c in cols})
    by_orbit = {Z: [c for c, col in enumerate(cols) if col[0] == Z] for Z in orbits}
    for Z in orbits:                                                        # (O)
        tag[0] = f'O {Z}'
        row([(c, 1) for c in by_orbit[Z]], 1, 1)
    for j in V:                                                             # (Y)
        O = [q for q in V if q != j]
        for i, q in itertools.combinations(O, 2):
            tag[0] = f'Ysym {j}:{i}-{q}'
            e1 = [(c, 1) for c in by_orbit[(j, i)] if q in cols[c][1]]
            e2 = [(c, -1) for c in by_orbit[(j, q)] if i in cols[c][1]]
            row(e1 + e2, 0, 0)
        for size in (3, 4, 5):                                             # subtour cuts (single 6-cycle)
            for S in itertools.combinations(O, size):
                tag[0] = f'Ysub {j}:{S}'
                ent = []
                for i, q in itertools.combinations(S, 2):
                    ent += [(c, 1) for c in by_orbit[(j, i)] if q in cols[c][1]]
                row(ent, 0, size - 1)
        for A in itertools.combinations(O, 3):
            if O[0] not in A:
                continue
            B = [q for q in O if q not in A]
            tag[0] = f'Ytri {j}:{A}'
            ent = []
            for S in (A, B):
                for i, q in itertools.combinations(S, 2):
                    ent += [(c, 1) for c in by_orbit[(j, i)] if q in cols[c][1]]
            row(ent, 0, 5)
    if ppe:
        for Q in PAIRS:
            tag[0] = f'PPE {sorted(Q)}'
            row([(c, 1) for c, col in enumerate(cols) if not (Q & set(col[0])) and col[2][Q] == 2], 2, 2)
            if typed:
                for a in sorted(Q):
                    tag[0] = f'PPEt {sorted(Q)} type {a}'
                    row([(c, 1) for c, col in enumerate(cols) if not (Q & set(col[0])) and col[2][Q] == 2
                         and col[3][Q] == a], 1, 1)
    if coupling:
        for Pp, Q in itertools.combinations(PAIRS, 2):
            if Pp & Q:
                continue
            tag[0] = f'E {sorted(Pp)}|{sorted(Q)}'
            ent = [(c, col[2][Q]) for c, col in enumerate(cols) if frozenset(col[0]) == Pp]
            ent += [(c, -col[2][Pp]) for c, col in enumerate(cols) if frozenset(col[0]) == Q]
            row(ent, 0, 0)
    if p2:
        for Q1, Q2 in itertools.combinations(PAIRS, 2):
            if len(Q1 & Q2) == 2:
                continue
            rhs = 16 if not (Q1 & Q2) else 12
            tag[0] = f'P2 {sorted(Q1)}|{sorted(Q2)}'
            ent = {}
            for c, col in enumerate(cols):
                vals = col[2]
                v = vals[Q1] * vals[Q2]
                if frozenset(col[0]) == Q1:
                    v += vals[Q2]
                if v:
                    ent[c] = ent.get(c, 0) + v
            row(list(ent.items()), rhs, rhs)
    A = coo_matrix(([d[2] for d in data], ([d[0] for d in data], [d[1] for d in data])), shape=(r, len(cols))).tocsr()
    return A, np.array(lo, float), np.array(hi, float), cols


def lp_feasible(A, lo, hi):
    eq = lo == hi
    from scipy.sparse import vstack
    res = linprog(np.zeros(A.shape[1]), A_ub=vstack([A[~eq], -A[~eq]]), b_ub=np.concatenate([hi[~eq], -lo[~eq]]),
                  A_eq=A[eq], b_eq=lo[eq], bounds=(0, 1), method='highs')
    return res.status == 0


def ip_feasible(A, lo, hi, limit=60):
    res = milp(c=np.zeros(A.shape[1]), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(A.shape[1]),
               bounds=Bounds(0, 1), options={'time_limit': limit})
    return {0: True, 2: False}.get(res.status, None)


if __name__ == '__main__':
    import time
    from ppe_typed_lp import SAMPLES
    t0 = time.time()
    PASCH = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 6), (1, 5, 6), (2, 3, 4, 5), (2, 4, 6)]
    for name, dec in [('PASCH (control, must be infeasible)', PASCH)] + [(f'sample {t}', d) for t, d in enumerate(SAMPLES)]:
        cols = columns(dec)
        A0, lo0, hi0, _ = build(dec, p2=False, cols=cols)
        A1, lo1, hi1, _ = build(dec, p2=True, cols=cols)
        print(f'{name}: {len(cols)} columns; without P2: LP {lp_feasible(A0, lo0, hi0)}, IP {ip_feasible(A0, lo0, hi0, 40)}'
              f' ; WITH P2: LP {lp_feasible(A1, lo1, hi1)}, IP {ip_feasible(A1, lo1, hi1, 60)}  [{time.time() - t0:.0f} s]',
              flush=True)
