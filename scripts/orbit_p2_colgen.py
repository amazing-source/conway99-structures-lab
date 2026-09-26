"""TEST (at most 5 minutes; tests the strength of a mathematical ingredient, proves nothing by itself).
C2's ORBIT-LEVEL second-order identity (invariant side) in the pure twisted frame:
    for all orbits U <= W:  sum_Z q_Z(U) q_Z(W) + (q_U(W) + q_W(U))/2 = 12 [U = W] + 4 - Sh(U,W),
q_Z(U) = neighbours of a label of Z among the 2 labels of U; Sh = 0 / 1 / 2 for disjoint / one common block / same
pair.  It is Q_out^2 + Q_out = 12 I + 4 J - Sh (eigenvalues 3, -4 on the 35-dim complement of the block space); the
pair-level (P2) is its compression.
Orbit-level column = (orbit Z, free neighbours p, exact local solution, and for every far pair of value 1 the orbit of
that pair which receives the neighbour; a far double hits both orbits of its pair; the square is the orbit D(Z)).
Rows: one column per orbit; free 6-cycles (symmetric edges, subtour cuts); PPE; orbit coupling q_U(W) = q_W(U) for
orbits of disjoint pairs; the 903 orbit identities.  LP by column generation, phase 1 (artificials of cost 1)."""
import itertools
import sys
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, hstack, identity, csr_matrix

sys.path.insert(0, 'scripts')
import config_model as cm
from zeros_budget import V
from local_system import frame

ORB = [(j, i) for j in V for i in V if i != j]          # orbit (free block j, other block i)
OI = {o: t for t, o in enumerate(ORB)}
PO = [frozenset(o) for o in ORB]


def sh(a, b):
    P, Q = PO[a], PO[b]
    return 2 if P == Q else (1 if P & Q else 0)


def build_rows():
    rows = []
    for Z in ORB:
        rows.append(('O', Z, 1, 1))
    for j in V:
        O = [q for q in V if q != j]
        for i, q in itertools.combinations(O, 2):
            rows.append(('Ysym', (j, i, q), 0, 0))
        for size in (3, 4, 5):
            for S in itertools.combinations(O, size):
                rows.append(('Ysub', (j, S), 0, size - 1))
    for Q in cm.PAIRS:
        rows.append(('PPE', Q, 2, 2))
    for a, b in itertools.combinations(range(42), 2):
        if not (PO[a] & PO[b]):
            rows.append(('C', (a, b), 0, 0))
    for a, b in itertools.combinations_with_replacement(range(42), 2):
        r = 12 * (a == b) + 4 - sh(a, b)
        rows.append(('P2o', (a, b), r, r))
    return rows


def base_columns(dec, RI):
    """Pair-level columns with: fixed sparse entries (O, Ysym, Ysub, PPE), base orbit vector v0, option list."""
    out = []
    for Z, p, vals in cm.columns(dec):
        j, i = Z
        cyc, k, R = frame(dec, j, i)
        O = [q for q in V if q != j]
        fixed = {RI[('O', Z)]: 1}
        for q in p:
            a, b = (i, q) if O.index(i) < O.index(q) else (q, i)
            r = RI[('Ysym', (j, a, b))]
            fixed[r] = fixed.get(r, 0) + (1 if a == i else -1)
        for size in (3, 4, 5):
            for S in itertools.combinations(O, size):
                if i in S:
                    c = sum(1 for q in p if q in S and O.index(i) < O.index(q))
                    if c:
                        fixed[RI[('Ysub', (j, S))]] = c
        v0 = np.zeros(42)
        opts = []
        for q in p:
            v0[OI[(j, q)]] += 1                 # free neighbour: orbit of {j,q} free at j
        v0[OI[(k, i)]] += 2                     # square: orbit of {i,k} free at k (paired via i)
        for Q, val in vals.items():
            if Q & {i, j}:
                continue
            a, b = sorted(Q)
            if val == 2:
                fixed[RI[('PPE', Q)]] = 1
                v0[OI[(a, b)]] += 1; v0[OI[(b, a)]] += 1
            elif val == 1:
                opts.append((OI[(a, b)], OI[(b, a)]))
        out.append((OI[Z], fixed, v0, opts))
    return out


def full_entries(zi, fixed, v, RI):
    ent = dict(fixed)
    nz = np.flatnonzero(v)
    for w in nz:
        if not (PO[w] & PO[zi]):
            key = ('C', (zi, w)) if zi < w else ('C', (w, zi))
            ent[RI[key]] = ent.get(RI[key], 0) + (v[w] if zi < w else -v[w])
    for x in range(len(nz)):
        for y in range(x, len(nz)):
            a, b = nz[x], nz[y]
            r = RI[('P2o', (a, b))]
            ent[r] = ent.get(r, 0) + v[a] * v[b]
    for w in nz:
        if w != zi:
            a, b = (zi, w) if zi < w else (w, zi)
            r = RI[('P2o', (a, b))]
            ent[r] = ent.get(r, 0) + 0.5 * v[w]
    return ent


def run(dec, time_limit=270, add_per_round=4000):
    t0 = time.time()
    rows = build_rows()
    RI = {(r[0], r[1]): t for t, r in enumerate(rows)}
    m = len(rows)
    lo = np.array([r[2] for r in rows], float)
    hi = np.array([r[3] for r in rows], float)
    eq = lo == hi
    base = base_columns(dec, RI)
    print(f'  {len(base)} pair-level columns, {sum(2 ** len(b[3]) for b in base)} orbit-level columns, {m} rows  '
          f'[{time.time() - t0:.0f} s]', flush=True)
    cols, seen = [], set()
    for bi, (zi, fixed, v0, opts) in enumerate(base):
        v = v0.copy()
        for o1, o2 in opts:
            v[o1] += 1
        cols.append(full_entries(zi, fixed, v, RI)); seen.add((bi, (0,) * len(opts)))
    it = 0
    ne = int(eq.sum()); ni = m - ne
    eq_idx = np.flatnonzero(eq); in_idx = np.flatnonzero(~eq)
    pos_eq = {r: t for t, r in enumerate(eq_idx)}
    pos_in = {r: t for t, r in enumerate(in_idx)}
    while True:
        it += 1
        n = len(cols)
        de, re_, ce, di, ri, ci = [], [], [], [], [], []
        for c, ent in enumerate(cols):
            for r, a in ent.items():
                if eq[r]:
                    re_.append(pos_eq[r]); ce.append(c); de.append(a)
                else:
                    ri.append(pos_in[r]); ci.append(c); di.append(a)
        Ae = coo_matrix((de, (re_, ce)), shape=(ne, n)).tocsr()
        Ai = coo_matrix((di, (ri, ci)), shape=(ni, n)).tocsr()
        # variables: x (n), sp (ne), sm (ne), su (ni)
        A_eq = hstack([Ae, identity(ne), -identity(ne), csr_matrix((ne, ni))]).tocsr()
        A_ub = hstack([Ai, csr_matrix((ni, 2 * ne)), -identity(ni)]).tocsr()
        cobj = np.concatenate([np.zeros(n), np.ones(2 * ne + ni)])
        res = linprog(cobj, A_ub=A_ub, b_ub=hi[in_idx], A_eq=A_eq, b_eq=lo[eq_idx],
                      bounds=[(0, 1)] * n + [(0, None)] * (2 * ne + ni), method='highs-ipm',
                      options={'time_limit': 120})
        if res.status != 0:
            return f'master status {res.status} ({res.message[:50]}) at round {it}'
        y = np.zeros(m)
        y[eq_idx] = res.eqlin.marginals
        y[in_idx] = res.ineqlin.marginals
        # pricing: score(c) = y . A_c ; improving if > 0
        Ymat = np.zeros((42, 42))
        yC = np.zeros((42, 42))
        for t, (kind, key, l, h) in enumerate(rows):
            if kind == 'P2o':
                a, b = key
                if a == b:
                    Ymat[a, a] += y[t]
                else:
                    Ymat[a, b] += y[t] / 2; Ymat[b, a] += y[t] / 2
            elif kind == 'C':
                a, b = key
                yC[a, b] = y[t]
        newc = []
        for bi, (zi, fixed, v0, opts) in enumerate(base):
            g = np.zeros(42)
            for w in range(42):
                if w == zi:
                    continue
                if not (PO[w] & PO[zi]):
                    g[w] += yC[zi, w] if zi < w else -yC[w, zi]
                g[w] += Ymat[zi, w]                       # the 1/2 (q_U(W) + q_W(U)) term, 2 x 1/2
            const = sum(y[r] * a for r, a in fixed.items())
            mo = len(opts)
            bits = np.array(list(itertools.product((0, 1), repeat=mo)), dtype=float).reshape(-1, mo)
            base1 = v0.copy()
            D = np.zeros((mo, 42))
            for t, (o1, o2) in enumerate(opts):
                base1[o1] += 1
                D[t, o2] += 1; D[t, o1] -= 1
            Vm = base1 + bits @ D
            sc = const + Vm @ g + np.einsum('ij,jk,ik->i', Vm, Ymat, Vm)
            best = int(np.argmax(sc))
            ch = tuple(int(x) for x in bits[best])
            if sc[best] > 1e-7 and (bi, ch) not in seen:
                newc.append((sc[best], bi, ch))
        newc.sort(reverse=True)
        print(f'  round {it}: {n} columns, phase-1 value {res.fun:.4f}, {len(newc)} improving  '
              f'[{time.time() - t0:.0f} s]', flush=True)
        if not newc:
            return f'CONVERGED: phase-1 value {res.fun:.4f} ({"LP INFEASIBLE" if res.fun > 1e-6 else "LP feasible"})'
        if time.time() - t0 > time_limit:
            return f'stopped at the time limit: phase-1 value {res.fun:.4f} (not converged)'
        for s, bi, ch in newc[:add_per_round]:
            zi, fixed, v0, opts = base[bi]
            v = v0.copy()
            for (o1, o2), c in zip(opts, ch):
                v[o2 if c else o1] += 1
            cols.append(full_entries(zi, fixed, v, RI)); seen.add((bi, ch))


if __name__ == '__main__':
    from certify_core import zero_triangle_classes
    cl = zero_triangle_classes()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    typ, dec, cnt = cl[n]
    print(f'class zero{n} {typ} x{cnt}', flush=True)
    print(run(dec), flush=True)
