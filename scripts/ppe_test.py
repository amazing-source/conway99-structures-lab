"""TEST (at most 5 minutes; it orients the maths and proves nothing).  Question: does the per-pair exact budget
PPE, together with the coupling e(P,Q), close the pure twisted decompositions whose two triangles intersect?
PPE comes from C2's 16.1 in the pure twisted case: every pair P' is doubled by exactly 2 orbits of pairs disjoint from
P', and receives exactly 16 in total, so it has exactly 6 zeros.
Model per decomposition, a MILP solved with scipy/HiGHS:
- the free 6-cycle at each block (binary edges, 2-regular, not two triangles);
- for each orbit X and each pair Q disjoint from P_X, a one-hot value phi_X(Q) in {0,1,2};
- the local system C0-C8;
- PPE (per-pair doubles = 2), column sums = 16;
- the coupling sum_{X in P} phi_X(Q) = sum_{W in Q} phi_W(P).
Controls:
- the triangle of triangles and Fano must be infeasible;
- the two decompositions with disjoint triangles must be infeasible;
- the same models without PPE and the coupling must be feasible for the survivors of the local run."""
import itertools
import sys
import time

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions, V
from local_system import frame

PAIRS = [frozenset(p) for p in itertools.combinations(V, 2)]
LABELS = []                                                       # label of each row of the last build()
VARS = {}                                                         # variable index -> key, for the last build()


def build(dec, ppe=True, coupling=True, colsum=None, typed=False):
    """typed=True adds the signed refinement (types of doubles): a double of X at {t,s} is of type t or s
    (its two neighbours share the t-vertex or the s-vertex); (L1) a y-double {k,r} is of type r; (L2) at a block t of
    R_x, the doubles incident to t are all of type t or all of another type; (P) each pair {a,b} is doubled by exactly
    one orbit of type a and one of type b."""
    if colsum is None:
        colsum = ppe
    var = {}

    def v(key):
        if key not in var:
            var[key] = len(var)
        return var[key]

    rows = []                                                     # (dict var->coef, lo, hi)
    LABELS.clear()
    cur = ['']

    def add(coefs, lo, hi):
        rows.append((coefs, lo, hi))
        LABELS.append(cur[0])

    # free 6-cycles
    for j in V:
        O = [r for r in V if r != j]
        cur[0] = f'cycle@{j}'
        for a in O:
            add({v(('c', j, frozenset((a, b)))): 1 for b in O if b != a}, 2, 2)
        for A in itertools.combinations(O, 3):
            if O[0] not in A:
                continue
            B = [r for r in O if r not in A]
            coefs = {v(('c', j, frozenset(e))): 1 for e in itertools.combinations(A, 2)}
            coefs.update({v(('c', j, frozenset(e))): 1 for e in itertools.combinations(B, 2)})
            add(coefs, 0, 5)
    tris = [set(c) for c in dec if len(c) == 3]

    def phi(X, Q):                                                # linear form of phi_X(Q)
        return {v(('t', X, Q, 1)): 1, v(('t', X, Q, 2)): 2}

    def n(X, r):
        j, i = X
        return v(('c', j, frozenset((i, r))))

    def lin(*parts):
        out = {}
        for p in parts:
            for key, c in p.items():
                out[key] = out.get(key, 0) + c
        return out

    orbits = [(j, i) for j in V for i in V if i != j]
    for X in orbits:
        j, i = X
        cyc, k, R = frame(dec, j, i)
        disj = [Q for Q in PAIRS if not (Q & {i, j})]
        cur[0] = f'onehot X=({j};{i})'
        for Q in disj:
            add({v(('t', X, Q, 0)): 1, v(('t', X, Q, 1)): 1, v(('t', X, Q, 2)): 1}, 1, 1)
        tag = f'X=({j};{i},k={k})'
        if len(cyc) == 3:
            cur[0] = f'C0 {tag}: n{k}=0'
            add({n(X, k): 1}, 0, 0)                                           # C0
        for t in tris:
            if t == set(cyc):
                continue
            if i in t:
                cur[0] = f'C1 {tag}: phi{sorted(t - {i})}=1'
                add(phi(X, frozenset(t - {i})), 1, 1)                         # C1
            elif j in t and k in t:
                (b,) = t - {j, k}
                cur[0] = f'C3 {tag}: n{k}+n{b}+y{b}=1'
                add(lin({n(X, k): 1, n(X, b): 1}, phi(X, frozenset((k, b)))), 1, 1)      # C3
            elif j in t:
                a, b = sorted(t - {j})
                cur[0] = f'C2 {tag}: n{a}+n{b}+phi{a}{b}=1'
                add(lin({n(X, a): 1, n(X, b): 1}, phi(X, frozenset((a, b)))), 1, 1)      # C2
            elif k in t:
                a, b = sorted(t - {k})
                cur[0] = f'C4 {tag}: y{a}+y{b}+phi{a}{b}=3'
                add(lin(phi(X, frozenset((k, a))), phi(X, frozenset((k, b))), phi(X, frozenset((a, b)))), 3, 3)  # C4
            else:
                cur[0] = f'C5 {tag}: sum over {sorted(t)}=3'
                add(lin(*[phi(X, frozenset(e)) for e in itertools.combinations(sorted(t), 2)]), 3, 3)          # C5
        cur[0] = f'C6 {tag}: n{k}+sum y=2'
        add(lin({n(X, k): 1}, *[phi(X, frozenset((k, r))) for r in R]), 2, 2)                  # C6
        if typed:
            for Q in disj:
                a, b = sorted(Q)
                cur[0] = f'type {tag} Q={sorted(Q)}'
                add({v(('tau', X, Q, a)): 1, v(('tau', X, Q, b)): 1, v(('t', X, Q, 2)): -1}, 0, 0)
            for r in R:
                cur[0] = f'L1 {tag} r={r}'
                add({v(('tau', X, frozenset((k, r)), k)): 1}, 0, 0)
            for t in R:
                inc = [frozenset((t, s)) for s in R if s != t] + [frozenset((k, t))]
                cur[0] = f'L2 {tag} t={t}'
                cA = {v(('tau', X, Q, t)): 1 for Q in inc}
                cA[v(('z', X, t))] = -2
                add(cA, -2, 0)                                                   # A_t <= 2 z
                cB = {v(('tau', X, Q, next(iter(Q - {t})))): 1 for Q in inc}
                cB[v(('z', X, t))] = 2
                add(cB, 0, 2)                                                    # B_t <= 2 (1 - z)
        for t in R:
            cur[0] = f'C8({t}) {tag}'
            add(lin({n(X, t): 1}, phi(X, frozenset((k, t))),
                    *[phi(X, frozenset((t, s))) for s in R if s != t]), 4, 4)                  # C8
    for Q in PAIRS:
        Xs = [X for X in orbits if not (Q & {X[0], X[1]})]
        assert len(Xs) == 20
        cur[0] = f'PPE Q={sorted(Q)}'
        if ppe:
            add({v(('t', X, Q, 2)): 1 for X in Xs}, 2, 2)
        if ppe and typed:
            for a in sorted(Q):
                cur[0] = f'PPEtyped Q={sorted(Q)} type={a}'
                add({v(('tau', X, Q, a)): 1 for X in Xs}, 1, 1)
        cur[0] = f'colsum Q={sorted(Q)}'
        if colsum:
            add(lin(*[phi(X, Q) for X in Xs]), 16, 16)
    if coupling:
        for P, Q in itertools.combinations(PAIRS, 2):
            if P & Q:
                continue
            XP = [X for X in orbits if frozenset(X) == P]
            XQ = [X for X in orbits if frozenset(X) == Q]
            cur[0] = f'coupling P={sorted(P)} Q={sorted(Q)}'
            coefs = lin(*[phi(X, Q) for X in XP])
            for key, c in lin(*[phi(W, P) for W in XQ]).items():
                coefs[key] = coefs.get(key, 0) - c
            add(coefs, 0, 0)
    VARS.clear()
    VARS.update({ix: key for key, ix in var.items()})
    A = lil_matrix((len(rows), len(var)))
    lo, hi = np.zeros(len(rows)), np.zeros(len(rows))
    for r, (coefs, l, h) in enumerate(rows):
        for key, c in coefs.items():
            A[r, key] = c
        lo[r], hi[r] = l, h
    return A.tocsr(), lo, hi, len(var)


def feasible(dec, **kw):
    A, lo, hi, nv = build(dec, **kw)
    res = milp(c=np.zeros(nv), constraints=LinearConstraint(A, lo, hi), integrality=np.ones(nv),
               bounds=Bounds(0, 1), options={'time_limit': 20})
    if res.status == 0:
        return True
    if res.status == 2:
        return False
    return None                                                   # time limit or other


def canon2(dec):
    """Canonical form of a decomposition with two triangles sharing a vertex (16 relabellings)."""
    tris = [set(c) for c in dec if len(c) == 3]
    (vv,) = tris[0] & tris[1]
    ef = sorted(set(V) - tris[0] - tris[1])
    best = None
    for T1, T2 in ((tris[0], tris[1]), (tris[1], tris[0])):
        for a in itertools.permutations(sorted(T1 - {vv})):
            for b in itertools.permutations(sorted(T2 - {vv})):
                for e in itertools.permutations(ef):
                    mp = {vv: 0, a[0]: 1, a[1]: 2, b[0]: 3, b[1]: 4, e[0]: 5, e[1]: 6}
                    form = tuple(sorted(tuple(sorted(tuple(sorted((mp[c[t]], mp[c[(t + 1) % len(c)]])))
                                                     for t in range(len(c)))) for c in dec))
                    if best is None or form < best[0]:
                        best = (form, [tuple(mp[x] for x in c) for c in dec])
    return best


if __name__ == '__main__':
    t0 = time.time()
    PASCH = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 6), (1, 5, 6), (2, 3, 4, 5), (2, 4, 6)]
    FANO = [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]
    print('control: triangle of triangles feasible?', feasible(PASCH), '(expected False)', flush=True)
    print('control: Fano feasible?', feasible(FANO), '(expected False)', flush=True)
    classes = {}
    disjoint_examples = {}
    for dec in decompositions((1 << 21) - 1, []):
        tris = [set(c) for c in dec if len(c) == 3]
        if len(tris) != 2:
            continue
        typ = tuple(sorted(len(c) for c in dec))
        if not (tris[0] & tris[1]):
            disjoint_examples.setdefault(typ, dec)
            continue
        form, rep = canon2(dec)
        classes.setdefault((typ, form), [rep, 0])[1] += 1
    print(f'{len(classes)} classes of decompositions with two intersecting triangles [{time.time() - t0:.0f} s]', flush=True)
    for typ, dec in disjoint_examples.items():
        print(f'control: disjoint triangles {typ} feasible?', feasible(dec), '(expected False)', flush=True)
    for (typ, form), (rep, cnt) in sorted(classes.items(), key=lambda kv: kv[0]):
        base = feasible(rep, ppe=False, coupling=False)
        full = feasible(rep)
        print(f'{typ} x{cnt}: {rep}\n    local only feasible: {base} ; with PPE + coupling feasible: {full}  '
              f'[{time.time() - t0:.0f} s]', flush=True)
