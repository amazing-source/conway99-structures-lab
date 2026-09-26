"""Independent verifier of the certificates written by certify_core.py.  No LP, no floating point: integers only.
For every certificate file:
  1. rebuild the columns (orbit, neighbour pair, exact local solution) from the decomposition, with the exact local
     enumerator (local_system.solve_core, itself cross-checked against doubles_budget3 on 1340 keys), and recompute the
     value table val(column, pair) with SEPARATE code (own frame function);
  2. rebuild the coefficient vector and the bounds of every certificate row FROM ITS LABEL, with separate code;
  3. check  sum_c min(0, w_c) u_c  >  sum_r max(y_r lo_r, y_r hi_r),  w = sum_r y_r A_r, u = branch box;
  4. check that the leaves form a complete branching tree (each internal node splits one orbit over ALL its pairs);
  5. check that the certified representatives cover all labelled decompositions of the mode (own canonical form)."""
import ast
import gzip
import itertools
import json
import os
import re
import sys

sys.path.insert(0, 'scripts')
from local_system import solve_core
from zeros_budget import decompositions

V7 = range(7)


def my_frame(dec, j, i):
    for c in dec:
        L = len(c)
        for t in range(L):
            if {c[t], c[(t + 1) % L]} == {i, j}:
                pos = c.index(i)
                nb = {c[(pos - 1) % L], c[(pos + 1) % L]}
                (k,) = nb - {j}
                return c, k
    raise ValueError


def my_columns(dec):
    cols = []
    for j in V7:
        for i in V7:
            if i == j:
                continue
            cyc, k = my_frame(dec, j, i)
            R = [r for r in V7 if r not in (i, j, k)]
            tris = [set(c) for c in dec if len(c) == 3 and set(c) != set(cyc)]
            for p in itertools.combinations([r for r in V7 if r not in (i, j)], 2):
                for phi, y in solve_core(tris, len(cyc) == 3, i, j, k, R, set(p)):
                    val = {}
                    for a, b in itertools.combinations(V7, 2):
                        Q = frozenset((a, b))
                        if Q == {i, j}:
                            v = 0
                        elif i not in Q and j not in Q:
                            v = y[b if a == k else a] if k in Q else phi[(a, b)]
                        elif j in Q:
                            v = 1 if (a if b == j else b) in p else 0
                        else:
                            v = 2 if (a if b == i else b) == k else 0
                        val[Q] = v
                    cols.append(((j, i), frozenset(p), val))
    return cols


def parse_pair(s):
    return frozenset(ast.literal_eval(s))


def row_from_label(label, cols):
    """Coefficients {column: int} and (lo, hi) rebuilt from the label."""
    kind, rest = label.split(' ', 1)
    coef = {}
    if kind == 'O':
        Z = ast.literal_eval(rest)
        for c, col in enumerate(cols):
            if col[0] == Z:
                coef[c] = 1
        return coef, 1, 1
    if kind == 'Ysym':
        j, rest2 = rest.split(':')
        i, q = map(int, rest2.split('-'))
        j = int(j)
        for c, col in enumerate(cols):
            if col[0] == (j, i) and q in col[1]:
                coef[c] = coef.get(c, 0) + 1
            if col[0] == (j, q) and i in col[1]:
                coef[c] = coef.get(c, 0) - 1
        return coef, 0, 0
    if kind in ('Ysub', 'Ytri'):
        j, rest2 = rest.split(':')
        j = int(j)
        S = ast.literal_eval(rest2)
        O = [q for q in V7 if q != j]
        groups = [S] if kind == 'Ysub' else [S, tuple(q for q in O if q not in S)]
        for G in groups:
            for i, q in itertools.combinations(G, 2):
                for c, col in enumerate(cols):
                    if col[0] == (j, i) and q in col[1]:
                        coef[c] = coef.get(c, 0) + 1
        return coef, 0, (len(S) - 1 if kind == 'Ysub' else 5)
    if kind == 'PPE':
        Q = parse_pair(rest)
        for c, col in enumerate(cols):
            if not (Q & set(col[0])) and col[2][Q] == 2:
                coef[c] = 1
        return coef, 2, 2
    if kind == 'E':
        P, Q = map(parse_pair, rest.split('|'))
        for c, col in enumerate(cols):
            if frozenset(col[0]) == P and col[2][Q]:
                coef[c] = coef.get(c, 0) + col[2][Q]
            if frozenset(col[0]) == Q and col[2][P]:
                coef[c] = coef.get(c, 0) - col[2][P]
        return coef, 0, 0
    if kind == 'P2':
        Q1, Q2 = map(parse_pair, rest.split('|'))
        rhs = 16 if not (Q1 & Q2) else 12
        for c, col in enumerate(cols):
            v = col[2][Q1] * col[2][Q2] + (col[2][Q2] if frozenset(col[0]) == Q1 else 0)
            if v:
                coef[c] = v
        return coef, rhs, rhs
    raise ValueError('unknown row label ' + label)


def check_certificate(cols, cert):
    fixed = {tuple(Z): frozenset(p) for Z, p in cert['fixings']}
    ub = [0 if (col[0] in fixed and col[1] != fixed[col[0]]) else 1 for col in cols]
    w = [0] * len(cols)
    right = 0
    for label, y in cert['rows']:
        coef, lo, hi = row_from_label(label, cols)
        for c, a in coef.items():
            w[c] += y * a
        right += max(y * lo, y * hi)
    left = sum(min(0, w[c]) * ub[c] for c in range(len(cols)))
    return left > right, left, right


def tree_complete(cols, paths, prefix=()):
    if prefix in paths:
        return True
    fixed = dict(prefix)
    # find the orbit branched on at this node: the next element of any path extending the prefix
    ext = [p for p in paths if p[:len(prefix)] == prefix and len(p) > len(prefix)]
    if not ext:
        return False
    Z = ext[0][len(prefix)][0]
    pairs = {col[1] for col in cols if col[0] == Z}
    return all(tree_complete(cols, paths, prefix + ((Z, p),)) for p in pairs)


def canon_any(dec):
    """Own canonical form (independent code).  One triangle: all 144 relabellings sending the triangle onto {0,1,2};
    no triangle: all relabellings sending a cycle of the rarest length onto (0,...,L-1), the rest in any order."""
    tri = [c for c in dec if len(c) == 3]
    maps = []
    if len(tri) == 1:
        T = list(tri[0]); rest = [v for v in V7 if v not in T]
        for a in itertools.permutations(T):
            for b in itertools.permutations(rest):
                maps.append(dict(zip(list(a) + list(b), V7)))
    else:
        lens = [len(c) for c in dec]
        L = min(set(lens), key=lambda l: (lens.count(l), -l))
        for C in [c for c in dec if len(c) == L]:
            for s in range(L):
                for d in (1, -1):
                    seq = [C[(s + d * t) % L] for t in range(L)]
                    rest = [v for v in V7 if v not in seq]
                    for b in itertools.permutations(rest):
                        maps.append(dict(zip(seq + list(b), V7)))
    best = None
    for mp in maps:
        form = tuple(sorted(tuple(sorted(tuple(sorted((mp[c[t]], mp[c[(t + 1) % len(c)]]))) for t in range(len(c))))
                            for c in dec))
        if best is None or form < best:
            best = form
    return best


if __name__ == '__main__':
    mode = sys.argv[1]
    files = sorted(f for f in os.listdir('logs/core_certs') if f.startswith(mode + '_'))
    ok_all, reps = True, []
    for fn in files:
        data = json.load(gzip.open(os.path.join('logs/core_certs', fn), 'rt', encoding='utf-8'))
        dec = [tuple(c) for c in data['decomposition']]
        reps.append((dec, data['labelled']))
        if data['status'] != 'certified':
            print(f'{fn}: status {data["status"]} -> NOT CERTIFIED'); ok_all = False; continue
        cols = my_columns(dec)
        good = True
        for cert in data['certificates']:
            ok, left, right = check_certificate(cols, cert)
            good &= ok
        paths = {tuple((tuple(Z), frozenset(p)) for Z, p in cert['fixings']) for cert in data['certificates']}
        complete = tree_complete(cols, paths)
        print(f'{fn}: {len(data["certificates"])} certificate(s): all valid {good}; branching tree complete {complete}')
        ok_all &= good and complete
    if len(sys.argv) > 2 and sys.argv[2] == 'coverage':
        target = {'one': 1, 'zero': 0}[mode]
        forms = {canon_any(d) for d, _ in reps}
        miss = 0
        total = 0
        for dec in decompositions((1 << 21) - 1, []):
            if sum(1 for c in dec if len(c) == 3) != target:
                continue
            total += 1
            if canon_any([tuple(c) for c in dec]) not in forms:
                miss += 1
        print(f'coverage: {total} labelled decompositions, {miss} not isomorphic to a certified representative')
        ok_all &= miss == 0
    print('ALL VERIFIED' if ok_all else 'NOT ALL VERIFIED')
