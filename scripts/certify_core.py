"""Certified computation (authorised by the user on 26/09, capped at about 1 h): the pure twisted core, i.e. the
decompositions with 1 triangle (15 classes) and without triangle, in the configuration model with C1's (P2) and the
subtour cuts.  ONLY EXACT CERTIFICATES COUNT.

For each class (and each branch when branching is needed) the program finds integer multipliers y, one per row of the
model, such that
    sum_c min(0, w_c) * u_c  >  sum_r max(y_r lo_r, y_r hi_r),    with  w = A^T y,
where x_c in [0, u_c] (u_c = 0 for columns removed by the branch), lo_r <= (A x)_r <= hi_r are the rows.  For any
feasible x, sum_r y_r (A x)_r = w^T x lies between the two sides, which is absurd.  The inequality is checked in exact
integer arithmetic.  Certificates are saved (logs/core_certs/*.json.gz) with the row LABELS, and an independent verifier
(verify_core_certs.py) rebuilds every row from its label with separate code.

Branching: if the LP is feasible, choose the orbit whose LP mass is most spread over its neighbour pairs, and split on
its pair (every pair with local solutions is a branch; together they cover all cases)."""
import gzip
import itertools
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import vstack, hstack, identity, csr_matrix

sys.path.insert(0, 'scripts')
import config_model as cm
from zeros_budget import decompositions, V

SCALE = 10 ** 9
OUT = 'logs/core_certs'


def lp_primal(A, lo, hi, ub):
    eq = lo == hi
    res = linprog(np.zeros(A.shape[1]), A_ub=vstack([A[~eq], -A[~eq]]), b_ub=np.concatenate([hi[~eq], -lo[~eq]]),
                  A_eq=A[eq], b_eq=lo[eq], bounds=list(zip(np.zeros(len(ub)), ub)), method='highs',
                  options={'time_limit': 300})
    return res


def farkas(A, lo, hi, ub):
    """LP for a certificate: variables yp, yn >= 0 (rows), t >= 0 (columns);
    t_c + (A^T (yp - yn))_c >= 0 ;  - sum_c u_c t_c - sum_r (hi_r yp_r - lo_r yn_r) >= 1 ; minimise L1."""
    m, n = A.shape
    AT = A.T.tocsr()
    # constraint rows: -(A^T yp - A^T yn + t) <= 0   (n rows),   sum u t + sum (hi yp - lo yn) <= -1   (1 row)
    top = hstack([-AT, AT, -identity(n)])
    bot = csr_matrix(np.concatenate([hi, -lo, ub]).reshape(1, -1))
    res = linprog(np.ones(2 * m + n), A_ub=vstack([top, bot]), b_ub=np.concatenate([np.zeros(n), [-1.0]]),
                  bounds=(0, None), method='highs', options={'time_limit': 300})
    if res.status != 0:
        return None
    y = res.x[:m] - res.x[m:2 * m]
    return y


def exact_check(A, lo, hi, ub, y_int):
    """Exact integer check of the certificate inequality.  A, lo, hi, ub integer; y_int list of Python ints."""
    A = A.tocsc()
    w = [0] * A.shape[1]
    Ac = A.tocoo()
    for r, c, v in zip(Ac.row.tolist(), Ac.col.tolist(), Ac.data.tolist()):
        w[c] += y_int[r] * int(v)
    left = sum(min(0, w[c]) * int(ub[c]) for c in range(len(w)))
    right = sum(max(y * int(l), y * int(h)) for y, l, h in zip(y_int, lo.tolist(), hi.tolist()))
    return left > right, left, right


def certify_node(dec, cols, fixings, labels_cache):
    A, lo, hi, _ = cm.build(dec, p2=True, cols=cols)
    labels = list(cm.LABELS)
    ub = np.ones(A.shape[1])
    for Z, p in fixings:
        for c, col in enumerate(cols):
            if col[0] == Z and col[1] != p:
                ub[c] = 0
    prim = lp_primal(A, lo, hi, ub)
    if prim.status == 0:
        return 'feasible', prim.x, None
    if prim.status != 2:
        return 'undecided', None, None
    y = farkas(A, lo, hi, ub)
    if y is None:
        return 'nocert', None, None
    for scale in (SCALE, SCALE * 1000):
        y_int = [int(round(v * scale)) for v in y]
        ok, left, right = exact_check(A, lo, hi, ub, y_int)
        if ok:
            cert = {'fixings': [[list(Z), sorted(p)] for Z, p in fixings],
                    'rows': [[labels[r], y_int[r]] for r in range(len(y_int)) if y_int[r] != 0],
                    'left': left, 'right': right}
            return 'certified', None, cert
    return 'rounding_failed', None, None


def branch_orbit(cols, x, fixings):
    fixed = {Z for Z, p in fixings}
    mass = {}
    for c, col in enumerate(cols):
        if col[0] in fixed:
            continue
        mass.setdefault(col[0], {}).setdefault(col[1], 0.0)
        mass[col[0]][col[1]] += x[c]
    best = min(mass, key=lambda Z: max(mass[Z].values()))
    pairs = sorted({col[1] for col in cols if col[0] == best}, key=lambda p: -mass[best].get(p, 0))
    return best, pairs


def certify_class(name, dec, deadline, max_depth=3):
    cols = cm.columns(dec)
    certs, stack, status = [], [[]], 'certified'
    nodes = 0
    while stack:
        if time.time() > deadline:
            return 'timeout', certs, nodes
        fixings = stack.pop()
        res, x, cert = certify_node(dec, cols, fixings, None)
        nodes += 1
        if res == 'certified':
            certs.append(cert)
        elif res == 'feasible':
            if len(fixings) >= max_depth:
                return 'depth_exhausted', certs, nodes
            Z, pairs = branch_orbit(cols, x, fixings)
            for p in pairs:
                stack.append(fixings + [(Z, p)])
        else:
            return res, certs, nodes
    return status, certs, nodes


def one_triangle_classes():
    from one_triangle_classes import canon1
    classes = {}
    for dec in decompositions((1 << 21) - 1, []):
        if sum(1 for c in dec if len(c) == 3) != 1:
            continue
        form, rep = canon1(dec)
        classes.setdefault((tuple(sorted(len(c) for c in dec)), form), [rep, 0])[1] += 1
    return [(typ, rep, cnt) for (typ, form), (rep, cnt) in sorted(classes.items())]


def canon0(dec):
    """Canonical form of a triangle-free decomposition: a cycle of the rarest length is mapped onto (0,1,...,L-1)
    in all 2L ways, the remaining blocks in all orders."""
    lens = [len(c) for c in dec]
    L = min(set(lens), key=lambda l: (lens.count(l), -l))
    best = None
    for C in [c for c in dec if len(c) == L]:
        for s in range(L):
            for rev in (False, True):
                seq = [C[(s + t) % L] for t in range(L)]
                if rev:
                    seq = seq[::-1]
                rest = [v for v in V if v not in seq]
                for perm in itertools.permutations(rest):
                    mp = dict(zip(seq + list(perm), range(7)))
                    form = tuple(sorted(tuple(sorted(tuple(sorted((mp[c[t]], mp[c[(t + 1) % len(c)]])))
                                                     for t in range(len(c)))) for c in dec))
                    if best is None or form < best[0]:
                        best = (form, [tuple(mp[x] for x in c) for c in dec])
    return best


def zero_triangle_classes():
    classes = {}
    for dec in decompositions((1 << 21) - 1, []):
        if any(len(c) == 3 for c in dec):
            continue
        form, rep = canon0(dec)
        classes.setdefault((tuple(sorted(len(c) for c in dec)), form), [rep, 0])[1] += 1
    return [(typ, rep, cnt) for (typ, form), (rep, cnt) in sorted(classes.items())]


if __name__ == '__main__':
    mode = sys.argv[1]
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 3000
    only = int(sys.argv[3]) if len(sys.argv) > 3 else None
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    cl = one_triangle_classes() if mode == 'one' else zero_triangle_classes()
    print(f'{mode}: {len(cl)} classes ({sum(c for _, _, c in cl)} labelled)  [{time.time() - t0:.0f} s]', flush=True)
    for n, (typ, rep, cnt) in enumerate(cl):
        if only is not None and n != only:
            continue
        deadline = t0 + budget
        if time.time() > deadline:
            print('BUDGET REACHED'); break
        st, certs, nodes = certify_class(f'{mode}{n}', rep, deadline)
        fn = f'{OUT}/{mode}_{n:02d}.json.gz'
        with gzip.open(fn, 'wt', encoding='utf-8') as f:
            json.dump({'mode': mode, 'index': n, 'type': typ, 'labelled': cnt, 'decomposition': [list(c) for c in rep],
                       'status': st, 'certificates': certs}, f)
        print(f'class {mode}{n} {typ} x{cnt}: {st.upper()}, {len(certs)} certificate(s), {nodes} node(s)  '
              f'[{time.time() - t0:.0f} s]', flush=True)
