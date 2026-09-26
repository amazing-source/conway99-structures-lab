"""Independent verifier of the ORBIT-LEVEL certificates written by orbit_p2_colgen.py (logs/orbit_certs/zero_NN.json.gz).
No LP, no floating point, no code shared with orbit_p2_colgen.py or config_model.py: integers only.
  1. Pair-level columns (orbit, free-neighbour pair, exact local solution, value of every block pair) come from
     verify_core_certs.my_columns (own frame function; exact local enumerator).  Each one is expanded HERE into
     orbit-level columns: the orbit vector v (number of neighbours of a vertex of the column's orbit in each of the 42
     orbits) is
       free neighbour q      -> orbit (j, q)   [pair {j,q}, free at j]        : 1
       square, pair {i,k}    -> orbit (k, i)   [free at k, paired via i]      : 2
       far pair {a,b}, value 2 -> orbits (a, b) and (b, a)                    : 1 each
       far pair {a,b}, value 1 -> orbit (a, b) OR orbit (b, a)                : 1   (every choice is a column)
  2. Every row is rebuilt FROM ITS LABEL (its coefficients and its bounds; the bounds stored in the file are only
     compared, never used):
       Ysym [j,a,b]  : [(j,a) has b as free neighbour] - [(j,b) has a]                     = 0
       Ysub [j,S]    : free edges inside S, each counted once                              in [0, |S| - 1]
       PPE  [a,b]    : columns with value 2 on the far pair {a,b}                          = 2
       C    [U,W]    : q_U(W) on the columns of U, - q_W(U) on the columns of W            = 0
       P2o  [U,W]    : DOUBLED orbit identity  2 sum_Z q_Z(U) q_Z(W) + q_U(W) + q_W(U)
                                                = 2 (12 [U = W] + 4 - |blocks(U) & blocks(W)|)
  3. Per-orbit Farkas inequality (the O rows, one unit per orbit, are used through the max):
       low  = sum_{rows except O} min(y lo, y hi)
       high = sum_Z max over ALL orbit-level columns c of orbit Z of sum_r y_r A_rc
     low > high is impossible for a solution x >= 0 with one unit per orbit.
  4. Coverage: every labelled triangle-free decomposition of K7 is isomorphic to a certified representative (own
     canonical form canon_any), and the labelled counts add up."""
import gzip
import itertools
import json
import os
import sys

sys.path.insert(0, 'scripts')
from verify_core_certs import my_columns, canon_any
from zeros_budget import decompositions

V7 = range(7)


def expand(col):
    """(orbit, p, fixed orbit vector, list of 2-option far pairs) for one pair-level column."""
    (j, i), p, val = col
    fixed = {}
    opts = []
    for Q, v in val.items():
        if v == 0:
            continue
        a, b = sorted(Q)
        if Q == frozenset((i, j)):
            raise ValueError('own pair with a nonzero value')
        if j in Q:                                     # free neighbour at j
            q = a if b == j else b
            assert v == 1 and q in p
            fixed[(j, q)] = fixed.get((j, q), 0) + 1
        elif i in Q:                                   # square: pair {i,k}, orbit free at k paired via i
            k = a if b == i else b
            assert v == 2
            fixed[(k, i)] = fixed.get((k, i), 0) + 2
        elif v == 2:
            fixed[(a, b)] = fixed.get((a, b), 0) + 1
            fixed[(b, a)] = fixed.get((b, a), 0) + 1
        else:
            assert v == 1
            opts.append(((a, b), (b, a)))
    return (j, i), p, fixed, opts


def main(n):
    fn = f'logs/orbit_certs/zero_{n:02d}.json.gz'
    data = json.load(gzip.open(fn, 'rt', encoding='utf-8'))
    dec = [tuple(c) for c in data['decomposition']]
    assert data['p2o_rows_doubled']
    ysym, ysub, yppe, yc_plus, yc_minus, yp = {}, {}, {}, {}, {}, {}
    low = 0
    bound_mismatch = 0
    for label, (lo_f, hi_f), y in data['rows']:
        kind, key = label
        if kind == 'O':
            continue
        if kind == 'Ysym':
            j, a, b = key
            ysym[(j, a, b)] = y; lo = hi = 0
        elif kind == 'Ysub':
            j, S = key
            ysub[(j, tuple(S))] = y; lo, hi = 0, len(S) - 1
        elif kind == 'PPE':
            yppe[frozenset(key)] = y; lo = hi = 2
        elif kind == 'C':
            U, W = tuple(key[0]), tuple(key[1])
            yc_plus.setdefault(U, {})[W] = yc_plus.get(U, {}).get(W, 0) + y
            yc_minus.setdefault(W, {})[U] = yc_minus.get(W, {}).get(U, 0) + y
            lo = hi = 0
        elif kind == 'P2o':
            U, W = tuple(key[0]), tuple(key[1])
            yp[frozenset((U, W))] = yp.get(frozenset((U, W)), 0) + y
            sh = len(set(U) & set(W))
            lo = hi = 2 * (12 * (U == W) + 4 - sh)
        else:
            raise ValueError(kind)
        bound_mismatch += (lo, hi) != (lo_f, hi_f)
        low += min(y * lo, y * hi)
    best = {}
    ncols = 0
    for col in my_columns(dec):
        Z, p, fixed, opts = expand(col)
        j, i = Z
        static = 0
        for q in p:                                    # Ysym
            a, b = sorted((i, q))
            static += ysym.get((j, a, b), 0) * (1 if a == i else -1)
        for (jj, S), y in ysub.items():                # Ysub
            if jj == j and i in S:
                static += y * sum(1 for q in p if q in S and q > i)
        for Q, y in yppe.items():                      # PPE
            if not (Q & {i, j}) and col[2][Q] == 2:
                static += y
        cp, cm_ = yc_plus.get(Z, {}), yc_minus.get(Z, {})
        for choice in itertools.product((0, 1), repeat=len(opts)):
            v = dict(fixed)
            for (o1, o2), c in zip(opts, choice):
                o = o2 if c else o1
                v[o] = v.get(o, 0) + 1
            w = static
            for W, x in v.items():                     # C rows
                w += x * (cp.get(W, 0) - cm_.get(W, 0))
            sup = list(v.items())
            for s in range(len(sup)):                  # P2o rows: 2 v_U v_W on every column
                U, xu = sup[s]
                w += 2 * yp.get(frozenset((U,)), 0) * xu * xu
                for t in range(s + 1, len(sup)):
                    W, xw = sup[t]
                    w += 2 * yp.get(frozenset((U, W)), 0) * xu * xw
            for W, x in v.items():                     # P2o rows: + q_Z(W) on the columns of Z (U = Z)
                if W != Z:
                    w += yp.get(frozenset((Z, W)), 0) * x
            ncols += 1
            if Z not in best or w > best[Z]:
                best[Z] = w
    high = sum(best.values())
    ok = low > high and bound_mismatch == 0 and len(best) == 42
    print(f'{fn}: {ncols} orbit-level columns, {len(best)} orbits; low = {low}, high = {high}, low - high = {low - high}; '
          f'bounds differing from the recomputed ones: {bound_mismatch}  ->  {"VALID" if ok else "NOT VALID"}', flush=True)
    return ok, dec, data['labelled']


if __name__ == '__main__':
    which = [int(a) for a in sys.argv[1:] if a.isdigit()]
    if not which:
        which = sorted(int(f[5:7]) for f in os.listdir('logs/orbit_certs') if f.startswith('zero_'))
    ok_all, reps = True, []
    for n in which:
        ok, dec, lab = main(n)
        ok_all &= ok
        reps.append((dec, lab))
    if 'coverage' in sys.argv:
        forms = {canon_any(d) for d, _ in reps}
        total = miss = 0
        for dec in decompositions((1 << 21) - 1, []):
            if any(len(c) == 3 for c in dec):
                continue
            total += 1
            miss += canon_any([tuple(c) for c in dec]) not in forms
        print(f'coverage: {total} labelled triangle-free decompositions, {miss} not isomorphic to a certified '
              f'representative; labelled counts of the certified classes add up to {sum(l for _, l in reps)}')
        ok_all &= miss == 0 and sum(l for _, l in reps) == total
    print('ALL VERIFIED' if ok_all else 'NOT ALL VERIFIED')
