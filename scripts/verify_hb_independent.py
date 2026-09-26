"""INDEPENDENT verification of C2's H-beta certificates (amazing-source/conway99-c2-experimental, logs/hb_certs),
written in the lab 'structures' from C2's DESCRIPTION of the model only (message of 26/09), without C2's code.
  1. columns rebuilt from the triplet list (i,j,k) of each record: z free at j, paired via i, third block k;
     free neighbours p = 2 orbits free at j, not sisters of z; typed rule (skip p if both are in {j,k});
     near values (0 on its pair, 2[R = {i,k}], #{y in p : pair(y) = R} for R through j); far values = all {0,1,2}
     weightings of the 10 far pairs with vertex sums deg(t) = 4 - near contacts at t;
  2. every certificate row rebuilt from its label ('one','sym','sub','sis','S','D','C','P2');
  3. rigorous check of  y^T b  outside  [LB, UB]  with x in [0,1] and slacks of the <= rows in [0, rhs]:
     y^T b in exact rationals; w = A^T y in float64 with an explicit a priori error bound
     |w_c - w^_c| <= (n+4) 2^-52 sum_r |y_r a_rc|;
  4. coverage: the 7 680 skeletons (60 orientations x 2^7 path orders) regenerated, canonical forms under D7
     recomputed, compared with the records (classes and multiplicities)."""
import gzip
import itertools
import json
import os
import random
import sys
from fractions import Fraction

import numpy as np

PAIRS = [frozenset(p) for p in itertools.combinations(range(7), 2)]
PIDX = {P: t for t, P in enumerate(PAIRS)}


def parse(lab):
    return eval(lab, {'__builtins__': {}}, {'frozenset': frozenset})


def weightings(far, deg):
    """All maps far-pair -> {0,1,2} with sum over pairs through t equal to deg[t]."""
    out = []
    far = list(far)

    def rec(t, w, rem):
        if t == len(far):
            if all(v == 0 for v in rem.values()):
                out.append(dict(w))
            return
        a, b = sorted(far[t])
        for v in (0, 1, 2):
            if rem[a] >= v and rem[b] >= v:
                rem[a] -= v; rem[b] -= v; w[far[t]] = v
                rec(t + 1, w, rem)
                rem[a] += v; rem[b] += v
        w.pop(far[t], None)
    rec(0, {}, dict(deg))
    return out


def columns(orbits):
    orbits = [tuple(z) for z in orbits]
    pair = {z: frozenset(z[:2]) for z in orbits}
    cols = []
    for z in orbits:
        i, j, k = z
        fr = [y for y in orbits if y != z and y[1] == j and pair[y] != pair[z]]
        others = [t for t in range(7) if t not in (i, j)]
        far = [frozenset(q) for q in itertools.combinations(others, 2)]
        for y1, y2 in itertools.combinations(fr, 2):
            if pair[y1] == frozenset((j, k)) and pair[y2] == frozenset((j, k)):
                continue
            near = {pair[z]: 0}
            for t in others:
                near[frozenset((i, t))] = 2 if t == k else 0
                near[frozenset((j, t))] = (pair[y1] == frozenset((j, t))) + (pair[y2] == frozenset((j, t)))
            deg = {t: 4 - near[frozenset((i, t))] - near[frozenset((j, t))] for t in others}
            if min(deg.values()) < 0:
                continue
            for w in weightings(far, deg):
                val = [0] * 21
                for P, v in near.items():
                    val[PIDX[P]] = v
                for P, v in w.items():
                    val[PIDX[P]] = v
                cols.append((z, frozenset((y1, y2)), val))
    return cols


def row(label, cols, zpair, VAL, Z, Pset):
    """(coefficient array over columns, rhs, is_le) rebuilt from the label."""
    t = parse(label)
    kind = t[0]
    n = len(cols)
    if kind == 'one':
        return np.array([1 if Z[c] == t[1] else 0 for c in range(n)], dtype=np.int64), 1, False
    if kind == 'sym':
        a, b = t[1], t[2]
        return np.array([(1 if (Z[c] == a and b in Pset[c]) else 0) - (1 if (Z[c] == b and a in Pset[c]) else 0)
                         for c in range(n)], dtype=np.int64), 0, False
    if kind == 'sub':
        S = t[1]
        coef = np.zeros(n, dtype=np.int64)
        for a in S:
            for b in S:
                if a != b:
                    coef += np.array([1 if (Z[c] == a and b in Pset[c]) else 0 for c in range(n)], dtype=np.int64)
        return coef, 2 * (len(S) - 1), True
    if kind == 'sis':
        a, b = t[1], t[2]
        return np.array([1 if Pset[c] == frozenset((a, b)) else 0 for c in range(n)], dtype=np.int64), 1, True
    if kind == 'S':
        return VAL[:, PIDX[t[1]]].copy(), 24, False
    if kind == 'D':
        return VAL[:, PIDX[t[1]]] ** 2, 32, False
    if kind == 'C':
        R, R2 = t[1], t[2]
        inR = np.array([zpair[c] == R for c in range(n)])
        inR2 = np.array([zpair[c] == R2 for c in range(n)])
        return inR * VAL[:, PIDX[R2]] - inR2 * VAL[:, PIDX[R]], 0, False
    if kind == 'P2':
        R, R2 = t[1], t[2]
        inR = np.array([zpair[c] == R for c in range(n)])
        return VAL[:, PIDX[R]] * VAL[:, PIDX[R2]] + inR * VAL[:, PIDX[R2]], 16 - 4 * (len(R & R2) > 0), False
    raise ValueError(label)


def verify_record(d):
    cols = columns(d['orbits'])
    n = len(cols)
    Z = [c[0] for c in cols]
    Pset = [c[1] for c in cols]
    zpair = [frozenset(c[0][:2]) for c in cols]
    VAL = np.array([c[2] for c in cols], dtype=np.int64)
    w_hat = np.zeros(n)
    w_abs = np.zeros(n)
    ytb = Fraction(0)
    le_hi = Fraction(0)
    le_lo = Fraction(0)
    nterms = 0
    for lab, ys in d['y']:
        y = Fraction(ys)
        coef, rhs, is_le = row(lab, cols, zpair, VAL, Z, Pset)
        yf = float(y)
        w_hat += yf * coef
        w_abs += abs(yf) * np.abs(coef)
        nterms += 1
        ytb += y * rhs
        if is_le:
            le_hi += max(0, y * rhs)
            le_lo += min(0, y * rhs)
    err = (nterms + 4) * 2.0 ** -52 * w_abs + 1e-300
    UB = float(np.maximum(0, w_hat + err).sum()) * (1 + 1e-10) + float(le_hi)
    LB = float(np.minimum(0, w_hat - err).sum()) * (1 + 1e-10) + float(le_lo)
    ok = float(ytb) > UB + 1e-6 or float(ytb) < LB - 1e-6
    return ok, n, float(ytb), UB, LB


def canon(triplets):
    best = None
    for s in (1, -1):
        for r in range(7):
            f = tuple(sorted(tuple((s * x + r) % 7 for x in t) for t in triplets))
            if best is None or f < best:
                best = f
    return best


def coverage(records):
    Tset = {frozenset((c, (c + 1) % 7)) for c in range(7)}
    chords = [frozenset(e) for e in itertools.combinations(range(7), 2) if frozenset(e) not in Tset]
    count = {}
    for bits in range(1 << len(chords)):
        indeg = [0] * 7
        inn = {c: [] for c in range(7)}
        for t, e in enumerate(chords):
            a, b = sorted(e)
            head, tail = (b, a) if bits >> t & 1 else (a, b)
            indeg[head] += 1
            inn[head].append(tail)
        if any(dd != 2 for dd in indeg):
            continue
        for order in range(128):
            trip = []
            for c in range(7):
                A, B = sorted(inn[c])
                seq = [(c - 1) % 7, A, B, (c + 1) % 7] if not (order >> c & 1) else [(c - 1) % 7, B, A, (c + 1) % 7]
                for x, y in zip(seq, seq[1:]):
                    trip += [(c, x, y), (c, y, x)]
            key = canon(trip)
            count[key] = count.get(key, 0) + 1
    rec = {canon([tuple(z) for z in d['orbits']]): d['mult'] for d in records}
    same = set(rec) == set(count)
    mult_ok = all(rec[k] == count[k] for k in rec if k in count)
    return len(count), sum(count.values()), same, mult_ok


if __name__ == '__main__':
    root = sys.argv[1]
    sample = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    records = []
    for fn in sorted(os.listdir(root)):
        if fn.startswith('hb_colall'):
            with gzip.open(os.path.join(root, fn), 'rt') as f:
                records += [json.loads(line) for line in f]
    ncl, tot, same, mult_ok = coverage(records)
    print(f'coverage: {ncl} classes regenerated, {tot} skeletons; same classes as the records: {same}; '
          f'multiplicities equal: {mult_ok}', flush=True)
    todo = records if not sample else random.Random(26).sample(records, sample)
    good = 0
    for d in todo:
        ok, n, ytb, UB, LB = verify_record(d)
        good += ok
        print(f"record {d['idx']} (mult {d['mult']}): {n} columns (record says {d['ncols']}); y^T b = {ytb:.3f}, "
              f"UB <= {UB:.3f}, LB >= {LB:.1f} -> {'VALID' if ok else 'NOT VALID'}", flush=True)
    print(f'{good}/{len(todo)} certificates valid (independent reconstruction)')
