"""Machine check (a check only, the proof is by hand) of the hand proof of family 4, the triangle of triangles, sent to
C2 on 25/09 and accepted by C2.  Own solver of the local system C0-C8 in the ORIGINAL block names, independent of the
canonical keys of doubles_budget3 (used only for a final cross-check).  Checked:
  (a) the Complement Lemma on every key of every block, and on every key of logs/doubles3_keys*.json;
  (b) X1, X2, X3 step by step, and that their exact minimum is 2;
  (c) per block (w = 3, p = 1, a = 0): the stated rules, the surviving cycles, the stated bound per orbit;
      every stated bound is <= the exact minimum and every excluded cycle is really infeasible;
  (d) the Z3 symmetry and the total 4 + 3*7 + 3*6 = 43."""
import itertools
import json
import sys

sys.path.insert(0, 'scripts')
from zeros_budget import V, CYC6
import doubles_budget3 as d3
from local_system import P, frame, solve_core, doubles

PASCH = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 6), (1, 5, 6), (2, 3, 4, 5), (2, 4, 6)]

_cache = {}


def solutions(dec, j, i, pair):
    key = (id(dec), j, i, frozenset(pair))
    if key not in _cache:
        cyc, k, R = frame(dec, j, i)
        tris = [set(c) for c in dec if len(c) == 3 and set(c) != set(cyc)]
        _cache[key] = solve_core(tris, len(cyc) == 3, i, j, k, R, pair)
    return _cache[key]


def exact(dec, j, i, pair):
    s = solutions(dec, j, i, pair)
    return None if not s else min(doubles(x) for x in s)


def cycles_at(j):
    others = [r for r in V if r != j]
    for order in CYC6:
        lab = [others[p] for p in order]
        yield lab, {lab[t]: frozenset((lab[t - 1], lab[(t + 1) % 6])) for t in range(6)}


def same_cycle(lab, target):
    m = len(target)
    rots = [tuple(target[(s + t) % m] for t in range(m)) for s in range(m)]
    rots += [tuple(reversed(r)) for r in rots]
    return tuple(lab) in rots


# ---------------------------------------------------------------- (a) Complement Lemma
def cl_check(tris, own_tri, i, j, k, R, pair, sols):
    n = {r: int(r in pair) for r in V}
    for t in tris:
        if k in t and i not in t and j not in t:
            cd = t - {k}
            a, b = sorted(set(R) - cd)
            N = n[a] + n[b]
            if N == 0:
                assert not sols, ('CL: N_ab = 0 but feasible', i, j, k, R, pair)
            for phi, y in sols:
                assert phi[P(a, b)] == 3 - N, ('CL fails', i, j, k, R, pair)


count = 0
for j in V:
    for i in V:
        if i == j:
            continue
        cyc, k, R = frame(PASCH, j, i)
        tris = [set(c) for c in PASCH if len(c) == 3 and set(c) != set(cyc)]
        for pair in itertools.combinations([r for r in V if r not in (i, j)], 2):
            cl_check(tris, len(cyc) == 3, i, j, k, R, set(pair), solutions(PASCH, j, i, pair))
            count += 1
print(f'(a) CL holds on all {count} (block, orbit, pair) keys of the triangle of triangles')
nkeys = 0
for fn in ('logs/doubles3_keys.json', 'logs/doubles3_keys_extra.json', 'logs/doubles3_keys_t3.json'):
    for tris_l, nbrs, val in json.load(open(fn)):
        tris_all = [set(t) for t in tris_l]
        own = {0, 1, 2} in tris_all
        tris = [t for t in tris_all if t != {0, 1, 2}]
        sols = solve_core(tris, own, 0, 1, 2, [3, 4, 5, 6], set(nbrs))
        mn = None if not sols else min(doubles(x) for x in sols)
        assert mn == val, ('own solver disagrees with doubles_budget3', tris_l, nbrs, mn, val)
        cl_check(tris, own, 0, 1, 2, [3, 4, 5, 6], set(nbrs), sols)
        nkeys += 1
print(f'    CL holds on all {nkeys} canonical keys of the exhaustive runs (and the own solver agrees with them)')

# ---------------------------------------------------------------- (b) X1, X2, X3
def check_extra(name, j, i, pair, k_exp, R_exp, steps):
    cyc, k, R = frame(PASCH, j, i)
    assert k == k_exp and R == R_exp, (name, k, R)
    sols = solutions(PASCH, j, i, pair)
    assert sols, name
    for phi, y in sols:
        f = lambda a, b: phi[P(a, b)]
        for label, cond in steps:
            assert cond(f, y), (name, label)
        assert doubles((phi, y)) >= 2
    mn = min(doubles(x) for x in sols)
    assert mn == 2, (name, mn)
    print(f'(b) {name} (block {j}, orbit {i}, pair {sorted(pair)}): {len(sols)} solutions, every step holds, '
          f'exact minimum {mn}')


check_extra('X1', 3, 5, (1, 6), 0, [1, 2, 4, 6], [
    ('CL phi46 = 2', lambda f, y: f(4, 6) == 2),
    ('C8(6): y6 = phi26 = 0', lambda f, y: y[6] == 0 and f(2, 6) == 0),
    ('C5: phi24 = 1', lambda f, y: f(2, 4) == 1),
    ('y1 = 0', lambda f, y: y[1] == 0),
    ('y2 + phi12 = 3', lambda f, y: y[2] + f(1, 2) == 3)])
check_extra('X2', 0, 3, (1, 4), 1, [2, 4, 5, 6], [
    ('C3: y2 = 0', lambda f, y: y[2] == 0),
    ('CL: phi24 = 2', lambda f, y: f(2, 4) == 2),
    ('C5: phi26 + phi46 = 1', lambda f, y: f(2, 6) + f(4, 6) == 1)])
check_extra('X3', 0, 5, (1, 6), 3, [1, 2, 4, 6], [
    ('C2: phi12 = 0', lambda f, y: f(1, 2) == 0),
    ('C8(1): y1 + phi14 = 2', lambda f, y: y[1] + f(1, 4) == 2),
    ('C5 + C8(6): phi24 = 1 + y6', lambda f, y: f(2, 4) == 1 + y[6]),
    ('phi14 = 2 + y1', lambda f, y: f(1, 4) == 2 + y[1]),
    ('y1 = 0, phi14 = 2', lambda f, y: y[1] == 0 and f(1, 4) == 2),
    ('(y2,y4,y6) in the 3 options', lambda f, y: (y[2], y[4], y[6]) in {(2, 0, 0), (1, 1, 0), (1, 0, 1)}),
    ('second double', lambda f, y: y[2] == 2 or f(2, 6) == 2 or f(2, 4) == 2)])


# ---------------------------------------------------------------- (c) blocks
def immediate_bound(j, i, pair):
    """Stated immediate bound: forced phi from C1 (1), C2 (1 - N), CL (3 - N); max(#2, n_k + #0) (lemma D)."""
    cyc, k, R = frame(PASCH, j, i)
    n = {r: int(r in pair) for r in V}
    forced = {}
    for c in PASCH:
        t = set(c)
        if len(c) != 3 or t == set(cyc):
            continue
        if i in t:
            forced[P(*(t - {i}))] = 1
        elif j in t and k not in t:
            a, b = sorted(t - {j})
            forced[P(a, b)] = 1 - n[a] - n[b]
        elif k in t and j not in t:
            a, b = sorted(set(R) - (t - {k}))
            forced[P(a, b)] = 3 - n[a] - n[b]
    return max(sum(v == 2 for v in forced.values()), n[k] + sum(v == 0 for v in forced.values()))


def check_block(j, rules, survivors, stated, extras, target):
    got = []
    for lab, nb in cycles_at(j):
        ok = all(rule(nb) for rule in rules)
        feasible = all(solutions(PASCH, j, i, nb[i]) for i in nb)
        if not ok:
            assert not feasible, ('a cycle excluded by the stated rules is feasible', j, lab)
            continue
        got.append(lab)
    assert len(got) == len(survivors), (j, len(got), got)
    for target_cyc, bounds in zip(survivors, stated):
        lab = next(l for l in got if same_cycle(l, target_cyc))
        nb = {target_cyc[t]: frozenset((target_cyc[t - 1], target_cyc[(t + 1) % 6])) for t in range(6)}
        for i, b in zip(target_cyc, bounds):
            imm = extras.get((i, nb[i]), immediate_bound(j, i, nb[i]))
            ex = exact(PASCH, j, i, nb[i])
            assert b == imm, ('stated bound differs from the rule', j, target_cyc, i, b, imm)
            assert ex is None or b <= ex, ('stated bound above the exact minimum', j, target_cyc, i)
        assert sum(bounds) >= target, (j, target_cyc)
    print(f'(c) block {j}: the stated rules leave exactly {len(got)} cycles, the listed ones; stated bounds = rules, '
          f'<= exact minima; each cycle >= {target}')


def has_nb_in(i, S):
    return lambda nb: len(nb[i] & set(S)) >= 1


def not_pair(i, S):
    return lambda nb: i not in nb or nb[i] != frozenset(S)


def not_adj(a, b):
    return lambda nb: b not in nb[a]


# block w = 3: CL rules; bounds N = 1 -> 1, N = 2 -> 0, X1
CL3 = {0: (1, 5), 1: (0, 5), 2: (0, 4), 4: (0, 2), 5: (4, 6), 6: (4, 5)}
rules3 = [has_nb_in(i, S) for i, S in CL3.items()]
surv3 = []
for lab, nb in cycles_at(3):
    if all(r(nb) for r in rules3):
        surv3.append(lab)
        z = sum(len(nb[i] & set(S)) == 2 for i, S in CL3.items())
        assert z <= 3
        if z == 3:
            assert same_cycle(lab, (0, 1, 5, 6, 4, 2)), lab
        bound = sum(len(nb[i] & set(S)) == 1 for i, S in CL3.items())
        if z == 3:
            assert nb[5] == frozenset((1, 6))
            bound += 1                                              # X1: orbit 5 gets 2 instead of 1
        assert bound >= 4, lab
        for i, S in CL3.items():
            b = int(len(nb[i] & set(S)) == 1) + int(z == 3 and i == 5)
            ex = exact(PASCH, 3, i, nb[i])
            assert ex is None or b <= ex
    else:
        assert not all(solutions(PASCH, 3, i, nb[i]) for i in nb)
print(f'(c) block 3: {len(surv3)} cycles pass the CL rules; z <= 3; z = 3 only on (0,1,5,6,4,2); every cycle >= 4 '
      f'(N = 1 doubles + X1); stated bounds <= exact minima')

rules1 = [not_adj(0, 2), not_adj(5, 6)] + [not_pair(i, (5, 6)) for i in (0, 2, 3, 4)] + \
         [not_pair(i, (0, 2)) for i in (3, 4, 5, 6)] + [has_nb_in(0, (3, 5)), has_nb_in(5, (0, 3))]
check_block(1, rules1, [(0, 3, 6, 2, 4, 5), (0, 3, 5, 2, 4, 6), (0, 4, 6, 2, 3, 5)],
            [(1, 2, 1, 1, 1, 1)] * 3, {}, 7)
rules0 = [not_adj(1, 2)] + [not_pair(i, (1, 2)) for i in (3, 4, 5, 6)] + \
         [has_nb_in(1, (3, 5)), has_nb_in(2, (3, 4)), has_nb_in(3, (2, 4)), has_nb_in(4, (2, 3))]
check_block(0, rules0,
            [(1, 3, 4, 2, 5, 6), (1, 3, 4, 2, 6, 5), (1, 4, 3, 2, 6, 5), (1, 5, 3, 2, 4, 6), (1, 5, 4, 2, 3, 6)],
            [(1, 2, 0, 1, 1, 1), (0, 2, 0, 1, 1, 2), (1, 1, 0, 1, 1, 2), (1, 2, 1, 0, 1, 1), (1, 1, 1, 0, 1, 2)],
            {(3, frozenset((1, 4))): 2, (5, frozenset((1, 6))): 2}, 6)

# ---------------------------------------------------------------- (d) symmetry, totals, cross-check
rho = {1: 6, 6: 2, 2: 1, 0: 5, 5: 4, 4: 0, 3: 3}
edges = lambda c: frozenset(frozenset((c[t], c[(t + 1) % len(c)])) for t in range(len(c)))
assert {edges(c) for c in PASCH} == {edges(tuple(rho[v] for v in c)) for c in PASCH}
block_min = {}
for j in V:
    best = None
    for lab, nb in cycles_at(j):
        vals = [exact(PASCH, j, i, nb[i]) for i in nb]
        if None not in vals and (best is None or sum(vals) < best):
            best = sum(vals)
    block_min[j] = best
assert block_min[1] == block_min[6] == block_min[2] and block_min[0] == block_min[5] == block_min[4]
print(f'(d) rho = (1 6 2)(0 5 4) is an automorphism; exact block minima {block_min}; exact total '
      f'{sum(block_min.values())}; hand total 4 + 3*7 + 3*6 = 43 > 42')
nx = 0
for j in V:
    for i in V:
        if i == j:
            continue
        cyc, k, R = frame(PASCH, j, i)
        tris = [tuple(sorted(c)) for c in PASCH if len(c) == 3]
        for pair in itertools.combinations([r for r in V if r not in (i, j)], 2):
            assert exact(PASCH, j, i, pair) == d3.local_min(d3.canon(tris, i, j, k, R, frozenset(pair)))
            nx += 1
print(f'    own solver = doubles_budget3 on all {nx} keys of the triangle of triangles')
print('ALL CHECKS PASSED')
