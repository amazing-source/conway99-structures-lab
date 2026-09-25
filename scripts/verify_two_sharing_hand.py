"""Machine check (a check only; the proofs are by hand) of the hand proofs for the pure twisted decompositions with
two triangles sharing a vertex.  Labels: v = 0, T1 = {0,1,2}, T2 = {0,3,4}, e = 5, f = 6.
Checked, with the exact local solver (local_system.py):
  (T) Lemma T on every canonical key of the exhaustive runs: if R_x = T u {r} with T a triangle, then
      y_r = 1 - n_k - n_r;
  (S) every local step of the hand proofs, for EVERY neighbour pair allowed by the facts invoked at that step.
Each step reads: orbit (j;i) with neighbours constrained by `allow`; `claim` holds on all solutions."""
import itertools
import json
import sys

sys.path.insert(0, 'scripts')
from zeros_budget import V
from local_system import P, frame, solve_core

# ---------------------------------------------------------------- (T) Lemma T on all canonical keys
cnt = 0
for fn in ('logs/doubles3_keys.json', 'logs/doubles3_keys_extra.json', 'logs/doubles3_keys_t3.json'):
    for tris_l, nbrs, val in json.load(open(fn)):
        tris_all = [set(t) for t in tris_l]
        own = {0, 1, 2} in tris_all
        tris = [t for t in tris_all if t != {0, 1, 2}]
        for T in tris:
            if T <= {3, 4, 5, 6}:
                (r,) = {3, 4, 5, 6} - T
                n = {q: int(q in nbrs) for q in V}
                for phi, y in solve_core(tris, own, 0, 1, 2, [3, 4, 5, 6], set(nbrs)):
                    assert y[r] == 1 - n[2] - n[r], ('Lemma T fails', tris_l, nbrs)
                cnt += 1
print(f'(T) Lemma T holds on all {cnt} (key, triangle in R) cases of the exhaustive runs')


def sols(dec, j, i, pair):
    cyc, k, R = frame(dec, j, i)
    tris = [set(c) for c in dec if len(c) == 3 and set(c) != set(cyc)]
    return solve_core(tris, len(cyc) == 3, i, j, k, R, set(pair))


def val(dec, j, i, Q, phi, y):
    """phi_X(Q) for the orbit (j;i): a phi value if Q is in R_x, a y value if Q = {k, r}."""
    cyc, k, R = frame(dec, j, i)
    Q = set(Q)
    if k in Q:
        (r,) = Q - {k}
        return y[r]
    return phi[P(*Q)]


def step(dec, j, i, k_exp, allow, claim, text):
    cyc, k, R = frame(dec, j, i)
    assert k == k_exp, (text, 'k is', k)
    others = [r for r in V if r not in (i, j)]
    npairs = 0
    for pair in itertools.combinations(others, 2):
        if not allow(set(pair)):
            continue
        for phi, y in sols(dec, j, i, pair):
            assert claim(lambda *Q: val(dec, j, i, Q, phi, y), set(pair)), (text, pair)
        npairs += 1
    assert npairs > 0, (text, 'no admissible pair')
    return npairs


def rng(dec, j, i, allow, Q):
    """(min, max) of phi_X(Q) over all solutions with admissible neighbour pairs."""
    vals = []
    others = [r for r in V if r not in (i, j)]
    for pair in itertools.combinations(others, 2):
        if allow(set(pair)):
            vals += [val(dec, j, i, Q, phi, y) for phi, y in sols(dec, j, i, pair)]
    return min(vals), max(vals)


ANY = lambda p: True
CLASSES = {
    '0 (3,3,4,4,7)': [(1, 2, 0), (1, 3, 2, 5, 0, 6, 4), (1, 5, 3, 6), (2, 6, 5, 4), (0, 3, 4)],
    '2 (3,3,4,5,6)b': [(1, 2, 0), (1, 5, 2, 6), (1, 3, 5, 0, 6, 4), (2, 3, 6, 5, 4), (0, 3, 4)],
    '1 (3,3,4,5,6)a': [(3, 4, 0), (3, 1, 4, 6, 0, 5), (3, 6, 1, 5, 2), (4, 5, 6, 2), (0, 1, 2)],
    '3 (3,3,4,5,6)c': [(1, 0, 2), (1, 5, 0, 6, 3), (1, 6, 5, 3, 2, 4), (0, 3, 4), (2, 5, 4, 6)],
    '4 (3,3,4,5,6)d': [(1, 0, 2), (1, 5, 0, 6, 3), (1, 6, 2, 3, 5, 4), (0, 3, 4), (2, 5, 6, 4)],
    '5 (3,3,4,5,6)e': [(1, 0, 2), (1, 5, 0, 6), (1, 3, 2, 5, 6, 4), (0, 3, 4), (2, 6, 3, 5, 4)],
    '6 (3,3,4,5,6)f': [(1, 0, 2), (1, 5, 0, 6), (1, 3, 5, 2, 6, 4), (0, 3, 4), (2, 3, 6, 5, 4)],
    '7 (3,3,5,5,5)': [(3, 4, 0), (3, 5, 4, 6, 2), (3, 6, 0, 5, 1), (4, 2, 5, 6, 1), (0, 2, 1)],
}

# ---------------------------------------------------------------- Lemma A (general), checked on every class
for name, dec in CLASSES.items():
    # orbits of {1,2}: phi(34) = 1 + y5 + y6 ; orbits of {3,4}: phi(12) = 1 + y5 + y6 (k = 0, own triangle)
    step(dec, 1, 2, 0, ANY, lambda f, p: f(3, 4) == 1 + f(0, 5) + f(0, 6), f'{name} A: Y')
    step(dec, 2, 1, 0, ANY, lambda f, p: f(3, 4) == 1 + f(0, 5) + f(0, 6), f'{name} A: Y\'')
    step(dec, 3, 4, 0, ANY, lambda f, p: f(1, 2) == 1 + f(0, 5) + f(0, 6), f'{name} A: Z')
    step(dec, 4, 3, 0, ANY, lambda f, p: f(1, 2) == 1 + f(0, 5) + f(0, 6), f'{name} A: Z\'')
    # orbits of {0,5}, {0,6}: the one paired via 0 has phi(12) = phi(34) = 1; the one free at 0 has
    # phi(12) = 1 - N12 and phi(34) = 1 - N34
    for s in (5, 6):
        cyc, k, R = frame(dec, s, 0)
        step(dec, s, 0, k, ANY, lambda f, p: f(1, 2) == 1 and f(3, 4) == 1, f'{name} A: ({s};0)')
        cyc, k, R = frame(dec, 0, s)
        step(dec, 0, s, k, ANY, lambda f, p: f(1, 2) == 1 - len(p & {1, 2}) and f(3, 4) == 1 - len(p & {3, 4}),
             f'{name} A: (0;{s})')
print('(A) Lemma A: its local identities hold in all 7 classes (the coupling step is by hand)')


def W_step(dec, name):
    """Corollary A1: W = (6;5), W' = (5;6): CL at the base pair of the triangle not containing k."""
    out = {}
    for j, i in ((6, 5), (5, 6)):
        cyc, k, R = frame(dec, j, i)
        base = (3, 4) if k in (1, 2) else (1, 2)
        step(dec, j, i, k, ANY, lambda f, p, b=base: f(*b) == 3 - len(p & set(b)), f'{name} W ({j};{i})')
        out[j] = set(base)                         # with PPE on the base pair: neighbours of i at block j = base
    return out


print('(S) class steps:')
for name, dec in CLASSES.items():
    nb = W_step(dec, name)                         # nb[6] = neighbours of 5 at block 6, nb[5] = neighbours of 6 at block 5
    print(f'   {name}: at block 6, 5 ~ {sorted(nb[6])}; at block 5, 6 ~ {sorted(nb[5])}')
print('   (per-class finishing steps below)')

# class 7: coupling {3,6}-{4,5}
d = CLASSES['7 (3,3,5,5,5)']
lo63, _ = rng(d, 6, 3, lambda p: 5 in p and 4 not in p, (4, 5))          # block 6: 3 ~ 5, 3 !~ 4
lo36, _ = rng(d, 3, 6, ANY, (4, 5))
_, hi45 = rng(d, 4, 5, ANY, (3, 6))
_, hi54 = rng(d, 5, 4, ANY, (3, 6))
assert lo63 + lo36 > hi45 + hi54, (lo63, lo36, hi45, hi54)
print(f'   class 7: E[36][45] >= {lo63}+{lo36} > {hi45}+{hi54} >= E[45][36]')

# class 0: coupling {3,5}-{4,6}; block 6: 5 ~ {1,2}; block 5: 6 ~ {3,4}
d = CLASSES['0 (3,3,4,4,7)']
lo64, _ = rng(d, 6, 4, lambda p: 5 not in p, (3, 5))
lo46, _ = rng(d, 4, 6, ANY, (3, 5))
_, hi53 = rng(d, 5, 3, lambda p: 6 in p, (4, 6))
_, hi35 = rng(d, 3, 5, ANY, (4, 6))
assert lo64 + lo46 > hi53 + hi35, (lo64, lo46, hi53, hi35)
print(f'   class 0: E[46][35] >= {lo64}+{lo46} > {hi53}+{hi35} >= E[35][46]')

# class 1: step (b) (3;5): phi12 = 2 - n4 + y6, so PPE{1,2} forces n4 = 1 at block 3; then coupling {3,6}-{4,5}
d = CLASSES['1 (3,3,4,5,6)a']
step(d, 3, 5, 0, ANY, lambda f, p: f(1, 2) == 2 - len(p & {4}) + f(0, 6), 'class 1 (3;5)')
lo63, _ = rng(d, 6, 3, lambda p: 5 not in p, (4, 5))                      # block 6: 5 ~ {1,2}
lo36, _ = rng(d, 3, 6, lambda p: not {4, 5} <= p, (4, 5))                 # block 3: 5 ~ 4, so 6 not ~ both
_, hi54 = rng(d, 5, 4, ANY, (3, 6))
_, hi45 = rng(d, 4, 5, ANY, (3, 6))
assert lo63 + lo36 > hi54 + hi45, (lo63, lo36, hi54, hi45)
print(f'   class 1: E[36][45] >= {lo63}+{lo36} > {hi54}+{hi45} >= E[45][36]')

# class 3: coupling {3,6}-{4,5}; block 6: 5 ~ {1,2}; block 5: 6 ~ {3,4}
d = CLASSES['3 (3,3,4,5,6)c']
lo63, _ = rng(d, 6, 3, lambda p: 5 not in p, (4, 5))
lo36, _ = rng(d, 3, 6, ANY, (4, 5))
_, hi54 = rng(d, 5, 4, lambda p: 6 in p, (3, 6))
_, hi45 = rng(d, 4, 5, ANY, (3, 6))
assert lo63 + lo36 > hi54 + hi45, (lo63, lo36, hi54, hi45)
print(f'   class 3: E[36][45] >= {lo63}+{lo36} > {hi54}+{hi45} >= E[45][36]')

# class 5: coupling {1,5}-{2,6}; block 6: 5 ~ {3,4}; block 5: 6 ~ {1,2}
d = CLASSES['5 (3,3,4,5,6)e']
_, hi51 = rng(d, 5, 1, lambda p: 6 in p, (2, 6))
_, hi15 = rng(d, 1, 5, ANY, (2, 6))
lo62, _ = rng(d, 6, 2, lambda p: 5 not in p, (1, 5))
lo26, _ = rng(d, 2, 6, ANY, (1, 5))
assert lo62 + lo26 > hi51 + hi15, (lo62, lo26, hi51, hi15)
print(f'   class 5: E[26][15] >= {lo62}+{lo26} > {hi51}+{hi15} >= E[15][26]')

# class 4: coupling {3,4}-{5,6} gives block 3: 4 ~ {5,6}; then coupling {3,6}-{4,5}
d = CLASSES['4 (3,3,4,5,6)d']
lo, hi = rng(d, 5, 6, ANY, (3, 4))                   # W' : y3 with 6 ~ {1,2} at block 5
lo_, hi_ = rng(d, 5, 6, lambda p: p == {1, 2}, (3, 4))
assert lo_ == hi_ == 1
step(d, 3, 4, 0, ANY, lambda f, p: f(5, 6) == 3 - len(p & {5, 6}), 'class 4 Z')
step(d, 4, 3, 0, ANY, lambda f, p: f(5, 6) == 3 - len(p & {5, 6}), "class 4 Z'")
lo36, _ = rng(d, 3, 6, lambda p: 4 in p and 5 not in p, (4, 5))
lo63, _ = rng(d, 6, 3, lambda p: 5 in p and 4 not in p, (4, 5))
_, hi54 = rng(d, 5, 4, ANY, (3, 6))
_, hi45 = rng(d, 4, 5, lambda p: 3 in p, (3, 6))
assert lo36 + lo63 > hi54 + hi45, (lo36, lo63, hi54, hi45)
print(f'   class 4: E[36][45] >= {lo36}+{lo63} > {hi54}+{hi45} >= E[45][36]  (after E[56][34] = 1+1 forces 4 ~ {{5,6}} at block 3)')

# class 6: saturations and block-0 forcing
d = CLASSES['6 (3,3,4,5,6)f']
assert rng(d, 6, 5, lambda p: p == {1, 2}, (3, 4)) == (1, 1) and rng(d, 5, 6, lambda p: p == {1, 2}, (3, 4)) == (1, 1)
step(d, 3, 4, 0, ANY, lambda f, p: f(5, 6) == 3 - len(p & {5, 6}), 'class 6 Z')
step(d, 4, 3, 0, ANY, lambda f, p: f(5, 6) == 3 - len(p & {5, 6}), "class 6 Z'")
blk0 = lambda s: (lambda p: len(p & {1, 2}) == 1 and len(p & {3, 4}) == 1)       # Lemma A at block 0
step(d, 0, 5, 1, blk0(5), lambda f, p: f(3, 6) + f(4, 6) == 3, 'class 6 (0;5): phi36 + phi46 = 3')
step(d, 0, 6, 1, blk0(6), lambda f, p: f(3, 5) + f(4, 5) == 3, 'class 6 (0;6): phi35 + phi45 = 3')
row0 = lambda p: len(p & {3, 4}) == 1 and len(p & {5, 6}) == 1                   # block 0 cycle structure
for i, k in ((1, 2), (2, 1)):
    step(d, 0, i, k, row0, lambda f, p: f(3, 6) == 2 or f(4, 5) == 2 or f(3, 5) == 2 or f(4, 6) == 2,
         f'class 6 (0;{i}): a double at 36, 45, 35 or 46')
    step(d, 0, i, k, row0, lambda f, p: not (f(3, 5) <= 1 and f(4, 6) <= 1) or f(3, 6) == 2 or f(4, 5) == 2,
         f'class 6 (0;{i}): with phi35, phi46 <= 1, a double at 36 or 45')
for (j, i, k, allow, Q, text) in ((5, 3, 1, lambda p: 6 not in p, (4, 6), '(5;3) phi46 = 2'),
                                  (3, 5, 2, lambda p: 4 in p and 6 not in p, (4, 6), '(3;5) phi46 = 2'),
                                  (6, 4, 1, lambda p: 5 not in p, (3, 5), '(6;4) phi35 = 2'),
                                  (4, 6, 2, lambda p: 3 in p and 5 not in p, (3, 5), '(4;6) phi35 = 2'),
                                  (5, 4, 2, lambda p: 6 not in p, (3, 6), '(5;4) phi36 = 2'),
                                  (6, 3, 2, lambda p: 5 not in p, (4, 5), '(6;3) phi45 = 2')):
    step(d, j, i, k, allow, lambda f, p, Q=Q: f(*Q) == 2, 'class 6 ' + text)
print('   class 6: {4,6} <- (5;3),(3;5); {3,5} <- (6;4),(4;6); {3,6} <- (5;4),(0;5); {4,5} <- (6;3),(0;6); '
      '(0;1) needs a third double at 36 or 45')
# class 2 (b): ten orbits forced to double on the four cross pairs {1,3},{1,4},{2,3},{2,4}; PPE allows 8
d = CLASSES['2 (3,3,4,5,6)b']
for j, i in ((1, 2), (2, 1), (3, 4), (4, 3)):            # Y, Y', Z, Z': CL at {5,6} (k = 0)
    step(d, j, i, 0, ANY, lambda f, p: f(5, 6) == 3 - len(p & {5, 6}), f'class 2 ({j};{i}) CL56')
assert rng(d, 6, 5, lambda p: p == {1, 2}, (3, 4)) == (1, 1) and rng(d, 5, 6, lambda p: p == {1, 2}, (3, 4)) == (1, 1)
# block 1 (2 ~ {5,6}) and block 2 (1 ~ {5,6}): the other neighbours of 5 and 6 are in {3,4}
for j, i, t in ((1, 5, 2), (1, 6, 2), (2, 5, 1), (2, 6, 1)):
    cyc, k, R = frame(d, j, i)
    step(d, j, i, k, lambda p, t=t: t in p, lambda f, p, t=t: p in ({t, 3}, {t, 4}), f'class 2 ({j};{i}) other nbr in 34')
for j, i, Q in ((1, 3, (2, 4)), (1, 4, (2, 3)), (2, 3, (1, 4)), (2, 4, (1, 3))):
    cyc, k, R = frame(d, j, i)
    step(d, j, i, k, lambda p: p in ({0, 5}, {0, 6}), lambda f, p, Q=Q: f(*Q) == 2, f'class 2 ({j};{i}) double at {Q}')
either = ((3, 6, lambda p: 4 in p, ((1, 4), (2, 4))), (4, 5, lambda p: 3 in p, ((1, 3), (2, 3))),
          (5, 1, lambda p: 6 in p, ((2, 3), (2, 4))), (5, 2, lambda p: 6 in p, ((1, 3), (1, 4))),
          (6, 1, lambda p: 5 in p, ((2, 3), (2, 4))), (6, 2, lambda p: 5 in p, ((1, 3), (1, 4))))
for j, i, allow, (Q1, Q2) in either:
    cyc, k, R = frame(d, j, i)
    step(d, j, i, k, allow, lambda f, p, Q1=Q1, Q2=Q2: f(*Q1) + f(*Q2) == 3, f'class 2 ({j};{i}) one double in {Q1},{Q2}')
print('   class 2: (1;3),(1;4),(2;3),(2;4) double at the opposite cross pair; (3;6),(4;5),(5;1),(5;2),(6;1),(6;2) '
      'have phi_a + phi_b = 3 on two cross pairs: 10 doubles on the four cross pairs > 8 = PPE')
print('ALL CHECKS PASSED')
