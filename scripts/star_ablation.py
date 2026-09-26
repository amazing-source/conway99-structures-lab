"""Ablation of the full-star obstruction at k = 22 (heuristic study, small exact searches): which ingredients are
essential?  Ingredients: lm = lambda/mu against the ball (Bn <= 2 - n); K = no pair of the shadow on a zero of MM^T;
overlap = n.n' <= 2 - a; diag = g(n,n) = 6.  Always kept: shadow shape, pairwise g in {2,3}, full-star capacities."""
import sys, time, random
sys.path.insert(0, 'scripts')
from star_shadows import build, ball_data, star_test, N, part

def rand_inv(rng):
    xs = list(range(12)); rng.shuffle(xs)
    f = [None] * 12
    for i in range(0, 12, 2):
        f[xs[i]], f[xs[i + 1]] = xs[i + 1], xs[i]
    return f

balls = [('Astra', [9,4,7,8,1,6,5,2,3,0,11,10], [7,2,1,9,11,8,10,0,5,3,6,4], [7,3,10,1,11,9,8,0,6,5,2,4], [(i + 1) % 12 for i in range(12)])]
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
nrand = int(sys.argv[2]) if len(sys.argv) > 2 else 3
while len(balls) < 1 + nrand:
    a, b, c = rand_inv(rng), rand_inv(rng), rand_inv(rng)
    hol = list(range(12)); rng.shuffle(hol)
    adj = build(a, b, c, hol)
    Bm, H, tau, rk = ball_data(adj)
    if rk != 35 or tau % 7 in (0, 1, 2, 4):
        continue
    B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
    if any(12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) < 0 for p in range(N) for q in range(N)):
        continue
    balls.append((f'rand{len(balls)}', a, b, c, hol))
variants = [('all', {}), ('-lm', {'use_lm': False}), ('-K', {'use_K': False}), ('-overlap', {'use_overlap': False}),
            ('-diag', {'use_diag': False}), ('-lm-K-overlap', {'use_lm': False, 'use_K': False, 'use_overlap': False})]
budget = float(sys.argv[3]) if len(sys.argv) > 3 else 280
t0 = time.time()
for name, a, b, c, hol in balls:
    adj = build(a, b, c, hol)
    groups = [v for v in range(1, 12) if v != a[0]]
    row = []
    for vname, kw in variants:
        if time.time() - t0 > budget:
            row.append(f'{vname}: (budget)'); continue
        t1 = time.time()
        r = star_test(adj, groups, full=True, **kw)
        row.append(f'{vname}: {"FEASIBLE" if r["solution"] is not None else "dead"} (min cand {min(r["counts"])}, {r["nodes"]} nodes, {time.time() - t1:.0f}s)')
    print(name, ' | '.join(row), flush=True)
