"""Ablation 'only_green' (lambda/mu against the ball, K-zeros and overlaps all dropped) on other balls:
python star_ablation_ball.py cyclic241 | randSEED"""
import sys, time, random
sys.path.insert(0, 'scripts')
from star_shadows import build, ball_data, star_test, N, part
from verify_astra_d5_cyclic import ball as cyc_ball
def rand_inv(rng):
    xs = list(range(12)); rng.shuffle(xs)
    f = [None] * 12
    for i in range(0, 12, 2):
        f[xs[i]], f[xs[i + 1]] = xs[i + 1], xs[i]
    return f
arg = sys.argv[1]
if arg == 'cyclic241':
    adj = cyc_ball(2, 4, 1); partner0 = 6
else:
    rng = random.Random(int(arg[4:]))
    while True:
        a, b, c = rand_inv(rng), rand_inv(rng), rand_inv(rng)
        hol = list(range(12)); rng.shuffle(hol)
        adj = build(a, b, c, hol)
        Bm, H, tau, rk = ball_data(adj)
        if rk != 35 or tau % 7 in (0, 1, 2, 4):
            continue
        B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
        if any(12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) < 0 for p in range(N) for q in range(N)):
            continue
        partner0 = a[0]; break
t0 = time.time()
r = star_test(adj, [v for v in range(1, 12) if v != partner0], full=True, use_lm=False, use_K=False, use_overlap=False)
print(f'{arg} only_green: {"FEASIBLE" if r["solution"] is not None else "dead"}; candidates {r["counts"]}; nodes {r["nodes"]}; {time.time() - t0:.0f}s', flush=True)
