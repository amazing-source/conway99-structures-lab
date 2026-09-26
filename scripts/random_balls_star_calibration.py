"""Calibration (heuristic, NOT a proof; <= 5 min): do random generic line balls with tau a non-square mod 7
(the only ones allowed at k = 22 by L23) all fail the full star test at u = A_0?  Random fixed-point-free
involutions a, b, c and a random holonomy; kept if K = 12I - B + 2J - B^2 - D >= 0 (ball consistent with
lambda/mu), rank_7(3I - B) = 35 and tau non-square."""
import random, sys, time
sys.path.insert(0, 'scripts')
from star_shadows import build, ball_data, star_test, N, part

def rand_inv(rng):
    xs = list(range(12)); rng.shuffle(xs)
    f = [None] * 12
    for i in range(0, 12, 2):
        f[xs[i]], f[xs[i + 1]] = xs[i + 1], xs[i]
    return f

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 2609)
budget = float(sys.argv[2]) if len(sys.argv) > 2 else 270
t0 = time.time(); tried = kept = killed = survived = 0
while time.time() - t0 < budget:
    a, b, c = rand_inv(rng), rand_inv(rng), rand_inv(rng)
    hol = list(range(12)); rng.shuffle(hol)
    adj = build(a, b, c, hol)
    tried += 1
    Bm, H, tau, rk = ball_data(adj)
    if rk != 35 or tau % 7 in (0, 1, 2, 4):
        continue
    B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
    if any(12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) < 0 for p in range(N) for q in range(N)):
        continue
    kept += 1
    res = star_test(adj, [v for v in range(1, 12) if v != a[0]], full=True)
    if res['solution'] is None:
        killed += 1
    else:
        survived += 1
        print('SURVIVOR at u = A_0:', a, b, c, hol, 'tau mod 7 =', tau % 7)
    print(f'ball {kept}: tau mod 7 = {tau % 7}, candidates {res["counts"]}, star solution: {res["solution"] is not None}, nodes {res["nodes"]}', flush=True)
print(f'random balls tried {tried}, admissible generic with tau non-square {kept}: star killed {killed}, star survived {survived} ({time.time() - t0:.0f} s)')
