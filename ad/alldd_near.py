"""All-dd skeleton + a random NEAR structure (mu_t, tau_t at every block t: R_t on the 12 orbits through t, no digon, no sister
edge, sisters with <= 1 common R_t-neighbour), then the pair model with exact near values and exact components.
Question: how often does fixing the near structure make the pair level infeasible?"""
import random, sys, time
from collections import Counter
from general_model import feasible, PAIRS, pr, N
from alldd import skeleton_from_pair_cycles, random_cycles
def rand_near(rng, orbs):
    near = {}
    for t in range(N):
        through = [i for i, o in enumerate(orbs) if t in o['pair']]
        sis = {}
        for i in through:
            for j in through:
                if i != j and orbs[i]['pair'] == orbs[j]['pair']: sis[i] = j
        for _ in range(2000):
            xs = through[:]; rng.shuffle(xs)
            m = {}; ok = True
            for a in range(0, 12, 2): m[xs[a]] = xs[a + 1]; m[xs[a + 1]] = xs[a]
            ys = through[:]; rng.shuffle(ys); p = {}
            for a in range(0, 12, 2): p[ys[a]] = ys[a + 1]; p[ys[a + 1]] = ys[a]
            if any(m[i] == p[i] or m[i] == sis[i] or p[i] == sis[i] for i in through): continue
            # 11.5: sisters have at most one common R_t neighbour
            if any(len({m[i], p[i]} & {m[sis[i]], p[sis[i]]}) > 1 for i in through): continue
            near[t] = {i: (m[i], p[i]) for i in through}
            break
        else:
            return None
    return near
if __name__ == '__main__':
  rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
  cnt = Counter(); t0 = time.time()
  for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 10):
      orbs = skeleton_from_pair_cycles(random_cycles(rng))
      near = rand_near(rng, orbs)
      if near is None: cnt['nogen'] += 1; continue
      r = feasible(orbs, near=near, time_limit=60)
      cnt['dead' if r is None else ('unknown' if isinstance(r, str) else 'ALIVE')] += 1
  print(dict(cnt), '%.0fs' % (time.time() - t0))
