"""All-dd + near structure with every G_t CONNECTED (balance vacuous): does the pair model still die?"""
import random, sys, time, numpy as np
from collections import Counter
from general_model import feasible, build, N
from alldd import skeleton_from_pair_cycles, random_cycles
from alldd_near import rand_near
def comps(orbs, near, t):
    par = {m: m for m in range(N) if m != t}
    def f(a):
        while par[a] != a: a = par[a]
        return a
    for i, nb in near[t].items():
        m = sum(orbs[i]['pair']) - t
        for j in nb:
            par[f(m)] = f(sum(orbs[j]['pair']) - t)
    g = {}
    for m in par: g.setdefault(f(m), set()).add(m)
    return sorted(len(s) for s in g.values())
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
stats = Counter(); t0 = time.time(); shapes = Counter()
tries = 0
while sum(stats.values()) < int(sys.argv[2]) if len(sys.argv) > 2 else 6:
    tries += 1
    orbs = skeleton_from_pair_cycles(random_cycles(rng))
    near = rand_near(rng, orbs)
    if near is None: continue
    sh = tuple(tuple(comps(orbs, near, t)) for t in range(N))
    for s in sh: shapes[s] += 1
    if any(s != (6,) for s in sh): continue
    r = feasible(orbs, near=near, time_limit=120)
    stats['dead' if r is None else ('unknown' if isinstance(r, str) else 'ALIVE')] += 1
    print(stats, flush=True)
print('component shapes seen (per block):', dict(shapes.most_common(6)), 'tries', tries, '%.0fs' % (time.time() - t0))
