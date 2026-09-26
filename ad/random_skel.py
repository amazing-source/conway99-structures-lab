"""Random full skeletons (Sigma = empty): a random perfect matching of the 42 orbits in which every matched pair of orbits
either shares exactly one block (one-block square, pivot = that block) or has disjoint supports (dis square); never sisters.
Measures how often the pair-level model (general_model) kills, by number of dis squares."""
import random, sys, time
from collections import Counter
from general_model import feasible, PAIRS, N
def rand_skel(rng, pdis):
    orbs = [(P, s) for P in PAIRS for s in (0, 1)]
    while True:
        rest = list(range(42)); rng.shuffle(rest); partner = {}
        ok = True
        while rest:
            i = rest.pop()
            want_dis = rng.random() < pdis
            cands = [j for j in rest if (not (set(orbs[i][0]) & set(orbs[j][0]))) == want_dis and orbs[j][0] != orbs[i][0]
                     and (want_dis or len(set(orbs[i][0]) & set(orbs[j][0])) == 1)]
            if not cands:
                cands = [j for j in rest if orbs[j][0] != orbs[i][0]]
                if not cands: ok = False; break
            j = rng.choice(cands); rest.remove(j); partner[i] = j; partner[j] = i
        if ok: break
    out = []
    for i, (P, s) in enumerate(orbs):
        Q = orbs[partner[i]][0]
        sh = set(P) & set(Q)
        out.append({'pair': P, 'pivot': (sh.pop() if sh else None), 'partner': partner[i]})
    return out
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
trials = int(sys.argv[2]) if len(sys.argv) > 2 else 40
stats = Counter(); t0 = time.time()
for t in range(trials):
    pdis = rng.choice([0.05, 0.15, 0.3])
    sk = rand_skel(rng, pdis)
    ndis = sum(1 for o in sk if o['pivot'] is None) // 2
    r = feasible(sk, time_limit=30)
    st = 'dead' if r is None else ('unknown' if isinstance(r, str) else 'ALIVE')
    stats[(min(ndis, 6), st)] += 1
print(dict(sorted(stats.items())), '%.0fs' % (time.time() - t0))
