import random, sys, time
from star_lemma import feasible_star
from general_model import feasible
from alldd import skeleton_from_pair_cycles, random_cycles
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
mode = sys.argv[3] if len(sys.argv) > 3 else 'lp'
for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3):
    orbs = skeleton_from_pair_cycles(random_cycles(rng))
    t0 = time.time()
    r = feasible_star(orbs, lp=(mode == 'lp'), time_limit=300)
    print('all-dd', t, 'with (*) [%s]:' % mode, 'DEAD' if r is None else r, '%.0fs' % (time.time() - t0), flush=True)
