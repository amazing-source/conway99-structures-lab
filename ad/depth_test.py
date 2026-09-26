"""How deep must a DFS over blocks go before the relaxed pair model (unassigned pivots = free partners) kills?
Take a class-9 leaf, keep the D-structure only at a subset of pivot blocks, free the others."""
import json, sys, time, itertools
from general_model import build, from_c2_key, feasible
from test_mix import LEAVES
from certify import certify
d = json.load(open(LEAVES)); job = d['class9_job42']
leaf = job['leaves'][int(sys.argv[1]) if len(sys.argv) > 1 else 0]
full = from_c2_key(leaf['triples'])
def partial(keep):
    orbs = []
    for o in full:
        o2 = dict(o)
        if o['pivot'] not in keep: o2['partner'] = None
        orbs.append(o2)
    # partners of kept orbits are kept orbits (same pivot), indices unchanged
    return orbs
for size in (int(sys.argv[2]) if len(sys.argv) > 2 else 1,):
    for keep in itertools.combinations(range(7), size):
        orbs = partial(set(keep))
        t0 = time.time()
        A, lo, hi, lb, ub, *_ = build(orbs)
        r = certify(A, lo, hi, lb, ub, maxnodes=60)
        print('keep', keep, r['status'], r.get('level', ''), r.get('nodes'), '%.1fs' % (time.time() - t0), flush=True)
