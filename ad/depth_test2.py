import json, sys, time, itertools
from general_model import feasible, from_c2_key
from test_mix import LEAVES
d = json.load(open(LEAVES)); job = d['class9_job42']
leaf = job['leaves'][0]
full = from_c2_key(leaf['triples'])
def partial(keep):
    return [dict(o, partner=(o['partner'] if o['pivot'] in keep else None)) for o in full]
size = int(sys.argv[1])
res = {}
for keep in itertools.combinations(range(7), size):
    t0 = time.time()
    r = feasible(partial(set(keep)), time_limit=60)
    st = 'dead' if r is None else ('unknown' if isinstance(r, str) else 'alive')
    print('keep', keep, st, '%.1fs' % (time.time() - t0), flush=True)
