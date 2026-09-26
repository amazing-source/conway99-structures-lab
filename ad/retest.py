"""Re-test pair-level survivors with the pair dichotomy over ALL D-structure choices."""
import sys, json, itertools
from pair_level import *
from run_all import dd_structures, dd_matching
def retest(rep):
    dd = [tuple(p) for p in rep['dd']]
    head = {pr(t, h): h for t, h in rep['arcs']}
    In = {c: [x for x in range(N) if x != c and head.get(pr(c, x)) == c] for c in range(N)}
    base = {}; free = []
    for k in range(N):
        d = len(In[k])
        if d in (2, 3): base[k] = [tuple(In[k])]
        elif d >= 4: free.append(k)
    opts = [two_regular_graphs(In[k]) for k in free]
    alive = []
    for ddc in dd_structures(dd):
        ddm = dd_matching(ddc)
        for choice in itertools.product(*opts):
            Ds = dict(base); Ds.update({k: c for k, c in zip(free, choice)})
            if feasible(In, dd, ddm, Ds) is None: continue
            r = feasible(In, dd, ddm, Ds, dichotomy=True)
            if r is not None: alive.append((ddc, {k: Ds[k] for k in free}))
    return alive
if __name__ == '__main__':
    reps = json.load(open(sys.argv[1]))
    for t in map(int, sys.argv[2].split(',')):
        al = retest(reps[t])
        print('rep', t, 'survivors with dichotomy:', len(al), al[:3]); sys.stdout.flush()
