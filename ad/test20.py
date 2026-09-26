import json
from pair_level import *
from run_all import dd_matching
reps = json.load(open('ad_reps_upto5.json'))
for t, ddc, extra in [(20, [((0, 1), (2, 3), (4, 5))], {6: [(0, 2, 4, 1, 3, 5)]})]:
    rep = reps[t]; dd = [tuple(p) for p in rep['dd']]
    head = {pr(a, h): h for a, h in rep['arcs']}
    In = {c: [x for x in range(N) if x != c and head.get(pr(c, x)) == c] for c in range(N)}
    Ds = {k: [tuple(In[k])] for k in range(N) if len(In[k]) in (2, 3)}; Ds.update(extra)
    r0 = feasible(In, dd, dd_matching(ddc), Ds)
    r1 = feasible(In, dd, dd_matching(ddc), Ds, dichotomy=True)
    print(t, 'plain:', 'feasible' if r0 is not None else 'infeasible', ' with dichotomy:', 'feasible' if r1 is not None else 'infeasible')
