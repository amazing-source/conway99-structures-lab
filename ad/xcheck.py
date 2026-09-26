import json
from general_model import *
from run_all import dd_matching
from pair_level import feasible as feas_old
reps = json.load(open('ad_reps_upto5.json'))
rep = reps[20]; dd = [tuple(p) for p in rep['dd']]
head = {pr(a, h): h for a, h in rep['arcs']}
In = {c: [x for x in range(N) if x != c and head.get(pr(c, x)) == c] for c in range(N)}
Ds = {k: [tuple(In[k])] for k in range(N) if len(In[k]) in (2, 3)}
Ds[6] = [(0, 2, 4, 1, 3, 5)]
ddm = dd_matching([((0, 1), (2, 3), (4, 5))])
orbs = from_ad(In, dd, Ds, ddm)
print('old  no dich:', feas_old(In, dd, ddm, Ds) is not None, ' with dich:', feas_old(In, dd, ddm, Ds, dichotomy=True) is not None)
print('new  no dich:', feasible(orbs, dichotomy=False) is not None, ' with dich:', feasible(orbs) is not None)
# AD2 rep 130 with 2+2+2
from pair_level import parse_survivors
C2 = r"C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/c2exp/logs/ad2_tail_survivors.txt"
R = [r for r in parse_survivors(C2) if r['rep'] == 130][0]; In = R['In']
Ds = {k: [tuple(In[k])] for k in range(N) if len(In[k]) in (2, 3)}; Ds[6] = [(0, 1), (2, 3), (4, 5)]
ddm = {((0, 1), 0): ((2, 3), 0), ((2, 3), 0): ((0, 1), 0), ((0, 1), 1): ((2, 3), 1), ((2, 3), 1): ((0, 1), 1)}
orbs = from_ad(In, [(0, 1), (2, 3)], Ds, ddm)
print('rep130 new no dich:', feasible(orbs, dichotomy=False) is not None, '; no balance, with dich:', feasible(orbs, balance=False) is not None)
