import sys, time
from pair_level import *
C2 = r"C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/c2exp/logs/ad2_tail_survivors.txt"
reps = parse_survivors(C2)
dd = [(0, 1), (2, 3)]
ddm = {((0, 1), 0): ((2, 3), 0), ((2, 3), 0): ((0, 1), 0), ((0, 1), 1): ((2, 3), 1), ((2, 3), 1): ((0, 1), 1)}
print(len(two_regular_graphs(range(6))), 'D-structures on 6 arms')
for R in reps:
    In = R['In']; t0 = time.time()
    base = {}
    Aopts = None
    for k in range(N):
        d = len(In[k])
        if d in (2, 3): base[k] = [tuple(In[k])]
        elif d == 6: A = k; Aopts = two_regular_graphs(In[k])
        elif d != 0: raise ValueError(d)
    alive = []
    for cyc in (Aopts or [None]):
        Ds = dict(base)
        if cyc is not None: Ds[A] = cyc
        r = feasible(In, dd, ddm, Ds)
        if r is not None: alive.append((cyc, r if isinstance(r, str) else 'feasible'))
    print('rep', R['rep'], 'A-structures alive:', len(alive), 'of', len(Aopts or [0]), ' %.1fs' % (time.time() - t0))
    for cyc, st in alive[:4]: print('   ', cyc, st)
    sys.stdout.flush()
