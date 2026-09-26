from single_orbit_proof import *
C2 = r"C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/c2exp/logs/ad2_tail_survivors.txt"
reps = parse_survivors(C2)
dd = [(0, 1), (2, 3)]
ddm = {((0, 1), 0): ((2, 3), 0), ((2, 3), 0): ((0, 1), 0), ((0, 1), 1): ((2, 3), 1), ((2, 3), 1): ((0, 1), 1)}
for R in reps:
    In = R['In']
    Ds = {k: [tuple(In[k])] for k in range(N) if len(In[k]) in (2, 3)}
    best = None
    for P in PAIRS:
        if 6 in P: continue
        for s in (0, 1):
            o = (P, s)
            eqs, lb, ub, why = facts_for_orbit(In, dd, ddm, Ds, o)
            eqs.append(({pr(6, m): 1 for m in range(6)}, 4, 'F1 at block 6'))
            if infeasible(eqs, lb, ub):
                mm = minimal(eqs, lb, ub)
                if best is None or len(mm) < len(best[1]): best = (o, mm, lb, ub, why)
    print('rep', R['rep'], 'In', In)
    show(best[0], best[1], best[2], best[3], best[4])
