"""'Step 1' at the full block of an AD residue (analogue of C1's 3K2 step 1): in the pair-level model with the D-structure of
the full block c fixed, look at phi_W(c b) for the far orbits W, and test whether the RIGID pattern is forced.
Fact (NOTE_S29): a far orbit W of a pair e1 e2 anchored at e1 is rigid at c iff phi_W(c e2) = 0 and phi_W(c b) = 1 for the
four blocks b outside {c, e1, e2} (shadow with one label in each of the 4 other blocks; phi_W(c e1) = 0 always).
Usage: python step1_fullblock.py reps.json rep_index c 'cycles' [test]"""
import sys, json, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from pair_level import N, PAIRS, pr
from pair_level_relaxed import build_relaxed

reps = json.load(open(sys.argv[1]))
t = int(sys.argv[2]); c = int(sys.argv[3]); cycles = [tuple(x) for x in json.loads(sys.argv[4])]
test = len(sys.argv) > 5 and sys.argv[5] == 'test'
rep = reps[t]; dd = [tuple(p) for p in rep['dd']]
head = {pr(a, h): h for a, h in rep['arcs']}
In = {k: [x for x in range(N) if x != k and head.get(pr(k, x)) == k] for k in range(N)}
assert len(In[c]) == 6, In[c]
A, lo, hi, lb, ub, nvar = build_relaxed(In, dd, True, {c: cycles})
orbits = [(p, s) for p in PAIRS for s in (0, 1)]
cnt = 0
phi = {}
for o in orbits:
    for P in PAIRS:
        phi[(o, P)] = cnt; cnt += 1


def solve(lb2, ub2, tl=60):
    res = milp(np.zeros(nvar), constraints=[LinearConstraint(A, lo, hi)], integrality=np.ones(nvar),
               bounds=Bounds(lb2, ub2), options={'time_limit': tl})
    return res
res = solve(lb, ub)
print('rep', t, 'full block', c, 'cycles', cycles, 'status', res.status, flush=True)
if res.status != 0:
    sys.exit()
far = [P for P in PAIRS if c not in P]
anch = {P: head.get(P) if P not in dd else None for P in far}
x = np.round(res.x).astype(int)
for P in far:
    for s in (0, 1):
        o = (P, s)
        prof = [int(x[phi[(o, pr(c, b))]]) if b != c else '-' for b in range(N)]
        print('  far orbit', P, s, 'anchor', anch[P], 'phi(c b), b = 0..6:', prof)
if test:
    t0 = time.time()
    for P in far:
        if anch[P] is None:
            continue
        e1 = anch[P]; e2 = P[0] + P[1] - e1
        for s in (0, 1):
            o = (P, s)
            verdicts = []
            # can phi(c e2) be >= 1 ?
            lb2, ub2 = lb.copy(), ub.copy(); lb2[phi[(o, pr(c, e2))]] = 1
            r = solve(lb2, ub2); verdicts.append(('e2>=1', {0: 'feas', 2: 'INF'}.get(r.status, 'unk')))
            for b in range(N):
                if b in (c, e1, e2):
                    continue
                for val in (0, 2):
                    lb2, ub2 = lb.copy(), ub.copy()
                    j = phi[(o, pr(c, b))]
                    if val == 0:
                        ub2[j] = 0
                    else:
                        lb2[j] = 2
                    r = solve(lb2, ub2); verdicts.append(('b%d=%d' % (b, val), {0: 'feas', 2: 'INF'}.get(r.status, 'unk')))
            print('  orbit', P, s, verdicts, '%.0fs' % (time.time() - t0), flush=True)
