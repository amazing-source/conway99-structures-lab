"""Test instance for the O_24 u-stage (structures 27/09): for one (y1, y2) class of C1's list, find ONE configuration of 7
signed octads with coverage 1 on A = T1 u T2 and 3 elsewhere, pairwise orthogonal, orthogonal to y1 and y2, even number of
minus signs (h = 2M in O_24).  CP-SAT is used only to produce a test input (C2 enumerates these exhaustively).
Then runs o24_ustage.ustage on it.  Usage: python o24_testinstance.py [class_index] [seconds]"""
import json, sys, time
import numpy as np
from ortools.sat.python import cp_model
from o24_ustage import C1_YPAIRS, golay_from_basis, sign_vectors, ustage

k = int(sys.argv[1]) if len(sys.argv) > 1 else 0
secs = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
D = json.load(open(C1_YPAIRS))
code = golay_from_basis(D['golay_basis'])
OCT = [W for W in code if bin(W).count('1') == 8]
orb = D['orbits'][k]
B1, B2 = np.array(orb['a1']), np.array(orb['a2'])
A = set(orb['T1']) | set(orb['T2'])
cov = [1 if i in A else 3 for i in range(24)]
assert int(B1 @ B2) == 0 and int(B1 @ B1) == 56 == int(B2 @ B2)

md = cp_model.CpModel()
sel = [[md.NewBoolVar('s%d_%d' % (c, o)) for o in range(759)] for c in range(7)]
pos = [[md.NewBoolVar('p%d_%d' % (c, i)) for i in range(24)] for c in range(7)]
neg = [[md.NewBoolVar('n%d_%d' % (c, i)) for i in range(24)] for c in range(7)]
z = [[md.NewIntVar(-1, 1, 'z%d_%d' % (c, i)) for i in range(24)] for c in range(7)]
idx = [md.NewIntVar(0, 758, 'idx%d' % c) for c in range(7)]
for c in range(7):
    md.AddExactlyOne(sel[c])
    md.Add(idx[c] == sum(o * sel[c][o] for o in range(759)))
    for i in range(24):
        md.Add(pos[c][i] + neg[c][i] == sum(sel[c][o] for o in range(759) if OCT[o] >> i & 1))
        md.Add(z[c][i] == pos[c][i] - neg[c][i])
    ev = md.NewIntVar(0, 4, 'ev%d' % c)
    md.Add(sum(neg[c]) == 2 * ev)
    md.Add(sum(int(B1[i]) * z[c][i] for i in range(24)) == 0)
    md.Add(sum(int(B2[i]) * z[c][i] for i in range(24)) == 0)
for c in range(6):
    md.Add(idx[c] <= idx[c + 1])
for i in range(24):
    md.Add(sum(pos[c][i] + neg[c][i] for c in range(7)) == cov[i])
for c in range(7):
    for d in range(c + 1, 7):
        prods = []
        for i in range(24):
            p = md.NewIntVar(-1, 1, 'q%d_%d_%d' % (c, d, i))
            md.AddMultiplicationEquality(p, [z[c][i], z[d][i]])
            prods.append(p)
        md.Add(sum(prods) == 0)
sv = cp_model.CpSolver()
sv.parameters.max_time_in_seconds = secs
sv.parameters.num_workers = 1
t0 = time.time()
st = sv.Solve(md)
print('class %d: CP-SAT status %s in %.1fs' % (k, sv.StatusName(st), time.time() - t0))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    M = np.array([[sv.Value(z[c][i]) for i in range(24)] for c in range(7)])
    assert (M @ M.T == 8 * np.eye(7, dtype=int)).all() and (M @ B1 == 0).all() and (M @ B2 == 0).all()
    assert (np.abs(M).sum(axis=0) == np.array(cov)).all()
    json.dump({'class': k, 'B1': B1.tolist(), 'B2': B2.tolist(), 'M': M.tolist()}, open('state/testinst_%d.json' % k, 'w'))
    SV = sign_vectors(code)
    r = ustage(B1, B2, M, SV)
    print({kk: v for kk, v in r.items() if kk not in ('lp_x',)} if 'proof' not in r else {kk: (v if kk != 'proof' else {'value': v.get('value'), 'den': v.get('den')}) for kk, v in r.items()})
