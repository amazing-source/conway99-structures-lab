"""Heuristic check: at a full block (6-cycle), are (i)+(ii) alone and (i)+(iii) alone easy (as the count predicts: ~1e24 and
~1e48 solutions), while (i)+(ii)+(iii) is not found?  Same CP-SAT model as fullblock_cpsat.py with parts switched off."""
import itertools, sys, time
from ortools.sat.python import cp_model
from fullblock_rigid import mates_from_D, K2K3, far_roots, pi, L

which = sys.argv[1]; tl = float(sys.argv[2])
cycles = [(0, 1, 2, 3, 4, 5)]
m = mates_from_D(cycles)
K2, K3 = K2K3(m)
far, reps = far_roots()
roots = [z for Z in reps for z in Z]
nonm = [e for e in itertools.combinations(L, 2) if m[e[0]] != e[1]]
md = cp_model.CpModel()
x = {(i, e): md.NewBoolVar('') for i in range(60) for e in nonm}
for i in range(60):
    md.AddExactlyOne(x[i, e] for e in nonm)
for e in nonm:
    md.AddExactlyOne(x[i, e] for i in range(60))
U = {}
for i in range(60):
    for b in L:
        U[i, b] = md.NewBoolVar('')
        md.Add(U[i, b] == sum(x[i, e] for e in nonm if b in e))
if 'ii' in which:
    for a in L:
        for b in L:
            md.Add(sum(U[i, b] for i in range(60) if a in roots[i]) == int(K2[a, b]))
if 'iii' in which:
    W = {}
    for Zi in range(30):
        i, j = 2 * Zi, 2 * Zi + 1
        for b in L:
            for d in L:
                w = md.NewBoolVar('')
                md.AddBoolAnd([U[i, b], U[j, d]]).OnlyEnforceIf(w)
                md.AddBoolOr([U[i, b].Not(), U[j, d].Not()]).OnlyEnforceIf(w.Not())
                W[Zi, b, d] = w
    for b in L:
        for d in L:
            md.Add(sum(W[Zi, b, d] + W[Zi, d, b] for Zi in range(30)) == int(K3[b, d]))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 2
t0 = time.time(); st = s.Solve(md)
print(which, s.StatusName(st), '%.1fs' % (time.time() - t0))
