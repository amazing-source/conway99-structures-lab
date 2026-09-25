"""Exact check of C1's certificate 3 (explicit data sent 25/09): triangle of classes {0,1,6}, T with 10 edges.
(a) lambda_max(T - (3I - P_d)) > 0 exactly (Sturm count of positive roots of the characteristic polynomial);
(b) local admissibility inside W = {x} ∪ Γ(x) ∪ Q: adjacent pairs have <= 1 common neighbour in W, non-adjacent <= 2."""
import itertools
import sympy as sp

labs = [(c, s) for c in (0, 1, 6) for s in '+-']
pairs = [((0, '+'), (1, '+')), ((0, '+'), (1, '-')), ((0, '-'), (1, '+')), ((0, '-'), (1, '-')),
         ((0, '+'), (6, '+')), ((0, '+'), (6, '-')), ((0, '-'), (6, '+')), ((0, '-'), (6, '-')),
         ((1, '+'), (6, '+')), ((1, '+'), (6, '-')), ((1, '-'), (6, '+')), ((1, '-'), (6, '-'))]
P = lambda s: tuple(tuple((int(t[0]), t[1])) for t in s)
def pr(a, b): return (a, b)
E = [(((0,'+'),(1,'+')), ((0,'+'),(6,'+'))), (((0,'+'),(1,'+')), ((1,'+'),(6,'+'))),
     (((0,'+'),(1,'-')), ((0,'-'),(6,'+'))), (((0,'+'),(1,'-')), ((1,'+'),(6,'+'))),
     (((0,'-'),(1,'-')), ((0,'+'),(6,'+'))), (((0,'-'),(1,'-')), ((1,'+'),(6,'-'))),
     (((0,'-'),(1,'-')), ((1,'-'),(6,'+'))), (((0,'+'),(6,'-')), ((1,'-'),(6,'+'))),
     (((0,'-'),(6,'+')), ((1,'+'),(6,'-'))), (((0,'-'),(6,'-')), ((1,'+'),(6,'+')))]
idx = {p: i for i, p in enumerate(pairs)}
T = sp.zeros(12, 12)
for a, b in E:
    T[idx[a], idx[b]] = T[idx[b], idx[a]] = 1
d = sp.Matrix([[1 if l in p else 0 for p in pairs] for l in labs])
Pd = d.T * (d * d.T).inv() * d
M = T - (3 * sp.eye(12) - Pd)
x = sp.symbols('x')
cp = sp.Poly(M.charpoly(x).as_expr(), x)
npos = sp.count_roots(cp, 0, sp.oo) - (1 if cp.eval(0) == 0 else 0)
print('(a) degree', cp.degree(), '; number of roots > 0 (exact, Sturm):', npos,
      '; largest root ~', max(sp.Poly(cp).nroots(), key=lambda r: sp.re(r)))
# (b) local admissibility in W, with all 14 labels (classes 0..6)
all_labs = [(c, s) for c in range(7) for s in '+-']
W = ['x'] + all_labs + pairs
def adj(u, v):
    if u == v: return 0
    if u == 'x' or v == 'x':
        o = v if u == 'x' else u
        return 1 if o in all_labs else 0
    ul, vl = u in all_labs, v in all_labs
    if ul and vl: return 1 if u[0] == v[0] else 0
    if ul or vl:
        l, p = (u, v) if ul else (v, u)
        return 1 if l in p else 0
    return T[idx[u], idx[v]]
bad = 0
for u, v in itertools.combinations(W, 2):
    common = sum(1 for z in W if z not in (u, v) and adj(u, z) and adj(v, z))
    if adj(u, v) and common > 1: bad += 1
    if not adj(u, v) and common > 2: bad += 1
print('(b) violations of lambda <= 1 / mu <= 2 inside W:', bad)
