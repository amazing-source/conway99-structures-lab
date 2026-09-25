"""Independent verification (own code, exact integers) of Astra's 23-vertex certificate (message relayed 25/09).

Rooted model at a vertex x of srg(99,14,1,2): labels 0..13, partner a <-> a xor 1; exterior vertices y_ab = {a,b},
a,b non-partners.  Branch (Astra section 4):
  * within the star of label 0 (and of label 1): the prescribed reflection Q of the Z7 witness f_(i,al)(j,be) = (2i-j, be):
    for labels 0 = (0,0) and 1 = (0,1), b = (j,be) -> (-j mod 7, be), i.e. b <-> 14 - b + 2*(b % 2);
  * between the stars of 0 and 1: y_{0,b} ~ y_{1,g(b)}, g = (2 3 4 ... 13);
  * y_{2,4} has exactly the neighbours y_{0,2}, y_{0,5}, y_{1,6}, y_{1,10} in S_0 ∪ S_1.
Checks:
  (A) algebra: N^T N = 11I + J - Pi; for a random symmetric H satisfying nothing, C N computed symbolically is not needed;
      instead we check the identity C N = 0 on the Im N part using the relations H N = 2J - N(I+Pi) only (done by hand
      in the session notes; here we re-check the diagonal C_yy = 20 and the formula of C entries);
  (B) w^T C_23 w with Astra's w  (expected -2);
  (C) the graph-level vector q: 1^T q, |q|^2, q^T A q (expected 0, 251994, 755998), all adjacencies in the support fixed;
  (D) consistency of the branch with the axioms that are checkable on these vertices: the Q-matchings are perfect and
      consistent with the witness at label 2 for the edge y_{2,4} ~ y_{0,2}; y_{2,4} has 2 neighbours containing 0 and 2
      containing 1 (HN rule); lambda/mu conditions among the 23 + x + labels that involve only fixed adjacencies."""
import itertools

L = range(14)
def partner(a): return a ^ 1
def Y(a, b): return (a, b) if a < b else (b, a)

def Q(b):  # reflection in the stars of labels 0 and 1
    return 14 - b + 2 * (b % 2)

assert all(Q(Q(b)) == b and Q(b) != b and 2 <= Q(b) <= 13 for b in range(2, 14))
g = {b: (b + 1 if b < 13 else 2) for b in range(2, 14)}

# fixed exterior-exterior edges of the branch
E = set()
for b in range(2, 14):
    E.add(frozenset((Y(0, b), Y(0, Q(b)))))
    E.add(frozenset((Y(1, b), Y(1, Q(b)))))
    E.add(frozenset((Y(0, b), Y(1, g[b]))))
y24 = Y(2, 4)
for u in (Y(0, 2), Y(0, 5), Y(1, 6), Y(1, 10)):
    E.add(frozenset((y24, u)))

verts = [Y(0, b) for b in range(2, 13)] + [Y(1, b) for b in range(2, 13)] + [y24]
w = [10, 1, -1, 4, -15, -12, -16, -9, -3, -1, 3,
     1, 3, 0, -3, -3, -17, -11, -17, -1, -5, -2,
     15]
assert len(verts) == 23 == len(w)

def H(u, v):
    return 1 if frozenset((u, v)) in E else 0

def Centry(u, v):
    inter = len(set(u) & set(v))
    part = len({partner(a) for a in u} & set(v))
    return 24 * (u == v) - 8 * H(u, v) - 3 * inter - part + 2

# (A) diagonal
assert all(Centry(u, u) == 20 for u in verts)
# N^T N = 11 I + J - Pi on all 84 vertices
allY = [Y(a, b) for a in L for b in L if a < b and b != partner(a)]
assert len(allY) == 84
for a in L:
    for b in L:
        val = sum(1 for y in allY if a in y and b in y)
        assert val == 11 * (a == b) + 1 - (b == partner(a)), (a, b, val)

# (B)
C23 = [[Centry(u, v) for v in verts] for u in verts]
wCw = sum(w[i] * C23[i][j] * w[j] for i in range(23) for j in range(23))
edge_sum = sum(w[i] * w[j] * H(verts[i], verts[j]) for i in range(23) for j in range(i + 1, 23))
print('(B) w^T C_23 w =', wCw, '| sum_{i<j} w_i w_j H_ij =', edge_sum, '| 25086 - 16*sum =', 25086 - 16 * edge_sum)

# (C) graph-level vector
W = sum(w)
t = {a: sum(w[i] for i, y in enumerate(verts) if a in y) for a in L}
q = {'x': -2 * W}
for a in L:
    q[('L', a)] = 3 * t[a] + t[partner(a)] - W
for i, y in enumerate(verts):
    q[y] = 8 * w[i]

def adj(u, v):
    if u == v: return 0
    if u == 'x' or v == 'x':
        o = v if u == 'x' else u
        return 1 if isinstance(o, tuple) and o[0] == 'L' else 0
    if u[0] == 'L' and v[0] == 'L':
        return 1 if v[1] == partner(u[1]) else 0
    if u[0] == 'L' or v[0] == 'L':
        lab_, ext = (u, v) if u[0] == 'L' else (v, u)
        return 1 if lab_[1] in ext else 0
    return H(u, v)

keys = list(q)
one_q = sum(q.values())
norm = sum(val * val for val in q.values())
qAq = sum(q[u] * q[v] * adj(u, v) for u in keys for v in keys)
print('(C) W =', W, '| 1^T q =', one_q, '| |q|^2 =', norm, '| q^T A q =', qAq, '| q^T A q - 3|q|^2 =', qAq - 3 * norm)

# (D) branch consistency checks on fixed data
cnt0 = sum(1 for u in verts if 0 in u and H(y24, u))
cnt1 = sum(1 for u in verts if 1 in u and H(y24, u))
print('(D) y_{2,4}: neighbours containing 0:', cnt0, ', containing 1:', cnt1, '(HN rule requires 2 and 2)')
# witness at label 2 = (1,0): f_2(j,be) = (2-j, be); y_{2,4}: 4 = (2,0) -> (0,0) = label 0: so y_{2,4} ~ y_{2,0} = y_{0,2}
j, be = 4 // 2, 4 % 2
f2 = 2 * ((2 - j) % 7) + be
print('    witness f_2(4) =', f2, '-> edge y_{2,4} ~ y_{0,2} consistent with the witness:', f2 == 0)
# lambda = 1 / mu = 2 among exterior vertices whose common neighbours are all inside the fixed support? Only partial:
# check no fixed edge creates a triangle of exterior vertices that is not a line of the structure (lambda = 1 violation)
bad = 0
V23 = verts
for u, v in itertools.combinations(V23, 2):
    if H(u, v):
        common_ext = [z for z in V23 if z not in (u, v) and H(u, z) and H(v, z)]
        common_lab = [a for a in L if a in u and a in v]
        if len(common_ext) + len(common_lab) > 1:
            bad += 1
print('    fixed adjacent pairs with more than one common neighbour among fixed data (lambda=1 violations):', bad)
