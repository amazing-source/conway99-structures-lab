"""Exact check of C1's claim: if G(L) has c components then rank_Q(A_B - 3I) = 38 - c and d(L) = 37 - rank_7(Q_BB) >= c,
with Q_BB = A_B - 3I + 3J on the 39-vertex ball T + N1(T).  Balls: Astra's (connected) and group-preserving ones."""
from fractions import Fraction

def ball(a, b, c, hol):
    n = 39                                               # 0,1,2 = T ; 3 + 12t + i = part t, index i
    A = [[0] * n for _ in range(n)]
    def e(p, q):
        A[p][q] = A[q][p] = 1
    e(0, 1); e(0, 2); e(1, 2)
    for t, f in enumerate([a, b, c]):
        for i in range(12):
            e(t, 3 + 12 * t + i)
            if f[i] > i:
                e(3 + 12 * t + i, 3 + 12 * t + f[i])
    for i in range(12):
        e(3 + i, 3 + 12 + i); e(3 + i, 3 + 24 + i); e(3 + 12 + i, 3 + 24 + hol[i])
    return A

def comps(a, b, c, hol):
    adj = {(t, i): set() for t in range(3) for i in range(12)}
    for t, f in enumerate([a, b, c]):
        for i in range(12):
            adj[t, i].add((t, f[i]))
    for i in range(12):
        for p, q in [((0, i), (1, i)), ((0, i), (2, i)), ((1, i), (2, hol[i]))]:
            adj[p].add(q); adj[q].add(p)
    seen, k = set(), 0
    for s in adj:
        if s not in seen:
            k += 1; st = [s]; seen.add(s)
            while st:
                p = st.pop()
                for q in adj[p]:
                    if q not in seen:
                        seen.add(q); st.append(q)
    return k

def rank_Q(M):
    M = [[Fraction(x) for x in r] for r in M]; r = 0
    for col in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][col] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][col] != 0:
                f = M[i][col] / M[r][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return r

def rank_p(M, p=7):
    M = [[x % p for x in r] for r in M]; r = 0
    for col in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][col]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][col], -1, p); M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col]:
                f = M[i][col]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r

def report(name, a, b, c, hol):
    A = ball(a, b, c, hol)
    n = len(A)
    A3 = [[A[i][j] - 3 * (i == j) for j in range(n)] for i in range(n)]
    Qbb = [[A3[i][j] + 3 for j in range(n)] for i in range(n)]
    cc = comps(a, b, c, hol)
    rq, r7 = rank_Q(A3), rank_p(Qbb)
    print(f'{name}: c = {cc}; rank_Q(A_B - 3I) = {rq} (38 - c = {38 - cc}); rank_7(Q_BB) = {r7}; d(L) = 37 - rank = {37 - r7}; d >= c: {37 - r7 >= cc}')

astra = ([9, 4, 7, 8, 1, 6, 5, 2, 3, 0, 11, 10], [7, 2, 1, 9, 11, 8, 10, 0, 5, 3, 6, 4], [7, 3, 10, 1, 11, 9, 8, 0, 6, 5, 2, 4],
         [(i + 1) % 12 for i in range(12)])
report("Astra's ball", *astra)
def grp(perm_in_group, size):
    out = []
    for g0 in range(0, 12, size):
        out += [g0 + perm_in_group[j] for j in range(size)]
    return out
# (4,4,4): three groups of 4 indices, all matchings and the holonomy preserve the groups
report('(4,4,4) ball', grp([1, 0, 3, 2], 4), grp([2, 3, 0, 1], 4), grp([3, 2, 1, 0], 4), grp([1, 2, 3, 0], 4))
# (6,6): two groups of 6
report('(6,6) ball', grp([1, 0, 3, 2, 5, 4], 6), grp([2, 4, 0, 5, 1, 3], 6), grp([3, 5, 4, 0, 2, 1], 6), grp([1, 2, 3, 4, 5, 0], 6))
# (4,8)
a48 = [1, 0, 3, 2] + [4 + x for x in [1, 0, 3, 2, 5, 4, 7, 6]]
b48 = [2, 3, 0, 1] + [4 + x for x in [2, 4, 0, 6, 1, 7, 3, 5]]
c48 = [3, 2, 1, 0] + [4 + x for x in [3, 5, 7, 0, 6, 1, 4, 2]]
h48 = [1, 2, 3, 0] + [4 + x for x in [1, 2, 3, 4, 5, 6, 7, 0]]
report('(4,8) ball', a48, b48, c48, h48)
