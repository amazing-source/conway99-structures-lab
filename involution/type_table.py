"""Independent recount (structures) of the block types of profile (2,4) at a block c, as in C1's large-T reduction:
pivot labels P = {l : mu(l) = tau(l)} (C2 notation: mu = mates, tau = pi o h); per-block counts |P cap lab(x)| sorted =
(1,1,1,1,2,2); K2, K3 >= 0; 'relevant for Sigma = empty' = no pivot pair {l, mu l} that is a pi-edge.
Classes of (mu, tau) under simultaneous conjugation by the hyperoctahedral group B6 (signed block permutations)."""
import itertools
L = list(range(12)); PI = [v ^ 1 for v in L]


def fpf(pts):
    if not pts:
        yield {}; return
    a = pts[0]
    for i in range(1, len(pts)):
        b = pts[i]; rest = pts[1:i] + pts[i + 1:]
        for m in fpf(rest):
            m = dict(m); m[a] = b; m[b] = a; yield m
ALL = [tuple(m[x] for x in L) for m in fpf(L)]
G = [tuple(2 * perm[i // 2] + ((i % 2) ^ ((fl >> (i // 2)) & 1)) for i in L) for perm in itertools.permutations(range(6)) for fl in range(64)]


def conj(t, g):
    gi = [0] * 12
    for i, x in enumerate(g): gi[x] = i
    return tuple(g[t[gi[x]]] for x in L)


def K_ok(mu, tau):
    for a in L:
        K2 = {b: 2 for b in L}
        for b in (a, PI[a], mu[a], tau[PI[a]]):
            K2[b] -= 1
        if min(K2.values()) < 0: return False
        K3 = {d: 2 for d in L}
        for d in (PI[a], tau[a], tau[mu[a]], mu[tau[a]]):
            K3[d] -= 1
        if min(K3.values()) < 0: return False
    return True
seen_mu, mu_reps = set(), []
for t in ALL:
    if t in seen_mu: continue
    seen_mu |= {conj(t, g) for g in G}; mu_reps.append(t)
import json, time
from collections import Counter
P7 = 7
SQ = {1, 2, 4}


def rank7(M):
    M = [[x % P7 for x in r] for r in M]; r = 0; n = len(M)
    for c in range(len(M[0])):
        p = next((i for i in range(r, n) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; inv = pow(M[r][c], 5, 7)
        M[r] = [x * inv % P7 for x in M[r]]
        for i in range(n):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % P7 for a, b in zip(M[i], M[r])]
        r += 1
    return r


def disc_reduced(M):
    """congruence diagonalisation mod 7 of a symmetric matrix: product of the nonzero pivots (disc of the nondegenerate part)"""
    M = [[x % P7 for x in r] for r in M]; n = len(M); d = 1; alive = list(range(n))
    while True:
        piv = next((i for i in alive if M[i][i]), None)
        if piv is None:
            pair = next(((i, j) for i in alive for j in alive if i < j and M[i][j]), None)
            if pair is None: break
            i, j = pair                                   # row/col i += row/col j  -> M[i][i] = 2 M[i][j] != 0
            for k in range(n): M[i][k] = (M[i][k] + M[j][k]) % P7
            for k in range(n): M[k][i] = (M[k][i] + M[k][j]) % P7
            piv = i
        a = M[piv][piv]; d = d * a % P7; inv = pow(a, 5, 7)
        for i in alive:
            if i != piv and M[i][piv]:
                f = M[i][piv] * inv % P7
                for k in range(n): M[i][k] = (M[i][k] - f * M[piv][k]) % P7
        for i in alive:
            if i != piv: M[i][piv] = 0
        for k in range(n):
            if k != piv: M[piv][k] = 0
        alive.remove(piv)
    return d


def ball(mu, tau):
    C = [[0] * 12 for _ in range(7)]
    for v in L:
        C[0][v] = 1; C[1 + v // 2][v] = 1 if v % 2 == 0 else -1
    H = [[0] * 19 for _ in range(19)]
    for i in range(7): H[i][i] = 4
    for i in range(7):
        for v in L: H[i][7 + v] = H[7 + v][i] = -C[i][v]
    for u in L:
        for v in L:
            H[7 + u][7 + v] = 3 * (u == v) - (mu[u] == v) + (tau[u] == v)
    return H


def union_type(mu, tau):
    seen, lens = set(), []
    for s in L:
        if s in seen: continue
        # walk alternating mu, tau
        cur, n, use_mu = s, 0, True
        while True:
            seen.add(cur); nxt = mu[cur] if use_mu else tau[cur]; n += 1; use_mu = not use_mu; cur = nxt
            if cur == s and use_mu: break
        lens.append(n)
    return tuple(sorted(lens))
t0 = time.time()
rows = []
for mu in mu_reps:
    stab = [g for g in G if conj(mu, g) == mu]
    seen = set()
    for tau in ALL:
        if tau in seen: continue
        orb = {conj(tau, g) for g in stab}; seen |= orb
        if not K_ok(mu, tau): continue
        P = [l for l in L if mu[l] == tau[l]]
        Q = [[(mu[u] == v) - (tau[u] == v) - 2 * (PI[u] == v) - (u == v) for v in L] for u in L]
        t = 12 - rank7(Q)
        H = ball(mu, tau)
        r = rank7(H); assert r == 18 - t, (r, t)
        dsc = disc_reduced(H)
        dim = 18 - t
        split = None if dim % 2 else ((dsc in SQ) == ((dim // 2) % 2 == 0))
        ut = union_type(mu, tau)
        pi_mu = any(mu[l] == PI[l] for l in L); pi_tau = any(tau[l] == PI[l] for l in L)
        rows.append({'pivots': len(P), 't': t, 'split': split, 'union': ut, 'pi_in_mu': pi_mu, 'pi_in_tau': pi_tau,
                     'relevant': all(mu[l] != PI[l] for l in P), 'orbit': len(orb)})
print('classes', len(rows), '%.0fs' % (time.time() - t0))
c = Counter((r['pivots'], r['t'], r['split']) for r in rows)
print('by (pivots, t, split):')
for k in sorted(c, key=lambda x: (x[0], x[1], str(x[2]))): print('  ', k, c[k])
free = [r for r in rows if r['pivots'] == 0]
cu = Counter(r['union'] for r in free)
print('free classes by union cycle type:', dict(cu))
sol = [r for r in free if all(x not in (2, 4) for x in r['union'])]
print('free, union without 2- or 4-cycles:', len(sol), Counter(r['union'] for r in sol), '| of which no pi-edge in mu:',
      sum(1 for r in sol if not r['pi_in_mu']), '| no pi-edge in mu nor tau:', sum(1 for r in sol if not r['pi_in_mu'] and not r['pi_in_tau']))
print('free by (t, split):', dict(Counter((r['t'], r['split']) for r in free)))
json.dump(rows, open('type_table.json', 'w'))
