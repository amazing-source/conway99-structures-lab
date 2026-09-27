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
tot = {'all': 0, 'relevant': 0}
by_mu = []
for mu in mu_reps:
    stab = [g for g in G if conj(mu, g) == mu]
    seen = set(); nall = nrel = 0
    for tau in ALL:
        if tau in seen: continue
        orb = {conj(tau, g) for g in stab}; seen |= orb
        P = [l for l in L if mu[l] == tau[l]]
        counts = sorted(sum(1 for l in P if l // 2 == x) for x in range(6))
        if counts != [1, 1, 1, 1, 2, 2]: continue
        if not K_ok(mu, tau): continue
        nall += 1
        if all(mu[l] != PI[l] for l in P): nrel += 1
    by_mu.append((nall, nrel)); tot['all'] += nall; tot['relevant'] += nrel
print('profile (2,4) classes: all', tot['all'], '| relevant (no pi-edge pivot pair)', tot['relevant'])
print('per mu-class (nonzero):', [x for x in by_mu if x[0]])
