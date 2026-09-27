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
from collections import Counter
cnt_all, cnt_rel = Counter(), Counter()
for mu in mu_reps:
    stab = [g for g in G if conj(mu, g) == mu]
    seen = set()
    for tau in ALL:
        if tau in seen: continue
        orb = {conj(tau, g) for g in stab}; seen |= orb
        if not K_ok(mu, tau): continue
        P = [l for l in L if mu[l] == tau[l]]
        cnt_all[len(P)] += 1
        if all(mu[l] != PI[l] for l in P): cnt_rel[len(P)] += 1
print('classes (K >= 0) by number of pivot labels:', dict(sorted(cnt_all.items())))
print('  of which relevant for Sigma = empty (no pi-edge pivot pair):', dict(sorted(cnt_rel.items())))
print('  total with pivots:', sum(v for k, v in cnt_all.items() if k > 0), '| >= 6 pivots:', sum(v for k, v in cnt_all.items() if k >= 6),
      '| relevant >= 6:', sum(v for k, v in cnt_rel.items() if k >= 6))
