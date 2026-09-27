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
from collections import Counter, defaultdict
t0 = time.time()


def canon(mu, tau):
    return min((conj(mu, g), conj(tau, g)) for g in G)
mine = {}
for mu in mu_reps:
    stab = [g for g in G if conj(mu, g) == mu]
    seen = set()
    for tau in ALL:
        if tau in seen: continue
        orb = {conj(tau, g) for g in stab}; seen |= orb
        if not K_ok(mu, tau): continue
        P = [l for l in L if mu[l] == tau[l]]
        if len(P) < 6: continue
        per = [sum(1 for l in P if l // 2 == x) for x in range(6)]
        prof = (per.count(2), per.count(1))
        rel = all(mu[l] != PI[l] for l in P)
        mine[canon(mu, tau)] = (prof, rel)
print('mine: %d classes with >= 6 pivots (%.0fs)' % (len(mine), time.time() - t0), flush=True)
C1 = json.load(open('C:/home/researcher/Desktop/c1_recherche/data/all_ge6.json'))
theirs = {}
for d in C1:
    mu, tau = tuple(d['mu']), tuple(d['tau'])
    theirs[canon(mu, tau)] = (tuple(d['anch_tw']), all(mu[l] != PI[l] for l in d['pivots']))
print('C1 : %d classes (%.0fs)' % (len(theirs), time.time() - t0), flush=True)
print('same set of classes:', set(mine) == set(theirs), '| only mine:', len(set(mine) - set(theirs)), '| only C1:', len(set(theirs) - set(mine)))
print('same profile + relevance for every class:', all(mine[k] == theirs[k] for k in mine if k in theirs))
cm = Counter(v[0] for v in mine.values()); cr = Counter(v[0] for v in mine.values() if v[1])
print('per profile (anch, tw): all', dict(sorted(cm.items(), reverse=True)))
print('                   relevant', dict(sorted(cr.items(), reverse=True)), '| total relevant', sum(cr.values()))
# the relevant files of C1 by tier
for fn in ('largeT_tier6_rel', 'largeT_tier12_rel', 'largeT_tier8_rel', 'tier_32_rel', 'tier_40_rel', 'allonblock_ge6_rel'):
    F = json.load(open('C:/home/researcher/Desktop/c1_recherche/data/%s.json' % fn))
    ks = [canon(tuple(d['mu']), tuple(d['tau'])) for d in F]
    ok = all(k in mine and mine[k][1] for k in ks) and len(set(ks)) == len(ks)
    profs = Counter(mine[k][0] for k in ks if k in mine)
    print('%-20s %3d classes, all relevant & distinct & in my list: %s, profiles %s' % (fn, len(F), ok, dict(profs)), flush=True)
