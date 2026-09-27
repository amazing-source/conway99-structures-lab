"""All 3308 block-type classes (mu, tau) under B6 with K2, K3 >= 0, in the SAME order as type_table.py / type_table.json
(same enumeration loop), with representatives: writes all_types_3308.json for the finish line (C1/C2, 27/09)."""
import itertools, json, time
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
        for b in (a, PI[a], mu[a], tau[PI[a]]): K2[b] -= 1
        if min(K2.values()) < 0: return False
        K3 = {d: 2 for d in L}
        for d in (PI[a], tau[a], tau[mu[a]], mu[tau[a]]): K3[d] -= 1
        if min(K3.values()) < 0: return False
    return True
t0 = time.time()
seen_mu, mu_reps = set(), []
for t in ALL:
    if t in seen_mu: continue
    seen_mu |= {conj(t, g) for g in G}; mu_reps.append(t)
rows = []
for mu in mu_reps:
    stab = [g for g in G if conj(mu, g) == mu]
    seen = set()
    for tau in ALL:
        if tau in seen: continue
        orb = {conj(tau, g) for g in stab}; seen |= orb
        if not K_ok(mu, tau): continue
        P = [l for l in L if mu[l] == tau[l]]
        rows.append({'id': len(rows), 'mu': list(mu), 'tau': list(tau), 'pivots': len(P),
                     'relevant': all(mu[l] != PI[l] for l in P), 'orbit': len(orb), 'stab_order': 46080 // len(orb) if False else None})
print('classes', len(rows), '%.0fs' % (time.time() - t0))
tt = json.load(open('type_table.json'))
assert len(tt) == len(rows) and all(a['pivots'] == b['pivots'] and a['orbit'] == b['orbit'] and a['relevant'] == b['relevant'] for a, b in zip(tt, rows)), 'order mismatch with type_table.json'
for r in rows:
    del r['stab_order']
json.dump(rows, open('all_types_3308.json', 'w'))
from collections import Counter
print('by pivots:', sorted(Counter(r['pivots'] for r in rows).items()), '| matches type_table.json order')
