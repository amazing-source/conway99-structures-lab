"""Before exporting the 156 classes (>= 6 pivots) without SEC: the flags --cuts/--digon/--path and --twist use HARD-CODED
right-hand sides (1 - [digon] for anchored non-sib Sp, 3 for twisted S3).  Recompute each from C1's general identity
slack(S, p) = 20 - 10|S| + Lambda(S) + C_p(S) (sat_local2.ident_list machinery) on every class, and list the flags that
actually add clauses.  Output: lemma_scope_156.json (C2's all_ge6.json order)."""
import json, itertools
import iota_design as D
import sat_local2 as S2
T = json.load(open('C:/home/researcher/Desktop/c2_experimental_ro/scripts/c1_tiers/data/all_ge6.json'))
rel = json.load(open('C:/home/researcher/Desktop/c1_recherche/data/allonblock_ge6_rel.json'))
relkey = {(tuple(d['mu']), tuple(d['tau'])): i for i, d in enumerate(rel)}
L = range(12)
def slack(mu, tau, S, p):
    K2, K3 = D.K2K3(mu, tau)
    S = sorted(S)
    lam = sum((mu[b] != d) + int(K3[b, d]) for b, d in itertools.combinations(S, 2)) + sum(int(K3[b, b]) // 2 for b in S)
    return 20 - 10 * len(S) + lam + sum(int(K2[2 * p, b]) + int(K2[2 * p + 1, b]) for b in S)
out, bad = [], []
for i, t in enumerate(T):
    mu, tau = t['mu'], t['tau']
    piedge = any(mu[l] == l ^ 1 or tau[l] == l ^ 1 for l in L)
    flags, notes = [], []
    dg, paths = S2.anchored_structures(mu, tau)
    anch = [p for p in range(6) if mu[2 * p] == tau[2 * p] and mu[2 * p + 1] == tau[2 * p + 1] and mu[2 * p] != 2 * p + 1]
    for p in anch:
        a, b = 2 * p, 2 * p + 1
        digon = mu[a] // 2 == mu[b] // 2
        s = slack(mu, tau, {a, b, mu[a], mu[b]}, p)
        if s != 1 - digon:
            bad.append((i, 'anchored Sp', p, s, 1 - digon))
    if anch:
        flags.append('--cuts')
    if dg:
        flags.append('--digon')
    if paths:
        flags.append('--path')
    tw = []
    for p in range(6):
        piv = [l for l in (2 * p, 2 * p + 1) if mu[l] == tau[l]]
        if len(piv) != 1:
            continue
        ps = piv[0]; po = ps ^ 1; qs = mu[ps]
        if qs // 2 == p:
            continue
        s = slack(mu, tau, {ps, po, qs}, p)
        tw.append(p)
        if s != 3:
            bad.append((i, 'twist S3', p, s, 3))
    if tw:
        flags.append('--twist')
    ids = S2.ident_list(mu, tau, 1)
    extra = [x for x in ids if not (x[0] in anch and set(x[1]) == {2 * x[0], 2 * x[0] + 1, mu[2 * x[0]], mu[2 * x[0] + 1]})]
    if extra:
        flags.append('--ident=1')
    out.append({'i': i, 'anch_tw': t['anch_tw'], 'pivots': t.get('n_pivot_labels'), 'relevant_index': relkey.get((tuple(mu), tuple(tau))),
                'pi_edge': piedge, 'flags_general': flags, 'n_new_idents': len(extra), 'twisted_blocks': tw, 'anchored_nonsib': anch})
json.dump(out, open('lemma_scope_156.json', 'w'), indent=1)
print('classes', len(out), '| relevant', sum(o['relevant_index'] is not None for o in out), '| pi-edge', sum(o['pi_edge'] for o in out))
print('hard-coded RHS mismatches:', len(bad), bad[:10])
from collections import Counter
print('flag sets:', Counter(' '.join(o['flags_general']) for o in out).most_common())
print('non-relevant without pi-edge:', [o['i'] for o in out if o['relevant_index'] is None and not o['pi_edge']])
