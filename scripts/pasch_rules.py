"""Human-readable rules per orbit for the PASCH hand proof: roles of the triangles (C0/C1/C2/C3/CL), the pairs kept after
immediate exclusion with their immediate double bound and the reason, and the surviving cycles with per-orbit bounds."""
import itertools, sys
sys.path.insert(0, 'scripts')
from pasch_greedy2 import block_data2
from pasch_cheapest import PASCH
from zeros_budget import V, CYC6

tris = [tuple(sorted(c)) for c in PASCH if len(c) == 3]
extras = {3: {(5, frozenset((1, 6))): 2}, 1: {}, 0: {(3, frozenset((1, 4))): 2, (5, frozenset((1, 6))): 2}}
for j in (3, 1, 0):
    cyc_of_edge = {}
    for c in PASCH:
        for t in range(len(c)):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % len(c)]))] = c
    others, exact, simple, imm = block_data2(PASCH, j)
    print(f'================ block {j}')
    for i in others:
        c = cyc_of_edge[frozenset((i, j))]; m = len(c); t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]; k = b if a == j else a
        R = [r for r in V if r not in (i, j, k)]
        roles = []
        if len(c) == 3: roles.append(f'C0: not {k}')
        for tr in tris:
            s = set(tr)
            if {i, j} <= s: continue
            if i in s: roles.append(f'C1 {tr}: phi{sorted(s - {i})} = 1')
            elif j in s and k in s: roles.append(f'C3 {tr}: not both {k},{sorted(s - {j, k})[0]}')
            elif j in s: roles.append(f'C2 {tr}: at most one of {sorted(s - {j})}, and then phi{sorted(s - {j})} = 1 - N')
            elif k in s: roles.append(f'CL {tr}: phi{sorted(set(R) - s)} = 3 - N{sorted(set(R) - s)} (so N >= 1)')
            else: roles.append(f'C5 {tr}')
        kept = [(sorted(pr), simple[(ii, pr)], exact[(ii, pr)]) for (ii, pr) in exact if ii == i and (ii, pr) not in imm]
        print(f' orbit {i}: k={k}, own {c}, R={R}; ' + '; '.join(roles))
        print(f'    kept pairs (immediate bound / exact): ' + ', '.join(f'{p}:{s}/{e}' for p, s, e in sorted(kept)))
    val = {key: (None if key in imm else simple[key]) for key in exact}
    val.update(extras[j])
    seen = set()
    for order in CYC6:
        lab = [others[p] for p in order]
        parts, tot, ok = [], 0, True
        for t in range(6):
            key = (lab[t], frozenset((lab[t - 1], lab[(t + 1) % 6])))
            v = val[key]
            if v is None: ok = False; break
            tot += v; parts.append(f'{lab[t]}{sorted(key[1])}={v}')
        if ok:
            print(f'   cycle {lab}: bound {tot}  [' + ', '.join(parts) + ']')
