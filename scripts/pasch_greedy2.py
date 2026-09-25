"""pasch_greedy with the Complement Lemma (CL) added to the immediate facts:
  CL: if k lies on a triangle {k,c,d} with c,d in R_x, then for {a,b} = R_x - {c,d}: phi(ab) = 3 - N_ab.
Immediate exclusions: C0, C2/C3 violations, and any contradiction between forced values of the same phi
(C1: 1 ; C2: 1 - N ; CL: 3 - N ; values must lie in {0,1,2}).
Immediate doubles bound: max(#phi forced to 2, n_k + #phi forced to 0)."""
import itertools, sys
sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V, CYC6
from pasch_cheapest import cycle_min, PASCH

def block_data2(dec, j):
    cyc_of_edge = {}
    for c in dec:
        for t in range(len(c)):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % len(c)]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]
    others = [r for r in V if r != j]
    exact, simple, imm = {}, {}, set()
    for i in others:
        c = cyc_of_edge[frozenset((i, j))]; m = len(c); t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]
        k = b if a == j else a
        R = [r for r in V if r not in (i, j, k)]
        own_tri = len(c) == 3
        for l, r in itertools.combinations([o for o in others if o != i], 2):
            pr = frozenset((l, r)); key = (i, pr)
            exact[key] = d3.local_min(d3.canon(tris, i, j, k, R, pr))
            n = {x: (1 if x in pr else 0) for x in V}
            bad = own_tri and n[k] == 1
            forced = {}
            def force(p, v):
                nonlocal bad
                p = frozenset(p)
                if v < 0 or v > 2 or (p in forced and forced[p] != v):
                    bad = True
                forced[p] = v
            for tr in tris:
                s = set(tr)
                if i in s and j not in s:
                    force(s - {i}, 1)                                     # C1
                elif j in s and i not in s:
                    rest = s - {j}
                    if sum(n[x] for x in rest) > 1: bad = True            # C2 / C3
                    elif k not in rest: force(rest, 1 - sum(n[x] for x in rest))   # C2 value
                if k in s and i not in s and j not in s:
                    cd = s - {k}; ab = set(R) - cd
                    force(ab, 3 - sum(n[x] for x in ab))                  # CL
            if bad:
                imm.add(key)
            twos = sum(1 for v in forced.values() if v == 2)
            zeros = sum(1 for v in forced.values() if v == 0)
            simple[key] = max(twos, n[k] + zeros)
    return others, exact, simple, imm

def run(j, target):
    others, exact, simple, imm = block_data2(PASCH, j)
    # soundness: immediate exclusions must be truly infeasible, simple bounds <= exact
    for key in exact:
        if key in imm: assert exact[key] is None, ('unsound exclusion', key)
        elif exact[key] is not None: assert simple[key] <= exact[key], ('unsound bound', key)
    val0 = {key: (None if key in imm else simple[key]) for key in exact}
    ex = [(key, v) for key, v in exact.items() if key not in imm and (v is None or v > simple[key])]
    val = dict(val0)
    for key, v in ex: val[key] = v
    kept = list(ex)
    for item in sorted(ex, key=lambda kv: (-(kv[1] if kv[1] is not None else 3), str(kv[0]))):
        key, v = item
        trial = dict(val); trial[key] = val0[key]
        mm = cycle_min(others, trial)
        if mm is None or mm >= target:
            val = trial; kept.remove(item)
    return len(kept), kept, cycle_min(others, val0)

if __name__ == '__main__':
    res = {}
    for name, j, targets in (('w', 3, range(0, 7)), ('p', 1, range(0, 9)), ('a', 0, range(0, 8))):
        for t in targets:
            c, kept, m0 = run(j, t)
            res[(name, t)] = (c, kept)
        print(f'{name} (block {j}): immediate-only minimum = {m0}; extras needed per target: '
              + ', '.join(f'{t}:{res[(name, t)][0]}' for t in targets), flush=True)
    best = None
    for tw in range(0, 7):
        for tp in range(0, 9):
            for ta in range(0, 8):
                if tw + 3 * tp + 3 * ta < 43: continue
                cs = [res[('w', tw)][0], res[('p', tp)][0], res[('a', ta)][0]]
                if best is None or sum(cs) < best[0]: best = (sum(cs), tw, tp, ta)
    print('cheapest:', best)
    for name, t in zip('wpa', best[1:]):
        print(f'  {name} >= {t}:', [(i, sorted(pr), 'EXCL' if v is None else v) for (i, pr), v in res[(name, t)][1]])
