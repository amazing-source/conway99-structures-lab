"""PASCH class (triangle of triangles): find, per block type, a SMALL set of hand-provable facts that suffices for a
target lower bound on the doubles of the block.  Facts are 'pair P excluded for orbit i' or 'pair P gives >= m doubles
for orbit i'.  Start from the full exact table (all facts), then greedily drop facts while the minimum over 6-cycles
stays >= target.  A dropped exclusion makes the pair allowed with 0 doubles; a dropped doubles fact sets it to 0.
Every retained fact is then re-checked against the exact local system (and will be derived by hand)."""
import itertools
import sys

sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V, CYC6

PASCH = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 6), (1, 5, 6), (2, 3, 4, 5), (2, 4, 6)]


def table(dec, j):
    cyc_of_edge = {}
    for c in dec:
        for t in range(len(c)):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % len(c)]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]

    def other_nb(i, j_):
        c = cyc_of_edge[frozenset((i, j_))]; m = len(c); t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]
        return b if a == j_ else a
    others = [r for r in V if r != j]
    tab, kk = {}, {}
    for i in others:
        k = other_nb(i, j)
        kk[i] = k
        Rx = [r for r in V if r not in (i, j, k)]
        for l, r in itertools.combinations([o for o in others if o != i], 2):
            tab[(i, frozenset((l, r)))] = d3.local_min(d3.canon(tris, i, j, k, Rx, frozenset((l, r))))
    return others, tab, kk


def cycle_min(others, val):
    best = None
    for order in CYC6:
        lab = [others[p] for p in order]
        tot = 0
        for t in range(6):
            v = val[(lab[t], frozenset((lab[t - 1], lab[(t + 1) % 6])))]
            if v is None:
                tot = None; break
            tot += v
        if tot is not None and (best is None or tot < best):
            best = tot
    return best


def reduce(j, target):
    others, tab, kk = table(PASCH, j)
    val = dict(tab)
    assert (cycle_min(others, val) or 99) >= target, 'target above the exact minimum'
    # try to drop high-value doubles facts first, then exclusions
    keys = sorted(tab, key=lambda key: (-(tab[key] if tab[key] is not None else 5), str(key)))
    for key in keys:
        old = val[key]
        if old == 0:
            continue
        val[key] = 0 if old is not None else 0          # weaken to 'allowed with 0 doubles'
        m = cycle_min(others, val)
        if m is None or m >= target:
            continue                                     # fact not needed
        val[key] = old                                   # needed: restore
    facts = {key: v for key, v in val.items() if v != 0}
    print(f'=== block {j} (k: {kk}), target {target}: {len(facts)} facts suffice; check min = {cycle_min(others, val)}')
    for (i, pr), v in sorted(facts.items(), key=lambda kv: (kv[0][0], sorted(kv[0][1]))):
        print(f'   orbit {i}: pair {sorted(pr)} -> {"EXCLUDED" if v is None else f">= {v} doubles"}')
    return facts


if __name__ == '__main__':
    reduce(1, 8)     # intersection block (type p)
    reduce(0, 7)     # private block (type a)
    reduce(3, 1)     # w block, just in case
