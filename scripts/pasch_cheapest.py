"""PASCH: cheapest hand proof of total doubles >= 43.
Immediate facts (free): exclusions by C0 / C2 / C3 (two neighbours in a triangle through j, or k when the own cycle
is a triangle, or k and the third vertex of the triangle {j,k,b}); doubles >= n_k + (zeros forced by C2 exports).
Extra facts (cost 1 each): non-immediate exclusions, and doubles values above the simple bound.
For each block type (w = 3, p = 1, a = 0; Z3 makes blocks of a type equivalent) and each target t, find (greedy then a
small exact search) a minimum set of extras giving min over 6-cycles >= t.  Then pick t_w + 3 t_p + 3 t_a >= 43
minimising cost_w + cost_p + cost_a."""
import itertools
import sys

sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V, CYC6

PASCH = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 6), (1, 5, 6), (2, 3, 4, 5), (2, 4, 6)]


def block_data(dec, j):
    cyc_of_edge = {}
    for c in dec:
        for t in range(len(c)):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % len(c)]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]
    others = [r for r in V if r != j]
    exact, simple, immediate_excl = {}, {}, set()
    for i in others:
        c = cyc_of_edge[frozenset((i, j))]; m = len(c); t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]
        k = b if a == j else a
        R = [r for r in V if r not in (i, j, k)]
        own_tri = len(c) == 3
        for l, r in itertools.combinations([o for o in others if o != i], 2):
            pr = frozenset((l, r))
            key = (i, pr)
            exact[key] = d3.local_min(d3.canon(tris, i, j, k, R, pr))
            n = {x: (1 if x in pr else 0) for x in V}
            imm = own_tri and n[k] == 1
            zeros = 0
            for tr in tris:
                s = set(tr)
                if j in s and i not in s:
                    rest = s - {j}
                    if sum(n[x] for x in rest) > 1:
                        imm = True                              # C2 or C3 violated
                    elif k not in rest and sum(n[x] for x in rest) == 1:
                        zeros += 1
            if imm:
                immediate_excl.add(key)
            simple[key] = n[k] + zeros
    return others, exact, simple, immediate_excl


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


def base_val(exact, simple, immediate_excl):
    val = {}
    for key in exact:
        val[key] = None if key in immediate_excl else simple[key]
    return val


def extras_list(exact, simple, immediate_excl):
    ex = []
    for key, v in exact.items():
        if key in immediate_excl:
            continue
        if v is None:
            ex.append((key, None))
        elif v > simple[key]:
            ex.append((key, v))
    return ex


def min_cost(j, target, max_size=7):
    others, exact, simple, imm = block_data(PASCH, j)
    val0 = base_val(exact, simple, imm)
    ex = extras_list(exact, simple, imm)
    m0 = cycle_min(others, val0)
    if m0 is not None and m0 >= target:
        return 0, []
    # exact search by increasing size (small)
    for size in range(1, max_size + 1):
        for comb in itertools.combinations(ex, size):
            val = dict(val0)
            for key, v in comb:
                val[key] = v
            m = cycle_min(others, val)
            if m is None or m >= target:
                return size, list(comb)
    return None, None


if __name__ == '__main__':
    res = {}
    for name, j, targets in (('w', 3, range(0, 7)), ('p', 1, range(0, 9)), ('a', 0, range(0, 8))):
        for t in targets:
            c, comb = min_cost(j, t)
            res[(name, t)] = (c, comb)
            print(f'{name} block {j}: target {t}: minimum extras {c}', flush=True)
    best = None
    for tw in range(0, 7):
        for tp in range(0, 9):
            for ta in range(0, 8):
                if tw + 3 * tp + 3 * ta < 43:
                    continue
                cs = [res[('w', tw)][0], res[('p', tp)][0], res[('a', ta)][0]]
                if None in cs:
                    continue
                cost = sum(cs)
                if best is None or cost < best[0]:
                    best = (cost, tw, tp, ta)
    print('cheapest:', best)
    cost, tw, tp, ta = best
    for name, t in (('w', tw), ('p', tp), ('a', ta)):
        print(f'  {name} >= {t}: extras {res[(name, t)][1]}')
