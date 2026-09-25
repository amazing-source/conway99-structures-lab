"""Print the local system (C0-C8, in the ORIGINAL block names) for given (block j, orbit i, neighbour pair), plus the
'simple bound' n_k + exported zeros, and the exact minimum (for cross-checking hand derivations)."""
import sys, itertools
sys.path.insert(0, 'scripts')
import doubles_budget3 as d3
from zeros_budget import V

def show(dec, j, i, pair):
    cyc_of_edge = {}
    for c in dec:
        for t in range(len(c)):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % len(c)]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]
    c = cyc_of_edge[frozenset((i, j))]; m = len(c); t = c.index(i)
    a, b = c[(t - 1) % m], c[(t + 1) % m]
    k = b if a == j else a
    R = [r for r in V if r not in (i, j, k)]
    n = {r: (1 if r in pair else 0) for r in V}
    lines = [f'block j={j}, orbit i={i}, k={k}, own cycle {c}, R={R}, neighbours {sorted(pair)}']
    own_tri = len(c) == 3
    if own_tri:
        lines.append(f'  C0: n_{k} = 0' + ('  VIOLATED' if n[k] else ''))
    for tr in tris:
        s = set(tr)
        if {i, j} <= s:
            continue
        if i in s:
            x, y = sorted(s - {i}); lines.append(f'  C1 {tr}: phi{x}{y} = 1')
        elif j in s and k in s:
            y = next(iter(s - {j, k})); lines.append(f'  C3 {tr}: n_{k} + n_{y} + y{y} = 1  ->  y{y} = {1 - n[k] - n[y]}')
        elif j in s:
            x, y = sorted(s - {j}); lines.append(f'  C2 {tr}: n_{x} + n_{y} + phi{x}{y} = 1  ->  phi{x}{y} = {1 - n[x] - n[y]}')
        elif k in s:
            x, y = sorted(s - {k}); lines.append(f'  C4 {tr}: y{x} + y{y} + phi{x}{y} = 3')
        else:
            x, y, z = sorted(s); lines.append(f'  C5 {tr}: phi{x}{y} + phi{x}{z} + phi{y}{z} = 3')
    lines.append(f'  C6: sum y = {2 - n[k]}')
    for tt in R:
        others = [s for s in R if s != tt]
        lines.append(f'  C8({tt}): {n[tt]} + y{tt} + ' + ' + '.join(f'phi{min(tt,s)}{max(tt,s)}' for s in others) + ' = 4')
    # simple bound
    zeros = set()
    for tr in tris:
        s = set(tr)
        if j in s and i not in s:
            rest = s - {j}
            if k not in rest and any(n[r] for r in rest):
                zeros.add(tuple(sorted(rest)))
    exact = d3.local_min(d3.canon(tris, i, j, k, R, frozenset(pair)))
    lines.append(f'  simple bound n_k + exported zeros = {n[k]} + {len(zeros)} {sorted(zeros)} ; exact minimum = {exact}')
    print('\n'.join(lines))

if __name__ == '__main__':
    PASCH = [(0, 1, 2), (0, 3, 1, 4), (0, 5, 3, 6), (1, 5, 6), (2, 3, 4, 5), (2, 4, 6)]
    reqs = [(1, 3, (0, 5)), (1, 5, (2, 4)), (1, 5, (0, 4)), (1, 0, (4, 5)), (1, 3, (2, 5))]
    for j, i, pr in reqs:
        show(PASCH, j, i, pr); print()
