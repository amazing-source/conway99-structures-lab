"""Pure twisted configurations of C2: the D-graph is a decomposition of E(K7) into cycles (length >= 3).
Step 1: count labelled decompositions (memo on the set of uncovered edges) and list the cycle-length types.
Step 2: isomorphism classes, by canonical form under S7 (only for the types we need)."""
import itertools, functools, sys
from collections import Counter

V = range(7)
EDGES = list(itertools.combinations(V, 2))
EID = {e: i for i, e in enumerate(EDGES)}
def eid(a, b): return EID[(a, b) if a < b else (b, a)]

def cycles_through(edge_mask, e0):
    """All cycles (as frozenset of edge ids + vertex tuple) using only edges in edge_mask and containing edge e0."""
    a, b = EDGES[e0]
    out = []
    # simple paths from b back to a avoiding e0, using available edges, length >= 2
    def dfs(path, used):
        u = path[-1]
        for w in V:
            if w in path[1:] or w == u: continue
            e = eid(u, w)
            if not (edge_mask >> e) & 1 or e in used: continue
            if w == a:
                if len(path) >= 3 and e != e0:
                    out.append((frozenset(used | {e, e0}), tuple(path)))
                continue
            dfs(path + [w], used | {e})
    dfs([a, b], set())
    # each cycle found twice? path starts a->b fixed direction, so each cycle containing e0 found once
    return out

FULL = (1 << 21) - 1
@functools.lru_cache(maxsize=None)
def count(mask):
    if mask == 0: return Counter({(): 1})
    e0 = (mask & -mask).bit_length() - 1
    tot = Counter()
    for es, cyc in cycles_through(mask, e0):
        sub = mask
        for e in es: sub &= ~(1 << e)
        for typ, c in count(sub).items():
            tot[tuple(sorted(typ + (len(cyc),)))] += c
    return tot

if __name__ == '__main__':
    res = count(FULL)
    print('labelled decompositions of K7 into cycles, by cycle-length type:')
    for typ, c in sorted(res.items()):
        print(f'   {typ}: {c}')
    print('total', sum(res.values()))
