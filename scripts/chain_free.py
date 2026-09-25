"""Which decompositions with >= 3 triangles contain NO chain T1 - T - T3 (T1, T3 through distinct vertices of T and
disjoint outside T), i.e. escape the hand lemma S8.2?  Report their types and triangle configurations, and whether the
exhaustive run killed them outright or via the budget."""
import itertools, sys
from collections import Counter
sys.path.insert(0, 'scripts')
from zeros_budget import decompositions
import doubles_budget3 as d3

def has_chain(tris):
    for T in tris:
        for T1, T3 in itertools.permutations([t for t in tris if t != T], 2):
            p, q = T1 & T, T3 & T
            if len(p) == 1 and len(q) == 1 and p != q and not (T1 & T3):
                return True
    return False

def config(tris):
    # pairwise intersection pattern
    inter = [len(a & b) for a, b in itertools.combinations(tris, 2)]
    common = set.intersection(*[set(t) for t in tris]) if tris else set()
    return (len(tris), tuple(sorted(Counter(inter).items())), len(common))

res = Counter()
killed = Counter()
for dec in decompositions((1 << 21) - 1, []):
    tris = [frozenset(c) for c in dec if len(c) == 3]
    if len(tris) < 3:
        continue
    typ = tuple(sorted(len(c) for c in dec))
    if has_chain(tris):
        res[(typ, 'chain')] += 1
        continue
    val = d3.analyse(dec)[0]
    res[(typ, 'chain-free', config(tris), 'outright' if val is None else f'budget {val}')] += 1
for k, v in sorted(res.items(), key=str):
    print(k, v)
