"""Quick check on the Fano-tight witness (its sigma satisfies KBLDP): is T(w) = 1 for every label w of X, for every
triad T_M and every ordered pair (Z -> X) of its blocks?  T(w) = number of far vertices of Z whose sigma-special
landing in X touches the label w.  (If it fails here, T = 1 is not a consequence of KBLDP.)"""
import io
import contextlib
import itertools
import sys
from collections import Counter
sys.path.insert(0, 'scripts')
with contextlib.redirect_stdout(io.StringIO()):
    import fano_cellcross_rank as F

dist = Counter()
bad = 0
for M in range(7):
    triad = [R for R in F.BLOCKS if F.third(*sorted(R)) == M]
    assert len(triad) == 3
    for Zb, Xb in itertools.permutations(triad, 2):
        T = Counter()
        for z in F.PAIRS:
            if F.blk(z) != Zb:
                continue
            for via in Zb:
                (e,) = [l for l in z if l[0] == via]
                (f,) = [l for l in z if l[0] != via]
                w = F.sigma[via][f] if e[1] == 0 else F.sinv[via][f]
                if w[0] in Xb:
                    T[w] += 1
        vals = tuple(sorted(T[(c, s)] for c in Xb for s in (0, 1)))
        dist[vals] += 1
        bad += vals != (1, 1, 1, 1)
print(f'witness: 7 triads x 6 ordered pairs = {sum(dist.values())}; T-profiles on the 4 labels of X: {dict(dist)}; '
      f'pairs with T not identically 1: {bad}')
