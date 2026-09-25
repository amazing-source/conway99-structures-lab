"""Check (exact, on real graphs) of the form of C1's (P2) used in C2's frame: for two DISJOINT class pairs R, R'
(no common class) at a vertex x of an srg(v,k,1,2) whose blocks have e(R) internal edges,
    sum_z cnt(z,R) cnt(z,R') = 32 - N(R,R')          (z over the far vertices; N counts ordered adjacent pairs),
and for two class pairs sharing ONE class:
    sum_z cnt(z,R) cnt(z,R') = 24 - N(R,R').
These are the off-diagonal (P2) identities, second order.  They are not implied by C2's 16.1 (the diagonal)."""
import itertools
import sys
import numpy as np
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
import graphs

for name, A in (('srg(9,4,1,2)', np.array(graphs.srg9(), dtype=np.int64)),
                ('BvLS(243,22,1,2)', np.array(graphs.bvls(), dtype=np.int64))):
    v = A.shape[0]
    for x in (0, 1, 7):
        Nx = [u for u in range(v) if A[x, u]]
        far = [u for u in range(v) if u != x and not A[x, u]]
        lab = {w: frozenset(u for u in Nx if A[w, u]) for w in far}
        classes = []
        for u in Nx:
            if not any(u in c for c in classes):
                (up,) = [w for w in Nx if w != u and A[u, w]]
                classes.append(frozenset((u, up)))
        blocks = {}
        for c1, c2 in itertools.combinations(range(len(classes)), 2):
            blocks[(c1, c2)] = [w for w in far if len(lab[w] & classes[c1]) == 1 and len(lab[w] & classes[c2]) == 1]
        cnt = {B: A[np.ix_(far, R)].sum(axis=1) for B, R in blocks.items()}
        checked = {0: 0, 1: 0}
        for B1, B2 in itertools.combinations(blocks, 2):
            common = len(set(B1) & set(B2))
            if common == 2:
                continue
            lhs = int(cnt[B1] @ cnt[B2])
            NRR = int(A[np.ix_(blocks[B1], blocks[B2])].sum())
            rhs = (32 if common == 0 else 24) - NRR
            assert lhs == rhs, (name, x, B1, B2, lhs, rhs)
            checked[common] += 1
    print(f'{name}: off-diagonal (P2) exact at 3 vertices: disjoint blocks -> 32 - N, one common class -> 24 - N '
          f'({checked[0]} + {checked[1]} block pairs at the last vertex)')
