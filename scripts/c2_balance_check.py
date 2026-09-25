"""D11 check of C2's balance lemma (C2 lab note 23 §8) in its GENERAL form, on the real BvLS graph (lambda=1, mu=2, k=22).
For a line l: N = points adjacent to exactly one point of l, F = points adjacent to none, B = A[N,F], G = A[N,N] (cubic).
Derivation: BB^T = (k-2)I - A_NN + 2J - D - A_NN^2 (D = same-part matrix);
for C a union of components of G: D 1_C = (|C|/3) 1_N (each component meets the 3 parts equally),
so BB^T 1_C = (k-14) 1_C + (5/3)|C| 1_N, sum_y n_C(y) = (k-4)|C|, sum_y n_C(y)^2 = (k-14)|C| + (5/3)|C|^2.
At k = 14 (Conway 99) this forces n_C(y) = |C|/6 for every y in F.  On BvLS we check every identity for every line and
every component, plus a negative control (a set that is NOT a union of components)."""
import sys
import numpy as np
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
import graphs

A = np.array(graphs.bvls(), dtype=np.int64); v = len(A); k = int(A[0].sum())
lines = {tuple(sorted((x, y, z))) for x in range(v) for y in np.nonzero(A[x])[0] for z in np.nonzero(A[x] * A[y])[0]}
ok_bbt = ok_part = ok_sums = 0; ncomp = 0; neg_break = neg_total = 0
rng = np.random.default_rng(1)
for L in lines:
    L = list(L)
    cnt = A[:, L].sum(axis=1)
    N = [y for y in range(v) if y not in L and cnt[y] == 1]
    F = [y for y in range(v) if y not in L and cnt[y] == 0]
    part = np.array([next(a for a in L if A[y, a]) for y in N])
    G = A[np.ix_(N, N)]; B = A[np.ix_(N, F)]
    D = (part[:, None] == part[None, :]).astype(np.int64)
    J = np.ones_like(D)
    ok_bbt += np.array_equal(B @ B.T, (k - 2) * np.eye(len(N), dtype=np.int64) - G + 2 * J - D - G @ G)
    # components of G
    seen = -np.ones(len(N), dtype=int); comps = []
    for i in range(len(N)):
        if seen[i] >= 0: continue
        stack = [i]; seen[i] = len(comps); cur = []
        while stack:
            u = stack.pop(); cur.append(u)
            for w in np.nonzero(G[u])[0]:
                if seen[w] < 0: seen[w] = len(comps); stack.append(w)
        comps.append(cur)
    for Cc in comps:
        ncomp += 1
        oneC = np.zeros(len(N), dtype=np.int64); oneC[Cc] = 1; c = len(Cc)
        ok_part += all((oneC[part == a]).sum() * 3 == c for a in L)
        n = B.T @ oneC
        ok_sums += (n.sum() == (k - 4) * c) and (3 * (n ** 2).sum() == 3 * (k - 14) * c + 5 * c * c)
    # negative control: random subset of N of size 12 that is not a union of components
    S = rng.choice(len(N), 12, replace=False); oneS = np.zeros(len(N), dtype=np.int64); oneS[S] = 1
    if not any(set(S) == set(Cc) for Cc in comps):
        neg_total += 1; n = B.T @ oneS; c = 12
        neg_break += not (3 * (n ** 2).sum() == 3 * (k - 14) * c + 5 * c * c)
print(f'BvLS: {len(lines)} lines, |N| = {len(N)}, |F| = {len(F)}')
print(f'  BB^T identity: {ok_bbt}/{len(lines)} lines')
print(f'  components meet the 3 parts equally: {ok_part}/{ncomp} components')
print(f'  sum n_C = (k-4)|C| and sum n_C^2 = (k-14)|C| + (5/3)|C|^2: {ok_sums}/{ncomp} components')
print(f'  negative control (random 12-sets, not unions of components): second-moment identity breaks on {neg_break}/{neg_total}')
