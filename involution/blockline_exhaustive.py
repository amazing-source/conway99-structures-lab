"""Independent check of C1's exhaustive claim: for EVERY pair (tau, iota) of fixed-point-free involutions of the 12 labels, the
19x19 odd block-line submatrix of S - 3I has rank_7 >= 14, i.e. dim ker Q_c <= 4.
  - iota = mate o rho is fixed-point-free by F0 (a root is never adjacent to its sigma-image), so only fpf iota occur;
  - tau is taken up to the centralizer of mate (signed permutations of the 6 blocks), which acts on (tau, iota) jointly and
    preserves the rank; its classes are the cycle types of the 2-factor mate u tau, i.e. the 11 partitions of 6.
Uses the matrix of blockline_rank.py (copied here so that importing does not run its sampling)."""
from collections import Counter
P = 7
labels = [(e, s) for e in range(6) for s in (0, 1)]
LI = {l: i for i, l in enumerate(labels)}
mate = {l: (l[0], 1 - l[1]) for l in labels}


def rank_mod(M, p=P):
    M = [[x % p for x in r] for r in M]
    r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


def matrix(tau, iota):
    n = 19
    M = [[0] * n for _ in range(n)]
    M[0][0] = -1
    for e in range(6):
        M[1 + e][1 + e] = -1
    for a, v in enumerate(labels):
        o = 7 + a
        M[0][o] = M[o][0] = 1
        e, s = v
        sg = 1 if s == 0 else -1
        M[1 + e][o] = M[o][1 + e] = sg
        M[o][7 + LI[tau[v]]] += 1
        M[o][7 + LI[iota[v]]] -= 1
    return [[M[i][j] - 3 * (i == j) for j in range(n)] for i in range(n)]


def fpf_involutions(pts):
    if not pts:
        yield {}
        return
    a = pts[0]
    for i in range(1, len(pts)):
        b = pts[i]
        rest = pts[1:i] + pts[i + 1:]
        for f in fpf_involutions(rest):
            g = dict(f)
            g[a] = b
            g[b] = a
            yield g


def partitions(n, mx=None):
    mx = n if mx is None else mx
    if n == 0:
        yield []
        return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k):
            yield [k] + rest


def tau_of_type(parts):
    """alternating cycles of mate u tau of lengths 2l: blocks b_1..b_l, tau((b_t,1)) = (b_{t+1},0) cyclically."""
    tau = {}
    b = 0
    for l in parts:
        blocks = list(range(b, b + l))
        b += l
        if l == 1:
            tau[(blocks[0], 0)] = (blocks[0], 1)
            tau[(blocks[0], 1)] = (blocks[0], 0)
            continue
        for t in range(l):
            x, y = (blocks[t], 1), (blocks[(t + 1) % l], 0)
            tau[x] = y
            tau[y] = x
    assert all(tau[tau[v]] == v and tau[v] != v for v in labels)
    return tau


if __name__ == '__main__':
    allinv = list(fpf_involutions(labels))
    assert len(allinv) == 10395
    total = Counter()
    for parts in partitions(6):
        tau = tau_of_type(parts)
        h = Counter(rank_mod(matrix(tau, iota)) for iota in allinv)
        total.update(h)
        print('tau type', parts, 'rank_7 histogram', dict(sorted(h.items())), flush=True)
    print('ALL (11 tau classes x 10395 iota):', dict(sorted(total.items())), '=> min rank', min(total), ', max dim ker Q_c =', 18 - min(total))
