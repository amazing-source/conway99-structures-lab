"""Per-block realisability test in the ROOTLESS envelope U = O24 (odd Leech lattice) - structures 27/09, route (C).
Sol (0836, ROOTLESS_SHADOWFRAME): in a rootless envelope, U = {a/sqrt8 : a_i = m mod 2, {i : a_i = m+2 mod 4} in C,
sum a = 0 mod 8} with C a doubly-even self-dual [24,12,8] code (= Golay, classical uniqueness); the 42 u are sign codewords
E/sqrt8 (E = +-1, negative set W in C), the 7 h are signed octads 2M/sqrt8 (M in {0,+-1}, weight 8).
One block c, gauge u_v = orbit of (c+, v) for its 12 labels v: then (h_c, u_v) = -C = -1 and, from S = [[-I, C],[C^t, T]]
restricted to star(c), (u_v, u_w) = -S = -[mu v = w] + [tau v = w] (v != w).  With h_c = 2*1_O/sqrt8 (WLOG, sign changes
by codewords):  |W_v & O| = 6 (two minority positions e_v on O),  3 h_c + sum_v u_v = 0  <=>  every point of O lies in
exactly 3 minority edges and every point off O lies in exactly 6 of the W_v;  |W_v ^ W_w| = 12 - 4 (u_v, u_w)*... i.e.
16 if w = mu v != tau v, 8 if w = tau v != mu v, 12 otherwise.  A type with no solution cannot occur at ANY block of a
frame in O24.  Usage: python o24_block_test.py [first last]"""
import json, sys, time
from itertools import combinations

wt = lambda x: bin(x).count('1')
g = sum(1 << b for b in (11, 9, 7, 6, 5, 1, 0))
WORDS = []
for m in range(4096):
    w = 0
    for j in range(12):
        if m >> j & 1:
            w ^= g << j
    if wt(w) % 2:
        w |= 1 << 23
    WORDS.append(w)
assert len(set(WORDS)) == 4096
CODE = set(WORDS)
OCT = [w for w in WORDS if wt(w) == 8]
assert len(OCT) == 759
O = OCT[0]
Opts = [i for i in range(24) if O >> i & 1]
DO = [w for w in WORDS if w & O == 0]                     # codewords disjoint from O: 32
assert len(DO) == 32
OPTS = []                                                 # (word W, edge e) with W & O = O \ e
for a, b in combinations(Opts, 2):
    e = (1 << a) | (1 << b)
    tgt = O & ~e
    base = next(w for w in WORDS if w & O == tgt)
    for d in DO:
        OPTS.append((base ^ d, e))
NOPT = len(OPTS)
assert NOPT == 28 * 32
# compatibility bitsets by required distance
COMP = {d: [0] * NOPT for d in (8, 12, 16)}
for i in range(NOPT):
    wi = OPTS[i][0]
    for j in range(NOPT):
        if i != j:
            dd = wt(wi ^ OPTS[j][0])
            if dd in COMP:
                COMP[dd][i] |= 1 << j
OFF = [i for i in range(24) if not O >> i & 1]


def solve(mu, tau):
    """backtracking: labels 0..11; returns a solution (list of option indices) or None"""
    L = list(range(12))
    need = {}
    for v in L:
        for w in L:
            if v != w:
                need[v, w] = 16 if (w == mu[v] and w != tau[v]) else 8 if (w == tau[v] and w != mu[v]) else 12
    order = [0]                                            # follow mu / tau links for early pruning
    while len(order) < 12:
        best = max((v for v in L if v not in order), key=lambda v: sum(need[v, w] != 12 for w in order))
        order.append(best)
    ALL = (1 << NOPT) - 1

    def rec(k, chosen, dom_cache):
        if k == 12:
            # column constraints: O points in exactly 3 edges, off-O points in exactly 6 words
            for p in Opts:
                if sum(1 for i in chosen.values() if OPTS[i][1] >> p & 1) != 3:
                    return None
            for p in OFF:
                if sum(1 for i in chosen.values() if OPTS[i][0] >> p & 1) != 6:
                    return None
            return dict(chosen)
        v = order[k]
        dom = ALL
        for w, i in chosen.items():
            dom &= COMP[need[v, w]][i]
            if not dom:
                return None
        # partial column pruning: edge multiplicity per O point <= 3
        cnt = {p: sum(1 for i in chosen.values() if OPTS[i][1] >> p & 1) for p in Opts}
        full = [p for p in Opts if cnt[p] >= 3]
        while dom:
            low = dom & -dom
            i = low.bit_length() - 1
            dom ^= low
            if any(OPTS[i][1] >> p & 1 for p in full):
                continue
            chosen[v] = i
            r = rec(k + 1, chosen, None)
            if r is not None:
                return r
            del chosen[v]
        return None
    # WLOG (Stab of the signed octad is transitive on the 896 options): first label -> option 0
    return rec(1, {order[0]: 0}, None)


if __name__ == '__main__':
    T = json.load(open('all_types_3308.json'))
    a = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    b = int(sys.argv[2]) if len(sys.argv) > 2 else len(T)
    t0 = time.time(); alive = []
    for t in T[a:b]:
        r = solve(t['mu'], t['tau'])
        if r is not None:
            alive.append(t['id'])
    print('types %d..%d: %d realisable in O24 at one block, %d impossible (%.0fs)' % (a, b - 1, len(alive), b - a - len(alive),
                                                                                   time.time() - t0), flush=True)
    print('realisable ids:', alive[:50])
