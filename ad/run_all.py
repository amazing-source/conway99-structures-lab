"""Pair-level test of every AD orientation class (C2's enumerator, degree + tail lemmas) over all D-structure choices:
A-blocks (d=6) and d=4 blocks: every 2-regular structure on In(k); d=2 digon, d=3 triangle; dd pairs: every 2-regular
multigraph (cycles >= 2) on the dd pairs using edges between disjoint pairs.  Prints survivors."""
import sys, json, itertools, time
from pair_level import *

def dd_structures(dd):
    dd = [tuple(p) for p in dd]
    out = []
    def disj(p, q): return not (set(p) & set(q))
    def rec(rem, cyc):
        if not rem: out.append(list(cyc)); return
        f = rem[0]; others = rem[1:]
        for L in range(1, len(others) + 1):
            for comb in itertools.combinations(others, L):
                for perm in itertools.permutations(comb):
                    if L >= 2 and perm[0] > perm[-1]: continue
                    c = (f,) + perm
                    if all(disj(c[i], c[(i + 1) % len(c)]) for i in range(len(c))) and (len(c) > 2 or disj(c[0], c[1])):
                        cyc.append(c); rec([x for x in others if x not in perm], cyc); cyc.pop()
    rec(dd, [])
    return out

def dd_matching(cycles):
    M = {}
    for c in cycles:
        L = len(c)
        if L == 2:
            a, b = c
            M[(a, 0)] = (b, 0); M[(b, 0)] = (a, 0); M[(a, 1)] = (b, 1); M[(b, 1)] = (a, 1)
        else:
            for i in range(L):
                a, b = c[i], c[(i + 1) % L]
                M[(a, 1)] = (b, 0); M[(b, 0)] = (a, 1)
    return M

def test_rep(rep, stop_first=True):
    dd = [tuple(p) for p in rep['dd']]
    head = {}
    for t, h in rep['arcs']:
        head[pr(t, h)] = h
    In = {c: [x for x in range(N) if x != c and head.get(pr(c, x)) == c] for c in range(N)}
    base = {}; free = []
    for k in range(N):
        d = len(In[k])
        if d in (2, 3): base[k] = [tuple(In[k])]
        elif d in (4, 5, 6): free.append(k)
        elif d == 1: return []
    opts = [two_regular_graphs(In[k]) for k in free]
    ddopts = dd_structures(dd)
    alive = []
    for ddc in ddopts:
        ddm = dd_matching(ddc)
        for choice in itertools.product(*opts):
            Ds = dict(base)
            for k, c in zip(free, choice): Ds[k] = c
            r = feasible(In, dd, ddm, Ds)
            if r is not None:
                alive.append({'ddc': ddc, 'Ds': {k: Ds[k] for k in free}, 'status': r if isinstance(r, str) else 'feasible'})
                if stop_first: return alive
    return alive

if __name__ == '__main__':
    reps = json.load(open(sys.argv[1]))
    lo_, hi_ = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, len(reps))
    t0 = time.time(); bynd = {}
    for t in range(lo_, min(hi_, len(reps))):
        rep = reps[t]; nd = len(rep['dd'])
        al = test_rep(rep)
        s = bynd.setdefault(nd, [0, 0]); s[0] += 1
        if al:
            s[1] += 1
            print('ALIVE rep', t, 'n_dd', nd, rep, al[0]); sys.stdout.flush()
        if t % 50 == 0: print('..', t, '%.0fs' % (time.time() - t0), bynd); sys.stdout.flush()
    print('done', bynd, '%.0fs' % (time.time() - t0))
