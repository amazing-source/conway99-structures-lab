"""Independent checks of c2exp/scripts/root_colmodel.py (no LP solved).
T1 column semantics on every near configuration x sampled far rows; T2 far_sets exhaustiveness vs an independent DFS;
T3 entries() vs row semantics recomputed from scratch; T4 scores() vs entries() for random duals."""
import sys, time, random
import numpy as np
sys.path.insert(0, r'C:\Users\RESEARCHER\AppData\Local\Temp\claude\C--Users-RESEARCHER-Desktop\c88c221a-efb9-4373-af95-1e500b2e81aa\scratchpad\c2exp\scripts')
import root_colmodel as R

t0 = time.time()
In = {0: [4], 1: [5], 2: [0], 3: [0, 1], 4: [1, 2], 5: [2, 3], 6: [0, 1, 2, 3, 4, 5]}     # M1 n=6 job 0 (note 40)
anch = {(min(c, x), max(c, x)) for c, xs in In.items() for x in xs}
tw = [pr for pr in R.PAIRS if pr not in anch]
print('twisted pairs', tw)
HOL = sys.argv[1] if len(sys.argv) > 1 else '0'
bits = [int(ch) for ch in HOL]
twe = {pr: (bits[t] if t < len(bits) else 0) for t, pr in enumerate(tw)}
try:
    g = R.RootModel(In, twe)
except ValueError as e:
    print('BUILD FAILED:', e)
    raise
print('built in %.0fs: near configs %d, columns %d, rows %d, far keys %d' % (
    time.time() - t0, sum(len(l) for l in g.near), g.ncols(), len(g.rname), len(R._farcache)))


def decode(x):
    a, b = R.PAIRS[x >> 2]
    return {a: -1 if (x >> 1) & 1 else 1, b: -1 if x & 1 else 1}


def orb(x):
    d = decode(x)
    (a, sa), (b, tb) = sorted(d.items())
    return (a, b, sa * tb)


def sigma(x):
    d = decode(x)
    (a, sa), (b, tb) = sorted(d.items())
    return R.root_id(a, b, 0 if sa < 0 else 1, 0 if tb < 0 else 1)


def orbit_index(x):
    return 2 * (x >> 2) + ((x >> 1 & 1) ^ (x & 1))


bad = []
rng = random.Random(1)
ncheck = 0
for Z in range(42):
    z = R.rep(Z)
    dz = decode(z)
    i = g.pivot[Z]
    j = [c for c in dz if c != i][0]
    for near, F in g.near[Z]:
        picks = [0, len(F) - 1] + [rng.randrange(len(F)) for _ in range(3)]
        for t in picks:
            rs = list(near) + [int(v) for v in F[t]]
            ncheck += 1
            S = set(rs)
            if len(S) != 12 or len(rs) != 12:
                bad.append(('size', Z, near)); continue
            if z in S or sigma(z) in S:
                bad.append(('self', Z, near)); continue
            ob = {}
            for x in rs:
                ob.setdefault(orb(x), []).append(x)
            twos = [o for o, l in ob.items() if len(l) == 2]
            if len(twos) != 1 or any(len(l) > 2 for l in ob.values()):
                bad.append(('orbits', Z, near)); continue
            Dorb = twos[0]
            if i not in Dorb[:2]:
                bad.append(('Dsupp', Z, near)); continue
            cnt = {}
            for x in rs:
                for c, s in decode(x).items():
                    cnt[(c, s)] = cnt.get((c, s), 0) + 1
            for c in range(7):
                for s in (1, -1):
                    want = 1 if c in (i, j) else 2
                    if cnt.get((c, s), 0) != want:
                        bad.append(('count', Z, near, c, s, cnt.get((c, s), 0)))
            for x in rs:
                if j in decode(x) and g.pivot[orbit_index(x)] == j:
                    bad.append(('pivot_j', Z, near, x))
            if g.pivot[orbit_index(ob[Dorb][0])] != i:
                bad.append(('Dpivot', Z, near))
print('T1 columns checked %d, problems %d %s' % (ncheck, len(bad), bad[:3]))


def my_far(fb, need):
    prs = [(a, b) for a in fb for b in fb if a < b]
    out = []
    rem = dict(need)
    chosen = []

    def opts(a, b):
        roots = [R.root_id(a, b, sa, tb) for sa in (0, 1) for tb in (0, 1)]
        res = [()] + [(x,) for x in roots]
        for x in roots:
            for y in roots:
                if x < y and orb(x) != orb(y):
                    res.append((x, y))
        return res

    def rec(p):
        if p == len(prs):
            if all(v == 0 for v in rem.values()):
                out.append(tuple(sorted(chosen)))
            return
        a, b = prs[p]
        for o in opts(a, b):
            ok = True
            for x in o:
                for c, s in decode(x).items():
                    k = (c, 0 if s > 0 else 1)
                    rem[k] -= 1
                    if rem[k] < 0:
                        ok = False
            if ok:
                chosen.extend(o)
                rec(p + 1)
                del chosen[len(chosen) - len(o):]
            for x in o:
                for c, s in decode(x).items():
                    rem[(c, 0 if s > 0 else 1)] += 1
    rec(0)
    return set(out)


keys = list(R._farcache.keys())
rng.shuffle(keys)
t2 = time.time()
mism = 0
nk = 0
for key in keys:
    if time.time() - t2 > 120:
        break
    fb, needt = key
    mine = my_far(fb, dict(needt))
    theirs = set(tuple(sorted(int(v) for v in row)) for row in R._farcache[key])
    nk += 1
    if mine != theirs:
        mism += 1
        print('  far mismatch', key, len(mine), len(theirs))
print('T2 far keys compared %d / %d, mismatches %d (%.0fs)' % (nk, len(keys), mism, time.time() - t2))


def my_entries(Z, rs):
    rs = set(rs)
    d = {}

    def add(name, v):
        r = g.rows[name]
        d[r] = d.get(r, 0) + v
    add(('one', Z), 1)
    for U in range(42):
        for W in range(U + 1, 42):
            if Z == U:
                add(('symP', U, W), int(R.rep(W) in rs)); add(('symM', U, W), int(sigma(R.rep(W)) in rs))
            if Z == W:
                add(('symP', U, W), -int(R.rep(U) in rs)); add(('symM', U, W), -int(sigma(R.rep(U)) in rs))
    for U in range(42):
        for W in range(U, 42):
            for sg in (1, -1):
                u = R.rep(U)
                v = R.rep(W) if sg > 0 else sigma(R.rep(W))
                ut = R.rep(U) if sg > 0 else sigma(R.rep(U))
                val = 2 * (int(u in rs and v in rs) + int(sigma(u) in rs and sigma(v) in rs))
                if Z == U:
                    val += int(v in rs)
                if Z == W:
                    val += int(ut in rs)
                if val:
                    add(('I', U, W, sg), val)
    for (c, S) in g.cuts:
        if Z in S:
            val = sum(int(R.rep(F) in rs) + int(sigma(R.rep(F)) in rs) for F in S if F != Z)
            if val:
                add(('cut', c, S), val)
    return {r: v for r, v in d.items() if v}


t3 = time.time()
diff = 0
n3 = 0
for Z in range(42):
    for near, F in rng.sample(g.near[Z], min(3, len(g.near[Z]))):
        rs = list(near) + [int(v) for v in F[rng.randrange(len(F))]]
        a = dict(g.entries(Z, frozenset(rs)))
        b = my_entries(Z, rs)
        n3 += 1
        if a != b:
            diff += 1
            print('  entries mismatch Z', Z, sorted(set(a.items()) ^ set(b.items()))[:6])
print('T3 entries compared %d, mismatches %d (%.0fs)' % (n3, diff, time.time() - t3))

y = np.array([rng.randint(-5, 5) for _ in range(len(g.rname))], dtype=np.int64)
Y = g.Ymat(y)
ff = g.ffpart(Y)
bad4 = 0
n4 = 0
for Z in rng.sample(range(42), 8):
    L = g.lin(y, Z)
    for (near, F, sc) in g.scores(y, Z, Y, L, ff)[:4]:
        for t in range(min(3, len(F))):
            rs = frozenset(near + tuple(int(v) for v in F[t]))
            val = sum(c * int(y[r]) for r, c in g.entries(Z, rs))
            n4 += 1
            if val != int(sc[t]):
                bad4 += 1
print('T4 score checks %d, mismatches %d; total %.0fs' % (n4, bad4, time.time() - t0))
