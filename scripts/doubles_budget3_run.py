"""doubles_budget2 + the mu = 2 condition with every element of every block of R_x (C8 below), with unit tests,
a fixed-relabelling cross-check of the canonicalisation, and an independent label-level solver (see labellevel_check).

Roles after relabelling: i = 0, j = 1, k = 2, R_x = {3,4,5,6}.  n_r = [r is a cycle-neighbour of X on the free 6-cycle
at j] (so exactly two of n_2..n_6 are 1).  Unknowns: phi(P) in {0,1,2} for the 6 pairs P in R_x, y_r in {0,1,2} for
r in R_x (y_r = number of labels of pair kr adjacent to x).  Constraints (derivations in NOTE / journal):
  C0  own triangle {0,1,2}:            n_2 = 0
  C1  triangle {0,a,b}, a,b in R_x:    phi(ab) = 1
  C2  triangle {1,a,b}, a,b in R_x:    n_a + n_b + phi(ab) = 1
  C3  triangle {1,2,b}, b in R_x:      n_2 + n_b + y_b = 1
  C4  triangle {2,a,b}, a,b in R_x:    y_a + y_b + phi(ab) = 3
  C5  triangle inside R_x:             phi over its three pairs = 3
  C6  (4 labels meet k):               sum_r y_r = 2 - n_2
  C8  (4 labels meet t, t in R_x):     n_t + y_t + sum_{s in R_x - t} phi(ts) = 4
  (C7: sum phi = 6 + n_2 follows from C6 + C8 and is asserted.)
doubles(x) >= #{P : phi(P) = 2} + #{r : y_r = 2}."""
import functools
import itertools
import sys

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions, CYC6, V

R = (3, 4, 5, 6)
PAIRS = [p for p in itertools.combinations(R, 2)]
PERMS4 = list(itertools.permutations(R))


KEYS = set()


def local_min(key):
    KEYS.add(key)
    return _local_min(key)


@functools.lru_cache(maxsize=None)
def _local_min(key):
    tris, nbrs = key
    n = {r: (1 if r in nbrs else 0) for r in range(7)}
    assert sum(n[r] for r in (2,) + R) == 2 and n[0] == 0 and n[1] == 0
    fixed, eq3, eqy, eqjk = {}, [], [], []
    for t in tris:
        s = set(t)
        if 0 in s and 1 in s:
            assert s == {0, 1, 2}, 'a triangle through edge ij must be the own cycle'
            if n[2]:
                return None
        elif 0 in s:
            assert 2 not in s, 'edge ik lies on the cycle through ij'
            a, b = sorted(s - {0})
            if fixed.get((a, b), 1) != 1:
                return None
            fixed[(a, b)] = 1
        elif 1 in s:
            a, b = sorted(s - {1})
            if 2 in (a, b):
                bb = a if b == 2 else b
                eqjk.append(bb)
                if n[2] + n[bb] > 1:
                    return None
            else:
                v = 1 - n[a] - n[b]
                if v < 0:
                    return None
                if fixed.get((a, b), v) != v:
                    return None
                fixed[(a, b)] = v
        elif 2 in s:
            a, b = sorted(s - {2})
            eqy.append((a, b))
        else:
            eq3.append([p for p in itertools.combinations(sorted(s), 2)])
    ysum = 2 - n[2]
    best = None
    free = [p for p in PAIRS if p not in fixed]
    for y in itertools.product(range(3), repeat=4):
        if sum(y) != ysum:
            continue
        yd = dict(zip(R, y))
        if any(yd[b] != 1 - n[2] - n[b] for b in eqjk):
            continue
        ydoub = sum(1 for v in y if v == 2)
        for vals in itertools.product(range(3), repeat=len(free)):
            phi = dict(fixed)
            phi.update(zip(free, vals))
            if any(sum(phi[p] for p in grp) != 3 for grp in eq3):
                continue
            if any(yd[a] + yd[b] + phi[(a, b)] != 3 for a, b in eqy):
                continue
            ok = True
            for t in R:
                if n[t] + yd[t] + sum(phi[tuple(sorted((t, s)))] for s in R if s != t) != 4:
                    ok = False
                    break
            if not ok:
                continue
            assert sum(phi.values()) == 6 + n[2]
            d = sum(1 for v in phi.values() if v == 2) + ydoub
            if best is None or d < best:
                best = d
    return best


def relabel_key(tris, i, j, k, Rx, nbrs, perm):
    mp = {i: 0, j: 1, k: 2}
    mp.update(dict(zip(Rx, perm)))
    return (tuple(sorted(tuple(sorted(mp[v] for v in t)) for t in tris)), tuple(sorted(mp[r] for r in nbrs)))


def canon(tris, i, j, k, Rx, nbrs):
    return min(relabel_key(tris, i, j, k, Rx, nbrs, p) for p in PERMS4)


def analyse(dec, canonical=True):
    cyc_of_edge = {}
    for c in dec:
        m = len(c)
        for t in range(m):
            cyc_of_edge[frozenset((c[t], c[(t + 1) % m]))] = c
    tris = [tuple(sorted(c)) for c in dec if len(c) == 3]

    def other_nb(i, j_):
        c = cyc_of_edge[frozenset((i, j_))]
        m = len(c)
        t = c.index(i)
        a, b = c[(t - 1) % m], c[(t + 1) % m]
        return b if a == j_ else a

    total, per_block = 0, []
    for j in V:
        others = [r for r in V if r != j]
        k_of = {i: other_nb(i, j) for i in others}
        cache = {}

        def lm(i, pair):
            if (i, pair) not in cache:
                k = k_of[i]
                Rx = [r for r in V if r not in (i, j, k)]
                key = canon(tris, i, j, k, Rx, pair) if canonical else relabel_key(tris, i, j, k, Rx, pair, R)
                cache[(i, pair)] = local_min(key)
            return cache[(i, pair)]

        best = None
        for order in CYC6:
            lab = [others[p] for p in order]
            tot = 0
            for t in range(6):
                i, l, rr = lab[t], lab[t - 1], lab[(t + 1) % 6]
                v = lm(i, frozenset((l, rr)))
                if v is None:
                    tot = None
                    break
                tot += v
            if tot is not None and (best is None or tot < best):
                best = tot
        if best is None:
            return None, per_block
        per_block.append(best)
        total += best
    return total, per_block


def unit_tests():
    assert local_min(((), (3, 4))) == 0
    assert local_min(((), (2, 3))) == 1
    fano = tuple(sorted([(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]))
    assert local_min((fano, (3, 4))) == 2
    assert local_min((fano, (3, 5))) is None          # two R_j-neighbours in the same triangle {1,3,5}
    assert local_min((fano, (2, 3))) is None          # own triangle forbids a neighbour in jk
    print('unit tests passed')


if __name__ == '__main__':
    unit_tests()
    mint = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    res = {}
    for dec in decompositions((1 << 21) - 1, []):
        nt = sum(1 for c in dec if len(c) == 3)
        if nt < mint:
            continue
        typ = tuple(sorted(len(c) for c in dec))
        res.setdefault(typ, []).append(analyse(dec)[0])
    print(f'local keys solved: {len(KEYS)}')
    import json
    json.dump([[list(map(list, k[0])), list(k[1]), local_min(k)] for k in sorted(KEYS)], open(sys.argv[2] if len(sys.argv) > 2 else 'logs/doubles3_keys.json', 'w'))
    for typ, zs in sorted(res.items(), key=lambda kv: (-kv[0].count(3), kv[0])):
        imp = sum(z is None for z in zs)
        vals = [z for z in zs if z is not None]
        killed = all(z is None or z > 42 for z in zs)
        print(f'{typ}: {len(zs)} labelled; impossible outright {imp}; lower bound: min {min(vals) if vals else "-"}, '
              f'max {max(vals) if vals else "-"}; KILLED: {killed}  (surviving labelled: '
              f'{sum(1 for z in zs if z is not None and z <= 42)})', flush=True)
