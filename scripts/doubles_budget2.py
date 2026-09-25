"""Exact local system for doubles in pure twisted configurations of C2 (all 21 pairs twisted, Sigma = empty, so n_s = 0).

Orbit X in pair ij, paired via i to D(X) in pair ik (k = other neighbour of i on the D-cycle through edge ij), free at j.
Label x = {a, b}, a in block i, b in block j.  R_x = [7] - {i, j, k} (4 blocks).  Facts (all from C2's lemmas):
  (F1) the only labels adjacent to x that meet block i are v, sigma v (in pair ik): phi_x(ir) = 0 for r != k;
  (F2) the labels adjacent to x that meet block j are exactly its two R_j-neighbours, i.e. its two neighbours on the free
       6-cycle at j (one pair jr each): phi_x(jr) = n_r := [r is a cycle-neighbour of X at j];
  (F3) exactly 4 labels adjacent to x meet block k (mu = 2 with (k,0) and with (k,1)): v, sigma v, the R_j-neighbour in
       pair jk if any (n_k), and 2 - n_k labels in pairs kr, r in R_x:  sum_r y_r = 2 - n_k,  y_r = phi_x(kr);
  (F4) hence the 12 labelled neighbours of x give  sum_{P in R_x} phi_x(P) = 6 + n_k   (omega(x) = n_k, C2 18.3);
  (F5) exact counts of complete D-triangles tau (C2 §25):  |N(x) ∩ Lab(tau)| = 2 / 1 / 3  if x's pair is in tau / meets
       tau in one block / is disjoint from tau.  Cases:
         tau = {i,j,k} (own):          n_k = 0;
         tau = {i,a,b}, a,b in R_x:    phi(ab) = 1           (phi(ia) = phi(ib) = 0 by F1; k in tau is impossible);
         tau = {j,a,b}, a,b in R_x:    n_a + n_b + phi(ab) = 1;
         tau = {j,k,b}, b in R_x:      n_k + n_b + y_b = 1;
         tau = {k,a,b}, a,b in R_x:    y_a + y_b + phi(ab) = 3;
         tau inside R_x:               phi over its 3 pairs = 3.
  doubles(x) >= #{P in R_x : phi(P) = 2} + #{r : y_r = 2}   (pairs jr carry <= 1, pairs ir carry 0 or D(X)).
Budget (C2 16.1 with n_s = 0): sum over the 84 labels of doubles = 84 exactly, i.e. sum over the 42 orbits = 42.
We minimise the local count for each orbit, then over the cyclic order of the free 6-cycle at each block (blocks are
independent; cross-constraints between different orbits are ignored, so the result is a valid LOWER bound).
If the bound exceeds 42 for every labelled decomposition of a type, that type is impossible.
If some block has no admissible cyclic order at all, the decomposition is impossible outright."""
import functools
import itertools
import sys

sys.path.insert(0, 'scripts')
from zeros_budget import decompositions, CYC6, V

PERMS4 = list(itertools.permutations(range(3, 7)))


@functools.lru_cache(maxsize=None)
def local_min(key):
    """key = (triangles as sorted tuples over {0..6}, nbrs as sorted tuple); roles: i=0, j=1, k=2, R_x={3,4,5,6}."""
    tris, nbrs = key
    n = {r: (1 if r in nbrs else 0) for r in range(7)}
    R = (3, 4, 5, 6)
    pairs = [tuple(p) for p in itertools.combinations(R, 2)]
    fixed = {}
    eq_phi3 = []      # lists of phi pairs summing to 3
    eq_y = []         # (a, b, pair_ab): y_a + y_b + phi(ab) = 3
    eq_jk = []        # b: n_k + n_b + y_b = 1  ->  y_b = 1 - n_k - n_b
    for t in tris:
        s = set(t)
        if s == {0, 1, 2}:
            if n[2]:
                return None
        elif 0 in s and 1 not in s:
            a, b = sorted(s - {0})
            if 2 in (a, b):
                return None                 # cannot happen (edge ik lies on the cycle through ij)
            if fixed.get((a, b), 1) != 1:
                return None
            fixed[(a, b)] = 1
        elif 1 in s and 0 not in s:
            a, b = sorted(s - {1})
            if 2 in (a, b):
                bb = b if a == 2 else a
                if n[2] + n[bb] > 1:
                    return None
                eq_jk.append(bb)
            else:
                v = 1 - n[a] - n[b]
                if v < 0:
                    return None
                if fixed.get((a, b), v) != v:
                    return None
                fixed[(a, b)] = v
        elif 0 not in s and 1 not in s:
            if 2 in s:
                a, b = sorted(s - {2})
                eq_y.append((a, b))
            else:
                eq_phi3.append([tuple(p) for p in itertools.combinations(sorted(s), 2)])
    ysum = 2 - n[2]
    target = 6 + n[2]
    best = None
    free = [p for p in pairs if p not in fixed]
    # enumerate y (compositions of ysum into 4 parts, each <= 2)
    for y in itertools.product(range(3), repeat=4):
        if sum(y) != ysum:
            continue
        yd = dict(zip(R, y))
        if any(yd[b] != 1 - n[2] - n[b] for b in eq_jk):
            continue
        ydoub = sum(1 for v in y if v == 2)
        for vals in itertools.product(range(3), repeat=len(free)):
            phi = dict(fixed)
            phi.update(zip(free, vals))
            if sum(phi.values()) != target:
                continue
            if any(sum(phi[p] for p in grp) != 3 for grp in eq_phi3):
                continue
            if any(yd[a] + yd[b] + phi[(a, b)] != 3 for a, b in eq_y):
                continue
            d = sum(1 for v in phi.values() if v == 2) + ydoub
            if best is None or d < best:
                best = d
    return best


def canon(tris, i, j, k, R, nbrs):
    best = None
    for perm in PERMS4:
        mp = {i: 0, j: 1, k: 2}
        mp.update(dict(zip(R, perm)))
        tt = tuple(sorted(tuple(sorted(mp[v] for v in t)) for t in tris))
        nb = tuple(sorted(mp[r] for r in nbrs))
        key = (tt, nb)
        if best is None or key < best:
            best = key
    return best


def analyse(dec):
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

    total = 0
    per_block = []
    for j in V:
        others = [r for r in V if r != j]
        k_of = {i: other_nb(i, j) for i in others}
        cache = {}

        def lm(i, pair):
            key = (i, pair)
            if key not in cache:
                k = k_of[i]
                R = [r for r in V if r not in (i, j, k)]
                cache[key] = local_min(canon(tris, i, j, k, R, pair))
            return cache[key]

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


if __name__ == '__main__':
    mint = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    res = {}
    for dec in decompositions((1 << 21) - 1, []):
        nt = sum(1 for c in dec if len(c) == 3)
        if nt < mint:
            continue
        typ = tuple(sorted(len(c) for c in dec))
        res.setdefault(typ, []).append(analyse(dec)[0])
    print(f'local keys solved: {local_min.cache_info().currsize}')
    for typ, zs in sorted(res.items(), key=lambda kv: (-kv[0].count(3), kv[0])):
        imp = sum(z is None for z in zs)
        vals = [z for z in zs if z is not None]
        killed = all(z is None or z > 42 for z in zs)
        print(f'{typ}: {len(zs)} labelled; impossible outright {imp}; lower bound on sum of doubles: '
              f'min {min(vals) if vals else "-"}, max {max(vals) if vals else "-"}; '
              f'KILLED: {killed}  (surviving labelled: {sum(1 for z in zs if z is not None and z <= 42)})', flush=True)
