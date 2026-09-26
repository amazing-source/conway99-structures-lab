"""Simulated annealing for the local system (i)-(iii) at a block c (any type (mu, tau), C2 notation), optionally with the
full-block extras (iv) unique neighbour and (v) D-cycle balance.  State: a bijection phi: far roots (60) -> non-mu pairs (60);
move = swap the images of two far roots; cost = |K2 - K2*|_1 + |K3 - K3*|_1 (+ penalties).  cost 0 = a local solution.
Usage: python sa_phi.py cfg.json seconds seed   (cfg: {"mu": [...], "tau": [...]} or {"cycles": [...]} for a full block)"""
import itertools, json, math, random, sys, time

L = list(range(12))
PI = [v ^ 1 for v in L]


def K2K3(mu, tau):
    K2 = [[2] * 12 for _ in L]
    K3 = [[2] * 12 for _ in L]
    for a in L:
        K2[a][a] -= 1; K2[a][PI[a]] -= 1; K2[a][mu[a]] -= 1; K2[a][tau[PI[a]]] -= 1
    for b in L:
        K3[b][PI[b]] -= 1; K3[b][tau[b]] -= 1; K3[b][tau[mu[b]]] -= 1; K3[b][mu[tau[b]]] -= 1
    return K2, K3


def setup(cfg):
    if 'cycles' in cfg:
        from fullblock_rigid import mates_from_D
        m = mates_from_D([tuple(c) for c in cfg['cycles']])
        mu, tau = m, m
        cycles = [tuple(c) for c in cfg['cycles']]
    else:
        mu, tau, cycles = cfg['mu'], cfg['tau'], []
    far = [z for z in itertools.combinations(L, 2) if z[0] // 2 != z[1] // 2]
    fidx = {z: i for i, z in enumerate(far)}
    sig = [fidx[tuple(sorted((PI[z[0]], PI[z[1]])))] for z in far]
    edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    assert len(far) == 60 and len(edges) == 60
    K2, K3 = K2K3(mu, tau)
    return mu, tau, cycles, far, sig, edges, K2, K3


def run(cfg, seconds, seed, extras=True, w3=1):
    rng = random.Random(seed)
    mu, tau, cycles, far, sig, edges, K2t, K3t = setup(cfg)
    n = 60
    labC = [set(2 * b + s for b in C for s in (0, 1)) for C in cycles]

    def pen_pair(z, e, e2):
        """penalty of the sigma-pair (z -> e, sigma z -> e2): (iv) and (v) and no cancellation (full block only)."""
        if not extras or not cycles:
            return 0
        p = 0
        zz = far[z]; se, se2 = set(edges[e]), set(edges[e2])
        p += len(se & se2)
        for a in zz:
            p += (a in se and PI[a] in se2) + (PI[a] in se and a in se2)
        for C, lc in zip(cycles, labC):
            p += abs(len({a // 2 for a in zz} & set(C)) + len(se & lc) + len(se2 & lc) - len(C))
        return p
    phi = list(range(n))
    rng.shuffle(phi)
    K2 = [[0] * 12 for _ in L]
    K3 = [[0] * 12 for _ in L]

    def contrib(z, sign):
        for a in far[z]:
            for b in edges[phi[z]]:
                K2[a][b] += sign
        for b in edges[phi[z]]:
            for d in edges[phi[sig[z]]]:
                K3[b][d] += sign
        for b in edges[phi[sig[z]]]:
            for d in edges[phi[z]]:
                K3[b][d] += sign
    # K3 counts ordered (z, sigma z) for all 60 z: sum_z [b in Phi z][d in Phi sigma z]; contrib adds both orders for the pair,
    # so we add it once per sigma-pair
    reps = [z for z in range(n) if z < sig[z]]
    for z in range(n):
        for a in far[z]:
            for b in edges[phi[z]]:
                K2[a][b] += 1
    for z in reps:
        for b in edges[phi[z]]:
            for d in edges[phi[sig[z]]]:
                K3[b][d] += 1
        for b in edges[phi[sig[z]]]:
            for d in edges[phi[z]]:
                K3[b][d] += 1

    def cost_all():
        c = sum(abs(K2[a][b] - K2t[a][b]) for a in L for b in L) + w3 * sum(abs(K3[b][d] - K3t[b][d]) for b in L for d in L)
        return c + sum(pen_pair(z, phi[z], phi[sig[z]]) for z in reps)
    cur = cost_all()
    best = cur
    t0 = time.time()
    T = 2.0
    it = 0
    while time.time() - t0 < seconds and cur > 0:
        it += 1
        z1 = rng.randrange(n); z2 = rng.randrange(n)
        if z1 == z2:
            continue
        # affected sigma-pairs
        prs = {min(z1, sig[z1]), min(z2, sig[z2])}
        touchedK2 = [(a, b) for z in (z1, z2) for a in far[z] for b in set(edges[phi[z1]]) | set(edges[phi[z2]])]
        touchedK3 = set()
        for p in prs:
            for zz in (p, sig[p]):
                for b in set(edges[phi[z1]]) | set(edges[phi[z2]]) | set(edges[phi[zz]]):
                    for d in L:
                        touchedK3.add((b, d)); touchedK3.add((d, b))
        before = sum(abs(K2[a][b] - K2t[a][b]) for a, b in set(touchedK2)) + \
            w3 * sum(abs(K3[b][d] - K3t[b][d]) for b, d in touchedK3) + sum(pen_pair(p, phi[p], phi[sig[p]]) for p in prs)
        # apply
        for p in prs:
            for b in edges[phi[p]]:
                for d in edges[phi[sig[p]]]:
                    K3[b][d] -= 1
            for b in edges[phi[sig[p]]]:
                for d in edges[phi[p]]:
                    K3[b][d] -= 1
        for z in (z1, z2):
            for a in far[z]:
                for b in edges[phi[z]]:
                    K2[a][b] -= 1
        phi[z1], phi[z2] = phi[z2], phi[z1]
        for z in (z1, z2):
            for a in far[z]:
                for b in edges[phi[z]]:
                    K2[a][b] += 1
        for p in prs:
            for b in edges[phi[p]]:
                for d in edges[phi[sig[p]]]:
                    K3[b][d] += 1
            for b in edges[phi[sig[p]]]:
                for d in edges[phi[p]]:
                    K3[b][d] += 1
        touchedK2b = touchedK2 + [(a, b) for z in (z1, z2) for a in far[z] for b in edges[phi[z]]]
        after = sum(abs(K2[a][b] - K2t[a][b]) for a, b in set(touchedK2b)) + \
            w3 * sum(abs(K3[b][d] - K3t[b][d]) for b, d in touchedK3) + sum(pen_pair(p, phi[p], phi[sig[p]]) for p in prs)
        before2 = before + sum(abs(K2[a][b] - K2t[a][b]) - abs(K2[a][b] - K2t[a][b]) for a, b in [])
        # recompute 'before' on the same K2 cells as 'after' (cells touched only after the swap had value unchanged before)
        extra = set(touchedK2b) - set(touchedK2)
        delta = after - before
        if extra:
            # those cells were not in 'before'; their old values: undo contribution of new images for the check
            pass
        if delta <= 0 or rng.random() < math.exp(-delta / T):
            cur += delta
            if cur < best:
                best = cur
        else:
            # undo
            for p in prs:
                for b in edges[phi[p]]:
                    for d in edges[phi[sig[p]]]:
                        K3[b][d] -= 1
                for b in edges[phi[sig[p]]]:
                    for d in edges[phi[p]]:
                        K3[b][d] -= 1
            for z in (z1, z2):
                for a in far[z]:
                    for b in edges[phi[z]]:
                        K2[a][b] -= 1
            phi[z1], phi[z2] = phi[z2], phi[z1]
            for z in (z1, z2):
                for a in far[z]:
                    for b in edges[phi[z]]:
                        K2[a][b] += 1
            for p in prs:
                for b in edges[phi[p]]:
                    for d in edges[phi[sig[p]]]:
                        K3[b][d] += 1
                for b in edges[phi[sig[p]]]:
                    for d in edges[phi[p]]:
                        K3[b][d] += 1
        if it % 20000 == 0:
            T = max(0.05, T * 0.97)
            cur = cost_all()
    final = cost_all()
    return final, best, it, phi, far, edges, sig


if __name__ == '__main__':
    cfg = json.load(open(sys.argv[1]))
    secs = float(sys.argv[2]) if len(sys.argv) > 2 else 60
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    final, best, it, phi, far, edges, sig = run(cfg, secs, seed)
    print('final cost', final, '| best seen', best, '| iterations', it, flush=True)
    if final == 0:
        sol = {str(list(far[z])): list(edges[phi[z]]) for z in range(60)}
        json.dump({'cfg': cfg, 'Phi': sol}, open(sys.argv[1].replace('.json', '_sa%d.json' % seed), 'w'), indent=1)
        print('solution saved', flush=True)
