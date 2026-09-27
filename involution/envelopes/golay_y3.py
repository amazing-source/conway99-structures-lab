"""id 1 (a1^8), step 3 (structures 27/09): pair classes (y1, y2) up to G1 and frame LP on their representatives.
For each y1 class k: the filtered y2 (golay_y2 logic), the stabiliser of +-y1 in G1 computed exactly by listing the
322 560 elements of Stab_M24(K) and keeping (pi, eps) with eps o pi (y1) = +-y1 (eps: an inner codeword on K^c, anything on K),
plus the free sign flips on K \\ supp(y1); then orbits of y2 under Stab(+-y1) and y2 -> -y2 (union-find on sorted keys).
Each representative pair goes to frame_lp (Golay-frame lattice, a-coordinates, scale 8).
Usage: python golay_y3.py k [k ...] [--nolp]"""
import itertools, json, os, sys, time, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from golay_family import WORDS, OM, M24, setwise_stabilizer, compose
from golay_y import lattice_data

T0 = time.time()


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def group_elements(gens):
    """all elements of <gens> (permutations of 24 points) by BFS closure"""
    idp = tuple(range(24))
    seen = {idp}
    frontier = [idp]
    while frontier:
        nf = []
        for p in frontier:
            for g in gens:
                q = tuple(g[p[i]] for i in range(24))
                if q not in seen:
                    seen.add(q)
                    nf.append(q)
        frontier = nf
    return np.array(sorted(seen), dtype=np.int8)


def y2_set(y1, Yn, Kl, Kcl, Wsign, EPS, nzK, BETA, N3):
    b1 = (BETA @ y1).astype(np.int64)
    n1 = np.abs(N3 @ y1)
    U1 = N3[n1 == 24].T.astype(np.int32)
    B32 = BETA.T.astype(np.int32)
    SKc = (Yn[:, Kcl] * y1[Kcl]) @ Wsign.T
    prodK = Yn[:, Kl] * y1[Kl]
    rootok = ~((np.abs(Yn[:, Kl]) == 4) & (np.abs(y1[Kl]) == 4)[None, :]).any(axis=1)
    out = []
    for s0 in range(0, len(Yn), 256):
        sl = np.nonzero(rootok[s0:s0 + 256])[0] + s0
        if not len(sl):
            continue
        sK = prodK[sl] @ EPS.T
        valid = ~((EPS[None, :, :] == -1) & ~nzK[sl][:, None, :]).any(axis=2)
        match = (SKc[sl][:, None, :] == -sK[:, :, None]) & valid[:, :, None]
        ii, ee, ww = np.nonzero(match)
        if not len(ii):
            continue
        Y2 = Yn[sl][ii].copy()
        Y2[:, Kl] *= EPS[ee]
        Y2[:, Kcl] *= Wsign[ww]
        ok = ~((np.abs(Y2[:, Kl]) == 4) & (np.abs(y1[Kl]) == 4)[None, :]).any(axis=1)
        Y2 = Y2[ok]
        if not len(Y2):
            continue
        bb = (Y2 @ B32).astype(np.int64)
        ok = ((bb ** 2 + b1[None, :] ** 2) <= 160).all(axis=1)
        if U1.shape[1]:
            ok &= ~(np.abs(Y2 @ U1) == 24).any(axis=1)
        out.append(Y2[ok].astype(np.int8))
    return np.concatenate(out) if out else np.zeros((0, 24), dtype=np.int8)


def keys_of(Y):
    """two int64 keys per int8 row (base 11 on 12 + 12 coordinates), for sorted lookup"""
    Z = (Y.astype(np.int64) + 5)
    w = 11 ** np.arange(12, dtype=np.int64)
    return Z[:, :12] @ w, Z[:, 12:] @ w


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    K, A, N3, BETA = lattice_data()
    Kc = OM ^ K
    Kl = [i for i in range(24) if K >> i & 1]
    Kcl = [i for i in range(24) if Kc >> i & 1]
    inner = [w & Kc for w in WORDS if w & K == 0]
    inner_set = set(inner)
    Wsign = np.array([[(-1) ** (w >> j & 1) for j in Kcl] for w in inner], dtype=np.int32)
    EPS = np.array(list(itertools.product((1, -1), repeat=8)), dtype=np.int32)
    Yn = np.load(os.path.join(HERE, 'state', 'u1_ynorm.npy')).astype(np.int32)
    nzK = (Yn[:, Kl] != 0)
    reps = json.load(open(os.path.join(HERE, 'state', 'u1_y1reps.json')))['reps']
    Sperm, _ = setwise_stabilizer(M24, K, samples=40)
    H = group_elements(Sperm[:8])
    log('Stab_M24(K) listed: %d elements' % len(H))
    assert len(H) == 322560
    P = None
    if '--nolp' not in sys.argv:
        import frame_lp as FL
        P = FL.prepare([list(map(int, r)) for r in A], scale=8, verbose=False)
    results = {}
    for k in [int(a) for a in args]:
        y1 = np.array(reps[k], dtype=np.int32)
        Y2 = y2_set(y1, Yn, Kl, Kcl, Wsign, EPS, nzK, BETA, N3)
        log('y1 class %d: %d filtered y2' % (k, len(Y2)))
        # stabiliser of +-y1: pi(y1)_{pi(i)} = y1_i ; need eps o pi(y1) = s y1, eps on K^c = inner codeword pattern
        img = np.zeros((len(H), 24), dtype=np.int32)
        rows = np.arange(len(H))[:, None]
        img[rows, H.astype(np.int64)] = y1[None, :]
        same_abs = (np.abs(img) == np.abs(y1)[None, :]).all(axis=1)
        stab = []
        for s in (1, -1):
            ratio = np.where(img[:, Kcl] * (s * y1[Kcl])[None, :] > 0, 0, 1)            # 1 where the sign must flip
            masks = (ratio * (1 << np.array(Kcl, dtype=np.int64))[None, :]).sum(axis=1)
            for idx in np.nonzero(same_abs)[0]:
                if int(masks[idx]) in inner_set:
                    # eps: flip K^c on the codeword mask, K signs to match s*y1 where y1 != 0
                    eps = np.ones(24, dtype=np.int32)
                    for j in Kcl:
                        if masks[idx] >> j & 1:
                            eps[j] = -1
                    for j in Kl:
                        if y1[j] != 0 and img[idx, j] * s * y1[j] < 0:
                            eps[j] = -1
                    stab.append((H[idx].astype(np.int64), eps))
        free = [j for j in Kl if y1[j] == 0]
        log('   Stab(+-y1): %d (pi, eps) elements, plus %d free K sign flips' % (len(stab), len(free)))
        # union-find on y2 under a generating sample of Stab and the free flips and y2 -> -y2
        def bkeys(Z):
            return np.ascontiguousarray((Z.astype(np.int16) + 6).astype(np.uint8)).view('S24').ravel()   # no zero bytes
        KY = bkeys(Y2)
        order = np.argsort(KY, kind='stable')
        SK = KY[order]

        def lookup(Z):
            q = bkeys(Z)
            pos = np.searchsorted(SK, q)
            pos = np.minimum(pos, len(SK) - 1)
            hit = SK[pos] == q
            res = np.where(hit, order[pos], -1)
            return res.astype(np.int64)
        rng = random.Random(k)
        gens = [stab[i] for i in rng.sample(range(len(stab)), min(12, len(stab)))]
        lab = np.arange(len(Y2), dtype=np.int64)
        maps = []
        for pi, eps in gens:
            Z = np.zeros_like(Y2)
            Z[:, pi] = Y2
            Z = (Z.astype(np.int32) * eps[None, :]).astype(np.int8)
            maps.append(lookup(Z))
        for j in free:
            Z = Y2.copy(); Z[:, j] = -Z[:, j]
            maps.append(lookup(Z))
        maps.append(lookup(-Y2))
        for f in maps:
            assert (f >= 0).all(), 'a stabiliser element leaves the y2 set'
        for it in range(100):
            old = lab.copy()
            for f in maps:
                lab = np.minimum(lab, lab[f])
                np.minimum.at(lab, f, lab.copy())
            for _ in range(5):
                lab = lab[lab]
            if (lab == old).all():
                break
        roots_ = np.unique(lab)
        log('   y2 orbits under Stab(+-y1) x {+-1}: %d' % len(roots_))
        res = {'y1': y1.tolist(), 'n_y2': int(len(Y2)), 'n_classes': int(len(roots_)), 'verdicts': {}}
        json.dump({'id': 1, 'frame': 'Golay A1^24 frame, a in Z^24, norm a.a/8 (golay_family.py)', 'scale': 8,
                   'basis': [[str(int(x)) for x in r] for r in A], 'y1_class': k,
                   'pairs': [[[str(int(v)) for v in y1], [str(int(v)) for v in Y2[r]]] for r in roots_],
                   'stats': {'n_y2_filtered': int(len(Y2)), 'stab_elements': len(stab), 'free_flips': len(free)}},
                  open(os.path.join(HERE, 'state', 'pairs_u1c%02d_001.json' % k), 'w'))
        if P is not None:
            import frame_lp as FL
            y1c = FL.coef_of(P, [int(v) for v in y1])
            for t, r in enumerate(roots_):
                y2c = FL.coef_of(P, [int(v) for v in Y2[r]])
                out = FL.frame_lp(P, y1c, y2c, save=os.path.join(HERE, 'frame_certs', 'u1_%02d_%06d.npz' % (k, t)))
                res['verdicts'][str(t)] = {'y2': Y2[r].tolist(), 'verdict': out['verdict'], 'phase1': out.get('phase1')}
                if (t + 1) % 50 == 0:
                    ex = sum(1 for v in res['verdicts'].values() if v['verdict'].startswith('EXCLUDED'))
                    log('   LP %d/%d, excluded %d' % (t + 1, len(roots_), ex))
            ex = sum(1 for v in res['verdicts'].values() if v['verdict'].startswith('EXCLUDED'))
            log('   y1 class %d: %d pair classes, %d EXCLUDED, %d not' % (k, len(roots_), ex, len(roots_) - ex))
        results[str(k)] = res
        json.dump(res, open(os.path.join(HERE, 'state', 'u1_classes_%02d.json' % k), 'w'))


if __name__ == '__main__':
    main()
