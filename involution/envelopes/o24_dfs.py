"""O_24 (id 0), second independent route (structures 27/09): for one (y1, y2) class of C1's list, enumerate the SETS of 7
signed octads that can be the 7 h of a frame, with the u-conditions used as pruning on partial configurations.
Scaled coordinates (sqrt8): h_c = 2 M_c (M_c = +-1 on an octad, even number of -1, taken up to sign = gauge), y_j = B_j,
u = E (sign vector of a Golay word, up to sign).
Exact necessary conditions (all proven, see o24_ustage.py):
  M.B1 = M.B2 = 0;  M_c.M_d = 0;  coverage sum_c |M_c| = 1 on A = {i : B1_i^2 + B2_i^2 = 10} and 3 elsewhere;
  every u: E.B_j = 0, E.M_c in {0,+-4}, exactly two +-4;
  for every pair {c,d}: one u of type "same sign" and one of type "opposite sign";
  for every c: the 12 u of its star (so >= 2 per other block).
Pruning on a partial set P (all monotone): each E is alive iff E.B = 0, no E.M = +-8 with M in P and at most 2 nonzero
products with P; for c, d in P both types must have an alive E nonzero exactly at c, d; for c in P at least 2 (7 - |P|)
alive E nonzero only at c; counts of alive E by number of nonzero products; each deficit point still coverable.
Enumeration of SETS without repetition: branch on the lowest point with a deficit, over its smallest-index covering member
(candidates of smaller index containing that point are excluded below the branch).
Leaves are handed to o24_ustage.ustage (LP with exact certificate).  Usage: python o24_dfs.py class_index [seconds]"""
import json, sys, time
import numpy as np
from o24_ustage import C1_YPAIRS, golay_from_basis, sign_vectors, ustage

T0 = time.time()


def run(k, budget, verbose=True, max_leaf_lp=20):
    D = json.load(open(C1_YPAIRS))
    code = golay_from_basis(D['golay_basis'])
    OCT = [W for W in code if bin(W).count('1') == 8]
    orb = D['orbits'][k]
    B1, B2 = np.array(orb['a1'], dtype=np.int64), np.array(orb['a2'], dtype=np.int64)
    A = [i for i in range(24) if B1[i] ** 2 + B2[i] ** 2 == 10]
    assert sorted(A) == sorted(orb['T1'] + orb['T2'])
    target = [1 if i in A else 3 for i in range(24)]
    SV = sign_vectors(code)
    SVp = SV[(SV @ B1 == 0) & (SV @ B2 == 0)]
    nE = len(SVp)
    # signed octads up to sign: first point of the octad positive, even number of minus signs
    Ms = []
    for W in OCT:
        pts = [i for i in range(24) if W >> i & 1]
        for s in range(128):
            signs = [1] + [-1 if s >> j & 1 else 1 for j in range(7)]
            if signs.count(-1) % 2:
                continue
            v = np.zeros(24, dtype=np.int64)
            v[pts] = signs
            if v @ B1 == 0 and v @ B2 == 0:
                Ms.append(v)
    Ms = np.array(Ms)
    PE = SVp @ Ms.T                                          # nE x nM, values in {0, +-4, +-8}
    star = (np.abs(PE) == 4).sum(axis=0)
    keep = star >= 12
    Ms, PE = Ms[keep], PE[:, keep]
    nM = len(Ms)
    if verbose:
        print('class %d: |A| = %d, %d sign words perp y (up to sign), %d signed octads perp y with >= 12 star candidates'
              % (k, len(A), nE, nM), flush=True)

    def bits(boolvec):
        return int.from_bytes(np.packbits(np.asarray(boolvec, dtype=np.uint8), bitorder='little').tobytes(), 'little')
    P4 = [bits(PE[:, m] == 4) for m in range(nM)]
    N4 = [bits(PE[:, m] == -4) for m in range(nM)]
    B8 = [bits(np.abs(PE[:, m]) == 8) for m in range(nM)]
    Mf = Ms.astype(np.float32)
    ORTH = []
    for a in range(0, nM, 512):                                  # exact: entries are small integers
        Gb = np.rint(Mf[a:a + 512] @ Mf.T).astype(np.int16)
        for row in Gb:
            ORTH.append(bits(row == 0))
    INC = [bits(Ms[:, p] != 0) for p in range(24)]
    if verbose:
        print('precomputation %.1fs' % (time.time() - T0), flush=True)
    ALLE = (1 << nE) - 1
    stats = {'nodes': 0, 'leaves': 0, 'leaf_verdicts': {}, 'per_depth': [0] * 8}
    leaves = []

    def check(chosen, nz1, nz2, dead):
        alive = ALLE & ~dead
        nz0 = alive & ~nz1 & ~nz2
        L = len(chosen)
        free = 7 - L
        if nz0.bit_count() < 2 * free * (free - 1) // 2 or nz1.bit_count() < 2 * L * free:
            return False
        for a in range(L):
            c = chosen[a]
            if ((P4[c] | N4[c]) & nz1).bit_count() < 2 * free:
                return False
            for b in range(a + 1, L):
                d = chosen[b]
                if not (((P4[c] & P4[d]) | (N4[c] & N4[d])) & nz2):
                    return False
                if not (((P4[c] & N4[d]) | (N4[c] & P4[d])) & nz2):
                    return False
        return True

    def rec(chosen, cov, avail, nz1, nz2, dead):
        stats['nodes'] += 1
        stats['per_depth'][len(chosen)] += 1
        if stats['nodes'] % 2000 == 0:
            now = time.time() - T0
            if now > budget:
                raise TimeoutError
            if now > stats.get('next_report', 0):
                stats['next_report'] = now + 30
                if verbose:
                    print('  [%.0fs] nodes %d per depth %s leaves %d' % (now, stats['nodes'], stats['per_depth'],
                                                                     stats['leaves']), flush=True)
        L = len(chosen)
        if L == 7:
            if cov == target:
                stats['leaves'] += 1
                leaves.append(list(chosen))
            return
        # deficit points must remain coverable
        p = None
        for q in range(24):
            dq = target[q] - cov[q]
            if dq:
                if (avail & INC[q]).bit_count() < dq:
                    return
                if p is None:
                    p = q
        if p is None:
            return
        opts = avail & INC[p]
        excl = 0
        while opts:
            low = opts & -opts
            o = low.bit_length() - 1
            opts ^= low
            nzo = P4[o] | N4[o]
            ndead = dead | B8[o] | (nz2 & nzo)
            nnz2 = (nz2 & ~nzo) | (nz1 & nzo)
            nnz1 = (nz1 & ~nzo) | (ALLE & ~nz1 & ~nz2 & ~dead & nzo)
            nnz1 &= ~ndead
            nnz2 &= ~ndead
            ncov = list(cov)
            for q in range(24):
                if Ms[o, q]:
                    ncov[q] += 1
            if any(ncov[q] > target[q] for q in range(24)):
                excl |= low
                continue
            nch = chosen + [o]
            if check(nch, nnz1, nnz2, ndead):
                navail = avail & ORTH[o] & ~excl & ~low
                for q in range(24):
                    if ncov[q] == target[q]:
                        navail &= ~INC[q]
                rec(nch, ncov, navail, nnz1, nnz2, ndead)
            excl |= low
    try:
        rec([], [0] * 24, (1 << nM) - 1, 0, 0, 0)
        complete = True
    except TimeoutError:
        complete = False
    json.dump({'class': k, 'complete': complete, 'leaves': [Ms[lf].tolist() for lf in leaves[:2000]]},
              open('state/o24_leaves_%03d.json' % k, 'w'))
    for lf in leaves[:max_leaf_lp]:
        r = ustage(B1, B2, Ms[lf], SV)
        v = r['verdict']
        stats['leaf_verdicts'][v] = stats['leaf_verdicts'].get(v, 0) + 1
    out = {'class': k, 'complete': complete, 'nodes': stats['nodes'], 'nodes_per_depth': stats['per_depth'],
           'survivors_at_7': stats['leaves'],
           'verdict': ('EXCLUDED' if complete and not leaves else 'OPEN' if not complete else
                       'EXCLUDED' if all(v.startswith('EXCLUDED') for v in stats['leaf_verdicts']) else 'SURVIVORS'),
           'leaves': stats['leaves'],
           'leaf_verdicts': stats['leaf_verdicts'], 'n_signed_octads': nM, 'n_sign_words': nE,
           'seconds': round(time.time() - T0, 1)}
    if verbose:
        print(out, flush=True)
    return out


if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 240.0
    run(k, budget)
