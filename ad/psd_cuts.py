"""Pair-level model + PSD cutting planes from the orbit Gram P = 6I - 2Q - Sh + J (PSD, = 14 x projection of rank 15).
Its compression to pair-constant vectors, P_Sym[P,P'] = 6 d - 2 T - 2|P n P'| + 2, with T[P,P'] = (phi_{P+}(P') + phi_{P-}(P'))/2
(P != P') and T[P,P] = n_s(P), must satisfy 0 <= P_Sym <= 14 I.  Each violated eigenvector v gives a linear cut
v^T P_Sym(phi) v >= 0 (or <= 14|v|^2).  Loop: solve ILP, cut, resolve."""
import numpy as np, itertools, time, sys
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp
from general_model import build, PAIRS, N
PI = {P: i for i, P in enumerate(PAIRS)}
def psym_from(x, orbs, phi, ns):
    by_pair = {P: [i for i, o in enumerate(orbs) if o['pair'] == P] for P in PAIRS}
    M = np.zeros((21, 21))
    for P in PAIRS:
        for P2 in PAIRS:
            a, b = PI[P], PI[P2]
            inter = len(set(P) & set(P2))
            if P == P2:
                nsv = x[ns[P]] if P in ns else 0.0
                M[a, b] = 4 - 2 * nsv
            else:
                T = sum(x[phi[(i, P2)]] for i in by_pair[P]) / 2.0
                M[a, b] = -2 * T - 2 * inter + 2
    return (M + M.T) / 2
def cut_row(v, orbs, phi, ns, nvar):
    """linear form in the variables equal to v^T P_Sym v, as (coef vector, constant)"""
    by_pair = {P: [i for i, o in enumerate(orbs) if o['pair'] == P] for P in PAIRS}
    coef = np.zeros(nvar); const = 0.0
    for P in PAIRS:
        for P2 in PAIRS:
            a, b = PI[P], PI[P2]; w = v[a] * v[b]
            if P == P2:
                const += 4 * w
                if P in ns: coef[ns[P]] += -2 * w
            else:
                const += (-2 * len(set(P) & set(P2)) + 2) * w
                # -2 T, T = (sum phi_{P,s}(P2) + sum phi_{P2,s}(P))/4 symmetrised
                for i in by_pair[P]: coef[phi[(i, P2)]] += -2 * w / 2.0
    return coef, const
def run(orbs, near=None, dichotomy=True, rounds=60, tl=60, verbose=False, integral=True):
    A, lo, hi, lb, ub, nvar, phi, ns = build(orbs, dichotomy=dichotomy, near=near)
    A = A.tocsr(); cuts = []; cuts_lo = []; cuts_hi = []
    for rnd in range(rounds):
        cons = [LinearConstraint(A, lo, hi)]
        if cuts:
            C = np.array(cuts); cons.append(LinearConstraint(C, np.array(cuts_lo), np.array(cuts_hi)))
        integ = np.ones(nvar) if integral else np.zeros(nvar)
        r = milp(np.zeros(nvar), constraints=cons, integrality=integ, bounds=Bounds(lb, ub), options={'time_limit': tl})
        if r.status == 2: return 'dead', rnd
        if r.status != 0: return 'unknown', rnd
        M = psym_from(r.x, orbs, phi, ns)
        ev, V = np.linalg.eigh(M)
        added = 0
        for k in range(len(ev)):
            if ev[k] < -1e-7:
                c, k0 = cut_row(V[:, k], orbs, phi, ns, nvar); cuts.append(c); cuts_lo.append(-k0); cuts_hi.append(np.inf); added += 1
            if ev[k] > 14 + 1e-7:
                c, k0 = cut_row(V[:, k], orbs, phi, ns, nvar); cuts.append(c); cuts_lo.append(-np.inf); cuts_hi.append(14 * V[:, k] @ V[:, k] - k0); added += 1
        if verbose: print('round', rnd, 'min eig %.3f max eig %.3f' % (ev[0], ev[-1]), 'cuts added', added, flush=True)
        if not added: return 'alive_psd', rnd
    return 'unknown_rounds', rounds
if __name__ == '__main__':
    import random
    from alldd import skeleton_from_pair_cycles, random_cycles
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3):
        orbs = skeleton_from_pair_cycles(random_cycles(rng))
        mode = sys.argv[3] if len(sys.argv) > 3 else 'lp'
        t0 = time.time(); res = run(orbs, verbose=True, dichotomy=(mode != 'lp_nodich'), integral=(mode == 'ilp'), rounds=200)
        print('all-dd trial', t, res, '%.0fs' % (time.time() - t0), flush=True)
