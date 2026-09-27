"""Parity + connectivity + trap test on an explicit envelope (structures 27/09).  All ingredients proven:
  (PAR) m(u) = r/8 - E(u) is EVEN for every u (Sol 27/09, checked by C2): L(u) = {chi_8 : |(chi,u)| = 3} is stable under
        chi -> chi - 2 sgn(chi.u) u and +-, orbits of size 4, and |L(u)| = 2 m(u).  E(x) = #{alpha > 0 : |(alpha,x)| = 1}.
  (CON) the 42 u are connected under |(u,u')| = 1 (C1's connectivity lemma, checked by structures).
  (P)   every w with |w|^2 < 8/7 and (w,u) in Z for all u lies in N-perp (dim 2).
Test: admissible = norm-3 x with m(x) even; connected components of the graph |(x,x')| = 1 on admissible vectors (up to
sign); a component can carry the u only if its span has rank >= 22 and (short vectors of its dual inside the span, norm
< 8/7) + (24 - rank) span at most 2 dimensions.  If no component qualifies: EXCLUDED.
Usage: python parity_trap.py id [id ...]   (results appended to parity_trap_results.json)"""
import json, os, sys, time
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build156 import Lattice, enum, rows_frac, hnf, gram_of, inv_exact, lll

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'parity_trap_results.json')


def dual_short(Bint, d):
    """rows Bint/d = a basis of a sublattice (rank k) of Q^24; returns (dimension spanned by the dual vectors of norm < 8/7,
    number of such vectors)"""
    Mr, Gs = lll(Bint, d)
    Gi = inv_exact(Gs)                                   # Gram of the dual basis (inside the span)
    den = 1
    for row in Gi:
        for x in row:
            den = den * x.denominator // np.gcd(den, x.denominator)
    Gd = [[int(x * den) for x in row] for row in Gi]    # integral after scaling by den
    bound = F(8, 7) * den
    sv = enum(Gd, float(bound) - 1e-9)                  # enum takes a float bound; exact filter below
    sv = [(z, v) for z, v in sv if v < bound]
    if not sv:
        return 0, 0
    Z = np.array([z for z, v in sv], dtype=float)
    return int(np.linalg.matrix_rank(Z)), 2 * len(sv)


def test(i):
    t0 = time.time()
    rec = json.load(open(os.path.join(HERE, 'U', '%03d.json' % i)))
    L = Lattice(rows_frac(rec['gens']))
    G = np.array(L.G, dtype=np.int64)
    sv = enum(L.G, 3)
    R = np.array([z for z, v in sv if v == 2], dtype=np.int64).reshape(-1, 24)
    X = np.array([z for z, v in sv if v == 3], dtype=np.int64)
    r = 2 * len(R)
    E = ((np.abs(X @ G @ R.T) == 1).sum(axis=1) if len(R) else np.zeros(len(X), dtype=np.int64))
    m = r // 8 - E
    assert r % 8 == 0 and (m >= 0).all()
    adm = X[m % 2 == 0]
    res = {'id': i, 'roots': rec['roots'], 'r': r, 'norm3': 2 * len(X), 'admissible': 2 * len(adm)}
    if len(adm) == 0:
        res['verdict'] = 'EXCLUDED (no norm-3 vector with m even)'
        return res
    P = adm @ G @ adm.T
    A = (np.abs(P) == 1)
    seen = np.zeros(len(adm), dtype=bool)
    comps = []
    for s in range(len(adm)):
        if seen[s]:
            continue
        stack, comp = [s], []
        seen[s] = True
        while stack:
            a = stack.pop()
            comp.append(a)
            for b in np.nonzero(A[a] & ~seen)[0]:
                seen[b] = True
                stack.append(b)
        comps.append(comp)
    info = []
    alive = False
    for comp in sorted(comps, key=len, reverse=True):
        rank = int(np.linalg.matrix_rank(adm[comp].astype(float)))
        entry = {'size': 2 * len(comp), 'rank': rank}
        if rank >= 22:
            coords = [L.coords(list(map(int, v))) for v in adm[comp]]
            Mi, d = [], 1
            for v in coords:
                for x in v:
                    d = d * x.denominator // np.gcd(d, x.denominator)
            Mi = [[int(x * d) for x in v] for v in coords]
            Bh = hnf(Mi)
            dim_short, nshort = dual_short(Bh, d)
            entry.update(dim_short=dim_short, short=nshort, codim=24 - rank)
            entry['ok'] = dim_short + 24 - rank <= 2
            alive |= entry['ok']
        else:
            entry['ok'] = False
        info.append(entry)
        if len(info) >= 6:
            break
    res['components'] = len(comps)
    res['largest'] = info
    res['verdict'] = 'survives' if alive else 'EXCLUDED (parity + connectivity + trap)'
    res['seconds'] = round(time.time() - t0, 1)
    return res


def main():
    done = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for a in sys.argv[1:]:
        i = int(a)
        if str(i) in done:
            continue
        res = test(i)
        done[str(i)] = res
        json.dump(done, open(OUT, 'w'), indent=0)
        print(res, flush=True)


if __name__ == '__main__':
    main()
