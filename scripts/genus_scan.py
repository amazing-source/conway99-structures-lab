"""S1 + S2 — general form of L22 (genus of the eigenlattices Lambda_r, Lambda_s of a strongly regular graph).

Clean class (everything determined except the p-adic glue ranks alpha_q and signs delta_q):
  integer eigenvalues k > r > s, d = r - s odd, every odd prime q | d has v_q(d) = 1 and q ∤ v.
For theta in {r, s} (theta' the other one, m = mult(theta), m' = mult(theta')):
  * q | d:           Lambda_theta ⊗ Z_q = U ⊥ qV, rank V = alpha_q (same alpha_q for both sides), 0 <= alpha_q <= min(m, m');
                     V' ≅ -V mod q, so delta'_q = delta_q + alpha_q [q ≡ 3 mod 4];  eps(U) eps(V) = (det / q^alpha | q).
  * q ∤ d, q | v:    the disc of 1^⊥ is cyclic, so exactly one of Lambda_r, Lambda_s carries Z/q^{v_q(v)} (S2);
                     it is Lambda_theta iff v_q(k - theta') < v_q(v); its component q^e<u> has eps = (-v/q^e | q).
  * 2-adic:          if v_2(k - theta') >= v_2(v), Lambda_theta is even and 2-unimodular: m even,
                     (-1)^{m/2} det ≡ 1 (mod 4), and the oddity formula
                     m + sum_{q|d} [alpha_q (q-1) + 4 delta_q] + sum_cyc [(q^e - 1) + 4 [e odd, eps = -1]] ≡ 0 (mod 8).
Output: parameter sets with no admissible alpha (would be non-existence), ceiling alpha = min(f,g) excluded, parity laws.
"""
import itertools
import math
import sys


def vp(x, p):
    x = abs(x)
    if x == 0:
        return 10 ** 9
    c = 0
    while x % p == 0:
        x //= p
        c += 1
    return c


def primes(x):
    x, out, p = abs(x), [], 2
    while p * p <= x:
        if x % p == 0:
            out.append(p)
            while x % p == 0:
                x //= p
        p += 1
    if x > 1:
        out.append(x)
    return out


def leg(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def feasible(v, k, l, mu):
    if k * (k - l - 1) != mu * (v - k - 1):
        return None
    D = (l - mu) ** 2 + 4 * (k - mu)
    sd = math.isqrt(D)
    if sd * sd != D or (l - mu + sd) % 2:
        return None
    r, s = (l - mu + sd) // 2, (l - mu - sd) // 2
    num_f, den = -k * (s + 1) * (k - s), (k + r * s) * (r - s)
    if den == 0 or num_f % den:
        return None
    f = num_f // den
    g = v - 1 - f
    if f <= 0 or g <= 0:
        return None
    if (r + 1) * (k + r + 2 * r * s) > (k + r) * (s + 1) ** 2 or (s + 1) * (k + s + 2 * r * s) > (k + s) * (r + 1) ** 2:
        return None
    if v > f * (f + 3) // 2 or v > g * (g + 3) // 2:
        return None
    return r, s, f, g


def s2_violations(v, k, r, s):
    """Primes q ∤ r-s, q | v for which NOT exactly one of E_r, E_s is q-integral."""
    bad = []
    for q in primes(v):
        if (r - s) % q == 0:
            continue
        ir = vp(k - s, q) >= vp(v, q)          # E_r q-integral
        is_ = vp(k - r, q) >= vp(v, q)         # E_s q-integral
        if ir == is_:
            bad.append((q, ir, is_))
    return bad


def side_data(v, k, theta, other, m):
    """Cyclic components of Lambda_theta at q ∤ d (list of (q, e, eps)) and whether Lambda_theta is 2-unimodular."""
    d = abs(theta - other)
    cyc = []
    for q in primes(v):
        if q == 2 or d % q == 0:
            continue
        if vp(k - other, q) < vp(v, q):
            e = vp(v, q)
            cyc.append((q, e, leg(-(v // q ** e), q)))
    two_unimod = vp(k - other, 2) >= vp(v, 2)
    return cyc, two_unimod


def side_ok(m, qs, alphas, deltas, cyc, two_unimod, mm):
    """Constraints on one side. mm = min(f,g) (only used for the forced sign when alpha = m)."""
    det_unit = 1
    for q, a in zip(qs, alphas):
        det_unit *= q ** a
    for q, e, _ in cyc:
        det_unit *= q ** e
    for q, a, dl in zip(qs, alphas, deltas):
        if a == 0 and dl:
            return False
        if a == m:                                   # U = 0: eps(V) = (det / q^a | q)
            rest = det_unit // q ** a
            if dl != (leg(rest, q) == -1):
                return False
    if not two_unimod:
        return True                                  # 2-adic part not determined: no oddity constraint
    if m % 2:
        return False
    if ((-1) ** (m // 2) * det_unit) % 4 != 1:
        return False
    tot = m
    for q, a, dl in zip(qs, alphas, deltas):
        tot += a * (q - 1) + 4 * dl
    for q, e, eps in cyc:
        tot += (q ** e - 1) + (4 if (e % 2 and eps == -1) else 0)
    return tot % 8 == 0


def analyse(v, k, l, mu):
    fe = feasible(v, k, l, mu)
    if not fe:
        return None
    r, s, f, g = fe
    d = r - s
    res = dict(params=(v, k, l, mu), eig=(r, s, f, g), d=d, s2=s2_violations(v, k, r, s))
    if d % 2 == 0:
        res['status'] = 'skip: d even'
        return res
    qs = primes(d)
    if any(vp(d, q) > 1 or v % q == 0 for q in qs):
        res['status'] = 'skip: q^2 | d or q | gcd(d, v)'
        return res
    mm = min(f, g)
    cyc_r, tu_r = side_data(v, k, r, s, f)
    cyc_s, tu_s = side_data(v, k, s, r, g)
    res.update(cyc_r=cyc_r, cyc_s=cyc_s, two_unimod=(tu_r, tu_s))
    # candidate alphas per q: boundaries and one representative per residue mod 8 in the interior
    cands = []
    for q in qs:
        c = {0, mm}
        for a in range(1, mm):
            if len(c) >= 10:
                break
            c.add(a)
        # make sure all residues mod 8 are represented (for large mm the loop above covers 1..8)
        cands.append(sorted(c))
    allowed = set()
    for alphas in itertools.product(*cands):
        for deltas in itertools.product((0, 1), repeat=len(qs)):
            deltas_s = deltas
            deltas_r = tuple((dl + a * (q % 4 == 3)) % 2 for q, a, dl in zip(qs, alphas, deltas))
            if side_ok(g, qs, alphas, deltas_s, cyc_s, tu_s, mm) and side_ok(f, qs, alphas, deltas_r, cyc_r, tu_r, mm):
                allowed.add(alphas)
                break
    res['qs'] = qs
    res['mm'] = mm
    res['allowed'] = allowed
    res['cands'] = cands

    def interior_residues(i):
        return sorted({a[i] % 8 for a in allowed if 0 < a[i] < mm})

    res['ceiling_excluded'] = [not any(a[i] == mm for a in allowed) for i in range(len(qs))]
    res['residues'] = [interior_residues(i) for i in range(len(qs))]
    res['status'] = 'EXCLUDED (no admissible alpha)' if not allowed else 'ok'
    return res


def scan(vmax):
    """Enumerate by eigenvalues: mu = k + r s, lambda = mu + r + s, v = 1 + k + k (k - lambda - 1) / mu."""
    out, seen = [], set()
    for s in range(-1, -vmax, -1):
        for r in range(0, vmax // max(1, -s) + 2):       # v > k >= 1 - r s, so -r s < vmax
            for k in range(max(r + 1, -r * s + 1), vmax):
                mu = k + r * s
                l = mu + r + s
                if mu < 1 or l < 0 or mu > k:
                    continue
                num = k * (k - l - 1)
                if num % mu:
                    continue
                v = 1 + k + num // mu
                if v > vmax:            # v is not monotone in k (mu = k + rs grows): keep scanning
                    continue
                if 2 * k > v - 1 or (v, k, l, mu) in seen:
                    continue
                seen.add((v, k, l, mu))
                res = analyse(v, k, l, mu)
                if res:
                    out.append(res)
    return out


if __name__ == '__main__':
    VMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1300
    print('srg(99,14,1,2):', {kk: vv for kk, vv in analyse(99, 14, 1, 2).items() if kk != 'allowed'})
    allres = scan(VMAX)
    clean = [x for x in allres if x['status'] in ('ok', 'EXCLUDED (no admissible alpha)')]
    s2bad = [x for x in allres if x['s2']]
    print(f'\nfeasible parameter sets (v <= {VMAX}, k <= (v-1)/2): {len(allres)}; clean class: {len(clean)}')
    print(f'S2 violations (not exactly one of E_r, E_s q-integral, q ∤ d, q | v): {len(s2bad)}')
    for x in s2bad[:40]:
        print('   srg%s eig %s  bad %s' % (x['params'], x['eig'], x['s2']))
    ex = [x for x in clean if x['status'].startswith('EXCLUDED')]
    print(f'\nclean class, NO admissible alpha: {len(ex)}')
    for x in ex:
        print('   srg%s eig r=%d^%d s=%d^%d  d=%d  cyc_r=%s cyc_s=%s 2-unimod=%s' %
              (x['params'], x['eig'][0], x['eig'][2], x['eig'][1], x['eig'][3], x['d'], x['cyc_r'], x['cyc_s'], x['two_unimod']))
    ceil = [x for x in clean if x['status'] == 'ok' and any(x['ceiling_excluded'])]
    print(f'\nclean class, ceiling alpha = min(f,g) excluded (L22 type): {len(ceil)}')
    for x in ceil:
        print('   srg%s eig r=%d^%d s=%d^%d  d=%d  q=%s  interior residues mod 8 allowed=%s' %
              (x['params'], x['eig'][0], x['eig'][2], x['eig'][1], x['eig'][3], x['d'], x['qs'], x['residues']))
