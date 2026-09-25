"""S5 — simplified criterion (no Jordan data needed):
srg with integer eigenvalues r > s, d = r - s ODD.  For theta with v_2(k - theta') >= v_2(v), Lambda_theta is even with odd
determinant D, so  (a) m_theta is even,  (b) (-1)^{m/2} D ≡ 1 (mod 4).
D = (q-parts, q | d) * prod_{q | v, q ∤ d, v_q(k - theta') < v_q(v)} q^{v_q(v)}.  The q | d parts are ≡ 1 (mod 4) when
q ≡ 1 (mod 4); if some q | d is ≡ 3 (mod 4) its exponent is unknown and (b) is not used (only (a))."""
import sys
sys.path.insert(0, 'scripts')
import genus_scan as gs

def crit(v, k, l, mu):
    fe = gs.feasible(v, k, l, mu)
    if not fe: return None
    r, s, f, g = fe
    d = r - s
    if d % 2 == 0: return None
    qd = gs.primes(d)
    det_free = all(q % 4 == 1 for q in qd)
    out = []
    for th, oth, m in ((r, s, f), (s, r, g)):
        if gs.vp(k - oth, 2) < gs.vp(v, 2): continue
        if m % 2:
            out.append((th, m, 'rank odd')); continue
        if not det_free: continue
        P = 1
        for q in gs.primes(v):
            if q == 2 or d % q == 0: continue
            if gs.vp(k - oth, q) < gs.vp(v, q): P *= q ** gs.vp(v, q)
        if ((-1) ** (m // 2) * P) % 4 != 1:
            out.append((th, m, f'disc mod 4: (-1)^{m//2} * {P} ≡ 3'))
    return (r, s, f, g, d, out)

if __name__ == '__main__':
    V = int(sys.argv[1]) if len(sys.argv) > 1 else 1300
    hits = []
    n_all = n_odd = 0
    for x in gs.scan.__wrapped__(V) if hasattr(gs.scan, '__wrapped__') else []:
        pass
    # reuse the eigenvalue enumeration of genus_scan but only for the parameter tuples
    seen = set()
    for s in range(-1, -V, -1):
        for r in range(0, V // max(1, -s) + 2):
            for k in range(max(r + 1, -r * s + 1), V):
                mu = k + r * s; l = mu + r + s
                if mu < 1 or l < 0 or mu > k: continue
                num = k * (k - l - 1)
                if num % mu: continue
                v = 1 + k + num // mu
                if v > V or 2 * k > v - 1 or (v, k, l, mu) in seen: continue
                seen.add((v, k, l, mu))
                c = crit(v, k, l, mu)
                if c is None: continue
                n_odd += 1
                if c[5]: hits.append(((v, k, l, mu), c))
    print(f'feasible sets with d odd (v <= {V}): {n_odd};  excluded by (a) or (b): {len(hits)}')
    for p, c in sorted(hits):
        r, s, f, g, d, out = c
        print(f'  srg{p}  r={r}^{f} s={s}^{g} d={d}  ' + '; '.join(f'Lambda_{th}(m={m}): {why}' for th, m, why in out))
