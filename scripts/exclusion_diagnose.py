"""Why does each excluded parameter set fail?  Level 1: rank parity (m odd with Lambda 2-unimodular).
Level 2: discriminant mod 4, (-1)^{m/2} det ≡ 1 (mod 4), independent of alpha when every q | d with q ≡ 3 (mod 4) is absent.
Level 3: only the oddity formula mod 8 kills it."""
import sys, itertools
sys.path.insert(0, 'scripts')
import genus_scan as gs

def diagnose(v, k, l, mu):
    res = gs.analyse(v, k, l, mu)
    r, s, f, g = res['eig']
    qs = res['qs']
    out = []
    for name, m, cyc, tu in (('r', f, res['cyc_r'], res['two_unimod'][0]), ('s', g, res['cyc_s'], res['two_unimod'][1])):
        if not tu:
            out.append(f'Lambda_{name}: 2-adic part not determined'); continue
        if m % 2:
            out.append(f'Lambda_{name}: rank {m} ODD but even 2-unimodular'); continue
        cycprod = 1
        for q, e, _ in cyc: cycprod *= q ** e
        free = [q for q in qs if q % 4 == 3]
        vals = set()
        for par in itertools.product((0, 1), repeat=len(free)):
            dd = cycprod
            for q, a in zip(free, par): dd *= q ** a
            vals.add(((-1) ** (m // 2) * dd) % 4)
        if vals == {3}:
            out.append(f'Lambda_{name}: disc mod 4 FAILS for every alpha (cyclic part {cyc}, (-1)^(m/2) = {(-1)**(m//2)})')
        elif 3 in vals:
            out.append(f'Lambda_{name}: disc mod 4 forces the parity of alpha at {free}')
        else:
            out.append(f'Lambda_{name}: disc mod 4 OK for every alpha')
    return res, out

if __name__ == '__main__':
    exc = [(430,39,8,3),(616,205,90,57),(209,52,15,12),(235,52,9,12),(1026,451,212,187),(494,85,12,15),
           (1148,481,210,195),(630,68,1,8),(871,348,137,140),(1248,203,22,35),(780,369,158,189),(690,318,121,168)]
    for p in exc:
        res, out = diagnose(*p)
        r, s, f, g = res['eig']
        print(f'srg{p}: r={r}^{f} s={s}^{g} d={res["d"]} q|d={res["qs"]}')
        for o in out: print('    ' + o)
