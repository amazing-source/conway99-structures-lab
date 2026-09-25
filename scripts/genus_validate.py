"""S1 validation (D11): on REAL graphs, compute Lambda_theta = Z^v ∩ ker(A - theta I) exactly and compare with the
predictions of genus_scan.py:
  * det Lambda_theta = prod_{q|d} q^alpha_q * prod_cyc q^e, 13-part elementary etc.;
  * Jordan symbols: at q | d, (0, m - alpha, eps_U), (1, alpha, eps_V); at a cyclic prime, (0, m - 1, .), (e, 1, eps_pred);
  * delta'_q = delta_q + alpha_q [q ≡ 3 mod 4];
  * the actual (alpha, delta) passes side_ok on both sides (never exclude a real graph);
  * alpha_q = rank_q(A - theta I) - [q ∤ k - theta] (L2 mechanism).
Graphs: Petersen, T(n) for odd n (d = n - 2), all complements implicitly (same lattices).
"""
import itertools
import sys
from fractions import Fraction as Fr

sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c1_recherche\code')
sys.path.insert(0, r'C:\Users\RESEARCHER\Desktop\c99_structures_lab\scripts')
import genus_tools as gt
from eigenlattice_check import int_kernel, rank_mod
import genus_scan as gs


def triangular(n):
    V = list(itertools.combinations(range(n), 2))
    return [[1 if (a != b and set(a) & set(b)) else 0 for b in V] for a in V]


def petersen():
    V = list(itertools.combinations(range(5), 2))
    return [[1 if not set(a) & set(b) else 0 for b in V] for a in V]


def lattice(A, theta):
    n = len(A)
    M = [[A[i][j] - (theta if i == j else 0) for j in range(n)] for i in range(n)]
    B = int_kernel(M)
    return [[sum(u[i] * w[i] for i in range(n)) for w in B] for u in B]


def check(name, A):
    v = len(A)
    k = sum(A[0])
    l = sum(A[0][j] * A[1][j] for j in range(v)) if A[0][1] else None
    # find lambda, mu from an adjacent and a non-adjacent pair
    adj = next(j for j in range(1, v) if A[0][j])
    non = next(j for j in range(1, v) if not A[0][j])
    lam = sum(A[0][t] * A[adj][t] for t in range(v))
    mu = sum(A[0][t] * A[non][t] for t in range(v))
    res = gs.analyse(v, k, lam, mu) if 2 * k <= v - 1 else gs.analyse(v, v - 1 - k, v - 2 - 2 * k + mu, v - 2 * k + lam)
    r, s, f, g = gs.feasible(v, k, lam, mu)
    d = r - s
    print(f'== {name}: srg({v},{k},{lam},{mu}), r={r}^{f}, s={s}^{g}, scan status: {res["status"]}')
    if not res['status'].startswith(('ok', 'EXCLUDED')):
        print('   (outside the clean class: only lattice facts printed)')
    data = {}
    for theta, other, m in ((r, s, f), (s, r, g)):
        G = lattice(A, theta)
        D = gt.det_frac(G)
        even = all(G[i][i] % 2 == 0 for i in range(len(G)))
        syms = {p: gt.jordan_symbol(G, p) for p in gt.odd_primes(D)}
        disc4 = ((-1) ** (m // 2) * D) % 4 if D % 2 else None
        print(f'   Lambda_{theta}: rank {len(G)}, det {D} = {gs.primes(int(D))}, even {even}, '
              f'(-1)^(m/2) det mod 4 = {disc4}, Jordan {syms}')
        data[theta] = (G, D, syms)
    if not res['status'].startswith(('ok', 'EXCLUDED')):
        return True
    ok = True
    qs = res['qs']
    alphas, deltas = [], {r: [], s: []}
    for q in qs:
        rk = rank_mod([[A[i][j] - (r if i == j else 0) for j in range(v)] for i in range(v)], q)
        a_pred = rk - (1 if (k - r) % q else 0)
        comps = {j: (nn, eps) for j, nn, eps in data[s][2].get(q, [])}
        a_true = comps.get(1, (0, 1))[0]
        if set(comps) - {0, 1}:
            print(f'   q={q}: UNEXPECTED Jordan scales {sorted(comps)}'); ok = False
        if a_true != a_pred:
            print(f'   q={q}: alpha from lattice {a_true} != rank formula {a_pred}'); ok = False
        alphas.append(a_true)
        for theta in (r, s):
            cq = {j: (nn, eps) for j, nn, eps in data[theta][2].get(q, [])}
            deltas[theta].append(1 if cq.get(1, (0, 1))[1] == -1 else 0)
        rel = (deltas[s][-1] + a_true * (q % 4 == 3)) % 2 == deltas[r][-1]
        print(f'   q={q}: alpha = {a_true} (rank formula {a_pred}), delta_r = {deltas[r][-1]}, delta_s = {deltas[s][-1]}, '
              f'delta relation {"OK" if rel else "FAILS"}')
        ok &= rel
    for theta, cyc in ((r, res['cyc_r']), (s, res['cyc_s'])):
        for q, e, eps in cyc:
            cq = {j: (nn, ee) for j, nn, ee in data[theta][2].get(q, [])}
            good = cq.get(e) == (1, eps) and set(cq) == {0, e}
            print(f'   cyclic prime {q} on Lambda_{theta}: predicted (e={e}, eps={eps}), lattice {data[theta][2].get(q)} -> {"OK" if good else "MISMATCH"}')
            ok &= good
    mm = res['mm']
    side_s = gs.side_ok(g, qs, tuple(alphas), tuple(deltas[s]), res['cyc_s'], res['two_unimod'][1], mm)
    side_r = gs.side_ok(f, qs, tuple(alphas), tuple(deltas[r]), res['cyc_r'], res['two_unimod'][0], mm)
    allowed = tuple(alphas) in res['allowed'] or any(
        all((a == b) or (0 < a < mm and 0 < b < mm and a % 8 == b % 8) for a, b in zip(tuple(alphas), t)) for t in res['allowed'])
    print(f'   actual (alpha, delta) passes side_ok: r {side_r}, s {side_s}; alpha in scan-allowed set: {allowed}')
    ok &= side_r and side_s and allowed
    print('   ->', 'CONSISTENT' if ok else '*** INCONSISTENT ***')
    return ok


if __name__ == '__main__':
    allok = True
    allok &= check('Petersen', petersen())
    for n in (5, 7, 9, 13, 15, 17, 19, 21):
        allok &= check(f'T({n})', triangular(n))
    print('\nALL CONSISTENT' if allok else '\n*** SOME INCONSISTENCY ***')
