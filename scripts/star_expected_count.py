"""Heuristic expected number of full-star solutions (independence heuristic): sum over patterns (apex set + perfect
matching of the rest) of prod_v #cand(v, role) * prod_pairs p_vw(value), with p_vw(value) the empirical fraction of
candidate pairs (right roles) having g = value and overlap ok.  Not a proof; a calibration of overdetermination."""
import sys, itertools, math
sys.path.insert(0, 'scripts')
from star_shadows import build, ball_data, green_mod7, N, part, P7
from verify_astra_d5_cyclic import ball as cyc_ball

def expected(name, adj, u=0):
    Bm, H, tau, rk = ball_data(adj)
    Lg = green_mod7(H)
    B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
    K = [[12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) for q in range(N)] for p in range(N)]
    partner = {p: next(q for q in adj[p] if part(q) == part(p)) for p in range(N)}
    uB = next(q for q in adj[u] if part(q) == 1); uC = next(q for q in adj[u] if part(q) == 2)
    prs = lambda t: [(12 * t + i, 12 * t + j) for i in range(12) for j in range(i + 1, 12) if partner[12 * t + i] != 12 * t + j]
    groups = [v for v in range(12) if v != u and v != partner[u]]
    cand = {}
    for v in groups:
        out = []
        for eb in prs(1):
            for ec in prs(2):
                S = (u, v) + eb + ec; Ss = set(S)
                if any(len(adj[p] & Ss) > 2 - (p in Ss) for p in range(N)): continue
                if any(K[p][q] == 0 for p, q in itertools.combinations(S, 2)): continue
                if sum(Lg[p][q] for p in S for q in S) % P7 != 6: continue
                out.append(S)
        cand[v] = out
    role = lambda S: 'apex' if (uB in S) != (uC in S) else ('bad' if (uB in S and uC in S) else 'm')
    byrole = {v: {'apex': [S for S in cand[v] if role(S) == 'apex'], 'm': [S for S in cand[v] if role(S) == 'm']} for v in groups}
    hv = {v: {r: [[sum(Lg[p][q] for p in S) % P7 for q in range(N)] for S in byrole[v][r]] for r in ('apex', 'm')} for v in groups}
    # empirical pair probabilities p[v,w][(rv,rw,val)]
    prob = {}
    for v, w in itertools.combinations(groups, 2):
        for rv in ('apex', 'm'):
            for rw in ('apex', 'm'):
                A, B_ = byrole[v][rv], byrole[w][rw]
                if not A or not B_:
                    continue
                cnt = {2: 0, 3: 0}; tot = len(A) * len(B_)
                for i, S in enumerate(A):
                    hS = hv[v][rv][i]; Ss = set(S)
                    for T in B_:
                        gg = sum(hS[q] for q in T) % P7
                        if gg in (2, 3) and len(Ss & set(T)) <= 2 - (3 - gg):
                            cnt[gg] += 1
                prob[v, w, rv, rw] = (cnt[2] / tot, cnt[3] / tot)
    # sum over patterns
    total = 0.0
    def pms(rest):
        if not rest:
            yield []
            return
        v = rest[0]
        for w in rest[1:]:
            for m in pms([x for x in rest if x not in (v, w)]):
                yield [(v, w)] + m
    for k in range(0, 11, 2):
        for ap in itertools.combinations(groups, k):
            rest = [v for v in groups if v not in ap]
            roles = {v: ('apex' if v in ap else 'm') for v in groups}
            base = 1.0
            ok = True
            for v in groups:
                n = len(byrole[v][roles[v]])
                if n == 0:
                    ok = False; break
                base *= n
            if not ok:
                continue
            for m in pms(rest):
                matched = {frozenset(e) for e in m}
                val = base
                for v, w in itertools.combinations(groups, 2):
                    key = (v, w, roles[v], roles[w])
                    if key not in prob:
                        val = 0; break
                    p2, p3 = prob[key]
                    val *= p2 if frozenset((v, w)) in matched else p3
                total += val
    print(f'{name}: candidates per group {[len(cand[v]) for v in groups]}; heuristic expected number of full-star solutions = {total:.3e}', flush=True)

a = [9,4,7,8,1,6,5,2,3,0,11,10]; b = [7,2,1,9,11,8,10,0,5,3,6,4]; c = [7,3,10,1,11,9,8,0,6,5,2,4]
expected('Astra', build(a, b, c, [(i + 1) % 12 for i in range(12)]))
expected('cyclic241', cyc_ball(2, 4, 1))
