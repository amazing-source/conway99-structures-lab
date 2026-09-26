"""Matchability structure of the star (all constraints): which groups can be apexes, which pairs of groups can be
matched (candidates with g = 2, both non-apex, overlap ok), and does a choice 'apexes + perfect matching of the rest'
exist at this pairwise level?  If not, the obstruction is a plain matching obstruction."""
import sys, itertools
sys.path.insert(0, 'scripts')
import star_shadows as ss
from star_shadows import build, ball_data, green_mod7, N, part, P7
from verify_astra_d5_cyclic import ball as cyc_ball

def analyse(name, adj, u=0):
    Bm, H, tau, rk = ball_data(adj)
    Lg = green_mod7(H)
    B2 = [[sum(Bm[p][r] * Bm[r][q] for r in range(N)) for q in range(N)] for p in range(N)]
    K = [[12 * (p == q) - Bm[p][q] + 2 - B2[p][q] - (part(p) == part(q)) for q in range(N)] for p in range(N)]
    partner = {p: next(q for q in adj[p] if part(q) == part(p)) for p in range(N)}
    uB = next(q for q in adj[u] if part(q) == 1); uC = next(q for q in adj[u] if part(q) == 2)
    pairs = lambda t: [(12 * t + i, 12 * t + j) for i in range(12) for j in range(i + 1, 12) if partner[12 * t + i] != 12 * t + j]
    groups = [v for v in range(12) if v != u and v != partner[u]]
    cand = {}
    for v in groups:
        out = []
        for eb in pairs(1):
            for ec in pairs(2):
                S = (u, v) + eb + ec; Ss = set(S)
                if any(len(adj[p] & Ss) > 2 - (p in Ss) for p in range(N)): continue
                if any(K[p][q] == 0 for p, q in itertools.combinations(S, 2)): continue
                if sum(Lg[p][q] for p in S for q in S) % P7 != 6: continue
                out.append(S)
        cand[v] = out
    apex_ok = {v: any((uB in S) != (uC in S) for S in cand[v]) for v in groups}
    nonapex = {v: [S for S in cand[v] if uB not in S and uC not in S] for v in groups}
    g = lambda S, T: sum(Lg[p][q] for p in S for q in T) % P7
    match_ok = {}
    for v, w in itertools.combinations(groups, 2):
        ok = False
        for S in nonapex[v]:
            for T in nonapex[w]:
                if g(S, T) == 2 and len(set(S) & set(T)) <= 1:
                    ok = True; break
            if ok: break
        match_ok[v, w] = ok
    # does some choice of apex set (even size, apex-capable groups) + perfect matching of the rest exist?
    def has_pm(rest):
        if not rest: return True
        v = rest[0]
        return any(match_ok[min(v, w), max(v, w)] and has_pm([x for x in rest if x not in (v, w)]) for w in rest[1:])
    choices = []
    for k in range(0, 11, 2):
        for ap in itertools.combinations([v for v in groups if apex_ok[v]], k):
            rest = [v for v in groups if v not in ap]
            if has_pm(rest):
                choices.append(ap)
    mg = [(v, w) for (v, w), ok in match_ok.items() if ok]
    print(f'{name}: apex-capable groups {[v for v in groups if apex_ok[v]]}; matchable pairs {len(mg)} of 45; '
          f'(apex set, perfect matching) choices at the pairwise level: {len(choices)}; e.g. {choices[:4]}', flush=True)

a = [9,4,7,8,1,6,5,2,3,0,11,10]; b = [7,2,1,9,11,8,10,0,5,3,6,4]; c = [7,3,10,1,11,9,8,0,6,5,2,4]
analyse('Astra', build(a, b, c, [(i + 1) % 12 for i in range(12)]))
analyse('cyclic241', cyc_ball(2, 4, 1))
