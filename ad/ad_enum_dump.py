"""All AD configurations (every pair anchored or doubly-disjoint 'dd'; Sigma = empty; no half pairs) at the orientation level,
filtered by proved lemmas only (note 41):
  - the dd pairs form a graph Delta on the 7 blocks admitting a 'dis graph' (2-regular multigraph on Delta's edges, edges between
    DISJOINT dd pairs);
  - degree lemmas (39.4 / 26.2): d_k not in {1,5}; d_k <= 3 if out_k >= 1; d_k = 4 only if out_k = 0 and dis_k = 2; d_k = 6 only
    if dis_k = 0 (then out_k = 0: an A-block); blocks with dis_k = 0: d_k in {0,2,3,6};
  - tail lemma 24.1 on the one-block squares forced by the orientation (digon at d = 2, triangle at d = 3, every A-arm has a
    square): component of the pivot inside F_k has >= 3 + [j in it] blocks => dead if |F_k| <= 2 or (|F_k| = 3 and j in F_k).
    (At d = 4 the D-graph on 4 pairs is a 4-cycle or two digons: every arm has two square-neighbours; the check is applied to
    the pairs that are forced in BOTH types only if they coincide, otherwise to the weaker 'every arm' form.)
Counts survivors up to S_7.  Usage: python ad_enum_all.py [max n_dd]"""
import itertools, sys, json, time
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
N = 7
ALLP = [(a, b) for a in range(N) for b in range(a + 1, N)]

def admits(dd):
    G = nx.Graph(); G.add_nodes_from(dd)
    for p, q in itertools.combinations(dd, 2):
        if not (set(p) & set(q)): G.add_edge(p, q)
    B = nx.Graph()
    for p in dd: B.add_node(('L', p)); B.add_node(('R', p))
    for p, q in G.edges(): B.add_edge(('L', p), ('R', q)); B.add_edge(('L', q), ('R', p))
    if not dd: return True
    M = nx.bipartite.maximum_matching(B, top_nodes=[('L', p) for p in dd])
    return len(M) // 2 == len(dd)

def delta_classes(m):
    seen = []
    for E in itertools.combinations(ALLP, m):
        G = nx.Graph(); G.add_nodes_from(range(N)); G.add_edges_from(E)
        key = tuple(sorted(d for _, d in G.degree()))
        if any(k == key and nx.is_isomorphic(G, H) for k, H, _ in seen): continue
        seen.append((key, G, E))
    return [E for _, _, E in seen if admits(list(E))]

def tail_dead(In, F, d):
    for i in range(N):
        arms = In[i]
        if len(arms) in (2, 3):
            pairs = itertools.permutations(arms, 2)
        elif len(arms) == 6:
            pairs = [(None, k) for k in arms]
        elif len(arms) == 4:
            pairs = [(None, k) for k in arms]
        else: continue
        for j, k in pairs:
            if len(F[k]) <= 2: return True
            if j is not None and len(F[k]) == 3 and j in F[k]: return True
    return False

DUMP = []
def run(maxdd):
    t0 = time.time(); summary = {}
    for m in range(2, maxdd + 1):
        classes = delta_classes(m); alive_total = 0; orient_total = 0
        for E in classes:
            dd = set(E); dis = {k: sum(1 for p in dd if k in p) for k in range(N)}
            rest = [p for p in ALLP if p not in dd]
            G = nx.Graph(); G.add_nodes_from(range(N)); G.add_edges_from(E)
            stab = [g for g in (dict(x) for x in GraphMatcher(G, G).isomorphisms_iter())]
            reps = set()
            # DFS over the anchored pairs with degree pruning
            d = {k: 0 for k in range(N)}; out = {k: 0 for k in range(N)}; heads = {}
            rem_at = {k: sum(1 for p in rest if k in p) for k in range(N)}
            def ok_final():
                for k in range(N):
                    if d[k] in (1, 5): return False
                    if out[k] >= 1 and d[k] > 3 and not (d[k] == 4 and out[k] == 0): return False
                    if d[k] == 4 and not (out[k] == 0 and dis[k] == 2): return False
                    if d[k] == 6 and dis[k] != 0: return False
                    if dis[k] == 0 and d[k] not in (0, 2, 3, 6): return False
                return True
            def prune(k):
                # partial: out >= 1 forbids d > 3 except the d=4/out=0 case (then out must stay 0)
                if out[k] >= 1 and d[k] > 3: return True
                return False
            def dfs(t):
                if t == len(rest):
                    if not ok_final(): return
                    In = {c: [x for x in range(N) if x != c and heads.get((min(c, x), max(c, x))) == c] for c in range(N)}
                    F = {c: {x for x in range(N) if x != c and x not in In[c]} for c in range(N)}
                    if tail_dead(In, F, d): return
                    arcs = tuple(sorted((p[0] + p[1] - h, h) for p, h in heads.items()))
                    key = min(tuple(sorted((g[a], g[b]) for a, b in arcs)) for g in stab)
                    reps.add(key); return
                a, b = rest[t]
                for h in (a, b):
                    o = a + b - h
                    d[h] += 1; out[o] += 1; heads[(a, b)] = h
                    if not (prune(h) or prune(o)): dfs(t + 1)
                    d[h] -= 1; out[o] -= 1; del heads[(a, b)]
            dfs(0)
            alive_total += len(reps); DUMP.extend({"dd": [list(p) for p in E], "arcs": [list(a) for a in r]} for r in reps)
        summary[m] = (len(classes), alive_total)
        print("n_dd = %d: Delta classes %d, AD orientation classes surviving degree + tail lemmas: %d  (%.0fs)" % (m, len(classes), alive_total, time.time() - t0)); sys.stdout.flush()
    return summary

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 6); json.dump(DUMP, open(sys.argv[2], "w"))
