"""Step 1' : near-design (= (i)+(iii)) together with the BLOCK-LEVEL counts (a relaxation of (ii)), as one CP-SAT model.
Variables: p[i,j] = the non-mu pairs i, j are iota-partners (Phi(z), Phi(sigma z) for one sigma-pair);
           a[i,j,k] = that iota-pair is the image of a sigma-pair of block pair k (15 block pairs, two sigma-pairs each).
Constraints: iota perfect matching with the K3 cross counts; sum_k a[i,j,k] = p[i,j]; each block pair gets exactly 2 iota-pairs;
  (ii-block) for each block r and label b:  sum_{k containing r} sum_{(i,j)} a[i,j,k] [b in e_i u e_j] = K2[r+,b] + K2[r-,b].
Infeasible  =>  (i)+(ii)+(iii) infeasible for this block type  =>  no block of this type in an involution.
Works for any block type (mu, tau) (C2 notation); cfg {"cycles": ...} for a full block, or {"mu": [...], "tau": [...]}.
Usage: python step1prime.py cfg.json time_limit workers"""
import itertools, json, sys, time
from ortools.sat.python import cp_model
import iota_design as D

L = list(range(12))


def build(mu, tau):
    K2, K3 = D.K2K3(mu, tau)
    edges = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    E = len(edges)
    BP = list(itertools.combinations(range(6), 2))
    md = cp_model.CpModel()
    P = {}
    for i in range(E):
        for j in range(i + 1, E):
            common = set(edges[i]) & set(edges[j])
            if any(K3[b, b] < 2 for b in common):
                continue
            P[i, j] = md.NewBoolVar('p%d_%d' % (i, j))
    for i in range(E):
        md.AddExactlyOne(P[min(i, j), max(i, j)] for j in range(E) if j != i and (min(i, j), max(i, j)) in P)
    acc = {}
    for (i, j), v in P.items():
        for b in edges[i]:
            for d in edges[j]:
                key = (min(b, d), max(b, d))
                acc.setdefault(key, []).append(v)
    for b in L:
        for d in L:
            if d < b:
                continue
            target = int(K3[b, d]) if b != d else int(K3[b, b]) // 2
            md.Add(sum(acc.get((b, d), [])) == target)
    A = {}
    for (i, j), v in P.items():
        lits = []
        for k in range(15):
            A[i, j, k] = md.NewBoolVar('a%d_%d_%d' % (i, j, k))
            lits.append(A[i, j, k])
        md.Add(sum(lits) == v)
    for k in range(15):
        md.Add(sum(A[i, j, k] for (i, j) in P) == 2)
    for r in range(6):
        ks = [k for k, bp in enumerate(BP) if r in bp]
        for b in L:
            terms = [A[i, j, k] for (i, j) in P if b in edges[i] or b in edges[j] for k in ks]
            md.Add(sum(terms) == int(K2[2 * r, b] + K2[2 * r + 1, b]))
    return md, P, A, edges


if __name__ == '__main__':
    cfg = json.load(open(sys.argv[1]))
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 120
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    if 'cycles' in cfg:
        from fullblock_rigid import mates_from_D
        m = mates_from_D([tuple(c) for c in cfg['cycles']])
        mu = tau = m
    else:
        mu, tau = cfg['mu'], cfg['tau']
    t0 = time.time()
    md, P, A, edges = build(mu, tau)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl
    s.parameters.num_workers = workers
    st = s.Solve(md)
    print(sys.argv[1], "step 1' (near-design + block-level (ii)):", s.StatusName(st), '%.1fs' % (time.time() - t0),
          '| vars', len(P) * 16, flush=True)
