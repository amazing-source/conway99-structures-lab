"""SYMMETRIC local solutions of the free star: Phi equivariant under a subgroup H of Stab(mu, tau) (signed block permutations).
sigma (flip every sign) is central in the signed permutation group, so H acts on the 30 sigma-pairs {z, sigma z}; an
H-equivariant Phi is determined by the ordered pair (Phi(z), Phi(sigma z)) on one representative per H-orbit of sigma-pairs,
invariant under the stabilizer of that representative.  (i), (ii), (iii) expand linearly.  Any solution is re-checked
explicitly.  Absence of symmetric solutions proves nothing about general ones."""
import itertools, json, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp
from full_local import K2, K3, sig, L


def signed_perms():
    for perm in itertools.permutations(range(6)):
        for signs in itertools.product((0, 1), repeat=6):
            yield {i: 2 * perm[i // 2] + ((i % 2) ^ signs[i // 2]) for i in L}


def stabilizer(mu, tau):
    return [g for g in signed_perms()
            if all(g[mu[i]] == mu[g[i]] for i in L) and all(g[tau[i]] == tau[g[i]] for i in L)]


def key(g):
    return tuple(g[i] for i in L)


def gen_subgroup(gens):
    ident = {i: i for i in L}
    H = {key(ident): ident}
    frontier = [ident]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                c = {i: g[a[i]] for i in L}
                if key(c) not in H:
                    H[key(c)] = c
                    new.append(c)
        frontier = new
    return list(H.values())


def act(g, e):
    return tuple(sorted((g[e[0]], g[e[1]])))


def sflip(z):
    return tuple(sorted((sig[z[0]], sig[z[1]])))


def check(mu, tau, Phi):
    k2, k3 = K2(mu, tau), K3(mu, tau)
    far = [e for e in itertools.combinations(L, 2) if sig[e[0]] != e[1]]
    if sorted(Phi) != sorted(far) or len(set(Phi.values())) != 60 or any(mu[e[0]] == e[1] for e in Phi.values()):
        return False
    M2 = np.zeros((12, 12), int)
    M3 = np.zeros((12, 12), int)
    for z in far:
        for a in z:
            for b in Phi[z]:
                M2[a, b] += 1
        for b in Phi[z]:
            for dd in Phi[sflip(z)]:
                M3[b, dd] += 1
    return bool((M2 == k2).all() and (M3 == k3).all())


def solve_equivariant(mu, tau, H, tl=120):
    k2, k3 = K2(mu, tau), K3(mu, tau)
    far = [e for e in itertools.combinations(L, 2) if sig[e[0]] != e[1]]
    nonmu = [e for e in itertools.combinations(L, 2) if mu[e[0]] != e[1]]
    spairs = []
    seen = set()
    for z in far:
        if z in seen:
            continue
        seen |= {z, sflip(z)}
        spairs.append((z, sflip(z)))
    # H-orbits of sigma-pairs (as unordered pairs; H commutes with sigma)
    def canon(p):
        return tuple(sorted(p))
    orbit_rep = {}
    reps = []
    for z, sz in spairs:
        c = canon((z, sz))
        if c in orbit_rep:
            continue
        for g in H:
            orbit_rep[canon((act(g, z), act(g, sz)))] = len(reps)
        reps.append((z, sz))
    var = []
    images = []
    for ri_, (z, sz) in enumerate(reps):
        stab = [g for g in H if act(g, z) == z]          # g fixing z also fixes sigma z
        swap = [g for g in H if act(g, z) == sz]          # g exchanging z and sigma z
        for e in nonmu:
            if not all(k2[a, b] > 0 for a in z for b in e):
                continue
            if not all(act(g, e) == e for g in stab):
                continue
            for e2 in nonmu:
                if e2 == e or not all(k2[a, b] > 0 for a in sz for b in e2):
                    continue
                if not all(act(g, e2) == e2 for g in stab):
                    continue
                if not all(act(g, e) == e2 and act(g, e2) == e for g in swap):
                    continue
                img = {}
                ok = True
                for g in H:
                    for w, f in ((act(g, z), act(g, e)), (act(g, sz), act(g, e2))):
                        if w in img and img[w] != f:
                            ok = False
                        img[w] = f
                if not ok or len(set(img.values())) != len(img):
                    continue
                var.append((ri_, e, e2))
                images.append(img)
    nv = len(var)
    if nv == 0:
        return 'INFEASIBLE (no admissible variable)', None, 0
    ri_l, ci_l, vv, vals = [], [], [], []
    r = 0

    def row(entries, val):
        nonlocal r
        for j, c in entries:
            if c:
                ri_l.append(r); ci_l.append(j); vv.append(c)
        vals.append(val); r += 1
    for ri_ in range(len(reps)):
        row([(j, 1) for j, v in enumerate(var) if v[0] == ri_], 1)
    for e0 in nonmu:
        row([(j, sum(1 for f in images[j].values() if f == e0)) for j in range(nv)], 1)
    for a in L:
        for b in L:
            row([(j, sum(1 for w, f in images[j].items() if a in w and b in f)) for j in range(nv)], int(k2[a, b]))
    for b in L:
        for dd in L:
            row([(j, sum(1 for w, f in images[j].items() if b in f and dd in images[j][sflip(w)])) for j in range(nv)], int(k3[b, dd]))
    A = sp.csr_matrix((np.array(vv, float), (np.array(ri_l), np.array(ci_l))), shape=(r, nv))
    v = np.array(vals, float)
    res = milp(np.zeros(nv), constraints=[LinearConstraint(A, v, v)], integrality=np.ones(nv), bounds=Bounds(0, 1),
               options={'time_limit': tl})
    if res.status == 0:
        Phi = {}
        for j in np.nonzero(res.x > 0.5)[0]:
            Phi.update(images[j])
        return ('FEASIBLE' if check(mu, tau, Phi) else 'FEASIBLE-but-check-failed'), Phi, nv
    return {2: 'INFEASIBLE'}.get(res.status, 'unknown(%d)' % res.status), None, nv


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\RESEARCHER\Desktop\c1_experimental\q7\A33_torsion_reps.json"
    d = json.load(open(src))
    if isinstance(d['mu'], list):
        mu = {i: d['mu'][i] for i in L}; tau = {i: d['tau'][i] for i in L}
    else:
        mu = {int(k): v for k, v in d['mu'].items()}; tau = {int(k): v for k, v in d['tau'].items()}
    G = stabilizer(mu, tau)
    print('|Stab(mu,tau)| =', len(G), flush=True)
    # try cyclic subgroups generated by each element (deduplicated), largest first, then the full group
    subs = {}
    for g in G:
        H = gen_subgroup([g])
        subs.setdefault(frozenset(key(h) for h in H), H)
    cand = sorted(subs.values(), key=len, reverse=True)
    tries = [G] + cand
    found = 0
    t0 = time.time()
    for H in tries:
        if len(H) == 1:
            continue
        st, Phi, nv = solve_equivariant(mu, tau, H, tl=60)
        print('|H| =', len(H), 'vars', nv, '->', st, '%.0fs' % (time.time() - t0), flush=True)
        if Phi and st == 'FEASIBLE':
            found += 1
            out = {'labels': 'i in 0..11: block i//2, sign + if i even; sigma(i) = i^1', 'mu': mu, 'tau': tau,
                   'H_order': len(H), 'Phi (far root -> non-mu pair)': {str(list(k)): list(v) for k, v in Phi.items()},
                   'satisfies': '(i)+(ii)+(iii), checked'}
            json.dump(out, open('symmetric_local_solution_%d.json' % found, 'w'), indent=1)
            print('  saved symmetric_local_solution_%d.json' % found, flush=True)
            if found >= 2:
                break
