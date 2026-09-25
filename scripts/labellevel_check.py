"""Independent label-level check of the local model of doubles_budget3.local_min.

Model x = {(0,0),(1,0)} (block i = 0, block j = 1), D(X) in pair {0,2} (k = 2): v = {(0,0),(2,0)}, sigma v = {(0,1),(2,1)}.
Everything is derived from lambda = 1, mu = 2 at the level of LABELS (no count-level shortcuts):
  * x has exactly 1 labelled neighbour containing each of (0,0), (0,1), (1,0), (1,1) (lambda = 1 with its own elements,
    mu = 2 with their mates, one common neighbour being the element itself), and exactly 2 containing each element of
    blocks 2..6 (mu = 2);
  * the neighbours meeting block 0 are v, sigma v (orbit paired via i); the neighbour containing (1,0) and the one
    containing (1,1) lie in the two pairs {1,l}, {1,r} given by the free 6-cycle at j (both assignments, all elements);
  * no two sigma-conjugate labels of an orbit other than D(X) are both adjacent to x (else that orbit would be D(X));
  * complete D-triangles tau: |N(x) ∩ Lab(tau)| = 2 / 1 / 3 (x's pair in tau / meets tau in one block / disjoint).
doubles(x) = number of pairs other than {0,2} carrying exactly 2 adjacent labels.  Returns the minimum, or None."""
import itertools

BLOCKS = range(7)


def lab(p, a, q, b):
    return ((p, a), (q, b)) if p < q else ((q, b), (p, a))


def conj(L):
    (p, a), (q, b) = L
    return ((p, 1 - a), (q, 1 - b))


def pair_of(L):
    return (L[0][0], L[1][0])


def solve(tris, nbrs):
    l, r = nbrs
    best = None
    v = lab(0, 0, 2, 0)
    sv = lab(0, 1, 2, 1)
    for swap in (0, 1):
        for e1 in (0, 1):
            for e2 in (0, 1):
                n1 = lab(1, 0, l, e1) if swap == 0 else lab(1, 1, l, e1)
                n2 = lab(1, 1, r, e2) if swap == 0 else lab(1, 0, r, e2)
                forced = [v, sv, n1, n2]
                need = {}
                for b in BLOCKS:
                    for e in (0, 1):
                        need[(b, e)] = 1 if b in (0, 1) else 2
                okf = True
                for L in forced:
                    for el in L:
                        need[el] -= 1
                        if need[el] < 0:
                            okf = False
                if not okf or any(need[(b, e)] != 0 for b in (0, 1) for e in (0, 1)):
                    continue
                cands = [lab(p, a, q, c) for p, q in itertools.combinations(range(2, 7), 2)
                         for a in (0, 1) for c in (0, 1)]
                cands = [L for L in cands if L not in forced]
                elements = [(b, e) for b in range(2, 7) for e in (0, 1)]
                chosen = []
                chosen_labels = []

                def evaluate(_):
                    nonlocal best
                    nb = set(forced) | set(chosen_labels)
                    if len(nb) != 12:
                        return
                    for L in nb:
                        if L not in (v, sv) and conj(L) in nb:
                            return
                    for t in tris:
                        s = set(t)
                        inter = len(s & {0, 1})
                        want = 2 if inter == 2 else (1 if inter == 1 else 3)
                        cnt = sum(1 for L in nb if set(pair_of(L)) <= s)
                        if cnt != want:
                            return
                    per_pair = {}
                    for L in nb:
                        per_pair[pair_of(L)] = per_pair.get(pair_of(L), 0) + 1
                    d = sum(1 for P, c in per_pair.items() if c == 2 and P != (0, 2))
                    if best is None or d < best:
                        best = d

                def dfs2():
                    el = next((x for x in elements if need[x] > 0), None)
                    if el is None:
                        evaluate(None)
                        return
                    last_idx = chosen[-1][1] if chosen and chosen[-1][0] == el else -1
                    for idx in range(last_idx + 1, len(cands)):
                        L = cands[idx]
                        if el not in L or L in chosen_labels:
                            continue
                        other = L[1] if L[0] == el else L[0]
                        if need[other] <= 0:
                            continue
                        need[el] -= 1
                        need[other] -= 1
                        chosen.append((el, idx))
                        chosen_labels.append(L)
                        dfs2()
                        chosen_labels.pop()
                        chosen.pop()
                        need[el] += 1
                        need[other] += 1

                dfs2()
    return best


if __name__ == '__main__':
    import sys
    import time
    sys.path.insert(0, 'scripts')
    import doubles_budget3 as d3
    fano = tuple(sorted([(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]))
    t0 = time.time()
    for key in [((), (3, 4)), ((), (2, 3)), (fano, (3, 4)), (fano, (3, 5)), (fano, (2, 3))]:
        print(key[1], 'tris' if key[0] else 'no tris', 'label-level:', solve(*key), ' count-level:', d3.local_min(key))
    print(f'{time.time() - t0:.1f} s')
