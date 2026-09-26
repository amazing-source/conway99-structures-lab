"""Block lines at a block c with NO one-block square through c (every orbit through c free at c: e.g. all pairs at c dd):
tau_c and iota_c share no edge (a common edge = a digon of R_c = a D-pair via c).  Histogram of (rank_7, local sign) over all
tau classes (with class sizes as weights: |orbit of tau under the centralizer|) and all such iota."""
from collections import Counter
from math import factorial
from blockline_exhaustive import labels, mate, rank_mod, matrix, fpf_involutions, partitions, tau_of_type
from blockline_Ablock import det_mod, chi
allinv = list(fpf_involutions(labels))
tot = Counter()
for parts in partitions(6):
    tau = tau_of_type(parts)
    h = Counter()
    for iota in allinv:
        if any(iota[v] == tau[v] for v in labels): continue
        M = matrix(tau, iota); r = rank_mod(M)
        s = chi(det_mod([row[1:] for row in M[1:]])) if r == 18 else None
        h[(r, s)] += 1
    tot.update(h)
    print(parts, dict(sorted(h.items(), key=str)), flush=True)
print('ALL:', dict(sorted(tot.items(), key=str)))
