"""Check of C1's pattern claim: whenever dim ker Q_c = 4 (block-line rank 14), eta_c (= iota here) has NO digon with mate
(iota(v) = mate(v) never happens), while tau has >= 2 digons (tau(v) = mate(v))."""
from collections import Counter
from blockline_exhaustive import labels, mate, rank_mod, matrix, fpf_involutions, partitions, tau_of_type
allinv = list(fpf_involutions(labels))
pat = Counter()
for parts in partitions(6):
    tau = tau_of_type(parts)
    for iota in allinv:
        if rank_mod(matrix(tau, iota)) == 14:
            dt = sum(1 for v in labels if tau[v] == mate[v]) // 2
            di = sum(1 for v in labels if iota[v] == mate[v]) // 2
            both = sum(1 for v in labels if tau[v] == mate[v] and iota[v] == mate[v]) // 2
            pat[(dt, di, both)] += 1
print('rank-14 cases by (#tau digons, #iota digons, #common):', dict(pat))
