"""At g- = 18 with a rank-18 block line E, the mod-7 form B = S - 3I is isometric (mod radicals) to B|_E, so
eps- = chi((-1)^g disc B) = chi(disc B|_E) must equal -1 (genus: [eps- = -1] = M - g/2 = 1 mod 2, M even).
Compute chi(det of the 18x18 principal submatrix without d_c) over random local configurations of rank 18."""
import random, sys
from collections import Counter
sys.path.insert(0, '.')
from blockline_rank import matrix, rand_involution, rank_mod, labels, mate
P = 7
def det_mod(M, p=P):
    M = [[x % p for x in r] for r in M]; n = len(M); d = 1
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c]), None)
        if piv is None: return 0
        if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
        d = d * M[c][c] % p; inv = pow(M[c][c], -1, p)
        for i in range(c + 1, n):
            if M[i][c]:
                f = M[i][c] * inv % p; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[c])]
    return d % p
chi = lambda a: 0 if a % P == 0 else (1 if pow(a, (P - 1) // 2, P) == 1 else -1)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
for label, fpf in [('iota fixed-point-free', True), ('iota with fixed points allowed', False)]:
    cnt = Counter()
    for _ in range(N):
        tau = rand_involution(rng, True); iota = rand_involution(rng, fpf)
        M = matrix(tau, iota)
        if rank_mod(M) != 18: continue
        sub = [row[1:] for row in M[1:]]            # drop d_c (the radical z_c has d_c-coefficient 3)
        cnt[chi(det_mod(sub))] += 1
    print(label, '-> chi(disc B|_E) over rank-18 configurations:', dict(cnt), flush=True)
print('tau = iota = mate (untwisted):', chi(det_mod([r[1:] for r in matrix(mate, mate)[1:]])))
