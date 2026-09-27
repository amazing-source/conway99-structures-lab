"""O_24, class 1 of C1's list: an explicit ONE-VECTOR certificate (energy trap), checked exactly (structures 27/09).
Scaled coordinates (a in Z^24, norm a.a/8).  g = s + 2 e_10, s = +-1 on {5,6,7,8,9,12,15,17}, so w = g/sqrt8 has |w|^2 = 12/8.
Claim: every norm-3 vector u of O_24 orthogonal to y1, y2 has |E.g| >= 4, i.e. (u.w)^2 >= 1/4.  The frame identity at w gives
sum_u (u.w)^2 + sum_h (h.w)^2 + sum_y (y.w)^2 = 7 |w|^2 = 21/2; the 42 u already give >= 42/4 = 21/2, hence y.w = 0 for both y.
But (B1.g, B2.g) != (0, 0).  Contradiction: class 1 has no frame."""
import json, sys
import numpy as np
from o24_ustage import C1_YPAIRS, golay_from_basis, sign_vectors
D = json.load(open(C1_YPAIRS)); code = golay_from_basis(D['golay_basis']); SV = sign_vectors(code)
o = D['orbits'][1]
B1, B2 = np.array(o['a1']), np.array(o['a2'])
g = np.zeros(24, dtype=np.int64)
for i, sg in zip([5, 6, 7, 8, 9, 12, 15, 17], [1, 1, -1, -1, -1, 1, -1, -1]):
    g[i] = sg
g[10] = 2
U = SV[(SV @ B1 == 0) & (SV @ B2 == 0)]
m = np.abs(U @ g).min()
print('u-candidates:', len(U), ' min |E.g| =', int(m), ' (B1.g, B2.g) =', (int(B1 @ g), int(B2 @ g)), ' |g|^2 =', int(g @ g))
ok = m >= 4 and (B1 @ g != 0 or B2 @ g != 0) and 42 * 16 >= 56 * int(g @ g)
print('certificate VALID' if ok else 'certificate INVALID')
