"""id 1 (a1^8), step 2 (structures 27/09): for each of the 14 y1 classes (state/u1_y1reps.json), all admissible y2 with
y1.y2 = 0 and the exact pair filters:
  roots (alpha = 4 e_j, j in K): (a1_j/2)^2 + (a2_j/2)^2 <= 5, i.e. not |a1_j| = |a2_j| = 4;
  shadow roots beta (a-norm 16): (beta.y1)^2 + (beta.y2)^2 <= 5/2, i.e. (a_b.a_1)^2 + (a_b.a_2)^2 <= 160;
  norm-3 u: never |a_u.a_1| = |a_u.a_2| = 24.
y2 ranges over all sign patterns of the sign-normalised admissible y (32 inner codewords on K^c x signs on K; signs on
K-zeros are not doubled).  Memory-lean: vectorised by chunks, y2 stored as int8.
Writes state/u1_pairs_KK.npy and prints counts.  Usage: python golay_y2.py k [k ...]"""
import itertools, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from golay_family import WORDS, OM
from golay_y import lattice_data

T0 = time.time()


def main():
    K, A, N3, BETA = lattice_data()
    Kc = OM ^ K
    Kl = [i for i in range(24) if K >> i & 1]
    Kcl = [i for i in range(24) if Kc >> i & 1]
    inner = [w & Kc for w in WORDS if w & K == 0]
    Wsign = np.array([[(-1) ** (w >> j & 1) for j in Kcl] for w in inner], dtype=np.int32)    # 32 x 16
    Yn = np.load(os.path.join(HERE, 'state', 'u1_ynorm.npy')).astype(np.int32)
    reps = json.load(open(os.path.join(HERE, 'state', 'u1_y1reps.json')))['reps']
    EPS = np.array(list(itertools.product((1, -1), repeat=8)), dtype=np.int32)                # 256 x 8
    B32, N32 = BETA.T.astype(np.int32), N3.T.astype(np.int32)
    nzK = (Yn[:, Kl] != 0)                                                                      # n x 8
    for k in [int(a) for a in sys.argv[1:] if not a.startswith('--')]:
        y1 = np.array(reps[k], dtype=np.int32)
        b1 = (BETA @ y1).astype(np.int64)
        n1 = np.abs(N3 @ y1)
        U1 = N3[n1 == 24].T.astype(np.int32)                       # 24 x (#u with |u.y1| = 3)
        SKc = (Yn[:, Kcl] * y1[Kcl]) @ Wsign.T                                                  # n x 32
        prodK = Yn[:, Kl] * y1[Kl]                                                             # n x 8
        out = []
        n_orth = 0
        count_only = '--count' in sys.argv
        rootok = ~((np.abs(Yn[:, Kl]) == 4) & (np.abs(y1[Kl]) == 4)[None, :]).any(axis=1)    # sign-independent
        n_kept = 0
        for s0 in range(0, len(Yn), 256):
            sl = np.nonzero(rootok[s0:s0 + 256])[0] + s0
            if not len(sl):
                continue
            sK = prodK[sl] @ EPS.T                                                              # c x 256
            valid = ~((EPS[None, :, :] == -1) & ~nzK[sl][:, None, :]).any(axis=2)              # no flips on K-zeros
            match = (SKc[sl][:, None, :] == -sK[:, :, None]) & valid[:, :, None]               # c x 256 x 32
            ii, ee, ww = np.nonzero(match)
            n_orth += len(ii)
            if not len(ii):
                continue
            Y2 = Yn[sl][ii].copy()
            Y2[:, Kl] *= EPS[ee]
            Y2[:, Kcl] *= Wsign[ww]
            # root pair filter
            ok = ~((np.abs(Y2[:, Kl]) == 4) & (np.abs(y1[Kl]) == 4)[None, :]).any(axis=1)
            Y2 = Y2[ok]
            if not len(Y2):
                continue
            ok = np.ones(len(Y2), dtype=bool)
            bb = (Y2 @ B32).astype(np.int64)
            ok &= ((bb ** 2 + b1[None, :] ** 2) <= 160).all(axis=1)
            del bb
            if U1.shape[1]:                                        # only the u with |u.y1| = 3 can violate
                ok &= ~(np.abs(Y2 @ U1) == 24).any(axis=1)
            n_kept += int(ok.sum())
            if not count_only:
                out.append(Y2[ok].astype(np.int8))
        if not count_only:
            Y2 = np.concatenate(out) if out else np.zeros((0, 24), dtype=np.int8)
            np.save(os.path.join(HERE, 'state', 'u1_pairs_%02d.npy' % k), Y2)
        print('[%5.0fs] y1 class %d: y2 orthogonal (with signs, root-pair ok) %d, after all pair filters %d'
              % (time.time() - T0, k, n_orth, n_kept), flush=True)


if __name__ == '__main__':
    main()
