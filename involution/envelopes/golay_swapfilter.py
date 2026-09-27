"""id 1 (a1^8): swap reduction of the pair files (structures 27/09).  A pair (y1, y2) with y1 in class k and y2 in class
l < k is equivalent (swap y1 <-> y2, then G1) to a pair with first vector in class l, which the class-l file already covers.
So in state/pairs_u1cKK_001.json keep only the y2 whose own G1-class (state/u1_ylabels.npy, same order as u1_ynorm.npy)
is >= k.  Class of a y2: sign-normalise it (K signs -> +, K^c signs -> canonical inner-codeword coset representative), then
look it up among the sign-normalised y (u1_ynorm.npy).  Writes state/pairs_u1cKK_001s.json (same format).
Usage: python golay_swapfilter.py k [k ...]"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from golay_family import WORDS, OM, u1_lattice


def main():
    K, _ = u1_lattice()
    Kc = OM ^ K
    Kl = [i for i in range(24) if K >> i & 1]
    Kcl = [i for i in range(24) if Kc >> i & 1]
    inner = [w & Kc for w in WORDS if w & K == 0]
    Yn = np.load(os.path.join(HERE, 'state', 'u1_ynorm.npy'))
    lab = np.load(os.path.join(HERE, 'state', 'u1_ylabels.npy'))
    index = {tuple(int(x) for x in v): int(lab[i]) for i, v in enumerate(Yn)}

    def normalise(a):
        a = list(a)
        for j in Kl:
            a[j] = abs(a[j])
        X = 0
        for j in Kcl:
            if a[j] % 4 == 3:
                X |= 1 << j
        best, bw = None, 0
        for w in inner:
            if best is None or (X ^ w) < best:
                best, bw = X ^ w, w
        for j in Kcl:
            if bw >> j & 1:
                a[j] = -a[j]
        return tuple(a)
    for k in [int(x) for x in sys.argv[1:]]:
        src = os.path.join(HERE, 'state', 'pairs_u1c%02d_001.json' % k)
        d = json.load(open(src))
        keep = []
        drop = 0
        for p in d['pairs']:
            y2 = [int(x) for x in p[1]]
            l = index.get(normalise(y2))
            if l is None:
                l = index.get(normalise([-x for x in y2]))
            assert l is not None, 'y2 not found among the admissible y'
            if l >= k:
                keep.append(p)
            else:
                drop += 1
        d['pairs'] = keep
        d['stats']['swap_filter'] = {'kept': len(keep), 'dropped_lower_class': drop}
        json.dump(d, open(os.path.join(HERE, 'state', 'pairs_u1c%02d_001s.json' % k), 'w'))
        print('class %d: kept %d, dropped %d (y2 in a lower class)' % (k, len(keep), drop), flush=True)


if __name__ == '__main__':
    main()
