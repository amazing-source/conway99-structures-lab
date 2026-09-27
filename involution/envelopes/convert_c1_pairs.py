"""Adapter (structures 27/09): C1's pair files (conway99-c1 data/finishC/ypairs/yp_NNN.json: integer coefficients in C1's basis
LLL(HNF(gens)) computed by C1's envelope_trap.lll / to_int_basis on U/NNN.json) -> the coordinate format read by
run_frame_lp.py ({id, pairs: [[y1 Q^24 strings], [y2 ...]], stats}).  Checks norms (7) and orthogonality after conversion.
Usage: python convert_c1_pairs.py yp_NNN.json [first last]  -> state/pairs_NNN[_a_b].json"""
import json, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
C1_CODE = os.environ.get('C1_CODE', 'C:/home/researcher/Desktop/c1_recherche/code')
sys.path.insert(0, C1_CODE)
from envelope_trap import to_int_basis, lll


def main():
    src = sys.argv[1]
    d = json.load(open(src))
    i = int(d['id'])
    gens = json.load(open(os.path.join(HERE, 'U', '%03d.json' % i)))['gens']
    B = lll(to_int_basis(gens))                                   # C1's basis (rational rows)
    pairs = d['pairs']
    a = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    b = int(sys.argv[3]) if len(sys.argv) > 3 else len(pairs)
    coords = lambda c: [sum(F(int(c[k])) * F(B[k][j]) for k in range(24)) for j in range(24)]
    dot = lambda u, v: sum(p * q for p, q in zip(u, v))
    out = []
    for y1c, y2c in pairs[a:b]:
        y1, y2 = coords(y1c), coords(y2c)
        assert dot(y1, y1) == 7 and dot(y2, y2) == 7 and dot(y1, y2) == 0
        out.append([[str(x) for x in y1], [str(x) for x in y2]])
    tag = '' if (a, b) == (0, len(pairs)) else '_%d_%d' % (a, b)
    dst = os.path.join(HERE, 'state', 'pairs_%03d%s.json' % (i, tag))
    json.dump({'id': i, 'pairs': out, 'stats': {'source': os.path.basename(src), 'range': [a, b], 'info': d.get('info')}},
              open(dst, 'w'))
    print('wrote', dst, len(out), 'pairs')


if __name__ == '__main__':
    main()
