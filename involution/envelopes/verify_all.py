"""Independent re-verification of the 156 envelopes in U/ (structures 27/09), from each basis alone:
det 1, odd, no norm-1 vector, root system and number of roots = table, and the table key: the two even neighbours of U
(U_0 = even sublattice, chi = characteristic vector, N_a = U_0 + Z chi/2, N_b = U_0 + Z(chi/2 + b), b odd) have the
table's Niemeier root systems (with the proven fix of id 150, see build156.TABLE_FIXES).  Resumable; results in
verify_all.json.  Usage: python verify_all.py [max_seconds]"""
import json, os, sys, time, collections
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build156 import (Lattice, enum, rtype, tkey, tname, parse_type, int_rows, gram_of, det_exact, inv_exact, rows_frac,
                      TABLE, TABLE_FIXES, save_json, load_json)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'verify_all.json')
T0 = time.time()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 270.0


def roots_type(L):
    sv = enum(L.G, 2)
    n1 = sum(1 for z, v in sv if v == 1)
    R = np.array([z for z, v in sv if v == 2], dtype=np.int64)
    return n1, 2 * len(R), (tkey(rtype(R, L.Gm)) if len(R) else ())


def check(rec, e):
    rows = rows_frac(rec['gens'])
    Mi, d = int_rows(rows)
    G = gram_of(Mi, d)
    res = {'det': det_exact(G), 'odd': any(G[i][i] % 2 for i in range(24))}
    U = Lattice(rows)
    n1, nr, rt = roots_type(U)
    res.update(norm1=n1, roots=nr, rtype=tname(rt))
    # characteristic vector and the two even neighbours
    Gi = inv_exact(U.G)
    diag = [U.G[i][i] % 2 for i in range(24)]
    c = [int(sum(Gi[i][j] * diag[j] for j in range(24))) % 2 for i in range(24)]
    assert all(sum(Gi[i][j] * diag[j] for j in range(24)).denominator == 1 for i in range(24))
    k0 = next(i for i in range(24) if diag[i])
    base = []
    for j in range(24):
        e_ = [0] * 24
        if j == k0:
            e_[j] = 2
        else:
            e_[j] = 1
            if diag[j]:
                e_[k0] = 1
        base.append(e_)
    half = [F(x, 2) for x in c]
    halfB = list(half)
    halfB[k0] += 1
    nbs = []
    for glue in (half, halfB):
        N = Lattice([U.coords(b) for b in base] + [U.coords(glue)])
        even = all(N.G[i][i] % 2 == 0 for i in range(24))
        _, nr2, t2 = roots_type(N)
        nbs.append((even, det_exact(N.G), nr2, t2))
    res['neighbours'] = [(tname(t) if t else 'Leech', ev, dt, nr2) for ev, dt, nr2, t in nbs]
    want = sorted([tkey(parse_type(e['nb1'])), tkey(parse_type(e['nb2']))])
    got = sorted([t for ev, dt, nr2, t in nbs])
    res['ok'] = (res['det'] == 1 and res['odd'] and n1 == 0 and nr == e['nroots'] and rt == tkey(parse_type(e['roots']))
                 and all(ev and dt == 1 for ev, dt, nr2, t in nbs) and got == want)
    return res


def main():
    table = {e['id']: e for e in json.load(open(TABLE))}
    for i, fx in TABLE_FIXES.items():
        table[i] = dict(table[i], **fx)
    done = load_json(OUT, {})
    for fn in sorted(os.listdir(os.path.join(HERE, 'U'))):
        if time.time() - T0 > BUDGET:
            break
        rec = json.load(open(os.path.join(HERE, 'U', fn)))
        if str(rec['id']) in done:
            continue
        r = check(rec, table[rec['id']])
        done[str(rec['id'])] = r
        save_json(OUT, done)
        print('[%5.0fs] id %3d %s %s' % (time.time() - T0, rec['id'], 'OK ' if r['ok'] else 'BAD', r), flush=True)
    bad = [k for k, v in done.items() if not v['ok']]
    print('checked %d/156, bad %s' % (len(done), bad))


if __name__ == '__main__':
    main()
