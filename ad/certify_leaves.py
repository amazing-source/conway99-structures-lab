"""Certify C2 mixture leaves with the pair-level model (general_model) and write one JSON line per leaf.
Usage: python certify_leaves.py <leaves.json> <job key> <out.jsonl> [lo hi]
  leaves.json: C2's format {job: {twE, In, A, leaves: [{triples, ...}]}}.  Leaves WITH an A-block: the D-structure at A is
  enumerated (130 options for d_A = 6); the leaf is dead iff every option is dead, and every option carries its own tree.
Verification: python certify_leaves.py --check <leaves.json> <job key> <out.jsonl>  (rebuilds every model, checks every tree
in exact rationals)."""
import json, sys, time
import numpy as np
from general_model import build, from_c2_key, pr
from pair_level import two_regular_graphs, cycle_matching
from certify import certify, check


def a_triples(A, cyc):
    return [(A, m, m2) for c in cyc for (m, s), (m2, s2) in cycle_matching(c).items()]


def models_for_leaf(job, leaf):
    A = job['A']
    In = {int(k): v for k, v in job['In'].items()}
    base = [tuple(t) for t in leaf['triples']]
    if A is None:
        yield None, base
    else:
        for cyc in two_regular_graphs(In[A]):
            yield [list(c) for c in cyc], base + a_triples(A, cyc)


def main():
    if sys.argv[1] == '--check':
        _, _, path, key, outp = sys.argv[:5]
        job = json.load(open(path))[key]
        ok = bad = 0
        for line in open(outp):
            rec = json.loads(line)
            if rec['status'] != 'dead':
                continue
            leaf = job['leaves'][rec['leaf']]
            opts = dict((json.dumps(c), t) for c, t in models_for_leaf(job, leaf))
            good = True
            for part in rec['parts']:
                trip = opts[json.dumps(part['A_structure'])]
                A_, lo, hi, lb, ub, *_ = build(from_c2_key(trip))
                good &= check(A_, lo, hi, lb, ub, part['tree'])
            good &= len(rec['parts']) == len(opts)
            ok += good
            bad += not good
        print('checked dead leaves: valid', ok, 'INVALID', bad)
        return
    path, key, outp = sys.argv[1:4]
    job = json.load(open(path))[key]
    lo_, hi_ = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (0, len(job['leaves']))
    with open(outp, 'a') as f:
        for t in range(lo_, min(hi_, len(job['leaves']))):
            t0 = time.time()
            leaf = job['leaves'][t]
            parts = []
            status = 'dead'
            for cyc, trip in models_for_leaf(job, leaf):
                A_, lo, hi, lb, ub, *_ = build(from_c2_key(trip))
                r = certify(A_, lo, hi, lb, ub)
                if r['status'] != 'dead':
                    status = r['status']
                    parts.append({'A_structure': cyc, 'status': r['status']})
                    break
                parts.append({'A_structure': cyc, 'level': r['level'], 'nodes': r['nodes'], 'tree': r['tree']})
            rec = {'leaf': t, 'index': leaf.get('index'), 'status': status,
                   'levels': sorted({p.get('level', p.get('status')) for p in parts}), 'parts': parts,
                   'seconds': round(time.time() - t0, 1)}
            f.write(json.dumps(rec) + '\n')
            f.flush()
            print(t, status, rec['levels'], rec['seconds'], 's', flush=True)


if __name__ == '__main__':
    main()
