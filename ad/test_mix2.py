"""Mixture leaves with the A-structure left FREE (partner variables) -- one ILP per leaf."""
import json, sys, time
from general_model import feasible, from_c2_key, N, pr
from test_mix import LEAVES, test_leaf

def orbs_freeA(triples, A, In_A):
    orbs = from_c2_key([tuple(t) for t in triples])
    for m in In_A:
        for s in (0, 1):
            orbs.append({'pair': pr(A, m), 'pivot': A, 'partner': None})
    return orbs

if __name__ == '__main__':
    d = json.load(open(LEAVES)); which = sys.argv[1]; job = d[which]
    A = job['A']; In = {int(k): v for k, v in job['In'].items()}
    lo_, hi_ = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, len(job['leaves']))
    t0 = time.time(); dead = 0; alive = []
    for t in range(lo_, min(hi_, len(job['leaves']))):
        leaf = job['leaves'][t]
        orbs = orbs_freeA(leaf['triples'], A, In[A]) if A is not None else from_c2_key(leaf['triples'])
        r = feasible(orbs)
        tag = 'opt1' if leaf.get('block0_option1') else ''
        st = 'ALIVE' if r is not None else 'dead'
        if r is None: dead += 1
        else: alive.append(t)
        print('leaf', t, tag, st, '(%.0fs)' % (time.time() - t0)); sys.stdout.flush()
    print('dead', dead, 'of', min(hi_, len(job['leaves'])) - lo_, 'alive', alive)
