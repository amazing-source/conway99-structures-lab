"""Pair-level test of C2's mixture leaves ('one block' skeletons, Sigma = empty) with the general model.
Leaves list only the orbits whose pivot is not the A-block; the D-structure at A (130 options) is enumerated."""
import json, sys, time
from general_model import feasible, from_c2_key, N, pr
from pair_level import two_regular_graphs, cycle_matching

LEAVES = r"C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/c2exp/logs/leaves_for_structures.json"


def a_triples(A, In_A, cycles):
    out = []
    for cyc in cycles:
        for (m, s), (m2, s2) in cycle_matching(cyc).items():
            out.append((A, m, m2))
    return out


def test_leaf(triples, A, In_A, dichotomy=True, balance=True, stop_first=True):
    if A is None:
        r = feasible(from_c2_key(triples), dichotomy=dichotomy, balance=balance)
        return [] if r is None else [('noA', r if isinstance(r, str) else 'feasible')]
    alive = []
    for cyc in two_regular_graphs(In_A):
        orbs = from_c2_key(list(map(tuple, triples)) + a_triples(A, In_A, cyc))
        r = feasible(orbs, dichotomy=dichotomy, balance=balance)
        if r is not None:
            alive.append((cyc, r if isinstance(r, str) else 'feasible'))
            if stop_first:
                break
    return alive


if __name__ == '__main__':
    d = json.load(open(LEAVES))
    which = sys.argv[1]
    job = d[which]
    A = job['A']
    In = {int(k): v for k, v in job['In'].items()}
    lo_, hi_ = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, len(job['leaves']))
    t0 = time.time()
    dead = 0
    for t in range(lo_, min(hi_, len(job['leaves']))):
        leaf = job['leaves'][t]
        al = test_leaf(leaf['triples'], A, In[A] if A is not None else None)
        tag = 'opt1' if leaf.get('block0_option1') else ''
        if al:
            print('leaf', t, leaf.get('index', ''), tag, 'ALIVE', al[0][0]); sys.stdout.flush()
        else:
            dead += 1
            print('leaf', t, leaf.get('index', ''), tag, 'dead  (%.0fs)' % (time.time() - t0)); sys.stdout.flush()
    print('dead', dead, 'of', min(hi_, len(job['leaves'])) - lo_)
