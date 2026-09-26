"""Orientation-level pre-filter for C2's one-block mixtures: pair types known (anchored via In, twisted via twE), every D-structure
free (general_model free partners), balance only on certain unions, plus C2's (*) for every (k,S) (exact at integers)."""
import json, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from general_model import build, PAIRS, pr, N
from star_generic import augment
SAMPLE = r"C:/home/researcher/AppData/Local/Temp/claude/C--Users-RESEARCHER-Desktop/c88c221a-efb9-4373-af95-1e500b2e81aa/scratchpad/c2exp/logs/orient_sample_for_structures.json"
def orbs_from_orientation(twE, In):
    tw = {pr(a, b) for a, b in twE}
    anchor = {}
    for c, xs in In.items():
        for x in xs: anchor[pr(int(c), x)] = int(c)
    orbs = []
    for P in PAIRS:
        if P in tw:
            orbs.append({'pair': P, 'pivot': P[0], 'partner': None})
            orbs.append({'pair': P, 'pivot': P[1], 'partner': None})
        else:
            c = anchor[P]
            orbs.append({'pair': P, 'pivot': c, 'partner': None})
            orbs.append({'pair': P, 'pivot': c, 'partner': None})
    return orbs
def test(orbs, star=True, tl=120, lp=False):
    A, lo, hi, lb, ub, nvar, phi, ns = build(orbs)
    if star:
        A, lo, hi, lb, ub, nv = augment(A, lo, hi, lb, ub, nvar, phi, {i: o['pair'] for i, o in enumerate(orbs)})
    else:
        nv = nvar
    integ = np.concatenate([np.ones(nvar), np.zeros(nv - nvar)]) if not lp else np.zeros(nv)
    res = milp(np.zeros(nv), constraints=[LinearConstraint(A, lo, hi)], integrality=integ, bounds=Bounds(lb, ub), options={'time_limit': tl})
    return {0: 'alive', 2: 'DEAD'}.get(res.status, 'unknown')
if __name__ == '__main__':
    d = json.load(open(SAMPLE))
    mode = sys.argv[1] if len(sys.argv) > 1 else 'lp'
    for ci, cls in d.items():
        for o in cls['orientations']:
            orbs = orbs_from_orientation(cls['twE'], o['In'])
            t0 = time.time()
            r0 = test(orbs, star=False, lp=(mode == 'lp'), tl=60)
            r1 = test(orbs, star=True, lp=(mode == 'lp'), tl=120) if r0 != 'DEAD' else 'DEAD'
            print('class', ci, 'idx', o['idx'], 'no(*):', r0, ' with(*):', r1, '%.0fs' % (time.time() - t0), flush=True)
