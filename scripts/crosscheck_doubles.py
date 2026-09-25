"""Cross-checks for doubles_budget3:
 (1) every local key used in the full run: label-level (independent) vs count-level.
     count None and label defined  -> ERROR (count model has a false constraint);
     both defined and label < count -> ERROR;  label > count or label None with count defined: fine (count is a relaxation).
 (2) canonicalisation: analyse(dec, canonical=True) == analyse(dec, canonical=False) on a sample of decompositions."""
import json, random, sys, time
sys.path.insert(0, 'scripts')
import labellevel_check as LL
import doubles_budget3 as d3
from zeros_budget import decompositions

keys = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'logs/doubles3_keys.json'))
t0 = time.time()
err, eq, gt, lnone = [], 0, 0, 0
for tris, nbrs, cval in keys:
    key = (tuple(tuple(t) for t in tris), tuple(nbrs))
    lval = LL.solve(key[0], key[1])
    if cval is None and lval is not None:
        err.append((key, cval, lval))
    elif cval is not None and lval is not None:
        if lval < cval:
            err.append((key, cval, lval))
        elif lval == cval:
            eq += 1
        else:
            gt += 1
    elif cval is not None and lval is None:
        lnone += 1
print(f'(1) {len(keys)} keys: equal {eq}, label > count {gt}, label infeasible but count feasible {lnone}, '
      f'ERRORS {len(err)}  [{time.time() - t0:.0f} s]')
for e in err[:10]:
    print('   ERROR', e)

random.seed(7)
decs = [d for d in decompositions((1 << 21) - 1, []) if sum(1 for c in d if len(c) == 3) >= 1]
sample = random.sample(decs, 300)
mism = 0
for dec in sample:
    a = d3.analyse(dec, canonical=True)[0]
    b = d3.analyse(dec, canonical=False)[0]
    mism += a != b
print(f'(2) canonical vs fixed relabelling on {len(sample)} random decompositions: mismatches {mism}')
