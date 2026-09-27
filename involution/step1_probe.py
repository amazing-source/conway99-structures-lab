"""Is the near-design (step 1 = (i)+(iii)) itself infeasible for a type?  CP-SAT with a long single call."""
import json, sys, time
import iota_design as D
cfg = json.load(open(sys.argv[1])); tl = float(sys.argv[2])
t0 = time.time()
st, sol, edges = D.step1(cfg['mu'], cfg['tau'], 1, tl)
print(sys.argv[1], 'step 1 (near-design, (i)+(iii)):', st, '%.0fs' % (time.time() - t0), flush=True)
