"""One ablation run: python star_ablation_one.py <variant> ; ball = Astra's original ball."""
import sys, time
sys.path.insert(0, 'scripts')
from star_shadows import build, star_test
V = {'all': {}, 'no_lm': {'use_lm': False}, 'no_K': {'use_K': False}, 'no_overlap': {'use_overlap': False},
     'no_diag': {'use_diag': False}, 'only_green': {'use_lm': False, 'use_K': False, 'use_overlap': False}}
a = [9,4,7,8,1,6,5,2,3,0,11,10]; b = [7,2,1,9,11,8,10,0,5,3,6,4]; c = [7,3,10,1,11,9,8,0,6,5,2,4]
adj = build(a, b, c, [(i + 1) % 12 for i in range(12)])
t0 = time.time()
r = star_test(adj, [v for v in range(1, 12) if v != a[0]], full=True, **V[sys.argv[1]])
print(f'{sys.argv[1]}: {"FEASIBLE" if r["solution"] is not None else "dead"}; candidates {r["counts"]}; nodes {r["nodes"]}; {time.time() - t0:.0f}s', flush=True)
