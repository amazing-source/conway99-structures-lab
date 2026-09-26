"""Minimal infeasible sub-stars (all constraints, full-star equality inside the subset): for each subset of the 10
groups of a given size, is there a choice of shadows with every member's matching degree equal to its capacity?"""
import sys, time, itertools
sys.path.insert(0, 'scripts')
from star_shadows import build, star_test
a = [9,4,7,8,1,6,5,2,3,0,11,10]; b = [7,2,1,9,11,8,10,0,5,3,6,4]; c = [7,3,10,1,11,9,8,0,6,5,2,4]
adj = build(a, b, c, [(i + 1) % 12 for i in range(12)])
groups = [v for v in range(1, 12) if v != a[0]]
size = int(sys.argv[1]); budget = float(sys.argv[2])
t0 = time.time(); dead = []; alive = 0; done = 0
for sub in itertools.combinations(groups, size):
    if time.time() - t0 > budget:
        break
    r = star_test(adj, list(sub), full=True)
    done += 1
    if r['solution'] is None:
        dead.append(sub)
    else:
        alive += 1
print(f'size {size}: tested {done} subsets in {time.time() - t0:.0f}s: feasible {alive}, infeasible {len(dead)}; first infeasible: {dead[:5]}', flush=True)
