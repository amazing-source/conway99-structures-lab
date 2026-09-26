"""TEST (at most 5 minutes).  Integer version of the configuration model with (P2) + subtour cuts on the one-triangle
classes whose LP relaxation is feasible (the 5 most symmetric classes)."""
import sys, time
sys.path.insert(0, 'scripts')
from zeros_budget import decompositions
from one_triangle_classes import canon1
from config_model import build, ip_feasible

t0 = time.time()
classes = {}
for dec in decompositions((1 << 21) - 1, []):
    if sum(1 for c in dec if len(c) == 3) != 1:
        continue
    form, rep = canon1(dec)
    classes.setdefault((tuple(sorted(len(c) for c in dec)), form), [rep, 0])[1] += 1
todo = [(k, v) for k, v in sorted(classes.items()) if v[1] in (630, 1260)]
print(f'{len(todo)} symmetric classes to test  [{time.time() - t0:.0f} s]', flush=True)
for (typ, form), (rep, cnt) in todo:
    left = 285 - (time.time() - t0)
    if left < 20:
        print('time budget reached'); break
    A, lo, hi, cols = build(rep, p2=True)
    r = ip_feasible(A, lo, hi, min(45, left - 10))
    print(f'  {typ} x{cnt}: IP with P2 {"INFEASIBLE" if r is False else ("FEASIBLE" if r else "undecided")}  '
          f'{rep}  [{time.time() - t0:.0f} s]', flush=True)
