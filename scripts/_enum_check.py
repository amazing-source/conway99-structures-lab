import sys; sys.path.insert(0, 'scripts')
import genus_scan as gs
V = 200
brute = set()
for v in range(5, V + 1):
    for k in range(2, (v - 1) // 2 + 1):
        for l in range(0, k):
            num = k * (k - l - 1)
            if num % (v - k - 1): continue
            mu = num // (v - k - 1)
            if mu < 1 or mu > k: continue
            if gs.feasible(v, k, l, mu): brute.add((v, k, l, mu))
fast = {x['params'] for x in gs.scan(V)}
print(len(brute), len(fast), "missing in fast:", sorted(brute - fast)[:20], "extra:", sorted(fast - brute)[:20])
