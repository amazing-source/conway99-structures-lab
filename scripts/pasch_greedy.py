"""Greedy version of pasch_cheapest (fast): for each block type and target, start from all extras, drop extras while
the minimum over cycles stays >= target (try expensive-looking ones first: doubles >= 2, then others)."""
import sys
sys.path.insert(0, 'scripts')
from pasch_cheapest import block_data, cycle_min, base_val, extras_list

def greedy(j, target):
    others, exact, simple, imm = block_data(__import__('pasch_cheapest').PASCH, j)
    val0 = base_val(exact, simple, imm)
    ex = extras_list(exact, simple, imm)
    val = dict(val0)
    for key, v in ex:
        val[key] = v
    m = cycle_min(others, val)
    if not (m is None or m >= target):
        return None, None
    kept = list(ex)
    for item in sorted(ex, key=lambda kv: (-(kv[1] if kv[1] is not None else 3), str(kv[0]))):
        key, v = item
        trial = dict(val)
        trial[key] = val0[key]
        mm = cycle_min(others, trial)
        if mm is None or mm >= target:
            val = trial
            kept.remove(item)
    return len(kept), kept

if __name__ == '__main__':
    res = {}
    for name, j, targets in (('w', 3, range(0, 7)), ('p', 1, range(0, 9)), ('a', 0, range(0, 8))):
        for t in targets:
            c, kept = greedy(j, t)
            res[(name, t)] = (c, kept)
            print(f'{name}: target {t}: extras {c}', flush=True)
    best = None
    for tw in range(0, 7):
        for tp in range(0, 9):
            for ta in range(0, 8):
                if tw + 3 * tp + 3 * ta < 43: continue
                cs = [res[('w', tw)][0], res[('p', tp)][0], res[('a', ta)][0]]
                if None in cs: continue
                if best is None or sum(cs) < best[0]:
                    best = (sum(cs), tw, tp, ta)
    print('cheapest (greedy):', best)
    for name, t in zip('wpa', best[1:]):
        print(f'  {name} >= {t}:')
        for (i, pr), v in res[(name, t)][1]:
            print(f'     orbit {i} pair {sorted(pr)}: {"EXCLUDED" if v is None else f">= {v}"}')
