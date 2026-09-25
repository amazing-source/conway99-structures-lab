"""PASCH: per block, (a) a minimal set of EXCLUSIONS whose consequence is that only the exact compatible orders survive;
(b) the doubles facts needed on those orders (per-pair exact minima, and the smallest values that still give the target)."""
import sys
sys.path.insert(0, 'scripts')
from pasch_minimal_facts import table, PASCH
from zeros_budget import CYC6

def surviving(others, excl):
    out = []
    for order in CYC6:
        lab = [others[p] for p in order]
        if all((lab[t], frozenset((lab[t - 1], lab[(t + 1) % 6]))) not in excl for t in range(6)):
            out.append(tuple(lab))
    return out

def canon_cycle(lab):
    n = len(lab); best = None
    for s in range(n):
        for seq in (lab[s:] + lab[:s], (lab[s::-1] + lab[:s:-1])):
            t = tuple(seq)
            if best is None or t < best: best = t
    return best

def run(j, target):
    others, tab, kk = table(PASCH, j)
    all_excl = {key for key, v in tab.items() if v is None}
    exact = {canon_cycle(list(o)) for o in surviving(others, all_excl)}
    excl = set(all_excl)
    for key in sorted(all_excl, key=str):
        trial = excl - {key}
        if {canon_cycle(list(o)) for o in surviving(others, trial)} == exact:
            excl = trial
    print(f'=== block {j} (k {kk}): exact compatible orders {sorted(exact)}; minimal exclusions needed: {len(excl)}')
    for i, pr in sorted(excl, key=lambda kp: (kp[0], sorted(kp[1]))):
        print(f'   orbit {i}: {sorted(pr)} excluded')
    needed = {}
    for cyc in sorted(exact):
        lab = list(cyc); tot = 0; parts = []
        for t in range(6):
            key = (lab[t], frozenset((lab[t - 1], lab[(t + 1) % 6])))
            tot += tab[key]; parts.append(f'{lab[t]}:{sorted(key[1])}={tab[key]}')
            needed[key] = tab[key]
        print(f'   order {cyc}: doubles {tot}  [{", ".join(parts)}]')
    print(f'   doubles facts used on these orders: {sum(1 for v in needed.values() if v)} nonzero')

if __name__ == '__main__':
    run(1, 8)
    run(0, 7)
