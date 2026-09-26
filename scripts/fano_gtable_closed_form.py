"""Arithmetic check of a HAND-DERIVED closed form (no solver): generic counts at a Fano-tight vertex.
For z in block {P,Q}, M = Phi(P,Q), C4 = the 4 points off the line PQM, matchings M_P, M_Q, M_M of C4 (pairs collinear
with P, Q, M); landing a (via P) in C4, b (via Q) in {Phi(P,a), Phi(Q,a)}; a', b' = M_M-partners of a, b:
  specials: PM, QM, Pa, Qb (1 each);
  g = 2 on both M_M pairs;  g = 1 - [a in pair] on M_P pairs;  g = 1 - [b in pair] on M_Q pairs;
  g(Ma) = [b = Phi(Q,a)], g(Ma') = [b = Phi(P,a)], g(Mb) = [b = Phi(P,a)], g(Mb') = [b = Phi(Q,a)].
Compared with the counts solved from the tight-family and slot equations (fano_cellcross_rank.counts_abstract)."""
import io
import itertools
import contextlib
import sys
sys.path.insert(0, 'scripts')
with contextlib.redirect_stdout(io.StringIO()):
    import fano_cellcross_rank as F

bad = total = 0
for P, Q in itertools.combinations(range(7), 2):
    M = F.third(P, Q)
    C4 = [x for x in range(7) if x not in (P, Q, M)]
    pairs = [frozenset(pr) for pr in itertools.combinations(C4, 2)]
    MP = [p for p in pairs if F.third(*p) == P]
    MQ = [p for p in pairs if F.third(*p) == Q]
    MM = [p for p in pairs if F.third(*p) == M]
    mm = lambda t: next(iter(next(p for p in MM if t in p) - {t}))
    for a in C4:
        for b in {F.third(P, a), F.third(Q, a)}:
            cf = {B: 0 for B in F.BLOCKS}
            for B in (frozenset((P, M)), frozenset((Q, M)), frozenset((P, a)), frozenset((Q, b))):
                cf[B] += 1
            for p in MM:
                cf[p] += 2
            for p in MP:
                cf[p] += 1 - (a in p)
            for p in MQ:
                cf[p] += 1 - (b in p)
            ap, bp = mm(a), mm(b)
            cf[frozenset((M, a))] += int(b == F.third(Q, a))
            cf[frozenset((M, ap))] += int(b == F.third(P, a))
            cf[frozenset((M, b))] += int(b == F.third(P, a))
            cf[frozenset((M, bp))] += int(b == F.third(Q, a))
            total += 1
            bad += [cf[B] for B in F.BLOCKS] != F.counts_abstract(P, Q, a, b)
print(f'closed form vs solved counts: {total} (block pair, landing) cases, {bad} differences')
