"""Structure of the (P2)/(PPE) part of a certificate, seen as a symmetric 21x21 matrix Y on pairs (P2 rows are the
off-diagonal entries of U^T U = M, PPE rows the diagonal ones).  Reports the spectrum of Y, of its projection on the
14-dim space W14 (block sums zero), and the weight of the other row families."""
import ast, gzip, json, sys, itertools
import numpy as np
sys.path.insert(0, 'scripts')
from zeros_budget import V

fn = sys.argv[1]
data = json.load(gzip.open(fn, 'rt', encoding='utf-8'))
cert = data['certificates'][0]
PAIRS = [frozenset(q) for q in itertools.combinations(V, 2)]
idx = {Q: t for t, Q in enumerate(PAIRS)}
Y = np.zeros((21, 21))
fam = {}
for label, y in cert['rows']:
    kind = label.split(' ', 1)[0]
    fam[kind] = fam.get(kind, 0) + abs(y)
    if kind == 'P2':
        Q1, Q2 = [frozenset(ast.literal_eval(s)) for s in label.split(' ', 1)[1].split('|')]
        Y[idx[Q1], idx[Q2]] += y / 2
        Y[idx[Q2], idx[Q1]] += y / 2
tot = sum(fam.values())
print('weight share by family:', {k: round(v / tot, 3) for k, v in sorted(fam.items(), key=lambda kv: -kv[1])})
Y /= np.abs(Y).max()
AT = np.array([[1.0 if len(P & Q) == 1 else 0.0 for Q in PAIRS] for P in PAIRS])
w, E = np.linalg.eigh(AT)
P14 = E[:, np.abs(w + 2) < 1e-9]
ev = np.linalg.eigvalsh(Y)
ev14 = np.linalg.eigvalsh(P14.T @ Y @ P14)
print('decomposition', data['decomposition'])
print('spectrum of Y (normalised):', np.round(ev, 3))
print('spectrum of Y on W14       :', np.round(ev14, 3))
