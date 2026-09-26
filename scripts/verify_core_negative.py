"""Negative control of verify_core_certs.check_certificate: a tampered certificate must be rejected."""
import gzip, json, sys, copy, random
sys.path.insert(0, 'scripts')
from verify_core_certs import my_columns, check_certificate

data = json.load(gzip.open('logs/core_certs/one_01.json.gz', 'rt', encoding='utf-8'))
dec = [tuple(c) for c in data['decomposition']]
cols = my_columns(dec)
cert = data['certificates'][0]
ok, l, r = check_certificate(cols, cert)
print('original certificate valid:', ok, f'(margin {l - r})')
rng = random.Random(1)
rej = 0
for trial in range(5):
    bad = copy.deepcopy(cert)
    idx = rng.randrange(len(bad['rows']))
    bad['rows'][idx][1] = -bad['rows'][idx][1] * 3 + 7           # corrupt one multiplier
    ok2, _, _ = check_certificate(cols, bad)
    rej += not ok2
bad = copy.deepcopy(cert)
bad['rows'] = [row for row in bad['rows'] if not row[0].startswith('P2')]   # drop all (P2) rows
ok3, _, _ = check_certificate(cols, bad)
print(f'corrupted multipliers rejected: {rej}/5 ; certificate without its (P2) rows valid: {ok3} (expected False)')
