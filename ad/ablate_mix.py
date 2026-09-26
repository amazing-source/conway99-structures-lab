import json
from test_mix import *
d = json.load(open(LEAVES)); job = d['class9_job42']
for t, leaf in enumerate(job['leaves'][:10]):
    res = {}
    for name, kw in [('full', {}), ('nodich', {'dichotomy': False}), ('nobal', {'balance': False}), ('neither', {'dichotomy': False, 'balance': False})]:
        res[name] = 'alive' if test_leaf(leaf['triples'], None, None, **kw) else 'dead'
    print(t, leaf['index'], res)
