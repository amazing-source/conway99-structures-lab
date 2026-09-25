"""D2 cross-check of the mod-4 discriminant criterion (disc4_scan.crit) against Brouwer's table (raw HTML, v <= 1300).
Every row carries the colour of its pair: green (#b0ffb0) = exists, yellow (#ffffb0) = unknown, red (#ffb0b0) = no.
Checks: (1) no existing graph (green) is excluded; (2) list excluded rows by status; (3) complement rows agree."""
import glob, re, sys
sys.path.insert(0, 'scripts')
from disc4_scan import crit

rows = []
for path in sorted(glob.glob('logs/brouwer/srgtab*-*.html')):
    html = open(path, encoding='utf-8', errors='replace').read()
    v_cur = None
    for m in re.finditer(r'<tr bgcolor="(#[0-9a-f]+)">(.*?)</tr>', html, re.S):
        colour, body = m.group(1), m.group(2)
        cells = [re.sub(r'<.*?>', '', c).replace('&nbsp;', '').strip() for c in re.findall(r'<td>(.*?)</td>', body, re.S)]
        if len(cells) < 5:
            continue
        mark, vv = cells[0], cells[1]
        if vv:
            v_cur = int(vv)
        try:
            k, l, mu = int(cells[2]), int(cells[3]), int(cells[4])
        except ValueError:
            continue
        status = {'#b0ffb0': 'exists', '#ffffb0': 'unknown', '#ffb0b0': 'no'}.get(colour, colour)
        comment = cells[7] if len(cells) > 7 else ''
        rows.append((v_cur, k, l, mu, status, bool(vv), comment))
print('rows parsed:', len(rows), ' by status:', {s: sum(1 for r in rows if r[4] == s) for s in ('exists', 'unknown', 'no')})
bad, excl = [], []
for v, k, l, mu, st, main, com in rows:
    c = crit(v, k, l, mu)
    if c and c[5]:
        excl.append((v, k, l, mu, st, main, com))
        if st == 'exists':
            bad.append((v, k, l, mu, com))
print('EXISTING graphs excluded (must be 0):', len(bad), bad[:10])
mains = [x for x in excl if x[5]]
comps = [x for x in excl if not x[5]]
print(f'excluded rows: {len(excl)} ({len(mains)} main, {len(comps)} complement rows)')
for st in ('no', 'unknown'):
    lst = [x for x in mains if x[4] == st]
    print(f'  status {st}: {len(lst)}')
    for x in lst:
        print(f'     srg{x[:4]}  {x[6]}')
