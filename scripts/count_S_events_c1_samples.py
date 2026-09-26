"""Answer to C1 (notes/125 sec. 7): in the tau/sigma data of the x19 certificate files (samples of notes/117-118),
count, for every class U and label v, the three S-events: sigma_U(v) = vbar, sigma_U(v) = tau_u(v), sigma_U(v) = tau_ubar(v)
(u = (c,0), ubar = (c,1); sigma[c] maps the label of the u-pair {u,v} to the label of the adjacent ubar-pair)."""
import json, glob, os, sys
lab = lambda s: (int(s.split('.')[0]), int(s.split('.')[1]))
mate = lambda l: (l[0], 1 - l[1])
files = sorted(glob.glob(os.path.join(sys.argv[1], 'x19*.json')))
for f in files:
    D = json.load(open(f))
    if 'tau' not in D or 'sigma' not in D:
        continue
    tau = {}
    for u, lst in D['tau'].items():
        m = {}
        for v, w in lst:
            m[lab(v)] = lab(w); m[lab(w)] = lab(v)
        tau[lab(u)] = m
    ev = {'vbar': 0, 'tau_u': 0, 'tau_ubar': 0}
    for c, lst in D['sigma'].items():
        c = int(c)
        u, ub = (c, 0), (c, 1)
        for v, w in lst:
            v, w = lab(v), lab(w)
            ev['vbar'] += (w == mate(v))
            ev['tau_u'] += (w == tau[u][v])
            ev['tau_ubar'] += (w == tau[ub][v])
    print(f'{os.path.basename(f):32s} events: sigma=vbar {ev["vbar"]:2d} | sigma=tau_u {ev["tau_u"]:2d} | sigma=tau_ubar {ev["tau_ubar"]:2d}')
