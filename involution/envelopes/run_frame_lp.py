"""Batch driver (structures 27/09): frame_lp on the (y1, y2) pairs produced by C1 (JSON {id, pairs: [[y1 coords], [y2 coords]],
stats}, coordinates of U/NNN.json).  Each infeasible pair gets an integer certificate saved with its basis
(frame_certs/NNN_PPPP.npz, checked by verify_frame_cert.py).  Results: frame_results/NNN.json (resumable).
Usage: python run_frame_lp.py pairs_file.json [pairs_file ...] [--jobs N] [--slices S]
  --slices S splits every file into S slices run in parallel (outputs frame_results/NNN_sA_B.json, merge with --merge).
  laptop: --jobs 1 (default) through run_capped; lab PC: --jobs up to the allowed threads."""
import json, os, sys, time
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')                     # one thread per worker (parallel slices)
from fractions import Fraction as F
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.makedirs(os.path.join(HERE, 'frame_certs'), exist_ok=True)
os.makedirs(os.path.join(HERE, 'frame_results'), exist_ok=True)


def run_file(task):
    from build156 import rows_frac
    from frame_lp import prepare, coef_of, frame_lp
    path, lo, hi = task if isinstance(task, tuple) else (task, None, None)
    data = json.load(open(path))
    i = int(data['id'])
    stem = os.path.splitext(os.path.basename(path))[0]                 # outputs keyed by the pairs file (no collisions)
    if lo is not None:
        data['pairs'] = data['pairs'][:hi]
    out_path = os.path.join(HERE, 'frame_results', '%s.json' % stem if lo is None else '%s_s%d_%d.json' % (stem, lo, hi))
    res = json.load(open(out_path)) if os.path.exists(out_path) else {'id': i, 'source': os.path.basename(path), 'pairs': {}}
    if 'basis' in data:                                                   # own frame (e.g. Golay a-coordinates, scale 8)
        P = prepare([[F(x) for x in r] for r in data['basis']], scale=int(data.get('scale', 1)), verbose=False)
        bpath = os.path.join(HERE, 'frame_certs', '%s_basis.npz' % stem)
        np.savez_compressed(bpath, basis=np.array([[str(x) for x in row] for row in P['M']]), den=P['d'],
                            scale=int(data.get('scale', 1)), own_frame=1)
    else:
        rec = json.load(open(os.path.join(HERE, 'U', '%03d.json' % i)))
        P = prepare(rows_frac(rec['gens']), verbose=False)
        bpath = os.path.join(HERE, 'frame_certs', '%03d_basis.npz' % i)      # one basis file per lattice
        np.savez_compressed(bpath, basis=np.array([[str(x) for x in row] for row in P['M']]), den=P['d'])
    t_last = time.time()
    for k, (c1, c2) in enumerate(data['pairs']):
        if str(k) in res['pairs'] or (lo is not None and k < lo):
            continue
        y1, y2 = coef_of(P, [F(x) for x in c1]), coef_of(P, [F(x) for x in c2])
        cert = os.path.join(HERE, 'frame_certs', '%s_%05d.npz' % (stem, k))
        r = frame_lp(P, y1, y2, save=cert)
        res['pairs'][str(k)] = r
        if time.time() - t_last > 20:
            json.dump(res, open(out_path, 'w'))
            t_last = time.time()
    json.dump(res, open(out_path, 'w'))
    ex = sum(1 for v in res['pairs'].values() if v['verdict'].startswith('EXCLUDED'))
    res['summary'] = {'pairs': len(res['pairs']), 'excluded': ex, 'not_excluded': len(res['pairs']) - ex}
    json.dump(res, open(out_path, 'w'))
    return i, res['summary']


def merge(stem):
    parts = [f for f in os.listdir(os.path.join(HERE, 'frame_results')) if f.startswith(stem + '_s')]
    res = {'stem': stem, 'pairs': {}}
    base = os.path.join(HERE, 'frame_results', '%s.json' % stem)
    if os.path.exists(base):
        res['pairs'].update(json.load(open(base))['pairs'])
    for f in parts:
        res['pairs'].update(json.load(open(os.path.join(HERE, 'frame_results', f)))['pairs'])
    ex = sum(1 for v in res['pairs'].values() if v['verdict'].startswith('EXCLUDED'))
    res['summary'] = {'pairs': len(res['pairs']), 'excluded': ex, 'not_excluded': len(res['pairs']) - ex,
                      'not_excluded_ids': sorted(int(k) for k, v in res['pairs'].items() if not v['verdict'].startswith('EXCLUDED'))}
    json.dump(res, open(base, 'w'))
    return res['summary']


def main():
    args = sys.argv[1:]
    jobs, slices = 1, 1
    if '--merge' in args:
        for a in args[args.index('--merge') + 1:]:
            print(a, merge(a))
        return
    if '--jobs' in args:
        k = args.index('--jobs'); jobs = int(args[k + 1]); del args[k:k + 2]
    if '--slices' in args:
        k = args.index('--slices'); slices = int(args[k + 1]); del args[k:k + 2]
    tasks = []
    for p in args:
        if slices == 1:
            tasks.append(p)
        else:
            n = len(json.load(open(p))['pairs'])
            cuts = [n * t // slices for t in range(slices + 1)]
            tasks += [(p, cuts[t], cuts[t + 1]) for t in range(slices)]
    t0 = time.time()
    if jobs == 1:
        for t in tasks:
            print(run_file(t), '%.0fs' % (time.time() - t0), flush=True)
    else:
        with Pool(jobs) as pool:
            for r in pool.imap_unordered(run_file, tasks):
                print(r, '%.0fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
