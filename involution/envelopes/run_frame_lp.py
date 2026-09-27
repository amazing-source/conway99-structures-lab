"""Batch driver (structures 27/09): frame_lp on the (y1, y2) pairs produced by C1 (JSON {id, pairs: [[y1 coords], [y2 coords]],
stats}, coordinates of U/NNN.json).  Each infeasible pair gets an integer certificate saved with its basis
(frame_certs/NNN_PPPP.npz, checked by verify_frame_cert.py).  Results: frame_results/NNN.json (resumable).
Usage: python run_frame_lp.py pairs_file.json [pairs_file ...] [--jobs N]
  laptop: --jobs 1 (default) through run_capped; lab PC: --jobs up to the allowed threads."""
import json, os, sys, time
from fractions import Fraction as F
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.makedirs(os.path.join(HERE, 'frame_certs'), exist_ok=True)
os.makedirs(os.path.join(HERE, 'frame_results'), exist_ok=True)


def run_file(path):
    from build156 import rows_frac
    from frame_lp import prepare, coef_of, frame_lp
    data = json.load(open(path))
    i = int(data['id'])
    out_path = os.path.join(HERE, 'frame_results', '%03d.json' % i)
    res = json.load(open(out_path)) if os.path.exists(out_path) else {'id': i, 'source': os.path.basename(path), 'pairs': {}}
    rec = json.load(open(os.path.join(HERE, 'U', '%03d.json' % i)))
    P = prepare(rows_frac(rec['gens']), verbose=False)
    M = np.array(P['M'], dtype=object)
    for k, (c1, c2) in enumerate(data['pairs']):
        if str(k) in res['pairs']:
            continue
        y1, y2 = coef_of(P, [F(x) for x in c1]), coef_of(P, [F(x) for x in c2])
        cert = os.path.join(HERE, 'frame_certs', '%03d_%05d.npz' % (i, k))
        r = frame_lp(P, y1, y2, save=cert)
        if os.path.exists(cert):
            z = dict(np.load(cert))
            np.savez_compressed(cert, **z, basis=np.array([[str(x) for x in row] for row in P['M']]), den=P['d'])
        res['pairs'][str(k)] = r
        json.dump(res, open(out_path, 'w'), indent=0)
    ex = sum(1 for v in res['pairs'].values() if v['verdict'].startswith('EXCLUDED'))
    res['summary'] = {'pairs': len(res['pairs']), 'excluded': ex, 'not_excluded': len(res['pairs']) - ex}
    json.dump(res, open(out_path, 'w'), indent=0)
    return i, res['summary']


def main():
    args = sys.argv[1:]
    jobs = 1
    if '--jobs' in args:
        k = args.index('--jobs'); jobs = int(args[k + 1]); del args[k:k + 2]
    t0 = time.time()
    if jobs == 1:
        for p in args:
            print(run_file(p), '%.0fs' % (time.time() - t0), flush=True)
    else:
        with Pool(jobs) as pool:
            for r in pool.imap_unordered(run_file, args):
                print(r, '%.0fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
