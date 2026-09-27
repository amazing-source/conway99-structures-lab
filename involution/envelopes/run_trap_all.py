"""Driver (structures 27/09): C1's proven trap criterion (conway99-c1 code/envelope_trap.py, trap_verdict) on the 156
explicit envelopes U/NNN.json.  Records the sha256 of the tool used.  Results merged into trap_results.json (resumable).
Usage: C1_CODE=<conway99-c1/code> python run_trap_all.py [--jobs N] [id ...]
  laptop: --jobs 1 only (one light process); lab PC: --jobs 10."""
import hashlib, json, os, sys, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
C1_CODE = os.environ.get('C1_CODE', 'C:/home/researcher/Desktop/c1_recherche/code')
sys.path.insert(0, C1_CODE)
OUT = os.path.join(HERE, 'trap_results.json')


def one(fn):
    from envelope_trap import trap_verdict
    rec = json.load(open(os.path.join(HERE, 'U', fn)))
    t = time.time()
    try:
        r = trap_verdict(rec['gens'])
    except Exception as ex:                       # report, never hide
        r = {'error': repr(ex)}
    r['seconds'] = round(time.time() - t, 1)
    return rec['id'], rec['roots'], r


def main():
    args = sys.argv[1:]
    jobs = 1
    if '--jobs' in args:
        k = args.index('--jobs')
        jobs = int(args[k + 1])
        del args[k:k + 2]
    want = {int(a) for a in args}
    sha = hashlib.sha256(open(os.path.join(C1_CODE, 'envelope_trap.py'), 'rb').read()).hexdigest()
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    fns = [f for f in sorted(os.listdir(os.path.join(HERE, 'U')))
           if (not want or int(f[:3]) in want) and str(int(f[:3])) not in res]
    print('envelope_trap.py sha256 %s; %d envelopes to run with %d jobs' % (sha, len(fns), jobs), flush=True)
    with Pool(jobs) as pool:
        for i, roots, r in pool.imap_unordered(one, fns):
            r['tool_sha256'] = sha
            res[str(i)] = dict(r, roots=roots)
            tmp = OUT + '.tmp'
            json.dump(res, open(tmp, 'w'), indent=0)
            os.replace(tmp, OUT)
            print('id %3d %-22s %s' % (i, roots or '(none)', {k: v for k, v in r.items() if k != 'tool_sha256'}), flush=True)
    ex = sorted(int(k) for k, v in res.items() if str(v.get('verdict', '')).startswith('EXCLUDED'))
    print('done %d/156; EXCLUDED by the trap: %d %s' % (len(res), len(ex), ex))


if __name__ == '__main__':
    main()
