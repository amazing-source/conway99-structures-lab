"""The two known members of the family srg(2m^2+1, 2m, 1, 2): m = 2 (H(2,3), 9 vertices) and m = 11 (Berlekamp-van
Lint-Seidel, coset graph of the perfect ternary Golay code, 243 vertices). Both are Cayley graphs on F_3^n."""
import itertools

def golay_parity_check():
    """Parity-check matrix (5 x 11) of the ternary Golay code: cyclic QR code of length 11 over F_3 with
    generator g(x) = x^5 + x^4 - x^3 + x^2 - 1; check polynomial h = (x^11 - 1)/g."""
    p = 3
    g = [2, 0, 1, 2, 1, 1]                                   # coefficients of 1, x, ..., x^5 (-1 = 2)
    # polynomial division (x^11 - 1) / g over F_3
    num = [p - 1] + [0] * 10 + [1]                           # -1 + x^11
    q = [0] * 7
    r = num[:]
    for d in range(11 - 5, -1, -1):
        c = r[d + 5] * pow(g[5], -1, p) % p
        q[d] = c
        for i in range(6):
            r[d + i] = (r[d + i] - c * g[i]) % p
    assert all(x == 0 for x in r), 'g does not divide x^11 - 1'
    h = q                                                    # degree 6
    # generator matrix rows: shifts of g ; code = multiples of g (dimension 6)
    G = [[0] * 11 for _ in range(6)]
    for s in range(6):
        for i, c in enumerate(g):
            G[s][(s + i) % 11] = c
    code = set()
    for coeffs in itertools.product(range(3), repeat=6):
        w = tuple(sum(coeffs[s] * G[s][j] for s in range(6)) % 3 for j in range(11))
        code.add(w)
    assert len(code) == 729
    dmin = min(sum(1 for x in w if x) for w in code if any(w))
    assert dmin == 5, dmin
    # parity-check matrix: rows = reversed shifts of h (dimension 5), checked against the code
    hr = h[::-1]
    Hm = [[0] * 11 for _ in range(5)]
    for s in range(5):
        for i, c in enumerate(hr):
            Hm[s][(s + i) % 11] = c
    assert all(sum(Hm[s][j] * w[j] for j in range(11)) % 3 == 0 for s in range(5) for w in code)
    return Hm, code

def cayley_F3(n, conn):
    """Cayley graph on F_3^n; vertices indexed 0..3^n-1 (base 3), connection set conn (list of n-tuples)."""
    V = list(itertools.product(range(3), repeat=n))
    idx = {v: i for i, v in enumerate(V)}
    adj = [set() for _ in V]
    for i, v in enumerate(V):
        for c in conn:
            w = tuple((a + b) % 3 for a, b in zip(v, c))
            adj[i].add(idx[w])
    return V, adj

def hamming_9():
    conn = [(1, 0), (2, 0), (0, 1), (0, 2)]
    return cayley_F3(2, conn)

def bvls_243():
    Hm, _ = golay_parity_check()
    cols = [tuple(Hm[s][j] for s in range(5)) for j in range(11)]
    conn = cols + [tuple((-x) % 3 for x in c) for c in cols]
    assert len(set(conn)) == 22
    return cayley_F3(5, conn)

def check_srg(adj, k, lam, mu):
    n = len(adj)
    for a in range(n):
        if len(adj[a]) != k:
            return False
        for b in range(a + 1, n):
            c = len(adj[a] & adj[b])
            if c != (lam if b in adj[a] else mu):
                return False
    return True

if __name__ == '__main__':
    for name, (V, adj), k in [('H(2,3)', hamming_9(), 4), ('BvLS', bvls_243(), 22)]:
        print(name, len(adj), 'srg(%d,%d,1,2):' % (len(adj), k), check_srg(adj, k, 1, 2))
