import sys
import time

GAMMA = 0.5772156649015329

# CONSTANTS


def harm(n):
    from math import log

    return log(n) + GAMMA + 1 / (2 * n)


def phi_q(q):
    from math import pi

    n = -(-(q - 2) // 2)
    return (4 / pi) * q + (2 * q / pi) * harm(n) + (1 - 2 / pi) * (q - 2) + 0.727


def pb_step3(q, m):
    return m**0.5 + phi_q(q) / q


def psi2(q):
    from math import pi

    p = q // 2
    h = harm(p - 1)
    if q % 2 == 0:
        return (q / pi) * (2 * h - 1 + 1 / p) + (1 - 2 / pi) * q / 2
    return (q / pi) * (2 * h - 1 + 2 / p) + (1 - 2 / pi) * (q / 2 + 1 / (2 * q))


def pb_chord(q, e0):
    from math import cos, pi

    sec = 1 / cos(pi * (e0 - (q - 1) / 2) / q)
    return ((4 / pi) * q + psi2(q) + q / 2 - sec / 2) / q


def gap_hp(q, e0, b, chord):
    import mpmath

    mpmath.mp.dps = 40
    q = mpmath.mpf(q)
    pi = mpmath.pi
    g = mpmath.euler
    if chord:
        p = mpmath.floor(q / 2)
        h = mpmath.log(p - 1) + g + 1 / (2 * (p - 1))
        if int(q) % 2 == 0:
            ps = (q / pi) * (2 * h - 1 + 1 / p) + (1 - 2 / pi) * q / 2
        else:
            ps = (q / pi) * (2 * h - 1 + 2 / p) + (1 - 2 / pi) * (q / 2 + 1 / (2 * q))
        sec = 1 / mpmath.cos(pi * (e0 - (q - 1) / 2) / q)
        pb = ((4 / pi) * q + ps + q / 2 - sec / 2) / q
    else:
        n = mpmath.ceil((q - 2) / 2)
        ph = (
            (4 / pi) * q
            + (2 * q / pi) * (mpmath.log(n) + g + 1 / (2 * n))
            + (1 - 2 / pi) * (q - 2)
            + mpmath.mpf("0.727")
        )
        pb = 1 + ph / q
    return (q - 1) * q ** (-mpmath.mpf(b)) - pb


def scan(chord, worst, lo, hi, b=0.8):
    from math import pi

    import numpy

    q = numpy.arange(lo, hi + 1, dtype=numpy.float64)
    if chord:
        p = numpy.floor(q / 2)
        h = numpy.log(p - 1) + GAMMA + 1 / (2 * (p - 1))
        ex = numpy.where(q % 2 == 0, 1 / p, 2 / p)
        fl = numpy.where(q % 2 == 0, q / 2, q / 2 + 1 / (2 * q))
        ps = (q / pi) * (2 * h - 1 + ex) + (1 - 2 / pi) * fl
        c = numpy.where(q % 2 == 1, 0.0, 0.5) if worst else (q - 1) / 2
        pb = ((4 / pi) * q + ps + q / 2 - 1 / numpy.cos(pi * c / q) / 2) / q
    else:
        n = numpy.ceil((q - 2) / 2)
        pb = (
            1
            + (
                (4 / pi) * q
                + (2 * q / pi) * (numpy.log(n) + GAMMA + 1 / (2 * n))
                + (1 - 2 / pi) * (q - 2)
                + 0.727
            )
            / q
        )
    gap = (q - 1) * q ** (-b) - pb
    ok = gap > 0
    i = numpy.nonzero(ok)[0]
    first = int(q[i[0]])
    return (
        first,
        int(ok.sum()),
        bool(ok[i[0] :].all()),
        float(gap[i[0] :].min()),
        int(q[i[0] + int(numpy.argmin(gap[i[0] :]))]),
    )


def wall():
    print(
        "THE REGION A WALL: the least base with PB < (base-1) base^(-4/5) at one excluded digit"
    )
    print(
        "form                          every e_0   e_0 in {0, base-1}   scan          held            up-set      least float gap above the wall"
    )
    hi = 200000
    for name, chord, lo in [
        ("step 3, PB_base(1)", False, 17),
        ("chord, PB'_base(1, e_0)", True, 36),
    ]:
        rows = [scan(chord, w, lo, hi) for w in (True, False)]
        print(
            "%-28s %10d %20d   %d..%d   %d %d   %s %s   %.3e at %d, %.3e at %d"
            % (
                name,
                rows[0][0],
                rows[1][0],
                lo,
                hi,
                rows[0][1],
                rows[1][1],
                rows[0][2],
                rows[1][2],
                rows[0][3],
                rows[0][4],
                rows[1][3],
                rows[1][4],
            )
        )
    print(
        "the worst e_0 is the middle digit at odd base and c = 1/2 at even base, the best e_0 in {0, base-1}"
    )
    print()
    print("THE MARGIN AT EACH WALL, 40 digits, gap = (base-1) base^(-4/5) - PB")
    for name, chord, q, worst in [
        ("step 3", False, 92317, False),
        ("chord, worst e_0", True, 39363, True),
        ("chord, e_0 = 0", True, 28352, False),
    ]:
        g0 = gap_hp(q - 1, (q - 1) // 2 if worst else 0, 0.8, chord)
        g1 = gap_hp(q, q // 2 if worst else 0, 0.8, chord)
        print(
            "%-18s base %6d gap %+.6e   base %6d gap %+.6e"
            % (name, q - 1, float(g0), q, float(g1))
        )
    print()
    print(
        "REGION A ASKS alpha_1 < 1/5, at each wall: 1/5 - alpha_1 = alpha - c - 4/5 = log(1 + gap/PB)/log base > 0, truncated down"
    )
    import mpmath

    for name, chord, q, worst in [
        ("step 3", False, 92317, False),
        ("chord, worst e_0", True, 39363, True),
        ("chord, e_0 = 0", True, 28352, False),
    ]:
        e0 = q // 2 if worst else 0
        g = gap_hp(q, e0, 0.8, chord)
        pb = (q - 1) * mpmath.mpf(q) ** (-mpmath.mpf(0.8)) - g
        v = mpmath.log(1 + g / pb) / mpmath.log(q)
        e = int(mpmath.floor(mpmath.log10(v)))
        mant = mpmath.floor(v / mpmath.mpf(10) ** e * 10**4) / 10**4
        print("%-18s base %6d  1/5 - alpha_1 >= %.4fe%03d" % (name, q, float(mant), e))

    print()
    print(
        "THE UNSHIFTED MASS region A pays: c_k = sum_(a mod base^k) |hat F_k(a/base^k)|, readings of c_k/c_(k-1) against 2(base - 1)"
    )
    import numpy

    for q, e0, k in [(33, 16, 4), (33, 0, 4), (17, 8, 5)]:
        F = [v for v in range(q) if v != e0]
        cs = []
        for kk in (k - 1, k):
            ind = numpy.zeros(q**kk)
            ind[strings(q, F, kk)] = 1.0
            cs.append(numpy.abs(numpy.fft.fft(ind)).sum())
        print(
            "  base %d, excluded %2d, k = %d: c_k/c_(k-1) = %.4f against 2(base - 1) = %d"
            % (q, e0, k, cs[1] / cs[0], 2 * (q - 1))
        )


# GRID


def mu_upto(N):
    import numpy

    mu = numpy.ones(N + 1, dtype=numpy.int8)
    mu[0] = 0
    sieve = numpy.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    primes = numpy.nonzero(sieve)[0]
    for p in primes:
        mu[p::p] = -mu[p::p]
    for p in primes:
        if p * p > N:
            break
        mu[p * p :: p * p] = 0
    return mu


def strings(q, F, L):
    import numpy

    out = numpy.zeros(1, dtype=numpy.int64)
    d = numpy.array(F, dtype=numpy.int64)
    for i in range(L):
        out = (out[:, None] + d[None, :] * (q**i)).ravel()
    out.sort()
    return out


def convergent(num, den, cap):
    import numpy

    n = num.astype(numpy.int64).copy()
    d = (
        numpy.full_like(n, den)
        if numpy.isscalar(den)
        else den.astype(numpy.int64).copy()
    )
    cap = numpy.broadcast_to(numpy.asarray(cap, dtype=numpy.int64), n.shape)
    h1 = numpy.ones_like(n)
    h2 = numpy.zeros_like(n)
    k1 = numpy.zeros_like(n)
    k2 = numpy.ones_like(n)
    live = d > 0
    c = numpy.zeros_like(n)
    c[live] = n[live] // d[live]
    h1, h2 = c * h1 + h2, h1
    k1, k2 = c * k1 + k2, k1
    n, d = d, n - c * d
    live = d > 0
    while live.any():
        c = numpy.zeros_like(n)
        c[live] = n[live] // d[live]
        hn = c * h1 + h2
        kn = c * k1 + k2
        go = live & (kn <= cap)
        h2 = numpy.where(go, h1, h2)
        k2 = numpy.where(go, k1, k2)
        h1 = numpy.where(go, hn, h1)
        k1 = numpy.where(go, kn, k1)
        nn = numpy.where(go, d, n)
        dd = numpy.where(go, n - c * d, 0)
        n, d = nn, dd
        live = go & (d > 0)
    return h1, k1


def smooth(d, q):
    from math import gcd

    if d < 1:
        return False
    while True:
        g = gcd(d, q)
        if g == 1:
            return d == 1
        while d % g == 0:
            d //= g


GAL_CAP = 10**5


def gmod(F, th):
    import numpy

    z = numpy.zeros(th.shape, dtype=complex)
    for v in F:
        z += numpy.exp(2j * numpy.pi * v * th)
    return numpy.abs(z)


def l1_norms(q, F, i1):
    import numpy

    V = q**i1
    t = numpy.arange(4 * V) / (4 * V)
    gs = [
        sum(numpy.exp(2j * numpy.pi * v * (q**j * t % 1.0)) for v in F)
        for j in range(i1)
    ]
    dg = [
        sum(
            2j * numpy.pi * v * numpy.exp(2j * numpy.pi * v * (q**j * t % 1.0))
            for v in F
        )
        for j in range(i1)
    ]
    f = numpy.ones(t.size, dtype=complex)
    for g in gs:
        f = f * g
    fp = numpy.zeros(t.size, dtype=complex)
    for j in range(i1):
        term = q**j * dg[j]
        for jj in range(i1):
            if jj != j:
                term = term * gs[jj]
        fp = fp + term
    return float(numpy.abs(f).mean()), float(numpy.abs(fp).mean())


def regions_one(q, e0, k, Z):
    from math import log, pi

    import numpy

    F = [v for v in range(q) if v != e0]
    fill = len(F)
    y = q**k
    D = numpy.zeros(y)
    D[strings(q, F, k)] = 1.0
    mu = mu_upto(y)[:y].astype(numpy.float64)
    hatF = numpy.conj(numpy.fft.fft(D))
    Sneg = numpy.fft.fft(mu)
    exact = float(mu[D > 0].sum())
    a = numpy.arange(y, dtype=numpy.int64)
    Q = int(y**0.6)
    l, d = convergent(a, y, Q)
    h = numpy.abs(a * d - l * y)
    assert (d <= Q).all() and (h * Q <= y).all()
    assert (numpy.gcd(l, d) == 1).all()
    y25 = y**0.4
    sm = numpy.array([smooth(int(v), q) for v in range(Q + 1)])
    A = d >= y25
    C = (~A) & (d < Z) & (h < Z)
    B = (~A) & (~C)
    C2 = C & sm[d]
    C1 = C & (~sm[d])
    term = hatF * Sneg / y
    parts = {
        nm: term[msk].sum().real
        for nm, msk in [("A", A), ("B", B), ("C1", C1), ("C2", C2)]
    }
    tot = sum(parts.values())
    print(
        "base %d, excluded %d, level %d, y = %d, Q = y^(3/5) = %d, Z = %d"
        % (q, e0, k, y, Q, Z)
    )
    print(
        "  exact sum of mu over the strings %d, the four regions sum to %.6f, difference %.2e"
        % (exact, tot, abs(tot - exact))
    )
    assert abs(tot - exact) < 1e-6 * y
    for nm in ("A", "B", "C1", "C2"):
        msk = {"A": A, "B": B, "C1": C1, "C2": C2}[nm]
        print(
            "  region %-2s  points %8d   contribution %+14.4f   = %+.6f fill^k"
            % (nm, msk.sum(), parts[nm], parts[nm] / fill**k)
        )
    assert C1.sum() <= 3 * Z * Z
    nZ = sum(1 for v in range(1, Z) if sm[v])
    assert C2.sum() <= 3 * Z * nZ
    print(
        "  counts: C1 %d <= 3 Z^2 = %d, C2 %d <= 3 Z N_Z = %d"
        % (C1.sum(), 3 * Z * Z, C2.sum(), 3 * Z * nZ)
    )
    dd = d[C2]
    assert (y % dd == 0).all()
    j = a[C2] - l[C2] * (y // dd)
    j = (j + y // 2) % y - y // 2
    assert (numpy.abs(j) * dd == h[C2]).all() and (numpy.abs(j) * dd < Z).all()
    print(
        "  C2 algebra: every d divides y and a = l y/d + j with |j| d = h < Z at all %d points"
        % C2.sum()
    )
    ib = numpy.nonzero(B & (h >= 1))[0]
    capq = (2 * y) // h[ib]
    l2, d2 = convergent(a[ib], y, capq)
    lo = d2 * 2 * h[ib] >= y
    hi = d2 * h[ib] <= 2 * y
    diff = l2 * d[ib] != l[ib] * d2
    assert lo.all() and hi.all() and diff.all()
    print(
        "  second approximation: y/(2h) <= d' <= 2y/h and l'/d' != l/d at all %d points of B with h >= 1"
        % ib.size
    )
    B3 = q * pb_step3(q, 1)
    a1 = log(B3 / fill) / log(q)
    af = numpy.abs(hatF)
    Dc = numpy.floor(numpy.log2(d)).astype(numpy.int64)
    Hc = numpy.where(
        h > 0, numpy.floor(numpy.log2(numpy.maximum(h, 1))).astype(numpy.int64) + 1, 0
    )
    key = Dc * 64 + Hc
    order = numpy.argsort(key, kind="stable")
    ks, starts = numpy.unique(key[order], return_index=True)
    sums = numpy.add.reduceat(af[order], starts)
    worst, tested, skipped = 0.0, 0, 0
    gal, top_ok, norms, galn = 0.0, True, {}, 0
    for kk, s, st in zip(ks, sums, starts):
        Dv = 2 ** int(kk // 64)
        hc = int(kk % 64)
        Hv = 0 if hc == 0 else 2 ** (hc - 1)
        V1 = 1
        while V1 < 4 * Dv * Dv:
            V1 *= q
        V2 = 1
        while V2 < 4 * Hv / Dv + 1:
            V2 *= q
        if V1 * V2 > y or 16 * Dv * Hv > y:
            skipped += 1
            continue
        tested += 1
        if Dv == 1 and Hv == 0:
            continue
        bound = (1 + pi * q) * fill**k * (V1 * V2) ** a1
        worst = max(worst, s / bound)
        mem = order[
            st : st
            + (
                len(order)
                if st == starts[-1]
                else starts[numpy.searchsorted(starts, st) + 1] - st
            )
        ]
        dm = d[mem]
        fk = (l[mem] % dm) * (Q + 1) + dm
        i1 = round(log(V1) / log(q))
        bot = numpy.ones(mem.size)
        for j in range(i1):
            bot = bot * gmod(F, (q**j * a[mem] % y) / y)
        o2 = numpy.argsort(fk, kind="stable")
        _, s2 = numpy.unique(fk[o2], return_index=True)
        lhs = numpy.maximum.reduceat(bot[o2], s2).sum()
        if V1 <= GAL_CAP:
            if i1 not in norms:
                norms[i1] = l1_norms(q, F, i1)
            rhs = 4 * Dv * Dv * norms[i1][0] + norms[i1][1]
            gal = max(gal, lhs / rhs)
            galn += 1
        grp = numpy.split(a[mem][o2] % V2, s2[1:])
        top_ok = top_ok and all(numpy.unique(g).size == g.size for g in grp)
    assert worst <= 1.0 and gal <= 1.0 and top_ok
    print(
        "  hybrid lemma, alpha_1 <= %.6f from the step 3 constant: %d classes (D, H) meet V_1 V_2 <= y and 16 D H <= y, %d do not"
        % (a1 + 5e-7, tested, skipped)
    )
    print(
        "  past the class of a = 0, where |hat F_k(0)| = fill^k: largest ratio of a class l^1 sum to the lemma %.6f;"
        % (worst + 5e-7)
    )
    print(
        "  the Gallagher step at the %d classes with V_1 <= %d: the fractions' largest |hat F_(i_1)| at their residues, summed,"
        % (galn, GAL_CAP)
    )
    print(
        "  against 4 D^2 ||f||_1 + ||f'||_1 with both norms read on 4 V_1 points, largest ratio %.6f; every fraction's residues"
        " distinct mod V_2" % (gal + 5e-7)
    )
    L = log(y)
    rb = numpy.abs(Sneg[ib]) / (y**0.8 + y * L**3 * numpy.maximum(d[ib], h[ib]) ** -0.5)
    ra = numpy.abs(Sneg[A]) / y**0.8
    print(
        "  readings, not bounds: max |S|/(y^(4/5) + y (log y)^3 max(d,h)^(-1/2)) on B %.3e, max |S|/y^(4/5) on A %.4f"
        % (rb.max(), ra.max())
    )
    print()


def lemma_a_prime(q, e0, k, Z, draws, seed):
    import random
    from math import cos, floor, gcd, log, pi, sin

    F = [v for v in range(q) if v != e0]
    fill = len(F)
    y = q**k
    rng = random.Random(seed)
    cp = 1 - (2 / fill) * (1 - cos(pi / (4 * q)))
    worst, n = -1e9, 0
    while n < draws:
        d = rng.randrange(2, Z)
        if smooth(d, q):
            continue
        l = rng.randrange(d)
        if gcd(l, d) != 1:
            continue
        j = rng.randrange(-Z, Z)
        a = (l * y) // d + j
        h = abs(a * d - l * y)
        if h >= Z:
            continue
        assert (h / (d * y)) < y ** (-2 / 3) / (4 * q * (q - 1))
        lf = 0.0
        for i in range(k):
            ph = (pow(q, i, y) * a % y) / y
            s = abs(sum(complex(cos(2 * pi * v * ph), sin(2 * pi * v * ph)) for v in F))
            lf += log(s / fill)
        md = max(1, floor(log(d / 2) / log(q) + 1e-12) + 1)
        rhs = floor(2 * k / (3 * md)) * log(cp)
        worst = max(worst, lf - rhs)
        n += 1
    assert worst <= 1e-9
    print(
        "perturbed Lemma A' at the grid point, base %d, excluded %d, level %d, d < %d with a prime outside the base, h < %d"
        % (q, e0, k, Z, Z)
    )
    print(
        "  %d seeded draws: log |hat F_k(a/y)|/fill^k minus floor(2k/(3 m_d)) log c' is at most %.4f, never above 0"
        % (draws, worst)
    )


def regions():
    t = time.time()
    print(
        "FALSIFICATION OF THE DISSECTION at small base and level, every grid point a mod y"
    )
    print(
        "A: d >= y^(2/5); B: d < y^(2/5) and max(d, h) >= Z; C1, C2: d < Z and h < Z, d with or without a prime outside the base"
    )
    print()
    regions_one(10, 5, 6, 16)
    regions_one(10, 0, 6, 16)
    regions_one(5, 2, 9, 12)
    lemma_a_prime(10, 5, 30, 60, 3000, 1009)
    lemma_a_prime(5, 0, 45, 40, 3000, 1013)
    print("runtime %.1f s" % (time.time() - t))


# BLOCKS


def blocks_of(x, q, F):
    ds = []
    v = x
    while v:
        ds.append(v % q)
        v //= q
    L = len(ds)
    out = []
    if 0 in F:
        if L >= 2:
            out.append((0, L - 1))
    else:
        for ell in range(1, L):
            out.append((0, ell))
    for j in range(L - 1, -1, -1):
        top = 0
        for i in range(L - 1, j, -1):
            top = top * q + ds[i]
        for f in F:
            if f < ds[j] and not (j == L - 1 and f == 0):
                out.append((top * q + f, j))
        if ds[j] not in F:
            break
    tail = all(v in F for v in ds)
    return out, tail, L


def blocks():
    import random

    import numpy

    t = time.time()
    print(
        "OFF THE POWERS OF THE BASE: S_F below x as blocks P base^k + D_k, at most fill + 1 per scale"
    )
    N = 2 * 10**6
    mu = mu_upto(N)
    rng = random.Random(7)
    for q, e0 in [(10, 0), (10, 9), (5, 0), (7, 3), (4, 1)]:
        F = [v for v in range(q) if v != e0]
        fill = len(F)
        memb = numpy.zeros(N + 1, dtype=bool)
        L = 1
        while q**L <= N:
            L += 1
        for ell in range(1, L + 1):
            s = strings(q, F, ell)
            s = s[s >= q ** (ell - 1)] if 0 in F else s
            memb[s[s <= N]] = True
        memb[0] = False
        cum = numpy.cumsum(numpy.where(memb, mu, 0).astype(numpy.int64))
        cnt = numpy.cumsum(memb.astype(numpy.int64))
        worst = 0
        for x in [rng.randrange(1, N) for _ in range(400)] + [
            q**e - 1 for e in range(1, L) if q**e <= N
        ]:
            bl, tail, Lx = blocks_of(x, q, F)
            per = {}
            tot = 0
            for P, kk in bl:
                per[kk] = per.get(kk, 0) + 1
                s = strings(q, F, kk) + P * q**kk
                tot += int(mu[s[s > 0]].sum())
            if tail:
                tot += int(mu[x])
            assert tot == cum[x], (q, e0, x)
            assert max(per.values(), default=0) <= fill + 1
            assert cnt[x] >= fill ** (Lx - 1) - 1
            worst = max(worst, max(per.values(), default=0))
        print(
            "  base %2d excluded %d: 400 random x below %d and every base^e - 1, the block sums meet M_F(x) exactly;"
            % (q, e0, N)
        )
        print(
            "           at most %d blocks at one scale against fill + 1 = %d, and A_F(x) >= fill^(L-1) - 1 at every x"
            % (worst, fill + 1)
        )
    print("runtime %.1f s" % (time.time() - t))


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "wall"
    {"wall": wall, "regions": regions, "blocks": blocks}[verb]()
