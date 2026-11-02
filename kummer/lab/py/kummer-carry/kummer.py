import math
import sys
import time
from fractions import Fraction
from functools import reduce
from itertools import pairwise

import numpy as np
from mpmath import iv
from sympy import isprime, nextprime, prevprime

# DIGITS


def digsum(n, p):
    n = np.array(n, dtype=np.int64)
    s = np.zeros_like(n)
    while (n > 0).any():
        s += n % p
        n //= p
    return s


def run(k, p, digits, start):
    c = np.full(np.shape(k), start, dtype=np.int64)
    ok = np.ones(np.shape(k), dtype=bool)
    kk = np.array(k, dtype=np.int64)
    for _ in range(digits):
        d = kk % p
        a0 = (3 * d + c <= p - 1) & (2 * d + c <= p - 1)
        a1 = (2 * d + c >= p) & (3 * d + c <= 2 * p - 1)
        ok &= a0 | a1
        c = np.where(a1, 1, 0)
        kk //= p
    return ok, c


def intervals(p):
    return {
        (0, 0): (0, (p - 1) // 3),
        (0, 1): ((p + 1) // 2, (2 * p - 1) // 3),
        (1, 0): (0, (p - 2) // 3),
        (1, 1): ((p - 1) // 2, (2 * p - 2) // 3),
    }


def table(p, a, b):
    out = {}
    for c1 in range(b):
        for c2 in range(a - b):
            row = []
            for d in range(p):
                x, y = b * d + c1, (a - b) * d + c2
                if x % p + y % p <= p - 1:
                    row.append((d, (x // p, y // p)))
            out[(c1, c2)] = row
    return out


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for q in range(2, int(n**0.5) + 1):
        if s[q]:
            s[q * q :: q] = False
    return np.nonzero(s)[0]


# KUMMER


def kummer():
    t0 = time.time()
    bad_sets = 0
    for p in [q for q in range(5, 3000) if isprime(q)]:
        T = table(p, 3, 1)
        I = intervals(p)
        for c in (0, 1):
            for c2 in (0, 1):
                ds = [d for d, s in T[(0, c)] if s == (0, c2)]
                lo, hi = I[(c, c2)]
                if ds != list(range(lo, hi + 1)):
                    bad_sets += 1
            if (
                len(T[(0, c)]) != (p + 1) // 2
                or (0, (0, 0)) not in T[(0, c)]
                or (1, (0, 0)) not in T[(0, c)]
            ):
                bad_sets += 1
    print(
        f"digit sets: the four intervals, row sums (p+1)/2 and the pair {{0,1}} into state 0, every prime 5..2999, mismatches {bad_sets}"
    )
    assert bad_sets == 0
    wit = 0
    ps = [q for q in range(5, 3000) if isprime(q)]
    for p in ps:
        u, v = (p + 1) // 2, (p - 1) // 3 if p % 6 == 1 else (2 * p - 1) // 3
        car = lambda k: (
            int(digsum(k, p) + digsum(2 * k, p) - digsum(3 * k, p)) // (p - 1)
        )
        if car(u) == 0 and car(v * p) == 0 and car(u + v * p) > 0:
            wit += 1
    print(
        f"not a digit design: u = (p+1)/2 and v p in K, u + v p outside K, with v = (p-1)/3 or (2p-1)/3, at {wit} of the {len(ps)} primes 5..2999"
    )
    assert wit == len(ps)
    for p in (5, 7, 11, 13):
        k = np.arange(p**5, dtype=np.int64)
        num = digsum(k, p) + digsum(2 * k, p) - digsum(3 * k, p)
        assert (num % (p - 1) == 0).all()
        ok, _ = run(k, p, 5, 0)
        bad = int(((num == 0) != ok).sum())
        print(
            f"p {p}: every k < p^5 ({p**5}), automaton against (s(k)+s(2k)-s(3k))/(p-1) = 0, mismatches {bad}, |K below p^5| {int(ok.sum())} = ((p+1)/2)^5 {((p + 1) // 2) ** 5}"
        )
        assert bad == 0 and int(ok.sum()) == ((p + 1) // 2) ** 5
    p = 101
    k = np.arange(p**3, dtype=np.int64)
    num = digsum(k, p) + digsum(2 * k, p) - digsum(3 * k, p)
    ok, _ = run(k, p, 3, 0)
    bad = int(((num == 0) != ok).sum())
    print(f"p 101: every k < p^3 directly, mismatches {bad}")
    assert bad == 0
    P2 = p * p
    r = np.arange(P2, dtype=np.int64)
    R = digsum(r, p) + digsum((2 * r) % P2, p) - digsum((3 * r) % P2, p)
    c2r, c3r = 2 * r // P2, 3 * r // P2
    okr, sr = run(r, p, 2, 0)
    q = np.arange(p**3, dtype=np.int64)
    sq = digsum(q, p)
    Q = {
        (u, v): sq + digsum(2 * q + u, p) - digsum(3 * q + v, p)
        for u in (0, 1)
        for v in (0, 1, 2)
    }
    okq = {s: run(q, p, 3, s)[0] for s in (0, 1)}
    pairs = bad = count = 0
    for key in set(zip(c2r.tolist(), c3r.tolist(), sr.tolist())):
        m = (c2r == key[0]) & (c3r == key[1]) & (sr == key[2])
        rv = set(zip(R[m].tolist(), okr[m].tolist()))
        qv = set(zip(Q[(key[0], key[1])].tolist(), okq[key[2]].tolist()))
        for x, ox in rv:
            for y, oy in qv:
                pairs += 1
                if (x + y) % (p - 1) != 0 or ((x + y) == 0) != (ox and oy):
                    bad += 1
        count += int(okr[m].sum()) * int(okq[key[2]].sum())
    print(
        f"p 101: every k < p^5 ({p**5}) through k = q p^2 + r, {pairs} value pairs over the carry classes, mismatches {bad}, |K below p^5| {count} = 51^5 {51**5}"
    )
    assert bad == 0 and count == 51**5
    print(f"kummer: {time.time() - t0:.1f} s")


# FAMILY


def limit_design(a, b):
    pts = sorted({Fraction(j, m) for m in (a, b, a - b) for j in range(m + 1)})
    iv_ = []
    for lo, hi in pairwise(pts):
        th = (lo + hi) / 2
        if (b * th) % 1 + ((a - b) * th) % 1 < 1:
            if iv_ and iv_[-1][1] == lo:
                iv_[-1] = (iv_[-1][0], hi)
            else:
                iv_.append((lo, hi))
    return iv_


def nu_limit(a, b, t=0.0):
    I = limit_design(a, b)
    per = reduce(math.lcm, [x.denominator for pair in I for x in pair])
    xs = np.arange(per) + t
    num = sum(
        np.exp(2j * np.pi * float(hi) * xs) - np.exp(2j * np.pi * float(lo) * xs)
        for lo, hi in I
    )
    return float(np.abs(num).mean()), I, per


def family():
    t0 = time.time()
    bad = 0
    fails = []
    for a in range(2, 8):
        for b in range(1, a):
            for p in [q for q in range(2, 38) if isprime(q)]:
                L = 2
                while p ** (L + 1) <= 200000:
                    L += 1
                k = np.arange(p**L, dtype=np.int64)
                v = digsum(b * k, p) + digsum((a - b) * k, p) - digsum(a * k, p)
                n = int((v == 0).sum())
                if p > a and n != ((p + 1) // 2) ** L:
                    bad += 1
                if p <= a and n != ((p + 1) // 2) ** L and p % 2 == 1:
                    fails.append((a, b, p, L, n, ((p + 1) // 2) ** L))
    print(
        f"count: |{{k < p^L : p does not divide C(ak,bk)}}| = ((p+1)/2)^L at every 2 <= a <= 7, 1 <= b < a, prime a < p <= 37, p^L <= 200000, mismatches {bad}"
    )
    assert bad == 0
    print(f"count at odd p <= a, first failures: {fails[:4]}")
    rows = bad_rows = 0
    pair_bad = 0
    least_fail = {}
    for a in range(2, 13):
        for b in range(1, a):
            for p in [q for q in range(a + 1, 400) if isprime(q)]:
                T = table(p, a, b)
                for st, row in T.items():
                    rows += 1
                    if len(row) != (p + 1) // 2:
                        bad_rows += 1
                    has = (0, (0, 0)) in row and (1, (0, 0)) in row
                    if p >= 2 * a - 1 and not has:
                        pair_bad += 1
                    if p < 2 * a - 1 and not has and (a, b) not in least_fail:
                        least_fail[(a, b)] = (p, st)
    print(
        f"rows: every carry state of every (a, b), 2 <= a <= 12, prime a < p < 400, {rows} rows, rows off (p+1)/2: {bad_rows}, rows missing the pair {{0,1}} into (0,0) at p >= 2a-1: {pair_bad}"
    )
    assert bad_rows == 0 and pair_bad == 0
    print(f"pair lost below 2a-1, first witnesses: {sorted(least_fail.items())[:5]}")
    worst = {}
    p = 10007
    for a in range(2, 8):
        for b in range(1, a):
            T = table(p, a, b)
            I = limit_design(a, b)
            lim = {d for d in range(p) if any(lo * p <= d < hi * p for lo, hi in I)}
            row0 = {d for d, s in T[(0, 0)]}
            dev = len(lim ^ row0)
            for st, row in T.items():
                lab = {(d, s) for d, s in row}
                lab0 = {(d, s) for d, s in T[(0, 0)]}
                dev = max(dev, len({d for d, s in lab ^ lab0}))
            worst[(a, b)] = dev
    print(
        f"p 10007: digits by which row (0,0) misses the limit design, or another row misses row (0,0) with its moves: max {max(worst.values())} over every (a, b) with a <= 7, {worst}"
    )
    for a in range(2, 8):
        for b in range(1, a // 2 + 1):
            nu, I, per = nu_limit(a, b)
            assert sum(hi - lo for lo, hi in I) == Fraction(1, 2)
            assert abs(nu - (3 - 4 * math.gcd(a, b) / a)) < 1e-12
            print(
                f"(a,b) = ({a},{b}): row-0 limit design {len(I)} intervals {[(str(lo), str(hi)) for lo, hi in I]}, period {per}, nu(0) = {nu:.6f}, (2/pi) nu = {2 * nu / math.pi:.6f}"
            )
    prod = []
    cells = 0
    for a in range(2, 7):
        for b in range(1, a):
            for p in [q for q in range(a + 1, 32) if isprime(q)]:
                k = np.arange(p**3, dtype=np.int64)
                ok = (digsum(b * k, p) + digsum((a - b) * k, p) - digsum(a * k, p)) == 0
                kk = k[ok]
                size = 1
                for j in range(3):
                    size *= len(set(((kk // p**j) % p).tolist()))
                cells += 1
                if size == len(kk):
                    prod.append((a, b, p))
    print(
        f"product sets below p^3 over every 2 <= a <= 6, 1 <= b < a, prime a < p <= 31 ({cells} cells): {sorted({(a, b) for a, b, p in prod})}, at {len(prod)} cells"
    )
    print(f"family: {time.time() - t0:.1f} s")


# CONSTANTS

E = lambda x: np.exp(2j * np.pi * x)
AF = lambda z: np.abs(E(z / 3) - 1 + E(2 * z / 3) - E(z / 2))
A1 = lambda z: np.abs(E(z / 6) - 1)
W = lambda z: 1 / z - np.log(z) / 6


def nu_h(a, t):
    nu = sum(a(r + t) for r in range(6)) / 6
    h = sum(a(r + t) * W(r + t) for r in range(1, 7)) + sum(
        a(r + 1 - t) * W(r + 1 - t) for r in range(6)
    )
    return nu, h


NUM = {"": ((1, 1, 3), (-1, 0, 1), (1, 2, 3), (-1, 1, 2)), "1": ((1, 1, 6), (-1, 0, 1))}


def w_iv(z):
    return 1 / z - iv.log(z) / 6


def a_iv(terms, z):
    re, im = iv.mpf(0), iv.mpf(0)
    for c, n, d in terms:
        x = 2 * iv.pi * n * z / d
        re += c * iv.cos(x)
        im += c * iv.sin(x)
    s = re**2 + im**2
    return iv.sqrt(iv.mpf([max(s.a, 0), s.b]))


def lip_sums():
    iv.prec = 120
    S, D = iv.mpf(0), iv.mpf(0)
    for m in range(1, 13):
        lo, hi = abs(w_iv(iv.mpf(m) / 2)), abs(w_iv(iv.mpf(m + 1) / 2))
        S += iv.mpf([max(lo.a, hi.a), max(lo.b, hi.b)])
        D += iv.mpf(4) / m**2 + iv.mpf(1) / (3 * m)
    return S, D


def lip_iv(la, amax):
    S, D = lip_sums()
    return la * S + amax * D


def up(x, d=6):
    return math.ceil(x * 10**d) / 10**d


def constants(verbose=True):
    N = 200000
    t = np.linspace(0, 0.5, N + 1)
    h = 0.5 / N
    out = {}
    for name, a, la, amax in (("", AF, 3 * iv.pi, 4), ("1", A1, iv.pi / 3, 2)):
        nu, H = nu_h(a, t)
        lh = float(lip_iv(la, amax).b)
        out["nu" + name] = up(nu.max() + float(la.b) * h / 2 + 1e-12)
        out["H" + name] = up(H.max() + lh * h / 2 + 1e-12)
        if verbose:
            print(
                f"nu{name}: grid max {nu.max():.9f} at t = {t[nu.argmax()]:.6f}, bound {out['nu' + name]}; H{name}: grid max {H.max():.9f} at t = {t[H.argmax()]:.6f}, Lipschitz <= {up(lh, 3)}, bound {out['H' + name]}; nu{name}(0) = {nu[0]:.9f}"
            )
    return out


def cover(f, lip, bound):
    cells, deep, todo = 0, 1, [(0, 1)]
    while todo:
        k, d = todo.pop()
        assert d < 60
        if (f(iv.mpf(2 * k + 1) / 2 ** (d + 1)) + lip / 2 ** (d + 1)).b <= bound.a:
            cells, deep = cells + 1, max(deep, d)
        else:
            todo += [(2 * k, d + 1), (2 * k + 1, d + 1)]
    return cells, deep + 1


def certify(K):
    S, D = lip_sums()
    Sc = iv.mpf(191944) / 45045 + iv.log(iv.mpf(858) / 7) / 6
    Dq = sum(Fraction(4, m * m) + Fraction(1, 3 * m) for m in range(1, 13))
    Di = iv.mpf(Dq.numerator) / Dq.denominator
    assert (
        S.a <= Sc.b
        and Sc.a <= S.b
        and D.a <= Di.b
        and Di.a <= D.b
        and Dq == Fraction(280246811, 38419920)
    )
    w4, w92 = w_iv(iv.mpf(4)), w_iv(iv.mpf(9) / 2)
    assert w4.a > 0 > w92.b and -w92.b > w4.b
    LH, LH1 = 3 * iv.pi * S + 4 * D, iv.pi / 3 * S + 2 * D
    print(
        f"Lipschitz: w(4) = {dn(w4.a, 5)}.. > 0 > w(9/2) = -{dn(-w92.b, 5)}.., S = 191944/45045 + log(858/7)/6 = {dn(S.a, 8)}.., S' = {Dq} = {dn(D.a, 8)}..; L_H = 3 pi S + 4 S' = {dn(LH.a, 7)}.. <= {up(float(LH.b), 3)}, L_H1 = (pi/3) S + 2 S' = {dn(LH1.a, 7)}.. <= {up(float(LH1.b), 3)}"
    )
    assert ((iv.sqrt(6) + iv.sqrt(2)) / 3).b <= iv.mpf(str(K["nu1"])).a
    for name, la, amax in (("", 3 * iv.pi, 4), ("1", iv.pi / 3, 2)):
        a = lambda z, terms=NUM[name]: a_iv(terms, z)
        nu = lambda t, a=a: sum((a(r + t) for r in range(6)), iv.mpf(0)) / 6
        H = lambda t, a=a: (
            sum((a(r + t) * w_iv(r + t) for r in range(1, 7)), iv.mpf(0))
            + sum((a(r + 1 - t) * w_iv(r + 1 - t) for r in range(6)), iv.mpf(0))
        )
        jobs = (
            (("H" + name, H, lip_iv(la, amax)),)
            if name
            else (("nu", nu, la), ("H", H, lip_iv(la, amax)))
        )
        for key, f, lip in jobs:
            cells, deep = cover(f, lip, iv.mpf(str(K[key])))
            print(
                f"{key} <= {K[key]} on [0, 1/2]: Lipschitz <= {up(float(lip.b), 3)}, interval cover closes in {cells} dyadic cells, finest half-width 2^-{deep}"
            )
    print(
        f"nu1 <= {K['nu1']} on [0, 1/2]: closed form, peak (sqrt 6 + sqrt 2)/3 at t = 1/2"
    )


# CERTIFICATE


def cert_iv(p, K):
    iv.prec = 120
    P = iv.mpf(p)
    pi = iv.pi
    L = iv.log(P / 2)
    fill = (P + 1) / 2
    nu, H, nu1, H1 = (iv.mpf(str(K[x])) for x in ("nu", "H", "nu1", "H1"))
    n1 = (P - 1) / 6 if p % 6 == 1 else (P + 1) / 6
    A = fill + (P / (2 * pi)) * (2 * nu * L + H) + 3 * (P - 1) / 2 + (1 - 2 / pi) * P
    B = n1 + (P / (2 * pi)) * (2 * nu1 * L + H1) + (P - 1) / 6 + (1 - 2 / pi) * P / 2
    C = 2 / iv.sin(pi / (2 * P))
    D = P if p % 6 == 1 else C
    lam = (A + D) / 2 + iv.sqrt(((A - D) / 2) ** 2 + B * C)
    return lam, fill, (lam - A) / C


def cert_float(p, K):
    p = np.asarray(p, dtype=float)
    L = np.log(p / 2)
    fill = (p + 1) / 2
    cls1 = np.round(p) % 6 == 1
    n1 = np.where(cls1, (p - 1) / 6, (p + 1) / 6)
    A = (
        fill
        + (p / (2 * np.pi)) * (2 * K["nu"] * L + K["H"])
        + 1.5 * (p - 1)
        + (1 - 2 / np.pi) * p
    )
    B = (
        n1
        + (p / (2 * np.pi)) * (2 * K["nu1"] * L + K["H1"])
        + (p - 1) / 6
        + (1 - 2 / np.pi) * p / 2
    )
    C = 2 / np.sin(np.pi / (2 * p))
    D = np.where(cls1, p, C)
    lam = (A + D) / 2 + np.sqrt(((A - D) / 2) ** 2 + B * C)
    return lam / fill


def c_inf_iv(K):
    iv.prec = 120
    pi = iv.pi
    nu, H, nu1, H1 = (iv.mpf(str(K[x])) for x in ("nu", "H", "nu1", "H1"))
    beta1 = 2 * pi * (iv.mpf(1) / 3 + (1 - 2 / pi) / 2)
    beta0 = 2 * pi * (3 - 6 / pi) - 2 * pi * (1 + 2 * (1 - 2 / pi) / 13) / 13
    assert ((H1 + beta1) * nu - (H + beta0) * nu1).b < 0
    return 4 + 2 * (1 - 2 / pi) + (H - 2 * nu * iv.log(2)) / pi + (
        8 / pi
    ) * nu1 / nu, 2 * nu / pi


def tail_gap(n, e, K):
    c, al = c_inf_iv(K)
    N = iv.mpf(n)
    return N ** (iv.mpf(1) / e) - al * iv.log(N) - c


def sci(x, way, d=4):
    ex = math.floor(math.log10(abs(x)))
    m = x / 10**ex
    m = (math.floor(m * 10**d) if way == "dn" else math.ceil(m * 10**d)) / 10**d
    return f"{m:.{d}f} * 10^{ex}"


def dn(x, d=7):
    return math.floor(float(x) * 10**d) / 10**d


def wall_bar(e, K, lo, hi, primes):
    ps = primes[(primes >= lo) & (primes <= hi)]
    g = ps ** (1 / e) - cert_float(ps, K)
    fail = ps[g <= 0]
    pm = int(fail.max())
    p0 = int(ps[ps > pm].min())
    _c, al = c_inf_iv(K)
    x = p0
    while not tail_gap(x, e, K).a > 0:
        x = int(x * 1.001) + 1
    lo_, hi_ = p0, x
    while hi_ - lo_ > 1:
        mid = (lo_ + hi_) // 2
        if tail_gap(mid, e, K).a > 0:
            hi_ = mid
        else:
            lo_ = mid
    ptail = hi_
    assert al.b * e < ptail ** (1 / e)
    lam, fill, _ = cert_iv(pm, K)
    gm = iv.mpf(pm) ** (iv.mpf(1) / e) - lam / fill
    assert gm.b < 0
    seg = (
        [int(q) for q in primes[(primes >= p0) & (primes <= ptail)]]
        if ptail <= primes[-1]
        else None
    )
    assert seg is not None
    worst = None
    least = None
    for q in seg:
        lam, fill, _ = cert_iv(q, K)
        gq = iv.mpf(q) ** (iv.mpf(1) / e) - lam / fill
        assert gq.a > 0
        if worst is None or gq.a < worst[1]:
            worst = (q, float(gq.a))
        mq = iv.mpf(1) / e - iv.log(lam / fill) / iv.log(iv.mpf(q))
        if least is None or mq.a < least[1]:
            least = (q, float(mq.a))
    lam, fill, eta = cert_iv(p0, K)
    marg = iv.mpf(1) / e - iv.log(lam / fill) / iv.log(iv.mpf(p0))
    return (
        pm,
        p0,
        ptail,
        len(seg),
        float(gm.b),
        worst,
        float(marg.a),
        float(eta.b),
        least,
    )


def wall():
    t0 = time.time()
    K = constants()
    certify(K)
    c, al = c_inf_iv(K)
    print(
        f"tail: lambda/fill <= (2 nu/pi) log p + c_inf at every prime p >= 13, 2 nu/pi <= {up(float(al.b), 7)}, c_inf <= {up(float(c.b), 6)}"
    )
    primes = primes_upto(16000000)
    for e, lo, hi in ((5, 5000000, 16000000), (4, 20000, 2000000)):
        pm, p0, pt, n, gm, worst, marg, eta, least = wall_bar(e, K, lo, hi, primes)
        print(
            f"bar 1/{e}: last failing prime {pm} (gap <= {sci(gm, 'up')}), P_0 = {p0}, closed form certified at 120 bits at all {n} primes {p0}..{pt}, least gap >= {sci(worst[1], 'dn')} at {worst[0]}, tail bound clears from {pt} and increases; 1/{e} - alpha_1 >= {sci(marg, 'dn')} at P_0 and >= {sci(least[1], 'dn')} at {least[0]}, the least over the certified primes; eta <= {up(eta, 6)}"
        )
    print(
        "alpha_1 of the certificate, upper ends, first prime above 10^e in each class p mod 12:"
    )
    for ex in range(3, 10):
        row = []
        for cl in (1, 5, 7, 11):
            q = nextprime(10**ex - 1)
            while q % 12 != cl:
                q = nextprime(q)
            lam, fill, _ = cert_iv(q, K)
            a1 = iv.log(lam / fill) / iv.log(iv.mpf(q))
            row.append(f"{q} ({q % 3},{q % 4}) {up(float(a1.b), 7):.7f}")
        print("  " + " | ".join(row))
    for cl in (1, 5):
        for target in (0.22, 0.18):
            lo, hi = 10**5, 10**10
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if (
                    math.log(float(cert_float(np.array([mid - (mid % 6) + cl]), K)[0]))
                    / math.log(mid)
                    > target
                ):
                    lo = mid
                else:
                    hi = mid
            q = hi
            while q % 6 != cl or not isprime(q):
                q += 1
            lam, fill, _ = cert_iv(q, K)
            a1 = iv.log(lam / fill) / iv.log(iv.mpf(q))
            qb = prevprime(q)
            while qb % 6 != cl:
                qb = prevprime(qb)
            lamb, fillb, _ = cert_iv(qb, K)
            a1b = iv.log(lamb / fillb) / iv.log(iv.mpf(qb))
            assert a1.b < target < a1b.a
            print(
                f"class p = {cl} mod 6: alpha_1 crosses {target} between the primes {qb} ({dn(float(a1b.a), 10)}) and {q} ({up(float(a1.b), 10)})"
            )
    print(f"wall: {time.time() - t0:.1f} s")


# TRANSFORMS


def dker(lo, hi, u):
    n = hi - lo + 1
    z = E(u)
    with np.errstate(all="ignore"):
        v = E(lo * u) * (E(n * u) - 1) / (z - 1)
    return np.where(np.abs(z - 1) < 1e-13, n + 0j, v)


def entries(p, u):
    I = intervals(p)
    D = {k: dker(*I[k], u) for k in I}
    return (
        D[(0, 0)] + D[(0, 1)],
        D[(0, 1)],
        D[(1, 0)] + D[(1, 1)] - D[(0, 0)] - D[(0, 1)],
        D[(1, 1)] - D[(0, 1)],
        D,
    )


def gsums(p, t):
    u = (t + np.arange(p)) / p
    f, d01, n10, n11, _ = entries(p, u)
    return np.abs(f).sum(), np.abs(d01).sum(), np.abs(n10).sum(), np.abs(n11).sum()


def member(n, p):
    L = 1
    while p**L <= int(np.max(n)):
        L += 1
    return run(n, p, L, 0)[0]


def check():
    t0 = time.time()
    K = constants(verbose=False)
    ts = np.linspace(0, 0.5, 201)
    nu, H = nu_h(AF, ts)
    nu1, H1 = nu_h(A1, ts)
    worst = [0, 0, 0, 0]
    for p in [q for q in range(13, 402) if isprime(q)] + [
        1009,
        1013,
        10007,
        10009,
        100003,
        100019,
    ]:
        lam, _, eta = cert_iv(p, K)
        lam, eta, fill = float(lam.b), float(eta.b), (p + 1) / 2
        n1 = (p - 1) // 6 if p % 6 == 1 else (p + 1) // 6
        for i, t in enumerate(ts):
            g, g01, c, d = gsums(p, t)
            gb = (
                fill
                + (p / (2 * math.pi)) * (2 * nu[i] * math.log(p / 2) + H[i])
                + 1.5 * (p - 1)
                + (1 - 2 / math.pi) * p
            )
            g1b = (
                n1
                + (p / (2 * math.pi)) * (2 * nu1[i] * math.log(p / 2) + H1[i])
                + (p - 1) / 6
                + (1 - 2 / math.pi) * p / 2
            )
            sig = (n1 * t) % 1
            cc = 2 * math.cos(math.pi * (sig - 0.5) / p) / math.sin(math.pi / (2 * p))
            assert abs(c - cc) <= 1e-9 * cc and (
                abs(d - p) <= 1e-9 * p if p % 6 == 1 else abs(d - cc) <= 1e-9 * cc
            )
            assert g < gb and g01 < g1b
            assert g + eta * c < lam and g01 + eta * d < lam * eta
            worst = [
                max(worst[0], g / gb),
                max(worst[1], g01 / g1b),
                max(worst[2], (g + eta * c) / lam),
                max(worst[3], (g01 + eta * d) / (lam * eta)),
            ]
    print(
        f"bounds at every prime 13..401 and six larger on 201 shifts of [0, 1/2]: G*/bound <= {up(worst[0])}, G01/bound <= {up(worst[1])}, C and D meet their closed forms, (T psi)_0/lambda <= {up(worst[2])}, (T psi)_1/(lambda eta) <= {up(worst[3])}"
    )
    rng = np.random.default_rng(1)
    err = 0.0
    for p in (13, 31):
        k = np.arange(p**3)
        kk = k[member(k, p)]
        for u in rng.random(20):
            direct = np.exp(2j * np.pi * kk * u).sum()
            Mprod = np.eye(2, dtype=complex)
            Nprod = np.eye(2, dtype=complex)
            for j in range(3):
                f, d01, n10, n11, D = entries(p, np.array([u * p**j]))
                Mprod = Mprod @ np.array(
                    [[D[(0, 0)][0], D[(0, 1)][0]], [D[(1, 0)][0], D[(1, 1)][0]]]
                )
                Nprod = Nprod @ np.array([[f[0], d01[0]], [n10[0], n11[0]]])
            err = max(
                err,
                abs(direct - Mprod[0].sum()) / abs(direct),
                abs(direct - Nprod[0, 0]) / abs(direct),
            )
    print(
        f"transfer identity at level 3, p 13 and 31, 20 points each: relative error <= {err:.1e}"
    )
    assert err < 1e-9
    worst = 0.0
    for p in (13, 17, 19, 23, 29, 31, 37, 41, 43):
        lam = float(cert_iv(p, K)[0].b)
        i = 1
        while p**i <= 2000000:
            n = p**i
            k = np.arange(n)
            ind = member(k, p).astype(float)
            for sg in (0.0, 0.25, 0.5, 0.8):
                S = np.abs(np.fft.fft(ind * np.exp(2j * np.pi * k * sg / n))).sum()
                worst = max(worst, S / lam**i)
            i += 1
    print(
        f"grid sums: Sigma_i(s)/lambda^i <= {up(worst)} at nine primes 13..43, four shifts, every level with p^i <= 2*10^6"
    )
    assert worst < 1
    print(f"check: {time.time() - t0:.1f} s")


# RATE


def power(p, N=100, it=25):
    tg = np.linspace(0, 1, N + 1)
    r = np.arange(p)
    U = [(tg[i] + r) / p for i in range(N + 1)]
    ent = (
        [tuple(np.abs(x) for x in entries(p, u)[:4]) for u in U] if p < 20000 else None
    )
    psi0 = np.ones(N + 1)
    psi1 = np.ones(N + 1)
    for _ in range(it):
        a0 = np.zeros(N + 1)
        a1 = np.zeros(N + 1)
        for i in range(N + 1):
            f, d01, n10, n11 = (
                ent[i] if ent else tuple(np.abs(x) for x in entries(p, U[i])[:4])
            )
            q0 = np.interp(U[i], tg, psi0)
            q1 = np.interp(U[i], tg, psi1)
            a0[i] = (q0 * f + q1 * n10).sum()
            a1[i] = (q0 * d01 + q1 * n11).sum()
        hi = np.max(np.maximum(a0 / psi0, a1 / psi1))
        lo = np.min(np.minimum(a0 / psi0, a1 / psi1))
        m = max(a0.max(), a1.max())
        psi0, psi1 = a0 / m, a1 / m
    return lo, hi


def rate():
    t0 = time.time()
    K = constants(verbose=False)
    lead = 10 / (3 * math.pi)
    print(
        "power iteration of T on 100 cells, readings: lambda_T/fill, minus (10/(3 pi)) log p, alpha_T, against the certificate"
    )
    for ex in (3, 4):
        for cl in (1, 5, 7, 11):
            q = nextprime(10**ex)
            while q % 12 != cl:
                q = nextprime(q)
            lo, hi = power(q)
            fill = (q + 1) / 2
            cf = float(cert_iv(q, K)[0].b) / fill
            print(
                f"  p {q} ({q % 3},{q % 4}): [{lo / fill:.4f}, {hi / fill:.4f}], offset {hi / fill - lead * math.log(q):.4f}, alpha_T {math.log(hi / fill) / math.log(q):.4f}; certificate {cf:.4f}, alpha_1 {math.log(cf) / math.log(q):.4f}"
            )
    q = 100003
    lo, hi = power(q, N=60, it=15)
    fill = (q + 1) / 2
    off_t = hi / fill - lead * math.log(q)
    print(
        f"  p {q} ({q % 3},{q % 4}): [{lo / fill:.4f}, {hi / fill:.4f}], offset {off_t:.4f}, alpha_T {math.log(hi / fill) / math.log(q):.4f}"
    )
    print(
        "exact unshifted grid sums by FFT, readings: Sigma_L/Sigma_(L-1)/fill for K and for the row-0 design F*, minus (10/(3 pi)) log p"
    )
    for p in (101, 103, 107, 109, 1009, 1013, 1019, 1051):
        L = 1
        while p ** (L + 1) <= 3000000:
            L += 1
        fill = (p + 1) / 2
        res = []
        xk = 0.0
        for name in ("K", "F*"):
            S = []
            for i in (L - 1, L):
                k = np.arange(p**i)
                if name == "K":
                    ind = member(k, p)
                else:
                    Fa = np.zeros(p, dtype=bool)
                    I = intervals(p)
                    Fa[I[(0, 0)][0] : I[(0, 0)][1] + 1] = True
                    Fa[I[(0, 1)][0] : I[(0, 1)][1] + 1] = True
                    ind = np.ones(len(k), dtype=bool)
                    kk = k.copy()
                    for _ in range(i):
                        ind &= Fa[kk % p]
                        kk //= p
                S.append(np.abs(np.fft.fft(ind.astype(float))).sum())
            x = S[1] / S[0] / fill
            xk = x if name == "K" else xk
            res.append(f"{name} {x:.4f} ({x - lead * math.log(p):+.4f})")
        print(f"  p {p} ({p % 3},{p % 4}) level {L}: " + ", ".join(res))
        off_k = xk - lead * math.log(p) if p > 1000 else None

    def meet(off):
        lo, hi = 1e3, 1e12
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            lo, hi = (mid, hi) if mid**0.2 < lead * math.log(mid) + off else (lo, mid)
        return hi

    print(
        f"readings carried to the bar 1/5, bounding nothing: (10/(3 pi)) log p + {off_t:.4f} meets p^(1/5) near {meet(off_t):.3e}, + {off_k:.4f} near {meet(off_k):.3e}"
    )
    print(f"rate: {time.time() - t0:.1f} s")


# METER


def count_k(x, p):
    fill = (p + 1) // 2
    M0 = np.array(
        [
            [
                len(range(lo, hi + 1))
                for lo, hi in (intervals(p)[(c, 0)], intervals(p)[(c, 1)])
            ]
            for c in (0, 1)
        ],
        dtype=object,
    )
    dig = []
    y = x
    while y:
        dig.append(y % p)
        y //= p
    L = len(dig)

    def acc(c, ds):
        ok, _ = run(np.array([sum(d * p**i for i, d in enumerate(ds))]), p, len(ds), c)
        return bool(ok[0])

    vec = [np.array([1, 0], dtype=object)]
    for _ in range(L):
        vec.append(vec[-1].dot(M0))
    tot = 0
    for j in range(L):
        for f in range(dig[j]):
            for c in (0, 1):
                if vec[j][c] and acc(c, [f] + dig[j + 1 :]):
                    tot += int(vec[j][c])
    assert sum(int(v) for v in vec[L]) == fill**L
    return tot + int(acc(0, dig)) - 1


def small():
    t0 = time.time()
    X = 10**8
    odd = np.ones(X // 2 + 1, dtype=bool)
    odd[0] = False
    for q in range(3, int(X**0.5) + 1, 2):
        if odd[q // 2]:
            odd[q * q // 2 :: q] = False
    pr = np.concatenate(([2], 2 * np.nonzero(odd)[0] + 1))
    pr = pr[pr <= X]
    for p in (5, 7, 11, 13):
        A = count_k(X, p)
        L = len(np.base_repr(X, p))
        ok, _ = run(pr, p, L, 0)
        psi = float(np.log(pr[ok].astype(float)).sum())
        for q in pr[pr <= 10**4]:
            q = int(q)
            pw = q * q
            while pw <= X:
                if run(np.array([pw]), p, L, 0)[0][0]:
                    psi += math.log(q)
                pw *= q
        kap = p / (p + 1)
        small_ok = count_k(10**6, p) == int(member(np.arange(1, 10**6 + 1), p).sum())
        assert small_ok
        print(
            f"p {p}: A_K(10^8) = {A}, psi_K(10^8)/((p/(p+1)) A_K(10^8)) = {psi / (kap * A):.4f}"
        )
    print(f"small: {time.time() - t0:.1f} s")


def meter():
    t0 = time.time()
    X = 10**7
    mu = np.ones(X + 1, dtype=np.int8)
    lam = np.zeros(X + 1)
    for q in primes_upto(X):
        q = int(q)
        mu[::q] *= -1
        if q * q <= X:
            mu[:: q * q] = 0
        pw = q
        while pw <= X:
            lam[pw] = math.log(q)
            pw *= q
    mu[0] = 0
    n = np.arange(X + 1)
    for p in (101, 1009, 10007):
        inK = member(n, p)
        inK[0] = False
        A = np.cumsum(inK)
        M = np.cumsum(np.where(inK, mu, 0))
        psi = np.cumsum(np.where(inK, lam, 0.0))
        kap = p / (p + 1)
        print(
            f"p {p}: A_K(10^7) = {A[-1]}, max abs M_K/A_K^(1/2) = {np.abs(M).max() / math.sqrt(A[-1]):.3f}, sum Lambda/(kappa A_K) = {psi[-1] / (kap * A[-1]):.4f}"
        )
    print(f"meter: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "kummer"
    {
        "kummer": kummer,
        "family": family,
        "wall": wall,
        "check": check,
        "rate": rate,
        "meter": meter,
        "small": small,
    }[verb]()
