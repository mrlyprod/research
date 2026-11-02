import math
import sys
import time
from fractions import Fraction

import numpy as np
from mpmath import iv, mp

GAMMA = 0.5772156649015329
ZETA3 = 1.2020569031595942
C_P = math.sqrt(2) - 4 / math.pi
K_2 = 2 * (4 / 3 + (36 / 35) * (7 * ZETA3 / 8 - 1))
KAP = 2 * (1 - 2 / math.pi)
TAIL = 3.4
CS = [
    Fraction(9, 10),
    Fraction(3, 4),
    Fraction(1, 2),
    Fraction(1, 3),
    Fraction(1, 4),
    Fraction(1, 5),
    Fraction(1, 10),
]

# CLOSED FORMS


def x0(N):
    return (N / math.pi) * (
        math.log(N + 3)
        + GAMMA
        + math.log(math.tan(3 * math.pi / 8 + math.pi / (4 * N)))
    ) + C_P * (N + 1) ** 2 / (8 * N)


def lam_odd(N):
    return (N - 1) / 2 + x0(N) + 0.5 / math.sin(math.pi / (2 * N))


def y_even(N):
    return (N / math.pi) * (math.log(N + 4) + GAMMA) + KAP * (N + 2) ** 2 / (8 * N)


def lam_even(N):
    return 3 * N / 4 + y_even(N)


def harm(n):
    return float(np.sum(1.0 / np.arange(1, n + 1))) if n else 0.0


def big_lam(L, N, h=None):
    h = harm(N // 2) if h is None else h
    return 2 * L + (2 * N / math.pi) * h + (1 - 2 / math.pi) * (N * N + 1) / (2 * N)


def iv_lamF(N):
    pi = iv.pi
    cp = iv.sqrt(2) - 4 / pi
    x = (N / pi) * (
        iv.log(N + 3) + iv.euler + iv.log(iv.tan(3 * pi / 8 + pi / (4 * N)))
    ) + cp * (N + 1) ** 2 / (8 * N)
    return (N + 1) / 2 + x + 1 / (2 * iv.sin(pi / (2 * N)))


def rate_odd(N):
    iv.prec = 120
    N = iv.mpf(N)
    return (iv_lamF(N) - 1) / ((N - 1) / 2)


def rate_even(N):
    iv.prec = 120
    N = iv.mpf(N)
    kap = 2 * (1 - 2 / iv.pi)
    return (
        3 * N / 4
        + (N / iv.pi) * (iv.log(N + 4) + iv.euler)
        + kap * (N + 2) ** 2 / (8 * N)
    ) / (N / 2)


def rate_unif(N, c):
    iv.prec = 120
    N = iv.mpf(N)
    k = (2 / iv.pi) * (iv.log(N / 2) + iv.euler + 1 / N) + (1 - 2 / iv.pi) * (
        N * N + 1
    ) / (2 * N * N)
    return 2 + k * iv.mpf(c.denominator) / c.numerator


def c_inf():
    iv.prec = 120
    pi = iv.pi
    return (
        1
        + 2 / pi
        + (2 / pi) * (iv.euler + iv.log(1 + iv.sqrt(2)))
        + (iv.sqrt(2) - 4 / pi) / 4
    )


def gap(rate, N, d):
    iv.prec = 120
    return iv.mpf(N) ** (iv.mpf(1) / d) - rate(N)


def down(x, d=5):
    v = mp.mpf(x.a.a)
    k = d - 1 - int(mp.floor(mp.log10(abs(v))))
    return mp.nstr(mp.floor(v * 10**k) / mp.mpf(10) ** k, d)


def up(x, d=5):
    v = mp.mpf(x.b.b)
    k = d - 1 - int(mp.floor(mp.log10(abs(v))))
    return mp.nstr(mp.ceil(v * 10**k) / mp.mpf(10) ** k, d)


def fup(x, d=6):
    return f"{math.ceil(x * 10**d) / 10**d:.{d}f}"


def fdown(x, d=6):
    return f"{math.floor(x * 10**d) / 10**d:.{d}f}"


def gdown(x, d=4):
    k = d - 1 - math.floor(math.log10(x))
    return f"{math.floor(x * 10**k) / 10**k:.{d - 1}e}"


def bisect(f, lo, hi, step):
    while hi - lo > step:
        mid = lo + ((hi - lo) // (2 * step)) * step
        if f(mid):
            hi = mid
        else:
            lo = mid
    return hi


# WALL


def mono_from(d, c=1.0, odd=False):
    N = 3
    f = (
        (lambda N: N ** (1 / d) * (N - 1) > (2 * d / math.pi) * (N + 1))
        if odd
        else (lambda N: N ** (1 / d) > 2 * d / (math.pi * c))
    )
    N = bisect(f, 2, 10**15, 1)
    iv.prec = 120
    test = (
        (lambda M: iv.mpf(M) ** (iv.mpf(1) / d) * (M - 1) - (2 * d / iv.pi) * (M + 1))
        if odd
        else (lambda M: iv.mpf(M) ** (iv.mpf(1) / d) - 2 * d / (iv.pi * c))
    )
    while test(N).a <= 0:
        N += 1
    return N


def wall_even(d):
    m = mono_from(d)
    f = lambda N: N ** (1 / d) > float(rate_even(N).b)
    W = bisect(f, m + m % 2, 10**9, 2)
    assert W >= m and gap(rate_even, W, d).a > 0 and gap(rate_even, W - 2, d).b < 0
    a1 = iv.mpf(1) / d - iv.log(rate_even(W)) / iv.log(W)
    print(
        f"  half interval at even base, bar 1/{d}: the closed form clears at every even base >= {W}, gap >= {down(gap(rate_even, W, d))} there, fails at {W - 2} with gap <= {up(gap(rate_even, W - 2, d))}; the gap increases from {m}; bar - alpha_1 >= {down(a1)} at the wall"
    )
    return W


def wall_odd(d):
    c = c_inf()
    tail = lambda N: (2 / iv.pi) * iv.log(N) + c + TAIL / iv.mpf(N)
    tgap = lambda N: iv.mpf(N) ** (iv.mpf(1) / d) * (N - 1) / (N + 1) - tail(N)
    m = max(mono_from(d, odd=True), 101)
    na = bisect(lambda N: tgap(N).a > 0, m, 10**9, 1)
    assert tgap(na).a > 0 and na > m
    top = na if na % 2 else na + 1
    N = top
    while gap(rate_odd, N - 2, d).a > 0:
        N -= 2
    for M in range(N, top + 1, 2):
        assert gap(rate_odd, M, d).a > 0
    a1 = iv.mpf(1) / d - iv.log(rate_odd(N)) / iv.log(N)
    print(
        f"  odd digits, bar 1/{d}: the tail (lambda/fill bound)(base+1)/(base-1) clears from {na}, increasing from {m}; the closed form (lambda - 1)/(fill - 1) certified at every odd base {N}..{top}, fails at {N - 2} with gap <= {up(gap(rate_odd, N - 2, d))}; bar - alpha_1 >= {down(a1)} at the wall {N}"
    )
    return N


def wall_unif(c):
    m = mono_from(5, float(c))
    assert gap(lambda N: rate_unif(N, c), m, 5).b < 0
    W = bisect(lambda N: N**0.2 > float(rate_unif(N, c).b), m, 10**15, 1)
    g = lambda N: gap(lambda M: rate_unif(M, c), N, 5)
    assert g(W).a > 0 and g(W - 1).b < 0
    return W, m, g(W), g(W - 1)


def wall():
    t0 = time.time()
    print(
        f"c_inf <= {up(c_inf(), 9)} (the half interval's tail constant, lambda/fill <= (2/pi) log base + c_inf + {TAIL}/base from base 101)"
    )
    wo5 = wall_odd(5)
    wall_odd(4)
    we5 = wall_even(5)
    wall_even(4)
    iv.prec = 120
    k = iv.mpf(3) / 2 + (2 / iv.pi) * iv.euler + (1 - 2 / iv.pi) / 2
    print(
        f"  half interval at even base: lambda_H/L - (2/pi) log base tends to 3/2 + (2/pi) gamma + (1 - 2/pi)/2 <= {up(k, 6)}"
    )
    for N in (100003, 10**6 + 3, 10**9 + 7):
        a = iv.log(rate_odd(N)) / iv.log(N)
        b = iv.log(rate_even(N + 1)) / iv.log(N + 1)
        print(
            f"  alpha_1 <= {up(a, 7)} for the odd digits at {N}, <= {up(b, 7)} for the half interval at {N + 1}"
        )
    print(
        "uniform certificate Lambda(L, base) = 2L + (2 base/pi) H_floor(base/2) + (1 - 2/pi)(base^2 + 1)/(2 base); B(c) = least base from which 2 + K(base)/c < base^(1/5)"
    )
    for c in CS:
        W, m, gw, gb = wall_unif(c)
        print(
            f"  c = {c}: B(c) = {W}, gap >= {down(gw)} there, fails at {W - 1} with gap <= {up(gb)}, increasing from {m}"
        )
    print(f"wall in {time.time() - t0:.1f}s")
    return wo5, we5


# CENSUS


def gsum(L, N, ts, chunk=3_000_000):
    r = np.arange(N, dtype=np.float64)
    out = np.empty(len(ts))
    step = max(1, chunk // N)
    for i in range(0, len(ts), step):
        u = (ts[i : i + step, None] + r[None, :]) / N
        den = np.sin(np.pi * u)
        num = np.abs(np.sin(np.pi * L * u))
        with np.errstate(divide="ignore", invalid="ignore"):
            out[i : i + step] = np.where(
                den > 1e-300, num / np.where(den > 1e-300, den, 1), L
            ).sum(1)
    return out


def readings(N, T):
    ts = np.linspace(0.0, 0.5, T)
    h = harm(N // 2)
    worst = (0.0, None)
    for L in range(2, N):
        g = gsum(L, N, ts)
        q = g.max() / big_lam(L, N, h)
        assert q < 1
        if q > worst[0]:
            worst = (q, (N, L, float(ts[g.argmax()])))
    return worst


def frac_grid(N, L, r, i, s, u=0):
    y = N**i
    a = np.arange(y, dtype=np.int64)
    out = np.ones(y)
    for j in range(i):
        th = ((pow(N, j, y) * a) % y) / y + float((Fraction(N) ** j * s) % 1)
        z = np.zeros(y, dtype=complex)
        for v in range(L):
            z += np.exp(2j * np.pi * (u + r * v) * th)
        out *= np.abs(z)
    return out


def interval_grid(N, L, i, s, sub=1):
    y = N**i
    a = np.arange(0, y, sub, dtype=np.int64)
    out = np.ones(len(a))
    for j in range(i):
        th = (((pow(N, j, y) * a) % y) / y + float((Fraction(N) ** j * s) % 1)) % 1.0
        den = np.sin(np.pi * th)
        num = np.abs(np.sin(np.pi * L * th))
        with np.errstate(divide="ignore", invalid="ignore"):
            out *= np.where(
                np.abs(den) > 1e-12,
                num / np.where(np.abs(den) > 1e-12, np.abs(den), 1),
                L,
            )
    return float(out.sum())


def smooth(r, N):
    g = 1
    while math.gcd(r // g, N) > 1:
        g *= math.gcd(r // g, N)
    return g


def census():
    t0 = time.time()
    print(
        "THE ONE-STEP BOUND: max over 201 shifts in [0, 1/2] of G_L(t) = sum_r |D_L((t+r)/base)| against Lambda(L, base), every length L = 2..base-1"
    )
    worst = (0.0, None)
    for N in range(3, 81):
        worst = max(worst, readings(N, 201))
    print(
        f"  every base 3..80: largest reading/Lambda <= {fup(worst[0])} at base {worst[1][0]}, L = {worst[1][1]}, t = {worst[1][2]:.3f}"
    )
    for N in (99, 100, 101, 102, 103, 104, 105, 106, 1000, 1001, 1002, 1003):
        w = readings(N, 101 if N < 1000 else 41)
        worst = max(worst, w)
        print(
            f"  base {N} ({N % 2} mod 2, {N % 3} mod 3, {N % 4} mod 4): largest reading/Lambda <= {fup(w[0])} at L = {w[1][1]}"
        )
    print(
        f"  overall worst cell: base {worst[1][0]}, L = {worst[1][1]}, reading/Lambda <= {fup(worst[0])}"
    )
    ts = np.linspace(0.0, 0.5, 201)
    print(
        "  slack: base, L, max G_L/L reading, Lambda/L, the sharp lambda/L where proved"
    )
    for N, L, sharp in (
        (10007, 5004, lam_odd(10007) + 1),
        (10007, 5003, lam_odd(10007)),
        (10006, 5003, lam_even(10006)),
        (10007, 3336, None),
        (10007, 2502, None),
        (10007, 10006, None),
    ):
        g = gsum(L, N, ts).max()
        sh = f"{sharp / L:.4f}" if sharp else "none"
        print(f"    {N} {L} {g / L:.4f} {big_lam(L, N) / L:.4f} {sh}")
    print(
        "THE CUT: alpha_1(L, base) = log_base(Lambda(L, base)/L), decreasing in L; the band within 0.02 of 1/5"
    )
    lo = min(
        math.log(big_lam(L, N) / L) / math.log(N)
        for N in list(range(3, 81)) + [99, 106, 1000, 1003]
        for L in (N - 1,)
    )
    print(
        f"  bases 3..80, 99..106, 1000..1003: smallest alpha_1 over every cell >= {fdown(lo)}, no cell within 0.02 of the cut"
    )
    for c in CS:
        W, *_ = wall_unif(c)
        h = math.log(W // 2) + GAMMA + 1 / (2 * (W // 2))
        q = (2 * W / math.pi) * h + (1 - 2 / math.pi) * (W * W + 1) / (2 * W)
        lc = math.ceil(c * W)
        a1 = math.log(2 + q / lc) / math.log(W)
        l_hi = min(W - 1, math.floor(q / (W**0.18 - 2)))
        l_lo = max(lc, math.ceil(q / (W**0.22 - 2)))
        print(
            f"  c = {c}, base B(c) = {W}: worst cell L = {lc}, alpha_1 = {a1:.7f}; within 0.02 of the cut: L = {l_lo}..{l_hi}, {max(0, l_hi - l_lo + 1)} of the {W - lc} covered lengths"
        )
    print(
        "PROGRESSIONS FROM THEIR OWN DIGITS: Sigma^P_i(s) against r_1 Sigma^I_i(r s) and the identity through gcd(r, base^i), every (r, L) with u the largest start"
    )
    err, rat, cnt = 0.0, 0.0, 0
    for N in range(3, 14):
        h = harm(N // 2)
        for r in range(1, N):
            for L in range(2, (N - 1) // r + 2):
                u = N - 1 - (L - 1) * r
                for i in range(1, 6):
                    if N**i * L > 400_000:
                        break
                    for s in (
                        Fraction(0),
                        Fraction(1, 2),
                        Fraction(1, 3),
                        Fraction(3819660113, 10**10),
                    ):
                        p = frac_grid(N, L, r, i, s, u).sum()
                        g = math.gcd(r, N**i)
                        rs = (r * s) % 1
                        ident = g * interval_grid(N, L, i, rs, g)
                        full = interval_grid(N, L, i, rs)
                        err = max(err, abs(p - ident) / p)
                        r1 = smooth(r, N)
                        assert p <= r1 * full * (1 + 1e-9) and full <= big_lam(
                            L, N, h
                        ) ** i * (1 + 1e-9)
                        rat = max(rat, p / (r1 * big_lam(L, N, h) ** i))
                        cnt += 1
    assert err < 1e-9
    print(
        f"  {cnt} cells at bases 3..13, levels with base^i L <= 4 * 10^5, 4 shifts: identity to {gdown(err) if err else 0}, largest Sigma^P/(r_1 Lambda^i) <= {fup(rat)}"
    )
    print(f"census in {time.time() - t0:.1f}s")


# CHECK


def odd_pieces(N, ts):
    n = (N + 1) // 2
    K = (n - 1) // 2
    J = (n - 2) // 2
    worst = [-1e9, -1e9]
    for t in ts:
        k = np.arange(1, K + 1, dtype=np.float64)
        le = np.sum(1 / (2 * np.sin(np.pi * (t + 2 * k) / (2 * N))))
        ho = np.sum(1 / (2 * np.sin(np.pi * (2 * k - t) / (2 * N))))
        bound = (N / math.pi) * (math.log(N + 4) - 0.0757) + C_P * N / 4
        worst[0] = max(worst[0], (le + ho) / bound)
        j = np.arange(J + 1, dtype=np.float64)
        q = lambda d: 1 / (2 * np.cos(np.pi * d / 2))
        res = np.sum(q((2 * j + 1 + t) / N) - q((2 * j + 1 - t) / N))
        worst[1] = max(worst[1], res)
    return worst


def check_odd(N, T):
    L = (N - 1) // 2
    ts = np.linspace(0.0, 0.5, T)
    G = gsum(L, N, ts)
    r = np.arange(N)
    S = np.array([np.abs(np.sin(np.pi * L * (t + r) / N)).sum() for t in ts])
    sig = np.mod(L * ts, 1.0)
    Sc = np.cos(np.pi * (sig - 0.5) / N) / np.sin(np.pi / (2 * N))
    assert np.max(np.abs(S - Sc) / Sc) < 1e-9
    a, b = np.sin(np.pi * ts / 2), np.cos(np.pi * ts / 2)
    X = x0(N)
    gb = L + (a + b) * X + (N / math.pi) * K_2 * ts**2
    lam = lam_odd(N)
    sn = np.sin(np.pi * ts)
    r1, r2 = np.max(G / gb), np.max((G + S / 2) / (lam * (1 + sn / 2)))
    assert r1 < 1 and r2 < 1 and K_2 * N / (2 * math.pi) <= L + np.min(S) / 2
    return r1, r2


def check_even(N, T):
    L = N // 2
    ts = np.linspace(0.0, 0.5, T)
    G = gsum(L, N, ts)
    r = np.arange(N)
    S = np.array([np.abs(np.sin(np.pi * L * (t + r) / N)).sum() for t in ts])
    a, b = np.sin(np.pi * ts / 2), np.cos(np.pi * ts / 2)
    assert np.max(np.abs(S - (N / 2) * (a + b))) < 1e-7 * N
    gb = L + (a + b) * y_even(N) + (N / math.pi) * K_2 * ts**2
    lam = lam_even(N)
    sn = np.sin(np.pi * ts)
    r1, r2 = np.max(G / gb), np.max((G + S / 2) / (lam * (1 + sn / 2)))
    assert r1 < 1 and r2 < 1
    return r1, r2


def check():
    t0 = time.time()
    w = [-1e9, -1e9]
    for N in list(range(3, 402, 2)) + [1001, 10001, 100001]:
        p = odd_pieces(N, np.linspace(0.0, 0.5, 101))
        w = [max(x, y) for x, y in zip(w, p)]
    assert w[0] < 1 and w[1] < math.sqrt(0.5) - 0.5
    print(
        f"odd digits, the two new pieces at every odd base 3..401 and 1001, 10001, 100001 on 101 shifts: unsigned aP sum over its bound <= {fup(w[0])}, aQ residue <= {fup(w[1])} against 1/sqrt 2 - 1/2"
    )
    bases_o = list(range(3, 402, 2)) + [1001, 3799, 10001, 30001, 94953, 100001]
    bases_e = list(range(4, 401, 2)) + [1000, 2414, 10000, 30000, 61270, 100000]
    for name, fn, bases in (
        ("odd digits (L = (base-1)/2)", check_odd, bases_o),
        ("half interval at even base (L = base/2)", check_even, bases_e),
    ):
        m1, m2 = 0.0, 0.0
        for N in bases:
            r1, r2 = fn(N, 801 if N < 1000 else 201)
            m1, m2 = max(m1, r1), max(m2, r2)
        print(
            f"{name}: {len(bases)} bases {bases[0]}..{bases[-1]}, max G/bound <= {fup(m1)}, max T phi/(lambda phi) <= {fup(m2)}"
        )
    print(
        "grid sums from the digits against (3/2) lambda^i, every level with base^i <= 2 * 10^5, shifts 0, 1/2, 1/3, 0.381966"
    )
    for name, lamf, bases, digits in (
        ("odd digits", lam_odd, range(3, 22, 2), lambda N: (1, 2, (N - 1) // 2)),
        (
            "half interval, even base",
            lam_even,
            range(4, 22, 2),
            lambda N: (0, 1, N // 2),
        ),
        ("odd digits, even base", lam_even, range(4, 22, 2), lambda N: (1, 2, N // 2)),
        ("even digits, even base", lam_even, range(4, 22, 2), lambda N: (0, 2, N // 2)),
    ):
        worst = 0.0
        for N in bases:
            u, r, L = digits(N)
            C = 1.5 * smooth(r, N)
            for i in range(1, 9):
                if N**i > 200_000:
                    break
                for s in (
                    Fraction(0),
                    Fraction(1, 2),
                    Fraction(1, 3),
                    Fraction(3819660113, 10**10),
                ):
                    v = frac_grid(N, L, r, i, s, u).sum()
                    worst = max(worst, v / (C * lamf(N) ** i))
        assert worst < 1
        print(f"  {name}: largest Sigma/(C lambda^i) <= {fup(worst)}, C = (3/2) r_1")
    print(f"check in {time.time() - t0:.1f}s")


# LOST


def lost():
    t0 = time.time()
    print(
        "THE DENOMINATORS A STEP LOSES: m_d = max over l prime to d of prod_(j<200) |D_L(r base^j l/d)|/L, every d = 2..60; lost class: the part of d prime to base divides r_0"
    )
    for N, r in (
        (7, 2),
        (9, 2),
        (15, 2),
        (11, 3),
        (10, 3),
        (8, 3),
        (13, 4),
        (12, 5),
        (25, 6),
        (7, 6),
    ):
        L = (N - 1) // r + 1
        r0 = r // smooth(r, N)
        pure, mixed, kept = [], [], []
        for d in range(2, 61):
            d0 = d // smooth(d, N)
            m = 0.0
            for l in range(1, d):
                if math.gcd(l, d) != 1:
                    continue
                p = 1.0
                for j in range(200):
                    x = (r * pow(N, j, d) * l % d) / d
                    p *= (
                        1.0
                        if x == 0
                        else abs(math.sin(math.pi * L * x) / math.sin(math.pi * x)) / L
                    )
                m = max(m, p)
            if r0 % d0:
                kept.append(m)
            elif d0 > 1 and d == d0:
                pure.append(m)
            elif d0 > 1:
                mixed.append(m)
        assert max(kept) < 1e-3 and all(v == 1.0 for v in pure)
        mx = (
            f", mixed g e: {len(mixed)} with m_d in [{min(mixed):.4f}, {max(mixed):.4f}]"
            if mixed
            else ""
        )
        print(
            f"  base {N}, step {r}, L = {L}, r_0 = {r0}: pure g | r_0, g > 1: {len(pure)} with m_d = 1 by construction{mx}; off the class: {len(kept)} with m_d <= {gdown(max(kept))}"
        )
    print(f"lost in {time.time() - t0:.1f}s")


# METER


def sieve(X):
    mu = np.ones(X + 1, dtype=np.int8)
    mu[0] = 0
    isp = np.ones(X + 1, dtype=bool)
    isp[:2] = False
    for p in range(2, int(X**0.5) + 1):
        if isp[p]:
            isp[p * p :: p] = False
    for p in np.nonzero(isp)[0]:
        mu[p::p] *= -1
        if p * p <= X:
            mu[p * p :: p * p] = 0
    return mu, isp


def member(X, N, digits):
    ok = np.ones(X + 1, dtype=bool)
    ln = np.zeros(X + 1, dtype=np.int8)
    y = np.arange(X + 1, dtype=np.int64)
    allowed = np.zeros(N, dtype=bool)
    allowed[list(digits)] = True
    while y.any():
        live = y > 0
        ok &= ~live | allowed[y % N]
        ln += live
        y //= N
    ok[0] = False
    return ok, ln


def phi(N):
    return sum(1 for f in range(1, N + 1) if math.gcd(f, N) == 1)


def kappa(N, digits):
    return N * sum(1 for f in digits if math.gcd(f, N) == 1) / (phi(N) * len(digits))


def meter():
    t0 = time.time()
    X = 10**7
    mu, isp = sieve(X)
    logk = np.log(np.maximum(np.arange(X + 1), 1))
    print(
        "sanity only, below 10^7: set, A(x), M(x), max |M|/A, max |M|/sqrt(A) over A >= 100, sum log p over the set against the predicted main term"
    )
    sets = (
        (
            "odd digits {1, 3, .., 99} at base 101",
            101,
            range(1, 100, 2),
            lambda ln, k: np.where(ln % 2 == 1, 2 * k, 0.0),
        ),
        (
            "half interval {0..49} at base 100",
            100,
            range(50),
            lambda ln, k: np.full(ln.shape, k),
        ),
        (
            "odd digits {1, 3, .., 99} at base 100",
            100,
            range(1, 100, 2),
            lambda ln, k: np.full(ln.shape, k),
        ),
        (
            "step 3 {1, 4, .., 97} at base 100",
            100,
            range(1, 100, 3),
            lambda ln, k: np.where(ln % 3 != 0, 1.5 * k, 0.0),
        ),
    )
    for name, N, digits, kap in sets:
        kp = kappa(N, list(digits))
        ok, ln = member(X, N, digits)
        A = np.cumsum(ok)
        M = np.cumsum(np.where(ok, mu, 0).astype(np.int64))
        sel = A >= 100
        th = np.sum(logk[ok & isp])
        main = np.sum(kap(ln[ok], kp))
        top = int(ln[ok & isp].max())
        n_top = int(np.sum(ok & isp & (ln == top)))
        print(
            f"  {name}, kappa_P = {kp:.6f}: {A[-1]} {M[-1]} {np.max(np.abs(M[sel]) / A[sel]):.4f} {np.max(np.abs(M[sel]) / np.sqrt(A[sel])):.3f}, primes {th:.1f} against {main:.1f}, ratio {th / main:.4f}; {n_top} primes of length {top}, abs(1 - ratio) sqrt of that {abs(1 - th / main) * math.sqrt(n_top):.2f}"
        )
    print(f"meter in {time.time() - t0:.1f}s")


# HEADROOM

RHOS = (Fraction(1, 2), Fraction(3, 4), Fraction(19, 20))
PRIMES = (10007, 100003, 1000003)
N0 = 1000


def a_env(N):
    return 2 / math.pi + 0.4 * math.pi / (N * N)


def fourier_cert(L, N, K0=20000):
    om = big_lam(L, N) - 2 * L
    k = np.arange(1, K0 + 1, dtype=np.int64)
    m = (k * L) % N
    with np.errstate(divide="ignore"):
        osc = np.where(
            m == 0, np.inf, 2 / (math.sin(math.pi / N) * np.sin(np.pi * m / N))
        )
    S = float(np.sum(np.minimum(om, osc) / (4.0 * k * k - 1))) + om / (2 * (2 * K0 + 1))
    return 2 * L + (2 / math.pi) * om + (4 / math.pi) * S


def pieces(c, rho1, n):
    return [(c + (rho1 - c) * j / n, c + (rho1 - c) * (j + 1) / n) for j in range(n)]


def envelope(a, b, K0=60):
    iv.prec = 120
    coef_k = 2 / iv.pi + (4 / iv.pi) / (2 * (2 * K0 + 1))
    coef_a = iv.mpf(0)
    for k in range(1, K0 + 1):
        w = 1 / iv.mpf(4 * k * k - 1)
        if math.ceil(k * a) <= math.floor(k * b):
            coef_k += (4 / iv.pi) * w
        else:
            sa = abs(iv.sin(iv.pi * iv.mpf(k * a.numerator) / a.denominator))
            sb = abs(iv.sin(iv.pi * iv.mpf(k * b.numerator) / b.denominator))
            coef_a += (4 / iv.pi) * w / iv.mpf(min(sa.a, sb.a))
    return coef_k, coef_a


def k_env(N):
    return (2 / iv.pi) * (iv.log(N / 2) + iv.euler + 1 / N) + (1 - 2 / iv.pi) * (
        N * N + 1
    ) / (2 * N * N)


def wall_fourier(c, rho1, n):
    iv.prec = 120
    envs = [(a, envelope(a, b)) for a, b in pieces(c, rho1, n)]
    coefs = [float(ck.b) / float(a) for a, (ck, ca) in envs] + [1 / float(rho1)]
    m = max(N0, mono_from(5, 1 / max(coefs)))
    a_hi = a_env(N0)

    def gaps(N):
        M = iv.mpf(N)
        kk = k_env(M)
        g = [
            M ** (iv.mpf(1) / 5)
            - 2
            - (ck * kk + ca * a_hi) / iv.mpf(a.numerator) * a.denominator
            for a, (ck, ca) in envs
        ]
        g.append(
            M ** (iv.mpf(1) / 5) - 2 - kk / iv.mpf(rho1.numerator) * rho1.denominator
        )
        return g

    ok = lambda N: all(g.a > 0 for g in gaps(N))
    W = bisect(ok, m, 10**12, 1)
    assert ok(W) and min(g.b for g in gaps(W - 1)) < 0
    return W, m


def abel_piece(N, L, ts, K0):
    r = np.arange(1, N - 1, dtype=np.float64)
    worst = 0.0
    for t in ts:
        th = (t + r) / N
        w = 1 / np.sin(np.pi * th)
        for k in range(1, K0 + 1):
            m = (k * L) % N
            if m == 0:
                continue
            C = abs(np.sum(np.cos(2 * np.pi * k * L * th) * w))
            worst = max(
                worst, C * math.sin(math.pi * m / N) * math.sin(math.pi / N) / 2
            )
    return worst


def read_wall(r, c):
    f = lambda N: r * (2 + float(k_env(iv.mpf(N)).b) / c) < N**0.2
    return bisect(f, 10**3, 10**9, 1)


def transfer_matrix(L, N, M):
    ts = np.linspace(0.0, 0.5, M)
    h = 0.5 / (M - 1)
    A = np.zeros((M, M))
    r = np.arange(N, dtype=np.float64)
    for j, t in enumerate(ts):
        u = (t + r) / N
        den = np.sin(np.pi * u)
        num = np.abs(np.sin(np.pi * L * u))
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(den > 1e-300, num / np.where(den > 1e-300, den, 1), L)
        x = np.minimum(u, 1 - u) / h
        i0 = np.minimum(np.floor(x).astype(np.int64), M - 2)
        f = x - i0
        A[j] = np.bincount(i0, weights=w * (1 - f), minlength=M) + np.bincount(
            i0 + 1, weights=w * f, minlength=M
        )
    return A


def perron(A):
    vals, vecs = np.linalg.eig(A)
    i = np.argmax(vals.real)
    v = np.abs(vecs[:, i].real)
    return float(vals[i].real), v / v[0]


def spectral():
    print(
        "THE SPECTRAL RADIUS: T_L on a 401-point grid of t in [0, 1/2], linear interpolation, (T phi)(t) = sum_r |D_L((t+r)/N)| phi((t+r)/N); rho(T_L) is the floor of every certificate T phi <= lambda phi"
    )
    print(
        "  c, N, L, rho/L, Lambda/L, rho/Lambda, max G/Lambda, alpha_1 = log_N(rho/L), the eigenvector's phi(1/2)/phi(0) - 1, the implied wall from rho/Lambda, and rho/L on 201 points where run"
    )
    for c in (Fraction(1, 3), Fraction(1, 2), Fraction(3, 4), Fraction(9, 10)):
        for N in PRIMES:
            L = round(c * N)
            lam = big_lam(L, N)
            A = transfer_matrix(L, N, 401)
            rho, v = perron(A)
            G = A.sum(1).max()
            coarse = (
                f"{perron(transfer_matrix(L, N, 201))[0] / L:.4f}" if N < 10**6 else "-"
            )
            assert rho <= G * (1 + 1e-9)
            print(
                f"    {c} {N} {L} {rho / L:.4f} {lam / L:.4f} {rho / lam:.6f} {G / lam:.6f} {math.log(rho / L) / math.log(N):.6f} {v[-1] - 1:.4f} {read_wall(rho / lam, float(c))} {coarse}"
            )
        if c == Fraction(1, 2):
            print(
                f"    the half interval's sharp certificate reads lambda/L = {(lam_odd(10007) + 1) / 5004:.4f} at base 10007, L = 5004, against rho/L above"
            )


def headroom():
    t0 = time.time()
    spectral()
    print(
        "THE PROOF'S ONE PIECE: |sum_(r=1)^(N-2) cos(2 pi k L (t+r)/N) csc(pi (t+r)/N)| against 2 csc(pi/N)/|sin(pi k L/N)|, k = 1..300, 41 shifts"
    )
    worst = 0.0
    for N, L in (
        (101, 51),
        (101, 76),
        (101, 100),
        (1003, 502),
        (1003, 753),
        (1003, 953),
    ):
        worst = max(worst, abel_piece(N, L, np.linspace(0.0, 0.5, 41), 300))
    assert worst < 1
    print(f"  six cells, largest sum over its bound <= {fup(worst)}")
    print(
        "THE FOURIER CERTIFICATE AGAINST THE ONE-STEP READING: max over 201 shifts of G_L(t) against lambda_F(L, base), every length 2..base-1 of every base 3..80"
    )
    ts = np.linspace(0.0, 0.5, 201)
    worst, below = 0.0, 0
    for N in range(3, 81):
        for L in range(2, N):
            lf = fourier_cert(L, N, 2000)
            q = gsum(L, N, ts).max() / lf
            assert q < 1
            worst = max(worst, q)
            below += lf < big_lam(L, N)
    print(
        f"  largest reading/lambda_F <= {fup(worst)}; lambda_F below Lambda in {below} of the {sum(N - 2 for N in range(3, 81))} cells"
    )
    print(
        "THE HEADROOM: max over 401 shifts t in [0, 1/2] of G_L(t) = sum_r |D_L((t + r)/N)|, L = rho N rounded, N prime; against Lambda(L, N) and against the Fourier certificate lambda_F(L, N)"
    )
    print(
        "  rho, N, L, t at the max, G/L, Lambda/L, lambda_F/L, G/Lambda, G/lambda_F, lambda_F/Lambda"
    )
    ts = np.linspace(0.0, 0.5, 401)
    for rho in RHOS:
        q = rho.denominator
        print(
            f"  rho = {rho}: the leading constant of the resonance at q = {q} is 1/(q sin(pi/(2q))) = {1 / (q * math.sin(math.pi / (2 * q))):.6f}, against 2/pi = {2 / math.pi:.6f}"
        )
        for N in PRIMES:
            L = round(rho * N)
            g = gsum(L, N, ts)
            G, t = g.max(), ts[g.argmax()]
            lam, lf = big_lam(L, N), fourier_cert(L, N)
            assert G < lf < 2 * lam
            print(
                f"    {rho} {N} {L} {t:.4f} {G / L:.4f} {lam / L:.4f} {lf / L:.4f} {fup(G / lam)} {fup(G / lf)} {fup(lf / lam)}"
            )
    print(
        "THE WALLS OF THE FOURIER CERTIFICATE: B_F(c), the least base from which every L >= c base clears 1/5 through min(Lambda, lambda_F), the envelope on [c, rho_1] in pieces and Lambda alone above rho_1"
    )
    for c in (Fraction(1, 2), Fraction(3, 4), Fraction(9, 10)):
        best = None
        for j in range(1, int(100 * (1 - c))):
            rho1 = c + Fraction(j, 100)
            W, m = wall_fourier(c, rho1, 8)
            if best is None or W < best[0]:
                best = (W, m, rho1)
        W, m, rho1 = best
        print(
            f"  c = {c}: B_F(c) = {W}, rho_1 = {rho1}, 8 pieces, certified at {W}, failing at {W - 1}, gaps increasing from {m}"
        )
    W1, *_ = wall_unif(Fraction(1, 1))
    g = gsum(W1 - 1, W1, np.array([0.5]))[0]
    print(
        f"  the floor of every envelope whose top piece is Lambda alone: 2 + K(base) < base^(1/5) first at {W1}; there the reading at L = base - 1, t = 1/2 is {g / big_lam(W1 - 1, W1):.4f} of Lambda"
    )
    for N in (11003, 12007):
        g = gsum(N - 1, N, np.linspace(0.0, 0.5, 401)).max()
        print(
            f"  the one-step reading itself at L = base - 1: base {N}, max G/L = {g / (N - 1):.4f} against base^(1/5) = {N**0.2:.4f}"
        )
    for r, name in (
        (0.688172, "the census ratio at rho = 1/2, base 1000003"),
        (1 / math.sqrt(2), "the limit constant 1/sqrt 2"),
    ):
        print(
            f"  read, not proved: with G <= {r:.6f} Lambda at every L >= base/2, {name}, the wall would be {read_wall(r, 0.5)}"
        )
    print(f"headroom in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "wall"
    {
        "wall": wall,
        "census": census,
        "check": check,
        "lost": lost,
        "meter": meter,
        "headroom": headroom,
    }[verb]()
