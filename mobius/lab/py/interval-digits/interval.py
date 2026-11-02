import math
import sys
import time

import numpy as np
from mpmath import iv, mp

GAMMA = 0.5772156649015329
ZETA3 = 1.2020569031595942
C_P = math.sqrt(2) - 4 / math.pi
K_2 = 2 * (4 / 3 + (36 / 35) * (7 * ZETA3 / 8 - 1))
TAIL = 3.4

# CLOSED FORMS


def x0_float(N):
    return (N / math.pi) * (
        math.log(N + 3)
        + GAMMA
        + math.log(math.tan(3 * math.pi / 8 + math.pi / (4 * N)))
    ) + C_P * (N + 1) ** 2 / (8 * N)


def lam_float(N):
    return (N + 1) / 2 + x0_float(N) + 0.5 / math.sin(math.pi / (2 * N))


def rate_iv(N):
    iv.prec = 120
    N = iv.mpf(N)
    pi = iv.pi
    cp = iv.sqrt(2) - 4 / pi
    x0 = (N / pi) * (
        iv.log(N + 3) + iv.euler + iv.log(iv.tan(3 * pi / 8 + pi / (4 * N)))
    ) + cp * (N + 1) ** 2 / (8 * N)
    n = (N + 1) / 2
    return (n + x0 + 1 / (2 * iv.sin(pi / (2 * N)))) / n


def c_inf_iv():
    iv.prec = 120
    pi = iv.pi
    return (
        1
        + 2 / pi
        + (2 / pi) * (iv.euler + iv.log(1 + iv.sqrt(2)))
        + (iv.sqrt(2) - 4 / pi) / 4
    )


def tail_gap_iv(N, d, c):
    iv.prec = 120
    N = iv.mpf(N)
    return N ** (iv.mpf(1) / d) - (2 / iv.pi) * iv.log(N) - c - TAIL / N


def exact_gap_iv(N, d):
    iv.prec = 120
    return iv.mpf(N) ** (iv.mpf(1) / d) - rate_iv(N)


def first_tail(e, c_hi, lo, hi):
    f = lambda N: N**e - (2 / math.pi) * math.log(N) - c_hi - TAIL / N
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return hi


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


def gdown(x, d=5):
    k = d - 1 - math.floor(math.log10(x))
    return f"{math.floor(x * 10**k) / 10**k:.{d - 1}e}"


def gup(x, d=3):
    k = d - 1 - math.floor(math.log10(x))
    return f"{math.ceil(x * 10**k) / 10**k:.{d - 1}e}"


def odd_up(N):
    return N if N % 2 else N + 1


# WALL


def wall_at(d, label, lo, hi, mono):
    e = 1 / d
    c = c_inf_iv()
    c_hi = float(c.b)
    na = first_tail(e, c_hi, lo, hi)
    while tail_gap_iv(na, d, c).a <= 0:
        na += 1
    assert tail_gap_iv(na - 1, d, c).a <= 0 or na - 1 < mono
    assert na >= max(mono, 101)
    N = odd_up(na)
    checked = 0
    while True:
        g = exact_gap_iv(N - 2, d)
        if g.a <= 0:
            break
        N -= 2
        checked += 1
    below = exact_gap_iv(N - 2, d)
    top = odd_up(na)
    n0 = N
    for M in range(n0, top + 1, 2):
        assert exact_gap_iv(M, d).a > 0
    margin = iv.mpf(1) / d - iv.log(rate_iv(n0)) / iv.log(n0)
    print(
        f"  {label}: tail bound holds for every base >= {na} (monotone from {mono}), closed form certified at every odd base {n0}..{top} ({(top - n0) // 2 + 1} bases), fails at {n0 - 2} with gap <= {up(below)}"
    )
    print(
        f"  {label}: wall {n0}, bar - alpha_1 >= {down(margin)} there, gap >= {down(exact_gap_iv(n0, d))}"
    )
    return n0


def maynard_alpha(q, consecutive_half):
    L = math.log(q)
    if consecutive_half:
        return math.log((2 + 2 / L) * (2 * q / (q + 1)) * L) / L
    return math.log((1 + 3 / L) * (q / (q - 1)) * L) / L


def maynard_cross(half):
    lo, hi = 10, 10**14
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if maynard_alpha(mid, half) < 0.2:
            hi = mid
        else:
            lo = mid
    return hi


def wall():
    t0 = time.time()
    c = c_inf_iv()
    print(
        f"c_inf = 1 + 2/pi + (2/pi)(gamma + log(1 + sqrt 2)) + (sqrt 2 - 4/pi)/4 <= {up(c, 9)}, tail {TAIL}/base from base 101"
    )
    worst = None
    for N in list(range(101, 3002, 2)) + [94939, 200001, 10**6 + 1, 10**8 + 1]:
        r = rate_iv(N)
        tl = (2 / iv.pi) * iv.log(N) + c + TAIL / N
        assert r.b < tl.a
        gap = tl - r
        worst = gap if worst is None or gap.a < worst.a else worst
    print(
        f"  direct: lambda/fill <= (2/pi) log base + c_inf + {TAIL}/base at every odd base 101..3001 and at 94939, 200001, 10^6 + 1, 10^8 + 1, smallest gap >= {down(worst)}"
    )
    wall_at(5, "bar 1/5", 1000, 10**7, 327)
    wall_at(4, "bar 1/4", 100, 10**7, 43)
    iv.prec = 120
    low = lambda N: (
        (2 * iv.mpf(N) / (iv.pi * (N + 1))) * iv.log(iv.mpf(N + 2) / 5)
        - iv.mpf(1) / 4
        - ((2 / iv.pi) * iv.log(N) - iv.mpf(131) / 100)
    )
    assert all(low(N).a > 0 for N in range(9, 10**4, 2))
    assert ((2 / iv.pi) * iv.log(9) - iv.mpf(131) / 100).a > 0 and (
        (2 / iv.pi) * iv.log(7) - iv.mpf(131) / 100
    ).b < 0
    tail = (
        iv.mpf(106) / 100
        - (2 / iv.pi) * iv.log(5)
        - (2 / iv.pi) * iv.log(iv.mpf(103) / 5) / 102
    )
    assert tail.a > 0
    print(
        f"  lower bound: (2 base/(pi (base+1))) log((base+2)/5) - 1/4 >= (2/pi) log base - 1.31 > 0 at every odd base 9..9999, and from 101 by 1.06 - (2/pi) log 5 - (2/pi) log((base+2)/5)/(base+1) >= {down(tail)}; the factor is negative at 7"
    )
    floor = lambda N: iv.mpf(N) ** (iv.mpf(1) / 5) - 2 * iv.mpf(N) / (N + 1)
    f27 = min(
        N for N in range(3, 200, 2) if all(floor(M).a > 0 for M in range(N, 200, 2))
    )
    assert floor(f27 - 2).b < 0
    print(
        f"  density floor: 1 - alpha_base < 1/5 at every odd base >= {f27}, fails at {f27 - 2}"
    )
    for N in (100003, 10**6 + 3, 10**9 + 7):
        r = rate_iv(N)
        a1 = iv.log(r) / iv.log(N)
        print(
            f"  base {N}: base^alpha_1 = lambda/fill <= {up(r, 8)}, (2/pi) log base = {2 / math.pi * math.log(N):.6f}, alpha_1 <= {up(a1, 7)}"
        )
    one = maynard_cross(False)
    half = maynard_cross(True)
    assert one == 1520573
    print(
        f"  Maynard 2022 alpha_q below 1/5: one missing digit from {one} (calibration), consecutive half interval from {half}"
    )
    print(f"wall in {time.time() - t0:.1f}s")


# CHECK


def grid_sums(N, ts):
    n = (N + 1) // 2
    r = np.arange(N, dtype=np.float64)
    G = np.empty(len(ts))
    S = np.empty(len(ts))
    step = max(1, 4_000_000 // N)
    for i in range(0, len(ts), step):
        t = ts[i : i + step, None]
        u = (t + r[None, :]) / N
        num = np.abs(np.sin(np.pi * n * u))
        den = np.sin(np.pi * u)
        with np.errstate(divide="ignore", invalid="ignore"):
            h = np.where(den > 0, num / np.where(den > 0, den, 1), n)
        G[i : i + step] = h.sum(1)
        S[i : i + step] = num.sum(1)
    return G, S


def check_base(N, T):
    n = (N + 1) // 2
    ts = np.linspace(0.0, 0.5, T)
    G, S = grid_sums(N, ts)
    sig = np.mod(n * ts, 1.0)
    Sc = np.cos(np.pi * (sig - 0.5) / N) / np.sin(np.pi / (2 * N))
    assert np.max(np.abs(S - Sc) / Sc) < 1e-9
    x0 = x0_float(N)
    lam = lam_float(N)
    sn = np.sin(np.pi * ts)
    gb = n + x0 + sn * (x0 / 2 + K_2 * N / (4 * math.pi))
    low = (N / math.pi) * math.log((N + 2) / 5) - n / 4
    r1 = np.max(G / gb)
    r2 = np.max((G + S / 2) / (lam * (1 + sn / 2)))
    r3 = np.min(G) / low if low > 0 else float("inf")
    assert r1 < 1 and r2 < 1 and r3 > 1
    assert K_2 * N / (2 * math.pi) <= n + np.min(S) / 2
    L = math.log(N)
    return (
        r1,
        r2,
        r3,
        G[0] / n - 2 / math.pi * L,
        np.max(G) / n - 2 * math.sqrt(2) / math.pi * L,
        ts[np.argmax(G)],
        np.min(G) / n - 2 / math.pi * L,
    )


def levels_mp(N, i, s):
    mp.dps = 40
    n = (N + 1) // 2
    tot = mp.mpf(0)
    Y = N**i
    for a in range(Y):
        v = mp.mpf(s) + mp.mpf(a) / Y
        p = mp.mpf(1)
        for j in range(i):
            w = v * N**j
            w = w - mp.floor(w)
            if w == 0:
                p *= n
            else:
                p *= abs(mp.sin(mp.pi * n * w) / mp.sin(mp.pi * w))
        tot += p
    return tot


def pieces(N, ts):
    n = (N + 1) // 2
    J = (n - 2) // 2
    K = (n - 1) // 2
    worst = [-1e9, -1e9, -1e9]
    for t in ts:
        m = 2 * np.arange(J + 1) + 1.0
        bp = np.sum(1 / (m + t) + 1 / (m - t)) - (math.log(N + 3) + GAMMA + K_2 * t * t)
        k = np.arange(1, K + 1, dtype=np.float64)
        le = np.sum(1 / (t + 2 * k)) if K else 0.0
        ho = 2 * k - t
        sg = np.sum(np.where(ho / N <= t, 1 / ho, -1 / ho)) if K else 0.0
        ap = le + sg - (math.log(N + 4) - 0.0757)
        dl = (t + 2 * k) / N
        dh = (2 * k - t) / N
        q = (
            np.sum(1 / (2 * np.cos(np.pi * dl / 2)))
            + np.sum(1 / (2 * np.cos(np.pi * dh / 2)))
            if K
            else 0.0
        )
        bq = q - (N / math.pi) * math.log(math.tan(3 * math.pi / 8 + math.pi / (4 * N)))
        worst = [max(worst[0], bp), max(worst[1], ap), max(worst[2], bq)]
    assert max(worst) < 0
    return worst


def check():
    t0 = time.time()
    w = [-1e9] * 3
    for N in list(range(3, 402, 2)) + [1001, 10001, 100001]:
        p = pieces(N, np.linspace(0.0, 0.5, 201))
        w = [max(a, b) for a, b in zip(w, p)]
    print(
        f"pieces at every odd base 3..401 and 1001, 10001, 100001 on 201 shifts, smallest margin under the bound: bP >= {gdown(-w[0])}, aP >= {gdown(-w[1])}, bQ >= {gdown(-w[2])}"
    )
    worst = [0, 0, 1e9]
    for N in range(3, 402, 2):
        r1, r2, r3, *_ = check_base(N, 801)
        worst = [max(worst[0], r1), max(worst[1], r2), min(worst[2], r3)]
    print(
        f"every odd base 3..401 at 801 shifts in [0, 1/2]: max G/bound <= {fup(worst[0])}, max T phi/(lambda phi) <= {fup(worst[1])}, min G/lower >= {fdown(worst[2])}"
    )
    print(
        "base, max G/bound, max T phi/(lambda phi), min G/lower, G(0)/fill - (2/pi) log base, max G/fill - (2 sqrt2/pi) log base, argmax t, min G/fill - (2/pi) log base"
    )
    for N in (1001, 4001, 10001, 30001, 100001):
        r1, r2, r3, g0, gm, tm, gmin = check_base(N, 401)
        print(
            f"  {N}: <= {fup(r1)} <= {fup(r2)} >= {fdown(r3)}, readings {g0:.6f} {gm:.6f} {tm:.4f} {gmin:.6f}"
        )
    print(f"  gamma' + 1 = {(2 / math.pi) * (GAMMA + math.log(8 / math.pi)) + 1:.6f}")
    print(
        "levels at 40 digits, every level 1..top: base, top level, sum at the last shift, (3/2) lambda^top, its ratio"
    )
    worst = 0.0
    for N, top in (
        (3, 8),
        (5, 5),
        (7, 4),
        (9, 4),
        (11, 3),
        (13, 3),
        (15, 3),
        (17, 3),
        (21, 3),
    ):
        lam = lam_float(N)
        for i in range(1, top + 1):
            for s in (0, 0.5, 1 / (2 * N), (math.sqrt(5) - 1) / 2):
                v = levels_mp(N, i, s)
                bound = 1.5 * lam**i
                worst = max(worst, float(v) / bound)
                assert v < bound
            if i == top:
                print(
                    f"  {N} {i} {mp.nstr(v, 10)} {fup(bound, 2)} {fup(float(v) / bound)}"
                )
    print(f"  largest sum/bound over all 4 shifts and levels 1..top: <= {fup(worst)}")
    err = 0.0
    for N in (5, 7, 9):
        n = (N + 1) // 2
        for s in (0, 0.5, 1 / (2 * N), (math.sqrt(5) - 1) / 2):
            t = math.fmod(N * N * s, 1.0)
            u = (t + np.arange(N)) / N
            h = np.array(
                [
                    abs(math.sin(math.pi * n * w) / math.sin(math.pi * w))
                    if w > 0
                    else n
                    for w in u
                ]
            )
            g, _ = grid_sums(N, u)
            two = float(np.sum(h * g))
            err = max(err, abs(two - float(levels_mp(N, 2, s))) / two)
    assert err < 1e-9
    print(
        f"  transfer identity sum_2(s) = (T^2 1)(base^2 s) at bases 5, 7, 9 and 4 shifts, largest relative error <= {gup(err)}"
    )
    print(f"check in {time.time() - t0:.1f}s")


# RATE


def power(N, M, it):
    n = (N + 1) // 2
    t = (np.arange(M) + 0.5) / M
    r = np.arange(N)
    phi = np.ones(M)
    lo = hi = 0.0
    for _ in range(it):
        new = np.zeros(M)
        for k in range(0, M, max(1, 2_000_000 // N)):
            u = (t[k : k + max(1, 2_000_000 // N), None] + r[None, :]) / N
            h = np.abs(np.sin(np.pi * n * u) / np.sin(np.pi * u))
            new[k : k + u.shape[0]] = (
                h * np.interp(u.ravel(), t, phi, period=1.0).reshape(u.shape)
            ).sum(1)
        q = new / phi
        lo, hi = q.min(), q.max()
        phi = new / new.max()
    return lo / n, hi / n, phi.min()


def rate():
    t0 = time.time()
    print(
        "power iteration on 1000 cells, readings: base, lambda/fill - (2/pi) log base low/high, min phi/max phi, base^(1/5) - lambda/fill"
    )
    for N in (101, 1001, 10001):
        lo, hi, pm = power(N, 1000, 30)
        print(
            f"  {N}: {lo - 2 / math.pi * math.log(N):.5f} {hi - 2 / math.pi * math.log(N):.5f} {pm:.4f} {N**0.2 - hi:.4f}"
        )
    print("power iteration on 300 cells near the route's own crossing")
    for N in (60001, 70001, 80001):
        lo, hi, pm = power(N, 300, 12)
        print(
            f"  {N}: {lo - 2 / math.pi * math.log(N):.5f} {hi - 2 / math.pi * math.log(N):.5f} {pm:.4f} {N**0.2 - hi:.4f}"
        )
    one = lambda N: grid_sums(N, np.array([0.5]))[0][0] / ((N + 1) // 2)
    for e in (0.25, 0.2):
        lo, hi = 1001, 2_000_001
        while hi - lo > 2:
            mid = odd_up((lo + hi) // 2)
            if one(mid) < mid**e:
                hi = mid
            else:
                lo = mid
        print(
            f"  one-step reading G(1/2)/fill < base^{e} first at odd base {hi}, G(1/2)/fill - (2 sqrt2/pi) log base = {one(hi) - 2 * math.sqrt(2) / math.pi * math.log(hi):.5f} there"
        )
    print(f"rate in {time.time() - t0:.1f}s")


# METER


def meter():
    t0 = time.time()
    X = 10**7
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
    k = np.arange(X + 1, dtype=np.int64)
    logk = np.log(np.maximum(k, 1))
    print(
        "sanity only: base, x, A(x), M(x), max |M|/A, max |M|/sqrt(A), primes: sum log p/(kappa A) at x, kappa = p/(p+1)"
    )
    for N in (101, 1009, 10007):
        h = (N - 1) // 2
        ok = np.ones(X + 1, dtype=bool)
        y = k.copy()
        while y.any():
            ok &= (y % N) <= h
            y //= N
        ok[0] = False
        A = np.cumsum(ok)
        M = np.cumsum(np.where(ok, mu, 0).astype(np.int64))
        sel = A >= 100
        th = np.sum(logk[ok & isp])
        kap = N / (N + 1)
        print(
            f"  {N}: {X} {A[-1]} {M[-1]} {np.max(np.abs(M[sel]) / A[sel]):.4f} {np.max(np.abs(M[sel]) / np.sqrt(A[sel])):.3f} {th / (kap * A[-1]):.4f}"
        )
    print(f"meter in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "wall"
    {"wall": wall, "check": check, "rate": rate, "meter": meter}[verb]()
