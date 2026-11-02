import math
import sys
import time
from fractions import Fraction

import numpy as np
from mpmath import iv, mp, mpf

# SET


def member(n, b):
    n = np.asarray(n, dtype=np.int64)
    ok = n >= 1
    prev = n % b
    m = n // b
    while (m > 0).any():
        d = m % b
        ok &= ~((m > 0) & (d == prev))
        prev = np.where(m > 0, d, prev)
        m = m // b
    return ok


def automaton(b, L):
    vals, tops, out = np.arange(b, dtype=np.int64), np.arange(b, dtype=np.int64), []
    out.append(vals[tops != 0])
    for l in range(1, L):
        nv, nt = [], []
        for d in range(b):
            keep = tops != d
            nv.append(vals[keep] + d * b**l)
            nt.append(np.full(keep.sum(), d, dtype=np.int64))
        vals, tops = np.concatenate(nv), np.concatenate(nt)
        out.append(vals[tops != 0])
    return np.sort(np.concatenate(out))


def strings(b, L):
    n = np.arange(b**L, dtype=np.int64)
    ok = np.ones(b**L, dtype=bool)
    prev, m = n % b, n // b
    for _ in range(1, L):
        d = m % b
        ok &= d != prev
        prev, m = d, m // b
    return ok, prev


def digits(x, b):
    out = []
    while x:
        out.append(x % b)
        x //= b
    return out


def blocks(x, b):
    ds = digits(x, b)
    L = len(ds)
    out = [(f, j, f) for j in range(L - 1) for f in range(1, b)]
    for j in range(L - 1, -1, -1):
        hi = ds[j + 1 :]
        if any(hi[i] == hi[i + 1] for i in range(len(hi) - 1)):
            break
        P = sum(d * b**i for i, d in enumerate(hi))
        for f in range(ds[j]):
            if (j == L - 1 and f == 0) or (j < L - 1 and f == ds[j + 1]):
                continue
            out.append((P * b + f, j, f))
    return out


def count(x, b):
    tot = sum(1 if j == 0 else (b - 1) ** j for _, j, _ in blocks(x, b))
    return tot + int(member(np.array([x]), b)[0])


def units(b):
    return [c for c in range(b) if math.gcd(c, b) == 1]


# SPECTRAL


def transform(b, L, u, w=None):
    ok, top = strings(b, L)
    n = np.nonzero(ok if w is None else ok & (w[top] != 0))[0]
    return np.exp(2j * np.pi * np.outer(u, n)).sum(axis=1)


def run_expansion(b, i, u):
    if i == 0:
        return np.ones_like(u, dtype=complex)
    tot = np.zeros_like(u, dtype=complex)
    for r in range(1, i + 1):
        rep = (b**r - 1) // (b - 1)
        S = sum(np.exp(2j * np.pi * c * rep * u) for c in range(b))
        tot += (-1) ** (r - 1) * S * run_expansion(b, i - r, (b**r * u) % 1.0)
    return tot


def grid_sum(ind, s):
    N = len(ind)
    return np.abs(np.fft.fft(ind * np.exp(2j * np.pi * np.arange(N) * s))).sum()


def absS(m, v):
    v = np.asarray(v, dtype=float) % 1.0
    den = np.abs(np.sin(np.pi * v))
    small = den < 1e-13
    return np.where(
        small, float(m), np.abs(np.sin(np.pi * m * v)) / np.where(small, 1.0, den)
    )


def matrix(b, u):
    M = np.zeros((b + 1, b + 1), dtype=complex)
    ph = np.exp(2j * np.pi * np.arange(b) * u)
    M[b, :b] = ph
    for c in range(b):
        M[c, :b] = ph
        M[c, c] = 0
    return M


# CERTIFICATE

SI_PI = "1.85193705198246617036105337015799136334580972898115490980477"
EULER = "0.57721566490153286060651209008240243104215933593992359880577"


def ivc(s):
    return iv.mpf([mpf(s) - mpf(10) ** -50, mpf(s) + mpf(10) ** -50])


def hstar(cells=20000):
    mp.dps = 40
    g = mpf(EULER)
    best = mp.pi + mp.sin(mp.pi / cells) * (
        1 / (1 - mpf(1) / cells) + g - mp.digamma(2 - mpf(1) / cells)
    )
    for k in range(1, cells // 2):
        t1, t2 = mpf(k) / cells, mpf(k + 1) / cells
        best = max(
            best,
            mp.sin(mp.pi * t1) / t1
            + mp.sin(mp.pi * t2)
            * (1 / (1 - t2) - mp.digamma(1 + t1) - mp.digamma(2 - t2)),
        )
    return best


def constants(hs):
    iv.prec = 120
    pi, g, si = iv.pi, ivc(EULER), ivc(SI_PI)
    c0 = (hs - 2 * iv.log(2)) / pi + (1 - 2 / pi) / 2
    ce = 1 + 2 / pi
    return pi, g, si, c0, ce


def Qform(b, y, C):
    pi, g, si, c0, ce = C
    b = iv.mpf(b)
    lb = iv.log(b)
    E = 2 * iv.log(b / 2) + 2 * g - 1 + 2 / (b - 1) + 2 / b
    K = E / pi + (1 - 2 / pi) * (iv.mpf(1) / 2 + 1 / b)
    P = (
        y
        * y
        * (
            (2 * si / pi + 2 * K / pi) / (1 - y)
            + (2 * (4 - pi) / pi) / (b - y)
            + (1 - 2 / pi) * K / (b * b - y)
        )
    )
    Ph = (
        2 * lb * y / (pi * (1 - y) ** 2) + c0 * y / (1 - y) + ce * (y / b) / (1 - y / b)
    )
    return P * (1 + Ph), Ph


def ystar(b, bar):
    b = iv.mpf(b)
    return iv.exp((1 - iv.mpf(bar)) * iv.log(b)) / (b - 1)


def ycert(b, C, tol=1e-10):
    lo, hi = 1e-6, 0.9
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if Qform(b, iv.mpf(mid), C)[0].b <= 1:
            lo = mid
        else:
            hi = mid
    return lo


def alpha_up(b, C):
    y = ycert(b, C)
    a = iv.log(iv.mpf(b) / ((b - 1) * iv.mpf(y))) / iv.log(iv.mpf(b))
    return float(a.b), y


def Ufull(N, hs):
    return (2 * math.log(N / 2 + 1) + hs) / math.pi + (1 - 2 / math.pi) * (0.5 + 1 / N)


def Rrun(b, r):
    M = b ** (r - 1)
    g = float(mpf(EULER))
    E = 2 * math.log(b / 2) + 2 * g - 1 + 2 / (b - 1) + 2 / b
    K = E / math.pi + (1 - 2 / math.pi) * (0.5 + 1 / b)
    return (2 * float(mpf(SI_PI)) + 2 * (4 - math.pi) / M) / math.pi + K * (
        2 / math.pi + (1 - 2 / math.pi) / M**2
    )


def tail(B1, bar, C):
    pi, g, si, c0, ce = C
    B = iv.mpf(B1)
    e1 = B / (B - 1)
    ymax = e1 * iv.exp(-iv.mpf(bar) * iv.log(B))
    k0 = (-2 * iv.log(2) + 2 * g - 1 + iv.mpf(401) / (100 * B)) / pi + (1 - 2 / pi) * (
        iv.mpf(1) / 2 + 1 / B
    )
    p1 = (2 / pi) * (2 / pi) / (1 - ymax) + (1 - 2 / pi) * (2 / pi) / (B * B - 1)
    p0 = (
        (2 * si / pi + 2 * k0 / pi) / (1 - ymax)
        + (2 * (4 - pi) / pi) / (B - 1)
        + (1 - 2 / pi) * k0 / (B * B - 1)
    )
    f1 = (2 / pi) / (1 - ymax) ** 2
    f0 = c0 / (1 - ymax) + ce / (B - 1)
    lb = iv.log(B)
    a = iv.mpf(bar)
    gB = (
        e1
        * e1
        * iv.exp(-2 * a * lb)
        * (p1 * lb + p0)
        * (1 + e1 * iv.exp(-a * lb) * (f1 * lb + f0))
    )
    turn = max(float((1 / (2 * a) - p0 / p1).b), float((1 / a - f0 / f1).b))
    return gB, turn


# VERBS


def count_verb():
    t0 = time.time()
    for b in (3, 4, 5, 10):
        X = b**6
        brute = np.nonzero(member(np.arange(X), b))[0]
        auto = automaton(b, 6)
        auto = auto[auto >= 1]
        assert np.array_equal(brute, auto)
        ok, _ = strings(b, 6)
        exact = [
            int(((brute >= b ** (L - 1)) & (brute < b**L)).sum()) for L in range(1, 7)
        ]
        assert exact == [(b - 1) ** L for L in range(1, 7)]
        assert ok.sum() == b * (b - 1) ** 5
        print(
            f"b {b}: automaton = brute force below b^6, {len(brute)} members, {ok.sum()} = b (b-1)^5 strings of length 6, (b-1)^L with L digits"
        )
    rng = np.random.default_rng(7)
    nx = 0
    for b in (3, 4, 5, 7, 10):
        for x in list(rng.integers(2, 200000, 80)) + [
            b**5 - 1,
            b**5,
            b**5 + 1,
            2 * b**4 + b,
        ]:
            x = int(x)
            got = []
            per = {}
            for P, j, f in blocks(x, b):
                ok, top = strings(b, j) if j > 0 else (np.array([True]), np.array([-1]))
                t = np.nonzero(ok & (top != f))[0]
                got.append(P * b**j + t)
                per[j] = per.get(j, 0) + 1
                assert len(t) == (b - 1) ** j
            got = np.concatenate(got)
            got = np.sort(
                np.concatenate([got, [x]] if member(np.array([x]), b)[0] else [got])
            )
            want = np.nonzero(member(np.arange(1, x + 1), b))[0] + 1
            assert np.array_equal(got, want) and max(per.values()) <= 2 * (b - 1)
            assert count(x, b) == len(want) and len(want) >= (b - 1) ** (
                len(digits(x, b)) - 1
            )
            nx += 1
    print(
        f"blocks: the split matches direct membership at {nx} values of x at b = 3, 4, 5, 7, 10, at most 2(b-1) blocks per scale, each of (b-1)^j members"
    )
    cells = 0
    for b in range(3, 13):
        U, phi = units(b), len(units(b))
        for j in range(1, 7):
            ok, top = strings(b, j)
            n = np.nonzero(ok)[0]
            u0 = np.isin(n % b, U)
            for f in range(b):
                c = int((u0 & (top[n] != f)).sum())
                assert b * c == phi * (b - 1) ** j + (-1) ** (j - 1) * (
                    phi - b * (f in U)
                )
                cells += 1
    print(
        f"units digit: b #(unit units digit, top != f) = phi(b)(b-1)^j + (-1)^(j-1)(phi(b) - b [f unit]) at all {cells} cells b = 3..12, j = 1..6, every f"
    )
    for b in range(3, 40):
        assert member(np.array([1, b]), b).all() and not member(np.array([b + 1]), b)[0]
    J = lambda b: np.ones((b, b), dtype=np.int64) - np.eye(b, dtype=np.int64)
    for b in range(3, 9):
        for n in range(9):
            A = np.linalg.matrix_power(J(b), n)
            want = ((b - 1) ** n - (-1) ** n) // b * np.ones((b, b), dtype=np.int64) + (
                -1
            ) ** n * np.eye(b, dtype=np.int64)
            assert np.array_equal(A, want)
    print(
        "witness 1, b in R_b and b + 1 not, at every b = 3..39; (J - I)^n closed form at b = 3..8, n = 0..8"
    )
    rng = np.random.default_rng(5)
    worst = 0
    for b in (4, 6, 8, 10):
        h, q = b // 2, b * b
        F = np.array(
            [a + b * be for a in range(h, b) for be in range(h)], dtype=np.int64
        )
        S = F.copy()
        for _ in range(2):
            S = (S[:, None] * q + F[None, :]).ravel()
        assert len(S) == h**6 and member(S, b).all()
        assert h in F and h + 1 in F and 0 not in F
        assert b * sum(math.gcd(int(f), q) == 1 for f in F) == len(units(b)) * h * h
        for i in (1, 2):
            t = rng.random(20)
            Fi = F.copy()
            for _ in range(i - 1):
                Fi = (Fi[:, None] * q + F[None, :]).ravel()
            lhs = np.abs(np.exp(2j * np.pi * np.outer(t, Fi)).sum(axis=1))
            Hi = np.arange(h ** (2 * i))
            Hi = sum(((Hi // h**m) % h) * b**m for m in range(2 * i))
            rhs = np.abs(np.exp(2j * np.pi * np.outer(t, Hi)).sum(axis=1))
            worst = max(worst, np.abs(lhs - rhs).max())
    print(
        f"embedded design at b = 4, 6, 8, 10: the {{a + b beta : b/2 <= a < b, beta < b/2}}-strings of three digits at base b^2 all lie in R_b, keep b/2 and b/2 + 1, kappa_F = 1 exactly, and abs(hat F_i) = abs(hat H_(2i)) at i = 1, 2 to {worst:.1e}"
    )
    print(f"count: {time.time() - t0:.1f} s")


def identity_verb():
    t0 = time.time()
    rng = np.random.default_rng(3)
    worst = 0
    for b, Lmax in ((3, 7), (4, 6), (5, 5), (7, 4), (10, 4)):
        u = rng.random(9)
        for i in range(1, Lmax + 1):
            worst = max(
                worst, np.abs(transform(b, i, u) - run_expansion(b, i, u)).max()
            )
            prod = np.array(
                [
                    np.linalg.multi_dot(
                        [np.eye(b + 1)[b]]
                        + [matrix(b, (b**k * x) % 1.0) for k in range(i)]
                        + [np.r_[np.ones(b), 0]]
                    )
                    for x in u
                ]
            )
            worst = max(worst, np.abs(prod - transform(b, i, u)).max())
    print(
        f"run expansion and matrix product against the direct transform: max error {worst:.2e}"
    )
    print(f"identity: {time.time() - t0:.1f} s")


def wall_verb():
    t0 = time.time()
    hs_raw = hstar()
    hs_up = math.ceil(float(hs_raw) * 10**6 + 1e-6) / 10**6
    print(
        f"h* = max of sin(pi t)(-psi(t) - psi(1-t)) on [0, 1/2] <= {float(hs_raw):.9f}, printed {hs_up:.6f}; value at t = 1/2: {float(2 * mpf(EULER) + 4 * mp.log(2)):.9f}"
    )
    C = constants(iv.mpf(hs_up))
    res = {}
    for bar in (0.2, 0.25):
        lo, hi = (1000, 2400) if bar == 0.2 else (150, 600)
        fails, passes = [], []
        for b in range(lo, hi + 1):
            q, _ = Qform(b, ystar(b, bar), C)
            if q.b < 1:
                passes.append(b)
            elif q.a > 1:
                fails.append(b)
            else:
                raise SystemExit(f"undecided at {b}")
        B0 = max(fails) + 1
        assert all(b in passes for b in range(B0, hi + 1))
        gB, turn = tail(hi, bar, C)
        assert gB.b < 1 and turn < math.log(hi)
        qB0, _ = Qform(B0, ystar(B0, bar), C)
        qF, _ = Qform(B0 - 1, ystar(B0 - 1, bar), C)
        a0, y0 = alpha_up(B0, C)
        res[bar] = B0
        print(
            f"bar {bar}: Q(y*) < 1 at every b = {B0}..{hi}, {hi - B0 + 1} bases, Q = {float(qB0.b):.7f} at {B0}; fails at {B0 - 1}, Q >= {float(qF.a):.7f}; tail g({hi}) <= {float(gB.b):.6f}, g decreasing from log b = {turn:.3f}; alpha_1 <= {a0:.7f} at {B0}, gap bar - alpha_1 >= {bar - a0:.4e}"
        )
    B0 = res[0.2]
    y0 = ycert(B0, C)
    q, Ph = Qform(B0, iv.mpf(y0), C)
    C1 = math.ceil(float((1 + Ph).b) * 10**4 + 1e-9) / 10**4
    Ce = math.ceil((C1 * y0 / (1 - y0) + 1) * 10**4 + 1e-9) / 10**4
    print(f"at b = {B0}: y = {y0:.8f}, C_1 = 1 + Phi(y) <= {C1:.4f}, C_e <= {Ce:.4f}")
    print("spread, alpha_1 upper end of the certificate:")
    for b in (
        3,
        4,
        5,
        7,
        10,
        31,
        100,
        1000,
        B0,
        10**4,
        10**5,
        10**6,
        10**7,
        10**8,
        10**9,
    ):
        a, y = alpha_up(b, C)
        print(
            f"  b {b}: alpha_1 <= {a:.7f}, y = {y:.6f}, z/(b-1) = {b / ((b - 1) * y):.5f}, Cauchy-Schwarz exponent 1 - log(b-1)/(2 log b) = {1 - math.log(b - 1) / (2 * math.log(b)):.4f}"
        )
    for cut, lo, hi in ((0.22, 300, 1000), (0.18, 2500, 5000)):
        ok = [Qform(b, ystar(b, cut), C)[0].b < 1 for b in range(lo, hi + 1)]
        first = lo + max(i for i, v in enumerate(ok) if not v) + 1
        assert all(ok[first - lo :])
        print(
            f"  alpha_1 < {cut} at every b = {first}..{hi} and not at {first - 1}: {alpha_up(first - 1, C)[0]:.7f} at {first - 1}, {alpha_up(first, C)[0]:.7f} at {first}"
        )
    print(f"wall: {time.time() - t0:.1f} s")


def check_verb():
    t0 = time.time()
    hs = 3.927267
    worst = 0
    for b in list(range(3, 41)) + [64, 101, 257, 1000, 1009]:
        for l in (1, 2):
            N = b**l
            if N > 200000:
                continue
            q = np.arange(N)
            best = max(absS(N, (t + q) / N).sum() for t in np.linspace(0, 0.5, 201)) / N
            worst = max(worst, best / Ufull(N, hs))
    print(
        f"full blocks: sup over 201 shifts of the grid sum of S_N over N = b^l, against U(N): at most {worst:.6f} of it"
    )
    worst = 0
    for b in list(range(3, 41)) + [64, 101, 257]:
        for r in (2, 3):
            M = b ** (r - 1)
            if b**r > 400000:
                continue
            j = np.arange(b)
            best = 0
            for t in np.linspace(0, 1, 101, endpoint=False):
                tm = (t + np.arange(M)) / M
                best = max(best, absS(b, (tm[:, None] + j[None, :]) / b).sum())
            worst = max(worst, best / b**r / Rrun(b, r))
    print(
        f"long runs: sup over 100 shifts of the grid sum of S_b over b^r points, against R_r: at most {worst:.6f} of it"
    )
    b = 1000
    q = np.arange(b)
    f1 = max(absS(b, (t + q) / b).sum() for t in np.linspace(0, 0.5, 201)) / b
    j = np.arange(b)
    r2 = (
        max(
            absS(b, ((t + np.arange(b))[:, None] / b + j[None, :]) / b).sum()
            for t in np.linspace(0, 1, 21, endpoint=False)
        )
        / b**2
    )
    print(
        f"at b = 1000: full block sup {f1:.4f} against U = {Ufull(b, hs):.4f}; long run sup {r2:.4f} against R_2 = {Rrun(b, 2):.4f}"
    )
    C = constants(iv.mpf(hs))
    worst1, worst2 = 0, 0
    for b, Lmax in ((3, 12), (4, 9), (5, 8), (7, 7), (10, 6)):
        y = ycert(b, C)
        z = b / y
        C1 = 1 + float(Qform(b, iv.mpf(y), C)[1].b) + 1e-9
        Ce = C1 * y / (1 - y) + 1
        for i in range(1, Lmax + 1):
            ok, top = strings(b, i)
            for s in (0.0, 0.37 / b**i, 0.5 / b**i, 0.81 / b**i):
                worst1 = max(worst1, grid_sum(ok.astype(float), s) / (C1 * z**i))
                for f in (0, b - 1):
                    worst2 = max(
                        worst2,
                        grid_sum((ok & (top != f)).astype(float), s)
                        / ((C1 + Ce) * z**i),
                    )
    print(
        f"grid sums by FFT at b = 3, 4, 5, 7, 10, every level with b^i <= 10^6, four shifts: Sigma_i(s) at most {worst1:.6f} of C_1 z^i, the blocks at most {worst2:.6f} of (C_1 + C_e) z^i"
    )
    rng = np.random.default_rng(11)
    worst = 0
    for b in (5, 10, 31, 101):
        th = 3 / (8 * b)
        sig = math.sin(3 * math.pi / 8) / math.sin(3 * math.pi / (8 * b))
        for _ in range(400):
            v = th + rng.random() * (1 - 2 * th)
            P2 = np.abs(matrix(b, v) @ matrix(b, rng.random()))
            kap = P2.sum(axis=1).max()
            worst = max(worst, kap / (b * min(sig, absS(b, v)) + 2 * b - 1))
        c2 = (b * sig + 2 * b - 1) / (b - 1) ** 2
        print(f"  b {b}: c_2 = (b sigma + 2b - 1)/(b-1)^2 = {c2:.6f}")
    print(
        f"two-step rows: kappa at most {worst:.6f} of b abs(S_b(v)) + 2b - 1 at 1600 random pairs with dist(v, Z) >= 3/(8b)"
    )
    cE = [
        (b * math.sin(3 * math.pi / 8) / math.sin(3 * math.pi / (8 * b)) + 2 * b - 1)
        / (b - 1) ** 2
        for b in range(3, 5001)
    ]
    first = 3 + max(i for i, v in enumerate(cE) if v >= 1) + 1
    assert all(v < 1 for v in cE[first - 3 :]) and all(
        cE[i + 1] < cE[i] for i in range(first - 3, len(cE) - 1)
    )
    print(
        f"c_E < 1 and decreasing at every b = {first}..5000, {cE[first - 4]:.6f} at {first - 1}; c_E = {cE[1153 - 3]:.6f} at 1153, limit (8/(3 pi)) sin(3 pi/8) = {8 / (3 * math.pi) * math.sin(3 * math.pi / 8):.6f}"
    )
    worst, cnt = 0, 0
    b, k = 31, 12
    c2 = (
        b * math.sin(3 * math.pi / 8) / math.sin(3 * math.pi / (8 * b)) + 2 * b - 1
    ) / (b - 1) ** 2
    den = 4 * b * (b - 1) * b ** math.ceil(2 * k / 3) * 1000
    for d in (2 * 3, 7, 9, 11, 13, 21, 62, 93, 99, 151, 961 * 3):
        md = max(1, math.floor(math.log(d / 2, b)) + 1)
        for l in range(1, d, max(1, d // 60)):
            if math.gcd(l, d) != 1 or d // math.gcd(d, b**k) == 1:
                continue
            eta = Fraction(int(rng.integers(-999, 1000)), den)
            u = Fraction(l, d) + eta
            v = np.eye(b + 1)[b]
            for j in range(k):
                v = v @ matrix(b, float((b**j * u) % 1))
            val = np.abs(v[:b] @ (np.ones(b) - np.eye(b)[0]))
            npairs = math.floor(2 * k / (3 * md)) // 2
            worst = max(worst, val / (b * (b - 1) ** (k - 1) * c2**npairs))
            cnt += 1
    print(
        f"C1 at b = 31, level 12, {cnt} reduced l/d over 11 moduli with a prime outside b, each plus a seeded rational perturbation: abs(hat W_k) at most {worst:.3e} of b (b-1)^(k-1) c_2^(pairs)"
    )
    print(f"check: {time.time() - t0:.1f} s")


def rate_verb():
    t0 = time.time()
    for b, Lmax in ((3, 13), (4, 10), (5, 8), (7, 7), (10, 6), (31, 4), (100, 3)):
        out, prev = [], None
        for L in range(1, Lmax + 1):
            ok, _ = strings(b, L)
            S0 = np.abs(np.fft.fft(ok.astype(float))).sum()
            if prev:
                out.append(S0 / prev / (b - 1))
            prev = S0
        print(
            f"b {b}: Sigma_L(0)/Sigma_(L-1)(0)/(b-1) at L = 2..{Lmax}: "
            + ", ".join(f"{x:.4f}" for x in out)
        )
    print(f"rate: {time.time() - t0:.1f} s")


def primes_upto(X):
    s = np.ones(X // 2 + 1, dtype=bool)
    s[0] = False
    for i in range(1, math.isqrt(X) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[p * p // 2 :: p] = False
    out = 2 * np.nonzero(s)[0] + 1
    return np.concatenate([[2], out[out <= X]])


def primes_verb():
    t0 = time.time()
    X = 10**8
    P = primes_upto(X)
    assert len(P) == 5761455
    small = P[P <= 10**4]
    for b in (10, 7, 3, 4, 5):
        pw = []
        for p in small:
            q = int(p) * int(p)
            while q <= X:
                pw.append((q, math.log(int(p))))
                q *= int(p)
        line = []
        for x in (10**6, 10**7, X):
            Px = P[P <= x]
            keep = member(Px, b)
            psi = np.log(Px[keep].astype(float)).sum() + sum(
                lg for q, lg in pw if q <= x and member(np.array([q]), b)[0]
            )
            A = count(x, b)
            line.append(f"{psi / A:.6f}")
        print(
            f"b {b}: psi_R/A at 10^6, 10^7, 10^8: "
            + ", ".join(line)
            + f"; A(10^8) = {count(X, b)}"
        )
    seg = 5 * 10**6
    tot, A, worst, mertens = 0, 0, 0.0, 0
    for lo in range(1, X + 1, seg):
        hi = min(lo + seg, X + 1)
        n = np.arange(lo, hi, dtype=np.int64)
        mu = np.ones(hi - lo, dtype=np.int8)
        rem = n.copy()
        for p in small:
            p = int(p)
            st = (-lo) % p
            mu[st::p] *= -1
            rem[st::p] //= p
            st2 = (-lo) % (p * p)
            mu[st2 :: p * p] = 0
        mu[rem > 1] *= -1
        mertens += int(mu.astype(np.int64).sum())
        keep = member(n, 10)
        Mc = tot + np.cumsum(np.where(keep, mu, 0).astype(np.int64))
        Ac = A + np.cumsum(keep)
        tot, A = int(Mc[-1]), int(Ac[-1])
        worst = max(worst, float(np.abs(Mc).max()))
    assert A == count(X, 10) and mertens == 1928
    print(
        f"b 10: M(10^8) = {tot}, A(10^8) = {A}, abs M(10^8)/A^(1/2) = {abs(tot) / math.sqrt(A):.3f}, max_(n <= 10^8) abs M(n)/A(10^8)^(1/2) = {worst / math.sqrt(A):.3f}"
    )
    print(f"primes: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "count"
    {
        "count": count_verb,
        "identity": identity_verb,
        "wall": wall_verb,
        "check": check_verb,
        "rate": rate_verb,
        "primes": primes_verb,
    }[verb]()
