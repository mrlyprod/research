import time
from fractions import Fraction as Q
from itertools import pairwise
from math import ceil, floor, gcd, log10, sqrt

import numpy as np

WINDOWS = 32
CODE = 7
LADDER = [(21, 23), (101, 103), (2321, 2323), (23231, 23233), (232321, 232323)]
GAPS = [(2321, 2323), (2320, 2322), (2321, 2322), (2321, 2325), (2319, 2325)]
LEVEL2 = [201, 1001, 2001, 200, 1000, 2000]
FINE = 1024
FINE_PAIR = (232321, 232323)

# DESIGNS


def up(x, digits=4):
    scale = 10 ** (digits - 1 - floor(log10(x)))
    return ceil(x * scale) / scale


def corners(code):
    return [[(code >> (2 * a + b)) & 1 for b in (0, 1)] for a in (0, 1)]


def walsh(code):
    t = corners(code)
    tau = [[1 - 2 * t[a][b] for b in (0, 1)] for a in (0, 1)]
    out = {}
    for sa in (0, 1):
        for sb in (0, 1):
            out[(sa, sb)] = Q(
                sum(
                    tau[a][b] * (-1) ** (sa * a + sb * b)
                    for a in (0, 1)
                    for b in (0, 1)
                ),
                4,
            )
    return out


def fill(code):
    return Q(code.bit_count(), 4)


# TRIANGLE WAVE


def dist2(x):
    r = x - 2 * (x // 2)
    return r if r <= 1 else 2 - r


def tri(x):
    return 1 - 2 * dist2(x)


def dist2_integral(t):
    whole = t // 2
    r = t - 2 * whole
    part = r * r / 2 if r <= 1 else 2 * r - r * r / 2 - 1
    return whole + part


def tri_mean(h, a, b):
    return 1 - 2 * (dist2_integral(h * b) - dist2_integral(h * a)) / (h * (b - a))


# LEVEL ONE


def law1(n, m, w):
    big = n * m * w
    cuts = np.unique(
        np.concatenate(
            [
                np.arange(0, big + 1, m * w),
                np.arange(0, big + 1, n * w),
                np.arange(0, big + 1, n * m),
            ]
        )
    )
    left, size = cuts[:-1], np.diff(cuts)
    state = (left // (m * w)) % 2 + 2 * ((left // (n * w)) % 2)
    window = left // (n * m)
    counts = np.bincount(
        window * 4 + state, weights=size.astype(np.float64), minlength=4 * w
    )
    return counts.reshape(w, 4) / float(n * m), np.rint(counts).astype(
        np.int64
    ).reshape(w, 4)


def kernel1(code):
    t = corners(code)
    k = np.zeros((4, 4))
    for sx in range(4):
        for sy in range(4):
            k[sx, sy] = float(t[sx & 1][sy & 1] != t[sx >> 1][sy >> 1])
    return k


def lag_law(rho):
    same, diff = (1 + rho) / 4, (1 - rho) / 4
    return [same, diff, diff, same]


def limit_law1(h, w):
    rows = []
    for i in range(w):
        rows.append([float(x) for x in lag_law(tri_mean(h, Q(i, w), Q(i + 1, w)))])
    return np.array(rows)


def walsh_limit(code, rs, rt):
    c = walsh(code)
    return (
        1
        - c[(0, 0)] ** 2
        - c[(1, 0)] ** 2 * rs
        - c[(0, 1)] ** 2 * rt
        - c[(1, 1)] ** 2 * rs * rt
    ) / 2


def kernel_limit(code, rs, rt):
    t = corners(code)
    ps, pt = lag_law(rs), lag_law(rt)
    return sum(
        ps[sx] * pt[sy]
        for sx in range(4)
        for sy in range(4)
        if t[sx & 1][sy & 1] != t[sx >> 1][sy >> 1]
    )


def code7_h(u, v):
    return (1 - abs(1 - 2 * u) * abs(1 - 2 * v)) / 2


def check_walsh():
    points = [Q(k, 12) for k in range(13)]
    bad = 0
    for code in range(16):
        for u in points:
            for v in points:
                for h in (1, 2, 3, 4, 6):
                    if walsh_limit(code, tri(h * u), tri(h * v)) != kernel_limit(
                        code, tri(h * u), tri(h * v)
                    ):
                        bad += 1
    classes = {}
    for code in range(16):
        c = walsh(code)
        key = (c[(1, 0)] ** 2, c[(0, 1)] ** 2, c[(1, 1)] ** 2)
        classes.setdefault(key, []).append(code)
    exact7 = all(
        walsh_limit(CODE, tri(2 * u), tri(2 * v)) == code7_h(u, v)
        for u in points
        for v in points
    )
    means = all(
        walsh_limit(code, tri_mean(h, Q(0), Q(1)), tri_mean(h, Q(0), Q(1)))
        == 2 * fill(code) * (1 - fill(code))
        for code in range(16)
        for h in range(1, 9)
    )
    print(
        f"walsh form against the lag kernel: {bad} mismatches over 16 codes, gaps 1 2 3 4 6, 169 points"
    )
    print(f"code 7 at gap 2 is (1 - |1-2u| |1-2v|)/2 on all 169 points: {exact7}")
    print(
        f"mean over the square is 2 p (1 - p) for all 16 codes and gaps 1..8: {means}"
    )
    for key, codes in sorted(classes.items()):
        print(
            f"  squared walsh (first, second, both) = {tuple(str(x) for x in key)}: codes {codes}"
        )


def direct2d(code, n, m, w):
    side = n * m
    i = np.arange(side)
    a, b = (i // m) % 2, (i // n) % 2
    t = np.array(corners(code))
    ink_n = t[a[:, None], a[None, :]]
    ink_m = t[b[:, None], b[None, :]]
    x = (ink_n != ink_m).astype(np.int64)
    step = side // w
    return x.reshape(w, step, w, step).sum(axis=(1, 3)), step * step


def check_direct():
    n, m, w = 21, 23, 7
    bad = 0
    for code in range(16):
        counts, cells = direct2d(code, n, m, w)
        _, exact = law1(n, m, w)
        k = kernel1(code).astype(np.int64)
        fact = exact @ k @ exact.T
        bad += int(np.sum(fact != counts * (n * m) ** 2 // cells))
    print(
        f"2d raster of side 483 against the factorised box mean, 16 codes, 7 by 7 windows: {bad} mismatches"
    )


def ladder():
    print(
        f"code 7, gap 2, {WINDOWS} by {WINDOWS} windows of side 1/{WINDOWS}, error = max |box mean - H(centre)|"
    )
    k = kernel1(CODE)
    us = (np.arange(WINDOWS) + 0.5) / WINDOWS
    hc = (1 - np.abs(1 - 2 * us)[:, None] * np.abs(1 - 2 * us)[None, :]) / 2
    lim = limit_law1(2, WINDOWS)
    exact_centre = np.max(np.abs(lim @ k @ lim.T - hc))
    print(f"  box average of H equals H at the box centre: {exact_centre:.1e}")
    for n, m in LADDER:
        t0 = time.time()
        law, _ = law1(n, m, WINDOWS)
        means = law @ k @ law.T
        err = np.max(np.abs(means - hc))
        bound = (20 * (m - n) + 36 * WINDOWS) / n
        print(
            f"  {n}/{m}: mean {means.mean():.6f}, error {err:.3e}, error N l {err * n / WINDOWS:.4f}, proved bound {up(bound)}, {time.time() - t0:.2f}s"
        )


def gaps():
    print(
        "other gaps and shared factors, code 7, error against the gap-h limit, global covariance of the two bits"
    )
    k = kernel1(CODE)
    for n, m in GAPS:
        h = m - n
        law, exact = law1(n, m, WINDOWS)
        lim = limit_law1(h, WINDOWS)
        err = np.max(np.abs(law @ k @ law.T - lim @ k @ lim.T))
        tot = exact.sum(axis=0)
        whole = n * m * WINDOWS
        p11, pa, pb = (
            Q(int(tot[3]), whole),
            Q(int(tot[1] + tot[3]), whole),
            Q(int(tot[2] + tot[3]), whole),
        )
        cov = p11 - pa * pb
        print(
            f"  {n}/{m}: gap {h}, gcd {gcd(n, m)}, error {err:.3e}, global covariance {cov}"
        )


def obstruction():
    print("local correlation against global covariance, code 7, 2321/2323")
    n, m = 2321, 2323
    law, exact = law1(n, m, WINDOWS)
    a = law[:, 1] + law[:, 3]
    b = law[:, 2] + law[:, 3]
    pearson = (law[:, 3] - a * b) / np.sqrt(a * (1 - a) * b * (1 - b))
    lim = np.array(
        [float(tri_mean(2, Q(i, WINDOWS), Q(i + 1, WINDOWS))) for i in range(WINDOWS)]
    )
    print(
        f"  local correlation runs {pearson.min():.4f} to {pearson.max():.4f}, max distance to the window mean of the triangle wave {np.max(np.abs(pearson - lim)):.2e}"
    )
    tot = exact.sum(axis=0)
    whole = n * m * WINDOWS
    p11, pa, pb = (
        Q(int(tot[3]), whole),
        Q(int(tot[1] + tot[3]), whole),
        Q(int(tot[2] + tot[3]), whole),
    )
    qn, qm = Q(n - 1, 2 * n), Q(m - 1, 2 * m)
    print(f"  global covariance {p11 - pa * pb}, marginals {pa == qn} {pb == qm}")
    fill_n, fill_m = 1 - qn * qn, 1 - qm * qm
    glob = fill_n + fill_m - 2 * fill_n * fill_m
    k = kernel1(CODE)
    joint = [Q(int(x), whole) for x in tot]
    direct = sum(
        joint[sx] * joint[sy] for sx in range(4) for sy in range(4) if k[sx, sy]
    )
    print(
        f"  global overlay mean from the window integrals {direct} = {float(direct):.9f}, independent value {glob}, equal {direct == glob}, limit 3/8"
    )


# OCTAGON


def octagon():
    kappa = sqrt(2) - 1
    c = (1 - kappa) / 2
    sag = (1 - sqrt(kappa)) ** 2 / sqrt(2)
    chord = sqrt(2) * (1 - kappa)
    print(
        f"regular chord octagon at kappa = sqrt 2 - 1, level c = 1 - 1/sqrt 2 = {c:.6f}"
    )
    print(
        f"  chord {chord:.6f}, sag {sag:.6f}, sag over chord {sag / chord:.6f} in centred units"
    )
    n, m = FINE_PAIR
    law, _ = law1(n, m, FINE)
    k = kernel1(CODE)
    chord_u = (1 - (1 + kappa) / 2) / 2
    arc_u = (1 - sqrt(kappa)) / 2
    out = []
    for u in (chord_u, arc_u):
        i = int(u * FINE)
        val = law[i] @ k @ law[i]
        lim = lag_law(tri_mean(2, Q(i, FINE), Q(i + 1, FINE)))
        lv = sum(
            float(lim[sx]) * float(lim[sy]) * k[sx, sy]
            for sx in range(4)
            for sy in range(4)
        )
        out.append((u, val, lv))
    print(f"  {n}/{m}, window side 1/{FINE} on the diagonal:")
    print(
        f"  chord midpoint u = {out[0][0]:.6f}: box mean {out[0][1]:.6f}, box average of H {out[0][2]:.6f}, H there 0.25"
    )
    print(
        f"  arc point u = {out[1][0]:.6f}: box mean {out[1][1]:.6f}, box average of H {out[1][2]:.6f}, level {c:.6f}"
    )


# LEVEL TWO


def law2(n, m, w):
    big = n * n * m * m * w
    cuts = np.unique(
        np.concatenate(
            [
                np.arange(0, big + 1, m * m * w),
                np.arange(0, big + 1, n * n * w),
                np.arange(0, big + 1, n * n * m * m),
            ]
        )
    )
    left, size = cuts[:-1], np.diff(cuts)
    fn, fm = left // (m * m * w), left // (n * n * w)
    a, c = (fn // n) % 2, (fn % n) % 2
    b, d = (fm // m) % 2, (fm % m) % 2
    state = a + 2 * b + 4 * c + 8 * d
    window = left // (n * n * m * m)
    counts = np.bincount(
        window * 16 + state, weights=size.astype(np.float64), minlength=16 * w
    )
    return counts.reshape(w, 16) / float(n * n * m * m)


def kernel2(code):
    t = corners(code)
    k = np.zeros((16, 16))
    for sx in range(16):
        for sy in range(16):
            ax, bx, cx, dx = sx & 1, (sx >> 1) & 1, (sx >> 2) & 1, (sx >> 3) & 1
            ay, by, cy, dy = sy & 1, (sy >> 1) & 1, (sy >> 2) & 1, (sy >> 3) & 1
            ink_n = t[ax][ay] & t[cx][cy]
            ink_m = t[bx][by] & t[dx][dy]
            k[sx, sy] = float(ink_n != ink_m)
    return k


def par(x):
    return int((x // 1) % 2)


def limit_law2(u, odd):
    s = 2 * u
    cuts = {Q(0), Q(1), Q(2)}
    for k in range(-4, 5):
        cuts.add(k - s)
    for floor_phi in (0, 1):
        for j in (0, 1, 2):
            for k in range(-12, 24):
                cuts.add(floor_phi + Q(k, 4) - s / 2 - Q(odd * j, 4))
    pts = sorted(x for x in cuts if 0 <= x <= 2)
    law = [Q(0)] * 16
    for lo, hi in pairwise(pts):
        if hi == lo:
            continue
        phi = (lo + hi) / 2
        f = phi - (phi // 1)
        j = (f + s) // 1
        lam = 4 * f + 2 * s + odd * j
        rho = tri(lam)
        a, b = par(phi), par(phi + s)
        wgt = (hi - lo) / 2
        for c in (0, 1):
            for d in (0, 1):
                law[a + 2 * b + 4 * c + 8 * d] += (
                    wgt * (1 + (1 if c == d else -1) * rho) / 4
                )
    return law


def simpson_law2(lo, hi, odd):
    pa, pm, pb = (
        limit_law2(lo, odd),
        limit_law2((lo + hi) / 2, odd),
        limit_law2(hi, odd),
    )
    return [(x + 4 * y + z) / 6 for x, y, z in zip(pa, pm, pb)]


def quadratic_check(odd):
    bad = 0
    for q in range(4):
        lo = Q(q, 4)
        xs = [lo + Q(k, 20) for k in range(6)]
        vals = [limit_law2(x, odd) for x in xs]
        for s in range(16):
            y = [v[s] for v in vals]
            d2 = [y[i] - 2 * y[i + 1] + y[i + 2] for i in range(4)]
            if len(set(d2)) != 1:
                bad += 1
    return bad


def h2(u, v, odd, k):
    pu, pv = limit_law2(u, odd), limit_law2(v, odd)
    return sum(pu[sx] * pv[sy] for sx in range(16) for sy in range(16) if k[sx, sy])


def level2():
    k = kernel2(CODE)
    print("level 2, code 7, gap 2")
    for odd in (1, 0):
        print(
            f"  limit law piecewise quadratic on each quarter, {'odd' if odd else 'even'} sides: {quadratic_check(odd)} breaches over 16 states"
        )
    lims = {}
    for odd in (1, 0):
        rows = [
            simpson_law2(Q(i, WINDOWS), Q(i + 1, WINDOWS), odd) for i in range(WINDOWS)
        ]
        lims[odd] = np.array([[float(x) for x in r] for r in rows])
    for n in LEVEL2:
        t0 = time.time()
        odd = n % 2
        law = law2(n, n + 2, WINDOWS)
        means = law @ k @ law.T
        lim = lims[odd] @ k @ lims[odd].T
        other = lims[1 - odd] @ k @ lims[1 - odd].T
        us = (np.arange(WINDOWS) + 0.5) / WINDOWS
        g = np.abs(1 - 2 * us)
        gg = g[:, None] * g[None, :]
        na = 2 * 9 / 16 - 2 * 9 / 16 * (0.5 + gg / 4)
        nb = 2 * 9 / 16 - 2 * (0.5 + gg / 4) ** 2
        print(
            f"  {n}/{n + 2}: error {np.max(np.abs(means - lim)):.3e}, against the other parity {np.max(np.abs(means - other)):.3e}, naive {np.max(np.abs(means - na)):.3e} and {np.max(np.abs(means - nb)):.3e}, {time.time() - t0:.2f}s"
        )
    for odd in (1, 0):
        whole = [
            sum(x)
            for x in zip(*[simpson_law2(Q(q, 4), Q(q + 1, 4), odd) for q in range(4)])
        ]
        whole = [x / 4 for x in whole]
        top = all(
            sum(whole[s] for s in range(16) if (s & 3) == ab) == Q(1, 4)
            for ab in range(4)
        )
        indep = all(whole[s] == Q(1, 16) for s in range(16))
        mean = sum(
            whole[sx] * whole[sy] for sx in range(16) for sy in range(16) if k[sx, sy]
        )
        print(
            f"  {'odd' if odd else 'even'} sides: global law uniform on 16 states {indep}, top pair uniform {top}, mean of the limit {mean} = {float(mean):.6f}, 2 p (1 - p) at p = 9/16 is {2 * Q(9, 16) * Q(7, 16)}"
        )
    closed_form(k)


def tee(y):
    r = y - (y // 1)
    return (1 if (y // 1) % 2 == 0 else -1) * r * (1 - r)


def weight(u):
    return -tee(4 * min(u, 1 - u)) / 2


def alpha(u):
    return (1 + tri(2 * u)) / 4


def gamma(u, odd):
    return (1 + 2 * weight(u)) / 4 if odd else Q(1, 4)


def beta(u):
    return alpha(u) / 4 + weight(u) / 8


def h2_closed(u, v, odd):
    return (
        Q(5, 8)
        - alpha(u) * alpha(v)
        - gamma(u, odd) * gamma(v, odd)
        - 2 * beta(u) * beta(v)
    )


def h2_naive(u, v):
    return Q(9, 16) - Q(9, 8) * alpha(u) * alpha(v)


def closed_form(k):
    pts = [Q(i, 96) for i in range(97)]
    for odd in (1, 0):
        bad = 0
        for u in pts:
            law = limit_law2(u, odd)
            a = sum(law[s] for s in range(16) if s & 3 == 3)
            c = sum(law[s] for s in range(16) if s >> 2 == 3)
            b = law[15]
            bad += int(a != alpha(u)) + int(c != gamma(u, odd)) + int(b != beta(u))
        grid = [Q(i, 12) for i in range(13)]
        bad2 = sum(
            int(h2(u, v, odd, k) != h2_closed(u, v, odd)) for u in grid for v in grid
        )
        print(
            f"  {'odd' if odd else 'even'} sides: alpha, beta, gamma against the lag law at 97 points {bad} mismatches, H2 closed form against the kernel at 169 points {bad2} mismatches"
        )
    pts = [Q(0), Q(1, 8), Q(1, 4)]
    for odd in (1, 0):
        mat = [[h2_closed(u, v, odd) for v in pts] for u in pts]
        det = (
            mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
            - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
            + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0])
        )
        p1, p2 = h2_closed(Q(0), Q(3, 8), odd), h2_closed(Q(1, 4), Q(1, 4), odd)
        print(
            f"  {'odd' if odd else 'even'}: H2 corner {h2_closed(Q(0), Q(0), odd)}, centre {h2_closed(Q(1, 2), Q(1, 2), odd)}, at (3/8, 3/8) {h2_closed(Q(3, 8), Q(3, 8), odd)}"
        )
        print(f"    3 by 3 minor at u, v in 0, 1/8, 1/4: {det}")
        print(
            f"    g(u) g(v) = 1/4 at both (0, 3/8) and (1/4, 1/4): H2 reads {p1} and {p2}"
        )
    fine = [Q(i, 256) for i in range(257)]
    parts = {odd: [(alpha(u), beta(u), gamma(u, odd)) for u in fine] for odd in (1, 0)}

    def at(p, q):
        return Q(5, 8) - p[0] * q[0] - p[2] * q[2] - 2 * p[1] * q[1]

    gap = max(
        abs(at(p, q) - at(r, t))
        for p, r in zip(parts[1], parts[0])
        for q, t in zip(parts[1], parts[0])
    )
    na = {
        odd: max(
            abs(at(p, q) - Q(9, 16) + Q(9, 8) * p[0] * q[0])
            for p in parts[odd]
            for q in parts[odd]
        )
        for odd in (1, 0)
    }
    print(
        f"  max |H2 odd - H2 even| on the 1/256 grid {gap} = {float(gap):.6f}, attained at (3/8, 3/8): {abs(h2_closed(Q(3, 8), Q(3, 8), 1) - h2_closed(Q(3, 8), Q(3, 8), 0)) == gap}"
    )
    print(
        f"  max |H2 - naive| on the 1/256 grid, a lower bound on the sup, truncated: odd {floor(float(na[1]) * 1e6) / 1e6:.6f}, even {floor(float(na[0]) * 1e6) / 1e6:.6f}"
    )
    grid = [Q(i, 64) + Q(1, 128) for i in range(64)]
    for odd in (1, 0):
        mat = np.array([[float(h2_closed(u, v, odd)) for v in grid] for u in grid])
        sv = np.linalg.svd(mat, compute_uv=False)
        rank = int(np.sum(sv > 1e-12 * sv[0]))
        print(
            f"  kernel rank on a 64 by 64 grid, {'odd' if odd else 'even'}: {rank}, singular values {', '.join(f'{x:.3e}' for x in sv[: rank + 1])}"
        )


def main():
    t0 = time.time()
    check_walsh()
    check_direct()
    ladder()
    gaps()
    obstruction()
    octagon()
    level2()
    print(f"total {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
