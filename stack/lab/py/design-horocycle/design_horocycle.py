import argparse
import math
import time

import numpy as np

N = 2**24
JMIN = 8.0
JMAX = 23.0
PER_OCTAVE = 16
SEED = 1729
DESIGNS = (
    ("base 3 {0,1}", 3, (0, 1)),
    ("base 5 {0,1,2}", 5, (0, 1, 2)),
    ("base 10 missing 7", 10, (0, 1, 2, 3, 4, 5, 6, 8, 9)),
)
T0 = time.time()


def clock(label):
    print(f"[{time.time() - T0:6.1f}s] {label}")


# ARITHMETIC


def mobius(n):
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    for p in np.flatnonzero(sieve).tolist():
        mu[p::p] *= -1
        if p <= n // p:
            mu[p * p :: p * p] = 0
    return mu


def totients(n):
    phi = np.arange(n + 1, dtype=np.int64)
    rem = np.arange(n + 1, dtype=np.int64)
    r = math.isqrt(n)
    small = np.ones(r + 1, dtype=bool)
    small[:2] = False
    for p in range(2, math.isqrt(r) + 1):
        if small[p]:
            small[p * p :: p] = False
    for p in np.flatnonzero(small).tolist():
        phi[p::p] -= phi[p::p] // p
        pk = p
        while pk <= n:
            rem[pk::pk] //= p
            pk *= p
    big = rem > 1
    phi[big] -= phi[big] // rem[big]
    return phi


def design(base, digits, n):
    dg = np.array(digits, dtype=np.int64)
    cur = dg[dg > 0]
    out = []
    while cur.size:
        cur = cur[cur <= n]
        out.append(cur)
        cur = cur[cur * base <= n]
        cur = (cur[:, None] * base + dg[None, :]).ravel()
    return np.sort(np.concatenate(out))


def convolve_id(n, ds, w):
    nu = np.zeros(n + 1, dtype=np.int64)
    d0 = math.isqrt(n)
    k0 = int(np.searchsorted(ds, d0, side="right"))
    for d, wd in zip(ds[:k0].tolist(), w[:k0].tolist()):
        nu[d::d] += wd * np.arange(1, n // d + 1, dtype=np.int64)
    big, wb = ds[k0:], w[k0:]
    for e in range(1, n // (d0 + 1) + 1):
        k = int(np.searchsorted(big, n // e, side="right"))
        if k == 0:
            break
        nu[big[:k] * e] += wb[:k] * e
    return nu


def brute_nu(c, members, mu):
    return sum(
        int(mu[e]) * (c // e) for e in range(1, c + 1) if c % e == 0 and e in members
    )


# WINDOWS

GL_X, GL_W = np.polynomial.legendre.leggauss(2000)
GL_U = 1.5 + 0.5 * GL_X
GL_WU = 0.5 * GL_W


def bump_cinf(u):
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u)
    ins = (u > 1.0) & (u < 2.0)
    v = u[ins]
    out[ins] = np.exp(4.0 - 1.0 / ((v - 1.0) * (2.0 - v)))
    return out


def bump_c2(u):
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u)
    ins = (u > 1.0) & (u < 2.0)
    v = u[ins]
    out[ins] = 64.0 * ((v - 1.0) * (2.0 - v)) ** 3
    return out


def moment(fn, s):
    return float(np.sum(GL_WU * fn(GL_U) * GL_U ** (s - 1.0)))


def balanced(fn):
    shift = moment(fn, 3.0) / moment(fn, 2.0)
    return lambda u: fn(u) * (np.asarray(u, dtype=float) - shift), shift


# BRIDGE

BX, BW = np.polynomial.legendre.leggauss(80)


def integrate(fn, lo, hi):
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    return half * float(np.dot(BW, fn(mid + half * BX)))


def g_transform(fn, t):
    if t >= 1:
        return 0.0
    hi = math.sqrt(1 / (t * t) - 1)
    lo = math.sqrt(1 / (2 * t * t) - 1) if t * t <= 0.5 else 0.0
    return 2 * integrate(lambda v: fn(1 / (t * t * (1 + v * v))), lo, hi)


def sieve_weight(n, ds, w):
    a = np.zeros(n + 1, dtype=np.int64)
    for d, wd in zip(ds.tolist(), w.tolist()):
        if d > n:
            break
        a[d::d] += wd
    return a


def horocycle_lattice(fn, h, a):
    q = math.isqrt(int(1 / h))
    total = 0.0
    pairs = 0
    arcs = 0
    for c in range(1, q + 1):
        s = h - c * c * h * h
        if s < 0:
            continue
        width = math.sqrt(s) / c
        inner = math.sqrt(max(h / 2 - c * c * h * h, 0.0)) / c
        for d in range(math.ceil(-c * (1 + width)), math.floor(c * width) + 1):
            wgt = int(a[math.gcd(c, abs(d))])
            if wgt == 0:
                continue
            centre = -d / c
            hit = 0
            for lo, hi in (
                (centre - width, centre - inner),
                (centre + inner, centre + width),
            ):
                lo, hi = max(0.0, lo), min(1.0, hi)
                if lo < hi:
                    total += wgt * integrate(
                        lambda t: fn(h / ((c * t + d) ** 2 + c * c * h * h)), lo, hi
                    )
                    hit += 1
            pairs += hit > 0
            arcs += hit
    return total, pairs, arcs


def theta_point(fn, z, a):
    x, y = z.real, z.imag
    total = 0.0
    d = 1
    while y / (d * d) >= 1:
        total += int(a[d]) * float(fn(y / (d * d)))
        d += 1
    for c in range(1, math.isqrt(int(1 / y)) + 2):
        base = math.floor(-c * x)
        for d in range(base - 2, base + 4):
            total += int(a[math.gcd(c, abs(d))]) * float(fn(y / abs(c * z + d) ** 2))
    return total


def bridge(mu, phi):
    print("BRIDGE")
    f = bump_cinf
    area_f = moment(f, -1.0)
    c_half = moment(f, -0.5)
    print(
        f"window exp(4 - 1/((u-1)(2-u))) on [1, 2]: F(-1) = int f w^-2 = {area_f:.15f}, F(-1/2) = int f w^-3/2 = {c_half:.15f}"
    )
    rng = np.random.default_rng(SEED)
    points = [complex(rng.uniform(0, 1), rng.uniform(0.003, 0.03)) for _ in range(5)]
    for name, base, digits in (DESIGNS[0], DESIGNS[2]):
        ds = design(base, digits, N)
        w = mu[ds].astype(np.int64)
        keep = w != 0
        dsq, wq = ds[keep], w[keep]
        a = sieve_weight(4096, dsq, wq)
        nu = convolve_id(4096, dsq, wq)
        members = set(ds[ds <= 300].tolist())
        bad = sum(int(nu[c]) != brute_nu(c, members, mu) for c in range(1, 301))
        print(
            f"{name}: nu_F by the two-regime sieve against the divisor definition at every c <= 300: {bad} mismatches"
        )
        print(
            f"{name}: k  h=4^-k  weighted pairs  arcs  lattice integral  h sum nu_F(c) g(c h^1/2)  gap"
        )
        worst = 0.0
        for k in range(2, 7):
            h = 4.0 ** (-k)
            lhs, pairs, arcs = horocycle_lattice(f, h, a)
            rhs = h * sum(
                int(nu[c]) * g_transform(f, c * math.sqrt(h))
                for c in range(1, 2**k + 1)
            )
            worst = max(worst, abs(lhs - rhs))
            print(
                f"{name}: {k}  {h:.3e}  {pairs}  {arcs}  {lhs:.15f}  {rhs:.15f}  {abs(lhs - rhs):.1e}"
            )
        print(f"{name}: worst gap {worst:.1e}")
        drift = 0.0
        for z in points:
            ref = theta_point(f, z, a)
            for img in (-1 / z, z + 1, z / (z + 1), (2 * z + 1) / (5 * z + 3)):
                drift = max(drift, abs(theta_point(f, img, a) - ref))
        print(
            f"{name}: Theta_F f at 5 seeded points and 4 images each under SL(2, Z): max |Theta(gamma z) - Theta(z)| = {drift:.1e}"
        )
        z1 = math.fsum((wq / dsq.astype(float)).tolist())
        z2 = math.fsum((wq / dsq.astype(float) ** 2).tolist())
        print(
            f"{name}: Z_F(1) partial sum to 2^24 = {z1:.10f}, Z_F(2) partial sum to 2^24 = {z2:.12f}"
        )
        point = -0.5 * z1 * c_half
        if base == 3:
            levels = math.floor(math.log(N, 3)) + 1
            top = design(3, digits, 3**levels - 1)
            top = top[top > N]
            s16 = math.fsum((1.0 / top.astype(float)).tolist())
            geo = sum((2.0 / 3.0) ** (L - 1) for L in range(levels + 1, 400))
            radius = sum(2.0**L / 3.0 ** (L - 1) for L in range(levels + 1, 400))
            assert s16 + geo <= radius
            lo, hi = z1 - radius, z1 + radius
            print(
                f"{name}: the elements above 2^24 are the {top.size} with {levels} digits above it, reciprocals summing to {s16:.6f}, and those of L > {levels} digits, 2^(L-1) of them each at least 3^(L-1), at most sum (2/3)^(L-1) = {geo:.6f}; tail <= {s16 + geo:.6f} <= {radius:.6f} = sum_(L > {levels}) 2^L/3^(L-1); Z_F(1) lies in [{lo:.6f}, {hi:.6f}]"
            )
            print(
                f"{name}: k  (C_F(h) - (pi/2) Z_F(2) F(-1))/h^(1/2)  against -(1/2) Z_F(1) F(-1/2) in [{-0.5 * c_half * hi:.6f}, {-0.5 * c_half * lo:.6f}], {point:.6f} at the partial sum"
            )
        else:
            print(
                f"{name}: k  (C_F(h) - (pi/2) Z_F(2) F(-1))/h^(1/2)  against -(1/2) Z_F(1) F(-1/2) = {point:.6f} at the uncertified partial sum"
            )
        nu_big = convolve_id(2**12, dsq, wq)
        for k in range(4, 13, 2):
            h = 4.0 ** (-k)
            y = math.sqrt(h)
            ch = h * sum(
                int(nu_big[c]) * g_transform(f, c * y) for c in range(1, 2**k + 1)
            )
            print(f"{name}: {k}  {(ch - math.pi / 2 * z2 * area_f) / y:.6f}")
    print(
        "full set control: k  (C(h) - (pi/2) F(-1)/zeta(2))/h^(1/2), the h^(1/2) term absent since 1/zeta(1) = 0"
    )
    for k in range(4, 13, 2):
        h = 4.0 ** (-k)
        y = math.sqrt(h)
        ch = h * sum(int(phi[c]) * g_transform(f, c * y) for c in range(1, 2**k + 1))
        print(f"full set: {k}  {(ch - math.pi / 2 * 6 / math.pi**2 * area_f) / y:.6f}")
    clock("bridge done")


# METER


def grid():
    j = np.arange(JMIN, JMAX + 1e-9, 1.0 / PER_OCTAVE)
    return j, 2.0 ** (-j)


def meter(arrays, y, fn):
    out = np.zeros((len(arrays), len(y)))
    for i, yy in enumerate(y):
        lo = math.ceil(1.0 / yy)
        hi = min(math.floor(2.0 / yy), N)
        win = fn(np.arange(lo, hi + 1, dtype=np.float64) * yy)
        for r, arr in enumerate(arrays):
            out[r, i] = (
                yy * yy * float(np.dot(arr[lo : hi + 1].astype(np.float64), win))
            )
    return out


def fit(x, v):
    a, b = np.polyfit(x, v, 1)
    return float(a), float(np.sqrt(np.mean((v - (a * x + b)) ** 2)))


def octave_rows(j, values, mass):
    oct_ = np.floor(j + 1e-9).astype(int)
    rows = []
    for k in sorted(set(oct_.tolist())):
        sel = oct_ == k
        if sel.sum() >= PER_OCTAVE // 2:
            rows.append(
                (
                    k,
                    float(np.sqrt(np.mean(values[sel] ** 2))),
                    float(np.exp(np.mean(np.log(mass[sel])))),
                )
            )
    return np.array(rows)


def slope_line(label, j, err, y, mass):
    rows = octave_rows(j, np.abs(err) / y**2, mass)
    x = np.log(rows[:, 2])
    v = np.log(rows[:, 1])
    half = len(rows) // 2
    a, r = fit(x, v)
    a_lo, _ = fit(x[:half], v[:half])
    a_hi, _ = fit(x[half:], v[half:])
    top = [
        float((v[i + 1] - v[i]) / (x[i + 1] - x[i]))
        for i in range(len(rows) - 4, len(rows) - 1)
    ]
    print(
        f"{label}: fitted slope {a:.4f} (residual {r:.3f}), lower {a_lo:.4f}, upper {a_hi:.4f}, top three octave-to-octave {', '.join(f'{t:.3f}' for t in top)}"
    )
    return a, rows


def phase_bins(label, j, err, y, mass, base, theta):
    resid = np.log(np.abs(err) / y**2 + 1e-300) - theta * np.log(mass)
    phase = (j * math.log(2) / math.log(base)) % 1.0
    bins = np.floor(phase * 8).astype(int)
    means = [float(np.mean(resid[bins == b])) for b in range(8)]
    spread = max(means) - min(means)
    print(
        f"{label}: base-phase bin means of log residual {', '.join(f'{m:+.2f}' for m in means)}, spread {spread:.2f}"
    )
    return spread


def run_meter(mu, phi):
    print("METER")
    j, y = grid()
    windows = []
    for wname, raw in (("cinf", bump_cinf), ("c2", bump_c2)):
        fn, shift = balanced(raw)
        f2 = float(np.sum(GL_WU * fn(GL_U) * GL_U))
        print(
            f"window {wname} balanced: f(u) = b(u)(u - {shift:.12f}), F(2) = {f2:.1e}"
        )
        windows.append((wname, fn))
    print(
        f"grid y = 2^-j, j in [{JMIN:.0f}, {JMAX:.0f}] step 1/{PER_OCTAVE}, {len(j)} samples, windows (1/y, 2/y] inside (2^8, 2^24], sieve to N = 2^24 = {N}"
    )
    rng = np.random.default_rng(SEED)
    summary = []
    for name, base, digits in DESIGNS:
        ds = design(base, digits, N)
        alpha = math.log(len(digits)) / math.log(base)
        w = mu[ds].astype(np.int64)
        keep = w != 0
        dsq, wq = ds[keep], w[keep]
        signs = rng.choice(np.array([-1, 1], dtype=np.int64), size=dsq.size)
        mass = np.searchsorted(ds, 1.0 / y, side="right").astype(float)
        top = int(
            np.searchsorted(ds, N, side="right")
            - np.searchsorted(ds, N // 2, side="right")
        )
        print(
            f"{name}: alpha = {alpha:.6f}, A_F(2^24) = {len(ds)}, squarefree {dsq.size}, elements in the top window (2^23, 2^24] {top}, A_F(1/y) from {int(mass[0])} to {int(mass[-1])}"
        )
        arrays = [
            convolve_id(N, dsq, wq),
            convolve_id(N, dsq, signs),
            convolve_id(N, dsq, np.abs(wq)),
        ]
        clock(f"{name}: three convolutions done")
        for wname, fn in windows:
            if wname == "cinf":
                direct = meter(
                    [arrays[0]], np.array([2.0**-12, 2.0**-14, 2.0**-16]), fn
                )[0]
                other = []
                for jj in (12, 14, 16):
                    yy = 2.0 ** (-jj)
                    sel = dsq <= 2.0 / yy
                    tot = 0.0
                    for d, wd in zip(dsq[sel].tolist(), wq[sel].tolist()):
                        e = np.arange(
                            1, math.floor(2.0 / (d * yy)) + 1, dtype=np.float64
                        )
                        tot += wd * float(np.dot(e, fn(e * d * yy)))
                    other.append(yy * yy * tot)
                print(
                    f"{name}: E_F at j = 12, 14, 16 by the convolution {', '.join(f'{v:+.6e}' for v in direct)}; by sum over d in S_F of mu(d) R(dy) {', '.join(f'{v:+.6e}' for v in other)}; max gap {max(abs(p - q) for p, q in zip(direct, other)):.1e}"
                )
            err = meter(arrays, y, fn)
            th_mu, rows = slope_line(f"{name} {wname} mu", j, err[0], y, mass)
            th_rand, _ = slope_line(f"{name} {wname} random signs", j, err[1], y, mass)
            th_abs, _ = slope_line(f"{name} {wname} |mu|", j, err[2], y, mass)
            if wname == "cinf":
                print(f"{name}: octave  A_F  rms |E_F|/y^2  ratio to A_F^(1/2)")
                for k, rms, m in rows:
                    print(
                        f"{name}: {int(k)}  {m:.1f}  {rms:.4e}  {rms / math.sqrt(m):.4f}"
                    )
                s_mu = phase_bins(f"{name} mu", j, err[0], y, mass, base, th_mu)
                s_rand = phase_bins(
                    f"{name} random signs", j, err[1], y, mass, base, th_rand
                )
                s_abs = phase_bins(f"{name} |mu|", j, err[2], y, mass, base, th_abs)
                ratio = rows[:, 1] / np.sqrt(rows[:, 2])
                point = np.abs(err[0]) / y**2 / np.sqrt(mass)
                flips = int(np.sum(np.sign(err[0][1:]) != np.sign(err[0][:-1])))
                summary.append(
                    (
                        name,
                        alpha,
                        th_mu,
                        th_rand,
                        th_abs,
                        s_mu,
                        s_rand,
                        s_abs,
                        float(ratio.min()),
                        float(ratio.max()),
                        float(rows[0, 2]),
                        float(rows[-1, 2]),
                        float(point.min()),
                        float(point.max()),
                        flips,
                    )
                )
        del arrays
        clock(f"{name} done")
    massf = 1.0 / y
    for wname, fn in windows:
        err = meter([phi], y, fn)
        slope_line(f"full set {wname} mu (nu = phi)", j, err[0], y, massf)
    clock("full set done")
    print("SUMMARY")
    for name, alpha, a, b, c, s1, s2, s3, r0, r1, m0, m1, p0, p1, flips in summary:
        print(
            f"{name}: alpha {alpha:.4f}, fitted slope mu {a:.4f}, random signs {b:.4f}, |mu| {c:.4f}; per-octave rms of |E_F|/y^2 over A_F^(1/2) in [{r0:.4f}, {r1:.4f}] while A_F^(1/2) grows by {math.sqrt(m1 / m0):.1f}; pointwise |E_F|/(y^2 A_F^(1/2)) in [{p0:.4f}, {p1:.4f}] with {flips} sign changes; base-phase spread mu {s1:.2f}, random signs {s2:.2f}, |mu| {s3:.2f}"
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "verb", choices=["bridge", "meter", "all"], nargs="?", default="all"
    )
    args = ap.parse_args()
    clock(f"sieves to {N}")
    mu = mobius(N)
    phi = totients(N)
    clock("sieves done")
    if args.verb in ("bridge", "all"):
        bridge(mu, phi)
    if args.verb in ("meter", "all"):
        run_meter(mu, phi)
    clock("done")


if __name__ == "__main__":
    main()
