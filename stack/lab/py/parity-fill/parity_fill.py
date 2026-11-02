import time
from fractions import Fraction
from math import gcd, sqrt

import numpy as np

N_SPIN = 55
ODDS_SPIN = list(range(1, N_SPIN + 1, 2))
L_SPIN = len(ODDS_SPIN)
GOLDEN = 137.507764
EXACT_TOP = 21
SPIN_R = 1024
SWEEP_R = 512
CHECK_R = 4096
LAYER_COUNTS = [4, 8, 14, 28]


def lcm(a, b):
    return a // gcd(a, b) * b


def odd_scales(top):
    return list(range(3, top + 1, 2))


def mass_table(scales):
    period = 1
    for n in scales:
        period = lcm(period, n)
    steps = [period // n for n in scales]
    counts = np.zeros(1 << len(scales), dtype=np.int64)
    block = 1 << 20
    for lo in range(0, period, block):
        j = np.arange(lo, min(lo + block, period), dtype=np.int64)
        code = np.zeros(j.shape, dtype=np.int64)
        for i, m in enumerate(steps):
            code |= ((j // m) & 1) << i
        counts += np.bincount(code, minlength=1 << len(scales))
    for i in range(len(scales)):
        bit = 1 << i
        for mask in range(1 << len(scales)):
            if not mask & bit:
                counts[mask] += counts[mask | bit]
    return period, counts


def fill_exact(period, counts, scales, top):
    bits = [i for i, n in enumerate(scales) if n <= top]
    total = 0
    for sub in range(1 << len(bits)):
        mask = 0
        size = 0
        for i, b in enumerate(bits):
            if sub >> i & 1:
                mask |= 1 << b
                size += 1
        total += (-2) ** size * int(counts[mask]) ** 2
    return Fraction(period * period - total, 2 * period * period)


def indicator(n, period):
    j = np.arange(period, dtype=np.int64)
    return ((j // (period // n)) & 1).astype(bool)


def fill_literal(top):
    scales = odd_scales(top)
    period = 1
    for n in scales:
        period = lcm(period, n)
    acc = np.zeros((period, period), dtype=bool)
    for n in scales:
        c = indicator(n, period)
        acc ^= c[:, None] & c[None, :]
    return Fraction(int(acc.sum()), period * period)


def fill_independent(top):
    prod = Fraction(1)
    for n in odd_scales(top):
        p = Fraction((n - 1) ** 2, 4 * n * n)
        prod *= 1 - 2 * p
    return (1 - prod) / 2


def fill_raster(top, r):
    u = (np.arange(r, dtype=np.float64) + 0.5) / r
    acc = np.zeros((r, r), dtype=bool)
    for n in odd_scales(top):
        c = parity(0.5 * n * u)
        acc ^= c[:, None] & c[None, :]
    return float(acc.mean())


def primes_up_to_count(k):
    out = []
    c = 2
    while len(out) < k:
        if all(c % d for d in range(2, int(sqrt(c)) + 1)):
            out.append(c)
        c += 1
    return out


def schedules():
    ps = primes_up_to_count(L_SPIN - 1)
    sched = {
        "unspun": [0.0] * L_SPIN,
        "degrees": [float(k + 1) for k in range(L_SPIN)],
        "primes": [0.0] + [float(p) for p in ps],
        "golden": [GOLDEN * (k + 1) for k in range(L_SPIN)],
    }
    for q in (2, 3, 4, 5):
        sched["eye 90/%d" % q] = [90.0 / q * k for k in range(L_SPIN)]
    return sched


def parity(z):
    return z - np.floor(z) >= 0.5


def spun_field(angles, r, layers=L_SPIN, dtype=np.float32):
    step = (np.arange(r, dtype=dtype) + dtype(0.5)) * dtype(1.0 / r)
    d = step - dtype(0.5)
    acc = np.zeros((r, r), dtype=bool)
    for k in range(layers):
        n = ODDS_SPIN[k]
        h = 0.5 * n
        t = np.float64(angles[k]) * np.pi / 180.0
        cx = dtype(h * np.cos(t))
        sx = dtype(h * np.sin(t))
        c0 = dtype(h * 0.5)
        a = c0 + cx * d[None, :] + sx * d[:, None]
        b = c0 - sx * d[None, :] + cx * d[:, None]
        acc ^= parity(a) & parity(b)
    return acc


def disc_mask(r):
    d = (np.arange(r, dtype=np.float64) + 0.5) / r - 0.5
    return d[:, None] ** 2 + d[None, :] ** 2 <= 0.25


def spun_fill(angles, r, layers=L_SPIN, disc=True):
    f = spun_field(angles, r, layers)
    return float(f[disc_mask(r)].mean()) if disc else float(f.mean())


def raster_row(top, r):
    u = (np.arange(r, dtype=np.float64) + 0.5) / r
    acc = np.zeros((r, r), dtype=bool)
    out = {}
    for n in odd_scales(top):
        c = parity(0.5 * n * u)
        acc ^= c[:, None] & c[None, :]
        out[n] = float(acc.mean())
    return out


def section_exact():
    scales = odd_scales(EXACT_TOP)
    period, counts = mass_table(scales)
    tops = list(range(3, EXACT_TOP + 1, 2))
    ras = raster_row(EXACT_TOP, CHECK_R)
    print(
        "exact fill of the parity fold, odd scales 3..N, 1D grid period %d, %d subsets"
        % (period, 1 << len(scales))
    )
    exact = {}
    for top in tops:
        f = fill_exact(period, counts, scales, top)
        exact[top] = f
        lit = fill_literal(top) if top <= 9 else None
        print(
            " N %2d L %2d fill %s = %.9f raster %.9f gap %.2e"
            % (top, top // 2 + 1, f, float(f), ras[top], abs(float(f) - ras[top])),
            ""
            if lit is None
            else ("literal 2D XOR %s %s" % (lit, "match" if lit == f else "MISMATCH")),
        )
    zero = sum(1 for k in range(1 << len(scales)) if counts[k] == 0)
    print(
        " joint 1D masses: %d of the %d subsets have m_S = 0, the smallest being {3,5,7}"
        % (zero, 1 << len(scales))
    )
    print(
        " joint 1D masses m_S of the first scales",
        " ".join(
            "m{%s} = %s"
            % (
                ",".join(str(scales[i]) for i in range(len(scales)) if k >> i & 1),
                Fraction(int(counts[k]), period),
            )
            for k in (1, 2, 3, 5, 7)
        ),
    )
    return exact


def section_independent(exact):
    print("exact against independent Bernoulli layers, p_n = ((n-1)/(2n))^2")
    for top in sorted(exact):
        f = exact[top]
        g = fill_independent(top)
        de = Fraction(1, 2) - f
        di = Fraction(1, 2) - g
        print(
            " N %2d exact %.9f independent %.9f difference %+.3e deviation %.6e against %.6e ratio %.6f"
            % (
                top,
                float(f),
                float(g),
                float(f - g),
                float(de),
                float(di),
                float(de / di),
            )
        )
    print(
        " decay of the exact deviation, successive ratios",
        " ".join(
            "%.6f"
            % float((Fraction(1, 2) - exact[t + 2]) / (Fraction(1, 2) - exact[t]))
            for t in sorted(exact)[:-1]
        ),
    )
    print(
        " decay of the independent deviation, successive ratios",
        " ".join(
            "%.6f"
            % float(
                (Fraction(1, 2) - fill_independent(t + 2))
                / (Fraction(1, 2) - fill_independent(t))
            )
            for t in sorted(exact)[:-1]
        ),
    )


def section_spun(exact):
    sched = schedules()
    print(
        "spun parity fill at N %d, %d layers, R %d, inscribed disc"
        % (N_SPIN, L_SPIN, SPIN_R)
    )
    rows = []
    for name, angles in sched.items():
        v = spun_fill(angles, SPIN_R)
        rows.append((name, v))
        print(" %-9s fill %.6f distance to 1/2 %.6f" % (name, v, abs(v - 0.5)))
    for name in ("unspun", "degrees", "eye 90/3"):
        print(
            " %-9s against resolution" % name,
            " ".join(
                "R %d %.6f" % (r, spun_fill(sched[name], r))
                for r in (256, 512, 1024, 2048)
            ),
        )
    best = min(rows, key=lambda r: abs(r[1] - 0.5))
    worst = max(rows, key=lambda r: abs(r[1] - 0.5))
    print(
        " closest to 1/2 %s at %.6f, furthest %s at %.6f"
        % (best[0], best[1], worst[0], worst[1])
    )
    print(
        "fill against layer count at R %d, exact and square raster on the unit square, spun readings on the disc"
        % SPIN_R
    )
    for lc in LAYER_COUNTS:
        top = 2 * lc - 1
        ex = "%.6f" % float(exact[top]) if top in exact else "-"
        print(
            " L %2d N %2d exact %10s square %.6f disc"
            % (lc, top, ex, spun_fill(sched["unspun"], SPIN_R, lc, disc=False)),
            " ".join(
                "%s %.6f" % (name, spun_fill(angles, SPIN_R, lc))
                for name, angles in sched.items()
                if name in ("unspun", "degrees", "primes", "golden")
            ),
        )


def section_sweep():
    print(
        "fixed increment sweep at N %d, R %d, disc fill against increment in whole degrees"
        % (N_SPIN, SWEEP_R)
    )
    rows = []
    for deg in range(91):
        angles = [float(deg) * k for k in range(L_SPIN)]
        rows.append((deg, spun_fill(angles, SWEEP_R)))
    for deg, v in rows:
        print(" increment %2d fill %.6f" % (deg, v))
    order = sorted(rows, key=lambda r: r[1])
    print(" lowest", " ".join("%d:%.6f" % r for r in order[:5]))
    print(" highest", " ".join("%d:%.6f" % r for r in order[-5:]))
    near = sorted(rows, key=lambda r: abs(r[1] - 0.5))
    print(" closest to 1/2", " ".join("%d:%.6f" % r for r in near[:5]))
    print(
        " eyes q = 1, 5, 3, 2 at increments 0, 18, 30, 45 and the quarter turn",
        " ".join("%d:%.6f" % (d, dict(rows)[d]) for d in (0, 18, 30, 45, 90)),
    )
    v = dict(rows)
    print(
        " mirror check fill(d) = fill(90 - d): %d of 46 pairs equal to the bit, largest gap %.3e, quarter-turn check fill(90) = fill(0): %s"
        % (
            sum(1 for d in range(46) if v[d] == v[90 - d]),
            max(abs(v[d] - v[90 - d]) for d in range(46)),
            v[0] == v[90],
        )
    )
    print(
        " spread of the 89 nonzero whole increments: min %.6f max %.6f, unspun %.6f"
        % (min(v[d] for d in range(1, 90)), max(v[d] for d in range(1, 90)), v[0])
    )


def main():
    t0 = time.perf_counter()
    exact = section_exact()
    section_independent(exact)
    section_spun(exact)
    section_sweep()
    print("seconds", round(time.perf_counter() - t0, 2))


main()
