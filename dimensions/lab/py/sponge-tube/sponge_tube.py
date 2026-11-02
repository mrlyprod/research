import math
import time
from fractions import Fraction as Fr

import numpy as np
from mpmath import iv, mp

mp.prec = 133

iv.prec = 133
D = iv.log(20) / iv.log(3)
WEIGHT = iv.mpf(27) / 20
CELLW = (3, 2, 3)
LEVELS = 400
DEEP_LEVELS = 7
TAIL_LEVELS = 80


def num(x):
    return iv.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else iv.mpf(x)


def asin(u):
    return iv.atan2(u, iv.sqrt(1 - u * u))


def a0(t, d):
    t = min(t, d)
    x, dd = num(t), num(d)
    if t == d:
        return iv.pi * dd * dd / 4
    return x / 2 * iv.sqrt(dd * dd - x * x) + dd * dd / 2 * asin(x / dd)


def a1(t, d):
    t = min(t, d)
    x, dd = num(t), num(d)
    if t == d:
        return dd**3 / 3
    y = dd * dd - x * x
    return (
        x
        * x
        * (3 * dd**4 - 3 * dd * dd * x * x + x**4)
        / (3 * (dd**3 + y * iv.sqrt(y)))
    )


def seg(t0, t1, alpha, beta, d):
    if t1 <= t0:
        return iv.mpf(0)
    return num(alpha) * (a0(t1, d) - a0(t0, d)) - num(beta) * (a1(t1, d) - a1(t0, d))


def hole_j(s, d):
    return 4 * seg(Fr(0), s / 2, s, 2, d)


def partial_j(s, c, d):
    left = seg(Fr(0), min(s / 2, c), s, 2, d)
    right = seg(max(Fr(0), s - c), s / 2, s, 2, d)
    if c <= s / 2:
        bottom = seg(Fr(0), c, c, 1, d)
    else:
        bottom = seg(Fr(0), s - c, c, 1, d) + seg(s - c, s / 2, s, 2, d)
    return left + right + 2 * bottom


def column_weight(i, digits):
    w = 1
    for _ in range(digits):
        w *= CELLW[i % 3]
        i //= 3
    return w


def columns_below(count, digits):
    if count >= 3**digits:
        return 8**digits
    total, prefix = 0, 1
    for d in range(digits - 1, -1, -1):
        digit = (count // 3**d) % 3
        total += prefix * sum(CELLW[:digit]) * 8**d
        prefix *= CELLW[digit]
    return total


def wall_half(delta):
    total = iv.mpf(0)
    for m in range(1, LEVELS + 1):
        total += 8 ** (m - 1) * hole_j(Fr(1, 3 ** (m + 1)), delta)
    total /= 2
    assert upper(total) <= lower(num(delta / 18))
    return total, num(delta * Fr(8, 9) ** LEVELS / 16)


def strip(delta):
    total = iv.mpf(0)
    for m in range(1, LEVELS + 1):
        s = Fr(1, 3 ** (m + 1))
        ratio = delta / s
        full = int((ratio - 2) // 3) + 1
        if full > 0:
            total += columns_below(min(full, 3 ** (m - 1)), m - 1) * hole_j(s, delta)
        cut = int((ratio - 1) // 3)
        if 0 <= cut < 3 ** (m - 1) and (3 * cut + 1) * s < delta < (3 * cut + 2) * s:
            total += column_weight(cut, m - 1) * partial_j(
                s, delta - (3 * cut + 1) * s, delta
            )
    assert upper(total) <= lower(num(delta * delta / 3))
    return total, num(delta * Fr(8, 9) ** LEVELS / 8)


def hole_rows(i, digits):
    rows = [0]
    for d in range(digits):
        digit = (i // 3**d) % 3
        allowed = (0, 2) if digit == 1 else (0, 1, 2)
        rows = [r + a * 3**d for r in rows for a in allowed]
    return sorted(rows)


def overlap(ints_a, ints_b):
    total, j = Fr(0), 0
    for lo, hi in ints_a:
        while j < len(ints_b) and ints_b[j][1] <= lo:
            j += 1
        k = j
        while k < len(ints_b) and ints_b[k][0] < hi:
            total += min(hi, ints_b[k][1]) - max(lo, ints_b[k][0])
            k += 1
    return total


def columns_meeting(lo, hi, s, digits):
    first = max(0, int((lo / s - 2) // 3))
    last = min(3**digits - 1, int((hi / s - 1) // 3) + 1)
    return [
        i
        for i in range(first, last + 1)
        if (3 * i + 2) * s > lo and (3 * i + 1) * s < hi
    ]


def deep_bound(delta):
    total = Fr(0)
    cache = {}
    for m in range(1, DEEP_LEVELS + 1):
        sm = Fr(1, 3 ** (m + 1))
        wm = sm * sm / (4 * delta)
        for n in range(1, DEEP_LEVELS + 1):
            sn = Fr(1, 3 ** (n + 1))
            wn = sn * sn / (4 * delta)
            for i in columns_meeting(delta - wn, delta, sm, m - 1):
                vlen = min(delta, (3 * i + 2) * sm) - max(delta - wn, (3 * i + 1) * sm)
                if (m, i) not in cache:
                    cache[(m, i)] = [
                        ((3 * r + 1) * sm, (3 * r + 2) * sm)
                        for r in hole_rows(i, m - 1)
                    ]
                for k in columns_meeting(delta - wm, delta, sn, n - 1):
                    ulen = min(delta, (3 * k + 2) * sn) - max(
                        delta - wm, (3 * k + 1) * sn
                    )
                    if (n, k) not in cache:
                        cache[(n, k)] = [
                            ((3 * r + 1) * sn, (3 * r + 2) * sn)
                            for r in hole_rows(k, n - 1)
                        ]
                    total += vlen * ulen * overlap(cache[(m, i)], cache[(n, k)])
    for m in range(1, TAIL_LEVELS + 1):
        for n in range(1, TAIL_LEVELS + 1):
            if max(m, n) <= DEEP_LEVELS:
                continue
            wm = min(Fr(1, 3 ** (2 * m + 2)) / (4 * delta), delta)
            wn = min(Fr(1, 3 ** (2 * n + 2)) / (4 * delta), delta)
            count = min(math.ceil(wn * 3**m) + 2, math.ceil(wm * 3**n) + 2)
            total += wm * wn * min(Fr(count, 9), Fr(1, 3))
    total += 2 * Fr(1, 3 ** (2 * TAIL_LEVELS + 6)) / (36 * delta * delta) / Fr(8, 9) / 8
    return total


def tube(delta):
    d = num(delta)
    half, half_tail = wall_half(delta)
    strip_v, strip_tail = strip(delta)
    value = iv.pi * d * d - 8 * iv.sqrt(2) * d**3 + 8 * d * d + 48 * (half - strip_v)
    return value - 48 * strip_tail - 24 * num(deep_bound(delta)), value + 48 * half_tail


def tube_upper(delta):
    s_sum = Fr(0)
    for m in range(1, LEVELS + 1):
        s = Fr(1, 3 ** (m + 1))
        s_sum += 8 ** (m - 1) * s * delta * min(s, 4 * delta)
    s_sum += delta * Fr(8, 9) ** LEVELS / 8
    return iv.pi * num(delta) ** 2 + 24 * num(s_sum)


def periodic(eps, depth):
    lo = hi = iv.mpf(20) / 27
    for l in range(depth + 1):
        tl, th = tube(eps / 3**l)
        lo += WEIGHT**l * tl
        hi += WEIGHT**l * th
    hi += WEIGHT ** (depth + 1) * tube_upper(eps / 3 ** (depth + 1)) * iv.mpf(20) / 9
    scale = iv.exp((D - 3) * iv.log(num(eps)))
    return lower(lo * scale), upper(hi * scale)


def carpet_distance(v, z, side, depth=40):
    v, z = v / side, z / side
    s = np.full(v.shape, side)
    out = np.zeros(v.shape)
    done = np.zeros(v.shape, dtype=bool)
    for _ in range(depth):
        s = s / 3
        dv, dz = np.floor(v * 3).astype(int), np.floor(z * 3).astype(int)
        hole = (dv == 1) & (dz == 1) & ~done
        rel = np.minimum.reduce([v * 3 - 1, 2 - v * 3, z * 3 - 1, 2 - z * 3])
        out[hole] = s[hole] * rel[hole]
        done |= hole
        v, z = v * 3 - dv, z * 3 - dz
    return out


def raster(delta, n):
    h = 1 / (3 * n)
    grid = (np.arange(n) + 0.5) * h
    u, v, z = np.meshgrid(grid, grid, grid, indexing="ij")
    dist = np.full(u.shape, np.inf)
    dvz, duz = carpet_distance(v, z, 1 / 3), carpet_distance(u, z, 1 / 3)
    for t, dd in ((u, dvz), (1 / 3 - u, dvz), (v, duz), (1 / 3 - v, duz)):
        dist = np.minimum(dist, np.sqrt(t * t + dd * dd))
    arm = dist
    e = np.minimum.reduce(
        [
            np.sqrt(np.minimum(u, 1 / 3 - u) ** 2 + np.minimum(v, 1 / 3 - v) ** 2),
            np.sqrt(np.minimum(u, 1 / 3 - u) ** 2 + np.minimum(z, 1 / 3 - z) ** 2),
            np.sqrt(np.minimum(v, 1 / 3 - v) ** 2 + np.minimum(z, 1 / 3 - z) ** 2),
        ]
    )
    r = np.sqrt(3) * h / 2
    inner = 6 * (arm <= delta - r).sum() + (e <= delta - r).sum()
    outer = 6 * (arm <= delta + r).sum() + (e <= delta + r).sum()
    return inner * h**3, outer * h**3


def lower(x):
    return mp.make_mpf(x._mpi_[0]) if hasattr(x, "_mpi_") else x


def upper(x):
    return mp.make_mpf(x._mpi_[1]) if hasattr(x, "_mpi_") else x


def down(x, k):
    return mp.mpf(int(mp.floor(lower(x) * 10**k))) / 10**k


def up(x, k):
    return mp.mpf(int(mp.ceil(upper(x) * 10**k))) / 10**k


def main():
    start = time.time()
    print(
        f"SPONGE: base 3, fill 20, D = log 20 / log 3 in [{mp.nstr(lower(D), 12)}, {mp.nstr(upper(D), 12)}], weight per level 27/20, constant term 20/27"
    )
    print(
        "TUBE: T(delta) = volume of F_delta inside the plus, delta in (0, 1/6], from the wall carpets and the centre cube edges, interval arithmetic at 133 bits"
    )
    for delta in (Fr(1, 8), Fr(1, 12)):
        lo, hi = tube(delta)
        lo, hi = lower(lo), upper(hi)
        rl, rh = raster(float(delta), 120)
        print(
            f"  delta = {delta}: closed form in [{mp.nstr(down(lo, 9), 10)}, {mp.nstr(up(hi, 9), 10)}], raster band [{rl:.5f}, {rh:.5f}] at 120^3 cells per cube"
        )
        assert rl <= float(lo) and float(hi) <= rh
    lo, hi = tube(Fr(1, 6))
    ident = (iv.pi + 8) / 36 - iv.sqrt(2) / 27
    assert lower(lo) <= upper(ident) and lower(ident) <= upper(hi)
    print(
        f"  delta = 1/6: A1 = V1 there, so T = (pi + 8)/36 - sqrt2/27 - 24 Deep = {mp.nstr(down(ident, 9), 10)} less at most {float(24 * deep_bound(Fr(1, 6))):.2e}; closed form in [{mp.nstr(down(lo, 9), 10)}, {mp.nstr(up(hi, 9), 10)}]"
    )
    print(
        "PERIODIC: p(eps) = eps^(D-3) (20/27 + sum_l (27/20)^l T(eps/3^l)), eps in (sqrt2/18, 1/6]"
    )
    bands = {}
    for eps in (Fr(1, 12), Fr(1, 8), Fr(1, 6)):
        lo, hi = periodic(eps, 40)
        bands[eps] = (lo, hi)
        print(
            f"  eps = {eps}: p in [{mp.nstr(down(lo, 6), 7)}, {mp.nstr(up(hi, 6), 7)}]"
        )
    gap = bands[Fr(1, 6)][0] - bands[Fr(1, 12)][1]
    gap2 = bands[Fr(1, 8)][0] - bands[Fr(1, 12)][1]
    assert gap > 0 and gap2 > 0
    print(
        f"  p(1/6) - p(1/12) >= {mp.nstr(down(gap, 6), 7)},  p(1/8) - p(1/12) >= {mp.nstr(down(gap2, 6), 7)}"
    )
    swing = gap / bands[Fr(1, 12)][1]
    print(
        f"  relative swing (p(1/6) - p(1/12)) / p(1/12) >= {mp.nstr(down(swing * 100, 4), 5)} %"
    )
    print(
        f"  doubly deep set, upper bound on its volume: {float(deep_bound(Fr(1, 6))):.3g} at delta = 1/6, {float(deep_bound(Fr(1, 8))):.3g} at delta = 1/8, {float(deep_bound(Fr(1, 12))):.3g} at delta = 1/12"
    )
    print(f"wall {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
