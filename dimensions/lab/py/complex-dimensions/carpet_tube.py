import math
import time
from fractions import Fraction

import numpy as np
from mpmath import mp, mpf
from scipy.ndimage import distance_transform_edt

mp.dps = 50
BASE = 3
FILL = 8
DIM = mp.log(FILL) / mp.log(BASE)
DIM_F = float(DIM)
LOW = (Fraction(1), Fraction(4, 5), Fraction(-4, 7))
HIGH = (Fraction(9, 8), Fraction(3, 10), Fraction(-1, 14))
SEAM = Fraction(379, 280)


def hole_tube(side, eps):
    return side * side if side <= 2 * eps else 4 * eps * side - 4 * eps * eps


def split_index(eps):
    m = 1
    while Fraction(1, BASE**m) > 2 * eps:
        m += 1
    return m


def direct_sum(eps):
    m0 = split_index(eps)
    ring = sum(
        FILL ** (m - 1) * hole_tube(Fraction(1, BASE**m), eps) for m in range(1, m0)
    )
    return ring + Fraction(FILL, BASE**2) ** (m0 - 1)


def closed_v(eps):
    m0 = split_index(eps)
    a = Fraction(4, 5) * eps * (Fraction(FILL, BASE) ** (m0 - 1) - 1)
    b = Fraction(4, 7) * eps * eps * (FILL ** (m0 - 1) - 1)
    return a - b + Fraction(FILL, BASE**2) ** (m0 - 1)


def branch(t):
    return LOW if t < Fraction(1, 2) else HIGH


def to_mpf(x):
    return mpf(x.numerator) / x.denominator if isinstance(x, Fraction) else mpf(x)


def profile(t):
    x = to_mpf(t)
    c0, c1, c2 = branch(t if isinstance(t, Fraction) else Fraction(float(x)))
    return x ** (DIM - 2) * (to_mpf(c0) + to_mpf(c1) * x + to_mpf(c2) * x * x)


def profile_float(t):
    c0, c1, c2 = branch(Fraction(t))
    return t ** (DIM_F - 2) * (float(c0) + float(c1) * t + float(c2) * t * t)


def stationary(coeffs):
    c0, c1, c2 = (to_mpf(c) for c in coeffs)
    a = -c2 * DIM
    b = -c1 * (DIM - 1)
    c = -c0 * (DIM - 2)
    disc = mp.sqrt(b * b - 4 * a * c)
    return sorted([(-b - disc) / (2 * a), (-b + disc) / (2 * a)])


def filled_cells(level):
    cells = np.ones((1, 1), dtype=bool)
    for _ in range(level):
        block = np.ones((BASE, BASE), dtype=bool)
        block[1, 1] = False
        cells = np.kron(cells, block)
    return cells


def hole_census(level):
    cells = filled_cells(level)
    side = BASE**level
    padded = np.pad(cells, 1, constant_values=True)
    counted = 0
    for m in range(1, level + 1):
        hole = BASE ** (level - m)
        parents = filled_cells(m - 1)
        origins = np.argwhere(parents) * (BASE * hole) + hole
        assert len(origins) == FILL ** (m - 1)
        for i, j in origins:
            assert not cells[i : i + hole, j : j + hole].any()
            ring = padded[i : i + hole + 2, j : j + hole + 2].copy()
            ring[1:-1, 1:-1] = True
            assert ring.all()
        counted += len(origins) * hole * hole
    assert counted == side * side - FILL**level
    return level, FILL**level, side * side - FILL**level


def raster_tube(level, cells_in):
    cells = filled_cells(level)
    side = BASE**level
    dist = distance_transform_edt(~cells)
    count = cells.sum() + ((dist <= cells_in + 1e-9) & ~cells).sum()
    return Fraction(int(count), side * side)


def window(ulo, uhi, samples):
    worst = 0.0
    lo, hi = math.inf, -math.inf
    for u in np.linspace(ulo, uhi, samples):
        eps = math.exp(-u)
        n = math.floor(-math.log(eps) / math.log(BASE))
        t = eps * BASE**n
        if t >= 1:
            n -= 1
            t = eps * BASE**n
        if t < 1 / BASE:
            n += 1
            t = eps * BASE**n
        m0 = split_index(Fraction(eps))
        ring = sum(
            FILL ** (m - 1) * float(hole_tube(Fraction(1, BASE**m), Fraction(eps)))
            for m in range(1, m0)
        )
        v = ring + (FILL / BASE**2) ** (m0 - 1)
        m = eps ** (DIM_F - 2) * v
        worst = max(worst, abs(m - profile_float(t)))
        lo, hi = min(lo, m), max(hi, m)
    return lo, hi, worst


def main():
    start = time.time()
    print(
        f"CARPET: base {BASE}, fill {FILL}, d = log {FILL} / log {BASE} = {mp.nstr(DIM, 12)}"
    )
    print(
        "HOLES: every hole is an open square whose boundary lies in the carpet, and the holes exhaust the complement"
    )
    for level in (3, 4, 5):
        lv, filled, empty = hole_census(level)
        print(
            f"  level {lv}: {filled} filled cells, {empty} empty cells in 8^(m-1) holes of side 3^-m, m = 1..{lv}, every hole ringed by filled cells"
        )
    print()
    print(
        "TUBE: V(eps) = sum_m 8^(m-1) h(3^-m, eps), h(s, eps) = s^2 if s <= 2 eps else 4 eps s - 4 eps^2, against the closed form, exact rationals"
    )
    grid = [
        Fraction(1, 3),
        Fraction(2, 5),
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(2, 3),
        Fraction(4, 5),
        Fraction(9, 10),
        Fraction(999, 1000),
    ]
    checks = 0
    for n in range(1, 13):
        for t in grid:
            eps = t / BASE**n
            assert direct_sum(eps) == closed_v(eps)
            checks += 1
    print(
        f"  {checks} pairs (n, t), n = 1..12, t in {{1/3, 2/5, 1/2, 3/5, 2/3, 4/5, 9/10, 999/1000}}: direct sum equals closed form exactly"
    )
    print()
    print("PROFILE: eps = 3^-n t, t in [1/3, 1), M(eps) = eps^(d-2) V(eps) -> G(t)")
    print(
        "  G(t) = t^(d-2) (1 + 4t/5 - 4t^2/7) on [1/3, 1/2),  G(t) = t^(d-2) (9/8 + 3t/10 - t^2/14) on [1/2, 1)"
    )
    half = Fraction(1, 2)
    left = sum(c * half**i for i, c in enumerate(LOW))
    right = sum(c * half**i for i, c in enumerate(HIGH))
    assert left == right == Fraction(44, 35)
    seam_low = (
        Fraction(1, FILL) * sum(c * Fraction(1, 3) ** i for i, c in enumerate(LOW)) * 9
    )
    seam_high = sum(HIGH)
    assert seam_low == seam_high == SEAM
    slope_low = sum(i * c * half ** (i - 1) for i, c in enumerate(LOW) if i)
    slope_high = sum(i * c * half ** (i - 1) for i, c in enumerate(HIGH) if i)
    assert slope_low == slope_high == Fraction(8, 35)
    print(
        f"  continuous at t = 1/2: both branches give 44/35 * 2^(2-d) = {mp.nstr(profile(half), 12)}, and C^1 there: both polynomial factors have slope 8/35"
    )
    print(f"  periodic: G(1/3) = G(1) = 379/280 = {mp.nstr(mpf(379) / 280, 12)}")
    tmax = [r for r in stationary(LOW) if Fraction(1, 3) <= Fraction(float(r)) < half]
    tmin = [r for r in stationary(HIGH) if half <= Fraction(float(r)) < 1]
    assert len(tmax) == 1 and len(tmin) == 1
    tmax, tmin = tmax[0], tmin[0]
    gmax, gmin = profile(tmax), profile(tmin)
    scan = [profile_float(t) for t in np.linspace(1 / 3, 1, 200001)]
    assert min(scan) >= float(gmin) - 1e-12 and max(scan) <= float(gmax) + 1e-12
    swing = (gmax - gmin) / gmin
    scale = mpf(10) ** 11
    print(
        f"  maximum {mp.nstr(gmax, 15)}, safe {mp.nstr(mp.ceil(gmax * scale) / scale, 12)} up, at t = {mp.nstr(tmax, 9)} (root of (4/7) d t^2 - (4/5)(d-1) t - (d-2) = 0 in [1/3, 1/2))"
    )
    print(
        f"  minimum {mp.nstr(gmin, 15)}, safe {mp.nstr(mp.floor(gmin * scale) / scale, 12)} down, at t = {mp.nstr(tmin, 9)} (root of (d/14) t^2 - (3/10)(d-1) t - (9/8)(d-2) = 0 in [1/2, 1))"
    )
    print(
        f"  swing (max - min) / min = {mp.nstr(swing * 100, 8)} %   scan of 200001 points stays inside [min, max]"
    )
    print(
        f"  Cantor design for scale: swing 3.53 %; the carpet swings {mp.nstr(swing * 100, 4)} %"
    )
    print()
    print("MEASURED: the level sum at eps = e^-u against G(t)")
    lo, hi, worst = window(50.0, 60.0, 2001)
    print(
        f"  u in [50, 60]: min {lo:.12f}   max {hi:.12f}   max |M(eps) - G(t)| = {worst:.1e}"
    )
    assert worst < 1e-9
    eps = math.exp(-60)
    print(
        f"  outer collar at u = 60: eps^(D-2) (4 eps + pi eps^2) = {eps ** (DIM_F - 2) * (4 * eps + math.pi * eps * eps):.1e}"
    )
    print()
    print(
        "RASTER: level-6 cells, Euclidean distance transform to the nearest filled cell centre, no hole lemma used"
    )
    level, k = 6, 21
    side = BASE**level
    raster = raster_tube(level, k)
    exact = closed_v(Fraction(k, side))
    print(
        f"  cells within {k} cell widths of a filled centre, plus the filled cells: {raster} = {float(raster):.12f}"
    )
    print(f"  exact V({k}/{side}) from the closed form: {exact} = {float(exact):.12f}")
    assert raster == exact
    print("  equal as rationals")
    print()
    print(f"wall {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
