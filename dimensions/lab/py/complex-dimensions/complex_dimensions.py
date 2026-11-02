import math
import time
from fractions import Fraction
from itertools import pairwise

import numpy as np
from scipy.special import gammaln

LATTICE = [
    (3, (0, 2)),
    (5, (0, 2, 4)),
    (15, (0, 4, 10, 14)),
    (15, (0, 2, 4, 10, 12, 14)),
]
RATIOS = (1 / 3, 1 / 5)
CONTROL_MAPS = ((1 / 3, 0.0), (1 / 5, 4 / 5))
CONTROL_SEED = Fraction(7, 15)
APERIODIC = ((3, (0, 2)), (5, (0, 4)))
POLE_RANGE = 40
REMAX, IMMAX = 3.0, 40.0
BOUNDARY = 200000
UMIN, UMAX, NU = math.log(1 / 0.03), math.log(1e6), 3000
SAFETY = 50.0
OMEGA_GRID = np.linspace(0.5, 30.0, 24001)
NBINS = 40
SAMPLES = 6000
WINDOWS = [6.0, 15.0, 24.0, 33.0, 42.0, 51.0, 60.0]
TAIL = 160
LN3, LN5 = math.log(3), math.log(5)


def label(base, digits):
    return f"base {base}, digits {{{','.join(str(d) for d in digits)}}}"


def gap_seeds(base, digits):
    ds = sorted(digits)
    return [Fraction(b - a - 1, base) for a, b in pairwise(ds) if b - a > 1]


def dimension(base, digits):
    return math.log(len(digits)) / math.log(base)


def solve_moran(ratios):
    d = 0.5
    for _ in range(200):
        f = sum(r**d for r in ratios) - 1
        fp = sum(math.log(r) * r**d for r in ratios)
        d -= f / fp
    return d


def lattice_poles(base, digits):
    k = len(digits)
    d = dimension(base, digits)
    omega = 2 * math.pi / math.log(base)
    s = d + 1j * omega * np.arange(-POLE_RANGE, POLE_RANGE + 1)
    denominator = np.abs(1 - k * np.exp(-s * math.log(base)))
    seeds = np.array([float(g) for g in gap_seeds(base, digits)])
    numerator = np.abs(np.exp(np.outer(s, np.log(seeds))).sum(axis=1))
    d1 = sum(gap_seeds(base, digits))
    print(f"  {label(base, digits)}")
    print(
        f"    d = {d:.6f} ({d:.9f})   omega = 2*pi/ln({base}) = {omega:.6f} ({omega:.9f})"
    )
    print(
        f"    {s.size} predicted poles m = -{POLE_RANGE}..{POLE_RANGE}: max |1 - k n^(-s)| = {denominator.max():.1e}"
    )
    print(
        f"    min |D(s)| = {numerator.min():.6f}   D(1) = {d1} = 1 - k/n: {d1 == 1 - Fraction(k, base)}"
    )


def moran(s):
    return 1 - sum(np.exp(s * math.log(r)) for r in RATIOS)


def moran_prime(s):
    return -sum(math.log(r) * np.exp(s * math.log(r)) for r in RATIOS)


def winding_number():
    corners = [
        -REMAX - 1j * IMMAX,
        REMAX - 1j * IMMAX,
        REMAX + 1j * IMMAX,
        -REMAX + 1j * IMMAX,
    ]
    steps = np.arange(BOUNDARY) / BOUNDARY
    path = np.concatenate(
        [a + (b - a) * steps for a, b in zip(corners, corners[1:] + corners[:1])]
    )
    phase = np.angle(moran(path))
    turns = np.diff(np.concatenate([phase, phase[:1]]))
    turns = (turns + np.pi) % (2 * np.pi) - np.pi
    return round(turns.sum() / (2 * np.pi))


def control_poles():
    re = np.linspace(-REMAX, REMAX, 61)
    im = np.linspace(-IMMAX, IMMAX, 121)
    s = (re[:, None] + 1j * im[None, :]).ravel()
    alive = np.ones(s.size, dtype=bool)
    with np.errstate(all="ignore"):
        for _ in range(200):
            s = s - moran(s) / moran_prime(s)
            alive &= np.isfinite(s) & (np.abs(s.real) <= 40) & (np.abs(s.imag) <= 4000)
        s = np.where(alive, s, np.nan)
        good = (
            alive
            & (np.abs(moran(s)) < 1e-13)
            & (np.abs(s.imag) <= IMMAX)
            & (np.abs(s.real) <= REMAX)
        )
    roots = []
    for z in sorted(s[good], key=lambda z: (z.imag, z.real)):
        if all(abs(z - q) > 1e-7 for q in roots):
            roots.append(z)
    roots = np.array(roots)
    gaps = np.diff(roots.imag)
    d = solve_moran(RATIOS)
    print(
        f"  two-ratio control, ratios 1/3 and 1/5, d = {d:.9f} solving 3^(-d) + 5^(-d) = 1"
    )
    print(
        f"    roots in Re [-{REMAX:.0f},{REMAX:.0f}], Im [-{IMMAX:.0f},{IMMAX:.0f}]: {roots.size}   winding number of 1 - 3^(-s) - 5^(-s) over the box: {winding_number()}"
    )
    print(
        f"    Re range {roots.real.min():.6f} .. {roots.real.max():.6f}   Im gaps {gaps.min():.6f} .. {gaps.max():.6f}"
    )
    for base in (3, 5, 15):
        omega = 2 * math.pi / math.log(base)
        ratio = roots.imag / omega
        print(
            f"    worst offset from spacing 2*pi/ln{base}: {np.abs(ratio - np.round(ratio)).max():.2f} of a step"
        )
    print(f"    ln3/ln5 = {LN3 / LN5:.12f}")


def compose(schedule):
    base, digits = 1, [0]
    for b, ds in schedule:
        digits = [x * b + d for x in digits for d in ds]
        base *= b
    return base, sorted(digits)


def integer_cover(schedule, levels):
    lefts, den = [0], 1
    for i in range(levels):
        b, ds = schedule[i % len(schedule)]
        lefts = [x * b + d for x in lefts for d in ds]
        den *= b
    return sorted(Fraction(x, den) for x in lefts)


def composition():
    for schedule in [((3, (0, 2)), (5, (0, 4))), ((3, (0, 2)), (5, (0, 2, 4)))]:
        base, digits = compose(schedule)
        print(
            f"  {' then '.join(label(b, ds) for b, ds in schedule)} -> {label(base, digits)}"
        )
    alternating = integer_cover(((3, (0, 2)), (5, (0, 4))), 8)
    product = integer_cover(((15, (0, 4, 10, 14)),), 4)
    print(
        f"    8 alternating levels and 4 base-15 levels: {len(alternating)} and {len(product)} intervals, identical as exact fractions: {alternating == product}"
    )


def level_cover(base, digits, floor):
    lefts = np.zeros(1)
    length = 1.0
    while length > floor:
        length /= base
        lefts = (lefts[:, None] + np.array(digits, dtype=float) * length).ravel()
    return np.sort(lefts), np.full(lefts.size, length)


def schedule_cover(schedule, floor):
    lefts = np.zeros(1)
    length = 1.0
    used = 0
    lnk, lnn = 0.0, 0.0
    while length > floor:
        base, digits = schedule[used % len(schedule)]
        used += 1
        length /= base
        lnk += math.log(len(digits))
        lnn += math.log(base)
        lefts = (lefts[:, None] + np.array(digits, dtype=float) * length).ravel()
    return np.sort(lefts), np.full(lefts.size, length), used, lnk / lnn


def ratio_cover(maps, floor):
    lefts, lens = np.zeros(1), np.ones(1)
    done_l, done_s = [], []
    while lefts.size:
        small = lens <= floor
        done_l.append(lefts[small])
        done_s.append(lens[small])
        lefts, lens = lefts[~small], lens[~small]
        if not lefts.size:
            break
        lefts = np.concatenate([lefts + c * lens for _, c in maps])
        lens = np.concatenate([r * lens for r, _ in maps])
    lefts, lens = np.concatenate(done_l), np.concatenate(done_s)
    order = np.argsort(lefts)
    return lefts[order], lens[order]


def thue_morse(n):
    return [i.bit_count() % 2 for i in range(n)]


def box_count(lefts, rights, eps):
    lo = np.floor(lefts / eps).astype(np.int64)
    hi = np.floor(rights / eps).astype(np.int64)
    return int((hi - lo + 1).sum() - np.count_nonzero(hi[:-1] == lo[1:]))


def detrended(lefts, rights, d):
    u = np.linspace(UMIN, UMAX, NU)
    g = np.array([math.log(box_count(lefts, rights, math.exp(-x))) for x in u]) - d * u
    slope, intercept = np.polyfit(u, g, 1)
    return u, g - (slope * u + intercept)


def power_at(w, dt, y):
    return abs(np.exp(-1j * w * dt) @ y) ** 2


def golden(f, lo, hi):
    phi = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(80):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = f(d)
    return (a + b) / 2


def periodogram_peak(u, g):
    j = np.arange(g.size)
    window = (
        0.42
        - 0.5 * np.cos(2 * np.pi * j / (g.size - 1))
        + 0.08 * np.cos(4 * np.pi * j / (g.size - 1))
    )
    y = g * window
    y = y - y.mean()
    dt = u - u[0]
    power = np.concatenate(
        [
            np.abs(np.exp(-1j * np.outer(chunk, dt)) @ y) ** 2
            for chunk in np.array_split(OMEGA_GRID, 48)
        ]
    )
    i = int(np.argmax(power))
    lo, hi = OMEGA_GRID[max(i - 1, 0)], OMEGA_GRID[min(i + 1, OMEGA_GRID.size - 1)]
    return golden(lambda w: -power_at(w, dt, y), lo, hi)


def folded_variance(u, g, period, nbins):
    phase = ((u - u[0]) / period) % 1.0
    index = np.minimum((phase * nbins).astype(int), nbins - 1)
    sums = np.bincount(index, g, nbins)
    counts = np.bincount(index, minlength=nbins)
    profile = np.where(counts > 0, sums / np.maximum(counts, 1), 0.0)
    return 1 - (g - profile[index]).var() / g.var()


def box_report(name, lefts, lens, d, base):
    u, g = detrended(lefts, lefts + lens, d)
    omega = periodogram_peak(u, g)
    fold = [folded_variance(u, g, math.log(p), NBINS) for p in (3, 5, 15)]
    print(f"  {name}: {lefts.size} intervals, d = {d:.6f}")
    if base:
        predicted = 2 * math.pi / math.log(base)
        err = abs(omega - predicted) / predicted * 100
        print(
            f"    predicted omega {predicted:.4f}   periodogram peak {omega:.4f}   error {err:.3f}%   within 1%: {err < 1}"
        )
    else:
        print(f"    no lattice prediction   periodogram peak {omega:.4f}")
    print(
        f"    folding variance explained at ln3 {fold[0]:.3f}   ln5 {fold[1]:.3f}   ln15 {fold[2]:.3f}"
    )


def box_counts():
    floor = math.exp(-UMAX) / SAFETY
    print(
        f"  u in [{UMIN:.4f},{UMAX:.4f}], {NU} points, covers refined below eps/{SAFETY:.0f}, {NBINS} folding bins"
    )
    for base, digits in LATTICE:
        lefts, lens = level_cover(base, digits, floor)
        box_report(label(base, digits), lefts, lens, dimension(base, digits), base)
    lefts, lens = ratio_cover(CONTROL_MAPS, floor)
    box_report("two-ratio control x/3, (x+4)/5", lefts, lens, solve_moran(RATIOS), None)
    schedule = [APERIODIC[bit] for bit in thue_morse(64)]
    lefts, lens, used, d = schedule_cover(schedule, floor)
    box_report(
        f"aperiodic control, Thue-Morse schedule of {label(3, (0, 2))} and {label(5, (0, 4))}, {used} levels",
        lefts,
        lens,
        d,
        None,
    )


def tube_digit(base, k, seeds, t):
    total = 0.0
    for g in seeds:
        a0 = 0 if g <= t else math.ceil(math.log(g / t) / math.log(base))
        total += t * (k**a0 - 1) / (k - 1) + g * (k / base) ** a0 * base / (base - k)
    return total


def brute_digit(base, digits, level, t):
    lefts = np.zeros(1, dtype=np.int64)
    for _ in range(level):
        lefts = (lefts[:, None] * base + np.array(sorted(digits))).ravel()
    den = base**level
    gaps = [Fraction(int(b - a - 1), den) for a, b in pairwise(lefts) if b - a > 1]
    cut = Fraction(t)
    covered = sum(min(g, cut) for g in gaps)
    return float(covered), lefts.size / den


A, B = np.meshgrid(np.arange(TAIL), np.arange(TAIL), indexing="ij")
LOG_COUNT = gammaln(A + B + 1) - gammaln(A + 1) - gammaln(B + 1)
LOG_WEIGHT = LOG_COUNT - A * LN3 - B * LN5
LOG_SCALE = A * LN3 + B * LN5


def tube_control(t):
    g = float(CONTROL_SEED)
    small = LOG_SCALE >= math.log(g / t)
    return t * np.exp(LOG_COUNT[~small]).sum() + g * np.exp(LOG_WEIGHT[small]).sum()


def content_curve(measure, d, ulo, uhi):
    u = np.linspace(ulo, uhi, SAMPLES)
    return u, np.array([measure(2 * math.exp(-x)) * math.exp(-x) ** (d - 1) for x in u])


def window_swings(measure, d):
    u, m = content_curve(measure, d, WINDOWS[0], WINDOWS[-1])
    swings = []
    for lo, hi in pairwise(WINDOWS):
        vals = m[(u >= lo) & (u <= hi)]
        swings.append((vals.max() - vals.min()) / vals.mean() * 100)
    return swings


def periodicity_defect(measure, d, period):
    worst = 0.0
    for u in np.linspace(40.0, 50.0, 600):
        a = measure(2 * math.exp(-u)) * math.exp(-u) ** (d - 1)
        b = measure(2 * math.exp(-u - period)) * math.exp(-u - period) ** (d - 1)
        worst = max(worst, abs(a - b) / abs(a))
    return worst


def swing_line(swings):
    return "    swing of M(eps) per window u in " + ", ".join(
        f"[{lo:.0f},{hi:.0f}] {s:.4f}%"
        for (lo, hi), s in zip(pairwise(WINDOWS), swings)
    )


def tube_lattice(base, digits):
    k = len(digits)
    d = dimension(base, digits)
    seeds = [float(g) for g in gap_seeds(base, digits)]
    measure = lambda t: tube_digit(base, k, seeds, t)
    print(
        f"  {label(base, digits)}: total gap length {sum(seeds) * base / (base - k):.12f}"
    )
    brute, tail = brute_digit(base, digits, 6, 0.01)
    closed = measure(0.01)
    print(
        f"    level 6, t = 0.01: brute {brute:.9f}   closed form {closed:.9f}   excess {closed - brute:.6e} <= cover length {tail:.6e}: {-1e-12 <= closed - brute <= tail * (1 + 1e-9) + 1e-12}"
    )
    print(swing_line(window_swings(measure, d)))
    defects = [periodicity_defect(measure, d, math.log(p)) for p in (3, 5, 15)]
    print(
        f"    max relative defect of M(eps) = M(eps/p) on u in [40,50]: p=3 {defects[0]:.1e}   p=5 {defects[1]:.1e}   p=15 {defects[2]:.1e}"
    )


def tube_two_ratio():
    d = solve_moran(RATIOS)
    total = CONTROL_SEED / (1 - sum(Fraction(1, round(1 / r)) for r in RATIOS))
    print(
        f"  two-ratio control, gap seed {CONTROL_SEED}, d = {d:.9f}, total gap length {total}, eps at u = 60: {math.exp(-60):.1e}"
    )
    swings = window_swings(tube_control, d)
    print(swing_line(swings))
    print(
        f"    swing decays {swings[0]:.2f}% -> {swings[-1]:.2f}%, monotone: {all(a > b for a, b in pairwise(swings))}"
    )


def cantor_profile():
    d = dimension(3, (0, 2))
    profile = lambda t: 2 ** (1 - d) * (t ** (d - 1) + t**d)
    tstar = (1 - d) / d
    lo, hi = profile(tstar), profile(1.0)
    seeds = [float(g) for g in gap_seeds(3, (0, 2))]
    _, m = content_curve(lambda t: tube_digit(3, 2, seeds, t), d, 50.0, 60.0)
    print("  Cantor limit profile 2^(1-d) (t^(d-1) + t^d), t in [1/3, 1)")
    print(
        f"    minimum {lo:.9f} at t* = (1-d)/d = {tstar:.6f}   maximum {hi:.9f} at the ends   swing {(hi - lo) / lo * 100:.2f}%"
    )
    print(
        f"    measured on u in [50,60]: min {m.min():.9f}   max {m.max():.9f}   min gap {abs(m.min() - lo):.1e}   max gap {abs(m.max() - hi):.1e}"
    )


def main():
    start = time.time()
    print(
        "POLES: zeros of 1 - k n^(-s) at s = d + 2 pi i m/ln(n), the numerator D(s) there, and the two-ratio Moran roots"
    )
    for base, digits in LATTICE:
        lattice_poles(base, digits)
    control_poles()
    print()
    print(
        "COMPOSITION: alternating one base per level multiplies into the product base"
    )
    composition()
    print()
    print(
        "BOX COUNT: g(u) = ln N(e^-u) - d u, Blackman periodogram on a direct DFT grid, and period folding"
    )
    box_counts()
    print()
    print(
        "TUBE: V(eps) = sum over gaps of min(gap, 2 eps) in closed form, M(eps) = eps^(d-1) V(eps)"
    )
    for base, digits in LATTICE:
        tube_lattice(base, digits)
    tube_two_ratio()
    cantor_profile()
    print()
    print(f"wall {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
