import math
import subprocess
import sys
import time
from bisect import bisect_right
from fractions import Fraction
from itertools import pairwise
from math import comb

import numpy as np

K = 5
HALF = Fraction(1, 2)
THIRD = Fraction(1, 3)
CHILDREN = (
    (HALF, (Fraction(0), Fraction(0))),
    (THIRD, (Fraction(2, 3), Fraction(0))),
    (THIRD, (Fraction(2, 3), Fraction(1, 3))),
    (THIRD, (Fraction(2, 3), Fraction(2, 3))),
    (THIRD, (Fraction(0), Fraction(2, 3))),
    (THIRD, (Fraction(1, 3), Fraction(2, 3))),
)
LN2, LN3 = math.log(2), math.log(3)
NUMER, DENOM = 65, 41
HEIGHT = 60.0
TERMS = 19
UMAX = 300.0
WINDOW = 10.0
WINDOW_STARTS = (10.0, 20.0, 40.0, 80.0, 160.0, 290.0)
SPEC_LO = 50.0
SPEC_STEP = 0.002
NBINS = 40
OMEGA_MAX = 600.0
PEAKS = 10
TEXTBOOK = 0.7675115443 + 45.55415979j
TOLERANCE = 1e-4
RIPPLE_BINS = 24
RIPPLE_PERIODS = 4
RIPPLE_CENTRES = (
    ((0.0, 0.0), 0.5, 2 / 3),
    ((1.0, 0.0), 1 / 3, 1 / 3),
    ((1.0, 1.0), 1 / 3, 1 / 3),
)
STEP = 2e-4
T0 = time.time()


def clock(label):
    print(f"[{time.time() - T0:6.1f}s] {label}")


def gp(script):
    return subprocess.run(
        ["gp", "-q"],
        input="default(realprecision, 60);\n" + script,
        capture_output=True,
        text=True,
        check=True,
    ).stdout


def number(text):
    return float(text.strip().replace(" E", "E"))


def dimension(k):
    return number(gp(f"print(solve(s = 0, 4, 2^-s + {k}*3^-s - 1))"))


def left_edge(k):
    return number(gp(f"print(solve(s = -2, 4, {k}*3^-s - 1 - 2^-s))"))


def overlap(a, b):
    (ra, (xa, ya)), (rb, (xb, yb)) = a, b
    return xa < xb + rb and xb < xa + ra and ya < yb + rb and yb < ya + ra


def letter_grid():
    grid = [["." for _ in range(6)] for _ in range(6)]
    for index, (r, (x, y)) in enumerate(CHILDREN):
        n = int(r * 6)
        for i in range(n):
            for j in range(n):
                grid[int(y * 6) + j][int(x * 6) + i] = "H" if index == 0 else str(index)
    return ["".join(row) for row in reversed(grid)]


def words(level):
    cells = [(Fraction(1), Fraction(0), Fraction(0), 0, 0)]
    for _ in range(level):
        cells = [
            (size * r, x + size * tx, y + size * ty, a + (r == HALF), b + (r == THIRD))
            for size, x, y, a, b in cells
            for r, (tx, ty) in CHILDREN
        ]
    return cells


def render(level):
    side = 6**level
    image = np.zeros((side, side), dtype=np.int32)
    for size, x, y, _, _ in words(level):
        n, i, j = int(size * side), int(x * side), int(y * side)
        image[j : j + n, i : i + n] += 1
    return image


def census(depth):
    counts = {}
    frontier = [(Fraction(1), Fraction(0), Fraction(0), 0, 0)]
    for _ in range(depth + 1):
        nxt = []
        for size, x, y, a, b in frontier:
            counts[(a, b)] = counts.get((a, b), 0) + 1
            if a + b < depth:
                nxt.extend(
                    (
                        size * r,
                        x + size * tx,
                        y + size * ty,
                        a + (r == HALF),
                        b + (r == THIRD),
                    )
                    for r, (tx, ty) in CHILDREN
                )
        frontier = nxt
    return counts


def patch():
    print("the letter on the 6-grid, H the 1/2 child, 1..5 the 1/3 children, . the gap")
    for row in letter_grid():
        print("  " + row)
    inside = all(
        0 <= x and x + r <= 1 and 0 <= y and y + r <= 1 for r, (x, y) in CHILDREN
    )
    disjoint = all(
        not overlap(CHILDREN[i], CHILDREN[j]) for i in range(6) for j in range(i + 1, 6)
    )
    print(
        f"children inside the closed unit square: {inside}   open images pairwise disjoint (open set condition): {disjoint}"
    )
    area = sum(r * r for r, _ in CHILDREN)
    print(f"area of the children {area} = 1/4 + {K}/9; gap {1 - area}")
    out = gp(
        f"D = solve(s = 0, 4, 2^-s + {K}*3^-s - 1); print(D); print(log(29)/log(6)); print(log(2)/log(3))"
    ).split()
    print(f"dimension D solving 2^-s + {K}*3^-s = 1: {out[0][:22]}")
    print(
        f"level render reading log(29)/log(6) = {out[1][:12]}   log 2/log 3 = {out[2][:12]}"
    )
    for level in range(1, 5):
        image = render(level)
        fill = int((image > 0).sum())
        print(
            f"  level {level}: side 6^{level} = {6**level}, fill {fill} = (9 + 4*{K})^{level}: {fill == (9 + 4 * K) ** level}   max cover {image.max()}"
        )
    tile = (render(1) > 0).astype(np.int32)
    kron = int((np.kron(tile, tile) != (render(2) > 0)).sum())
    print(
        f"level 2 render against the Kronecker square of the level 1 tile: {kron} of {36 * 36} cells differ"
    )
    depth = 7
    counts = census(depth)
    good = all(counts[(a, b)] == comb(a + b, a) * K**b for (a, b) in counts)
    print(
        f"cell census to word length {depth}: {len(counts)} sizes 2^-a 3^-b, every count C(a+b, a) {K}^b: {good}"
    )
    clock("patch")


def moran(s, k):
    return 1 - np.exp(-s * LN2) - k * np.exp(-s * LN3)


def moran_prime(s, k):
    return LN2 * np.exp(-s * LN2) + k * LN3 * np.exp(-s * LN3)


def seeds(k, height):
    out = gp(
        f'r = polroots(1 - x^{DENOM} - {k}*x^{NUMER}); for(j = 1, #r, z = -{DENOM}*log(r[j])/log(2); print(real(z), " ", imag(z)))'
    )
    rows = np.array(
        [[float(t) for t in line.split()] for line in out.strip().splitlines()]
    )
    z = rows[:, 0] + 1j * rows[:, 1]
    period = 2 * math.pi * DENOM / LN2
    z = np.concatenate([z + 1j * period * m for m in (-1, 0, 1)])
    return z[np.abs(z.imag) <= height + 2]


def newton(z, k):
    with np.errstate(all="ignore"):
        for _ in range(60):
            z = z - moran(z, k) / moran_prime(z, k)
    return z


def roots_in_box(k, lo, hi, height):
    re = np.linspace(lo, hi, 9)
    im = np.linspace(-height, height, int(8 * height) + 1)
    grid = (re[:, None] + 1j * im[None, :]).ravel()
    lattice = seeds(k, height)
    z = newton(np.concatenate([lattice, grid]), k)
    with np.errstate(all="ignore"):
        good = (
            np.isfinite(z)
            & (np.abs(moran(z, k)) < 1e-12)
            & (z.real >= lo)
            & (z.real <= hi)
            & (np.abs(z.imag) <= height)
        )
    found = []
    for w in sorted(z[good], key=lambda w: (w.imag, w.real)):
        if all(abs(w - q) > 1e-8 for q in found):
            found.append(w)
    polished = newton(lattice, k)
    with np.errstate(all="ignore"):
        hit = (
            np.isfinite(polished)
            & (np.abs(moran(polished, k)) < 1e-12)
            & (np.abs(polished.imag) <= height)
            & (polished.real >= lo)
            & (polished.real <= hi)
        )
    from_lattice = sum(
        1 for w in found if np.min(np.abs(polished[hit] - w), initial=1.0) < 1e-8
    )
    return np.array(found), from_lattice


def edge_floor(k, lo, hi, height):
    x = np.linspace(lo, hi, 4001)
    return float(np.abs(moran(x + 1j * height, k)).min())


def winding(k, lo, hi, height, step):
    corners = [lo - 1j * height, hi - 1j * height, hi + 1j * height, lo + 1j * height]
    lipschitz = LN2 * 2 ** (-lo) + k * LN3 * 3 ** (-lo)
    turns, margin = 0.0, math.inf
    for a, b in zip(corners, corners[1:] + corners[:1]):
        n = math.ceil(abs(b - a) / step)
        path = a + (b - a) * np.arange(n + 1) / n
        values = moran(path, k)
        margin = min(
            margin, float((np.abs(values[:-1]) / (lipschitz * abs(b - a) / n)).min())
        )
        turns += float(np.angle(values[1:] / values[:-1]).sum())
    return round(turns / (2 * math.pi)), margin


def box(k):
    d, edge = dimension(k), left_edge(k)
    lo, hi = edge - 0.25, d + 0.25
    height = HEIGHT
    while edge_floor(k, lo, hi, height) < 0.02:
        height += 0.5
    found, from_lattice = roots_in_box(k, lo, hi, height)
    count, margin = winding(k, lo, hi, height, STEP)
    return d, edge, lo, hi, height, found, from_lattice, count, margin


def offsets(found):
    upper = found[found.imag > 1e-9]
    out = []
    for base in (2, 3, 6):
        omega = 2 * math.pi / math.log(base)
        ratio = upper.imag / omega
        out.append(
            math.floor(1000 * float(np.abs(ratio - np.round(ratio)).max())) / 1000
        )
    return out


def convergents():
    script = "\n".join(
        [
            f"k = {K}; D = solve(s = 0, 4, 2^-s + k*3^-s - 1); P = 2^-D; Q = k*3^-D; mu = P*log(2) + Q*log(3);",
            "mf(z) = 1 - 2^-z - k*3^-z;",
            "mg(z) = log(2)*2^-z + k*log(3)*3^-z;",
            f"c = contfrac(log(3)/log(2), , {TERMS}); p0 = 1; q0 = 0; p1 = c[1]; q1 = 1;",
            'for(j = 2, #c - 1, p2 = c[j]*p1 + p0; q2 = c[j]*q1 + q0; p0 = p1; q0 = q1; p1 = p2; q1 = q2; w = D + I*2*Pi*q1/log(2); for(i = 1, 80, w = w - mf(w)/mg(w)); th = 2*Pi*(q1*log(3)/log(2) - p1); pr = P*Q*log(2)^2*th^2/(2*mu^3); if(q1 > 1, print(q1, " ", p1, " ", imag(w), " ", D - real(w), " ", pr, " ", abs(mf(w)))))',
        ]
    )
    rows = [line.replace(" E", "E").split() for line in gp(script).strip().splitlines()]
    print(
        f"roots near D + 2 pi i q/ln 2 for the convergents p/q of log2(3), k = {K}: D - Re w against P Q (ln 2)^2 theta^2/(2 f'(D)^3), theta = 2 pi (q log2 3 - p)"
    )
    for q, p, im, gap, pred, res in rows:
        print(
            f"  q {int(q):>7}  p {int(p):>7}  Im w {float(im):>14.6f}  D - Re w {float(gap):.6e}  law {float(pred):.6e}  ratio {float(gap) / float(pred):.6f}  |f(w)| {float(res):.0e}"
        )


def poles():
    print(
        f"complex dimensions: the zeros of 1 - 2^-s - k 3^-s in the box Re [D_l - 0.25, D + 0.25], Im [-T, T], T >= {HEIGHT}"
    )
    print(
        f"seeds from PARI polroots of 1 - x^{DENOM} - k x^{NUMER} ({NUMER}/{DENOM} a convergent of log2 3), then Newton; the count certified by the winding number"
    )
    print(
        "   k        D        D_l      T   roots  winding  margin  from-lattice  Re min      Re max     offset 2pi/ln2 ln3 ln6, floored"
    )
    worst = 0.0
    for k in range(1, 28):
        d, edge, _, _, height, found, from_lattice, count, margin = box(k)
        o = offsets(found)
        print(
            f"  {k:>2}  {d:.6f}  {edge:+.6f}  {height:5.1f}  {found.size:>4}  {count:>5}  {margin:7.1f}  {from_lattice:>5}       {found.real.min():+.6f}  {found.real.max():+.6f}   {o[0]:.3f} {o[1]:.3f} {o[2]:.3f}"
        )
        worst = max(worst, found.real.min() - edge)
        if k == K:
            keep = found
        if k == 1:
            textbook = found[np.argmin(np.abs(found - TEXTBOOK))]
    print(
        f"largest gap from D_l to the lowest real part over k = 1..27, rounded up: {math.ceil(worst * 10000) / 10000:.4f}"
    )
    z = complex(TEXTBOOK)
    for _ in range(60):
        z -= (1 - 2**-z - 2 ** (-485 * z / 306)) / (
            LN2 * 2**-z + 485 / 306 * LN2 * 2 ** (-485 * z / 306)
        )
    print(
        f"k = 1, the 2-3 nonlattice equation: the printed root .7675115443 + 45.55415979 i; the root of its lattice approximant 1 - 2^-s - 2^(-485 s/306) {z.real:.10f} + {z.imag:.8f} i, at {abs(z - TEXTBOOK):.1e}; the true root {textbook.real:.10f} + {textbook.imag:.8f} i, at {abs(textbook - TEXTBOOK):.1e}"
    )
    print(f"the first complex dimensions of the patch, k = {K}, Im >= 0:")
    for w in keep[keep.imag > -1e-9]:
        print(f"  {w.real:+.6f} {w.imag:+.6f} i")
    clock("poles")
    convergents()
    clock("convergents")


def breakpoints(k, top):
    rows = []
    for b in range(int(top / LN3) + 1):
        for a in range(int((top - b * LN3) / LN2) + 1):
            rows.append((a * LN2 + b * LN3, a, b))
    rows.sort()
    exact = [2**a * 3**b for _, a, b in rows]
    ordered = all(x < y for x, y in pairwise(exact))
    weights = [comb(a + b, a) * k**b for _, a, b in rows]
    totals, running = [], 0
    for w in weights:
        running += w
        totals.append(running)
    return np.array([u for u, _, _ in rows]), exact, totals, ordered


def cells_at_least(exact, totals, x):
    i = bisect_right(exact, x)
    return totals[i - 1] if i else 0


def renewal_check(k, exact, totals):
    return all(
        totals[j]
        == 1
        + cells_at_least(exact, totals, Fraction(v, 2))
        + k * cells_at_least(exact, totals, Fraction(v, 3))
        for j, v in enumerate(exact)
    )


def normalized(u_breaks, log_totals, d, grid):
    i = np.searchsorted(u_breaks, grid, side="right") - 1
    return np.exp(log_totals[i] - d * grid)


def fold(u, g, period, nbins):
    phase = np.floor((u % period) / period * nbins).astype(int)
    means = np.bincount(phase, weights=g, minlength=nbins) / np.maximum(
        np.bincount(phase, minlength=nbins), 1
    )
    centred = g - g.mean()
    return 1 - float(((g - means[phase]) ** 2).sum() / (centred**2).sum())


def spectrum(u, g, top):
    y = (g - g.mean()) * np.blackman(g.size)
    pad = 1 << 21
    power = np.abs(np.fft.rfft(y, pad)) ** 2
    omega = 2 * np.pi * np.fft.rfftfreq(pad, u[1] - u[0])
    keep = (omega > 2.0) & (omega < OMEGA_MAX)
    w, pw = omega[keep], power[keep]
    peaks = np.where((pw[1:-1] > pw[:-2]) & (pw[1:-1] >= pw[2:]))[0] + 1
    peaks = peaks[np.argsort(pw[peaks])[::-1][:top]]
    amplitude = np.sqrt(pw[peaks]) / (np.blackman(g.size).sum() / 2)
    return w[peaks], amplitude


def roots_to(k, height):
    found, _ = roots_in_box(k, left_edge(k) - 0.25, dimension(k) + 0.25, height)
    return found[found.imag > 0]


def count():
    d = dimension(K)
    p, q = 2**-d, K * 3**-d
    mu = p * LN2 + q * LN3
    limit = 1 / (d * mu)
    u_breaks, exact, totals, ordered = breakpoints(K, UMAX)
    print(
        f"cells of size 2^-a 3^-b >= r = e^-U, U <= {UMAX:.0f}: {len(exact)} sizes, sorted exactly by 2^a 3^b: {ordered}"
    )
    print(
        f"N(e^-U) = sum of C(a+b, a) {K}^b; N({UMAX:.0f}) has {len(str(totals[-1]))} digits"
    )
    clock("breakpoints")
    print(
        f"renewal identity N(U) = 1 + N(U - ln 2) + {K} N(U - ln 3) exact at every breakpoint: {renewal_check(K, exact, totals)}"
    )
    clock("renewal")
    log_totals = np.array([math.log(t) for t in totals])
    print(
        f"limit N(r) r^D -> C = 1/(D f'(D)), f'(D) = 2^-D ln 2 + {K} 3^-D ln 3 = {mu:.9f}: {limit:.9f}"
    )
    print(
        "   window U     min F/C     max F/C     mean F/C - 1    swing      swing sqrt(U)   carpet swing"
    )
    dc = math.log(8) / LN3
    for start in WINDOW_STARTS:
        grid = np.linspace(start, start + WINDOW, 200001)
        f = normalized(u_breaks, log_totals, d, grid) / limit
        m = np.floor(grid / LN3)
        carpet = (8 ** (m + 1) - 1) / 7 * np.exp(-dc * grid)
        cswing = (carpet.max() - carpet.min()) / carpet.mean()
        swing = f.max() - f.min()
        print(
            f"  [{start:5.0f},{start + WINDOW:5.0f}]  {f.min():.6f}   {f.max():.6f}   {f.mean() - 1:+.6f}      {swing:.6f}   {swing * math.sqrt(start):.6f}       {cswing:.6f}"
        )
    clock("windows")
    u = np.arange(SPEC_LO, UMAX, SPEC_STEP)
    g = np.log(normalized(u_breaks, log_totals, d, u))
    m = np.floor(u / LN3)
    gc = np.log((8 ** (m + 1) - 1) / 7) - dc * u
    print(
        f"folded variance of g(u) = ln N(e^-u) - D u on [{SPEC_LO:.0f}, {UMAX:.0f}], {NBINS} bins:   ln 2    ln 3    ln 6"
    )
    print(
        f"  patch                                               {fold(u, g, LN2, NBINS):.3f}   {fold(u, g, LN3, NBINS):.3f}   {fold(u, g, math.log(6), NBINS):.3f}"
    )
    print(
        f"  carpet control                                      {fold(u, gc, LN2, NBINS):.3f}   {fold(u, gc, LN3, NBINS):.3f}   {fold(u, gc, math.log(6), NBINS):.3f}"
    )
    roots = roots_to(K, OMEGA_MAX + 5)
    print(
        f"Blackman periodogram of g on [{SPEC_LO:.0f}, {UMAX:.0f}], the {PEAKS} highest peaks in ({2.0}, {OMEGA_MAX:.0f}), each against the nearest complex dimension:"
    )
    w, a = spectrum(u, g, PEAKS)
    taper = np.blackman(u.size)
    for wi, ai in sorted(zip(w, a)):
        near = roots[np.argmin(np.abs(roots.imag - wi))]
        decay = float((taper * np.exp(-(d - near.real) * u)).sum() / taper.sum())
        predicted = 2 * decay / (abs(near * moran_prime(near, K)) * limit)
        print(
            f"  peak {wi:8.3f}  amplitude {ai:.3e}  root {near.real:.6f} + {near.imag:.3f} i  D - Re {d - near.real:.2e}  offset {wi - near.imag:+.3f}  residue amplitude {predicted:.3e}"
        )
    wc, _ = spectrum(u, gc, 6)
    print(
        "carpet control peaks, as multiples of 2 pi/ln 3:",
        " ".join(f"{x / (2 * math.pi / LN3):.3f}" for x in sorted(wc)),
    )
    clock("count")


def ball_mass(cx, cy, radius, d):
    ratios = np.array([float(r) for r, _ in CHILDREN])
    shifts = np.array([[float(x), float(y)] for _, (x, y) in CHILDREN])
    weights = ratios**d
    x, y, size, mass = np.zeros(1), np.zeros(1), np.ones(1), np.ones(1)
    inside, r2 = 0.0, radius * radius
    while size.size:
        far = (
            np.maximum(np.abs(x - cx), np.abs(x + size - cx)) ** 2
            + np.maximum(np.abs(y - cy), np.abs(y + size - cy)) ** 2
        )
        nx = np.clip(cx, x, x + size) - cx
        ny = np.clip(cy, y, y + size) - cy
        near = nx * nx + ny * ny
        full = far <= r2 * (1 - 1e-12)
        inside += float(mass[full].sum())
        open_ = ~full & (near < r2 * (1 + 1e-12))
        x, y, size, mass = x[open_], y[open_], size[open_], mass[open_]
        if size.size == 0 or size.max() < TOLERANCE * radius:
            return inside, inside + float(mass.sum())
        x = (x[:, None] + size[:, None] * shifts[None, :, 0]).ravel()
        y = (y[:, None] + size[:, None] * shifts[None, :, 1]).ravel()
        mass = (mass[:, None] * weights[None, :]).ravel()
        size = (size[:, None] * ratios[None, :]).ravel()
    return inside, inside


def ripple():
    d = dimension(K)
    print(
        f"the spin mass M(r) = mu(B(c, r)) of the natural measure, weights 2^-D and 3^-D, D = {d:.9f}, enclosed to relative {TOLERANCE:g} by cells"
    )
    print(
        "   centre    map ratio  window   identity M(r) = ratio^D M(r/ratio)   ripple swing   bar     swing/bar  fold ln 2   fold ln 3   fold ln 6"
    )
    for (cx, cy), ratio, window in RIPPLE_CENTRES:
        period = math.log(1 / ratio)
        steps = (
            np.arange(RIPPLE_PERIODS * RIPPLE_BINS + RIPPLE_BINS) * period / RIPPLE_BINS
        )
        radii = window * 0.999 * np.exp(-steps)
        bounds = np.array([ball_mass(cx, cy, r, d) for r in radii])
        mid = bounds.mean(axis=1)
        bar = float(np.max(np.log(bounds[:, 1] / bounds[:, 0])))
        lo, hi = bounds[RIPPLE_BINS:, 0], bounds[RIPPLE_BINS:, 1]
        scaled_lo, scaled_hi = (
            ratio**d * bounds[:-RIPPLE_BINS, 0],
            ratio**d * bounds[:-RIPPLE_BINS, 1],
        )
        meets = bool(
            np.all((lo <= scaled_hi * (1 + 1e-9)) & (scaled_lo <= hi * (1 + 1e-9)))
        )
        g = np.log(mid) - d * np.log(radii)
        logr = np.log(radii)
        folds = [fold(logr, g, math.log(b), RIPPLE_BINS) for b in (2, 3, 6)]
        print(
            f"   ({cx},{cy})   1/{round(1 / ratio)}       {window:.4f}   bands meet at all {lo.size} shifted radii: {meets}          {g.max() - g.min():.5f}   {bar:.1e}   {math.floor((g.max() - g.min()) / bar):>6}     {folds[0]:.3f}       {folds[1]:.3f}       {folds[2]:.3f}"
        )
    clock("ripple")


VERBS = {"patch": patch, "poles": poles, "count": count, "ripple": ripple}

if __name__ == "__main__":
    VERBS[sys.argv[1] if len(sys.argv) > 1 else "patch"]()
