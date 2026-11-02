import time
from fractions import Fraction
from math import atan2, gcd, sqrt
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent

N = 55
ODDS = list(range(1, N + 1, 2))
L = len(ODDS)
GOLDEN = 137.507764
RESOLUTIONS = [256, 512, 1024, 2048]
LAYER_COUNTS = [4, 8, 14, 28]
FIG_R = 1024
PEAK_COUNT = 20
PEAK_RADIUS = 0.01
ZOOM_R = 1024
ZOOM_W = 0.02


def primes_up_to_count(k):
    out = []
    c = 2
    while len(out) < k:
        if all(c % d for d in range(2, int(sqrt(c)) + 1)):
            out.append(c)
        c += 1
    return out


def gaussian_angle(n):
    best = None
    for b in range(1, int(sqrt(n)) + 1):
        a2 = n - b * b
        a = round(sqrt(a2))
        if a * a == a2 and a >= b >= 1:
            best = (a, b)
    if best is None:
        return 0.0
    return np.degrees(atan2(best[1], best[0]))


def schedules():
    ps = primes_up_to_count(L - 1)
    return {
        "unspun": [0.0] * L,
        "degrees": [float(k + 1) for k in range(L)],
        "golden": [GOLDEN * (k + 1) for k in range(L)],
        "primes": [0.0] + [float(p) for p in ps],
        "random": list(np.random.default_rng(1).uniform(0.0, 360.0, L)),
        "gaussian": [gaussian_angle(n) for n in ODDS],
    }


def parity(z):
    return z - np.floor(z) >= 0.5


def stack(angles, r, layers=L, dtype=np.float32, box=(0.0, 0.0, 1.0)):
    x0, y0, w = box
    step = (np.arange(r, dtype=dtype) + dtype(0.5)) * dtype(w / r)
    dx = step + dtype(x0 - 0.5)
    dy = step + dtype(y0 - 0.5)
    acc = np.zeros((r, r), dtype=np.uint8)
    for k in range(layers):
        n = ODDS[k]
        h = 0.5 * n
        t = np.float64(angles[k]) * np.pi / 180.0
        cx = dtype(h * np.cos(t))
        sx = dtype(h * np.sin(t))
        c0 = dtype(h * 0.5)
        a = c0 + cx * dx[None, :] + sx * dy[:, None]
        b = c0 - sx * dx[None, :] + cx * dy[:, None]
        acc += parity(a) & parity(b)
    return (np.float32(layers) - acc) / np.float32(layers)


def disc_mask(r):
    d = (np.arange(r, dtype=np.float64) + 0.5) / r - 0.5
    return d[:, None] ** 2 + d[None, :] ** 2 <= 0.25


def maximum_filter3(f):
    p = np.pad(f, 1, mode="constant", constant_values=-np.inf)
    r = f.shape[0]
    m = f.copy()
    for dy in (0, 1, 2):
        for dx in (0, 1, 2):
            np.maximum(m, p[dy : dy + r, dx : dx + r], out=m)
    return m


def plateau(field, val, iy, ix, wide):
    r = field.shape[0]
    y0, y1 = max(0, iy - wide), min(r, iy + wide + 1)
    x0, x1 = max(0, ix - wide), min(r, ix + wide + 1)
    w = field[y0:y1, x0:x1] == val
    seed = np.zeros_like(w)
    seed[iy - y0, ix - x0] = True
    while True:
        g = seed.copy()
        g[1:] |= seed[:-1]
        g[:-1] |= seed[1:]
        g[:, 1:] |= seed[:, :-1]
        g[:, :-1] |= seed[:, 1:]
        g &= w
        if g.sum() == seed.sum():
            return np.nonzero(seed)[0] + y0, np.nonzero(seed)[1] + x0
        seed = g


def peaks(field, mask, count=PEAK_COUNT, radius=PEAK_RADIUS):
    r = field.shape[0]
    cand = (field >= maximum_filter3(field)) & mask
    ys, xs = np.nonzero(cand)
    vs = field[ys, xs]
    order = np.lexsort((xs, ys, -vs))
    ys, xs, vs = ys[order], xs[order], vs[order]
    py = (ys + 0.5) / r
    px = (xs + 0.5) / r
    alive = np.ones(len(vs), dtype=bool)
    wide = int(radius * r) + 1
    out = []
    for _ in range(count):
        idx = int(np.argmax(alive))
        if not alive[idx]:
            break
        cy, cx = plateau(field, vs[idx], ys[idx], xs[idx], wide)
        out.append(
            (
                float((cx + 0.5).mean() / r),
                float((cy + 0.5).mean() / r),
                float(vs[idx]),
                float(px[idx]),
                float(py[idx]),
                len(cx),
            )
        )
        alive &= ~((px - px[idx]) ** 2 + (py - py[idx]) ** 2 <= radius**2)
    return out


def drift(a, b, r):
    pa = np.array([[p[0], p[1]] for p in a])
    pb = np.array([[p[0], p[1]] for p in b])
    dd = np.sqrt(((pa[:, None, :] - pb[None, :, :]) ** 2).sum(2)).min(1)
    return float(np.median(dd)), float(dd.max()), int((dd * r <= 1.0).sum())


def exact_variance(layers):
    s = [2 * k + 1 for k in range(layers)]
    tot = Fraction(0)
    for m in s:
        for n in s:
            d = gcd(m, n)
            tot += Fraction(
                (d * d - 1) * (2 * (m - 1) * (n - 1) + d * d - 1), 16 * m * m * n * n
            )
    return tot / (layers * layers)


def diagonal_maximum():
    cuts = sorted(
        {Fraction(k, n) for n in ODDS for k in range(1, n + 1)}
        | {Fraction(0), Fraction(1)}
    )
    best = 0
    rows = []
    for i in range(len(cuts) - 1):
        u = (cuts[i] + cuts[i + 1]) / 2
        sel = frozenset(
            n for n in ODDS if (n * u).numerator // (n * u).denominator % 2 == 1
        )
        if len(sel) > best:
            best, rows = len(sel), []
        if len(sel) == best:
            rows.append((cuts[i], cuts[i + 1], sel))
    cells = [(a, b, c, d) for (a, b, s1) in rows for (c, d, s2) in rows if s1 == s2]
    area = sum(((b - a) * (d - c) for (a, b, c, d) in cells), Fraction(0))
    return best, rows, cells, area, len(cuts)


def dark_layers(angles, x, y):
    out = []
    for k in range(L):
        n = ODDS[k]
        t = np.float64(angles[k]) * np.pi / 180.0
        u = 0.5 + (x - 0.5) * np.cos(t) + (y - 0.5) * np.sin(t)
        v = 0.5 - (x - 0.5) * np.sin(t) + (y - 0.5) * np.cos(t)
        if parity(0.5 * n * u) and parity(0.5 * n * v):
            out.append(n)
    return out


def centre_value(field, r):
    i = r // 2
    return float(field[i - 1 : i + 1, i - 1 : i + 1].mean())


def save_png(field, path, side):
    r = field.shape[0]
    b = r // side
    q = field.reshape(side, b, side, b).mean((1, 3))
    q = np.round(np.round(q * L) * (255.0 / L)).astype(np.uint8)
    Image.fromarray(q, mode="L").save(path, optimize=True)
    return path.stat().st_size


def contact_sheet(fields, path, side):
    sheet = np.zeros((2 * side, 2 * side), dtype=np.uint8)
    for i, f in enumerate(fields):
        r = f.shape[0]
        b = r // side
        q = f.reshape(side, b, side, b).mean((1, 3))
        q = np.round(np.round(q * L) * (255.0 / L)).astype(np.uint8)
        sheet[
            side * (i // 2) : side * (i // 2) + side,
            side * (i % 2) : side * (i % 2) + side,
        ] = q
    Image.fromarray(sheet, mode="L").save(path, optimize=True)
    return path.stat().st_size


def frac(v):
    return "%d/%d" % (round(v * L), L)


def report(name, angles):
    print("schedule", name, "angles", " ".join("%.6g" % a for a in angles))
    tops = {}
    for r in RESOLUTIONS:
        f = stack(angles, r)
        msk = disc_mask(r)
        ink = 1.0 - f
        d = ink[msk]
        pk = peaks(ink, msk)
        tops[r] = pk
        top = pk[0][2]
        area = float((d == d.max()).sum()) / (r * r)
        print(
            " R %4d peak %s locus %.6g mean %.6f rms %.6f centre %s top"
            % (
                r,
                frac(top),
                area,
                float(f[msk].mean()),
                float(d.std()),
                frac(1.0 - centre_value(f, r)),
            ),
            " ".join(
                "(%.5f,%.5f,%s,%dpx)" % (q[0], q[1], frac(q[2]), q[5]) for q in pk[:3]
            ),
        )
    for r in RESOLUTIONS[:-1]:
        med, mx, hit = drift(tops[r], tops[2 * r], r)
        print(
            " drift %4d->%4d median %.6f (%.2f px) max %.6f (%.2f px) within 1 px %d/%d"
            % (r, 2 * r, med, med * r, mx, mx * r, hit, PEAK_COUNT)
        )
    x, y, v, sx, sy, area = tops[RESOLUTIONS[-1]][0]
    z = stack(angles, ZOOM_R, box=(x - ZOOM_W / 2, y - ZOOM_W / 2, ZOOM_W))
    zi = 1.0 - z
    lev = float((zi == zi.max()).mean())
    print(
        " zoom window %g wide at effective R %d: peak %s, cell at the top peak %.6g"
        % (ZOOM_W, int(ZOOM_R / ZOOM_W), frac(float(zi.max())), lev * ZOOM_W * ZOOM_W)
    )
    dk = dark_layers(angles, sx, sy)
    print(
        " top peak (%.5f, %.5f) at %s, plateau %d px at R %d, %d layers dark there"
        % (x, y, frac(v), area, RESOLUTIONS[-1], len(dk)),
        " ".join(str(n) for n in dk),
    )
    return tops


def main():
    t0 = time.perf_counter()
    sched = schedules()
    print(
        "layers",
        L,
        "odd scales 1..%d" % N,
        "peak radius",
        PEAK_RADIUS,
        "peak count",
        PEAK_COUNT,
    )
    print(
        "centre cell inradius 1/%d in every layer at every angle, scales dark at the centre"
        % (2 * N),
        " ".join(str(n) for n in dark_layers([0.0] * L, 0.5, 0.5)),
    )
    for name, angles in sched.items():
        report(name, angles)

    c, rows, cells, area, cuts = diagonal_maximum()
    print(
        "unspun ink maximum by exact scan of %d breakpoints on the diagonal: %d/%d on %d intervals"
        % (cuts, c, L, len(rows)),
        " ".join("[%s, %s)" % (a, b) for a, b, _ in rows),
    )
    print(
        "unspun maximum locus: %d cells, total area %s = %.6g, one cell %s = %.6g, scales"
        % (
            len(cells),
            area,
            float(area),
            (rows[0][1] - rows[0][0]) ** 2,
            float((rows[0][1] - rows[0][0]) ** 2),
        ),
        " ".join(str(n) for n in sorted(rows[0][2])),
    )

    print("fade at R %d, c = rms * sqrt(layers)" % FIG_R)
    msk = disc_mask(FIG_R)
    for name, angles in sched.items():
        row = []
        for lc in LAYER_COUNTS:
            f = stack(angles, FIG_R, lc)
            row.append((float(f.std()) * sqrt(lc), float(f[msk].std()) * sqrt(lc)))
        print(
            " %-9s square" % name,
            " ".join("%.6f" % a for a, b in row),
            " disc",
            " ".join("%.6f" % b for a, b in row),
        )
    print(
        " exact   square",
        " ".join("%.6f" % sqrt(float(exact_variance(lc)) * lc) for lc in LAYER_COUNTS),
    )

    f32 = stack(sched["golden"], 512)
    f64 = stack(sched["golden"], 512, dtype=np.float64)
    print(
        "float32 against float64 at R 512 golden: max abs diff %.6g, pixels differing %d of %d"
        % (float(np.abs(f32 - f64).max()), int((f32 != f64).sum()), 512 * 512)
    )

    figs = {}
    for name, fname in (
        ("unspun", "stack-unspun.png"),
        ("degrees", "stack-degrees.png"),
        ("primes", "stack-primes.png"),
        ("gaussian", "stack-gaussian.png"),
    ):
        f = stack(sched[name], FIG_R)
        f = np.where(disc_mask(FIG_R), f, 1.0).astype(np.float32)
        figs[name] = f
        print("figure", fname, "bytes", save_png(f, HERE / fname, 512))
    print(
        "figure stack-sheet.png bytes",
        contact_sheet(list(figs.values()), HERE / "stack-sheet.png", 256),
    )
    print("seconds", round(time.perf_counter() - t0, 2))


main()
