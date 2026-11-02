import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from mpmath import asin, exp, log, mp, mpf, pi, quad, sqrt
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sponge_tube as st

LEVEL = 4
SIDE = 3**LEVEL
REMOVED = np.array(
    [(1, 1, 1), (0, 1, 1), (2, 1, 1), (1, 0, 1), (1, 2, 1), (1, 1, 0), (1, 1, 2)]
)
PLUS = 7 / 27
SQRT2 = np.sqrt(2)

# THE LEVEL-4 PREFRACTAL


def prefractal():
    idx = np.arange(SIDE)
    digits = np.stack([(idx // 3**k) % 3 for k in range(LEVEL)], 1)
    one = digits == 1
    i, j, k = np.meshgrid(idx, idx, idx, indexing="ij")
    count = one[i].astype(int) + one[j].astype(int) + one[k].astype(int)
    keep = (count <= 1).all(axis=-1)
    return np.stack([i[keep], j[keep], k[keep]], 1)


class Oracle:
    def __init__(self, cubes):
        h = 1.0 / SIDE
        self.lo = cubes * h
        self.hi = self.lo + h
        self.tree = cKDTree(self.lo + h / 2)
        self.halfdiag = np.sqrt(3) * h / 2

    def dist(self, points, k=300):
        d, ii = self.tree.query(points, k=k)
        q = points[:, None, :]
        g = np.maximum(np.maximum(self.lo[ii] - q, q - self.hi[ii]), 0)
        best = np.sqrt((g * g).sum(-1)).min(1)
        complete = d[:, -1] - self.halfdiag > best
        return best, complete


def wall_faces(cubes):
    n = SIDE
    faces = []
    for a in range(3):
        for arm in (0, 2):
            for b in ((a + 1) % 3, (a + 2) % 3):
                c = 3 - a - b
                for side in (0, 1):
                    plane = 1 / 3 if side == 0 else 2 / 3
                    sel = (
                        (cubes[:, b] == n // 3 - 1)
                        if side == 0
                        else (cubes[:, b] == 2 * n // 3)
                    )
                    sel &= (cubes[:, a] >= arm * n // 3) & (
                        cubes[:, a] < (arm + 1) * n // 3
                    )
                    sel &= (cubes[:, c] >= n // 3) & (cubes[:, c] < 2 * n // 3)
                    faces.append((b, plane, sel))
    return faces


def dist_to_faces(points, oracle, faces, own):
    c = np.floor(points * 3).astype(int)
    centre = (c == 1).all(1)
    best = np.full(len(points), np.inf)
    for b, plane, sel in faces:
        lo, hi = oracle.lo[sel], oracle.hi[sel]
        others = [k for k in range(3) if k != b]
        q = points[:, None, :]
        g = np.maximum(
            np.maximum(
                lo[None, :, others] - q[:, :, others],
                q[:, :, others] - hi[None, :, others],
            ),
            0,
        )
        dd = np.sqrt((g * g).sum(-1) + (points[:, None, b] - plane) ** 2).min(1)
        if own:
            inside = np.ones(len(points), bool)
            for k in others:
                inside &= (points[:, k] >= lo[:, k].min() - 1e-12) & (
                    points[:, k] <= hi[:, k].max() + 1e-12
                )
            inside &= np.abs(points[:, b] - plane) <= 1 / 3 + 1e-12
            best = np.where(inside & ~centre, np.minimum(best, dd), best)
        else:
            best = np.minimum(best, dd)
    if own:
        m = np.minimum(points[centre] - 1 / 3, 2 / 3 - points[centre])
        best[centre] = np.minimum.reduce(
            [
                np.hypot(m[:, 0], m[:, 1]),
                np.hypot(m[:, 0], m[:, 2]),
                np.hypot(m[:, 1], m[:, 2]),
            ]
        )
    return best


# THE REDUCED DISTANCE OF THE LEMMAS


def carpet_dist(v, z, side, depth=40):
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


def reduced(points):
    c = np.floor(points * 3).astype(int)
    centre = (c == 1).all(1)
    out = np.full(len(points), np.inf)
    for a in range(3):
        arm = (c[:, a] != 1) & (c[:, (a + 1) % 3] == 1) & (c[:, (a + 2) % 3] == 1)
        along = points[arm, a] - c[arm, a] / 3
        for b in ((a + 1) % 3, (a + 2) % 3):
            other = (a + 2) % 3 if b == (a + 1) % 3 else (a + 1) % 3
            u = points[arm, b] - 1 / 3
            dd = carpet_dist(points[arm, other] - 1 / 3, along, 1 / 3)
            out[arm] = np.minimum(out[arm], np.sqrt(u * u + dd * dd))
            out[arm] = np.minimum(out[arm], np.sqrt((1 / 3 - u) ** 2 + dd * dd))
    m = np.minimum(points[centre] - 1 / 3, 2 / 3 - points[centre])
    out[centre] = np.minimum.reduce(
        [
            np.hypot(m[:, 0], m[:, 1]),
            np.hypot(m[:, 0], m[:, 2]),
            np.hypot(m[:, 1], m[:, 2]),
        ]
    )
    return out


def sample_plus(rng, n):
    base = REMOVED[rng.integers(0, 7, n)] / 3.0
    return base + rng.random((n, 3)) / 3.0


def three_sigma(hits, n, scale):
    p = hits / n
    return p * scale, 3 * np.sqrt(p * (1 - p) / n) * scale


# THE ORACLE VERB


def verb_oracle():
    t0 = time.time()
    rng = np.random.default_rng(7)
    cubes = prefractal()
    oracle = Oracle(cubes)
    faces = wall_faces(cubes)
    print(
        f"ORACLE: level-{LEVEL} prefractal, {len(cubes)} cubes, {sum(int(f[2].sum()) for f in faces)} cube faces on the 24 walls"
    )
    points = sample_plus(rng, 20000)
    near = sample_plus(rng, 4000)
    axis = rng.integers(0, 3, 4000)
    near[np.arange(4000), axis] = np.where(
        rng.random(4000) < 0.5, 1 / 3, 2 / 3
    ) + rng.normal(0, 1e-3, 4000)
    near = np.clip(near, 0, 1)
    c = np.minimum(np.floor(near * 3).astype(int), 2)
    near = near[(c == 1).sum(1) >= 2]
    points = np.vstack([points, near])
    d_f4, complete = oracle.dist(points)
    d_walls = dist_to_faces(points, oracle, faces, own=False)
    d_own = dist_to_faces(points, oracle, faces, own=True)
    print(
        f"  {len(points)} points of the plus, {len(near)} of them within about 1e-3 of a face; candidate set complete for every point: {bool(complete.all())}"
    )
    print(
        f"  max |dist(x, F_4) - dist(x, 24 level-4 wall carpets)| = {np.abs(d_f4 - d_walls).max():.1e}"
    )
    print(
        f"  max |dist(x, F_4) - own four walls or centre edges at level 4| = {np.abs(d_f4 - d_own).max():.1e}"
    )
    d_red = reduced(points)
    print(
        f"  reduced distance of the lemmas (infinite carpets) minus dist(x, F_4): min {(d_red - d_f4).min():.1e}, max {(d_red - d_f4).max():.2e}, against sqrt(2)/{3**LEVEL} = {SQRT2 / 3**LEVEL:.2e}"
    )
    centre = oracle.dist(np.array([[0.5, 0.5, 0.5]]))[0][0]
    print(
        f"  dist(centre, F_4) = {centre:.17f}, sqrt(2)/6 = {SQRT2 / 6:.17f}, difference {abs(centre - SQRT2 / 6):.1e}"
    )
    print(f"  wall {time.time() - t0:.1f} s")


# THE MONTE CARLO VERB


def verb_montecarlo():
    t0 = time.time()
    rng = np.random.default_rng(7)
    cubes = prefractal()
    oracle = Oracle(cubes)
    print(
        "MONTECARLO: T(delta) from the reduced distance of the lemmas, and a lemma-free bracket from the level-4 prefractal"
    )
    n, chunk = 1_000_000, 100_000
    margin = SQRT2 / (6 * SIDE)
    low = {Fr(1, 8): 0, Fr(1, 12): 0}
    high = {Fr(1, 8): 0, Fr(1, 12): 0}
    for _ in range(n // chunk):
        d_f4, complete = oracle.dist(sample_plus(rng, chunk), k=200)
        assert complete.all()
        for d in low:
            high[d] += int((d_f4 <= float(d)).sum())
            low[d] += int((d_f4 <= float(d) - margin).sum())
    for d, lo in low.items():
        print(
            f"  lemma-free, {n} points against F_4: T({d}) in about [{lo / n * PLUS:.5f}, {high[d] / n * PLUS:.5f}], the count within delta an upper bound on T since F lies in F_4, the count within delta - sqrt(2)/486 a lower bound since every point of a kept cube is within sqrt(2)/486 of F, Monte Carlo 3 sigma {3 * np.sqrt(0.25 / n) * PLUS:.5f}"
        )
    n, chunk = 40_000_000, 2_000_000
    hits = {Fr(1, 12): 0, Fr(1, 8): 0, Fr(1, 6): 0}
    for _ in range(n // chunk):
        d_red = reduced(sample_plus(rng, chunk))
        for d in hits:
            hits[d] += int((d_red <= float(d)).sum())
    for d, h in hits.items():
        value, err = three_sigma(h, n, PLUS)
        print(
            f"  reduced distance, {n} points: T({d}) = {value:.6f} +- {err:.6f} (3 sigma)"
        )
    delta = Fr(1, 12)
    d = float(delta)
    a1 = st.lower(st.strip(delta)[0])
    n, chunk = 20_000_000, 5_000_000
    hits = 0
    for _ in range(n // chunk):
        u, v, z = rng.random(chunk) * d, rng.random(chunk) * d, rng.random(chunk) / 3
        du, dv = carpet_dist(v, z, 1 / 3), carpet_dist(u, z, 1 / 3)
        hits += int(((u * u + du * du <= d * d) & (v * v + dv * dv <= d * d)).sum())
    value, err = three_sigma(hits, n, d * d / 3)
    print(
        f"  overlap of the two wall tubes, {n} points of the column: V2({delta}) = {value:.8f} +- {err:.1e}, against 2 A1 - delta^2/3 = {float(2 * a1 - mpf(d) ** 2 / 3):.8f} with A1 from the generator"
    )
    print(f"  wall {time.time() - t0:.1f} s")


# THE SEEDED VERB


def verb_seeded():
    t0 = time.time()
    rng = np.random.default_rng(20260921)
    print(
        "SEEDED: Deep sampled on its confinement box [delta - w, delta]^2 x [0, 1/3], w = (1/81)/(4 delta), and V2 on the column, seed 20260921"
    )
    for delta in (Fr(1, 6), Fr(1, 8)):
        d = float(delta)
        w = (1 / 81) / (4 * d)
        bound = float(st.deep_bound(delta))
        n, chunk = 10_000_000, 1_000_000
        deep_hits = 0
        v2_hits = 0
        for _ in range(n // chunk):
            u, v, z = (
                d - w * rng.random(chunk),
                d - w * rng.random(chunk),
                rng.random(chunk) / 3,
            )
            du, dv = carpet_dist(v, z, 1 / 3), carpet_dist(u, z, 1 / 3)
            deep_hits += int(
                ((u * u + du * du > d * d) & (v * v + dv * dv > d * d)).sum()
            )
            u, v, z = (
                d * rng.random(chunk),
                d * rng.random(chunk),
                rng.random(chunk) / 3,
            )
            du, dv = carpet_dist(v, z, 1 / 3), carpet_dist(u, z, 1 / 3)
            v2_hits += int(
                ((u * u + du * du <= d * d) & (v * v + dv * dv <= d * d)).sum()
            )
        deep, deep_err = three_sigma(deep_hits, n, w * w / 3)
        v2, v2_err = three_sigma(v2_hits, n, d * d / 3)
        a1 = st.lower(st.strip(delta)[0])
        ident = float(2 * a1 - mpf(d) ** 2 / 3) + deep
        print(
            f"  delta = {delta}: Deep = {deep:.4e} +- {deep_err:.1e} (3 sigma, {deep_hits} hits of {n}), bound {bound:.4e}, ratio {bound / deep:.1f}"
        )
        print(
            f"  delta = {delta}: V2 = {v2:.7f} +- {v2_err:.1e} (3 sigma), against 2 A1 - delta^2/3 + Deep = {ident:.7f}, difference {v2 - ident:.1e}"
        )
    print(f"  wall {time.time() - t0:.1f} s")


# THE RECOMPUTE VERB, OWN CLOSED FORMS AT 60 DIGITS


def R(x):
    return mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mpf(x)


def a0(tau, d):
    b = sqrt(d * d - tau * tau)
    return (tau * b + d * d * asin(tau / d)) / 2


def a1(tau, d):
    b = sqrt(d * d - tau * tau)
    return tau * tau * (d * d + d * b + b * b) / (3 * (d + b))


def a1_naive(tau, d):
    return (d**3 - (d * d - tau * tau) ** mpf(1.5)) / 3


def hole(s, d, prim=a1):
    tau = R(min(d, s / 2))
    s, d = R(s), R(d)
    return 4 * s * a0(tau, d) - 8 * prim(tau, d)


def cut_hole(s, c, d, prim=a1):
    if c <= s / 2:
        tau = R(min(c, d))
        s, c, d = R(s), R(c), R(d)
        return (s + 2 * c) * a0(tau, d) - 4 * prim(tau, d)
    t1, t2 = R(min(s - c, d)), R(min(s / 2, d))
    s, c, d = R(s), R(c), R(d)
    return (
        (s + 2 * c) * a0(t1, d)
        - 4 * prim(t1, d)
        + 4 * s * (a0(t2, d) - a0(t1, d))
        - 8 * (prim(t2, d) - prim(t1, d))
    )


def column_holes(i, m):
    n = 1
    for _ in range(m - 1):
        n *= 2 if i % 3 == 1 else 3
        i //= 3
    return n


def columns_below(count, digits):
    if digits == 0:
        return 1 if count >= 1 else 0
    q, rem = divmod(count, 3 ** (digits - 1))
    if q >= 3:
        return 8**digits
    w = (3, 2, 3)
    return sum(w[:q]) * 8 ** (digits - 1) + w[q] * columns_below(rem, digits - 1)


def wall_half(d, levels, prim=a1):
    return (
        sum(
            8 ** (m - 1) * hole(Fr(1, 3 ** (m + 1)), d, prim)
            for m in range(1, levels + 1)
        )
        / 2
    )


def strip(d, levels, prim=a1):
    total = mpf(0)
    for m in range(1, levels + 1):
        s = Fr(1, 3 ** (m + 1))
        count = (d / s - 2) // 3 + 1 if d / s >= 2 else 0
        count = min(count, 3 ** (m - 1))
        if count > 0:
            total += columns_below(count, m - 1) * hole(s, d, prim)
        i = (d / s - 1) // 3
        if 0 <= i < 3 ** (m - 1) and (3 * i + 1) * s < d < (3 * i + 2) * s:
            total += column_holes(i, m) * cut_hole(s, d - (3 * i + 1) * s, d, prim)
    return total


def tube_open(d, levels):
    x = R(d)
    return (
        (pi + 8) * x**2
        - 8 * sqrt(2) * x**3
        + 48 * (wall_half(d, levels) - strip(d, levels))
    )


def tube_upper(d):
    x = R(d)
    total = mpf(0)
    for m in range(1, 200):
        s = Fr(1, 3 ** (m + 1))
        total += 8 ** (m - 1) * R(s) * x * R(min(s, 4 * d))
    return pi * x**2 + 24 * (total + x * R(Fr(8, 9) ** 199 / 8))


def hole_tables(levels):
    cells = 3**levels
    tables = {}
    for m in range(1, levels + 1):
        span = 3 ** (levels - m)
        table = np.zeros((3 ** (m - 1), cells))
        for i in range(3 ** (m - 1)):
            rows = [0]
            for k in range(m - 1):
                allowed = (0, 2) if (i // 3**k) % 3 == 1 else (0, 1, 2)
                rows = [r + a * 3**k for r in rows for a in allowed]
            for r in rows:
                table[i, (3 * r + 1) * span : (3 * r + 2) * span] = 1.0
        tables[m] = table
    return tables


def deep_bound(delta, tables, levels=7, tail_levels=80):
    delta = float(delta)
    cells = 3**levels
    boxes = 0.0
    for m in range(1, levels + 1):
        sm = 3.0 ** (-(m + 1))
        wm = sm * sm / (4 * delta)
        for n in range(1, levels + 1):
            sn = 3.0 ** (-(n + 1))
            wn = sn * sn / (4 * delta)
            i = np.arange(3 ** (m - 1))
            vlen = np.clip(
                np.minimum((3 * i + 2) * sm, delta)
                - np.maximum((3 * i + 1) * sm, delta - wn),
                0,
                None,
            )
            k = np.arange(3 ** (n - 1))
            ulen = np.clip(
                np.minimum((3 * k + 2) * sn, delta)
                - np.maximum((3 * k + 1) * sn, delta - wm),
                0,
                None,
            )
            boxes += float((vlen @ tables[m]) @ (ulen @ tables[n])) / cells / 3
    tail = 0.0
    for m in range(1, tail_levels + 1):
        for n in range(1, tail_levels + 1):
            if max(m, n) <= levels:
                continue
            wm = min(3.0 ** (-(2 * m + 2)) / (4 * delta), delta)
            wn = min(3.0 ** (-(2 * n + 2)) / (4 * delta), delta)
            count = min(np.ceil(wn * 3.0**m) + 2, np.ceil(wm * 3.0**n) + 2)
            tail += wm * wn * min(count / 9, 1 / 3)
    tail += 3.0 ** (-(2 * tail_levels + 4)) / (1536 * delta * delta)
    return boxes + tail


def verb_recompute():
    t0 = time.time()
    mp.dps = 60
    levels = 200
    print(
        "RECOMPUTE: the tube and the bands from the closed forms of the paper alone, 60 decimal digits, 200 hole levels"
    )
    for m in range(1, 7):
        digits = m - 1
        brute = [
            sum(
                1
                for r in range(3**digits)
                if all(
                    not ((i // 3**k) % 3 == 1 and (r // 3**k) % 3 == 1)
                    for k in range(digits)
                )
            )
            for i in range(3**digits)
        ]
        assert brute == [column_holes(i, m) for i in range(3**digits)], m
        assert all(
            columns_below(count, digits) == sum(brute[:count])
            for count in range(3**digits + 1)
        ), m
    print("  column counts of the strip validated by brute force to level 6")
    for s, d in ((Fr(1, 27), Fr(1, 12)), (Fr(1, 9), Fr(1, 6))):
        q = quad(
            lambda t: 4 * (R(s) - 2 * t) * sqrt(R(d) ** 2 - t**2),
            [0, min(R(d), R(s) / 2)],
        )
        print(
            f"  J({s}, {d}) closed form minus quadrature: {mp.nstr(hole(s, d) - q, 3)}"
        )
    mp.dps = 20
    for s, c, d in (
        (Fr(1, 9), Fr(1, 72), Fr(1, 8)),
        (Fr(1, 9), Fr(7, 108), Fr(1, 6)),
        (Fr(1, 81), Fr(1, 200), Fr(1, 12)),
    ):
        s_, c_, d_ = R(s), R(c), R(d)

        def inner(v):
            def f(z):
                r = min(v, s_ - v, z, s_ - z)
                return sqrt(d_**2 - r**2) if r < d_ else mpf(0)

            return quad(f, [0, min(v, s_ - v), s_ / 2, s_ - min(v, s_ - v), s_])

        q = quad(inner, [0, min(c_, s_ / 2), c_] if c_ > s_ / 2 else [0, c_])
        print(
            f"  cut hole (s, c, delta) = ({s}, {c}, {d}) closed form minus 2D quadrature: {mp.nstr(cut_hole(s, c, d) - q, 3)}"
        )
    mp.dps = 60
    tail_v1 = lambda d: R(d) * R(Fr(8, 9) ** levels) / 18
    tail_a1 = lambda d: R(d) * R(Fr(8, 9) ** levels) / 9
    tables = hole_tables(7)
    for d in (Fr(1, 12), Fr(1, 8), Fr(1, 6)):
        v1, s1 = wall_half(d, levels), strip(d, levels)
        t = tube_open(d, levels)
        deep = deep_bound(d, tables)
        print(
            f"  delta = {d}: V1 = {mp.nstr(v1, 15)}, A1 = {mp.nstr(s1, 15)}, Deep bound {deep:.4e}, T in [{mp.nstr(t - 48 * tail_a1(d) - 24 * mpf(deep), 12)}, {mp.nstr(t + 48 * tail_v1(d), 12)}]"
        )
    weight = mpf(27) / 20
    dim = log(20) / log(3)
    for eps in (Fr(1, 12), Fr(1, 8), Fr(1, 6)):
        lo = hi = mpf(20) / 27
        for l in range(41):
            d = eps / 3**l
            t = tube_open(d, levels)
            deep = deep_bound(d, tables)
            lo += weight**l * (t - 48 * tail_a1(d) - 24 * mpf(deep))
            hi += weight**l * (t + 48 * tail_v1(d))
        series_tail = weight**41 * tube_upper(eps / 3**41) * mpf(20) / 9
        hi += series_tail
        scale = exp((dim - 3) * log(R(eps)))
        print(
            f"  eps = {eps}: p in [{mp.nstr(lo * scale, 12)}, {mp.nstr(hi * scale, 12)}], series tail {mp.nstr(series_tail, 3)}"
        )
    mp.dps = 200
    d = Fr(1, 12)
    v1, s1 = wall_half(d, levels, a1_naive), strip(d, levels, a1_naive)
    mp.dps = 60
    print(
        f"  delta = 1/12 at 200 digits with a1 = (delta^3 - (delta^2 - t^2)^(3/2))/3: V1 = {mp.nstr(v1, 15)}, A1 = {mp.nstr(s1, 15)}"
    )
    n = 3000
    g = (np.arange(n) + 0.5) / (3 * n)
    v, z = np.meshgrid(g, g, indexing="ij")
    f = np.sqrt(np.maximum(float(d) ** 2 - carpet_dist(v, z, 1 / 3) ** 2, 0))
    cell = (1 / (3 * n)) ** 2
    print(
        f"  delta = 1/12 on a {n}^2 midpoint grid of the carpet distance: V1 = {f.sum() * cell / 2:.9f}, A1 = {f[v <= float(d)].sum() * cell:.9f}"
    )
    print(f"  wall {time.time() - t0:.1f} s")


def main():
    verbs = {
        "oracle": verb_oracle,
        "montecarlo": verb_montecarlo,
        "seeded": verb_seeded,
        "recompute": verb_recompute,
    }
    want = sys.argv[1:] or list(verbs)
    for v in want:
        verbs[v]()


if __name__ == "__main__":
    main()
