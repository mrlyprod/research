import argparse
import cmath
import json
import math
import os
import time

LN2 = math.log(2.0)
PERIOD = 2.0 * math.pi / LN2
ROUNDING = 1e-14
HERE = os.path.dirname(os.path.abspath(__file__))
CONTROL = os.path.join(HERE, "control.json")

# ENGINE


def roots(poly):
    n = len(poly) - 1
    while n > 0 and abs(poly[n]) < 1e-14:
        n -= 1
    c = [poly[i] / poly[n] for i in range(n + 1)]
    z = [(0.4 + 0.9j) ** i for i in range(n)]
    for _ in range(400):
        move = 0.0
        for i in range(n):
            num = sum(c[j] * z[i] ** j for j in range(n + 1))
            den = 1.0 + 0j
            for j in range(n):
                if j != i:
                    den *= z[i] - z[j]
            step = num / den
            z[i] -= step
            move = max(move, abs(step))
        if move < 1e-15:
            break
    return z


class Rule:
    def __init__(self, width, code, peel):
        self.width, self.code, self.peel = width, code, peel
        self.size = 1 << (width - 1)
        mask = self.size - 1
        self.T = [[0.0] * self.size for _ in range(self.size)]
        self.G1 = [[0.0] * self.size for _ in range(self.size)]
        for u in range(self.size):
            for a in (0, 1):
                w = (u << 1) | a
                if (code >> w) & 1:
                    self.T[w & mask][u] += 1.0
                    self.G1[w & mask][u] += float(a)
        self.faddeev()
        self.degree = max(i for i in range(self.size + 1) if abs(self.det[i]) > 0.5)
        raw = (
            [1.0 / x for x in roots(self.det[: self.degree + 1])] if self.degree else []
        )
        self.eigen = [
            complex(z.real, 0.0) if abs(z.imag) < 1e-9 * abs(z) else z for z in raw
        ]
        self.rho = max((abs(z) for z in self.eigen), default=0.0)
        self.alpha = math.log(self.rho) / LN2 if self.rho else float("-inf")
        self.low = [n for n in range(1, 1 << (peel - 1)) if self.accepts(n)]
        self.mid = [[] for _ in range(self.size)]
        for n in range(1 << (peel - 1), 1 << peel):
            if self.accepts(n):
                self.mid[n & mask].append(n)
        self.loglow = [math.log(n) for n in self.low]
        self.logmid = [[math.log(n) for n in row] for row in self.mid]
        self.counts = [0] * 64
        for j in range(1, 64):
            self.counts[j] = self.words(j)

    def accepts(self, n):
        bits = bin(n)[2:]
        for i in range(len(bits) - self.width + 1):
            if not (self.code >> int(bits[i : i + self.width], 2)) & 1:
                return False
        return True

    def words(self, j):
        if j < self.width:
            return 1 << (j - 1)
        live = [1.0 if u >> (self.width - 2) else 0.0 for u in range(self.size)]
        for _ in range(j - self.width + 1):
            live = [
                sum(self.T[v][u] * live[u] for u in range(self.size))
                for v in range(self.size)
            ]
        return sum(live)

    def faddeev(self):
        n = self.size
        eye = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        m = [row[:] for row in eye]
        self.adj, c = [m], [1.0]
        for k in range(1, n + 1):
            tm = [
                [sum(self.T[i][t] * m[t][j] for t in range(n)) for j in range(n)]
                for i in range(n)
            ]
            ck = -sum(tm[i][i] for i in range(n)) / k
            c.append(ck)
            m = [
                [tm[i][j] + (ck if i == j else 0.0) for j in range(n)] for i in range(n)
            ]
            if k < n:
                self.adj.append(m)
        self.charpoly = list(reversed(c))
        self.det = c[:]

    def deter(self, x):
        return sum(self.det[i] * x**i for i in range(self.size + 1))

    def apply_adj(self, x, v):
        out = [0j] * self.size
        p = 1.0 + 0j
        for k in range(self.size):
            for i in range(self.size):
                out[i] += p * sum(self.adj[k][i][j] * v[j] for j in range(self.size))
            p *= x
        return out

    def adj_abs(self, x):
        a = abs(x)
        return [
            [
                sum(abs(self.adj[k][i][j]) * a**k for k in range(self.size))
                for j in range(self.size)
            ]
            for i in range(self.size)
        ]

    def tail(self, sigma):
        if sigma <= self.alpha + 1e-9:
            return float("inf")
        total = 0.0
        for j in range(self.peel, 64):
            total += self.counts[j] * math.exp(-(j - 1) * sigma * LN2)
        return total

    def combs(self, m):
        out = []
        for z in self.eigen:
            re = math.log(abs(z)) / LN2 - m
            out.append((re, cmath.phase(z) / LN2))
        return out


def two_pow(w):
    return cmath.exp(-w * LN2)


def ladder(rule, s, deep, shift, cut):
    n = rule.size
    levels = max(1, math.ceil(shift - s.real))
    rows = levels + cut + 2
    g = [[0j] * n for _ in range(rows)]
    bound = [[0.0] * n for _ in range(rows)]
    for j in range(levels, rows):
        t = rule.tail(s.real + j)
        bound[j] = [t] * n
    num, bnum = None, None
    for j in range(levels - 1, -1, -1):
        w = s + j
        sig, aw = w.real, abs(w)
        acc = [sum(cmath.exp(-w * L) for L in rule.logmid[u]) for u in range(n)]
        bacc = [0.0] * n
        c, cb = 1.0 + 0j, 1.0
        for l in range(1, cut + 1):
            c = c * (-w - (l - 1)) / l
            cb = cb * (aw + l - 1) / l
            weight = c * two_pow(w + l)
            wb = cb * math.exp(-(sig + l) * LN2)
            for i in range(n):
                acc[i] += weight * sum(rule.G1[i][t] * g[j + l][t] for t in range(n))
                bacc[i] += wb * sum(rule.G1[i][t] * bound[j + l][t] for t in range(n))
        ratio = ((aw + cut + 1) / (cut + 2)) * 2.0 ** (-rule.peel)
        if ratio >= 0.5:
            raise ValueError(f"the l-cut does not close at s = {s}, ratio {ratio}")
        rest = cb * (aw + cut) / (cut + 1) * math.exp(-(sig + cut + 1) * LN2)
        rest *= rule.tail(sig + cut + 1) / (1.0 - ratio)
        for i in range(n):
            bacc[i] += rest * sum(rule.G1[i][t] for t in range(n))
        if j == 0:
            num, bnum = acc, bacc
            if not deep:
                break
        x = two_pow(w)
        det = rule.deter(x)
        g[j] = [z / det for z in rule.apply_adj(x, acc)]
        ia = rule.adj_abs(x)
        bound[j] = [
            sum(ia[i][t] * bacc[t] for t in range(n)) / abs(det) for i in range(n)
        ]
    return g[0], bound[0], num, bnum


class Engine:
    def __init__(self, rule, shift, cut, add=()):
        self.rule, self.shift, self.cut = rule, shift, cut
        self.add = tuple(sorted(add))
        self.logadd = [math.log(n) for n in self.add]

    def poly_low(self, s):
        return sum(cmath.exp(-s * L) for L in self.rule.loglow) + sum(
            cmath.exp(-s * L) for L in self.logadd
        )

    def poly_abs(self, sigma):
        return sum(math.exp(-sigma * L) for L in self.rule.loglow) + sum(
            math.exp(-sigma * L) for L in self.logadd
        )

    def cofactor(self, s):
        _, _, num, bnum = ladder(self.rule, s, False, self.shift, self.cut)
        x = two_pow(s)
        det = self.rule.deter(x)
        top = sum(self.rule.apply_adj(x, num))
        ia = self.rule.adj_abs(x)
        n = self.rule.size
        lift = [sum(ia[i][j] for i in range(n)) for j in range(n)]
        bound = sum(lift[j] * bnum[j] for j in range(n))
        scale = abs(det) * self.poly_abs(s.real) + sum(
            lift[j] * abs(num[j]) for j in range(n)
        )
        return det * self.poly_low(s) + top, bound + ROUNDING * scale

    def zeta(self, s):
        g, b, _, _ = ladder(self.rule, s, True, self.shift, self.cut)
        value = self.poly_low(s) + sum(g)
        return value, sum(b) + ROUNDING * (
            self.poly_abs(s.real) + sum(abs(z) for z in g)
        )

    def residue(self, s0):
        value, bound = self.cofactor(s0)
        x = two_pow(s0)
        slope = sum(
            i * self.rule.det[i] * x ** (i - 1) for i in range(1, self.rule.size + 1)
        )
        den = -x * LN2 * slope
        return value / den, bound / abs(den)


# BOX


def geometry(rule, left, right, height):
    combs, poles = [], []
    for re, off in rule.combs(0):
        pts = [
            complex(re, off + PERIOD * j) for j in range(-2, int(height / PERIOD) + 3)
        ]
        combs.append((re, [p for p in pts if 0.0 < p.imag < height]))
        m = 1
        while re - m > left:
            poles += [complex(re - m, p.imag) for p in pts if 0.0 < p.imag < height]
            m += 1
    combs.sort(key=lambda c: -c[0])
    lines = sorted(
        {round(c[0], 12) for c in combs} | {round(p.real, 12) for p in poles}
    )
    lines = [v for v in lines if left < v < right]
    edges = (
        [left]
        + [0.5 * (lines[i] + lines[i + 1]) for i in range(len(lines) - 1)]
        + [right]
    )
    wide = []
    for i in range(len(edges) - 1):
        span = edges[i + 1] - edges[i]
        parts = max(1, math.ceil(span / 0.75))
        wide += [edges[i] + span * t / parts for t in range(parts)]
    edges = wide + [right]
    heights = sorted(
        {round(p.imag, 9) for c in combs for p in c[1]}
        | {round(p.imag, 9) for p in poles}
    )
    rows = [0.02] + [
        0.5 * (heights[i] + heights[i + 1]) for i in range(len(heights) - 1)
    ]
    rows = [r for r in rows if r < height] + [height]
    return combs, poles, edges, rows


def show(z, w=12):
    return f"{z.real:.{w}f}{z.imag:+.{w}f}i"


class Census:
    def __init__(self, engine, combs, poles, edges, rows, seed):
        self.engine, self.combs, self.poles = engine, combs, poles
        self.edges, self.rows, self.seed = edges, rows, seed
        self.cache, self.worst, self.bound, self.calls = {}, 0.0, 0.0, 0
        self.probes = 0

    def value(self, s):
        key = (round(s.real, 12), round(s.imag, 12))
        hit = self.cache.get(key)
        if hit is None:
            hit, b = self.engine.cofactor(s)
            self.cache[key] = hit
            self.bound = max(self.bound, b)
            self.calls += 1
        return hit

    def segment(self, a, b, fa, fb, depth):
        step = cmath.phase(fb / fa)
        if abs(step) <= 1.0 or depth >= 18:
            self.worst = max(self.worst, abs(step))
            return step
        m = 0.5 * (a + b)
        fm = self.value(m)
        return self.segment(a, m, fa, fm, depth + 1) + self.segment(
            m, b, fm, fb, depth + 1
        )

    def side(self, a, b):
        steps = max(4, math.ceil(abs(b - a) / self.seed))
        total, prev, fprev = 0.0, a, self.value(a)
        for i in range(1, steps + 1):
            nxt = a + (b - a) * i / steps
            fnext = self.value(nxt)
            total += self.segment(prev, nxt, fprev, fnext, 0)
            prev, fprev = nxt, fnext
        return total

    def winding(self, box):
        x0, x1, y0, y1 = box
        corner = [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)]
        return sum(self.side(corner[i], corner[(i + 1) % 4]) for i in range(4)) / (
            2.0 * math.pi
        )

    def count(self, box):
        x0, x1, y0, y1 = box
        turn = self.winding(box)
        held = sum(1 for p in self.poles if x0 < p.real < x1 and y0 < p.imag < y1)
        return round(turn + held), turn, held

    def split(self, box):
        x0, x1, y0, y1 = box
        marks = self.poles + [p for c in self.combs for p in c[1]]
        if x1 - x0 >= y1 - y0:
            cut = 0.5 * (x0 + x1)
            while any(abs(p.real - cut) < 1e-4 for p in marks):
                cut += 1.7e-3
            return (x0, cut, y0, y1), (cut, x1, y0, y1)
        cut = 0.5 * (y0 + y1)
        while any(abs(p.imag - cut) < 1e-4 for p in marks):
            cut += 1.7e-3
        return (x0, x1, y0, cut), (x0, x1, cut, y1)

    def hunt(self, box, n, depth):
        if n <= 0:
            return []
        x0, x1, y0, y1 = box
        if (x1 - x0 < 0.01 and y1 - y0 < 0.01) or depth > 44:
            return [complex(0.5 * (x0 + x1), 0.5 * (y0 + y1))] * n
        left, right = self.split(box)
        nl = self.count(left)[0]
        return self.hunt(left, nl, depth + 1) + self.hunt(right, n - nl, depth + 1)

    def polish(self, z):
        a, b = z, z + 1e-4
        fa, fb = self.value(a), self.value(b)
        for _ in range(60):
            if abs(fb - fa) < 1e-300:
                break
            c = b - fb * (b - a) / (fb - fa)
            if abs(c - b) < 1e-13:
                b = c
                break
            a, fa, b, fb = b, fb, c, self.value(c)
        return b, abs(self.value(b))

    def circle(self, s0, eps, n=48):
        total, bound = 0j, 0.0
        for i in range(n):
            u = cmath.exp(2j * math.pi * i / n)
            v, b = self.engine.cofactor(s0 + eps * u)
            total += v * u
            bound = max(bound, b)
            self.probes += 1
        return eps * total / n, eps * bound

    def arc(self, c, rad, a, b, fa, fb, depth):
        step = cmath.phase(fb / fa)
        if abs(step) <= 1.0 or depth >= 18:
            self.worst = max(self.worst, abs(step))
            return step
        m = 0.5 * (a + b)
        fm = self.value(c + rad * cmath.exp(1j * m))
        return self.arc(c, rad, a, m, fa, fm, depth + 1) + self.arc(
            c, rad, m, b, fm, fb, depth + 1
        )

    def ring(self, c, rad, n=40):
        f0 = self.value(c + rad)
        total, prev, fprev = 0.0, 0.0, f0
        for i in range(1, n + 1):
            th = 2.0 * math.pi * i / n
            fnext = f0 if i == n else self.value(c + rad * cmath.exp(1j * th))
            total += self.arc(c, rad, prev, th, fprev, fnext, 0)
            prev, fprev = th, fnext
        return total / (2.0 * math.pi)

    def disc(self, c, rad):
        held = sum(1 for p in self.poles if abs(p - c) < rad)
        return round(self.ring(c, rad)) + held

    def guarded(self, c, rad, guard):
        near = min((abs(abs(p - c) - rad) for p in self.poles), default=float("inf"))
        return self.disc(c, rad - guard), self.disc(c, rad + guard), near

    def regular(self, s0, eps, n=64):
        total = 0j
        for i in range(n):
            u = eps * cmath.exp(2j * math.pi * i / n)
            total += self.value(s0 + u) / self.engine.rule.deter(two_pow(s0 + u))
        return total / n


# VERBS


def census(args):
    start = time.perf_counter()
    rule = Rule(args.width, args.code, args.peel)
    engine = Engine(rule, args.shift, args.cut)
    combs, poles, edges, rows = geometry(rule, args.left, args.right, args.height)
    height = rows[-1]
    run = Census(engine, combs, poles, edges, rows, args.seed)
    print(
        f"code {args.code}, k = {args.width}, D = 1, base 2, states {rule.size},"
        f" peel {args.peel}, shift {args.shift}, cut {args.cut}"
    )
    print(f"det(I - x T) coefficients {[round(c, 12) for c in rule.det]}")
    print(
        f"eigenvalues {' '.join(show(z, 12) for z in rule.eigen)}, alpha {rule.alpha:.12f}"
    )
    for i, (re, pts) in enumerate(combs):
        print(
            f"comb {i} Re s = {re:.12f} teeth {len(pts)} first"
            f" {pts[0].imag if pts else float('nan'):.12f} spacing {PERIOD:.12f}"
        )
    print(
        f"box Re s in [{args.left}, {args.right}], Im s in [{rows[0]}, {height:.12f}]"
    )
    print(
        f"cofactor poles in the box {len(poles)}"
        f" on Re s = {sorted({round(p.real, 9) for p in poles})}"
    )
    ranked = sorted(poles, key=lambda q: q.imag)
    for p in ranked:
        r1, b1 = run.circle(p, args.ring)
        r2, _ = run.circle(p, 0.4 * args.ring)
        lift, _ = engine.cofactor(p + 1e-5)
        print(
            f"  pole Im {p.imag:10.6f} residue {show(r1, 12)} bound {b1:.2e},"
            f" at radius {0.4 * args.ring:.3f} gap {abs(r1 - r2):.2e},"
            f" simple to {abs(1e-5 * lift - r1) / abs(r1):.2e}"
        )
    if len(ranked) > 1:
        blank = complex(ranked[0].real, 0.5 * (ranked[0].imag + ranked[1].imag))
        rb, _ = run.circle(blank, args.ring)
        print(
            f"  blank point on the same line Im {blank.imag:10.6f}"
            f" circle mean {abs(rb):.2e}"
        )
    cells, total = [], 0
    for i in range(len(edges) - 1):
        for j in range(len(rows) - 1):
            box = (edges[i], edges[i + 1], rows[j], rows[j + 1])
            n, turn, held = run.count(box)
            cells.append((box, n))
            total += n
            if n or held:
                print(
                    f"cell Re [{box[0]:7.4f},{box[1]:7.4f}] Im [{box[2]:9.5f},{box[3]:9.5f}]"
                    f" winding {turn:+9.5f} poles {held} zeros {n}"
                )
    print(f"zeros in the box {total}, cells {len(cells)}")
    zeros = []
    for box, n in cells:
        for z in run.hunt(box, n, 0):
            zeros.append(run.polish(z))
    zeros.sort(key=lambda p: p[0].imag)
    print(
        f"located {len(zeros)} of {total},"
        f" largest residual {max((r for _, r in zeros), default=0.0):.3e}"
    )
    for z, _ in zeros:
        gaps = " ".join(
            f"c{i} {min((abs(z - t) for t in pts), default=float('inf')):.9f}"
            for i, (_, pts) in enumerate(combs)
        )
        print(f"  zero {show(z)} {gaps}")
    on_edge = [
        z
        for z, _ in zeros
        if min(
            abs(z.real - args.left),
            abs(z.real - args.right),
            abs(z.imag - rows[0]),
            abs(z.imag - height),
        )
        < 0.02
    ]
    on_pole = [
        p
        for p in poles
        if min(
            abs(p.real - args.left),
            abs(p.real - args.right),
            abs(p.imag - rows[0]),
            abs(p.imag - height),
        )
        < 0.02
    ]
    print(
        f"zeros within 0.02 of a box edge {len(on_edge)} {[show(z, 9) for z in on_edge]}"
    )
    print(f"poles within 0.02 of a box edge {len(on_pole)}")
    inner = [
        min(
            [abs(z.real - e) for e in edges[1:-1]]
            + [abs(z.imag - r) for r in rows[1:-1]]
            + [float("inf")]
        )
        for z, _ in zeros
    ]
    tight = min(range(len(zeros)), key=lambda i: inner[i])
    print(
        f"tightest clearance to an internal cell edge {inner[tight]:.6f}"
        f" at {show(zeros[tight][0])}"
    )
    held = set()
    for i, (re, pts) in enumerate(combs):
        near = [min(abs(z - t) for z, _ in zeros) for t in pts]
        for t in pts:
            for z, _ in zeros:
                if abs(z - t) < args.rho:
                    held.add((round(z.real, 9), round(z.imag, 9)))
        print(
            f"comb {i} Re s = {re:.12f} teeth {len(pts)}"
            f" occupied at radius {args.rho} {sum(1 for d in near if d < args.rho)}"
        )
        print(f"  distances {' '.join(f'{d:.9f}' for d in near)}")
        for t in pts:
            r, rb = engine.residue(t)
            reg = run.regular(t, args.eps)
            u1 = -r / reg
            best = min(zeros, key=lambda p: abs(p[0] - t))
            disc = sum(1 for z, _ in zeros if abs(z - t) < args.rho)
            print(
                f"  tooth Im {t.imag:10.6f} r {show(r, 9)} bound {rb:.2e}"
                f" R {show(reg, 9)} u1 {show(u1, 9)} abs {abs(u1):.9f}"
                f" zero at {abs(best[0] - t):.9f}"
                f" modulus miss {abs(abs(best[0] - t) - abs(u1)):.9f}"
                f" vector miss {abs(best[0] - t - u1):.9f}"
                f" zeros in the disc {disc}"
            )
    family = [z for z, _ in zeros if (round(z.real, 9), round(z.imag, 9)) not in held]
    print(f"tooth zeros {len(zeros) - len(family)}, second family {len(family)}")
    if family:
        print(
            f"second family Re s in [{min(z.real for z in family):.12f},"
            f" {max(z.real for z in family):.12f}],"
            f" mean {sum(z.real for z in family) / len(family):.12f}"
        )
    for i, (re, _) in enumerate(combs):
        hug = [z for z, _ in zeros if abs(z.real - re) < args.hug]
        print(
            f"zeros within {args.hug} of comb {i}'s line: {len(hug)}"
            f" {[f'{show(z, 9)} at {abs(z.real - re):.9f}' for z in hug]}"
        )
    axis = 0.5 * (combs[0][0] + combs[-1][0])
    mates = sum(
        1
        for i, (z, _) in enumerate(zeros)
        for j, (w, _) in enumerate(zeros)
        if j != i
        and abs(w.imag - z.imag) < 0.05
        and abs(w.real + z.real - 2.0 * axis) < 0.05
    )
    selfish = [z for z, _ in zeros if abs(2.0 * (z.real - axis)) < 0.05]
    near = min(abs(2.0 * (z.real - axis)) for z, _ in zeros)
    print(
        f"reflection partners other than the zero itself about Re s = {axis:.12f},"
        f" the axis of the comb lines: {mates} of {len(zeros)}"
    )
    print(
        f"zeros whose own reflection lands within 0.05 of them {len(selfish)},"
        f" nearest self-reflection {near:.12f} {[show(z, 12) for z in selfish]}"
    )
    guard, gbound = engine.zeta(complex(args.right, 0.0))
    print(
        f"zeta_W({args.right}) = {guard.real:.12f} with a_min = 1, bound {gbound:.3e}"
    )
    print(
        f"cuts: Im s in [0, {rows[0]}] uncounted, Im s > {height:.12f} uncounted,"
        f" Re s < {args.left} uncounted"
    )
    print(
        f"largest surviving phase step {run.worst:.6f} rad, largest bound {run.bound:.3e},"
        f" evaluations {run.calls} on the contour and the hunt,"
        f" {run.probes} off it on the residue circles, runtime {time.perf_counter() - start:.2f} s"
    )


# CLASSES


def window_maps(k):
    tables = set()
    for flip in (0, 1):
        for rev in (False, True):
            t = []
            for w in range(1 << k):
                ds = [((w >> (k - 1 - j)) & 1) ^ flip for j in range(k)]
                if rev:
                    ds = ds[::-1]
                v = 0
                for c in ds:
                    v = (v << 1) | c
                t.append(v)
            tables.add(tuple(t))
    return sorted(tables)


def classes(k):
    nbits = 1 << k
    maps = []
    for t in window_maps(k):
        arr = [0] * (1 << nbits)
        for code in range(1, 1 << nbits):
            low = code & -code
            arr[code] = arr[code ^ low] | (1 << t[low.bit_length() - 1])
        maps.append(arr)
    seen = bytearray(1 << nbits)
    out = []
    for c in range(1 << nbits):
        if seen[c]:
            continue
        orb = {m[c] for m in maps}
        for x in orb:
            seen[x] = 1
        out.append((c, len(orb)))
    return out


def spectrum(rule):
    order = sorted(rule.eigen, key=lambda z: -abs(z))
    mods = []
    for z in order:
        if not any(abs(abs(z) - m) < 1e-9 for m in mods):
            mods.append(abs(z))
    gap = min(
        (
            abs(order[i] - order[j])
            for i in range(len(order))
            for j in range(i + 1, len(order))
        ),
        default=float("inf"),
    )
    return order, mods, gap


# TEETH


def flank(rule, walls):
    left, _right, low, high = walls
    out = []
    for re, off in rule.combs(0):
        pts = [off + PERIOD * j for j in range(-2, int(high / PERIOD) + 3)]
        m = 1
        while re - m > left - 1.0:
            out += [complex(re - m, y) for y in pts if low - 1.0 < y < high + 1.0]
            m += 1
    return out


def wall(p, walls):
    left, right, low, high = walls
    dx = max(left - p.real, p.real - right, 0.0)
    dy = max(low - p.imag, p.imag - high, 0.0)
    if dx == 0.0 and dy == 0.0:
        return min(p.real - left, right - p.real, p.imag - low, high - p.imag)
    return math.hypot(dx, dy)


def sweep(engine, rule, args):
    combs, poles, edges, rows = geometry(rule, args.left, args.right, args.height)
    run = Census(engine, combs, poles, edges, rows, args.seed)
    total = 0
    cells = []
    for i in range(len(edges) - 1):
        for j in range(len(rows) - 1):
            box = (edges[i], edges[i + 1], rows[j], rows[j + 1])
            n = run.count(box)[0]
            cells.append((box, n))
            total += n
    zeros = []
    for box, n in cells:
        for z in run.hunt(box, n, 0):
            zeros.append(run.polish(z))
    zeros.sort(key=lambda p: p[0].imag)
    return combs, poles, edges, rows, cells, total, zeros, run


def teeth(args):
    start = time.perf_counter()
    reps = classes(args.width)
    print(
        f"D = 1, base 2, width {args.width}, classes under G_(1,k) {len(reps)},"
        f" peel {args.peel}, shift {args.shift}, cut {args.cut}"
    )
    print(
        f"box Re s in [{args.left}, {args.right}], Im s in [0.02, {args.height}],"
        f" contour seed {args.seed}, occupancy radius {args.rho}"
    )
    picked, skipped = [], []
    for code, size in reps:
        rule = Rule(args.width, code, args.peel)
        if not rule.eigen:
            continue
        order, mods, gap = spectrum(rule)
        if len(mods) != 2:
            continue
        if gap < 1e-9:
            skipped.append(code)
            continue
        picked.append((code, size, rule, order, mods))
    picked.sort(key=lambda r: -r[4][1] / r[4][0])
    print(
        f"two-line classes {len(picked)} of {len(reps)},"
        f" dropped for a repeated eigenvalue {len(skipped)} {skipped}"
    )
    table, law, guard = [], [], []
    for code, size, rule, order, mods in picked:
        engine = Engine(rule, args.shift, args.cut)
        combs, poles, _edges, rows, cells, total, zeros, run = sweep(engine, rule, args)
        height = rows[-1]
        print(
            f"code {code}, k = {args.width}, class size {size}, states {rule.size},"
            f" det {[round(c, 9) for c in rule.det]}"
        )
        print(
            f"  eigenvalues {' '.join(show(z, 9) for z in order)}, rho {mods[0]:.12f},"
            f" abs(l2)/rho {mods[1] / mods[0]:.12f}"
        )
        print(
            f"  cells {len(cells)} poles {len(poles)} zeros {total} located {len(zeros)}"
            f" largest residual {max((r for _, r in zeros), default=0.0):.3e}"
            f" phase step {run.worst:.6f} bound {run.bound:.3e} calls {run.calls}"
        )
        held, row = set(), []
        for j, m in enumerate(mods):
            pts = sorted(
                [p for re, ps in combs if abs(2.0**re - m) < 1e-9 for p in ps],
                key=lambda p: p.imag,
            )
            args_of = sorted(
                {round(cmath.phase(z), 9) for z in order if abs(abs(z) - m) < 1e-9}
            )
            near = [
                min((abs(z - t) for z, _ in zeros), default=float("inf")) for t in pts
            ]
            for t, d in zip(pts, near):
                if d < args.rho:
                    for z, _ in zeros:
                        if abs(z - t) < args.rho:
                            held.add((round(z.real, 9), round(z.imag, 9)))
            hit = sum(1 for d in near if d < args.rho)
            print(
                f"  line {j} Re s = {math.log(m) / LN2:.12f} arg"
                f" {' '.join(f'{a:+.9f}' for a in args_of)} teeth {len(pts)}"
                f" occupied {hit} least distance {min(near, default=float('nan')):.9f}"
            )
            print(f"    distances {' '.join(f'{d:.9f}' for d in near)}")
            for t, d in zip(pts, near):
                r, rb = engine.residue(t)
                reg = run.regular(t, args.eps)
                u1 = -r / reg
                best = min(zeros, key=lambda p: abs(p[0] - t), default=None)
                vec = abs(best[0] - t - u1) if best else float("inf")
                disc = sum(1 for z, _ in zeros if abs(z - t) < args.rho)
                law.append((code, j, abs(u1), d, d < args.rho, vec, disc))
                print(
                    f"    tooth Im {t.imag:10.6f} r {show(r, 9)} bound {rb:.2e}"
                    f" R {show(reg, 9)} u1 {show(u1, 9)} abs {abs(u1):.9f}"
                    f" zero at {d:.9f} modulus miss {abs(d - abs(u1)):.9f}"
                    f" vector miss {vec:.9f} zeros in the disc {disc}"
                )
            row.append((len(pts), hit, args_of))
        off = [z for z, _ in zeros if (round(z.real, 9), round(z.imag, 9)) not in held]
        edge = [
            z
            for z, _ in zeros
            if min(
                abs(z.real - args.left),
                abs(z.real - args.right),
                abs(z.imag - rows[0]),
                abs(z.imag - height),
            )
            < 0.02
        ]
        walls = (args.left, args.right, rows[0], height)
        reach = [wall(p, walls) for p in flank(rule, walls)]
        on_pole = [d for d in reach if d < 0.02]
        clear = min(reach + [float("inf")])
        guard.append((code, len(edge), len(on_pole), clear))
        print(
            f"  zeros off every tooth {len(off)} of {total},"
            f" within 0.02 of a contour {len(edge)} {[show(z, 9) for z in edge]},"
            f" poles within 0.02 of a contour {len(on_pole)},"
            f" least pole clearance {clear:.12f}"
        )
        reach = mods[1] and math.log(mods[1]) / LN2 - args.rho
        print(
            f"  second line disc reaches Re s = {reach:.12f},"
            f" inside the box {reach > args.left}"
        )
        table.append((code, size, mods[1] / mods[0], row, total, len(off)))
    print("table rule class ratio arg2 teeth1 occupied1 teeth2 occupied2 zeros off")
    for code, size, ratio, row, total, off in table:
        a2 = " ".join(f"{a:+.6f}" for a in row[1][2])
        print(
            f"  {code:4d} {size:3d} {ratio:.9f} [{a2}]"
            f" {row[0][0]:3d} {row[0][1]:3d} {row[1][0]:3d} {row[1][1]:3d} {total:3d} {off:3d}"
        )
    full = [r for r in table if r[3][1][0] and r[3][1][1] == r[3][1][0]]
    empty = [r for r in table if r[3][1][0] and r[3][1][1] == 0]
    print(
        f"second line full {len(full)} {[r[0] for r in full]},"
        f" empty {len(empty)} {[r[0] for r in empty]},"
        f" partial {len(table) - len(full) - len(empty)}"
    )
    inside = [t for t in law if t[2] < args.eps]
    outside = [t for t in law if t[2] >= args.rho]
    print(
        f"teeth {len(law)}, occupied {sum(1 for t in law if t[4])};"
        f" first-order prediction abs(u1) below {args.eps}, inside the disc that builds R:"
        f" {len(inside)} teeth, occupied {sum(1 for t in inside if t[4])}"
    )
    print(
        f"prediction abs(u1) at or above {args.rho}: {len(outside)} teeth,"
        f" occupied {sum(1 for t in outside if t[4])}"
    )
    band = [t for t in law if args.eps <= t[2] < args.rho]
    print(
        f"prediction in [{args.eps}, {args.rho}): {len(band)} teeth,"
        f" occupied {sum(1 for t in band if t[4])} {[(t[0], round(t[2], 6), round(t[3], 6)) for t in band]}"
    )
    worst = max(law, key=lambda t: abs(t[3] - t[2]), default=None)
    print(
        f"largest modulus miss of the first-order law {abs(worst[3] - worst[2]):.9f}"
        f" at code {worst[0]} line {worst[1]}, abs(u1) {worst[2]:.9f} zero at {worst[3]:.9f}"
    )
    vworst = max(law, key=lambda t: t[5], default=None)
    print(
        f"largest vector miss of the first-order law {vworst[5]:.9f}"
        f" at code {vworst[0]} line {vworst[1]}, abs(u1) {vworst[2]:.9f} zero at {vworst[3]:.9f}"
    )
    doubles = [t for t in law if t[6] > 1]
    print(
        f"teeth holding more than one zero inside {args.rho} {len(doubles)}"
        f" {[(t[0], t[1], t[6]) for t in doubles]}, zeros inside a tooth disc"
        f" {sum(t[6] for t in law)} against occupied teeth"
        f" {sum(1 for t in law if t[4])}"
    )
    tight = min(guard, key=lambda g: g[3], default=None)
    print(
        f"contour guard: zeros within 0.02 of a contour {sum(g[1] for g in guard)},"
        f" poles within 0.02 of a contour {sum(g[2] for g in guard)},"
        f" least pole clearance {tight[3]:.12f} at code {tight[0]}"
    )
    print(f"runtime {time.perf_counter() - start:.2f} s")


# BRIDGE


def bridge(args):
    start = time.perf_counter()
    a, b = Rule(2, 7, args.peel), Rule(3, 55, args.peel)
    top = 1 << args.span
    extra = [n for n in range(1, top) if b.accepts(n) and not a.accepts(n)]
    lost = [n for n in range(1, top) if a.accepts(n) and not b.accepts(n)]
    print(
        f"code 7 at k = 2 against code 55 at k = 3 over 1 .. {top - 1}:"
        f" in 55 and not 7 {extra}, in 7 and not 55 {lost}"
    )
    print(
        f"det(I - x T) code 7 {[round(c, 12) for c in a.det]},"
        f" code 55 {[round(c, 12) for c in b.det]}, states {a.size} against {b.size}"
    )
    ea, eb = Engine(a, args.shift, args.cut), Engine(b, args.shift, args.cut)
    worst, over = 0.0, 0
    alpha = math.log((1.0 + math.sqrt(5.0)) / 2.0) / LN2
    for s in [
        complex(2.0, 0.0),
        complex(0.8, 0.0),
        complex(-0.95, 20.0),
        complex(alpha, PERIOD),
        complex(-alpha, 0.5 * PERIOD),
        complex(-alpha, 1.5 * PERIOD),
        complex(-0.442302243578, 4.612546440182),
    ]:
        va, ba = ea.cofactor(s)
        vb, bb = eb.cofactor(s)
        want = b.deter(two_pow(s)) * cmath.exp(-s * math.log(3.0))
        gap = abs(vb - va - want)
        worst = max(worst, gap)
        over += gap > ba + bb
        print(
            f"  s {show(s, 6)} Z7 {show(va, 12)} Z55 {show(vb, 12)}"
            f" gap {gap:.3e} bound {ba + bb:.3e} {'in' if gap <= ba + bb else 'OUT'}"
        )
    print(
        f"largest gap {worst:.3e}, outside its bound {over} of 7,"
        f" runtime {time.perf_counter() - start:.2f} s"
    )


def control(args):
    rule = Rule(2, 7, args.peel)
    engine = Engine(rule, args.shift, args.cut)
    with open(CONTROL) as f:
        data = json.load(f)
    print(
        f"control {data['source']}, dps {data['dps']}, peel {data['peel']},"
        f" shift {data['shift']}, cut {data['cut']}"
    )
    worst, over = 0.0, 0
    for row in data["rows"]:
        s = complex(float(row["re"]), float(row["im"]))
        got, bound = (
            engine.cofactor(s)
            if row["kind"] == "Z"
            else engine.zeta(s)
            if row["kind"] == "zeta"
            else engine.residue(s)
        )
        want = complex(float(row["value"][0]), float(row["value"][1]))
        gap = abs(got - want)
        worst = max(worst, gap)
        over += gap > bound
        print(
            f"  {row['kind']:8s} {show(s, 6)} value {show(got, 12)} gap {gap:.3e}"
            f" bound {bound:.3e} {'in' if gap <= bound else 'OUT'}"
        )
    print(f"largest gap {worst:.3e}, outside its bound {over} of {len(data['rows'])}")


# DIAL


def lines_of(combs):
    out = []
    for re, pts in combs:
        key = round(re, 9)
        for i, (r, ps) in enumerate(out):
            if abs(r - key) < 1e-9:
                out[i] = (r, ps + list(pts))
                break
        else:
            out.append((key, list(pts)))
    out.sort(key=lambda c: -c[0])
    return [(r, sorted(ps, key=lambda p: p.imag)) for r, ps in out]


def exact_min(reg, t, cand, cap):
    import numpy as np

    half = len(cand) // 2
    left, right = cand[:half], cand[half:]

    def table(part):
        sums = np.array([0.0 + 0.0j])
        size = np.array([0])
        for n in part:
            v = cmath.exp(-t * math.log(n))
            sums = np.concatenate([sums, sums + v])
            size = np.concatenate([size, size + 1])
        return sums, size

    ls, lk = table(left)
    rs, rk = table(right)
    best, pick = float("inf"), None
    for j in range(ls.size):
        room = cap - int(lk[j])
        if room < 0:
            continue
        ok = np.flatnonzero(rk <= room)
        d = np.abs(reg + ls[j] + rs[ok])
        w = int(np.argmin(d))
        if float(d[w]) < best:
            best = float(d[w])
            idx = int(ok[w])
            pick = [left[e] for e in range(len(left)) if j >> e & 1] + [
                right[e] for e in range(len(right)) if idx >> e & 1
            ]
    return best, sorted(pick)


def dial(args):
    start = time.perf_counter()
    rule = Rule(args.width, args.code, args.peel)
    base = Engine(rule, args.shift, args.cut)
    combs, poles, edges, rows = geometry(rule, args.left, args.right, args.height)
    lines = lines_of(combs)
    height = rows[-1]
    cand = [n for n in range(2, args.top + 1) if not rule.accepts(n)]
    print(
        f"code {args.code}, k = {args.width}, D = 1, base 2, states {rule.size},"
        f" peel {args.peel}, shift {args.shift}, cut {args.cut}"
    )
    print(f"det(I - x T) coefficients {[round(c, 12) for c in rule.det]}")
    print(
        f"eigenvalues {' '.join(show(z, 12) for z in rule.eigen)}, alpha {rule.alpha:.12f}"
    )
    print(
        f"box Im s in [{rows[0]}, {height:.12f}], Re s in [{args.left}, {args.right}],"
        f" occupancy radius {args.rho}, ring samples 40, edge guard {args.guard},"
        f" R radius {args.eps}"
    )
    for i, (re, pts) in enumerate(lines):
        print(
            f"line {i} Re s = {re:.12f} teeth {len(pts)}"
            f" Im {' '.join(f'{p.imag:.6f}' for p in pts)}"
        )
    walls = (args.left, args.right, rows[0], height)
    wide = flank(rule, walls)
    reach = min(t.real for _, pts in lines for t in pts) - args.rho - args.guard
    outer = max([p.real for p in wide if p.real <= args.left] + [float("-inf")])
    print(
        f"cofactor poles in the box {len(poles)}"
        f" on Re s = {sorted({round(p.real, 9) for p in poles})}, padded pole list"
        f" {len(wide)} on Re s = {sorted({round(p.real, 9) for p in wide})}"
    )
    print(
        f"deepest disc reach Re s = {reach:.12f}, nearest pole line left of the box"
        f" {outer:.12f}, clear {reach > outer}"
    )
    print(
        f"candidates, the integers of 2 .. {args.top} outside S_W, {len(cand)}: {cand}"
    )
    ref = Census(base, combs, wide, edges, rows, args.seed)
    print(
        "falsification, the perturbation is entire so every pole and residue must hold"
    )
    worst_res, worst_val, worst_off, worst_add = 0.0, 0.0, 0.0, 0.0
    for probe in ([cand[0]], cand[:3], [cand[-1]]):
        eng = Engine(rule, args.shift, args.cut, probe)
        run = Census(eng, combs, wide, edges, rows, args.seed)
        for p in sorted(poles, key=lambda q: q.imag):
            r0, _ = ref.circle(p, args.ring)
            r1, _ = run.circle(p, args.ring)
            worst_res = max(worst_res, abs(r1 - r0))
        gaps = []
        for _, pts in lines:
            for t in pts:
                v0, _ = base.cofactor(t)
                v1, _ = eng.cofactor(t)
                gaps.append(abs(v1 - v0))
                worst_val = max(worst_val, abs(v1 - v0))
        off, add = 0.0, 0.0
        for _, pts in lines:
            for t in pts:
                s0 = t + args.rho
                part = rule.deter(two_pow(s0)) * sum(
                    cmath.exp(-s0 * math.log(n)) for n in probe
                )
                v0, _ = base.cofactor(s0)
                v1, _ = eng.cofactor(s0)
                off = max(off, abs(v1 - v0 - part))
                add = max(add, abs(part))
        worst_off, worst_add = max(worst_off, off), max(worst_add, add)
        print(
            f"  F {probe} poles {len(poles)} largest residue gap"
            f" {worst_res:.3e}, teeth {len(gaps)} largest Z gap {max(gaps):.3e},"
            f" off-tooth identity largest miss {off:.3e} against an added part of"
            f" {add:.6f}"
        )
    print(
        f"largest residue gap over every probe {worst_res:.3e},"
        f" largest tooth value gap {worst_val:.3e}, largest off-tooth identity miss"
        f" {worst_off:.3e} against an added part of up to {worst_add:.6f}"
    )
    print(
        "the residue and tooth probes are blind to the size of an entire addition,"
        " the circle mean annihilating it and the determinant vanishing at a tooth;"
        " the off-tooth identity is the one that measures it"
    )
    rows_out = []
    for i, (re, pts) in enumerate(lines):
        for t in pts:
            r, rb = base.residue(t)
            reg = ref.regular(t, args.eps)
            u1 = -r / reg
            lo, hi, near = ref.guarded(t, args.rho, args.guard)
            rows_out.append((i, t, r, reg, u1, lo, hi, near))
            print(
                f"  tooth line {i} Im {t.imag:10.6f} r {show(r, 9)} bound {rb:.2e}"
                f" R {show(reg, 9)} u1 {show(u1, 9)} abs {abs(u1):.9f}"
                f" zeros in the disc {lo}/{hi} pole clearance {near:.6f}"
            )
    print(
        f"baseline occupied {sum(1 for q in rows_out if q[5] > 0)} of {len(rows_out)} teeth,"
        f" edge guard splits {sum(1 for q in rows_out if q[5] != q[6])}"
    )
    grid, split = {}, 0
    for n in cand:
        eng = Engine(rule, args.shift, args.cut, [n])
        run = Census(eng, combs, wide, edges, rows, args.seed)
        for i, t, r, reg, u1, lo, hi, near in rows_out:
            a, b, _ = run.guarded(t, args.rho, args.guard)
            split += a != b
            grid[(round(t.imag, 6), n)] = (a, b)
    agree, cells, flips = 0, 0, []
    perline = {}
    print("grid, one row per tooth, one column per candidate in the printed order")
    for i, t, r, reg, u1, lo, hi, near in rows_out:
        meas, pred, bad = "", "", []
        for n in cand:
            a, b = grid[(round(t.imag, 6), n)]
            p = abs(-r / (reg + cmath.exp(-t * math.log(n)))) < args.rho
            meas += "1" if a else "0"
            pred += "1" if p else "0"
            cells += 1
            hit, tot, one, und, hitb, oneb = perline.get(i, (0, 0, 0, 0, 0, 0))
            perline[i] = (
                hit + (bool(a) == p),
                tot + 1,
                one + bool(a),
                und + (not a and b > 0),
                hitb + (bool(b) == p),
                oneb + bool(b),
            )
            if bool(a) == p:
                agree += 1
            else:
                bad.append(n)
        got = [n for n, c in zip(cand, meas) if c == "1"]
        miss = [n for n, c in zip(cand, meas) if c == "0"]
        print(
            f"  tooth line {i} Im {t.imag:10.6f} base {1 if lo else 0}"
            f" measured {meas} occupied {len(got)} of {len(cand)}"
        )
        print(f"    predicted {pred} mismatch {bad}")
        seam = [
            n
            for n in cand
            if grid[(round(t.imag, 6), n)][0] != grid[(round(t.imag, 6), n)][1]
        ]
        print(
            f"    last candidate reading empty {max(miss) if miss else None},"
            f" empty candidates {miss}"
        )
        print(
            f"    candidates whose disc holds a zero within {args.guard} of the"
            f" occupancy circle {seam}"
        )
        out_got = [n for n in cand if grid[(round(t.imag, 6), n)][1]]
        out_miss = [n for n in cand if not grid[(round(t.imag, 6), n)][1]]
        print(
            f"    outer reading occupied {len(out_got)} of {len(cand)},"
            f" smallest singleton {out_got[0] if out_got else None},"
            f" last candidate reading empty"
            f" {max(out_miss) if out_miss else None}"
        )
        first = got[0] if got else None
        if first is None:
            print(f"    no singleton at or below {args.top} occupies it")
        else:
            ph = (t.imag * math.log(first) / (2.0 * math.pi)) % 1.0
            print(
                f"    smallest singleton {{{first}}} phase {ph:.6f}"
                f" predicted abs(u1) {abs(-r / (reg + cmath.exp(-t * math.log(first)))):.9f}"
                f" on the seam {first in seam}"
            )
        if lo and len(got) < len(cand):
            flips.append((round(t.imag, 6), miss))
    print(
        f"grid cells {cells}, first-order law agrees {agree},"
        f" disagrees {cells - agree}, edge guard splits {split}"
    )
    for i in sorted(perline):
        hit, tot, one, und, hitb, oneb = perline[i]
        print(
            f"  line {i} cells {tot} measured occupied {one},"
            f" first-order law agrees {hit}, the constant occupied predictor agrees {one}"
        )
        print(
            f"    occupancy undetermined {und}, outer reading occupied {oneb},"
            f" first-order law agrees {hitb},"
            f" the constant occupied predictor agrees {oneb}"
        )
    print(f"occupied teeth emptied by a singleton {len(flips)} {flips}")
    print("greedy against the regular part, the F that drives R + P_F to zero")
    for i, t, r, reg, u1, lo, hi, near in rows_out:
        if not lo:
            continue
        pool, pick, done = list(cand), [], None
        for _ in range(args.deep):
            best = min(
                pool,
                key=lambda n: abs(
                    reg + sum(cmath.exp(-t * math.log(m)) for m in pick + [n])
                ),
            )
            pick.append(best)
            pool.remove(best)
            eng = Engine(rule, args.shift, args.cut, pick)
            run = Census(eng, combs, wide, edges, rows, args.seed)
            a, b, _ = run.guarded(t, args.rho, args.guard)
            pv = abs(reg + sum(cmath.exp(-t * math.log(m)) for m in pick))
            print(
                f"  tooth line {i} Im {t.imag:10.6f} F {pick} abs(R + P_F) {pv:.9f}"
                f" predicted abs(u1) {abs(r) / pv:.9f} zeros in the disc {a}/{b}"
            )
            if a == 0 and b == 0:
                done = list(pick)
                break
        print(
            f"  tooth line {i} Im {t.imag:10.6f} emptied by {done}"
            f" size {len(done) if done else None}"
        )
        if args.deep and len(cand) <= 24:
            floor, best = exact_min(reg, t, cand, args.deep)
            eng = Engine(rule, args.shift, args.cut, best)
            run = Census(eng, combs, wide, edges, rows, args.seed)
            a, b, _ = run.guarded(t, args.rho, args.guard)
            print(
                f"  tooth line {i} Im {t.imag:10.6f} exact minimiser over the"
                f" {sum(math.comb(len(cand), e) for e in range(args.deep + 1))}"
                f" subsets of size at most {args.deep} {best}"
                f" abs(R + P_F) {floor:.9f} predicted abs(u1)"
                f" {abs(r) / floor:.9f} emptying threshold abs(r)/rho"
                f" {abs(r) / args.rho:.9f} zeros in the disc {a}/{b}"
            )
    print(f"runtime {time.perf_counter() - start:.2f} s")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=2)
    parser.add_argument("--code", type=int, default=7)
    parser.add_argument("--peel", type=int, default=8)
    parser.add_argument("--shift", type=float, default=12.0)
    parser.add_argument("--cut", type=int, default=24)
    parser.add_argument("--left", type=float, default=-0.95)
    parser.add_argument("--right", type=float, default=2.0)
    parser.add_argument("--height", type=float, default=43.1)
    parser.add_argument("--seed", type=float, default=0.05)
    parser.add_argument("--rho", type=float, default=0.45)
    parser.add_argument("--eps", type=float, default=0.3)
    parser.add_argument("--hug", type=float, default=0.05)
    parser.add_argument("--ring", type=float, default=0.05)
    parser.add_argument("--span", type=int, default=18)
    parser.add_argument("--top", type=int, default=40)
    parser.add_argument("--deep", type=int, default=8)
    parser.add_argument("--guard", type=float, default=0.02)
    parser.add_argument(
        "verb", choices=["census", "control", "teeth", "bridge", "dial"]
    )
    args = parser.parse_args()
    {
        "census": census,
        "control": control,
        "teeth": teeth,
        "bridge": bridge,
        "dial": dial,
    }[args.verb](args)


main()
