import argparse
import importlib.util
import math
import os
import time
from collections import Counter
from fractions import Fraction

import numpy as np
from mpmath import iv, mp

HERE = os.path.dirname(os.path.abspath(__file__))
JP_LO = "0.53128050627720514162446864736847178549305910901839877988839780392752953564383134591810957018118523987"
JP_HI = "0.53128050627720514162446864736847178549305910901839877988839780392752953564383134591810957018118523989"
PV = {
    (1, 3): "0.454489077661828743845",
    (2, 3): "0.337436780806063636304",
    (1, 2, 3): "0.705660908028738230607",
}
S1_NOTE = ("0.457015235231", "6.958882679527")
WALKS = 4 * 10**6


def sibling():
    path = os.path.join(HERE, "..", "question-mark", "question.py")
    spec = importlib.util.spec_from_file_location("question", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ivq(q):
    return iv.mpf(q.numerator) / iv.mpf(q.denominator)


def up(x):
    return mp.mpf(x._mpi_[1])


def low(x):
    return mp.mpf(x._mpi_[0])


def exact(x):
    man, ex = x.man, x.exp
    return Fraction(int(man)) * (Fraction(2) ** int(ex))


def fixed(q, d):
    sign = "-" if q < 0 else ""
    q = abs(q)
    whole = (
        q.numerator // q.denominator * 10**d
        + (q - q.numerator // q.denominator) * 10**d
    )
    n = int(whole)
    text = str(n).rjust(d + 1, "0")
    return f"{sign}{text[:-d]}.{text[-d:]}"


def bracket(lo, hi, d=None):
    lo, hi = exact(lo), exact(hi)
    if d is None:
        d = 1
        while Fraction(1, 10 ** (d + 1)) > (hi - lo) / 10 or d < 3:
            d += 1
    down = Fraction(math.floor(lo * 10**d), 10**d)
    upq = Fraction(math.ceil(hi * 10**d), 10**d)
    return f"[{fixed(down, d)}, {fixed(upq, d)}]"


def ceil_sig(x, k):
    q = x if isinstance(x, Fraction) else exact(mp.mpf(x))
    assert q > 0
    e = int(mp.floor(mp.log10(mp.mpf(q.numerator) / q.denominator)))
    while q >= Fraction(10) ** (e + 1):
        e += 1
    while q < Fraction(10) ** e:
        e -= 1
    n = math.ceil(q * Fraction(10) ** (k - 1 - e))
    if n >= 10**k:
        n, e = math.ceil(Fraction(n, 10)), e + 1
    t = str(n)
    if -4 <= e <= 5:
        d = max(0, k - 1 - e)
        return fixed(Fraction(n) * Fraction(10) ** (e - k + 1), d)
    return f"{t[0]}.{t[1:]}e{e}" if k > 1 else f"{t}e{e}"


def trunc(x, digits):
    return math.floor(x * 10**digits) / 10**digits


# GRAPH


def closed_walks(edges, period):
    states = sorted({u for u, _, _ in edges} | {v for _, v, _ in edges})
    count = Counter()
    for u0 in states:
        st = np.array([u0])
        mat = np.array([[1, 0, 0, 1]], dtype=np.int64)
        for n in range(1, period + 1):
            ns, nm = [], []
            for u, v, a in edges:
                pick = st == u
                if pick.any():
                    x = mat[pick]
                    nm.append(
                        np.stack(
                            [
                                x[:, 1],
                                x[:, 0] + a * x[:, 1],
                                x[:, 3],
                                x[:, 2] + a * x[:, 3],
                            ],
                            axis=1,
                        )
                    )
                    ns.append(np.full(len(x), v))
            st, mat = np.concatenate(ns), np.concatenate(nm)
            assert len(st) < WALKS and mat.max() < 2**60
            shut = st == u0
            t, c = np.unique(mat[shut, 0] + mat[shut, 3], return_counts=True)
            for x, y in zip(t, c):
                count[(n, int(x))] += int(y)
    return sorted(count.items())


def minimize(states, edges):
    out = {u: {} for u in states}
    for u, v, a in edges:
        out[u][a] = v
    part = {u: tuple(sorted(out[u])) for u in states}
    while True:
        new = {
            u: (part[u], tuple((a, part[out[u][a]]) for a in sorted(out[u])))
            for u in states
        }
        if len(set(new.values())) == len(set(part.values())):
            break
        part = new
    idx = {b: i for i, b in enumerate(sorted(set(part.values()), key=str))}
    return len(idx), sorted(
        {(idx[part[u]], idx[part[out[u][a]]], a) for u in states for a in out[u]}
    )


def period_of(n, edges):
    level = {0: 0}
    todo = [0]
    while todo:
        u = todo.pop()
        for x, y, _ in edges:
            if x == u and y not in level:
                level[y] = level[u] + 1
                todo.append(y)
    assert len(level) == n
    g = 0
    for x, y, _ in edges:
        g = math.gcd(g, level[x] + 1 - level[y])
    return g


def pieces(q, code, k):
    n = 1 << (k - 1)
    edges, tails = q.graph_of(code, k)
    comp = q.components(n, edges, tails)
    out = []
    for c in sorted(set(comp.values()), key=sorted):
        e = [(u, v, a) for (u, v), ls in edges.items() for a in ls if u in c and v in c]
        cof = any(u in c and v in c for u, v in tails.items())
        if cof:
            out.append(("tail", len(c), None))
        elif len(e) > len(c):
            s, me = minimize(sorted(c), e)
            out.append(("finite", len(c), (s, me)))
    return out


# DISC


def disc_numbers(edges, c, r, beta):
    labels = sorted({a for _, _, a in edges})
    ends = [abs(1 / (c - r + a) - c) for a in labels] + [
        abs(c - 1 / (c + r + a)) for a in labels
    ]
    return labels, max(ends) / r


def choose_disc(edges, s, states, period):
    th = np.linspace(0, 2 * np.pi, 721)
    best = None
    for c in np.arange(0.30, 0.92, 0.04):
        for r in np.arange(0.30, 1.30, 0.05):
            if c - r <= -0.9:
                continue
            labels, h = disc_numbers(edges, c, r, None)
            if h >= 0.95:
                continue
            z = c + r * np.exp(1j * th)
            for beta in np.arange(0.6, 4.01, 0.2):
                if c - r <= -beta + 0.05:
                    continue
                w = {
                    a: np.exp(
                        (2 * s * (np.log(z + beta) - np.log(beta * (z + a) + 1))).real
                    ).max()
                    for a in labels
                }
                sig = sum(
                    sum(w[a] for u, v, a in edges if u == x and v == y) ** 2
                    for x in range(states)
                    for y in range(states)
                )
                cc = math.sqrt(sig / (1 - h * h))
                m = np.arange(1, 60)
                e = np.concatenate([[1.0], np.cumprod(cc * h ** (m - 1) / (1 - h**m))])
                f = e
                for _ in range(states - 1):
                    f = np.convolve(f, e)[:60]
                tail = f[period + 1 :].sum()
                if best is None or tail < best[0]:
                    best = (tail, round(c, 2), round(r, 2), round(beta, 2))
    return tuple(Fraction(str(x)) for x in best[1:])


def weight_sup(labels, c, r, beta, sig, tau, pad, arcs):
    out = {}
    for a in labels:
        top = None
        for j in range(arcs):
            th = iv.mpf([2 * j, 2 * j + 2]) * iv.pi / arcs
            x = ivq(c) + ivq(r) * iv.cos(th)
            y = ivq(r) * iv.sin(th)
            if beta is None:
                u = x + a
                assert low(u) > 0
                re = -iv.log(u * u + y * y) / 2
                im = -iv.atan2(y, u)
            else:
                u1 = x + ivq(beta)
                u2 = ivq(beta) * x + ivq(beta * a + 1)
                v2 = ivq(beta) * y
                assert beta > 0 and low(u1) > 0 and low(u2) > 0
                re = (iv.log(u1 * u1 + y * y) - iv.log(u2 * u2 + v2 * v2)) / 2
                im = iv.atan2(y, u1) - iv.atan2(v2, u2)
            val = (2 * (sig * re - tau * im) + 2 * pad * iv.sqrt(re * re + im * im)).b
            top = val if top is None or val > top else top
        out[a] = iv.exp(iv.mpf(top))
    return out


def euler_tail(cc, h, states, period, radius):
    x = cc * radius
    e = [iv.mpf(1)]
    m = 0
    while True:
        m += 1
        e.append(e[-1] * x * h ** (m - 1) / (1 - h**m))
        if m > period + 2 and (x * h**m / (1 - h ** (m + 1))).b < 0.5:
            break
    rest = e[-1]
    f = e
    for _ in range(states - 1):
        f = [
            sum(
                (
                    f[i] * e[j - i]
                    for i in range(max(0, j - len(e) + 1), min(j, len(f) - 1) + 1)
                ),
                iv.mpf(0),
            )
            for j in range(len(f) + len(e) - 1)
        ]
    total = sum(e, iv.mpf(0)) + rest
    return up(sum(f[period + 1 :], iv.mpf(0)) + states * rest * total ** (states - 1))


def walk_cap(edges, states, period):
    adj = np.zeros((states, states))
    for u, v, _ in edges:
        adj[u, v] += 1
    growth = max(abs(np.linalg.eigvals(adj)))
    return min(period, int(math.log(WALKS / 2) / math.log(growth)))


class Certificate:
    def __init__(self, edges, states, period, point):
        period = walk_cap(edges, states, period)
        self.edges, self.states, self.period = edges, states, period
        self.labels = sorted({a for _, _, a in edges})
        self.c, self.r, self.beta = choose_disc(edges, point, states, period)
        assert self.c - self.r > -min(self.labels) and self.c - self.r > -self.beta
        self.hq = (
            max(
                max(
                    abs(1 / (self.c - self.r + a) - self.c),
                    abs(self.c - 1 / (self.c + self.r + a)),
                )
                for a in self.labels
            )
            / self.r
        )
        assert self.hq < 1
        self.h = ivq(self.hq)
        self.groups = closed_walks(edges, period)
        self.iv_rows = []
        self.mp_rows = []
        for (n, t), m in self.groups:
            e = (-1) ** n
            mu = (t + iv.sqrt(iv.mpf(t * t - 4 * e))) / 2
            self.iv_rows.append((n, iv.log(mu), m / (1 - e / (mu * mu))))
            mum = (t + mp.sqrt(t * t - 4 * e)) / 2
            self.mp_rows.append((n, mp.log(mum), m / (1 - e / (mum * mum))))

    def constant(self, sig, tau, pad, arcs):
        w = weight_sup(self.labels, self.c, self.r, self.beta, sig, tau, pad, arcs)
        sig2 = iv.mpf(0)
        for x in range(self.states):
            for y in range(self.states):
                sig2 += (
                    sum(
                        (w[a] for u, v, a in self.edges if u == x and v == y), iv.mpf(0)
                    )
                    ** 2
                )
        return iv.sqrt(sig2 / (1 - self.h * self.h))

    def tail(self, sig, tau, pad, radius=1, arcs=720):
        return euler_tail(
            self.constant(sig, tau, pad, arcs), self.h, self.states, self.period, radius
        )

    def coefficients(self, s, ctx, deriv=False):
        rows = self.iv_rows if ctx is iv else self.mp_rows
        zero = ctx.mpc(0) if isinstance(s, (ctx.mpc,)) else ctx.mpf(0)
        tr = [zero] * (self.period + 1)
        dtr = [zero] * (self.period + 1)
        for n, lg, f in rows:
            x = f * ctx.exp(-2 * s * lg)
            tr[n] += x
            if deriv:
                dtr[n] += -2 * lg * x
        d, dd = [zero + 1], [zero]
        for n in range(1, self.period + 1):
            d.append(-sum((tr[k] * d[n - k] for k in range(1, n + 1)), zero) / n)
            if deriv:
                dd.append(
                    -sum(
                        (
                            dtr[k] * d[n - k] + tr[k] * dd[n - k]
                            for k in range(1, n + 1)
                        ),
                        zero,
                    )
                    / n
                )
        return d, dd

    def newton(self, s, steps=8):
        for _ in range(steps):
            d, dd = self.coefficients(s, mp, True)
            s = s - sum(d) / sum(dd)
        return s

    def winding(self, d, radius, bound, arcs=360):
        def horner(z):
            v = d[-1]
            for coef in reversed(d[:-1]):
                v = v * z + coef
            return v

        planes = []
        pts = []
        for j in range(arcs):
            th = iv.mpf([2 * j, 2 * j + 2]) * iv.pi / arcs
            v = horner(radius * iv.mpc(iv.cos(th), iv.sin(th)))
            tests = [
                (0, 1, low(v.real) - bound > 0),
                (0, -1, up(v.real) + bound < 0),
                (1, 1, low(v.imag) - bound > 0),
                (1, -1, up(v.imag) + bound < 0),
            ]
            plane = next((a, b) for a, b, ok in tests if ok)
            planes.append(plane)
            t = iv.mpf(2 * j) * iv.pi / arcs
            u = horner(radius * iv.mpc(iv.cos(t), iv.sin(t)))
            pts.append(complex(float(u.real.mid), float(u.imag.mid)))
        turn = 0.0
        for j in range(arcs):
            a, b = planes[j]
            for p in (pts[j], pts[(j + 1) % arcs]):
                assert b * (p.real if a == 0 else p.imag) > 0
            q = pts[(j + 1) % arcs] / pts[j]
            turn += math.atan2(q.imag, q.real)
        count = round(turn / (2 * math.pi))
        assert abs(turn / (2 * math.pi) - count) < 0.01
        return count


def top_zero(cert):
    f = lambda x: sum(cert.coefficients(mp.mpf(x), mp)[0])
    hi = 1.0
    fh = f(hi)
    assert fh > 0
    while True:
        lo = hi - 0.01
        assert lo > 0
        fl = f(lo)
        if fl < 0:
            break
        hi, fh = lo, fl
    for _ in range(30):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return cert.newton(mp.mpf((lo + hi) / 2), 4)


def certify_real(cert, guess, width_floor):
    s = cert.newton(mp.mpf(guess), 4) if guess is not None else top_zero(cert)
    cc = cert.constant(iv.mpf(s), iv.mpf(0), iv.mpf("1e-6"), 720)
    t1 = euler_tail(cc, cert.h, cert.states, cert.period, 1)
    t2 = euler_tail(cc, cert.h, cert.states, cert.period, 2)
    _, dd = cert.coefficients(s, mp, True)
    r = mp.mpf(mp.nstr(max(10 * t1 / abs(sum(dd)), width_floor), 2))
    assert r < mp.mpf("1e-7")
    vals = []
    for x in (s - r, s + r):
        d, _ = cert.coefficients(iv.mpf(x), iv)
        assert cert.winding(d, 2, t2) == cert.period_graph
        vals.append(sum(d, iv.mpf(0)))
    assert (vals[0] + t1).b < 0 < (vals[1] - t1).a, (vals, t1)
    return s - r, s + r, t1


# ZERO


def zero(period, arcs):
    t0 = time.time()
    mp.dps = iv.dps = 90
    cert = Certificate(
        [(0, 0, 1), (0, 0, 2)], 1, period, complex(0.457015235231, 6.958882679527)
    )
    cert.period_graph = 1
    h = float(cert.hq)
    print(
        f"zero: A = {{1,2}}, periods 1..{period}, {len(cert.groups)} (period, trace) classes, disc centre {cert.c} radius {cert.r}, conjugation beta {cert.beta}, contraction h = {cert.hq} = {h:.6f}"
    )
    s0 = cert.newton(mp.mpc(*S1_NOTE))
    S0 = iv.mpc(s0.real, s0.imag)
    kk = cert.constant(S0.real, S0.imag, iv.mpf(0), arcs) / cert.h
    tail0 = cert.tail(S0.real, S0.imag, iv.mpf(0), arcs=arcs)
    f0, fp0 = [sum(x, iv.mpc(0)) for x in cert.coefficients(S0, iv, True)]
    r = mp.mpf(mp.nstr(10 * tail0 / abs(mp.mpc(fp0.real.mid, fp0.imag.mid)), 2))
    rho2 = mp.mpf("0.001")
    tail1 = up(
        iv.mpf(cert.tail(S0.real, S0.imag, iv.mpf(rho2) + 4 * r, arcs=arcs)) / rho2
    )
    box = iv.mpf([-r, r])
    loose = S0 + iv.mpc(box, box)
    re_text = bracket(low(loose.real), up(loose.real))
    im_text = bracket(low(loose.imag), up(loose.imag))
    X = iv.mpc(iv.mpf(re_text[1:-1].split(", ")), iv.mpf(im_text[1:-1].split(", ")))
    assert up(X.real) - low(X.real) < 3 * r and up(X.imag) - low(X.imag) < 3 * r
    half_re, half_im = [
        (Fraction(t[1:-1].split(", ")[1]) - Fraction(t[1:-1].split(", ")[0])) / 2
        for t in (re_text, im_text)
    ]
    _, fpX = [sum(x, iv.mpc(0)) for x in cert.coefficients(X, iv, True)]
    y = 1 / mp.mpc(fp0.real.mid, fp0.imag.mid)
    Y = iv.mpc(y.real, y.imag)
    b0 = iv.mpf([-tail0, tail0])
    b1 = iv.mpf([-tail1, tail1])
    K = S0 - Y * (f0 + iv.mpc(b0, b0)) + (1 - Y * (fpX + iv.mpc(b1, b1))) * (X - S0)
    inside = (
        K.real.a > X.real.a
        and K.real.b < X.real.b
        and K.imag.a > X.imag.a
        and K.imag.b < X.imag.b
    )
    assert inside
    print(
        f"zero: weight constant K <= {ceil_sig(up(kk), 6)} on {arcs} arcs, Euler tail past period {period} at s_1 <= {ceil_sig(tail0, 3)}, derivative tail on the box <= {ceil_sig(tail1, 3)}"
    )
    print(
        f"zero: |D_{period}(s_0)| = {mp.nstr(abs(mp.mpc(f0.real.mid, f0.imag.mid)), 3)}, D_{period}'(s_0) = {mp.nstr(mp.mpc(fp0.real.mid, fp0.imag.mid), 12)}, radius before rounding {mp.nstr(r, 2)}, printed box half-widths {ceil_sig(half_re, 2)} (Re) and {ceil_sig(half_im, 2)} (Im)"
    )
    print(
        f"zero: Krawczyk image, on the printed box, inside it: {inside}; exactly one zero of D in"
    )
    print(f"zero:   Re s in {re_text}")
    print(f"zero:   Im s in {im_text}")
    re_t, im_t = trunc(float(X.real.a), 12), trunc(float(X.imag.a), 12)
    note = (float(S1_NOTE[0]), float(S1_NOTE[1]))
    rounds = (
        abs(note[0] - float(s0.real)) < 5e-13 and abs(note[1] - float(s0.imag)) < 5e-13
    )
    assert rounds
    print(
        f"zero: the note reads {S1_NOTE[0]} + {S1_NOTE[1]} i, the box rounded to twelve digits: {rounds}"
    )
    jc, jr = (
        Fraction("0.758687144013554292899790137015621955739"),
        Fraction("0.957589818521375342814351002388265920293"),
    )
    jh = ivq(
        max(
            max(abs(1 / (jc - jr + a) - jc), abs(jc - 1 / (jc + jr + a)))
            for a in (1, 2)
        )
        / jr
    )
    for sig, tau, name in (
        (iv.mpf(JP_LO), iv.mpf(0), "dim E_2"),
        (S0.real, S0.imag, "s_1"),
    ):
        ks = []
        for beta in (None, Fraction(33, 20), cert.beta):
            ws = weight_sup([1, 2], jc, jr, beta, sig, tau, iv.mpf(0), arcs)
            ks.append(ceil_sig(up((ws[1] + ws[2]) / (jh * iv.sqrt(1 - jh * jh))), 7))
        print(
            f"zero: on the Jenkinson-Pollicott disc, h = {mp.nstr(up(jh), 6)}, upper bounds on the constant K at {name}, rounded up: {ks[0]} unconjugated, {ks[1]} at beta = 33/20, {ks[2]} at beta = {cert.beta}"
        )
    print(f"zero: {time.time() - t0:.1f} s")
    return X, re_t, im_t


def control(period, arcs):
    t0 = time.time()
    mp.dps = iv.dps = 90
    rows = [((1, 2), JP_LO, JP_HI, "Jenkinson-Pollicott 2018 Theorem 1")]
    rows += [
        (a, PV[a], PV[a], "Pollicott-Vytnova 2022 Table 3, +- 1e-20")
        for a in ((1, 3), (2, 3), (1, 2, 3))
    ]
    for alpha, lo, hi, src in rows:
        cert = Certificate([(0, 0, a) for a in alpha], 1, period, 0.5)
        cert.period_graph = 1
        guess = mp.mpf(lo)
        slo, shi, tail = certify_real(cert, guess, mp.mpf("1e-60"))
        if src.startswith("Pollicott"):
            ok = (
                slo < mp.mpf(lo) + mp.mpf("1e-20")
                and mp.mpf(lo) - mp.mpf("1e-20") < shi
            )
        else:
            ok = slo < mp.mpf(lo) and mp.mpf(hi) < shi
        print(
            f"control: A = {{{','.join(map(str, alpha))}}}, periods 1..{cert.period}, disc ({cert.c}, {cert.r}, beta {cert.beta}), h = {float(cert.hq):.6f}, tail {mp.nstr(tail, 3)}: dim in {bracket(slo, shi)}, {src} inside: {ok}"
        )
        assert ok
    print(f"control: {time.time() - t0:.1f} s")


# ORPHANS


def orphans(period, modes):
    t0 = time.time()
    mp.dps = iv.dps = 60
    q = sibling()
    x, w = q.nodes(modes)
    B = q.branch_table(x, w, 2000)
    rowsq = q.taylor_rows(x, w, 4)
    _, _, _, live = q.symmetric_orphans(4)
    order = [11892] + [c for c in live if c != 11892]
    print(
        "| code | states | minimal | period | P | h at most | tail at most | dim_CF in | collocation |"
    )
    print("|---|---|---|---|---|---|---|---|---|")
    walls = []
    kinds = {}
    for k, code in [(4, c) for c in order] + [(3, 54)]:
        parts = pieces(q, code, k)
        if any(p[0] == "tail" for p in parts):
            walls.append(code)
            continue
        colloc, _ = q.graph_zero(code, k, x, B, rowsq, 2000)
        found = []
        seen = set()
        for _, n, (s, me) in parts:
            if (s, tuple(me)) in seen:
                continue
            seen.add((s, tuple(me)))
            cert = Certificate(me, s, period, colloc)
            cert.period_graph = period_of(s, me)
            try:
                slo, shi, tail = certify_real(cert, colloc, mp.mpf("1e-45"))
            except AssertionError:
                slo, shi, tail = certify_real(cert, None, mp.mpf("1e-45"))
            found.append((slo, shi, tail, n, s, cert))
        slo, shi, tail, n, s, cert = max(found, key=lambda f: f[0])
        assert all(f[1] < slo for f in found if f[5] is not cert)
        reading = Fraction(math.floor(Fraction(colloc) * 10**12), 10**12)
        cut = [Fraction(math.floor(exact(x) * 10**12), 10**12) for x in (slo, shi)]
        kind = (
            "truncates"
            if cut[0] == cut[1] == reading
            else ("contains" if exact(slo) <= reading <= exact(shi) else None)
        )
        ok = kind is not None
        kinds[kind] = kinds.get(kind, 0) + 1
        print(
            f"| {code} | {n} | {s} | {cert.period_graph} | {cert.period} | {fixed(Fraction(math.ceil(cert.hq * 10**4), 10**4), 4)} | {ceil_sig(tail, 2)} | {bracket(slo, shi)} | {trunc(colloc, 12):.12f} |"
        )
        assert ok, (code, slo, shi, colloc)
    print(
        f"orphans: the printed twelve-digit reading is the bracket truncated at twelve digits on {kinds.get('truncates', 0)} rows and lies inside the bracket on {kinds.get('contains', 0)}"
    )
    print(
        f"orphans: a cofinite edge on a cycle, outside the periodic-point determinant: {walls}"
    )
    print(f"orphans: {time.time() - t0:.1f} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("verb", choices=["zero", "control", "orphans", "all"])
    ap.add_argument("--period", type=int, default=18)
    ap.add_argument("--arcs", type=int, default=2000)
    ap.add_argument("--modes", type=int, default=40)
    a = ap.parse_args()
    if a.verb in ("zero", "all"):
        zero(a.period, a.arcs)
    if a.verb in ("control", "all"):
        control(a.period, a.arcs)
    if a.verb in ("orphans", "all"):
        orphans(a.period + 4, a.modes)


if __name__ == "__main__":
    main()
