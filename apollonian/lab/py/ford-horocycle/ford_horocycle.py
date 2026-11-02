import argparse
import math
import sys
import time
from fractions import Fraction

import numpy as np

DELTA = {2: 0.5312805062772051416, 3: 0.705660908028}

# TOTIENTS


def totients(n):
    phi = np.arange(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if phi[p] == p:
            phi[p::p] -= phi[p::p] // p
    return phi


# QUOTIENTS


def quotients(a, b):
    out = []
    while b:
        q, r = divmod(a, b)
        out.append(q)
        a, b = b, r
    return out[1:]


def member(a, b, m):
    qs = quotients(a, b)
    if not qs:
        return True
    alt = qs[:-1] + [qs[-1] - 1, 1]
    return max(qs) <= m or max(alt) <= m


def rule(a, b, m):
    qs = quotients(a, b)
    if not qs:
        return True
    return max(qs[:-1], default=1) <= m and qs[-1] <= m + 1


# CROSSINGS


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def mobius_image(a, ap, b, d, z):
    x, y = z
    num = (Fraction(a) * x + ap, Fraction(a) * y)
    den = (Fraction(b) * x + d, Fraction(b) * y)
    n2 = den[0] ** 2 + den[1] ** 2
    re = (num[0] * den[0] + num[1] * den[1]) / n2
    im = (num[1] * den[0] - num[0] * den[1]) / n2
    return re, im


def crossings(bound, heights):
    phi = totients(bound)
    checked = 0
    for b in range(1, bound + 1):
        for a in range(b):
            if math.gcd(a, b) != 1:
                continue
            _g, s, t = egcd(a, b)
            ap, d = -t, s
            assert a * d - ap * b == 1
            centre = (Fraction(a, b), Fraction(1, 2 * b * b))
            r2 = Fraction(1, 4 * b**4)
            for x in (
                Fraction(0),
                Fraction(1, 3),
                Fraction(-7, 2),
                Fraction(5),
                Fraction(-d, b),
            ):
                re, im = mobius_image(a, ap, b, d, (x, Fraction(1)))
                assert (re - centre[0]) ** 2 + (im - centre[1]) ** 2 == r2
                assert im <= Fraction(1, b * b)
                if x == Fraction(-d, b):
                    assert im == Fraction(1, b * b) and re == Fraction(a, b)
                checked += 1
    print(
        f"horocycle images checked exactly: {checked} points on {int(phi[1 : bound + 1].sum())} Ford circles, b <= {bound}, every image on the Ford circle over a/b, top point 1/b^2 at x = -d/b"
    )
    print(
        "Q  height  closed-disk count  sum_{b<=Q} phi  open-disk count  sum_{b<Q} phi  count at 1/(2Q^2)  sum_{b<=floor(Q sqrt2)} phi"
    )
    for Q in heights:
        h = Fraction(1, Q * Q)
        closed = sum(
            int(phi[b]) for b in range(1, bound + 1) if Fraction(1, b * b) >= h
        )
        opened = sum(int(phi[b]) for b in range(1, bound + 1) if Fraction(1, b * b) > h)
        half = sum(
            int(phi[b]) for b in range(1, bound + 1) if Fraction(1, b * b) >= h / 2
        )
        mQ = int(phi[1 : Q + 1].sum())
        mQ1 = int(phi[1:Q].sum())
        m2 = int(phi[1 : math.isqrt(2 * Q * Q) + 1].sum())
        assert closed == mQ and opened == mQ1 and half == m2
        print(f"{Q}  1/{Q * Q}  {closed}  {mQ}  {opened}  {mQ1}  {half}  {m2}")


# BRIDGE


def bump_smooth(u):
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u)
    inside = (u > 1) & (u < 2)
    w = u[inside]
    out[inside] = np.exp(4 - 1 / ((w - 1) * (2 - w)))
    return out


def bump_c2(u):
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u)
    inside = (u > 1) & (u < 2)
    w = u[inside]
    out[inside] = 64 * (w - 1) ** 3 * (2 - w) ** 3
    return out


def gl(nodes):
    x, w = np.polynomial.legendre.leggauss(nodes)
    return x, w


def integrate(f, lo, hi, x, w):
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    return half * np.dot(w, f(mid + half * x))


def g_transform(f, t, x, w):
    if t >= 1:
        return 0.0
    hi = math.sqrt(1 / (t * t) - 1)
    lo = math.sqrt(1 / (2 * t * t) - 1) if t * t <= 0.5 else 0.0
    return 2 * integrate(lambda v: f(1 / (t * t * (1 + v * v))), lo, hi, x, w)


def horocycle_integral(f, h, x, w):
    Q = math.isqrt(int(1 / h))
    total = 0.0
    bumps = 0
    for c in range(1, Q + 1):
        s = h - c * c * h * h
        if s < 0:
            continue
        width = math.sqrt(s) / c
        inner = math.sqrt(max(h / 2 - c * c * h * h, 0.0)) / c
        dlo = math.ceil(-c * (1 + width))
        dhi = math.floor(c * width)
        for d in range(dlo, dhi + 1):
            if math.gcd(c, abs(d)) != 1:
                continue
            centre = -d / c
            pieces = (
                (centre - width, centre - inner),
                (centre + inner, centre + width),
            )
            hit = False
            for lo, hi in pieces:
                lo, hi = max(0.0, lo), min(1.0, hi)
                if lo >= hi:
                    continue
                total += integrate(
                    lambda t: f(h / ((c * t + d) ** 2 + c * c * h * h)), lo, hi, x, w
                )
                hit = True
            bumps += hit
    return total, bumps


def bridge(kmax, nodes):
    x, w = gl(nodes)
    for name, f in (("smooth bump", bump_smooth), ("C2 bump", bump_c2)):
        main = 3 / math.pi * integrate(lambda u: f(u) / (u * u), 1, 2, x, w)
        print(f"{name}: main term (3/pi) int f(w) w^-2 dw = {main:.15f}")
        print(
            "k  h=4^-k  Q  bumps  horocycle integral  h S_g(h^1/2)  abs gap  error E  E/h^(3/4)"
        )
        Qmax = 2**kmax
        phi = totients(Qmax)
        for k in range(2, kmax + 1):
            h = 4.0 ** (-k)
            Q = 2**k
            t0 = time.time()
            lhs, bumps = horocycle_integral(f, h, x, w)
            rhs = h * sum(
                int(phi[c]) * g_transform(f, c * math.sqrt(h), x, w)
                for c in range(1, Q + 1)
            )
            err = lhs - main
            print(
                f"{k}  {h:.3e}  {Q}  {bumps}  {lhs:.15f}  {rhs:.15f}  {abs(lhs - rhs):.1e}  {err:+.3e}  {err / h**0.75:+.4f}  ({time.time() - t0:.1f}s)"
            )


# CENSUS


def census_walk(m, qmax):
    hist = np.zeros(qmax + 1, dtype=np.int64)
    hist[1] = 1
    stack = [(1, 2, 0, 1, 1, 1, 0, 1)]
    while stack:
        p, q, p1, q1, p2, q2, last, run = stack.pop()
        hist[q] += 1
        ql = q + q1
        if ql <= qmax:
            r = run + 1 if last == 0 else 1
            if r <= m:
                stack.append((p + p1, ql, p1, q1, p, q, 0, r))
        qr = q + q2
        if qr <= qmax:
            r = run + 1 if last == 1 else 1
            if r <= m:
                stack.append((p + p2, qr, p, q, p2, q2, 1, r))
    return np.cumsum(hist)


def census_nodes(m, qmax):
    out = []
    stack = [(1, 2, 0, 1, 1, 1, 0, 1)]
    while stack:
        node = stack.pop()
        p, q, p1, q1, p2, q2, last, run = node
        out.append(node)
        ql, qr = q + q1, q + q2
        if ql <= qmax:
            r = run + 1 if last == 0 else 1
            if r <= m:
                stack.append((p + p1, ql, p1, q1, p, q, 0, r))
        if qr <= qmax:
            r = run + 1 if last == 1 else 1
            if r <= m:
                stack.append((p + p2, qr, p, q, p2, q2, 1, r))
    return out


def census(jmax2, jmax3, jcontrol, jcheck):
    qc = 2**jcontrol
    t0 = time.time()
    cum = census_walk(10**9, qc)
    phi = totients(qc)
    tot = np.cumsum(phi)
    assert all(int(cum[2**j]) == int(tot[2**j]) for j in range(jcontrol + 1))
    print(
        f"control A = N to Q = 2^{jcontrol}: the walk's count is sum_(b <= Q) phi(b) at every Q = 2^j, {int(cum[qc])} at the top ({time.time() - t0:.1f}s)"
    )
    qk = 2**jcheck
    for m in (2, 3):
        walked = {(p, q) for p, q, *_ in census_nodes(m, qk)} | {(0, 1)}
        brute = {
            (a, b)
            for b in range(1, qk + 1)
            for a in range(b)
            if math.gcd(a, b) == 1 and member(a, b, m)
        }
        ruled = {
            (a, b)
            for b in range(1, qk + 1)
            for a in range(b)
            if math.gcd(a, b) == 1 and rule(a, b, m)
        }
        assert walked == brute == ruled
        both = one = 0
        for p, q, p1, q1, p2, q2, last, run in census_nodes(m, qk):
            qs = quotients(p, q)
            an = qs[-1]
            assert run == an - 1
            older, younger = ((p1, q1), (p2, q2)) if last == 0 else ((p2, q2), (p1, q1))
            assert older[1] < younger[1] or (p, q) == (1, 2)
            with_younger = member(p + younger[0], q + younger[1], m)
            with_older = member(p + older[0], q + older[1], m)
            assert with_younger
            assert with_older == (an <= m)
            if with_older:
                both += 1
            else:
                one += 1
        print(
            f"m = {m}, b <= 2^{jcheck}: walk, either-expansion test and the rule a_1..a_(n-1) <= m, a_n <= m+1 agree on {len(walked)} fractions; mediant with the younger neighbour stays in {both + one} of {both + one}, with the older in {both}, exactly the nodes with a_n <= m, and refused in {one}, exactly a_n = m + 1"
        )
    for m, jmax in ((2, jmax2), (3, jmax3)):
        qmax = 2**jmax
        t0 = time.time()
        cum = census_walk(m, qmax)
        elapsed = time.time() - t0
        d = 2 * DELTA[m]
        print(
            f"m = {m}, A = {{1..{m}}}, 2 delta = {d:.13f}, walk to Q = 2^{jmax}: {int(cum[qmax])} fractions in [0, 1) ({elapsed:.1f}s)"
        )
        print(
            "j  Q  N(Q)  Q^(2 delta)  N/Q^(2 delta)  octave exponent log2(N(2Q)/N(Q))  octave exponent minus 2 delta  two-octave exponent log4(N(4Q)/N(Q))  min and max of N/Q^(2 delta) on the quarter octaves of [j, j+1)"
        )
        for j in range(jmax + 1):
            Q = 2**j
            n = int(cum[Q])
            power = Q**d
            one = f"{math.log2(int(cum[2 * Q]) / n):.4f}" if j + 1 <= jmax else "-"
            dev = f"{math.log2(int(cum[2 * Q]) / n) - d:+.4f}" if j + 1 <= jmax else "-"
            two = f"{math.log2(int(cum[4 * Q]) / n) / 2:.4f}" if j + 2 <= jmax else "-"
            band = [
                int(cum[round(2 ** (j + k / 4))]) / 2 ** ((j + k / 4) * d)
                for k in range(4)
                if j + k / 4 <= jmax
            ]
            print(
                f"{j}  {Q}  {n}  {power:.3f}  {n / power:.5f}  {one}  {dev}  {two}  {min(band):.5f} {max(band):.5f}"
            )


# MAIN


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("verb", choices=["crossings", "bridge", "census", "all"])
    ap.add_argument("--bound", type=int, default=60)
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--nodes", type=int, default=200)
    ap.add_argument("--jmax2", type=int, default=24)
    ap.add_argument("--jmax3", type=int, default=18)
    ap.add_argument("--jcontrol", type=int, default=10)
    ap.add_argument("--jcheck", type=int, default=9)
    args = ap.parse_args()
    t0 = time.time()
    if args.verb in ("crossings", "all"):
        crossings(args.bound, [1, 2, 3, 10, 32, 42])
    if args.verb in ("bridge", "all"):
        bridge(args.kmax, args.nodes)
    if args.verb in ("census", "all"):
        census(args.jmax2, args.jmax3, args.jcontrol, args.jcheck)
    print(f"total {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
