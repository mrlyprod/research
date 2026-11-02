from fractions import Fraction

import mpmath as mp
import numpy as np

TOP = 11
MC_TOP = 7
MC_SAMPLES = 10**7
BATCH = 10**6


def wallis_table(top):
    table = {}
    even = Fraction(16, 3)
    odd = Fraction(3, 8)
    table[2] = ("r/Pi^2", even)
    table[3] = ("rational", odd)
    k = 1
    while 2 * k + 2 <= top or 2 * k + 3 <= top:
        even *= Fraction(4 * k * (k + 1), (2 * k + 1) * (2 * k + 3))
        odd *= Fraction((2 * k + 1) * (2 * k + 3), (2 * k + 2) * (2 * k + 4))
        if 2 * k + 2 <= top:
            table[2 * k + 2] = ("r/Pi^2", even)
        if 2 * k + 3 <= top:
            table[2 * k + 3] = ("rational", odd)
        k += 1
    return table


def value(form, rational):
    if form == "r/Pi^2":
        return mp.mpf(rational.numerator) / mp.mpf(rational.denominator) / mp.pi**2
    return mp.mpf(rational.numerator) / mp.mpf(rational.denominator)


def half_ball(rng, count, dim):
    pts = rng.standard_normal((count, dim))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    pts *= rng.random((count, 1)) ** (1.0 / dim)
    pts[:, dim - 1] = np.abs(pts[:, dim - 1])
    return pts


def flat_face_rate(rng, dim, samples):
    hits = 0
    done = 0
    while done < samples:
        take = min(BATCH, samples - done)
        p = half_ball(rng, take, dim)
        q = half_ball(rng, take, dim)
        gap = p[:, dim - 1] - q[:, dim - 1]
        good = gap != 0.0
        step = np.where(good, p[:, dim - 1] / np.where(good, gap, 1.0), 0.0)
        cross = p[:, : dim - 1] + step[:, None] * (q[:, : dim - 1] - p[:, : dim - 1])
        hits += int(np.count_nonzero(good & (np.linalg.norm(cross, axis=1) <= 1.0)))
        done += take
    return hits / samples


def main():
    mp.mp.dps = 30
    table = wallis_table(TOP)
    print("VERSION L: two uniform points, the line through them meets the flat face")
    print("  seeds f(2) = 16/(3 Pi^2) and f(3) = 3/8")
    print("  even step f(2k+2)/f(2k) = 4k(k+1)/((2k+1)(2k+3))")
    print("  odd step  f(2k+3)/f(2k+1) = (2k+1)(2k+3)/((2k+2)(2k+4))")
    print(f"  exact rationals to d = {TOP}")
    for d in sorted(table):
        form, r = table[d]
        v = value(form, r)
        shown = f"({r})/Pi^2" if form == "r/Pi^2" else f"{r}"
        print(f"  d = {d:>2}  {shown:>26}  {form:>9}  {mp.nstr(v, 12)}")

    print("  page decimals")
    for d, digits in ((2, 7), (3, 3), (4, 7), (5, 6), (6, 7), (7, 10)):
        form, r = table[d]
        print(f"  d = {d}  {mp.nstr(value(form, r), digits)}")

    rng = np.random.default_rng(20250828)
    print(f"  random-point check, {MC_SAMPLES} samples per dimension")
    for d in range(2, MC_TOP + 1):
        form, r = table[d]
        exact = value(form, r)
        est = mp.mpf(flat_face_rate(rng, d, MC_SAMPLES))
        se = mp.sqrt(est * (1 - est) / MC_SAMPLES)
        print(
            f"  d = {d}  estimate {mp.nstr(est, 8)}  exact {mp.nstr(exact, 8)}"
            f"  deviation {mp.nstr(est - exact, 3)}  one sigma {mp.nstr(se, 3)}"
        )


main()
