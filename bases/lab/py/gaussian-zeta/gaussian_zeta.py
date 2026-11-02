import math
import time
from fractions import Fraction

import numpy as np
from mpmath import bernoulli, catalan, mp, mpf, zeta

DIGITS = 90
CORNERS = [(0, 0), (0, 1), (1, 0), (1, 1)]
DIAGONAL = [(0, 0), (1, 1)]
MATE = [(0, 1), (1, 0)]
EXACT_TO = 80
LIMITS = [49, 199, 999]
SHELL_MAX = 1000000
RADIUS = 6000
STATED = {
    1: "1/24",
    2: "1/4",
    3: "7/24",
    4: "1/4",
    5: "7/24",
    6: "1/2",
    7: "13/24",
    8: "1/8",
    9: "1/6",
    10: "3/8",
    11: "5/12",
    12: "3/8",
    13: "5/12",
    14: "5/8",
    15: "2/3",
}


def code_of(filled):
    return sum(1 << k for k, cell in enumerate(CORNERS) if cell in filled)


def fill(filled, n):
    table = np.zeros((2, 2), dtype=bool)
    for i, j in filled:
        table[i, j] = True
    parity = np.arange(n) % 2
    return int(np.count_nonzero(table[parity[:, None], parity[None, :]]))


def fluctuation(filled, n):
    return Fraction(fill(filled, n), n * n) - Fraction(len(filled), 4)


def pi_machin():
    def arctan_inv(x):
        total = mpf(0)
        term = mpf(1) / x
        k = 0
        while term > mpf(10) ** (-DIGITS - 10):
            total += term / (2 * k + 1) * (1 if k % 2 == 0 else -1)
            term /= x * x
            k += 1
        return total

    return 16 * arctan_inv(5) - 4 * arctan_inv(239)


def odd_zeta(s, m=600, k=14):
    total = sum(mpf(1) / mpf(2 * j - 1) ** s for j in range(1, m))
    a = mpf(2 * m - 1)
    total += a ** (1 - s) / (2 * (s - 1)) + a ** (-s) / 2
    for step in range(1, k + 1):
        poch = mpf(1)
        for i in range(2 * step - 1):
            poch *= s + i
        derivative = -(mpf(2) ** (2 * step - 1)) * poch * a ** (-(s + 2 * step - 1))
        total -= bernoulli(2 * step) / mp.factorial(2 * step) * derivative
    return total


def lambda_part():
    print("DIAGONAL DESIGN, FILLED CORNERS (0,0) AND (1,1), CODE", code_of(DIAGONAL))
    even_ok = True
    odd_ok = True
    for n in range(1, EXACT_TO + 1):
        value = fluctuation(DIAGONAL, n)
        if n % 2 == 0:
            even_ok = even_ok and value == 0
        else:
            odd_ok = odd_ok and value == Fraction(1, 2 * n * n)
        if n <= 11:
            print(f"  n={n:2d}  fill {fill(DIAGONAL, n):4d}  fluctuation {value}")
    print(f"  even n <= {EXACT_TO} fluctuation exactly 0: {even_ok}")
    print(f"  odd n <= {EXACT_TO} fluctuation exactly 1/(2n^2): {odd_ok}")
    forms_ok = all(
        fill(DIAGONAL, 2 * m) == 2 * m * m
        and fill(DIAGONAL, 2 * m - 1) == m * m + (m - 1) ** 2
        for m in range(1, EXACT_TO // 2 + 1)
    )
    print(f"  fill(2m) = 2m^2 and fill(2m-1) = m^2 + (m-1)^2: {forms_ok}")
    mate_ok = all(
        fluctuation(MATE, n) == (0 if n % 2 == 0 else Fraction(-1, 2 * n * n))
        for n in range(1, 41)
    )
    print(
        f"  orbit mate code {code_of(MATE)} fluctuation -1/(2n^2) on odd n <= 40: {mate_ok}"
    )
    print()

    print("EXACT RATIONAL REDUCTION")
    z2 = Fraction(1, 2) * (1 - Fraction(1, 16)) / 90
    z4 = Fraction(1, 2) * (1 - Fraction(1, 64)) / 945
    print(f"  (1/2)(1-2^-4)/90  = {z2}   target 1/192   {z2 == Fraction(1, 192)}")
    print(f"  (1/2)(1-2^-6)/945 = {z4}   target 1/1920  {z4 == Fraction(1, 1920)}")
    print()

    print(f"CONSTANTS AT {DIGITS} DIGITS")
    mp.dps = DIGITS
    machin = pi_machin()
    print(f"  pi machin          {mp.nstr(machin, 60)}")
    print(f"  pi library         {mp.nstr(mp.pi, 60)}")
    print(f"  pi difference      {mp.nstr(abs(machin - mp.pi), 3)}")
    for s in (4, 6):
        direct = odd_zeta(s)
        identity = (1 - mpf(2) ** (-s)) * zeta(s)
        print(f"  lambda({s}) two routes differ {mp.nstr(abs(direct - identity), 3)}")
    targets = {}
    cases = ((2, 4, 192, "0.50733901580"), (4, 6, 1920, "0.50072353832"))
    for s, power, denom, prefix in cases:
        left = odd_zeta(s + 2) / 2
        right = machin**power / denom
        targets[s] = right
        text = mp.nstr(left, 50)
        print(f"  Z({s}) = lambda({s + 2})/2  {text}")
        print(f"  pi^{power}/{denom:<5}          {mp.nstr(right, 50)}")
        gap = mp.nstr(abs(left - right), 3)
        print(f"  difference {gap}  page prefix {prefix}: {text.startswith(prefix)}")
    print()

    print("SERIES SUMMED FROM COUNTED FILLS, NO CLOSED FORM")
    counted = [fluctuation(DIAGONAL, n) for n in range(1, LIMITS[-1] + 1)]
    for s in (2, 4):
        total = mpf(0)
        for n, value in enumerate(counted, start=1):
            if value:
                total += mpf(value.numerator) / value.denominator / mpf(n) ** s
            if n in LIMITS:
                gap = abs(total - targets[s])
                print(
                    f"  s={s}  n <= {n:4d}  partial {mp.nstr(total, 16)}  gap {mp.nstr(gap, 2)}"
                )
    print()


def shell_counts(nmax):
    bound = math.isqrt(nmax)
    axis = np.arange(-bound, bound + 1)
    norm = axis[:, None] ** 2 + axis[None, :] ** 2
    cls = (np.abs(axis)[:, None] % 2) * 2 + (np.abs(axis)[None, :] % 2)
    keep = norm <= nmax
    return [np.bincount(norm[keep & (cls == c)], minlength=nmax + 1) for c in range(4)]


def chi4_excess(nmax):
    out = np.zeros(nmax + 1, dtype=np.int64)
    for d in range(1, nmax + 1, 2):
        out[d::d] += 1 if d % 4 == 1 else -1
    return out


def q_poly(a, t):
    return a[0] * t * t + a[3] * t * (1 - t) + (a[1] + a[2]) * (1 - t) / 2


def indicators(code):
    return [(code >> k) & 1 for k in range(4)]


def lattice_class_sums(radius):
    axis = np.arange(-radius, radius + 1, dtype=np.float64)
    limit = float(radius) ** 2
    out = [0.0] * 4
    for ri in range(2):
        for rj in range(2):
            rows = axis[np.abs(axis) % 2 == ri]
            cols = axis[np.abs(axis) % 2 == rj]
            squares = cols * cols
            pieces = []
            for i in rows:
                norm = i * i + squares
                norm = norm[(norm <= limit) & (norm > 0)]
                pieces.append(float(np.sum(1.0 / (norm * norm))))
            out[ri * 2 + rj] = math.fsum(pieces)
    return out


def gaussian_part():
    print(f"SHELL IDENTITIES ON Z^2, n <= {SHELL_MAX}")
    ee, eo, oe, oo = shell_counts(SHELL_MAX)
    r2 = ee + eo + oe + oo
    print(
        f"  r2(n) = 4(d1(n) - d3(n)): {bool(np.all(r2[1:] == 4 * chi4_excess(SHELL_MAX)[1:]))}"
    )
    want = np.zeros(SHELL_MAX + 1, dtype=np.int64)
    want[4::4] = r2[1 : SHELL_MAX // 4 + 1]
    print(f"  S_ee(4m) = r2(m), zero elsewhere: {bool(np.all(ee[1:] == want[1:]))}")
    want = np.zeros(SHELL_MAX + 1, dtype=np.int64)
    want[2::2] = (eo + oe)[1 : SHELL_MAX // 2 + 1]
    print(f"  S_oo(2m) = S_mix(m), zero elsewhere: {bool(np.all(oo[1:] == want[1:]))}")
    print(f"  S_eo = S_oe: {bool(np.all(eo == oe))}")
    print(
        "  so S_mix = (1-t) S and Q_c(t) = a_ee t^2 + a_oo t(1-t) + (a_eo + a_oe)(1-t)/2"
    )
    print()

    print("THE FIFTEEN NONEMPTY DESIGNS AT s = 2, Z_c(2) / (pi^2 G)")
    mp.dps = 50
    for code in range(1, 16):
        value = q_poly(indicators(code), Fraction(1, 4)) * Fraction(2, 3)
        ok = str(value) == STATED[code]
        print(f"  code {code:2d}  {value!s:>5}  page {STATED[code]:>5}  {ok}")
    z7 = mpf(13) / 24 * mp.pi**2 * catalan
    print(
        f"  code 7 (Sierpinski) 13 pi^2 G / 24 = {mp.nstr(z7, 20)}  page 4.8967847822"
    )
    print()

    print(f"TRUNCATED LATTICE SUMS, RADIUS {RADIUS}, TAIL pi/(4 R^2) PER CLASS")
    sums = lattice_class_sums(RADIUS)
    tail = math.pi / (4.0 * RADIUS * RADIUS)
    base = 4.0 * (math.pi**2 / 6.0) * float(catalan)
    worst = 0.0
    for code in range(1, 16):
        a = indicators(code)
        measured = sum(sums[k] for k in range(4) if a[k]) + tail * sum(a)
        predicted = base * float(q_poly(a, Fraction(1, 4)))
        gap = abs(measured - predicted)
        worst = max(worst, gap)
        print(
            f"  code {code:2d}  predicted {predicted:.12f}  measured {measured:.12f}  gap {gap:.1e}"
        )
    print(f"  worst gap {worst:.1e}")
    print()


def main():
    start = time.time()
    lambda_part()
    gaussian_part()
    print(
        f"DOMAIN  fluctuation n <= {EXACT_TO}, series n <= {LIMITS[-1]}, "
        f"shells n <= {SHELL_MAX}, radius {RADIUS}"
    )
    print(f"WALL {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
