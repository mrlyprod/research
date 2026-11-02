from fractions import Fraction
from math import cos, gcd, log, pi, sqrt

import numpy as np

C = 2.0 * sqrt(2.0) / pi


def sine_coefficient(n):
    return C / n if n % 2 == 1 else 0.0


def dilate_coefficient(n, k):
    total = 0.0
    for j in range(n):
        total += (
            ((-1) ** j) * (cos(k * pi * j / n) - cos(k * pi * (j + 1) / n)) / (k * pi)
        )
    return sqrt(2.0) * total


def check_dilation_shift(scales, harmonics):
    worst = 0.0
    for n in range(1, scales + 1):
        for k in range(1, harmonics + 1):
            want = sine_coefficient(k // n) if k % n == 0 else 0.0
            worst = max(worst, abs(dilate_coefficient(n, k) - want))
    print(
        "shift: <s(nx), e_k> equals a_{k/n} when n divides k and 0 otherwise, max deviation "
        "%.2e over n <= %d and k <= %d, so the dilate s(nx) carries the Dirichlet series n^-s S(s)"
        % (worst, scales, harmonics)
    )


def grid_integral(m, n):
    lcm = m * n // gcd(m, n)
    total = 0
    for k in range(lcm):
        a = (m * (2 * k + 1)) // (2 * lcm)
        b = (n * (2 * k + 1)) // (2 * lcm)
        total += 1 if (a + b) % 2 == 0 else -1
    return Fraction(total, lcm)


def symbol_entry(m, n):
    g = gcd(m, n)
    if (m // g) % 2 == 1 and (n // g) % 2 == 1:
        return Fraction(g * g, m * n)
    return Fraction(0)


def symbol_series(m, n, terms):
    g = gcd(m, n)
    mm, nn = m // g, n // g
    total = 0.0
    for j in range(1, terms + 1):
        total += sine_coefficient(j * nn) * sine_coefficient(j * mm)
    return total


def check_gram_two_ways(cap, terms):
    pairs = 0
    worst = 0.0
    for m in range(1, cap + 1):
        for n in range(m, cap + 1):
            exact = grid_integral(m, n)
            assert exact == symbol_entry(m, n), "gram entry (%d,%d)" % (m, n)
            g = gcd(m, n)
            slack = C * C / ((m // g) * (n // g) * 2.0 * terms)
            gap = abs(symbol_series(m, n, terms) - float(exact))
            assert gap <= slack + 1e-12, "symbol series (%d,%d)" % (m, n)
            worst = max(worst, gap)
            pairs += 1
    print(
        "gram: the L^2(0,1) integral of s(mx)s(nx) on the lcm grid equals the symbol sum "
        "sum_j a_{jn'} a_{jm'} at all %d pairs m <= n <= %d, exact in rationals, reading %s at "
        "(3,5), %s at (3,9) and %s at (2,3)"
        % (pairs, cap, grid_integral(3, 5), grid_integral(3, 9), grid_integral(2, 3))
    )
    print(
        "gram: the truncated symbol sum at %d terms sits within its own tail bound of the closed "
        "form at every pair, largest gap %.2e" % (terms, worst)
    )


def check_blocks(cap):
    off = 0
    for m in range(1, cap + 1):
        for n in range(1, cap + 1):
            e = symbol_entry(m, n)
            am, an = (m & -m).bit_length() - 1, (n & -n).bit_length() - 1
            if am != an:
                assert e == 0, "block leak (%d,%d)" % (m, n)
                off += 1
            else:
                mo, no = m >> am, n >> an
                assert e == Fraction(gcd(mo, no) ** 2, mo * no), (
                    "block value (%d,%d)" % (m, n)
                )
    print(
        "blocks: the Gram entry vanishes whenever v_2(m) differs from v_2(n), %d ordered pairs to "
        "%d, and the block at v_2 = a is the odd Gram itself, so the full dilation system's Gram "
        "is a countable direct sum of copies of the odd one" % (off, cap)
    )


def exact_det(rows):
    a = [row[:] for row in rows]
    n = len(a)
    det = Fraction(1)
    for i in range(n):
        p = next(r for r in range(i, n) if a[r][i] != 0)
        if p != i:
            a[i], a[p] = a[p], a[i]
            det = -det
        det *= a[i][i]
        inv = Fraction(1) / a[i][i]
        for r in range(i + 1, n):
            f = a[r][i] * inv
            if f:
                for c in range(i, n):
                    a[r][c] -= f * a[i][c]
    return det


def smith_product(odds):
    want = Fraction(1)
    for k in odds:
        num, den, t, p = 1, 1, k, 3
        if t % 2 == 0:
            t //= 2
        while p * p <= t:
            if t % p == 0:
                num *= p * p - 1
                den *= p * p
                while t % p == 0:
                    t //= p
            p += 2
        if t > 1:
            num *= t * t - 1
            den *= t * t
        want *= Fraction(num, den)
    return want


def check_determinant(cap):
    for k in range(1, cap + 1):
        odds = list(range(1, 2 * k, 2))
        rows = [[Fraction(gcd(a, b) ** 2, a * b) for b in odds] for a in odds]
        assert exact_det(rows) == smith_product(odds), "det at K=%d" % k
    odds = list(range(1, 26, 2))
    print(
        "det: the odd Gram determinant equals prod over odd k of prod over p | k of (1 - p^-2) at "
        "every K = 1..%d, the value over the %d odd scales n <= 25 being %s, the number "
        "lab/py/stack-levels prints" % (cap, len(odds), smith_product(odds))
    )


def check_inverse(cap, primecap):
    ok = 0
    for n in range(1, cap + 1):
        total = Fraction(0)
        for d in range(1, n + 1):
            if n % d == 0 and d % 2 == 1 and (n // d) % 2 == 1:
                mu, t, p, sq = 1, d, 3, False
                while p * p <= t:
                    if t % p == 0:
                        t //= p
                        if t % p == 0:
                            sq = True
                            break
                        mu = -mu
                    p += 2
                if sq:
                    continue
                if t > 1:
                    mu = -mu
                total += Fraction(mu, d) * Fraction(1, n // d)
        assert total == (1 if n == 1 else 0), "inverse at n=%d" % n
        ok += 1
    sieve = bytearray([1]) * (primecap + 1)
    sieve[0] = sieve[1] = 0
    for p in range(2, int(primecap**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = bytearray(len(sieve[p * p :: p]))
    tail = sum(C / p for p in range(3, primecap + 1, 2) if sieve[p])
    print(
        "inverse: the Dirichlet inverse of a_n/a_1 is mu(n)/n on the odd n and 0 on the even, "
        "checked at every n <= %d, so 1/S(s) = 1/(a_1 (1 - 2^-1-s) zeta(1 + s)) and the Mobius "
        "square wave, coefficients mu(n) a_n, has symbol a_1^2/S(s)" % ok
    )
    print(
        "inverse: a_n/a_1 and its inverse are totally multiplicative, and the Riesz test sum over "
        "the primes reads %.4f already at p <= %d, a divergent sum by Mertens, so neither system "
        "is a Riesz basis; the verdict is invariant under scaling phi by a_1"
        % (tail, primecap)
    )


def spectrum(caps):
    print(
        "spectrum: K, top odd scale N, lambda_max, lambda_min, condition number, "
        "lambda_max/(log N)^2, lambda_max/(log log N)^2"
    )
    for k in caps:
        odds = np.arange(1, 2 * k, 2, dtype=np.float64)
        g = np.gcd(
            odds.astype(np.int64)[:, None], odds.astype(np.int64)[None, :]
        ).astype(np.float64)
        mat = g * g / (odds[:, None] * odds[None, :])
        ev = np.linalg.eigvalsh(mat)
        top, bot = ev[-1], ev[0]
        n = 2 * k - 1
        print(
            "spectrum: %5d %6d %10.5f %12.3e %12.3e %8.4f %8.4f"
            % (k, n, top, bot, top / bot, top / (log(n) ** 2), top / (log(log(n)) ** 2))
        )


def main():
    check_dilation_shift(8, 60)
    check_gram_two_ways(12, 20000)
    check_blocks(64)
    check_determinant(13)
    check_inverse(400, 100000)
    spectrum([25, 50, 100, 150, 200])


if __name__ == "__main__":
    main()
