import cmath
import math
import os
from fractions import Fraction

import numpy as np
from PIL import Image

N_LINEAR = 30
N_QUAD = 60
B_MAX = 12
DRIFTS = [
    Fraction(1, 2),
    Fraction(1, 3),
    Fraction(1, 4),
    Fraction(2, 5),
    Fraction(1, 6),
    Fraction(3, 8),
    Fraction(5, 9),
]
LINEAR_DRIFT = Fraction(2, 7)
GAUSS_PRIMES = [5, 7, 11, 13]
ORIGIN_D = 60
WEYL_NS = [1000, 10000]
PAIR_MAX = 40
SHARE_MAX = 12
SHARE_D = 16
WIDE_B = 20
WIDE_D = 16
WIDE_PAIR = 120
DISTINCT_N = 25
RENDER_N = 30
PANEL = 500
DIGITS = 50

# EXACT IRRATIONAL DRIFTS


def sqrt_fraction(k):
    scale = 10**DIGITS
    return Fraction(math.isqrt(k * scale * scale), scale)


def silver():
    return sqrt_fraction(2) - 1


def golden():
    return (sqrt_fraction(5) - 1) / 2


# THE LEAN


def phase_linear(n, drift):
    return n * drift


def phase_quadratic(n, drift):
    return n * n * drift


def lit(n, x, theta):
    return (n * x - theta).denominator == 1


def brightness_literal(top, x, drift, phase):
    return sum(1 for n in range(1, top + 1) if lit(n, x, phase(n, drift)))


def reduced_points(top_b):
    out = [(0, 1)]
    for b in range(2, top_b + 1):
        for a in range(1, b):
            if math.gcd(a, b) == 1:
                out.append((a, b))
    return out


# THE QUADRATIC CONGRUENCE


def valuation(m, p):
    e = 0
    while m % p == 0:
        m //= p
        e += 1
    return e


def prime_factors(m):
    out = []
    q = m
    p = 2
    while p * p <= q:
        if q % p == 0:
            out.append(p)
            while q % p == 0:
                q //= p
        p += 1
    if q > 1:
        out.append(q)
    return out


def local_classes(p, a, b, c, d):
    beta = valuation(b, p)
    delta = valuation(d, p)
    modulus = p ** max(beta, (delta + 1) // 2)
    if beta >= 1 and beta <= delta < 2 * beta:
        k = 2 * beta - delta
        unit = (a * (d // p**delta) * pow(c * (b // p**beta), -1, p**k)) % p**k
        return modulus, [0, (p ** (delta - beta) * unit) % modulus]
    return modulus, [0]


def crt_merge(state, modulus, classes):
    period, residues = state
    step = period * modulus
    inv = pow(period % modulus, -1, modulus)
    out = []
    for r in residues:
        for s in classes:
            out.append((r + period * ((s - r) * inv % modulus)) % step)
    return step, sorted(out)


def solution_classes(a, b, c, d):
    state = (1, [0])
    for p in prime_factors(b * d):
        state = crt_merge(state, *local_classes(p, a, b, c, d))
    return state


def count_class(top, period, residue):
    if residue == 0:
        return top // period
    return 0 if residue > top else (top - residue) // period + 1


def brightness_form(top, a, b, c, d):
    period, residues = solution_classes(a, b, c, d)
    return sum(count_class(top, period, r) for r in residues)


def first_lit(a, b, c, d):
    period, residues = solution_classes(a, b, c, d)
    return min(period if r == 0 else r for r in residues)


def square_root_kernel(d):
    out = 1
    for p in prime_factors(d):
        out *= p ** ((valuation(d, p) + 1) // 2)
    return out


def middle_primes(b, d):
    return [
        p
        for p in prime_factors(b)
        if valuation(b, p) <= valuation(d, p) < 2 * valuation(b, p)
    ]


def period_form(b, d):
    k = square_root_kernel(d)
    return b * k // math.gcd(b, k)


def halving_rule(b, d):
    period = period_form(b, d)
    beta = valuation(b, 2)
    return period // 2 if beta >= 1 and valuation(d, 2) == 2 * beta - 1 else period


# THE SOLUTION SET READ LITERALLY


def literal_mask(a, b, c, d):
    m = b * d
    n = np.arange(m, dtype=np.int64)
    return (n * (a * d - b * c * n)) % m == 0


def divisors(m):
    out = []
    k = 1
    while k * k <= m:
        if m % k == 0:
            out.append(k)
            out.append(m // k)
        k += 1
    return sorted(set(out))


def minimal_period(mask):
    m = mask.size
    count = int(mask.sum())
    for q in divisors(m):
        if count % (m // q):
            continue
        if np.array_equal(mask, np.roll(mask, q)):
            return q
    return m


def formula_set(a, b, c, d):
    period, residues = solution_classes(a, b, c, d)
    return {r + period * k for r in residues for k in range(b * d // period)}


def shares_point(m, n, drift):
    return (m * n // math.gcd(m, n) * (m - n) * drift).denominator == 1


def layer_points(n, drift):
    theta = (n * n * drift) % 1
    return {(k + theta) / n for k in range(n)}


# GAUSS SUMS


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def phase_census(c, p):
    out = {}
    for n in range(p):
        j = c * n * n % p
        out[j] = out.get(j, 0) + 1
    return out


def gauss_sum(coefficient, m):
    return sum(
        cmath.exp(2j * math.pi * (coefficient * n * n % m) / m) for n in range(m)
    )


def epsilon(p):
    return 1 if p % 4 == 1 else 1j


def quadratic_sum(a_coeff, b_coeff, m):
    return sum(
        cmath.exp(2j * math.pi * ((a_coeff * n * n + b_coeff * n) % m) / m)
        for n in range(m)
    )


# SECTIONS


def linear_lean():
    print("LINEAR LEAN  t_n = delta, phase n delta, N =", N_LINEAR)
    drift = LINEAR_DRIFT
    bad = 0
    for a, b in reduced_points(N_LINEAR):
        x = drift + Fraction(a, b)
        if brightness_literal(N_LINEAR, x, drift, phase_linear) != N_LINEAR // b:
            bad += 1
    off = [drift + Fraction(1, q) for q in (31, 37, 41, 43)]
    blank = sum(
        1 for x in off if brightness_literal(N_LINEAR, x, drift, phase_linear) != 0
    )
    print(
        "  delta =",
        drift,
        " nodes tested",
        len(reduced_points(N_LINEAR)),
        " brightness != floor(N/b):",
        bad,
        " nonzero off the translate:",
        blank,
    )
    for name, value in (
        ("sqrt 2 - 1", math.sqrt(2) - 1),
        ("phi - 1", (math.sqrt(5) - 1) / 2),
    ):
        bad = 0
        for a, b in reduced_points(N_LINEAR):
            x = value + a / b
            hit = sum(
                1
                for n in range(1, N_LINEAR + 1)
                if abs(n * x - n * value - round(n * x - n * value)) < 1e-12
            )
            if hit != N_LINEAR // b:
                bad += 1
        closest = min(
            min(
                abs(n * (a / b) - n * value - round(n * (a / b) - n * value))
                for n in range(1, N_LINEAR + 1)
            )
            for a, b in reduced_points(N_LINEAR)[1:]
        )
        print(
            "  delta = %s  brightness != floor(N/b) at 1e-12: %d  closest approach at a rational x: %.6f"
            % (name, bad, closest)
        )


def quadratic_lean():
    print(
        "QUADRATIC LEAN  t_n = n delta, phase n^2 delta, N =",
        N_QUAD,
        " all reduced a/b with b <=",
        B_MAX,
    )
    points = reduced_points(B_MAX)
    for drift in DRIFTS:
        c, d = drift.numerator, drift.denominator
        bad = sum(
            1
            for a, b in points
            if brightness_literal(N_QUAD, Fraction(a, b), drift, phase_quadratic)
            != brightness_form(N_QUAD, a, b, c, d)
        )
        set_bad = 0
        period_bad = 0
        halved = 0
        for a, b in points:
            mask = literal_mask(a, b, c, d)
            if formula_set(a, b, c, d) != {int(v) for v in np.nonzero(mask)[0]}:
                set_bad += 1
            least = minimal_period(mask)
            if least != halving_rule(b, d):
                period_bad += 1
            if least < period_form(b, d):
                halved += 1
        top = max(brightness_form(N_QUAD, a, b, c, d) for a, b in points)
        best = [
            Fraction(a, b)
            for a, b in points
            if brightness_form(N_QUAD, a, b, c, d) == top
        ]
        origin = brightness_form(N_QUAD, 0, 1, c, d)
        print(
            "  delta = %s  points %d  brightness mismatches %d  solution-set mismatches %d  minimal-period mismatches %d  minimal period below lcm(b, d*) at %d points"
            % (drift, len(points), bad, set_bad, period_bad, halved)
        )
        print(
            "    top brightness %d at %s  origin %d  origin is top: %s"
            % (
                top,
                ", ".join(str(v) for v in best),
                origin,
                "yes" if origin == top else "no",
            )
        )


def minimal_period_law():
    print(
        "THE MINIMAL PERIOD  lcm(b, d*) is a period; the least one halves it exactly when v_2(b) >= 1 and v_2(d) = 2 v_2(b) - 1"
    )
    points = reduced_points(WIDE_B)
    drifts = [
        Fraction(c, e)
        for e in range(1, WIDE_B + 1)
        for c in range(1, e + 1)
        if math.gcd(c, e) == 1
    ]
    tuples = 0
    halved = 0
    bad = 0
    for drift in drifts:
        c, d = drift.numerator, drift.denominator
        for a, b in points:
            tuples += 1
            least = minimal_period(literal_mask(a, b, c, d))
            if least < period_form(b, d):
                halved += 1
            if least != halving_rule(b, d):
                bad += 1
    print(
        "  tuples %d over b, d <= %d  minimal period below lcm(b, d*) at %d  breaches of the rule %d"
        % (tuples, WIDE_B, halved, bad)
    )
    for a, b, c, d in ((1, 2, 1, 2), (1, 6, 1, 2), (1, 4, 3, 8)):
        mask = literal_mask(a, b, c, d)
        print(
            "    x = %d/%d  delta = %d/%d  lcm(b, d*) = %d  classes %s  minimal period %d"
            % (
                a,
                b,
                c,
                d,
                period_form(b, d),
                sorted(solution_classes(a, b, c, d)[1]),
                minimal_period(mask),
            )
        )


def lit_set():
    print("THE LIT SET  first layer lighting a/b, and its dependence on the numerator")
    points = reduced_points(B_MAX)
    for drift in DRIFTS:
        c, d = drift.numerator, drift.denominator
        bad = 0
        for a, b in points:
            reach = first_lit(a, b, c, d)
            literal = next(
                (
                    n
                    for n in range(1, N_QUAD + 1)
                    if lit(n, Fraction(a, b), phase_quadratic(n, drift))
                ),
                None,
            )
            if literal != (reach if reach <= N_QUAD else None):
                bad += 1
        generic = sum(
            1 for a, b in points if first_lit(a, b, c, d) == period_form(b, d)
        )
        print(
            "  delta = %s  first-lit mismatches %d  points first lit at lcm(b, d*) itself %d of %d"
            % (drift, bad, generic, len(points))
        )
    print("  numerator dependence, delta = 1/4:")
    for a in (1, 3):
        period, residues = solution_classes(a, 4, 1, 4)
        print(
            "    x = %d/4  classes mod %d: %s  first lit n = %d  B_61 = %d"
            % (
                a,
                period,
                residues,
                first_lit(a, 4, 1, 4),
                brightness_form(61, a, 4, 1, 4),
            )
        )
    split = [
        (a, b, Fraction(c, d))
        for b in range(2, B_MAX + 1)
        for d in (2, 4, 8, 9)
        for c in range(1, d)
        if math.gcd(c, d) == 1
        for a in range(1, b)
        if math.gcd(a, b) == 1 and first_lit(a, b, c, d) != first_lit(1, b, c, d)
    ]
    print(
        "  reduced (a/b, delta) pairs with b <= %d and d in 2, 4, 8, 9 whose first lit layer moves with the numerator: %d"
        % (B_MAX, len(split))
    )


def zeta_em(s, terms=400000):
    k = np.arange(1, terms + 1, dtype=np.float64)
    total = float(np.sum(k ** (-float(s))))
    m = float(terms)
    return total + m ** (1 - s) / (s - 1) - 0.5 * m ** (-s) + s * m ** (-s - 1) / 12.0


def inner_roots(limit):
    inner = np.ones(limit + 1, dtype=np.int64)
    k = 2
    while k * k <= limit:
        inner[k * k :: k * k] = k
        k += 1
    return inner


def origin_law():
    print("THE ORIGIN  x = 0, brightness floor(N/d*) with d* the least k with d | k^2")
    bad = 0
    for d in range(1, ORIGIN_D + 1):
        for c in range(1, d + 1):
            if math.gcd(c, d) != 1:
                continue
            drift = Fraction(c, d)
            if brightness_literal(
                N_QUAD, Fraction(0), drift, phase_quadratic
            ) != N_QUAD // square_root_kernel(d):
                bad += 1
    kernels = [square_root_kernel(d) for d in range(1, 31)]
    print("  drifts tested to d =", ORIGIN_D, " mismatches", bad)
    print("  d* for d = 1..30:", ", ".join(str(v) for v in kernels))
    limit = 10**6
    d = np.arange(1, limit + 1, dtype=np.float64)
    star = d / inner_roots(limit)[1:].astype(np.float64)
    breaches = sum(
        1 for k in range(1, 1001) if int(star[k - 1]) != square_root_kernel(k)
    )
    print("  d/A000188(d) against the local form to d = 1000: breaches", breaches)
    for s in (1, 2):
        partials = [
            float(np.sum(1.0 / (star[:cut] * d[:cut] ** s)))
            for cut in (10**4, 10**5, 10**6)
        ]
        target = zeta_em(2 * s + 1) * zeta_em(s + 1) / zeta_em(2 * s + 2)
        print(
            "  s = %d  partial sums %s  zeta(%d) zeta(%d)/zeta(%d) = %.6f"
            % (
                s,
                ", ".join("%.6f" % v for v in partials),
                2 * s + 1,
                s + 1,
                2 * s + 2,
                target,
            )
        )


def gauss_sums():
    print(
        "GAUSS SUMS  the phase census of the quadratic lean at prime drift denominator, and the centred twist S(c,p) = (c/p) S(1,p)"
    )
    for p in GAUSS_PRIMES:
        residue = next(c for c in range(1, p) if legendre(c, p) == 1)
        nonresidue = next(c for c in range(1, p) if legendre(c, p) == -1)
        for c in (residue, nonresidue):
            census = phase_census(c, p)
            inverse = pow(c, -1, p)
            bad = sum(
                1
                for j in range(p)
                if census.get(j, 0) != (1 if j == 0 else 1 + legendre(j * inverse, p))
            )
            value = gauss_sum(c, p)
            claim = legendre(c, p) * epsilon(p) * math.sqrt(p)
            twist = sum(
                1
                for j in range(p)
                if census.get(j, 0) - 1
                != legendre(c, p) * (phase_census(1, p).get(j, 0) - 1)
            )
            print(
                "  p = %2d  c = %2d  (c/p) = %2d  distinct phases %d  census breaches %d  centred twist breaches %d  S(c,p) = %.6f%+.6fi  (c/p) eps_p sqrt p = %.6f%+.6fi"
                % (
                    p,
                    c,
                    legendre(c, p),
                    len(census),
                    bad,
                    twist,
                    value.real,
                    value.imag,
                    claim.real,
                    claim.imag,
                )
            )
    print(
        "  the count of solutions as a Fourier sum, R = (1/m) sum_h sum_n e(h Q(n)/m)"
    )
    for a, b, c, d in (
        (1, 1, 1, 5),
        (1, 2, 1, 3),
        (1, 4, 1, 4),
        (1, 6, 1, 6),
        (5, 12, 3, 8),
    ):
        m = b * d
        total = sum(quadratic_sum(h * b * c, -h * a * d, m) for h in range(m)) / m
        period, residues = solution_classes(a, b, c, d)
        exact = len(residues) * (m // period)
        print(
            "  a/b = %d/%d  delta = %d/%d  m = %2d  Fourier count %.6f%+.6fi  classes mod %d %s  R = %d"
            % (a, b, c, d, m, total.real, total.imag, period, residues, exact)
        )


def irrational_lean():
    print(
        "IRRATIONAL DRIFT  no two layers share a point, brightness at most 1 everywhere"
    )
    for name, drift in (("sqrt 2 - 1", silver()), ("phi - 1", golden())):
        worst = None
        for n in range(1, PAIR_MAX + 1):
            for m in range(1, n):
                span = m * n // math.gcd(m, n)
                t = drift * (m - n) * span
                frac = t - math.floor(t)
                gap = Fraction(1, span) * min(frac, 1 - frac)
                if worst is None or gap < worst[0]:
                    worst = (gap, m, n)
        print(
            "  delta = %s  closest two layers to N = %d come: %.3e at (m, n) = (%d, %d)"
            % (name, PAIR_MAX, float(worst[0]), worst[1], worst[2])
        )
        row = []
        for top in WEYL_NS:
            q = drift.denominator
            p = drift.numerator
            values = sorted((n * n * p) % q for n in range(1, top + 1))
            star = 0.0
            for i, r in enumerate(values, start=1):
                u = r / q
                star = max(star, i / top - u, u - (i - 1) / top)
            row.append(
                "N = %5d  D*_N = %.6f  1/sqrt N = %.6f  ratio %.4f"
                % (top, star, 1 / math.sqrt(top), star * math.sqrt(top))
            )
        for line in row:
            print("   ", line)


# ADVERSARIAL PASS


def sharing_law():
    print(
        "SHARING  layers m and n share a lit point iff lcm(m, n)(m - n) delta is an integer"
    )
    drifts = [
        Fraction(c, d)
        for d in range(1, SHARE_D + 1)
        for c in range(1, d + 1)
        if math.gcd(c, d) == 1
    ]
    pairs = 0
    sharing = 0
    bad = 0
    loose = 0
    slack = 0
    witness = None
    for drift in drifts:
        for n in range(2, SHARE_MAX + 1):
            for m in range(1, n):
                pairs += 1
                truth = len(layer_points(m, drift) & layer_points(n, drift)) > 0
                if truth:
                    sharing += 1
                if truth != shares_point(m, n, drift):
                    bad += 1
                if not truth and (Fraction(n) * (m - n) * drift).denominator == 1:
                    loose += 1
                if not truth and (Fraction(m * n) * (m - n) * drift).denominator == 1:
                    slack += 1
                    if witness is None:
                        witness = (m, n, drift)
    print(
        "  pairs %d over m < n <= %d and every reduced c/d with d <= %d  sharing %d  criterion failures %d"
        % (pairs, SHARE_MAX, SHARE_D, sharing, bad)
    )
    print(
        "  the rational reading is vacuous at rational drift, every one of the %d non-sharing pairs having n (m - n) delta rational; n (m - n) delta in Z holds at %d of them and m n (m - n) delta in Z at %d, first at (m, n, delta) = (%d, %d, %s), so both are weaker than the criterion"
        % (pairs - sharing, loose, slack, witness[0], witness[1], witness[2])
    )


def adversarial():
    print(
        "ADVERSARIAL  the closed form against literal stacking on the widest domain this study runs"
    )
    points = reduced_points(WIDE_B)
    drifts = [
        Fraction(c, d)
        for d in range(1, WIDE_D + 1)
        for c in range(1, d + 1)
        if math.gcd(c, d) == 1
    ]
    bad = 0
    tested = 0
    for drift in drifts:
        c, d = drift.numerator, drift.denominator
        for a, b in points:
            tested += 1
            if brightness_literal(
                N_QUAD, Fraction(a, b), drift, phase_quadratic
            ) != brightness_form(N_QUAD, a, b, c, d):
                bad += 1
    print(
        "  points %d  drifts %d  pairs %d  brightness mismatches at N = %d: %d"
        % (len(points), len(drifts), tested, N_QUAD, bad)
    )
    for name, drift in (("sqrt 2 - 1", silver()), ("phi - 1", golden())):
        seen = set()
        for n in range(1, DISTINCT_N + 1):
            theta = (n * n * drift) % 1
            for k in range(n):
                seen.add((k + theta) / n)
        worst = None
        for n in range(1, WIDE_PAIR + 1):
            for m in range(1, n):
                span = m * n // math.gcd(m, n)
                t = drift * (m - n) * span
                frac = t - math.floor(t)
                gap = Fraction(1, span) * min(frac, 1 - frac)
                if worst is None or gap < worst[0]:
                    worst = (gap, m, n)
        print(
            "  delta = %s  lit points of layers 1..%d: %d distinct of %d drawn  closest pair to N = %d: %.3e at (%d, %d)"
            % (
                name,
                DISTINCT_N,
                len(seen),
                DISTINCT_N * (DISTINCT_N + 1) // 2,
                WIDE_PAIR,
                float(worst[0]),
                worst[1],
                worst[2],
            )
        )


# RENDER


def lean_axis(top, drift, resolution):
    counts = np.zeros(resolution, dtype=np.int64)
    for n in range(1, top + 1):
        theta = float((n * n * drift) % 1)
        for k in range(n):
            counts[int((k + theta) / n * resolution) % resolution] += 1
    return counts


def panel_grey(drift):
    axis = lean_axis(RENDER_N, drift, PANEL)
    field = np.outer(axis, axis)
    levels = np.unique(field)
    rank = np.searchsorted(levels, field).astype(np.float64)
    grey = np.full(field.shape, 255.0)
    lit_mask = field > 0
    grey[lit_mask] = 190 - 190 * (rank[lit_mask] - 1) / max(levels.size - 2, 1)
    return grey.astype(np.uint8), int(field.max())


def render():
    left, peak_left = panel_grey(Fraction(1, 5))
    right, peak_right = panel_grey(silver())
    sheet = np.full((PANEL, 2 * PANEL + 8), 255, dtype=np.uint8)
    sheet[:, :PANEL] = left
    sheet[:, PANEL + 8 :] = right
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "leaning-stack.png")
    Image.fromarray(sheet, mode="L").save(path, optimize=True)
    print(
        "RENDER  leaning-stack.png  N = %d  panels delta = 1/5 and delta = sqrt 2 - 1  peaks %d and %d  %d bytes"
        % (RENDER_N, peak_left, peak_right, os.path.getsize(path))
    )


def main():
    linear_lean()
    quadratic_lean()
    minimal_period_law()
    lit_set()
    origin_law()
    gauss_sums()
    sharing_law()
    adversarial()
    irrational_lean()
    render()


main()
