from fractions import Fraction
from math import gcd, isqrt, lcm

from sympy import Catalan as SCatalan
from sympy import N as N_
from sympy import (
    Poly,
    Rational,
    Symbol,
    cos,
    cyclotomic_poly,
    minimal_polynomial,
    sin,
    sympify,
)
from sympy import pi as spi
from sympy import sqrt as ssqrt
from sympy import zeta as szeta

pi = spi

X = Symbol("x")
CYC = 360
PHI360 = Poly(cyclotomic_poly(CYC, X), X, domain="QQ")
SCHEDULE_N = 30
NEAR_BOX = 6
STACK_N = 50
COUNT_NS = [50, 200, 800]
SPOT_DEGREES = [0, 7, 30, 45, 60, 90, 97, 113, 143, 180, 270, 271]
BASE_BOX = 10
ODD_N = 55
INCREMENTS = [
    (Fraction(0), "0"),
    (Fraction(18), "18"),
    (Fraction(45, 2), "22.5"),
    (Fraction(30), "30"),
    (Fraction(45), "45"),
    (Fraction(60), "60"),
    (Fraction(135, 2), "67.5"),
    (Fraction(36), "36"),
    (Fraction(15), "15"),
    (Fraction(45, 4), "11.25"),
    (Fraction(10), "10"),
    (90 * (ssqrt(2) - 1), "90 (sqrt2 - 1)"),
]
BASE_DEPTH = 8


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gdivides(w, z):
    p = gmul(z, gconj(w))
    n = norm(w)
    return p[0] % n == 0 and p[1] % n == 0


def gquo(z, w):
    p = gmul(z, gconj(w))
    n = norm(w)
    return (p[0] // n, p[1] // n)


def gnearest(z, w):
    p = gmul(z, gconj(w))
    n = norm(w)
    return ((2 * p[0] + n) // (2 * n), (2 * p[1] + n) // (2 * n))


def ggcd(z, w):
    while w != (0, 0):
        q = gnearest(z, w)
        z, w = w, (z[0] - gmul(q, w)[0], z[1] - gmul(q, w)[1])
    return z


def canon(z):
    if z == (0, 0):
        return z
    for c in (z, (-z[1], z[0]), (-z[0], -z[1]), (z[1], -z[0])):
        if c[0] > 0 and c[1] >= 0:
            return c
    return z


def classes_up_to(n):
    out = []
    for a in range(1, isqrt(n) + 1):
        for b in range(isqrt(n - a * a) + 1):
            out.append((a, b))
    return sorted(out, key=lambda z: (norm(z), z))


def circle_classes_jacobi(t):
    s = 0
    j = 0
    while 4 * j + 1 <= t:
        s += t // (4 * j + 1) - t // (4 * j + 3)
        j += 1
    return s


def circle_classes_direct(t):
    return len(classes_up_to(t))


def gauss_primes_dividing(z):
    n = norm(z)
    out = []
    m = n
    q = 2
    while q * q <= m:
        if m % q == 0:
            while m % q == 0:
                m //= q
            out.append(q)
        q += 1
    if m > 1:
        out.append(m)
    primes = []
    for q in out:
        if q == 2:
            primes.append((1, 1))
        elif q % 4 == 3:
            primes.append((q, 0))
        else:
            r = 2
            a = None
            while a is None:
                if pow(r, (q - 1) // 2, q) == q - 1:
                    a = pow(r, (q - 1) // 4, q)
                r += 1
            p = ggcd((q, 0), (a, 1))
            primes.append(canon(p))
            primes.append(canon(gconj(p)))
    return [p for p in primes if gdivides(p, z)]


def gauss_totient(z):
    v = norm(z)
    if v == 1:
        return 1
    for p in gauss_primes_dividing(z):
        v = v // norm(p) * (norm(p) - 1)
    return v


def totient_sum(n_max):
    return sum(gauss_totient(z) for z in classes_up_to(n_max))


def cyc_reduce(terms):
    acc = {}
    for e, c in terms:
        k = (e % CYC,)
        acc[k] = acc.get(k, 0) + c
    acc = {k: v for k, v in acc.items() if v}
    if not acc:
        return Poly(0, X, domain="QQ")
    return Poly(acc, X, domain="QQ").rem(PHI360)


RATIONAL_CACHE = {}


def rational_cos_sin(d):
    d = d % CYC
    if d not in RATIONAL_CACHE:
        c = cyc_reduce([(d, 1), (-d, 1)])
        s = cyc_reduce([(90 + d, -1), (90 - d, 1)])
        RATIONAL_CACHE[d] = (c.degree() <= 0, s.degree() <= 0)
    return RATIONAL_CACHE[d]


def rational_angle_degrees():
    return [d for d in range(CYC) if all(rational_cos_sin(d))]


def spot_check_degrees():
    out = []
    for d in SPOT_DEGREES:
        c = minimal_polynomial(cos(pi * Rational(d, 180)), X, polys=True).degree()
        s = minimal_polynomial(sin(pi * Rational(d, 180)), X, polys=True).degree()
        out.append((d, c == 1, s == 1, rational_cos_sin(d)))
    return out


def primes(k):
    out = []
    n = 2
    while len(out) < k:
        if all(n % p for p in out if p * p <= n):
            out.append(n)
        n += 1
    return out


def prime_schedule(k):
    return list(zip(range(1, k + 1), [0] + primes(k - 1)))


def dead_spin_pairs(schedule):
    hits = []
    for i in range(len(schedule)):
        for j in range(i + 1, len(schedule)):
            m, a = schedule[i]
            n, b = schedule[j]
            if all(rational_cos_sin((a - b) % CYC)):
                hits.append((m, a, n, b))
    return hits


def shared_witness(m, a, n, b):
    g = gcd(m, n)
    q = ((a - b) % CYC) // 90
    v = [(1, 0), (0, 1), (-1, 0), (0, -1)][q]
    return (Fraction(v[0], g), Fraction(v[1], g))


def rot(d):
    from math import cos as fc
    from math import radians
    from math import sin as fs

    return fc(radians(d)), fs(radians(d))


def near_miss(schedule, box):
    best = None
    for i in range(len(schedule)):
        for j in range(i + 1, len(schedule)):
            m, a = schedule[i]
            n, b = schedule[j]
            if all(rational_cos_sin((a - b) % CYC)):
                continue
            c, s = rot(a - b)
            for u in range(-box, box + 1):
                for v in range(-box, box + 1):
                    if u == 0 and v == 0:
                        continue
                    x = (c * u - s * v) * n / m
                    y = (s * u + c * v) * n / m
                    dx = abs(x - round(x))
                    dy = abs(y - round(y))
                    d = max(dx, dy)
                    if best is None or d < best[0]:
                        best = (d, m, a, n, b, u, v)
    return best


def layer_nodes(z):
    a, b = z
    n = norm(z)
    out = set()
    span = a + b + 1
    for m in range(-span, span + 1):
        for k in range(-span, span + 1):
            x = Fraction(m * a + k * b, n)
            y = Fraction(k * a - m * b, n)
            out.add((x - int(x // 1), y - int(y // 1)))
    return out


def literal_stack(n_max):
    hits = {}
    for z in classes_up_to(n_max):
        for p in layer_nodes(z):
            hits[p] = hits.get(p, 0) + 1
    return hits


def reduced_denominator(node):
    x, y = node
    d = lcm(x.denominator, y.denominator)
    u = (x.numerator * (d // x.denominator), y.numerator * (d // y.denominator))
    g = ggcd(u, (d, 0))
    return canon(gquo((d, 0), g))


def closed_brightness(n_max, d):
    return circle_classes_jacobi(n_max // norm(d))


def pythagorean_hits(bound, reach):
    seen = set()
    for a in range(-reach, reach + 1):
        for b in range(-reach, reach + 1):
            if a == 0 and b == 0:
                continue
            n = a * a + b * b
            seen.add((Fraction(a * a - b * b, n), Fraction(2 * a * b, n)))
    direct = set()
    for r in range(1, bound + 1):
        for p in range(-r, r + 1):
            q2 = r * r - p * p
            s = isqrt(q2)
            if s * s == q2:
                for t in (s, -s):
                    direct.add((Fraction(p, r), Fraction(t, r)))
    return len(direct), direct <= seen


def base_c_overlap(c, box):
    best = None
    count = 0
    for m in range(-box, box + 1):
        for k in range(-box, box + 1):
            if m == 0 and k == 0:
                continue
            x = c[0] * m - c[1] * k
            y = c[1] * m + c[0] * k
            d = max(abs(x - round(x)), abs(y - round(y)))
            if d < 1e-12:
                count += 1
            if best is None or d < best:
                best = d
    return count, best


def base_layer(c, k):
    z = (1, 0)
    for _ in range(k):
        z = gmul(z, c)
    return layer_nodes(canon(z))


def unit_square_shares(g, a):
    c, s = rot(a)
    count = 0
    for k in range(-4 * g - 4, 4 * g + 5):
        for l in range(-4 * g - 4, 4 * g + 5):
            if k == 0 and l == 0:
                continue
            x = (c * k - s * l) / g
            y = (s * k + c * l) / g
            if 1e-9 < x < 1 - 1e-9 and 1e-9 < y < 1 - 1e-9:
                count += 1
    return count


def is_whole_turn(ratio, d):
    v = d * ratio
    if isinstance(v, Fraction):
        return v.denominator == 1
    return bool(sympify(v).is_integer)


def increment_period(ratio, layers):
    for q in range(1, layers + 1):
        if is_whole_turn(ratio, q):
            return q
    return None


def increment_classes(theta_deg, layers):
    ratio = theta_deg / 90
    q = increment_period(ratio, layers)
    if q is None:
        sizes = [1] * layers
    else:
        sizes = [len(range(r, layers, q)) for r in range(q)]
    closed = sum(c * (c - 1) // 2 for c in sizes)
    pairwise = sum(
        1
        for j in range(layers)
        for k in range(j + 1, layers)
        if is_whole_turn(ratio, j - k)
    )
    niven = None
    if isinstance(theta_deg, Fraction) and theta_deg.denominator == 1:
        niven = sum(
            1
            for j in range(layers)
            for k in range(j + 1, layers)
            if all(rational_cos_sin(int(theta_deg) * (j - k)))
        )
    return q, sizes, closed, pairwise, niven


def share_count_spread(g):
    return sorted({unit_square_shares(g, a) for a in range(1, 90)})


def base_depth_check(c, depth):
    layers = [base_layer(c, k) for k in range(depth + 1)]
    nested = all(layers[k] <= layers[k + 1] for k in range(depth))
    bad = 0
    for p in layers[depth]:
        d = min(k for k in range(depth + 1) if p in layers[k])
        b = sum(1 for k in range(depth + 1) if p in layers[k])
        if b != depth + 1 - d:
            bad += 1
    return nested, len(layers[depth]), bad


def main():
    rats = rational_angle_degrees()
    print("rational rotation degrees", *rats)
    print("rational rotation count", len(rats))
    spots = spot_check_degrees()
    print("spot degrees", len(spots))
    print(
        "spot minpoly agrees with cyclotomic", all((c, s) == r for _, c, s, r in spots)
    )

    hits, covered = pythagorean_hits(60, 12)
    print("rational unit-circle points denominator <= 60", hits)
    print("all are w^2/N(w), Gaussian w in box 12", covered)

    sched = prime_schedule(SCHEDULE_N)
    print("schedule layers", len(sched))
    print("schedule angles", *[a for _, a in sched])
    pairs = dead_spin_pairs(sched)
    print("schedule pairs", len(sched) * (len(sched) - 1) // 2)
    print("sharing pairs", len(pairs))
    for m, a, n, b in pairs:
        g = gcd(m, n)
        ok = (a - b) % 90 == 0 and n % g == 0 and m % g == 0
        w = shared_witness(m, a, n, b)
        print(
            "share",
            m,
            a,
            n,
            b,
            "gcd",
            g,
            "exact",
            ok,
            "witness",
            str(w[0]),
            str(w[1]),
            "open square nodes",
            unit_square_shares(g, a),
        )
    print(
        "all sharing pairs congruent mod 90",
        all((a - b) % 90 == 0 for _, a, _, b in pairs),
    )
    print("shared node window", "open unit square, origin excluded")
    print("shared lattice density per unit area", "g^2")
    for g in (2, 3, 7):
        print("share counts over whole degrees 1..89 at g", g, *share_count_spread(g))
    odds = len(range(1, ODD_N + 1, 2))
    print("increment schedule odd scales 1 to", ODD_N, "layers", odds)
    for theta, name in INCREMENTS:
        q, sizes, closed, pairwise, niven = increment_classes(theta, odds)
        print(
            "increment",
            name,
            "classes",
            "none" if q is None else q,
            "sizes",
            *(["all 1"] if q is None else sizes),
            "pairs",
            closed,
            "pairwise equal",
            closed == pairwise,
            "niven equal",
            "skipped" if niven is None else closed == niven,
        )
    d, m, a, n, b, u, v = near_miss(sched, NEAR_BOX)
    print("near miss box", NEAR_BOX)
    print("near miss min distance", "%.6f" % d, "at", m, a, n, b, u, v)

    stack = literal_stack(STACK_N)
    print("stack norm bound", STACK_N)
    print("stack layers", len(classes_up_to(STACK_N)))
    print("stack nodes", len(stack))
    print("gaussian totient sum", totient_sum(STACK_N))
    print("nodes equal totient sum", len(stack) == totient_sum(STACK_N))
    bad = 0
    for p, b in stack.items():
        if closed_brightness(STACK_N, reduced_denominator(p)) != b:
            bad += 1
    print("brightness comparisons", len(stack))
    print("brightness mismatches", bad)
    print(
        "brightest node", max(stack.values()), "at 0", stack[(Fraction(0), Fraction(0))]
    )
    print(
        "circle count jacobi equals direct",
        all(circle_classes_jacobi(t) == circle_classes_direct(t) for t in range(401)),
    )
    print("circle counts t=1..12", *[circle_classes_jacobi(t) for t in range(1, 13)])
    for n in COUNT_NS:
        s = totient_sum(n)
        print("node count N", n, "=", s, "ratio N^2", "%.6f" % (s / (n * n)))
    print(
        "pi / (8 zeta(2) Catalan)", "%.6f" % float(N_(spi / (8 * szeta(2) * SCatalan)))
    )

    for c, name, depth in (((1, 1), "1+i", BASE_DEPTH), ((2, 1), "2+i", 4)):
        nested, size, bad = base_depth_check(c, depth)
        print(
            "base",
            name,
            "depth",
            depth,
            "nested",
            nested,
            "deepest layer nodes",
            size,
            "address mismatches",
            bad,
        )
    for c, name in (
        ((1.0, 1.0), "1+i"),
        ((1.5, 0.5), "3/2+i/2"),
        (
            (
                2**0.5 * __import__("math").cos(1.0),
                2**0.5 * __import__("math").sin(1.0),
            ),
            "sqrt2 e^i",
        ),
    ):
        count, best = base_c_overlap(c, BASE_BOX)
        print(
            "base",
            name,
            "box",
            BASE_BOX,
            "overlaps",
            count,
            "min distance",
            "%.6g" % best,
        )


main()
