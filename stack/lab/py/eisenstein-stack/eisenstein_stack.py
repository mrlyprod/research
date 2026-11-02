from fractions import Fraction
from math import isqrt, lcm, log
from pathlib import Path

from sympy import N as N_
from sympy import (
    Poly,
    Rational,
    Symbol,
    cos,
    cyclotomic_poly,
    minimal_polynomial,
    sin,
    sqrt,
)
from sympy import pi as spi
from sympy import zeta as szeta

HERE = Path(__file__).resolve().parent

X = Symbol("x")
CYC = 360
PHI360 = Poly(cyclotomic_poly(CYC, X), X, domain="QQ")
SPOT_DEGREES = [0, 30, 45, 60, 90, 120, 137, 180, 240, 270, 300, 359]
STACK_N = 50
COUNT_NS = [50, 200, 800, 1200, 2000, 3200]
CIRCLE_T = 400
CSL_M = 100
ROT_BOUND = 60
ROT_REACH = 20
FIG_N = 60
FIG_W = 900
FIG_H = 600
FIG_OUT = HERE / "eisenstein-stack.png"
ROOT3 = 3.0**0.5

SURF = (255, 255, 255)
BLUE = (0, 140, 255)
GRAY = (186, 186, 191)
INK = (0, 0, 0)


def enorm(z):
    a, b = z
    return a * a - a * b + b * b


def emul(z, w):
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c - b * d)


def econj(z):
    a, b = z
    return (a - b, -b)


def edivides(w, z):
    p = emul(z, econj(w))
    n = enorm(w)
    return p[0] % n == 0 and p[1] % n == 0


def equo(z, w):
    p = emul(z, econj(w))
    n = enorm(w)
    return (p[0] // n, p[1] // n)


def enearest(z, w):
    p = emul(z, econj(w))
    n = enorm(w)
    return ((2 * p[0] + n) // (2 * n), (2 * p[1] + n) // (2 * n))


def egcd(z, w):
    while w != (0, 0):
        q = enearest(z, w)
        t = emul(q, w)
        z, w = w, (z[0] - t[0], z[1] - t[1])
    return z


def associates(z):
    out = [z]
    for _ in range(2):
        a, b = out[-1]
        out.append((-b, a - b))
    return out + [(-a, -b) for a, b in out]


def ecanon(z):
    if z == (0, 0):
        return z
    for c in associates(z):
        if c[0] > 0 and 0 <= c[1] < c[0]:
            return c
    return z


def classes_up_to(n):
    out = set()
    r = isqrt(4 * n // 3) + 2
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            if (a, b) == (0, 0):
                continue
            if enorm((a, b)) <= n:
                out.add(ecanon((a, b)))
    return sorted(out, key=lambda z: (enorm(z), z))


def hex_classes_closed(t):
    s = 0
    j = 0
    while 3 * j + 1 <= t:
        s += t // (3 * j + 1) - t // (3 * j + 2)
        j += 1
    return s


def hex_classes_direct(t):
    return len(classes_up_to(t))


def eisenstein_primes_dividing(z):
    n = enorm(z)
    rational = []
    m = n
    q = 2
    while q * q <= m:
        if m % q == 0:
            while m % q == 0:
                m //= q
            rational.append(q)
        q += 1
    if m > 1:
        rational.append(m)
    primes = []
    for q in rational:
        if q == 3:
            primes.append(ecanon((1, -1)))
        elif q % 3 == 2:
            primes.append((q, 0))
        else:
            r = 2
            root = None
            while root is None:
                c = pow(r, (q - 1) // 3, q)
                if (c * c + c + 1) % q == 0:
                    root = c
                r += 1
            p = egcd((q, 0), (root, -1))
            primes.append(ecanon(p))
            primes.append(ecanon(econj(p)))
    return [p for p in primes if edivides(p, z)]


def eisenstein_totient(z):
    v = enorm(z)
    if v == 1:
        return 1
    for p in eisenstein_primes_dividing(z):
        v = v // enorm(p) * (enorm(p) - 1)
    return v


def totient_sum(n_max):
    return sum(eisenstein_totient(z) for z in classes_up_to(n_max))


def layer_nodes(z):
    n = enorm(z)
    c = econj(z)
    out = set()
    for m in range(n):
        for k in range(n):
            p, q = emul((m, k), c)
            out.add((Fraction(p % n, n), Fraction(q % n, n)))
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
    g = egcd(u, (d, 0))
    return ecanon(equo((d, 0), g))


def closed_brightness(n_max, d):
    return hex_classes_closed(n_max // enorm(d))


def cyc_reduce(terms):
    acc = {}
    for e, c in terms:
        k = (e % CYC,)
        acc[k] = acc.get(k, 0) + c
    acc = {k: v for k, v in acc.items() if v}
    if not acc:
        return Poly(0, X, domain="QQ")
    return Poly(acc, X, domain="QQ").rem(PHI360)


def field_cos_sin(d):
    c = cyc_reduce([(d, 1), (-d, 1)])
    s = cyc_reduce([(d + 120, 1), (d + 240, -1), (-d + 120, -1), (-d + 240, 1)])
    return c.degree() <= 0, s.degree() <= 0


def field_rotation_degrees():
    return [d for d in range(CYC) if all(field_cos_sin(d))]


def spot_check_degrees():
    out = []
    for d in SPOT_DEGREES:
        c = minimal_polynomial(cos(spi * Rational(d, 180)), X, polys=True).degree()
        s = minimal_polynomial(
            sin(spi * Rational(d, 180)) / sqrt(3), X, polys=True
        ).degree()
        out.append((d, c == 1, s == 1, field_cos_sin(d)))
    return out


def rotation_hits(bound, reach):
    seen = set()
    for a in range(-reach, reach + 1):
        for b in range(-reach, reach + 1):
            if (a, b) == (0, 0):
                continue
            n = enorm((a, b))
            seen.add(
                (
                    Fraction(2 * a * a - 2 * a * b - b * b, 2 * n),
                    Fraction(2 * a * b - b * b, 2 * n),
                )
            )
    direct = set()
    for r in range(1, bound + 1):
        for p in range(-r, r + 1):
            rest = r * r - p * p
            if rest < 0 or rest % 3:
                continue
            s = isqrt(rest // 3)
            if 3 * s * s == rest:
                for t in (s, -s):
                    direct.add((Fraction(p, r), Fraction(t, r)))
    return len(direct), direct <= seen


def chi3(n):
    return (0, 1, -1)[n % 3]


def dirichlet_mul(u, v, m):
    out = [0] * (m + 1)
    for i in range(1, m + 1):
        if u[i] == 0:
            continue
        for j in range(1, m // i + 1):
            out[i * j] += u[i] * v[j]
    return out


def csl_zeta_ratio(m):
    a = [0] * (m + 1)
    for d in range(1, m + 1):
        c = chi3(d)
        if c:
            for n in range(d, m + 1, d):
                a[n] += c
    inv_square = [0] * (m + 1)
    for k in range(1, isqrt(m) + 1):
        f = k
        mu = 1
        q = 2
        while q * q <= f:
            if f % q == 0:
                f //= q
                if f % q == 0:
                    mu = 0
                    break
                mu = -mu
            q += 1
        if mu and f > 1:
            mu = -mu
        inv_square[k * k] = mu
    inv_three = [0] * (m + 1)
    e = 1
    j = 0
    while 3**j <= m:
        inv_three[3**j] = e
        e = -e
        j += 1
    return dirichlet_mul(dirichlet_mul(a, inv_square, m), inv_three, m)


def csl_euler_product(m):
    out = [0] * (m + 1)
    out[1] = 1
    for q in range(2, m + 1):
        if q % 3 != 1 or any(q % r == 0 for r in range(2, isqrt(q) + 1)):
            continue
        local = [0] * (m + 1)
        local[1] = 1
        k = q
        while k <= m:
            local[k] = 2
            k *= q
        out = dirichlet_mul(out, local, m)
    return out


def draw(n_max):
    from PIL import Image, ImageDraw

    stack = literal_stack(n_max)
    top = max(stack.values())
    ss = 3
    w, h = FIG_W * ss, FIG_H * ss
    img = Image.new("RGB", (w, h), SURF)
    pen = ImageDraw.Draw(img)
    unit = 560.0 * ss
    ox = (w - 1.5 * unit) / 2 + unit / 2
    oy = (h + unit * ROOT3 / 2) / 2
    corners = [(0, 0), (1, 0), (1, 1), (0, 1)]
    place = lambda x, y: (ox + unit * (x - y / 2), oy - unit * y * ROOT3 / 2)
    pen.line([place(*c) for c in corners] + [place(0, 0)], fill=GRAY, width=2 * ss)
    for node, b in sorted(stack.items(), key=lambda kv: kv[1]):
        px, py = place(float(node[0]), float(node[1]))
        r = (1.4 + 7.0 * (b / top) ** 0.5) * ss
        pen.ellipse([px - r, py - r, px + r, py + r], fill=BLUE if b < top else INK)
    small = img.resize((FIG_W, FIG_H), Image.LANCZOS)
    small.convert("P", palette=Image.ADAPTIVE, colors=48).save(FIG_OUT, optimize=True)
    return len(stack), top


def main():
    print("ring Z[omega], omega^2 = -1 - omega, norm a^2 - ab + b^2, units 6")
    print(
        "coordinates in the basis 1, omega; fundamental domain the unit square of that chart"
    )

    rots = field_rotation_degrees()
    print("field rotation degrees", *rots)
    print("field rotation count", len(rots))
    spots = spot_check_degrees()
    print("spot degrees", len(spots))
    print(
        "spot minpoly agrees with cyclotomic", all((c, s) == r for _, c, s, r in spots)
    )
    hits, covered = rotation_hits(ROT_BOUND, ROT_REACH)
    print("rational rotations denominator <=", ROT_BOUND, hits)
    print("all are w/conj(w), Eisenstein w in box", ROT_REACH, covered)

    print(
        "hex circle count closed equals direct",
        all(
            hex_classes_closed(t) == hex_classes_direct(t) for t in range(CIRCLE_T + 1)
        ),
    )
    print("hex circle counts t=1..12", *[hex_classes_closed(t) for t in range(1, 13)])

    stack = literal_stack(STACK_N)
    layers = classes_up_to(STACK_N)
    print("stack norm bound", STACK_N)
    print("stack layers", len(layers))
    print("stack nodes", len(stack))
    print("eisenstein totient sum", totient_sum(STACK_N))
    print("nodes equal totient sum", len(stack) == totient_sum(STACK_N))
    bad = 0
    for node, b in stack.items():
        if closed_brightness(STACK_N, reduced_denominator(node)) != b:
            bad += 1
    print("brightness comparisons", len(stack))
    print("brightness mismatches", bad)
    print(
        "origin brightness",
        stack[(Fraction(0), Fraction(0))],
        "max",
        max(stack.values()),
    )

    lval = (szeta(2, Rational(1, 3)) - szeta(2, Rational(2, 3))) / 9
    z2 = spi**2 / 6
    rho = spi / (3 * sqrt(3))
    const = rho / (2 * z2 * lval)
    print("zeta(2)", "%.15f" % float(N_(z2, 30)))
    print("L(2, chi_-3)", "%.15f" % float(N_(lval, 30)))
    print("zeta_K(2) = zeta(2) L(2, chi_-3)", "%.15f" % float(N_(z2 * lval, 30)))
    print("residue of zeta_K at 1 = pi/(3 sqrt 3)", "%.15f" % float(N_(rho, 30)))
    c = float(N_(const, 30))
    print("constant pi/(6 sqrt3 zeta(2) L(2, chi_-3))", "%.15f" % c)
    for n in COUNT_NS:
        s = totient_sum(n)
        r = s / (n * n)
        print(
            "node count N",
            n,
            "=",
            s,
            "ratio N^2",
            "%.6f" % r,
            "deviation",
            "%+.6f" % (r - c),
            "scaled by N/log N",
            "%+.3f" % ((r - c) * n / log(n)),
        )

    ratio = csl_zeta_ratio(CSL_M)
    euler = csl_euler_product(CSL_M)
    print("csl series bound", CSL_M)
    print("csl zeta ratio equals euler product", ratio == euler)
    print(
        "csl nonzero coefficients",
        *[(n, ratio[n]) for n in range(1, CSL_M + 1) if ratio[n]],
    )

    nodes, top = draw(FIG_N)
    print("figure norm bound", FIG_N, "nodes", nodes, "brightest", top)
    print("figure bytes", FIG_OUT.stat().st_size)


main()
