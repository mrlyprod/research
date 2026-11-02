import math
import random
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np

PREC = 60
EULER_TERMS = 100
BERNOULLI = [
    Fraction(1, 6),
    Fraction(-1, 30),
    Fraction(1, 42),
    Fraction(-1, 30),
    Fraction(5, 66),
    Fraction(-691, 2730),
    Fraction(7, 6),
    Fraction(-3617, 510),
    Fraction(43867, 798),
    Fraction(-174611, 330),
]

PAIRS = 200000
CUTOFF = 4000000
SIDE = 3000
TRIALS = 400000
RADIUS = 500
SEED = 20260725


def atan_inv(x):
    base = Decimal(x)
    term = Decimal(1) / base
    total = term
    square = base * base
    eps = Decimal(10) ** -58
    k = 1
    sign = -1
    while True:
        term = term / square
        add = Decimal(sign) * term / Decimal(2 * k + 1)
        if abs(add) < eps:
            break
        total += add
        sign = -sign
        k += 1
    return total


def machin_pi():
    return 16 * atan_inv(5) - 4 * atan_inv(239)


def trigamma(x):
    total = Fraction(0)
    for n in range(EULER_TERMS):
        shifted = n + x
        total += 1 / (shifted * shifted)
    m = EULER_TERMS + x
    total += 1 / m
    total += Fraction(1, 2) / (m * m)
    for index, bern in enumerate(BERNOULLI):
        total += bern / m ** (2 * (index + 1) + 1)
    return total


def dec(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


def contfrac(value, terms):
    out = []
    x = value
    for _ in range(terms):
        a = int(x)
        out.append(a)
        frac = x - a
        if frac == 0:
            break
        x = Decimal(1) / frac
    return out


def denominators(quotients):
    out = []
    previous, current = 0, 1
    for index, a in enumerate(quotients):
        if index == 0:
            out.append(current)
            continue
        previous, current = current, a * current + previous
        out.append(current)
    return out


def cut(value, width):
    return str(value)[:width]


def precision_block():
    getcontext().prec = PREC
    pi = machin_pi()
    pi2 = pi * pi
    psi13 = trigamma(Fraction(1, 3))
    psi23 = trigamma(Fraction(2, 3))
    psi14 = trigamma(Fraction(1, 4))
    psi34 = trigamma(Fraction(3, 4))
    l3 = dec(psi13 - psi23) / 9
    catalan = dec(psi14 - psi34) / 16
    z2 = pi2 / 6
    zk3 = z2 * l3
    zki = z2 * catalan
    sqrt3 = Decimal(3).sqrt()
    reflection = dec(psi13 + psi23) - 4 * pi2 / 3
    cl2 = l3 * 3 * sqrt3 / 4
    r1 = l3 * sqrt3 / pi2
    r2 = l3 / pi2
    print("pi           =", cut(pi, 45))
    print("zeta(2)      =", cut(z2, 45))
    print("6/pi^2       =", cut(6 / pi2, 45))
    print("L(2,chi_-3)  =", cut(l3, 45))
    print("Catalan G    =", cut(catalan, 45))
    print("zeta_K3(2)   =", cut(zk3, 45))
    print("1/zeta_K3(2) =", cut(1 / zk3, 45))
    print("zeta_Qi(2)   =", cut(zki, 45))
    print("1/zeta_Qi(2) =", cut(1 / zki, 45))
    print("psi1(1/3)+psi1(2/3)-4pi^2/3 =", reflection)
    print("Cl2(pi/3) implied =", cut(cl2, 30))
    for label, value, terms in (
        ("L3*sqrt3/pi^2", r1, 22),
        ("L3/pi^2      ", r2, 22),
        ("G/pi^2       ", catalan / pi2, 22),
        ("zeta(2)/pi^2 ", z2 / pi2, 8),
    ):
        quotients = contfrac(value, terms)
        print(f"{label} = {cut(value, 30)} CF: {quotients}")
        print(f"{label}   last convergent denominator: {denominators(quotients)[-1]}")
    return float(zk3), float(zki), float(1 / zk3)


def dirichlet(pairs):
    m = np.arange(pairs, dtype=np.float64)
    odd = 3 * m + 1
    even = 3 * m + 2
    return math.fsum(1.0 / (odd * odd) - 1.0 / (even * even))


def lattice_radius(cutoff):
    return int(2.2 * math.sqrt(cutoff))


def lattice_sums(cutoff):
    r = lattice_radius(cutoff)
    b = np.arange(-r, r + 1, dtype=np.int64)
    hexagonal = 0.0
    square = 0.0
    for a in range(-r, r + 1):
        hexnorm = a * a - a * b + b * b
        sqnorm = a * a + b * b
        take = hexnorm[(hexnorm > 0) & (hexnorm <= cutoff)].astype(np.float64)
        hexagonal += float(np.sum(1.0 / (take * take)))
        take = sqnorm[(sqnorm > 0) & (sqnorm <= cutoff)].astype(np.float64)
        square += float(np.sum(1.0 / (take * take)))
    return hexagonal / 6.0, square / 4.0


def visible_pairs(side):
    axis = np.arange(1, side + 1, dtype=np.int64)
    count = 0
    for a in range(1, side + 1):
        count += int(np.count_nonzero(np.gcd(a, axis) == 1))
    return count


def emul(z, w):
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c - b * d)


def enorm(z):
    a, b = z
    return a * a - a * b + b * b


def edivmod(z, w):
    c, d = w
    n = enorm(w)
    p = emul(z, (c - d, -d))
    q = (round(p[0] / n), round(p[1] / n))
    s = emul(q, w)
    return q, (z[0] - s[0], z[1] - s[1])


def egcd(z, w):
    while enorm(w) != 0:
        z, w = w, edivmod(z, w)[1]
    return z


def seeded_sieve(trials, radius):
    rng = random.Random(SEED)
    coprime = 0
    for _ in range(trials):
        z = (rng.randint(-radius, radius), rng.randint(-radius, radius))
        w = (rng.randint(-radius, radius), rng.randint(-radius, radius))
        if enorm(z) == 0 or enorm(w) == 0:
            continue
        if enorm(egcd(z, w)) == 1:
            coprime += 1
    return coprime


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            r -= r // p
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        r -= r // m
    return r


def fresh(base, level):
    return sum(1 for k in range(1, base**level + 1) if k % base)


def sieve_block(zk3, zki, inv_zk3):
    print(f"L(2,chi_-3) direct Dirichlet sum ({PAIRS} pairs) = {dirichlet(PAIRS):.12f}")
    hexagonal, square = lattice_sums(CUTOFF)
    hextail = math.pi / (3 * math.sqrt(3)) / CUTOFF
    sqtail = math.pi / 4 / CUTOFF
    print(
        f"zeta_K3(2) lattice sum (norm<={CUTOFF}) = {hexagonal:.10f}"
        f"  (+tail~{hextail:.2e} -> {hexagonal + hextail:.10f})"
    )
    print(f"zeta_K3(2) expected zeta(2)*L       = {zk3:.10f}")
    print(
        f"zeta_Qi(2) lattice sum (norm<={CUTOFF}) = {square:.10f}"
        f"  (+tail~{sqtail:.2e} -> {square + sqtail:.10f})"
    )
    print(f"zeta_Qi(2) expected zeta(2)*G       = {zki:.10f}")
    seen = visible_pairs(SIDE)
    print(
        f"square-lattice visible density N={SIDE}: {seen / SIDE**2:.6f}"
        f"  vs 6/pi^2 = {6 / math.pi**2:.6f}"
    )
    coprime = seeded_sieve(TRIALS, RADIUS)
    print(
        f"Eisenstein sieve: {coprime}/{TRIALS} coprime = {coprime / TRIALS:.5f}"
        f"  vs 1/zeta_K3(2) = {inv_zk3:.5f}"
    )
    assert enorm(emul((3, 2), (1, 5))) == enorm((3, 2)) * enorm((1, 5))
    print(
        f"gcd(6,9) in Z[omega] has norm {enorm(egcd((6, 0), (9, 0)))} (expect 9 = N(3))"
    )
    for base in (2, 3, 5):
        for level in range(1, 5):
            new = fresh(base, level)
            assert new == base**level - base ** (level - 1) == phi(base**level)
    print(
        "fresh nodes per level = b^L - b^(L-1) = phi(b^L)"
        " for b=2,3,5, L=1..4: OK; fraction (b-1)/b"
    )
    print(
        f"b=6,L=2: new={fresh(6, 2)}, b^L-b^(L-1)={6**2 - 6}, phi(b^L)={phi(36)}"
        " (phi identity needs prime b)"
    )


def main():
    print(
        f"domain: decimal prec {PREC}, Euler-Maclaurin {EULER_TERMS} terms through B_20"
    )
    print(
        f"domain: Dirichlet {PAIRS} pairs, lattice norm <= {CUTOFF}"
        f" (radius {lattice_radius(CUTOFF)}), census [1,{SIDE}]^2,"
        f" sieve {TRIALS} pairs in [-{RADIUS},{RADIUS}]^2"
    )
    zk3, zki, inv_zk3 = precision_block()
    sieve_block(zk3, zki, inv_zk3)


main()
