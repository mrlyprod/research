from decimal import Decimal, getcontext
from fractions import Fraction
from math import comb, gcd, isqrt

getcontext().prec = 50

BOUNDS = [50, 200, 800, 1200, 2000, 3200, 12800, 51200, 102400]
SET_BOUNDS = [50, 200, 800]
PUBLISHED = {
    "Q(i)": {50: 672, 200: 10608, 800: 168088},
    "Q(sqrt -3)": {
        50: 630,
        200: 9606,
        800: 151020,
        1200: 337026,
        2000: 945486,
        3200: 2419950,
    },
}
FIELDS = [
    {"name": "Q(i)", "ram": 2, "mod": 4, "w": 4, "disc": 4, "D": -4},
    {"name": "Q(sqrt -3)", "ram": 3, "mod": 3, "w": 6, "disc": 3, "D": -3},
]
DISCS = [-4, -3, -20, -23]
EULER_M = 60
EULER_J = 10


def bernoulli(m):
    b = [Fraction(0)] * (m + 1)
    b[0] = Fraction(1)
    for n in range(1, m + 1):
        s = sum(comb(n + 1, k) * b[k] for k in range(n))
        b[n] = -s / (n + 1)
    return b


def dec(f):
    return Decimal(f.numerator) / Decimal(f.denominator)


def hurwitz2(a, bern):
    base = dec(a)
    s = sum(1 / ((base + k) * (base + k)) for k in range(EULER_M))
    x = base + EULER_M
    s += 1 / x + 1 / (2 * x * x)
    for j in range(1, EULER_J + 1):
        s += dec(bern[2 * j]) / x ** (2 * j + 1)
    return s


def arctan_inv(n):
    total = Decimal(0)
    term = Decimal(1) / Decimal(n)
    sq = Decimal(n) * Decimal(n)
    k = 0
    while True:
        t = term / (2 * k + 1)
        if t < Decimal(10) ** -46:
            return total
        total += t if k % 2 == 0 else -t
        term /= sq
        k += 1


def dpi():
    return 16 * arctan_inv(5) - 4 * arctan_inv(239)


def jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def kron(d, n):
    if n == 0:
        return 0
    r = 1
    while n % 2 == 0:
        if d % 2 == 0:
            return 0
        n //= 2
        r *= 1 if d % 8 in (1, 7) else -1
    return r * jacobi(d % n, n) if n > 1 else r


def units(d):
    return 4 if d == -4 else (6 if d == -3 else 2)


def class_number(d):
    q = -d
    return units(d) * -sum(kron(d, a) * a for a in range(1, q)) // (2 * q)


def lvalue(d, bern):
    q = -d
    s = sum(kron(d, a) * hurwitz2(Fraction(a, q), bern) for a in range(1, q))
    return s / (q * q)


def spf_sieve(n):
    spf = list(range(n + 1))
    for p in range(2, isqrt(n) + 1):
        if spf[p] == p:
            for m in range(p * p, n + 1, p):
                if spf[m] == m:
                    spf[m] = p
    return spf


def factor(n, spf):
    out = []
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        out.append((p, e))
    return out


def kind(fld, p):
    if p == fld["ram"]:
        return "r"
    return "s" if p % fld["mod"] == 1 else "i"


def mulc(fld, u, v):
    x, y = u
    c, d = v
    if fld["ram"] == 2:
        return (x * c - y * d, x * d + y * c)
    return (x * c - y * d, x * d + y * c - y * d)


def conjc(fld, z):
    a, b = z
    return (a, -b) if fld["ram"] == 2 else (a - b, -b)


def nrm(fld, z):
    a, b = z
    return a * a + b * b if fld["ram"] == 2 else a * a - a * b + b * b


def classes(fld, m):
    if fld["ram"] == 2:
        return [
            (a, b) for a in range(1, isqrt(m) + 1) for b in range(isqrt(m - a * a) + 1)
        ]
    r = isqrt(4 * m // 3) + 1
    out = []
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            n = a * a - a * b + b * b
            if 0 < n <= m:
                z = (a, b)
                best = z
                for _ in range(3):
                    z = (-z[1], z[0] - z[1])
                    best = min(best, z, (-z[0], -z[1]))
                if best == (a, b):
                    out.append((a, b))
    return out


def totient(fld, z, n, spf):
    v = n
    a, b = z
    for p, e in factor(n, spf):
        k = kind(fld, p)
        if k == "r":
            v = v // p * (p - 1)
        elif k == "i":
            v = v // (p * p) * (p * p - 1)
        else:
            v = v // p * (p - 1)
            if a % p == 0 and b % p == 0:
                v = v // p * (p - 1)
    return v


def mobius(fld, z, n, spf):
    v = 1
    a, b = z
    for p, e in factor(n, spf):
        k = kind(fld, p)
        if k == "r":
            if e != 1:
                return 0
            v = -v
        elif k == "i":
            if e != 2:
                return 0
            v = -v
        else:
            if e == 1:
                v = -v
            elif e == 2 and a % p == 0 and b % p == 0:
                pass
            else:
                return 0
    return v


def tables(fld, m, spf):
    cls = classes(fld, m)
    phi = [0] * (m + 1)
    nsum = [0] * (m + 1)
    mob = [0] * (m + 1)
    for z in cls:
        n = nrm(fld, z)
        phi[n] += totient(fld, z, n, spf)
        nsum[n] += n
        mob[n] += mobius(fld, z, n, spf)
    for n in range(1, m + 1):
        phi[n] += phi[n - 1]
        nsum[n] += nsum[n - 1]
    return cls, phi, nsum, mob


def norm_arrays(d, m):
    chi = [kron(d, n) for n in range(m + 1)]
    ideals = [0] * (m + 1)
    for k in range(1, m + 1):
        if chi[k]:
            for n in range(k, m + 1, k):
                ideals[n] += chi[k]
    prime = [True] * (m + 1)
    prime[0] = prime[1] = False
    for k in range(2, isqrt(m) + 1):
        if prime[k]:
            for n in range(k * k, m + 1, k):
                prime[n] = False
    mu = [1] * (m + 1)
    sq = [True] * (m + 1)
    for k in range(2, m + 1):
        if prime[k]:
            for n in range(k, m + 1, k):
                mu[n] = -mu[n]
            for n in range(k * k, m + 1, k * k):
                sq[n] = False
    mk = [0] * (m + 1)
    for k in range(1, m + 1):
        if sq[k] and mu[k]:
            for j in range(1, m // k + 1):
                if sq[j] and chi[j]:
                    mk[k * j] += mu[k] * mu[j] * chi[j]
    tt = [0] * (m + 1)
    for n in range(1, m + 1):
        tt[n] = tt[n - 1] + n * ideals[n]
    return mk, tt


def norm_totient_sum(n, mk, tt):
    return sum(mk[k] * tt[n // k] for k in range(1, n + 1) if mk[k])


def totient_sum(phi, m):
    return phi[m]


def convolution_sum(m, nsum, mob):
    return sum(mob[n] * nsum[m // n] for n in range(1, m + 1) if mob[n])


def hnf(r1, r2):
    a1, b1 = r1
    a2, b2 = r2
    x, y, d = 1, 0, a1
    u, v, e = 0, 1, a2
    while e:
        q = d // e
        d, e = e, d - q * e
        x, u = u, x - q * u
        y, v = v, y - q * v
    if d < 0:
        d, x, y = -d, -x, -y
    return d, x * b1 + y * b2, abs((a2 // d) * b1 - (a1 // d) * b2)


def farey_set(fld, m):
    seen = set()
    for z in classes(fld, m):
        n = nrm(fld, z)
        r1 = mulc(fld, (1, 0), z)
        r2 = mulc(fld, (0, 1), z)
        h11, _h12, h22 = hnf(r1, r2)
        cj = conjc(fld, z)
        for x in range(h11):
            for y in range(h22):
                u, v = mulc(fld, (x, y), cj)
                u %= n
                v %= n
                g = gcd(gcd(u, v), n)
                seen.add((u // g, v // g, n // g))
    return len(seen)


def d12(x):
    return str(+x.quantize(Decimal("1.000000000000")))


def main():
    bern = bernoulli(2 * EULER_J)
    pi = dpi()
    z2 = pi * pi / 6
    lv = {d: lvalue(d, bern) for d in DISCS}
    print("pi", d12(pi))
    print("zeta(2)", d12(z2))
    print("L(2, chi_-4) Catalan", d12(lv[-4]))
    print("L(2, chi_-3)", d12(lv[-3]))
    const = {}
    for fld in FIELDS:
        dk = fld["disc"]
        zk = z2 * lv[fld["D"]]
        rho = 2 * pi / (fld["w"] * Decimal(dk).sqrt())
        c = rho / (2 * zk)
        const[fld["name"]] = c
        print("field", fld["name"], "w", fld["w"], "|D|", dk)
        print("  zeta_K(2)", d12(zk))
        print("  rho_K", d12(rho))
        print("  c = rho_K/(2 zeta_K(2))", d12(c))
        print(
            "  Sayous c_K = pi/(sqrt|D| zeta_K(2))",
            d12(pi / (Decimal(dk).sqrt() * zk)),
            "= w c",
            d12(fld["w"] * c),
        )
    spf = spf_sieve(max(BOUNDS))
    byclass = {}
    for fld in FIELDS:
        name = fld["name"]
        c = const[name]
        cls, phi, nsum, mob = tables(fld, max(BOUNDS), spf)
        byclass[fld["D"]] = phi
        print("field", name, "classes to", max(BOUNDS), len(cls))
        for m in BOUNDS:
            s = totient_sum(phi, m)
            dn = Decimal(m)
            main_term = c * dn * dn
            ratio = Decimal(s) / main_term
            dev = Decimal(s) - main_term
            mark = PUBLISHED[name].get(m)
            tag = (
                "published " + str(mark) + " match " + str(mark == s) if mark else "new"
            )
            print(
                "  N",
                m,
                "count",
                s,
                "count/(c N^2)",
                d12(ratio),
                "dev/N^1.5",
                d12(dev / (dn * dn.sqrt())),
                "dev/(N ln N)",
                d12(dev / (dn * dn.ln())),
                tag,
            )
        bad = [
            m for m in BOUNDS if convolution_sum(m, nsum, mob) != totient_sum(phi, m)
        ]
        print("  convolution identity mismatches", len(bad))
        for m in SET_BOUNDS:
            k = farey_set(fld, m)
            print(
                "  Farey set size at N",
                m,
                k,
                "equals totient sum",
                k == totient_sum(phi, m),
            )
    mx = max(BOUNDS)
    for d in DISCS:
        h = class_number(d)
        w = units(d)
        zk = z2 * lv[d]
        rho = 2 * pi * h / (w * Decimal(-d).sqrt())
        c = rho / (2 * zk)
        mk, tt = norm_arrays(d, mx)
        print(
            "discriminant",
            d,
            "h",
            h,
            "w",
            w,
            "L(2, chi_D)",
            d12(lv[d]),
            "rho_K",
            d12(rho),
            "c",
            d12(c),
        )
        for m in BOUNDS:
            s = norm_totient_sum(m, mk, tt)
            dn = Decimal(m)
            main_term = c * dn * dn
            tail = (
                " class route " + str(s == totient_sum(byclass[d], m))
                if d in byclass
                else ""
            )
            print(
                "  N",
                m,
                "count",
                s,
                "count/(c N^2)",
                d12(s / main_term),
                "dev/N^1.5",
                d12((s - main_term) / (dn * dn.sqrt())) + tail,
            )


main()
