import math
from fractions import Fraction

import numpy as np

N_CONV = 60
N_SELECT = [30, 61, 200, 501]
HARMONIC_NS = [1000, 4000, 16000]
S_VALUES = [3, 4, 5]
PRIME_LS = [5, 10, 100, 1000]
ODD_LS = [10, 100, 1000, 4000]
SQUAREFREE_LS = [10, 50, 200, 1000]
DAVENPORT_NS = [10**3, 10**4, 10**5]
DAVENPORT_X = [
    Fraction(1, 7),
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(3, 11),
]

# SIEVES


def mobius_table(n):
    mu = np.ones(n + 1, dtype=np.int64)
    composite = np.zeros(n + 1, dtype=bool)
    for p in range(2, n + 1):
        if composite[p]:
            continue
        composite[p::p] = True
        mu[p::p] *= -1
        square = p * p
        if square <= n:
            mu[square::square] = 0
    mu[0] = 0
    return mu


def primes_upto(n):
    flag = np.ones(n + 1, dtype=bool)
    flag[:2] = False
    for p in range(2, int(n**0.5) + 1):
        if flag[p]:
            flag[p * p :: p] = False
    return np.nonzero(flag)[0]


def first_primes(count):
    limit = 64
    while True:
        found = primes_upto(limit)
        if found.size >= count:
            return [int(p) for p in found[:count]]
        limit *= 2


def prime_powers_upto(n):
    out = {}
    for p in primes_upto(n):
        p = int(p)
        power = p
        exponent = 1
        while power <= n:
            out[power] = (p, exponent)
            power *= p
            exponent += 1
    return out


def totients(n):
    phi = np.arange(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if phi[p] == p:
            phi[p::p] -= phi[p::p] // p
    return phi


# DIRICHLET CONVOLUTION


def convolve(u, v, n):
    out = [Fraction(0)] * (n + 1)
    for k in range(1, n + 1):
        if u[k] == 0:
            continue
        for m in range(1, n // k + 1):
            out[k * m] += u[k] * v[m]
    return out


def literal_stack(w, n):
    nodes = {}
    for scale in range(1, n + 1):
        if w[scale] == 0:
            continue
        for j in range(1, scale + 1):
            g = math.gcd(j, scale)
            key = (j // g, scale // g)
            nodes[key] = nodes.get(key, Fraction(0)) + w[scale]
    return nodes


def literal_double_stack(u, v, n):
    nodes = {}
    for k in range(1, n + 1):
        if u[k] == 0:
            continue
        for m in range(1, n // k + 1):
            weight = u[k] * v[m]
            if weight == 0:
                continue
            scale = k * m
            for j in range(1, scale + 1):
                g = math.gcd(j, scale)
                key = (j // g, scale // g)
                nodes[key] = nodes.get(key, Fraction(0)) + weight
    return nodes


def nested_stack(u, v, n):
    nodes = {}
    for k in range(1, n + 1):
        if u[k] == 0:
            continue
        for m in range(1, n // k + 1):
            weight = u[k] * v[m]
            if weight == 0:
                continue
            for i in range(1, k + 1):
                for j in range(1, m + 1):
                    num = (i - 1) * m + j
                    den = k * m
                    g = math.gcd(num, den)
                    key = (num // g, den // g)
                    nodes[key] = nodes.get(key, Fraction(0)) + weight
    return nodes


def rectangular_weights(u, v, outer, inner):
    out = {}
    for k in range(1, outer + 1):
        for m in range(1, inner + 1):
            out[k * m] = out.get(k * m, Fraction(0)) + u[k] * v[m]
    return out


def node_formula(w, n, b):
    return sum((w[k * b] for k in range(1, n // b + 1)), Fraction(0))


def convolution_check():
    n = N_CONV
    mu = mobius_table(n)
    one = [Fraction(0)] + [Fraction(1)] * n
    mob = [Fraction(0)] + [Fraction(int(mu[k])) for k in range(1, n + 1)]
    inv = [Fraction(0)] + [Fraction(1, k) for k in range(1, n + 1)]
    pairs = [
        ("1 * 1", one, one),
        ("1 * mu", one, mob),
        ("mu * mu", mob, mob),
        ("1 * n^-1", one, inv),
    ]
    print("CONVOLUTION  hyperbolic cut k*m <= N, N =", n)
    for name, u, v in pairs:
        w = convolve(u, v, n)
        left = literal_double_stack(u, v, n)
        right = literal_stack(w, n)
        nested = nested_stack(u, v, n)
        keys = set(left) | set(right) | set(nested)
        bad = sum(
            1
            for key in keys
            if left.get(key, Fraction(0)) != right.get(key, Fraction(0))
        )
        geo = sum(
            1
            for key in keys
            if nested.get(key, Fraction(0)) != right.get(key, Fraction(0))
        )
        form = sum(
            1
            for (a, b) in keys
            if right.get((a, b), Fraction(0)) != node_formula(w, n, b)
        )
        rect = rectangular_weights(u, v, 24, 40)
        cut = sum(1 for m in range(1, 25) if rect.get(m, Fraction(0)) != w[m])
        above = sum(1 for m in range(25, 41) if rect.get(m, Fraction(0)) != w[m])
        print(
            " ",
            name,
            "nodes",
            len(keys),
            "scaled-copy mismatches",
            geo,
            "hyperbolic mismatches",
            bad,
            "formula mismatches",
            form,
        )
        print(
            "   ",
            name,
            "rectangular cut 24 x 40: mismatches at m <= 24:",
            cut,
            " at 25 <= m <= 40:",
            above,
        )
    e = convolve(one, mob, n)
    print(
        "  1 * mu identity  e(1) =",
        e[1],
        " nonzero beyond 1:",
        sum(1 for k in range(2, n + 1) if e[k] != 0),
    )
    d = convolve(one, one, n)
    print("  1 * 1 is d(n), first ten", [int(d[k]) for k in range(1, 11)])
    print("  scale weight W(M) = (u*v)(M) for every M <= N under the hyperbolic cut")


# SELECTIONS


def selection_literal(weight, n, b):
    return sum(1 for k in range(1, n // b + 1) if weight(k * b))


def even_form(n, b):
    return n // (b if b % 2 == 0 else 2 * b)


def odd_form(n, b):
    return 0 if b % 2 == 0 else (n // b + 1) // 2


def prime_form(n, b, prime_set, pi_n):
    if b == 1:
        return pi_n
    return 1 if b in prime_set else 0


def coprime_count(y, b_divisors_mu):
    return sum(m * (y // e) for e, m in b_divisors_mu)


def squarefree_form(n, b, mu):
    if mu[b] == 0:
        return 0
    divisors = [(e, int(mu[e])) for e in range(1, b + 1) if b % e == 0 and mu[e] != 0]
    x = n // b
    total = 0
    d = 1
    while d * d <= x:
        if mu[d] != 0 and math.gcd(d, b) == 1:
            total += int(mu[d]) * coprime_count(x // (d * d), divisors)
        d += 1
    return total


def prime_power_form(n, b, powers):
    if b == 1:
        return sum(int(math.log(n, p) + 1e-9) for p in primes_upto(n))
    if b not in powers:
        return 0
    p, i = powers[b]
    top = int(math.log(n, p) + 1e-9)
    return max(0, top - i + 1)


def selection_closed_forms():
    print("SELECTIONS  closed form against literal stacking, all b <= N")
    for n in N_SELECT:
        mu = mobius_table(n)
        prime_set = {int(p) for p in primes_upto(n)}
        powers = prime_powers_upto(n)
        pi_n = len(prime_set)
        rows = []
        checks = [
            ("evens", lambda m: m % 2 == 0, lambda b: even_form(n, b)),
            ("odds", lambda m: m % 2 == 1, lambda b: odd_form(n, b)),
            (
                "primes",
                lambda m: m in prime_set,
                lambda b: prime_form(n, b, prime_set, pi_n),
            ),
            ("squarefree", lambda m: mu[m] != 0, lambda b: squarefree_form(n, b, mu)),
            (
                "prime powers",
                lambda m: m in powers,
                lambda b: prime_power_form(n, b, powers),
            ),
        ]
        for name, weight, form in checks:
            bad = sum(
                1 for b in range(1, n + 1) if selection_literal(weight, n, b) != form(b)
            )
            rows.append("%s %d" % (name, bad))
        print("  N =", n, " mismatches:", ", ".join(rows))
    n = N_SELECT[-1]
    prime_set = {int(p) for p in primes_upto(n)}
    lit = sum(
        1
        for b in range(1, n + 1)
        if selection_literal(lambda m: m in prime_set, n, b) > 0
    )
    print(
        "  primes-only line stack at N =",
        n,
        "lights",
        lit,
        "denominators: b = 1 and the",
        len(prime_set),
        "primes",
    )


# HARMONIC STACK


def zeta_em(s, terms=200000):
    k = np.arange(1, terms + 1, dtype=np.float64)
    total = float(np.sum(k ** (-float(s))))
    m = float(terms)
    return total + m ** (1 - s) / (s - 1) - 0.5 * m ** (-s) + s * m ** (-s - 1) / 12.0


def harmonic_stack():
    print("HARMONIC STACK  weights n^-s, node a/b reads b^-s H_s(floor(N/b))")
    for s in S_VALUES:
        target = zeta_em(s - 1)
        row = []
        for n in HARMONIC_NS:
            k = np.arange(1, n + 1, dtype=np.float64)
            partial = np.concatenate((np.zeros(1), np.cumsum(k ** (-float(s)))))
            phi = totients(n)
            b = np.arange(1, n + 1, dtype=np.float64)
            mass = float(
                np.sum(phi[1:] * b ** (-float(s)) * partial[n // np.arange(1, n + 1)])
            )
            row.append("N=%d %.6f" % (n, mass))
        print(
            "  s =", s, " zeta(s-1) = %.6f" % target, " total node mass", ", ".join(row)
        )
    n = HARMONIC_NS[-1]
    s = S_VALUES[0]
    k = np.arange(1, n + 1, dtype=np.float64)
    partial = np.concatenate((np.zeros(1), np.cumsum(k ** (-float(s)))))
    z = zeta_em(s)
    for b in (1, 2, 5, 17):
        node = float(partial[n // b]) * b ** (-float(s))
        literal = float(
            np.sum((b * np.arange(1, n // b + 1, dtype=np.float64)) ** (-float(s)))
        )
        print(
            "  s = %d  b = %2d  node %.9f  literal %.9f  zeta(s)/b^s %.9f"
            % (s, b, node, literal, z / b**s)
        )


# CARPET LAYERS


def sign_at(n, i, length):
    return 1 if ((n * i) // length) % 2 == 0 else -1


def mean_sign(m):
    return Fraction(sum(sign_at(m, i, m) for i in range(m)), m)


def pair_integral(m, n):
    length = m * n // math.gcd(m, n)
    return Fraction(
        sum(sign_at(m, i, length) * sign_at(n, i, length) for i in range(length)),
        length,
    )


def carpet_exact(m, n):
    am, an = mean_sign(m), mean_sign(n)
    a = pair_integral(m, n)
    mu_m, mu_n = (1 - am) / 2, (1 - an) / 2
    nu = (1 - am - an + a) / 4
    return nu * nu - (mu_m * mu_n) ** 2


def carpet_var(n):
    q = Fraction(n - 1, 2 * n)
    return q**2 - q**4


def carpet_cov(m, n):
    d = math.gcd(m, n)
    return Fraction(
        (d * d - 1) * (2 * (m - 1) * (n - 1) + d * d - 1), 16 * m * m * n * n
    )


def pair_sum(layers):
    a = np.array(layers, dtype=np.int64)
    m = a.astype(np.float64)
    var = float(np.sum((m - 1) ** 2 / (4 * m**2) - (m - 1) ** 4 / (16 * m**4)))
    off = 0.0
    hits = 0
    step = 512
    for start in range(0, a.size, step):
        block = a[start : start + step]
        d = np.gcd.outer(block, a).astype(np.float64)
        x = block.astype(np.float64)[:, None]
        y = m[None, :]
        cov = (d * d - 1) * (2 * (x - 1) * (y - 1) + d * d - 1) / (16 * x * x * y * y)
        cov[np.arange(block.size), np.arange(start, start + block.size)] = 0.0
        off += float(np.sum(cov))
        hits += int(np.count_nonzero(d > 1))
    return var, off, (hits - a.size) // 2


def prime_carpet_variance():
    print("CARPET LAYER LAW  closed form against exact rational integration")
    for m, n in ((3, 5), (5, 7), (3, 9), (15, 21), (5, 15)):
        print(
            "  (%d,%d) exact %s closed %s equal %s"
            % (
                m,
                n,
                carpet_exact(m, n),
                carpet_cov(m, n),
                carpet_exact(m, n) == carpet_cov(m, n),
            )
        )
    twelve = first_primes(13)[1:]
    bad = sum(
        1
        for i, p in enumerate(twelve)
        for q in twelve[i + 1 :]
        if carpet_exact(p, q) != 0
    )
    print("  odd prime pairs among", twelve, "with nonzero exact covariance:", bad)
    print("PRIME CARPET STACK  L*Var of the L-layer mean, odd primes")
    for count in PRIME_LS:
        ps = first_primes(count + 1)[1:]
        total = sum((carpet_var(p) for p in ps), Fraction(0))
        var, off, sharing = pair_sum(ps)
        exact = str(total / count) if count <= 5 else "-"
        print(
            "  L = %4d  pairs sharing a factor %d  covariance sum %.1f  L*Var = %.10f  independent %.10f  ratio %.10f"
            % (
                count,
                sharing,
                off,
                float(total / count),
                var / count,
                (var + off) / var,
            )
        )
        if exact != "-":
            print("    exact L*Var = %s" % exact)
    print(
        "  limit 3/16 = %.10f  c = sqrt(3)/4 = %.7f  c factor over independent exactly 1 at every L"
        % (3 / 16, math.sqrt(3) / 4)
    )
    bad = sum(
        1
        for p in first_primes(1001)[1:]
        if 16 * carpet_var(p) * p**4 != 3 * p**4 - 4 * p**3 - 2 * p**2 + 4 * p - 1
    )
    print(
        "  16 p^4 Var_p = 3p^4 - 4p^3 - 2p^2 + 4p - 1 breaches over the first 1000 odd primes:",
        bad,
    )
    print("ODD CARPET STACK  the tree's full odd stack for comparison")
    for count in ODD_LS:
        layers = list(range(3, 2 * count + 3, 2))
        var, off, _ = pair_sum(layers)
        actual = (var + off) / count
        indep = var / count
        print(
            "  L = %4d  L*Var = %.7f  independent %.7f  ratio %.6f  c factor %.6f"
            % (count, actual, indep, actual / indep, math.sqrt(actual / indep))
        )


def squarefree_carpet_variance():
    print("SQUAREFREE ODD STACK  the rival selection")
    mu = mobius_table(20000)
    pool = [n for n in range(3, 20001, 2) if mu[n] != 0]
    print(
        "  Cov(15,21) = %s = %.10f  Pearson r = %.10f"
        % (
            carpet_cov(15, 21),
            float(carpet_cov(15, 21)),
            float(
                carpet_cov(15, 21) / (carpet_var(15) * carpet_var(21)) ** Fraction(1, 2)
            ),
        )
    )
    for count in SQUAREFREE_LS:
        layers = pool[:count]
        var, off, sharing = pair_sum(layers)
        actual = (var + off) / count
        indep = var / count
        print(
            "  L = %4d  pairs sharing a factor %6d  L*Var = %.7f  independent %.7f  ratio %.6f"
            % (count, sharing, actual, indep, actual / indep)
        )


# DAVENPORT


def davenport():
    print("DAVENPORT  sum mu(n)/n ((nx)) against -sin(2 pi x)/pi")
    top = DAVENPORT_NS[-1]
    mu = mobius_table(top)
    n = np.arange(1, top + 1, dtype=np.int64)
    weight = mu[1:].astype(np.float64) / n.astype(np.float64)
    errors = {}
    for x in DAVENPORT_X:
        a, q = x.numerator, x.denominator
        r = (n * a) % q
        frac = r.astype(np.float64) / q
        saw = np.where(r == 0, 0.0, frac - 0.5)
        target = -math.sin(2 * math.pi * float(x)) / math.pi
        row = []
        for cut in DAVENPORT_NS:
            row.append("%d %.6f" % (cut, float(np.sum(weight[:cut] * saw[:cut]))))
        tail = float(np.sum(weight * frac))
        errors[x] = [
            abs(float(np.sum(weight[:cut] * saw[:cut])) - target)
            for cut in DAVENPORT_NS
        ]
        print(
            "  x = %s  target %.6f  partial sums %s  {nx} form %.6f"
            % (x, target, ", ".join(row), tail)
        )
    for i, cut in enumerate(DAVENPORT_NS):
        print(
            "  cut n <= %6d  max error over the five points %.2e"
            % (cut, max(e[i] for e in errors.values()))
        )
    print(
        "  sum mu(n)/n to %d = %.8f  the {nx} form needs this zero, which is PNT-equivalent"
        % (top, float(np.sum(weight)))
    )


def main():
    convolution_check()
    selection_closed_forms()
    harmonic_stack()
    prime_carpet_variance()
    squarefree_carpet_variance()
    davenport()


main()
