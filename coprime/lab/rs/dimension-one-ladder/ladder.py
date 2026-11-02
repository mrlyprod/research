from fractions import Fraction
from math import comb, log

import mpmath
import numpy as np
import sympy as sp

DIGITS = ((0, 0), (1, 0), (0, 1))
LOG3 = log(3)
KAPPA = 3 - log(5) / LOG3
WALL = 2 / (3 + log(5) / LOG3)


def log3(x):
    return log(x) / LOG3


def beta(kappa):
    return kappa / (2 * kappa + 2 - KAPPA)


def delta_counts(K):
    types = [
        (b, c, comb(K, b) * comb(K - b, c))
        for b in range(K + 1)
        for c in range(K + 1 - b)
    ]
    out = {}
    for b, c, m in types:
        for B, C, M in types:
            out[(b - B, c - C)] = out.get((b - B, c - C), 0) + m * M
    return out


def carry_matrix(K):
    radius = (K - 1) // 2
    states = [
        (i, j) for i in range(-radius, radius + 1) for j in range(-radius, radius + 1)
    ]
    index = {s: i for i, s in enumerate(states)}
    M = [[0] * len(states) for _ in states]
    for s in states:
        for (d1, d2), w in delta_counts(K).items():
            z1, z2 = s[0] + d1, s[1] + d2
            if z1 % 3 == 0 and z2 % 3 == 0:
                M[index[s]][index[(z1 // 3, z2 // 3)]] += w
    return M, index[(0, 0)]


def energies(K, levels):
    M, zero = carry_matrix(K)
    v = [0] * len(M)
    v[zero] = 1
    out = [1]
    for _ in range(levels):
        v = [sum(v[i] * M[i][j] for i in range(len(v))) for j in range(len(v))]
        out.append(v[zero])
    return out


def gasket(a):
    pts = [(0, 0)]
    for l in range(a):
        pts = [(x + dx * 3**l, y + dy * 3**l) for x, y in pts for dx, dy in DIGITS]
    return pts


def energy_direct(K, a):
    pts = gasket(a)
    top = K * (3**a - 1) + 1
    sums = np.zeros((top, top), dtype=np.int64)
    sums[0, 0] = 1
    for _ in range(K):
        fresh = np.zeros_like(sums)
        for x, y in pts:
            fresh[x:, y:] += sums[: top - x, : top - y]
        sums = fresh
    return sum(int(v) * int(v) for v in sums.flat if v)


print(
    "kappa = 3 - log_3 5 = %.6f, beta_0^(4) = %.6f, wall 2/(3 + log_3 5) = %.12f"
    % (KAPPA, beta(KAPPA), WALL)
)

print("energies E_2K(G_a), carry matrix exact integers, a = 0..6")
rows = {}
for K in (2, 3, 4, 5):
    rows[K] = energies(K, 6)
    print("2K = %d:" % (2 * K), rows[K])
print("E_4 = 15^a for a <= 6:", all(rows[2][a] == 15**a for a in range(7)))
for K in (2, 3, 4, 5):
    top = 4 if K < 5 else 3
    agree = all(energy_direct(K, a) == rows[K][a] for a in range(1, top + 1))
    print("2K = %d direct convolution agrees at a = 1..%d: %s" % (2 * K, top, agree))

print("characteristic polynomials and Perron roots")
x = sp.symbols("x")
mpmath.mp.dps = 30
perron = {}
for K in (3, 4, 5):
    M, _ = carry_matrix(K)
    poly = sp.Matrix(M).charpoly(x).as_expr()
    print("2K = %d states %d factors: %s" % (2 * K, len(M), sp.factor(poly)))
    perron[K] = max(
        mpmath.mpf(str(r))
        for f, _ in sp.factor_list(poly)[1]
        for r in sp.Poly(f, x).nroots(n=25)
        if r.is_real
    )
lam6 = 57 + 6 * mpmath.sqrt(46)
lam8 = 456 + 3 * mpmath.sqrt(11017)
lam10 = perron[5]
print(
    "lambda_6 = 57 + 6 sqrt 46 = %s, matches largest factor root %s"
    % (mpmath.nstr(lam6, 12), abs(lam6 - perron[3]) < 1e-15)
)
print(
    "lambda_8 = 456 + 3 sqrt 11017 = %s, matches largest factor root %s"
    % (mpmath.nstr(lam8, 12), abs(lam8 - perron[4]) < 1e-15)
)
print("lambda_10 largest root of the quartic factor = %s" % mpmath.nstr(lam10, 13))


def reachable(A, start):
    seen, stack = {start}, [start]
    while stack:
        i = stack.pop()
        for j in range(len(A)):
            if A[i][j] and j not in seen:
                seen.add(j)
                stack.append(j)
    return seen


def scc_zero(M, zero):
    T = [[M[j][i] for j in range(len(M))] for i in range(len(M))]
    S = sorted(reachable(M, zero) & reachable(T, zero))
    sub = [[M[i][j] for j in S] for i in S]
    return S, sub, all(len(reachable(sub, i)) == len(S) for i in range(len(S)))


print(
    "energy cap: strongly connected component of the zero state, Perron root, E_2K(G_a) <= lambda_2K^a"
)
for K in (2, 3, 4, 5):
    M, zero = carry_matrix(K)
    radius = (K - 1) // 2
    spread = max(max(abs(d1), abs(d2)) for d1, d2 in delta_counts(K))
    closed = (radius + spread) // 3 <= radius
    S, sub, irreducible = scc_zero(M, zero)
    root = max(
        mpmath.mpf(str(r))
        for f, _ in sp.factor_list(sp.Matrix(sub).charpoly(x).as_expr())[1]
        for r in sp.Poly(f, x).nroots(n=25)
        if r.is_real
    )
    lam = mpmath.mpf(15) if K == 2 else perron[K]
    print(
        "2K = %2d box radius %d digit spread %d closed %s scc %2d of %2d irreducible %s self loop %s perron equals lambda_2K %s cap holds a = 0..6 %s ratio at a = 6 %s"
        % (
            2 * K,
            radius,
            spread,
            closed,
            len(S),
            len(M),
            irreducible,
            M[zero][zero] > 0,
            abs(root - lam) < mpmath.mpf(10) ** -18,
            all(mpmath.mpf(rows[K][a]) <= lam**a for a in range(7)),
            mpmath.nstr(mpmath.mpf(rows[K][6]) / lam**6, 6),
        )
    )

QUARTIC = sp.Poly(x**4 - 7833 * x**3 + 7916949 * x**2 - 850684437 * x + 13054946580, x)
LO, HI = sp.Rational(66641136625, 10**7), sp.Rational(66641136626, 10**7)
print(
    "lambda_10 certified by Sturm on the exact quartic: roots above %s %d, roots in the bracket %d, width %s"
    % (
        HI,
        QUARTIC.count_roots(HI, sp.oo),
        QUARTIC.count_roots(LO, HI),
        sp.nsimplify(HI - LO),
    )
)
KAP10 = 10 - mpmath.log(mpmath.mpf(HI.p) / HI.q) / mpmath.log(3)


def truncate(v, d):
    q = int(mpmath.floor(v * mpmath.mpf(10) ** d))
    return "%d.%0*d" % (q // 10**d, d, q % 10**d)


print(
    "certified kappa_10 >= %s, beta_0^(10) >= %s, both truncated down at 12 places"
    % (truncate(KAP10, 12), truncate(beta(KAP10), 12))
)
print("safe short edge, truncated down at 7 places: %s" % truncate(beta(KAP10), 7))

print("ladder rungs: moments, kappa_2K, beta_0^(2K)")
kappas = {}
for K, lam in ((2, 15), (3, lam6), (4, lam8), (5, lam10)):
    kappas[K] = 2 * K - mpmath.log(lam) / mpmath.log(3)
    print(
        "2K = %2d kappa %s beta %s"
        % (2 * K, mpmath.nstr(kappas[K], 13), mpmath.nstr(beta(kappas[K]), 12))
    )
print(
    "tenth rung above eighth by %s" % mpmath.nstr(beta(kappas[5]) - beta(kappas[4]), 9)
)
print(
    "Holder block exponents 2K = 6, 8, 10:",
    ", ".join(str(Fraction(2, 1) / (1 - Fraction(1, 2 * K))) for K in (3, 4, 5)),
)

print("rows above the tenth, floating Perron roots of exact matrices")
for K in range(6, 11):
    M, _ = carry_matrix(K)
    lam = max(np.linalg.eigvals(np.array(M, dtype=float)).real)
    kappa = 2 * K - log3(lam)
    print(
        "2K = %2d states %3d lambda %.6f beta %.9f wall - beta %.2e"
        % (2 * K, len(M), lam, beta(kappa), WALL - beta(kappa))
    )


def primes(lo, hi):
    return list(sp.primerange(lo, hi + 1))


def fourier(p, n):
    t = np.arange(p)
    e = np.exp(2j * np.pi * t / p)
    W = np.abs(1 + e[:, None] + e[None, :]) / 3
    F = np.ones((p, p))
    out = []
    for l in range(n):
        I = (t * pow(3, l, p)) % p
        F = F * W[np.ix_(I, I)]
        out.append(F.copy())
    return out


print("moment identities on (Z/p)^2, primes 5..199")
checks, worst = 0, 0.0
for p in primes(5, 199):
    F = fourier(p, 6)
    for a in range(1, 6):
        targets = []
        if 3**a <= p:
            targets.append((2, p * p * 3.0 ** (-a)))
        if 2 * 3**a <= p:
            targets.append((4, p * p * (5 / 27) ** a))
        for K in (3, 4, 5):
            if K * (3**a - 1) < p and a < len(rows[K]):
                targets.append((2 * K, p * p * rows[K][a] / 3 ** (2 * K * a)))
        for power, target in targets:
            worst = max(worst, abs((F[a - 1] ** power).sum() - target) / target)
            checks += 1
print("identities checked %d, worst relative error %.1e" % (checks, worst))


def half(p):
    a = 0
    while 2 * 3 ** (a + 1) <= p:
        a += 1
    return a


KAPS = {K: float(kappas[K]) for K in (3, 4, 5)}
EDGE10 = KAPS[5] / (2 - KAPPA + 2 * KAPS[5])


def master10(p, n):
    a = half(p)
    b = min(a - 1, n - 2 * a)
    if n < 2 * a or b < 0:
        return None
    return p * p * 3.0 ** (-((8 + KAPPA) * a + KAPS[5] * b) / 10)


def master(p, n):
    a = half(p)
    bounds = [float(p * p)]
    if 3 ** (n // 2) <= p:
        bounds.append(p * p * 3.0 ** (-(n // 2)))
    if n >= 4 * a:
        bounds.append(p * p * 3.0 ** (-KAPPA * a))
    if n >= 3 * a:
        bounds.append(p * p * 3.0 ** (-(1 + KAPPA) * a / 2))
    if n >= 2 * a:
        b4 = min(a, n - 2 * a)
        bounds.append(p * p * 3.0 ** (-a / 2 - KAPPA * (a + b4) / 4))
        b = min(a - 1, n - 2 * a)
        if b >= 0:
            for K in (3, 4, 5):
                bounds.append(
                    p * p * 3.0 ** (-((2 * K - 2 + KAPPA) * a + KAPS[K] * b) / (2 * K))
                )
                if n >= 2 * K * b:
                    bounds.append(p * p * 3.0 ** (-KAPS[K] * b))
    return min(bounds)


print(
    "master bound against exact L_n(p), primes 5..199, n = 2..24, the min over all orders and the order-10 block alone"
)
solo, solo_worst, solo_where, solo_cases = 0, 0.0, None, 0
violations, worst, where = 0, 0.0, None
for p in primes(5, 199):
    for n, F in enumerate(fourier(p, 24), start=1):
        if n < 2:
            continue
        only = master10(p, n)
        if only is not None:
            solo_cases += 1
            r10 = F.sum() / only
            if r10 > 1 + 1e-9:
                solo += 1
            if r10 > solo_worst:
                solo_worst, solo_where = r10, (p, n)
        ratio = F.sum() / master(p, n)
        if ratio > 1 + 1e-9:
            violations += 1
        if ratio > worst:
            worst, where = ratio, (p, n)
print(
    "min over all orders: violations %d, worst ratio L/bound %.4f at (p, n) = %s"
    % (violations, worst, where)
)
print(
    "order-10 block alone: %d cases, violations %d, worst ratio L/bound %.4f at (p, n) = %s"
    % (solo_cases, solo, solo_worst, solo_where)
)


def regime(p, n):
    a = half(p)
    if 2 * a > n:
        return None
    return "I" if n >= 4 * a else "II" if n >= 3 * a else "III"


def tail_bound(p, n):
    a = half(p)
    r = regime(p, n)
    if r == "I":
        return p * p * (5 / 27) ** a
    if r == "II":
        return p * p * (5 / 81) ** (a / 2)
    return p * p * 3.0 ** (-a / 2) * (5 / 27) ** ((n - a) / 4)


print(
    "dyadic tail: exponents (kappa - 1)/8 = %.4f, (2 + kappa)/4 = %.4f"
    % ((KAPPA - 1) / 8, (2 + KAPPA) / 4)
)
print(
    "tail constants 16/(kappa - 1) = %.2f, 10 (1 + kappa)/(kappa - 1) 2^((1 - kappa)/2) = %.2f, 8/(2 + kappa) = %.2f"
    % (
        16 / (KAPPA - 1),
        10 * (1 + KAPPA) / (KAPPA - 1) * 2 ** ((1 - KAPPA) / 2),
        8 / (2 + KAPPA),
    )
)
uncovered = [
    (p, n)
    for n in (8, 10, 12, 14, 16, 20, 24)
    for p in primes(5, 199)
    if p <= 3 ** (EDGE10 * n) and regime(p, n) is None
]
print(
    "regime cover below the order-10 edge 3^(%.12f n) for n in 8..24: uncovered %d"
    % (EDGE10, len(uncovered))
)
fails, cases = 0, 0
for n in (8, 10, 12):
    for p in primes(5, 199):
        if regime(p, n) is None:
            continue
        cases += 1
        R = fourier(p, n)[-1].sum() - 1
        if R > tail_bound(p, n) * (1 + 1e-9):
            fails += 1
print("per-prime tail bounds at n = 8, 10, 12: %d cases, %d failures" % (cases, fails))

LAMBDA = {K: 2 - KAPPA + 2 * KAPS[K] for K in (3, 4, 5)}
print(
    "master-bound constants per order: Lambda_2K = 2 - kappa + 2 kappa_2K, edge kappa_2K/Lambda_2K, decay Lambda_2K/2K, geometric constant 2/(1 - 3^(-Lambda_2K/2K))"
)
for K in (3, 4, 5):
    print(
        "2K = %2d Lambda %.9f edge %.12f decay %.9f constant %.4f"
        % (
            2 * K,
            LAMBDA[K],
            KAPS[K] / LAMBDA[K],
            LAMBDA[K] / (2 * K),
            2 / (1 - 3.0 ** (-LAMBDA[K] / (2 * K))),
        )
    )


def gain10(a, n):
    return (KAPS[5] * n - LAMBDA[5] * a) / 10


def four_term(n, z, eta):
    return (
        2 / z
        + 35 * z ** (1 - KAPPA)
        + 40 * 3.0 ** (-(KAPPA - 1) * n / 8)
        + 6 * 3.0 ** (-(LAMBDA[5] / 10) * eta * n)
    )


edge10 = EDGE10
bad, ratio = 0, 0.0
for eta in (0.001, 0.01, 0.05, 0.1):
    for n in range(6, 401):
        lo, hi = n // 3 + 1, int((edge10 - eta) * n)
        if hi < lo:
            continue
        block = [gain10(a, n) for a in range(lo, hi + 1)]
        total = sum(2 * 3.0 ** (-g) for g in block)
        cap = 6 * 3.0 ** (-(LAMBDA[5] / 10) * eta * n)
        if min(block) <= 0 or total > cap:
            bad += 1
        ratio = max(ratio, total / cap)
print(
    "order-10 main range, primes with 3a > n: gains positive and geometric cap holds at eta in 0.001..0.1, n = 6..400, failures %d, worst sum/cap %.4f"
    % (bad, ratio)
)


def energy_sum(p, n):
    t = np.arange(p)
    e = np.exp(2j * np.pi * t / p)
    W = np.abs(1 + e[:, None] + e[None, :]) / 3
    F = np.ones((p, p))
    for l in range(n):
        I = (t * pow(3, l, p)) % p
        F = F * W[np.ix_(I, I)]
    return F.sum()


print("four-term dyadic bound against the exact prime sum, eta = 0.02")
for n in (6, 8, 10, 12):
    for z in (5, 11):
        top = 3.0 ** ((edge10 - 0.02) * n)
        total = sum(energy_sum(p, n) / (p * p) for p in primes(z + 1, int(top)))
        print(
            "n = %2d z = %2d exact %.6f bound %.6f holds %s"
            % (n, z, total, four_term(n, z, 0.02), total <= four_term(n, z, 0.02))
        )
