from fractions import Fraction
from math import gcd, isqrt

import numpy as np
from sympy import Catalan, Poly, Symbol, cyclotomic_poly
from sympy import zeta as szeta

X = Symbol("x")
ZETA_K2 = float(szeta(2) * Catalan)
ZETA_2 = float(szeta(2))
CHECK_N = 50
CHECK_LAMBDA = 20
RESIDUE_N = 200
SAYOUS_TS = [2, 3, 4, 5, 6]
T2_NS = [20, 50]
T2_CUT = 200000
CLASSICAL_Q = 40
CLASSICAL_QS = [125, 250, 500, 1000, 2000, 4000, 8000]
METER_NS = [100, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000]

# GAUSSIAN ARITHMETIC


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gsub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def gnearest(z, w):
    p = gmul(z, gconj(w))
    n = norm(w)
    return ((2 * p[0] + n) // (2 * n), (2 * p[1] + n) // (2 * n))


def ggcd(z, w):
    while w != (0, 0):
        q = gnearest(z, w)
        z, w = w, gsub(z, gmul(q, w))
    return z


def gdivides(w, z):
    p = gmul(z, gconj(w))
    n = norm(w)
    return p[0] % n == 0 and p[1] % n == 0


def gquo(z, w):
    p = gmul(z, gconj(w))
    n = norm(w)
    return (p[0] // n, p[1] // n)


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


def elements_up_to(n):
    r = isqrt(n)
    out = []
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            if 0 < a * a + b * b <= n:
                out.append((a, b))
    return sorted(out, key=lambda z: (norm(z), z))


def is_rational_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def is_prime_class(z):
    n = norm(z)
    if is_rational_prime(n):
        return True
    r = isqrt(n)
    return r * r == n and r % 4 == 3 and is_rational_prime(r)


def factor_class(z):
    out = []
    cur = canon(z)
    while norm(cur) > 1:
        for p in classes_up_to(norm(cur)):
            if norm(p) > 1 and is_prime_class(p) and gdivides(p, cur):
                out.append(p)
                cur = canon(gquo(cur, p))
                break
        else:
            raise ValueError
    return out


def mobius_class(z):
    f = factor_class(z)
    return 0 if len(set(f)) != len(f) else (-1) ** len(f)


def totient_class(z):
    v = norm(z)
    if v == 1:
        return 1
    for p in set(factor_class(z)):
        v = v // norm(p) * (norm(p) - 1)
    return v


# RESIDUES MOD A GAUSSIAN INTEGER


def residue_box(d):
    g = abs(gcd(d[0], d[1]))
    return g, norm(d) // g


def residues(d):
    g, h = residue_box(d)
    return [(x, y) for x in range(g) for y in range(h)]


def check_residues(n_max):
    bad_count = 0
    bad_distinct = 0
    bad_totient = 0
    for d in classes_up_to(n_max):
        n = norm(d)
        r = residues(d)
        if len(r) != n:
            bad_count += 1
        reps = {tuple(v % n for v in gmul(u, gconj(d))) for u in r}
        if len(reps) != n:
            bad_distinct += 1
        if len(coprime_residues(d)) != totient_class(d):
            bad_totient += 1
    return len(classes_up_to(n_max)), bad_count, bad_distinct, bad_totient


def coprime_residues(d):
    return [u for u in residues(d) if norm(ggcd(u, d)) == 1]


# EXACT SUMS OF ROOTS OF UNITY


def roots_of_unity_sum(counts, n, phi):
    deg = len(phi) - 1
    c = list(counts)
    for k in range(n - 1, deg - 1, -1):
        f = c[k]
        if f:
            for j in range(deg + 1):
                c[k - j] -= f * phi[j]
    if any(c[1:deg]):
        return None
    return c[0]


# THEOREM 1


def ramanujan_exact(d, lam, phi_cache):
    n = norm(d)
    counts = [0] * n
    for u in coprime_residues(d):
        k = gmul(gmul(lam, u), gconj(d))[0] % n
        counts[k] += 1
    if n not in phi_cache:
        phi_cache[n] = [int(c) for c in Poly(cyclotomic_poly(n, X), X).all_coeffs()]
    return roots_of_unity_sum(counts, n, phi_cache[n])


def ramanujan_formula(d, lam, mob):
    tot = 0
    for e in classes_up_to(norm(d)):
        if gdivides(e, d) and (lam == (0, 0) or gdivides(e, lam)):
            tot += mob[canon(gquo(d, e))] * norm(e)
    return tot


def mertens_classes(t, mob):
    return sum(v for z, v in mob.items() if norm(z) <= t)


def sum_formula(n_max, lam, mob):
    tot = 0
    for e in classes_up_to(n_max if lam == (0, 0) else min(n_max, norm(lam))):
        if lam == (0, 0) or gdivides(e, lam):
            tot += norm(e) * mertens_classes(n_max // norm(e), mob)
    return tot


def node_set(n_max):
    out = []
    for d in classes_up_to(n_max):
        n = norm(d)
        for u in coprime_residues(d):
            w = gmul(u, gconj(d))
            out.append((Fraction(w[0] % n, n), Fraction(w[1] % n, n)))
    return out


def torus_point(w, n):
    x, y = w[0] % n, w[1] % n
    g = gcd(gcd(x, y), n)
    return (x // g, y // g, n // g)


def sayous_set(t):
    out = set()
    for q in elements_up_to(t * t):
        n = norm(q)
        cj = gconj(q)
        for x in range(n):
            for y in range(n):
                out.add(torus_point(gmul((x, y), cj), n))
    return out


def node_points(n_max):
    out = set()
    for x, y in node_set(n_max):
        n = x.denominator * y.denominator // gcd(x.denominator, y.denominator)
        out.add(
            torus_point(
                (
                    x.numerator * (n // x.denominator),
                    y.numerator * (n // y.denominator),
                ),
                n,
            )
        )
    return out


def check_sayous(ts):
    rows = []
    for t in ts:
        s = sayous_set(t)
        g = node_points(t * t)
        rows.append((t, t * t, len(s), len(g), s == g))
    return rows


def literal_sum(nodes, lam):
    p, q = lam
    x = np.array([float(a) for a, _ in nodes])
    y = np.array([float(b) for _, b in nodes])
    return complex(np.sum(np.exp(2j * np.pi * (p * x - q * y))))


def check_theorem_1(n_max, lam_bound):
    mob = {z: mobius_class(z) for z in classes_up_to(n_max)}
    phi_cache = {}
    nodes = node_set(n_max)
    bad_ram = 0
    bad_sum = 0
    bad_lit = 0
    worst = 0.0
    lams = elements_up_to(lam_bound)
    for lam in lams:
        tot = 0
        for d in classes_up_to(n_max):
            ex = ramanujan_exact(d, lam, phi_cache)
            fo = ramanujan_formula(d, lam, mob)
            if ex is None or ex != fo:
                bad_ram += 1
            tot += fo
        if tot != sum_formula(n_max, lam, mob):
            bad_sum += 1
        lit = literal_sum(nodes, lam)
        err = abs(lit - tot)
        worst = max(worst, err)
        if err > 1e-9:
            bad_lit += 1
    return len(nodes), len(lams), bad_ram, bad_sum, bad_lit, worst


# THEOREM 2


def gauss_weights(n_max, mob):
    out = []
    for e in classes_up_to(n_max):
        w = norm(e) * mertens_classes(n_max // norm(e), mob)
        if w:
            out.append((e, w))
    return out


def franel_exact(n_max, mob):
    cls = classes_up_to(n_max)
    mert = {e: mertens_classes(n_max // norm(e), mob) for e in cls}
    tot = Fraction(0)
    for a in cls:
        if mert[a] == 0:
            continue
        for b in cls:
            if mert[b] == 0:
                continue
            g = norm(canon(ggcd(a, b)))
            tot += Fraction(g * g * mert[a] * mert[b], norm(a) * norm(b))
    return tot


def lambda_side(n_max, cut, mob, a_arr):
    r = isqrt(cut)
    side = 2 * r + 1
    acc = np.zeros((side, side), dtype=np.int64)
    for e, w in gauss_weights(n_max, mob):
        lim = cut // norm(e)
        rr = isqrt(lim)
        s, t = np.meshgrid(
            np.arange(-rr, rr + 1), np.arange(-rr, rr + 1), indexing="ij"
        )
        keep = (s * s + t * t <= lim) & ((s != 0) | (t != 0))
        s, t = s[keep], t[keep]
        xs = e[0] * s - e[1] * t
        ys = e[0] * t + e[1] * s
        np.add.at(acc, (xs + r, ys + r), w)
    ii, jj = np.meshgrid(np.arange(-r, r + 1), np.arange(-r, r + 1), indexing="ij")
    nn = ii * ii + jj * jj
    keep = (nn > 0) & (nn <= cut)
    vals = acc[keep].astype(np.float64)
    nrm = nn[keep].astype(np.float64)
    lhs = float(np.sum(vals * vals / (nrm * nrm)))
    weight = float(np.sum(1.0 / (nrm * nrm)))
    return lhs, weight


def check_theorem_2(n_max, cut, mob, a_arr):
    m = sum(totient_class(d) for d in classes_up_to(n_max))
    exact = franel_exact(n_max, mob)
    rhs = 4.0 * ZETA_K2 * float(exact)
    lhs, weight = lambda_side(n_max, cut, mob, a_arr)
    tail = m * m * (4.0 * ZETA_K2 - weight)
    return m, exact, rhs, lhs, rhs - lhs, tail


# THE SIEVES


def chi4(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


def sieve(n_max, gaussian):
    a = np.zeros(n_max + 1, dtype=np.int64)
    if gaussian:
        for d in range(1, n_max + 1, 2):
            a[d::d] += chi4(d)
    else:
        a[1:] = 1
    mob = np.zeros(n_max + 1, dtype=np.int64)
    mob[1] = 1
    for n in range(1, n_max // 2 + 1):
        v = mob[n]
        if v:
            mob[2 * n :: n] -= a[2 : n_max // n + 1] * v
    j2 = np.zeros(n_max + 1, dtype=np.int64)
    ph = np.zeros(n_max + 1, dtype=np.int64)
    for d in range(1, n_max + 1):
        if a[d]:
            j2[d::d] += (d * d * a[d]) * mob[1 : n_max // d + 1]
            ph[d::d] += (d * a[d]) * mob[1 : n_max // d + 1]
    mg = np.zeros(n_max + 1, dtype=np.int64)
    mg[1:] = np.cumsum(mob[1:])
    return a, mob, mg, j2, np.cumsum(ph)


def inner_h(t, a, mg):
    ms = np.arange(1, t + 1)
    return float(np.sum(a[1 : t + 1] * mg[t // ms] / ms))


def franel_form(n, a, mg, j2):
    ns = np.arange(1, n + 1)
    ts = n // ns
    cache = {int(t): inner_h(int(t), a, mg) for t in np.unique(ts)}
    hv = np.array([cache[int(t)] for t in ts])
    return float(np.sum(j2[1 : n + 1] * hv * hv / (ns.astype(np.float64) ** 2)))


# THE CLASSICAL CONTROL


def farey_delta_square(q):
    nodes = sorted({Fraction(a, b) for b in range(1, q + 1) for a in range(1, b + 1)})
    m = len(nodes)
    return m, sum((r - Fraction(j + 1, m)) ** 2 for j, r in enumerate(nodes))


def classical_exact(q):
    mu = {1: 1}
    for n in range(2, q + 1):
        mu[n] = -sum(mu[d] for d in range(1, n) if n % d == 0)
    mert = {n: sum(mu[k] for k in range(1, q // n + 1)) for n in range(1, q + 1)}
    tot = Fraction(0)
    for x in range(1, q + 1):
        for y in range(1, q + 1):
            g = gcd(x, y)
            tot += Fraction(g * g * mert[x] * mert[y], x * y)
    return tot


# THE RUN


def readout(x, a, mg):
    ms = np.arange(1, x + 1)
    return int(np.sum(a[1 : x + 1] * mg[x // ms]))


def main():
    nodes, lams, bad_ram, bad_sum, bad_lit, worst = check_theorem_1(
        CHECK_N, CHECK_LAMBDA
    )
    ncls, bc, bd, bt = check_residues(RESIDUE_N)
    print(
        f"RESIDUE SYSTEMS to norm bound {RESIDUE_N}: {ncls} classes, {bc} wrong sizes, {bd} collisions mod d, {bt} totient mismatches"
    )
    print("SAYOUS SET")
    for t, n, ns, ng, same in check_sayous(SAYOUS_TS):
        print(
            f"  T = {t}, N = T^2 = {n}: literal G_T has {ns} points, node set has {ng}, equal is {same}"
        )
    print(f"zeta_K(2) = zeta(2) * Catalan = {ZETA_K2:.6f}")
    print("THEOREM 1")
    print(
        f"  nodes at N = {CHECK_N}: {nodes}; lambda with N(lambda) <= {CHECK_LAMBDA}: {lams}"
    )
    print(
        f"  exact Ramanujan mismatches: {bad_ram} of {lams * len(classes_up_to(CHECK_N))}"
    )
    print(f"  Mertens-form mismatches: {bad_sum} of {lams}")
    print(
        f"  literal node-sum mismatches at 1e-9: {bad_lit} of {lams}, worst {worst:.3e}"
    )

    a_big, _mob_big, mg_big, j2_big, ph_big = sieve(max(METER_NS), True)
    print("THEOREM 2")
    for n in T2_NS:
        mob = {z: mobius_class(z) for z in classes_up_to(n)}
        m, exact, rhs, lhs, diff, tail = check_theorem_2(n, T2_CUT, mob, a_big)
        ok = 0.0 <= diff <= tail
        print(f"  N = {n}: m = {m}, F(N) = {exact} = {float(exact):.6f}")
        print(
            f"    identity {rhs:.6f}, truncated Fourier side {lhs:.6f}, gap {diff:.6f}, tail bound {tail:.6f}, inside {ok}"
        )

    print("THE ZERO MODE")
    for n in T2_NS + [200]:
        mob = {z: mobius_class(z) for z in classes_up_to(n)}
        m = sum(totient_class(d) for d in classes_up_to(n))
        print(
            f"  N = {n}: sum Phi(d) = {m}, sum N(e) M_G(N/N(e)) = {sum_formula(n, (0, 0), mob)}"
        )
    print(
        f"  sum_(N(a) <= x) M_G(x/N(a)) over x = 1..2000: always 1 is {all(readout(x, a_big, mg_big) == 1 for x in range(1, 2001))}"
    )

    print("CLASSICAL CONTROL")
    m_cl, s2 = farey_delta_square(CLASSICAL_Q)
    c_q = classical_exact(CLASSICAL_Q)
    print(f"  Q = {CLASSICAL_Q}: m = {m_cl}, sum delta^2 = {float(s2):.10f}")
    print(f"  C(Q) - 1 = 12 m sum delta^2 exactly: {c_q - 1 == 12 * m_cl * s2}")
    a_cl, _mob_cl, mg_cl, j2_cl, ph_cl = sieve(max(METER_NS), False)
    print(
        f"  sieve C(Q) = {franel_form(CLASSICAL_Q, a_cl, mg_cl, j2_cl):.9f} against exact {float(c_q):.9f}"
    )

    print("  Q | Phi(Q) | S2 * Q from the identity")
    for q in CLASSICAL_QS:
        cq = franel_form(q, a_cl, mg_cl, j2_cl)
        print(
            f"  {q} | {int(ph_cl[q])} | {(cq - 1.0) * q / (12.0 * int(ph_cl[q])):.4f}"
        )

    print("THE METER")
    print(
        "  N | m | F(N) | F(N)/N | D_2(N)^2 * N^3 | slope | M_G(N) | M_G(N)^2 / (zeta_K(2) F(N)) | classical C(N)/N"
    )
    prev = None
    for n in METER_NS:
        f = franel_form(n, a_big, mg_big, j2_big)
        m = int(ph_big[n])
        d2sq = 4.0 * ZETA_K2 * f / (m * m)
        c = franel_form(n, a_cl, mg_cl, j2_cl)
        slope = (
            "-" if prev is None else f"{np.log(f / prev[1]) / np.log(n / prev[0]):.4f}"
        )
        mg = int(mg_big[n])
        print(
            f"  {n} | {m} | {f:.4f} | {f / n:.6f} | {d2sq * n**3:.4f} | {slope} | {mg} | {mg * mg / (ZETA_K2 * f):.6f} | {c / n:.6f}"
        )
        prev = (n, f)


if __name__ == "__main__":
    main()
