import sys
from math import gcd

# DIGITS


def dmask(n, q):
    m = 0
    while n:
        m |= 1 << (n % q)
        n //= q
    return m


def in_set(n, q, F):
    return n >= 1 and (dmask(n, q) & ~F) == 0


def repunit(q, L):
    return (q**L - 1) // (q - 1)


def show(F, q):
    return "".join(str(d) for d in range(q) if (F >> d) & 1)


# WALL


def theorem_witness(q, F):
    full = (1 << q) - 1
    if not (F >> 1) & 1:
        return ("unit", None, None)
    if F == full:
        return ("full", None, None)
    missing = [c for c in range(q) if not (F >> c) & 1]
    high = [c for c in missing if c >= 2]
    if high:
        c = min(high)
        return ("repunit", repunit(q, c), repunit(q, c + 1))
    if q % 2 == 1:
        return ("odd", 2, (q * q + 1) // 2)
    return ("even", q * q - 1, q * q + 1)


def breaks(q, F, m, n):
    if gcd(m, n) != 1 or m < 2 or n < 2:
        return False
    lhs = 1 if in_set(m * n, q, F) else 0
    rhs = (1 if in_set(m, q, F) else 0) * (1 if in_set(n, q, F) else 0)
    return lhs != rhs


def least_witness(q, F, pairs, masks):
    keep = ~F
    for m, n in pairs:
        a = 1 if (masks[m] & keep) == 0 else 0
        b = 1 if (masks[n] & keep) == 0 else 0
        c = 1 if (masks[m * n] & keep) == 0 else 0
        if c != a * b:
            return (m, n)
    return None


def coprime_pairs(bound):
    out = []
    for m in range(2, bound):
        for n in range(m + 1, bound // m + 1):
            if gcd(m, n) == 1:
                out.append((m, n))
    out.sort(key=lambda p: (p[0] * p[1], p[0]))
    return out


def wall(qmax=12, bound=4000):
    pairs = coprime_pairs(bound)
    print("WALL: is 1_(S_F) multiplicative")
    print("q  sets  unit  full  witnessed  kinds")
    total = {"unit": 0, "full": 0, "repunit": 0, "odd": 0, "even": 0}
    rows = []
    for q in range(2, qmax + 1):
        masks = [0] * (bound + 1)
        for n in range(1, bound + 1):
            masks[n] = dmask(n, q)
        kinds = {"unit": 0, "full": 0, "repunit": 0, "odd": 0, "even": 0}
        unfound = 0
        for F in range(1, 1 << q):
            kind, m, n = theorem_witness(q, F)
            kinds[kind] += 1
            total[kind] += 1
            if kind in ("unit", "full"):
                continue
            assert breaks(q, F, m, n), (q, F, kind, m, n)
            lw = least_witness(q, F, pairs, masks)
            if lw is None:
                unfound += 1
                print("  no witness below %d at q %d F {%s}" % (bound, q, show(F, q)))
            else:
                rows.append((q, F, kind, m, n, lw[0], lw[1]))
        print(
            "%2d %5d %5d %5d %10d  repunit %d odd %d even %d unfound %d"
            % (
                q,
                (1 << q) - 1,
                kinds["unit"],
                kinds["full"],
                kinds["repunit"] + kinds["odd"] + kinds["even"],
                kinds["repunit"],
                kinds["odd"],
                kinds["even"],
                unfound,
            )
        )
    hardest = max(rows, key=lambda r: r[5] * r[6])
    print(
        "hardest least witness: q %d F {%s} pair (%d, %d) product %d"
        % (
            hardest[0],
            show(hardest[1], hardest[0]),
            hardest[5],
            hardest[6],
            hardest[5] * hardest[6],
        )
    )
    print(
        "totals unit %d full %d repunit %d odd %d even %d"
        % (total["unit"], total["full"], total["repunit"], total["odd"], total["even"])
    )
    print()
    print("least witness, every F with 1 in F and F not full, q <= 5")
    print("q F         kind     theorem pair        least pair   product")
    for q, F, kind, m, n, a, b in rows:
        if q > 5:
            continue
        if not (F >> 1) & 1:
            continue
        print(
            "%d {%-8s} %-8s (%d, %d)%s(%d, %d) %s= %d"
            % (
                q,
                show(F, q),
                kind,
                m,
                n,
                " " * max(1, 16 - len("(%d, %d)" % (m, n))),
                a,
                b,
                " " * max(1, 12 - len("(%d, %d)" % (a, b))),
                a * b,
            )
        )


VERBS = {"wall": wall}

# LERCH


def lerch_table(s, terms):
    from mpmath import factorial, gamma, zeta

    return ([zeta(s - j) / factorial(j) for j in range(terms)], gamma(1 - s))


def lerch(s, t, tab):
    from mpmath import j as I
    from mpmath import mp, mpf, pi

    zs, g = tab
    t = mpf(t)
    if t <= mpf(1) / 2:
        mu = -2 * pi * I * t
        head = g * (2 * pi * I * t) ** (s - 1)
    else:
        u = 1 - t
        mu = 2 * pi * I * u
        head = g * (-2 * pi * I * u) ** (s - 1)
    acc = mp.mpc(0)
    p = mp.mpc(1)
    for c in zs:
        acc += c * p
        p *= mu
    return head + acc


def gl_nodes(n):
    import numpy
    from mpmath import mpf

    x, w = numpy.polynomial.legendre.leggauss(n)
    return [mpf(float(a)) for a in x], [mpf(float(a)) for a in w]


def design_sum(q, F, L, s):
    from mpmath import mp

    digs = [d for d in range(q) if (F >> d) & 1]
    vals = [0]
    for i in range(L):
        p = q**i
        vals = [v + d * p for v in vals for d in digs]
    return mp.fsum([mp.mpf(n) ** (-s) for n in vals if n >= 1]), len(vals)


def gseq(q, F, L, t):
    from mpmath import j as I
    from mpmath import mp, pi

    out = mp.mpc(1)
    for i in range(L):
        x = t * (q**i)
        x = x - mp.floor(x)
        w = mp.exp(2 * pi * I * x)
        p = mp.mpc(1)
        acc = mp.mpc(0)
        for d in range(q):
            if (F >> d) & 1:
                acc += p
            p *= w
        out *= acc
    return out


def cell_integral(q, F, L, s, tab, nodes, weights, lo, hi):
    from mpmath import mp

    half = (hi - lo) / 2
    mid = (hi + lo) / 2
    acc = mp.mpc(0)
    for x, w in zip(nodes, weights):
        t = mid + half * x
        acc += w * gseq(q, F, L, t) * lerch(s, t, tab)
    return acc * half


def position_integral(q, F, L, s, npts=14, terms=80, grade=64):
    from mpmath import mp

    tab = lerch_table(s, terms)
    nodes, weights = gl_nodes(npts)
    N = q**L
    h = mp.mpf(1) / N
    total = mp.mpc(0)
    for a in range(1, N - 1):
        total += cell_integral(q, F, L, s, tab, nodes, weights, a * h, (a + 1) * h)
    for j in range(grade):
        lo = h / (2 ** (j + 1))
        hi = h / (2**j)
        total += cell_integral(q, F, L, s, tab, nodes, weights, lo, hi)
        total += cell_integral(q, F, L, s, tab, nodes, weights, 1 - hi, 1 - lo)
    return total


def position(qmax=None):
    from mpmath import exp, mp, mpc, pi, polylog
    from mpmath import j as I

    mp.dps = 25
    print("LERCH EVALUATOR against mpmath polylog")
    for s in [mpc(3, 3) + mpc("0.3"), mpc("2.7", "1.9")]:
        tab = lerch_table(s, 80)
        for t in ["0.03125", "0.25", "0.4", "0.6", "0.9"]:
            a = lerch(s, mp.mpf(t), tab)
            b = polylog(s, exp(-2 * pi * I * mp.mpf(t)))
            print("  s %s t %s  err %s" % (mp.nstr(s, 6), t, mp.nstr(abs(a - b), 3)))
    print()
    print("POSITION PRODUCT: sum over D_L of n^(-s) against the integral")
    print(
        "q F        L  k^L    s              direct                     integral                   err"
    )
    cases = [(10, 3), (10, 4), (3, 3), (3, 4), (3, 5)]
    for q, L in cases:
        F = ((1 << q) - 1) & ~(1 << (q - 1)) if q == 10 else 0b011
        for s in [mpc("3.3", 0), mpc("2.7", "1.9")]:
            direct, kL = design_sum(q, F, L, s)
            integ = position_integral(q, F, L, s)
            print(
                "%2d %-8s %d %6d %-14s %-26s %-26s %s"
                % (
                    q,
                    show(F, q),
                    L,
                    kL,
                    mp.nstr(s, 5),
                    mp.nstr(direct, 14),
                    mp.nstr(integ, 14),
                    mp.nstr(abs(direct - integ), 3),
                )
            )


VERBS["position"] = position

# PAIR


def mobius_sieve(N):
    mu = [1] * (N + 1)
    primes = []
    comp = [False] * (N + 1)
    for i in range(2, N + 1):
        if not comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > N:
                break
            comp[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    return mu


def pair_coeff(n, q, F, mu, divs, masks):
    keep = ~F
    tot = 0
    for d in divs[n]:
        e = n // d
        if (masks[d] & keep) == 0 and (masks[e] & keep) == 0:
            tot += mu[e]
    return tot


def pair(N=4000):
    mu = mobius_sieve(N)
    divs = [[] for _ in range(N + 1)]
    for d in range(1, N + 1):
        for m in range(d, N + 1, d):
            divs[m].append(d)
    print("ZETA_F TIMES M_F IS NOT 1")
    print("coefficient c(n) = sum over a b = n with a, b in S_F of mu(b)")
    print("q F          k  first n > 1 with c(n) nonzero   c(n)")
    cases = []
    for q in range(2, 9):
        for F in range(1, 1 << q):
            if (F >> 1) & 1:
                cases.append((q, F))
    for f in [((1 << 10) - 1) & ~(1 << 9), ((1 << 10) - 1) & ~1, (1 << 10) - 1]:
        cases.append((10, f))
    shown = {
        (3, 0b011),
        (3, 0b101),
        (3, 0b110),
        (4, 0b0011),
        (5, 0b00011),
        (10, ((1 << 10) - 1) & ~(1 << 9)),
        (10, ((1 << 10) - 1) & ~1),
        (3, 0b111),
        (10, (1 << 10) - 1),
    }
    masks_by_q = {}
    firsts = []
    for q, F in cases:
        if q not in masks_by_q:
            masks_by_q[q] = [0] + [dmask(n, q) for n in range(1, N + 1)]
        masks = masks_by_q[q]
        assert pair_coeff(1, q, F, mu, divs, masks) == 1
        hit = None
        for n in range(2, N + 1):
            c = pair_coeff(n, q, F, mu, divs, masks)
            if c != 0:
                hit = (n, c)
                break
        firsts.append((q, F, hit))
        if (q, F) in shown:
            print(
                "%2d %-10s %2d  %-30s %s"
                % (
                    q,
                    show(F, q),
                    F.bit_count(),
                    "none below %d" % N if hit is None else str(hit[0]),
                    "-" if hit is None else str(hit[1]),
                )
            )
    full = [(q, F, h) for (q, F, h) in firsts if F == (1 << q) - 1]
    part = [(q, F, h) for (q, F, h) in firsts if F != (1 << q) - 1]
    print(
        "full digit sets tested %d, all with c(n) = 0 for 1 < n <= %d: %s"
        % (len(full), N, all(h is None for (_, _, h) in full))
    )
    print(
        "proper sets tested %d, all with a nonzero c(n): %s, largest first witness %d"
        % (
            len(part),
            all(h is not None for (_, _, h) in part),
            max(h[0] for (_, _, h) in part),
        )
    )


VERBS["pair"] = pair

# WORD


def mu_int(n):
    r = 1
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0:
                return 0
            r = -r
        d += 1
    if m > 1:
        r = -r
    return r


def lyndon(k, L):
    t = 0
    for d in range(1, L + 1):
        if L % d == 0:
            t += mu_int(d) * k ** (L // d)
    return t // L


def series_mul(a, b, n):
    out = [0] * n
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if i + j >= n:
                break
            out[i + j] += x * y
    return out


def binom(n, r):
    t = 1
    for i in range(r):
        t = t * (n - i) // (i + 1)
    return t


def word(order=17):
    print("WORD EULER PRODUCT: 1/(1 - k u) = prod over L of (1 - u^L)^(-c_k(L))")
    print("c_k(L) is the Lyndon count (1/L) sum over d dividing L of mu(d) k^(L/d)")
    print("k   c_k(1..10)")
    for k in [2, 3, 4, 9, 10]:
        print("%2d  %s" % (k, " ".join(str(lyndon(k, L)) for L in range(1, 11))))
    print()
    print("k   product to u^%d equals 1/(1 - k u)" % (order - 1))
    for k in [2, 3, 4, 9, 10]:
        acc = [0] * order
        acc[0] = 1
        for L in range(1, order):
            c = lyndon(k, L)
            f = [0] * order
            j = 0
            while L * j < order:
                f[L * j] = binom(c + j - 1, j)
                j += 1
            acc = series_mul(acc, f, order)
        print("%2d  %s" % (k, acc == [k**i for i in range(order)]))
    print()
    print("word Mobius: 1 - k u has coefficients 1, -k and zero after, so the word")
    print("Mertens is 1 at norm 1 and 1 - k at every norm above, bounded in the norm")
    print("word zeta 1/(1 - k q^(-s)) has no zero and simple poles exactly at")
    print("s = alpha + 2 pi i m / log q with alpha = log_q k, the design pole lattice")


VERBS["word"] = word

# CHARACTERS


def prime_factors(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def primitive_root(p):
    if p == 2:
        return 1
    fs = list(prime_factors(p - 1))
    g = 2
    while True:
        if all(pow(g, (p - 1) // f, p) != 1 for f in fs):
            return g
        g += 1


def crt_lift(x, m, Q):
    o = Q // m
    if o == 1:
        return x % Q
    inv = pow(o % m, -1, m)
    return (1 + o * ((x - 1) * inv % m)) % Q


def cyclic_parts(Q):
    parts = []
    for p, e in prime_factors(Q).items():
        m = p**e
        if p == 2:
            if e == 1:
                continue
            if e == 2:
                parts.append((crt_lift(3, 4, Q), 2))
            else:
                parts.append((crt_lift(m - 1, m, Q), 2))
                parts.append((crt_lift(5, m, Q), 1 << (e - 2)))
        else:
            g = primitive_root(p)
            if e > 1 and pow(g, p - 1, p * p) == 1:
                g += p
            parts.append((crt_lift(g, m, Q), (p - 1) * p ** (e - 1)))
    return parts


def char_table(Q):
    import cmath

    parts = cyclic_parts(Q)
    orders = [n for _, n in parts]
    logs = {}

    def walk(i, cur, exps):
        if i == len(parts):
            logs[cur] = tuple(exps)
            return
        g, n = parts[i]
        x = cur
        for a in range(n):
            walk(i + 1, x, exps + [a])
            x = x * g % Q

    walk(0, 1 % Q, [])
    idx = []

    def tuples(i, cur):
        if i == len(orders):
            idx.append(tuple(cur))
            return
        for a in range(orders[i]):
            tuples(i + 1, cur + [a])

    tuples(0, [])
    chars = []
    for a in idx:
        tab = [0j] * Q
        for r, l in logs.items():
            ph = sum(a[i] * l[i] / orders[i] for i in range(len(orders)))
            tab[r] = cmath.exp(2j * cmath.pi * ph)
        chars.append(tab)
    return chars


# FIBRE


def fibre(N=3000):
    import cmath
    from math import gcd as G

    mu = mobius_sieve(N)
    print("LERCH-MOBIUS AT A RATIONAL: M(s, a/Q) as inverse L-functions")
    print("M(s, a/Q) = sum over g dividing Q of mu(g) g^(-s) (1/phi(Q/g))")
    print("            times sum over chi mod Q/g of Gauss(chi, a) A_chi(s)")
    print(
        "A_chi(s) = 1/L(s, chi) times prod over p dividing Q not Q/g of (1 - chi(p) p^(-s))^(-1)"
    )
    print(
        "Q    a  divisors used  characters  max coefficient error to n = %d  Euler step"
        % N
    )
    for Q, a in [
        (3, 1),
        (9, 1),
        (9, 2),
        (27, 5),
        (5, 2),
        (10, 3),
        (10, 7),
        (100, 21),
        (7, 1),
        (8, 3),
        (12, 5),
    ]:
        lhs = [0j] * (N + 1)
        for n in range(1, N + 1):
            lhs[n] = mu[n] * cmath.exp(-2j * cmath.pi * n * a / Q)
        rhs = [0j] * (N + 1)
        gs = [g for g in range(1, Q + 1) if Q % g == 0 and mu[g] != 0]
        nchars = 0
        euler_err = 0.0
        for g in gs:
            Qp = Q // g
            chars = char_table(Qp) if Qp > 1 else [[1j * 0 + 1]]
            nchars += len(chars)
            phi = len([r for r in range(Qp) if G(r, Qp) == 1]) if Qp > 1 else 1
            extra = [p for p in prime_factors(Q) if Qp % p != 0]
            for chi in chars:
                cv = (lambda m: chi[m % Qp]) if Qp > 1 else (lambda m: 1.0 + 0j)
                gauss = (
                    sum(
                        (cv(r).conjugate()) * cmath.exp(-2j * cmath.pi * r * a / Qp)
                        for r in range(Qp)
                        if G(r, Qp) == 1
                    )
                    if Qp > 1
                    else 1.0 + 0j
                )
                A = [0j] * (N // g + 1)
                for m in range(1, N // g + 1):
                    if G(m, Q) == 1:
                        A[m] = mu[m] * cv(m)
                for m in range(1, N // g + 1):
                    if A[m] != 0j:
                        rhs[g * m] += mu[g] * gauss * A[m] / phi
                B = [0j] * (N // g + 1)
                for m in range(1, N // g + 1):
                    B[m] = mu[m] * cv(m)
                S = [0j] * (N // g + 1)
                S[1] = 1.0 + 0j
                for p in extra:
                    T = list(S)
                    pk = p
                    while pk <= N // g:
                        for m in range(1, N // g // pk + 1):
                            T[m * pk] += S[m] * cv(pk)
                        pk *= p
                    S = T
                C = [0j] * (N // g + 1)
                for m in range(1, N // g + 1):
                    if B[m] != 0j:
                        for l in range(1, N // g // m + 1):
                            if S[l] != 0j:
                                C[m * l] += B[m] * S[l]
                euler_err = max(
                    euler_err, max(abs(C[m] - A[m]) for m in range(1, N // g + 1))
                )
        err = max(abs(lhs[n] - rhs[n]) for n in range(1, N + 1))
        print(
            "%-4d %-2d %-14s %-11d %-33s %s"
            % (
                Q,
                a,
                ",".join(str(g) for g in gs),
                nchars,
                "%.3e" % err,
                "%.3e" % euler_err,
            )
        )
    print()
    print(
        "t = 0 fibre: M(s, 0) = 1/zeta(s), so the classical RH is one fibre of the family"
    )


VERBS["fibre"] = fibre

# BEURLING


def primes_to(X):
    sieve = bytearray([1]) * (X + 1)
    sieve[0] = sieve[1] = 0
    i = 2
    while i * i <= X:
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
        i += 1
    return [i for i in range(2, X + 1) if sieve[i]]


def generated(P, X):
    out = [(1, 1)]

    def go(i, val, sign, sqf):
        for j in range(i, len(P)):
            p = P[j]
            if val * p > X:
                break
            v = val * p
            out.append((v, -sign if sqf else 0))
            go(j + 1, v, -sign, sqf)
            w = v * p
            while w <= X:
                out.append((w, 0))
                go(j + 1, w, 0, False)
                w *= p

    go(0, 1, 1, True)
    out.sort()
    return out


def beurling(X=10**6):
    from math import log

    print("BEURLING SYSTEM ON A DESIGN: N_B is the free semigroup on the primes of S_F")
    print(
        "N_B is not S_F. mu_B is the restriction of mu, M_B(x) = sum of mu over N_B below x"
    )
    allp = primes_to(X)
    cases = [
        (10, ((1 << 10) - 1) & ~(1 << 9), "base 10 missing 9"),
        (3, 0b011, "base 3 {0,1}"),
        (3, 0b101, "base 3 {0,2}"),
        (10, (1 << 10) - 1, "base 10 full, control"),
    ]
    for q, F, name in cases:
        k = F.bit_count()
        alpha = log(k) / log(q)
        P = [p for p in allp if in_set(p, q, F)]
        gen = generated(P, X)
        print()
        print(
            "%s  k %d  alpha %.6f  alpha/2 %.6f  primes in S_F below %d: %d"
            % (name, k, alpha, alpha / 2, X, len(P))
        )
        print(
            "  x            pi_F(x)   N_B(x)     M_B(x)    max abs M_B   log max / log x  log max / log N_B  log_q x"
        )
        idx = 0
        run = 0
        mx = 0
        rows = []
        checks = []
        j = 2
        while True:
            for r in range(4):
                x = round(q ** (j + r / 4.0))
                if x > X:
                    break
                if x >= 10:
                    checks.append(x)
            if q**j > X:
                break
            j += 1
        checks = sorted({c for c in checks if c <= X})
        pi_at = 0
        pidx = 0
        for x in checks:
            while idx < len(gen) and gen[idx][0] <= x:
                run += gen[idx][1]
                mx = max(mx, abs(run))
                idx += 1
            while pidx < len(P) and P[pidx] <= x:
                pidx += 1
            pi_at = pidx
            nf = idx
            e = log(mx) / log(x) if mx > 0 else float("nan")
            e2 = log(mx) / log(nf) if mx > 0 and nf > 1 else float("nan")
            rows.append((log(x) / log(q), e, e2))
            print(
                "  %-12d %-9d %-10d %-9d %-13d %-16.6f %-18.6f %.4f"
                % (x, pi_at, nf, run, mx, e, e2, log(x) / log(q))
            )
        print(
            "  exponent log max / log x by residue class of log_q x, last four checkpoints each"
        )
        for r in [0.0, 0.25, 0.5, 0.75]:
            cls = [
                row
                for row in rows
                if abs((row[0] % 1.0) - r) < 0.02 or abs((row[0] % 1.0) - r - 1) < 0.02
            ]
            if cls:
                print(
                    "    r = %.2f  %s" % (r, "  ".join("%.4f" % c[1] for c in cls[-4:]))
                )


VERBS["beurling"] = beurling

# DUAL


def dual():
    from mpmath import exp, gamma, mp, mpc, mpf, pi, zeta
    from mpmath import j as I

    mp.dps = 30
    print(
        "REFLECTION: Z(s, t) from the Hurwitz formula, DLMF 25.13.3 solved for F(-t, s)"
    )
    print("Z(s, t) = (2 pi)^s Gamma(1-s) / (2 pi i) times")
    print("          [ e^(pi i s/2) zeta(1-s, t) - e^(-pi i s/2) zeta(1-s, 1-t) ]")
    print(
        "s              t         series value               reflection value           err"
    )
    for s in [mpc("3.3", 0), mpc("2.7", "1.9"), mpc("0.6", "4.1")]:
        tab = lerch_table(s, 120)
        for t in ["0.125", "0.3", "0.5", "0.77"]:
            t = mpf(t)
            a = lerch(s, t, tab)
            b = ((2 * pi) ** s * gamma(1 - s) / (2 * pi * I)) * (
                exp(pi * I * s / 2) * zeta(1 - s, t)
                - exp(-pi * I * s / 2) * zeta(1 - s, 1 - t)
            )
            print(
                "%-14s %-9s %-26s %-26s %s"
                % (
                    mp.nstr(s, 5),
                    mp.nstr(t, 4),
                    mp.nstr(a, 14),
                    mp.nstr(b, 14),
                    mp.nstr(abs(a - b), 3),
                )
            )
    print()
    print("the reflection moves the kernel and not the design measure, so the position")
    print("identity returns a dual integral of the same G_L against zeta(1-s, t) and")
    print("zeta(1-s, 1-t), never a relation between zeta_F(s) and zeta_F(1-s)")


VERBS["dual"] = dual

if __name__ == "__main__":
    v = sys.argv[1] if len(sys.argv) > 1 else "wall"
    VERBS[v]()
