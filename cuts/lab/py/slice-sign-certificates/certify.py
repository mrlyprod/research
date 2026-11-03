import sys
import time
from itertools import pairwise

import sympy
from mpmath import iv, mp

iv.dps = 80
mp.dps = 80
PI = iv.pi
TOP = 70


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def symbol(q):
    m = (q - 1) // 2
    return [0 if d == m else 1 for d in range(q)]


class Machine:
    def __init__(self, q, D, power=None):
        self.q, self.D, self.m = q, D, (q - 1) // 2
        a = symbol(q)
        if power is None:
            power = [1]
            for _ in range(D - 1):
                power = polymul(power, a)
        tip = a[:]
        tip[self.m] += D
        self.P = polymul(power, tip)
        self.fill = sum(self.P)
        assert self.fill == (q - 1) ** (D - 1) * (q - 1 + D)
        assert self.P == self.P[::-1] and len(self.P) == 2 * D * self.m + 1
        self.r = (D - 1) // 2
        n = self.r + 1
        self.v = [q - 1 if c % q == 0 else -1 for c in range(n)]
        e = self.entry
        self.W = [
            [e(0, c) if k == 0 else e(k, c) + e(-k, c) for k in range(n)]
            for c in range(n)
        ]
        self.even = [
            [e(k, 0) if c == 0 else e(k, c) + e(k, -c) for c in range(n)]
            for k in range(n)
        ]
        for c in range(-self.r, self.r + 1):
            assert all(e(k, c) == 0 for k in (self.r + 1, -self.r - 1))

    def entry(self, new, old):
        i = old + self.D * self.m - self.q * new
        return self.P[i] if 0 <= i < len(self.P) else 0

    def step(self, x):
        return [sum(w * y for w, y in zip(row, x)) for row in self.W]


def machines(q, dims):
    a = symbol(q)
    power, have = [1], 0
    for D in dims:
        while have < D - 1:
            power = polymul(power, a)
            have += 1
        yield Machine(q, D, power)


def depth(M, K):
    beta = [1] * (M.r + 1)
    for _ in range(K):
        beta = M.step(beta)
    nxt = M.step(beta)
    assert min(beta) > 0
    bad = [c for c in range(M.r + 1) if not M.q * nxt[c] < M.fill * beta[c]]
    return bad, beta, nxt


def kmin(M, cap):
    beta = [1] * (M.r + 1)
    for K in range(cap + 1):
        nxt = M.step(beta)
        assert min(beta) > 0
        if all(M.q * b1 < M.fill * b0 for b0, b1 in zip(beta, nxt)):
            return K
        beta = nxt
    raise SystemExit(f"no certificate at q = {M.q}, D = {M.D} up to depth {cap}")


def row_depth(M, cap):
    row, prefix = M.v[:], []
    for t in range(1, cap + 1):
        prefix.append(row[0])
        row = M.step(row)
        if min(row) >= 0:
            assert min(prefix) >= 0
            return t
    raise SystemExit(f"no row certificate at q = {M.q}, D = {M.D}")


def tower(M, top):
    row, out = M.v[:], []
    for _ in range(top + 1):
        out.append(row[0])
        row = M.step(row)
    return out


def bareiss(A):
    A = [row[:] for row in A]
    n, sign, prev = len(A), 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if A[i][k]), None)
            if swap is None:
                return 0
            A[k], A[swap] = A[swap], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[-1][-1]


def detmod(A, p):
    A = [[x % p for x in row] for row in A]
    n, det = len(A), 1
    for k in range(n):
        piv = next((i for i in range(k, n) if A[i][k]), None)
        if piv is None:
            return 0
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            det = -det
        det = det * A[k][k] % p
        inv = pow(A[k][k], -1, p)
        for i in range(k + 1, n):
            f = A[i][k] * inv % p
            if f:
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[k])]
    return det % p


def v2det(A, bits=256):
    prec = bits
    A = [[x % (1 << prec) for x in row] for row in A]
    n, total = len(A), 0
    for k in range(n):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                x = A[i][j]
                if x:
                    e = (x & -x).bit_length() - 1
                    if best is None or e < best[0]:
                        best = (e, i, j)
        assert best is not None
        e, i, j = best
        A[k], A[i] = A[i], A[k]
        for row in A:
            row[k], row[j] = row[j], row[k]
        total += e
        prec -= e
        mod = 1 << prec
        inv = pow(A[k][k] >> e, -1, mod)
        for i in range(k + 1, n):
            if A[i][k]:
                f = ((A[i][k] >> e) * inv) % mod
                A[i] = [(x - f * y) % mod for x, y in zip(A[i], A[k])]
        A = [[x % mod for x in row] for row in A]
    assert total < prec
    return total


def v2(x):
    return (x & -x).bit_length() - 1


def charpoly(A):
    n = len(A)
    c = [1] + [0] * n
    Mk = [[0] * n for _ in range(n)]
    for k in range(1, n + 1):
        Mk = [
            [sum(A[i][t] * Mk[t][j] for t in range(n)) for j in range(n)]
            for i in range(n)
        ]
        for i in range(n):
            Mk[i][i] += c[k - 1]
        tr = sum(A[i][t] * Mk[t][i] for i in range(n) for t in range(n))
        assert tr % k == 0
        c[k] = -tr // k
    return c


def above(M):
    x = sympy.symbols("x")
    chi = sympy.Poly(charpoly(M.even), x)
    t = sympy.Rational(M.fill, M.q)
    assert chi.eval(t) != 0
    return chi.count_roots(t, None)


def pencil(M):
    n = M.r + 1
    return [
        [(M.fill if i == j else 0) - M.q * M.even[i][j] for j in range(n)]
        for i in range(n)
    ]


def lt(x, y):
    return (x < y) is True


def up(x):
    return mp.make_mpf(x._mpi_[1])


def lo(x):
    return mp.make_mpf(x._mpi_[0])


def show(x, digits=8):
    return mp.nstr(up(x), digits)


def base3_constants():
    tg = {i: 2 * iv.cos(PI / 3**i) for i in range(2, TOP + 1)}
    hg = {i: 2 * iv.cos(2 * PI / 3**i) for i in range(2, TOP + 1)}
    ell = {i: iv.log(tg[i] / hg[i]) for i in range(2, TOP + 1)}
    for i in range(3, TOP + 1):
        assert lt(ell[i], iv.mpf("2.04") * (PI / 3**i) ** 2)
    tail = iv.mpf("2.04") * PI**2 / (8 * iv.mpf(9) ** TOP)
    logR = sum(ell.values()) + iv.mpf([0, up(tail)])
    R = iv.exp(logR)
    exit3 = iv.cos(7 * PI / 27) / iv.cos(2 * PI / 27)
    half3 = iv.cos(8 * PI / 27) / iv.cos(2 * PI / 27)
    true = {
        "0.7528157": (tg[2] / hg[2]) * half3,
        "0.7052518": exit3,
        "0.6137010": half3,
        "0.2266816": iv.cos(4 * PI / 9) / iv.cos(2 * PI / 9),
        "0.8900159": iv.sqrt(3) / (2 * iv.cos(2 * PI / 27)),
        "0.66966": (tg[2] / hg[2])
        * (tg[3] / hg[3])
        * iv.cos(PI / 3 - PI / 81)
        / iv.cos(2 * PI / 81),
        "0.6419": R * iv.cos(PI / 3 - PI / 243) / iv.cos(2 * PI / 243),
    }
    for printed, value in true.items():
        assert lt(value, iv.mpf(printed)), printed
    return tg, hg, ell, logR, R, true


def base3_tails(tg, hg, ell):
    worst_l, worst_s = mp.mpf(0), mp.mpf(0)
    tail_l = iv.mpf("2.04") * PI**2 / (8 * iv.mpf(9) ** TOP)
    tail_s = 5 * PI**2 / (8 * iv.mpf(9) ** TOP)
    for K in range(1, 8):
        LT = sum(ell[i] for i in range(K + 2, TOP + 1)) + iv.mpf([0, up(tail_l)])
        ST = sum(4 - tg[i] - hg[i] for i in range(K + 2, TOP + 1)) + iv.mpf(
            [0, up(tail_s)]
        )
        a, b = LT * iv.mpf(9) ** K, ST * iv.mpf(9) ** K
        assert lt(a, iv.mpf("0.3579")) and lt(b, iv.mpf("0.6919"))
        worst_l, worst_s = max(worst_l, up(a)), max(worst_s, up(b))
    return worst_l, worst_s


def s_of(D, tg, hg):
    total = sum(
        iv.log((D - tg[i]) / (D - 2)) - iv.log((D + hg[i]) / (D + 2))
        for i in range(2, TOP + 1)
    )
    tail = 5 * PI**2 / (8 * (D - 2) * iv.mpf(9) ** TOP)
    return total + iv.mpf([0, up(tail)])


def envelope(D, K):
    a, b, c, s = (
        iv.mpf(x) for x in ("0.7528157", "0.7052518", "0.2266816", "0.8900159")
    )
    return (4 * (K - 1) * (a ** (D - 1) + b ** (D - 1)) + 2 * c ** (D - 1)) * iv.exp(
        2 * (K + 1) * s ** (D - 1)
    )


def verb_base3():
    start = time.time()
    tg, hg, ell, logR, R, true = base3_constants()
    print(
        f"R = {show(R, 12)}, log R / 4 = {show(logR / 4, 12)}, log R in [{mp.nstr(lo(logR), 12)}, {mp.nstr(up(logR), 12)}]"
    )
    for printed, value in true.items():
        print(f"  constant {printed}: true value at most {show(value, 10)}")
    worst_l, worst_s = base3_tails(tg, hg, ell)
    print(
        f"Lemma 7.4 tails at K = 1..7: 9^K LT(K) <= {mp.nstr(worst_l, 5)} <= 0.3579, 9^K (D-2) ST(D, K) <= {mp.nstr(worst_s, 5)} <= 0.6919"
    )
    rows, widest, window_min = [], mp.mpf(0), None
    for D in range(38, 401, 2):
        sD = s_of(D, tg, hg)
        A = iv.log(iv.mpf(D + 2) / (D - 2))
        X = ((D - 1) * logR + sD) / A
        K1 = int(mp.floor(up(X))) + 2
        assert lt(X + 1, iv.mpf(K1))
        widest = max(widest, (up(X) - lo(X)) / lo(X))
        E = envelope(D, K1)
        delta = PI**2 * (D - 2) ** 2 / (2 * iv.mpf(3) ** (2 * K1 + 2))
        margin = iv.mpf(4) / (D + 2) - delta
        assert lt(E, margin), D
        assert lt(2 * E, iv.mpf(8 * D) / (D + 2) ** 2), D
        L0 = int(mp.floor(lo(X)))
        L0 = L0 if L0 % 2 else L0 - 1
        T = (iv.mpf("1.9") * (D - 1) + iv.mpf(7) / (D - 2)) / iv.mpf(9) ** L0
        assert lt(T, iv.mpf("1e-70")) and lt(T, E / 1000)
        w = (-iv.log(1 - E) + T) / A
        assert lt(w * A, iv.mpf("1.1") * E)
        gap = iv.log((1 - 2 * E) * (iv.mpf(D + 2) / (D - 2)) ** 2)
        assert lt(T, gap)
        if window_min is None or lo(gap) < window_min[0]:
            window_min = (lo(gap), D)
        rows.append((D, K1, E, margin, lo(margin / E)))
    assert len(rows) == 182
    for a, b in pairwise(rows):
        assert lt(b[2], a[2])
    E38 = rows[0][2]
    assert lt(E38, iv.mpf("0.0897"))
    print(
        f"(*) E_(K_1)(D) < 4/(D+2) - delta(D) certified at all {len(rows)} even D = 38..400, K_1 asserted >= K* + 1 at each"
    )
    for D, K1, E, margin, slack in rows:
        if D in (38, 40, 42, 60, 100, 200, 400):
            print(
                f"  D = {D}: K_1 = {K1}, E <= {show(E, 6)}, slack factor >= {mp.nstr(slack, 6)}"
            )
    print(
        f"  widest relative enclosure of K*: {mp.nstr(widest, 3)}; E at D = 38 at most {show(E38, 5)}, E decreasing on the range"
    )
    print(
        f"  2E < 8D/(D+2)^2, T(D) < 1e-70, T(D) < E/1000 and w A < 1.1 E at every D; window margin log((1-2E)((D+2)/(D-2))^2) at least {mp.nstr(window_min[0], 4)} (D = {window_min[1]})"
    )
    D = 400
    assert lt(logR, iv.mpf("0.2275"))
    bound = D * ((D - 1) * iv.mpf("0.2275") + iv.mpf("0.69") / (D - 2)) / 4 + 3
    assert lt(bound, iv.mpf(D * D) / 16) and rows[-1][1] + 1 <= D * D // 16
    a, b, c, s = (
        iv.mpf(x) for x in ("0.7528157", "0.7052518", "0.2266816", "0.8900159")
    )
    F = (
        (D + 2)
        * ((iv.mpf(D * D) / 4) * (a ** (D - 1) + b ** (D - 1)) + 2 * c ** (D - 1))
        * iv.exp((iv.mpf(D * D) / 8) * s ** (D - 1))
    )
    G = (D + 2) * PI**2 * D**2 / (2 * iv.mpf(9) ** (D + 1))
    assert lt(F + G, iv.mpf("1.02e-42"))
    print(
        f"majorant for D >= 400 with K_1 = floor(K*) + 2: F(400) <= {show(F, 4)}, (D+2) delta <= {show(G, 3)}, F(400) + G(400) < 1.02e-42, so (D+2)(E + delta) < 2 beyond 400"
    )
    certs = {}
    for M in machines(3, range(2, 37, 2)):
        certs[M.D] = kmin(M, 200)
    assert certs[12] == 8 and certs[30] == 50
    print(
        "exact Collatz-Wielandt certificates, base 3, even D = 2..36, K_min:",
        " ".join(f"{D}:{K}" for D, K in certs.items()),
    )
    print(f"base3 done in {time.time() - start:.1f} s")


def g5(x):
    c = iv.cos(x)
    return 4 * c * c + 2 * c - 2


def verb_base5():
    start = time.time()
    hat = {i: g5(2 * PI / 5**i) for i in range(2, 41)}
    a = iv.mpf(5) / 4 * (2 * PI) ** 2
    tail = a / (24 * iv.mpf(25) ** 40)
    C = {}
    for j in (2, 3):
        prod = iv.mpf(1)
        for i in range(j + 1, 41):
            prod *= 4 / hat[i]
        C[j] = prod / (1 - iv.mpf([0, up(tail)]))
    assert lt(C[2], iv.mpf("1.0033008")) and lt(C[3], iv.mpf("1.0001317"))
    vals = {
        n: abs(g5(2 * PI * n / 25)) for n in range(25) if n % 5 and n not in (1, 24)
    }
    top = max(vals, key=lambda n: up(vals[n]))
    assert top in (2, 23) and all(
        lt(vals[n], vals[2]) for n in vals if n not in (2, 23)
    )
    assert lt(vals[2], iv.mpf("2.8242670")) and lt(iv.mpf("3.6897796"), hat[2])
    first = vals[2] / hat[2]
    assert lt(first, iv.mpf("0.7654303")) and lt(first * C[2], iv.mpf("0.7679580"))
    third = (1 + 12 * PI / 125) / (4 - 5 * (2 * PI / 125) ** 2) * C[3]
    fourth = (1 + 12 * PI / 625) / (4 - 5 * (2 * PI / 625) ** 2) * C[3]
    assert lt(third, iv.mpf("0.3264722")) and lt(fourth, iv.mpf("0.2651481"))
    print(
        f"Lemma 6.4: C_2 <= {show(C[2])}, C_3 <= {show(C[3])}; branch j = 2 <= {show(first * C[2])} (n = +-2, |g| <= {show(vals[2])}, hat g_2 >= {mp.nstr(lo(hat[2]), 8)}), j = 3 <= {show(third)}, j >= 4 <= {show(fourth)}; all below 0.768"
    )
    mp.dps = 30
    hatf = {
        i: 2 * mp.cos(2 * mp.pi / 5**i) + 2 * mp.cos(4 * mp.pi / 5**i)
        for i in range(2, 7)
    }
    for K in (2, 3, 4):
        N = 5 ** (K + 1)
        best, arg, branch = mp.mpf(0), None, {}
        for n in range(N):
            if n % 5 == 0 or n in (1, N - 1):
                continue
            ratio = mp.mpf(1)
            for i in range(2, K + 2):
                x = 2 * mp.pi * n / 5**i
                ratio *= abs(2 * mp.cos(x) + 2 * mp.cos(2 * x)) / hatf[i]
            j = next(i for i in range(2, K + 2) if n % 5**i not in (1, 5**i - 1))
            branch[min(j, 4)] = max(branch.get(min(j, 4), 0), ratio)
            if ratio > best:
                best, arg = ratio, n
        assert arg in (2, N - 2)
        assert (
            branch[2] < 0.7679580
            and branch.get(3, 0) < 0.3264722
            and branch.get(4, 0) < 0.2651481
        )
        print(
            f"  exact maximum of Q_K(n)/Q_K(1) at K = {K}: {mp.nstr(best, 6)} at n = +-2; branch maxima {', '.join(f'j{k}: {mp.nstr(v, 5)}' for k, v in sorted(branch.items()))}"
        )
    mp.dps = 80
    r = iv.mpf("0.768")

    def kappa(D, K):
        return (iv.mpf(D + 4) / (D + hat[2])) * (iv.mpf(D + 4) / (D + hat[3])) ** (
            K - 1
        )

    def lhs(D, K):
        return (2 * 5**K - 1) * r ** (D - 1) * kappa(D, K)

    def rhs(D, K):
        return iv.cos(PI * (D - 2) / 5 ** (K + 1))

    left, right = [], []
    for D in range(18, 33, 2):
        assert 5**3 >= 4 * (D - 2) and lt(lhs(D, 2), rhs(D, 2))
        left.append(mp.nstr(up(lhs(D, 2)), 4))
        right.append(mp.nstr(lo(rhs(D, 2)), 4))
    assert lt(rhs(16, 2), lhs(16, 2))
    print(
        f"Proposition 6.5 at K = 2, D = 18..32: left {', '.join(left)}; right {', '.join(right)}; at D = 16 it fails, {mp.nstr(up(lhs(16, 2)), 4)} against {mp.nstr(lo(rhs(16, 2)), 4)}"
    )

    def K_of(D):
        return next(K for K in range(1, 40) if 5 ** (K + 1) >= 4 * (D - 2))

    for D in range(34, 601, 2):
        K = K_of(D)
        assert K >= 2 and lt(lhs(D, K), rhs(D, K)), D
    ranges = 0
    for K in range(2, 10):
        first = max(34, 5**K // 4 + 3 + (5**K // 4 + 3) % 2)
        last = min(10**6, (5 ** (K + 1) // 4 + 2) // 2 * 2)
        assert K_of(first) == K == K_of(last) and K_of(first - 2) < K or first == 34
        assert lt(lhs(first, K), rhs(last, K)), K
        ranges += 1
    B = (1 + (4 - hat[2]) / (34 + hat[2])) * iv.exp(
        (4 - hat[3]) * (iv.log(136) / iv.log(5) - 1) / (34 + hat[3])
    )
    slope = (34 + hat[3]) / (34 * iv.log(5)), iv.log(136) / iv.log(5) - 1
    h34 = 34 * r**33
    goal = iv.cos(PI / 4) / (8 * B)
    assert (
        lt(B, iv.mpf("1.0089188"))
        and lt(h34, iv.mpf("0.0056028"))
        and lt(iv.mpf("0.0876069"), goal)
    )
    assert lt(slope[0], slope[1])
    print(
        f"Proposition 6.5 holds at every even D = 34..600 and on all {ranges} depth ranges K = 2..9 to D = 10^6 (left side at the range start below right side at its end)"
    )
    print(
        f"proof of Theorem 6.1 for D >= 34: B(34) <= {show(B)}, slope terms {show(slope[0], 4)} < {mp.nstr(lo(slope[1]), 5)}, h(34) <= {show(h34)}, cos(pi/4)/(8 B(34)) >= {mp.nstr(lo(goal), 7)}, 1/|ln r| = {mp.nstr(-1 / mp.log(mp.mpf('0.768')), 5)}"
    )
    dims = sorted(set(range(2, 401, 2)) | {314, 316})
    staircase, deaths, rowdepth = {}, {}, {}
    for M in machines(5, dims):
        D = M.D
        if D <= 64:
            staircase[D] = kmin(M, 4)
        rowdepth[D] = row_depth(M, 6)
        e = next(k for k in range(1, 9) if 5**k >= 2 * D - 3)
        assert rowdepth[D] == max(1, e - 1), D
        if 66 <= D <= 314:
            assert not depth(M, 3)[0], D
        if 16 <= D <= 320:
            bad, beta, nxt = depth(M, 2)
            assert (not bad) == (D <= 64), D
            if D == 66:
                c = bad[0]
                deficit = mp.mpf(M.fill * beta[c] - 5 * nxt[c]) / (M.fill * beta[c])
                deaths[2] = (D, bad, deficit)
        if D in (14, 16):
            bad, _, _ = depth(M, 1)
            assert (not bad) == (D == 14)
            if D == 16:
                deaths[1] = (D, bad, None)
        if D in (314, 316):
            bad, _, _ = depth(M, 3)
            assert (not bad) == (D == 314)
            if D == 316:
                deaths[3] = (D, bad, None)
                assert not depth(M, 4)[0]
        if D == 8:
            left8 = perron(M)
        if D == 40:
            row = M.step(M.v)
            assert row[19] > 0
    assert all(
        staircase[D] == (0 if D == 2 else 1 if D <= 14 else 2) for D in staircase
    )
    print(
        "exact certificates, base 5: K_min = 0 at D = 2, 1 at D = 4..14, 2 at D = 16..64 (every even D)"
    )
    print(
        f"  depth 1 fails first at D = 16, carries {deaths[1][1]}; depth 2 holds at every even D = 16..64 and fails at every even D = 66..320, first at D = 66, carries {deaths[2][1]}, relative deficit {mp.nstr(deaths[2][2], 3)}; depth 3 holds at every even D = 66..314, fails at 316 at carries {deaths[3][1]}, depth 4 holds at 316"
    )
    print(
        "  minimal row-certificate depth t(D) = max(1, ceil(log_5(2D-3)) - 1) at every even D = 2..400; (v^T M)_19 > 0 at D = 40"
    )
    print(
        f"  D = 8: right Perron mass p = {mp.nstr(left8[0], 6)} > 1/5, left-vector mass {mp.nstr(left8[1], 4)} < 1/5, 5 rho - fill + (D-1)(5p - 1) = {mp.nstr(left8[2], 3)}"
    )
    print(f"base5 done in {time.time() - start:.1f} s")


def perron(M):
    mp.dps = 40
    S = range(-M.r, M.r + 1)
    A = mp.matrix([[M.entry(new, old) for old in S] for new in S])
    E, ER = mp.eig(A)
    k = max(range(len(E)), key=lambda i: mp.re(E[i]))
    rho = mp.re(E[k])
    w = [mp.re(ER[i, k]) for i in range(len(S))]
    E2, EL = mp.eig(A.T)
    k2 = max(range(len(E2)), key=lambda i: mp.re(E2[i]))
    y = [mp.re(EL[i, k2]) for i in range(len(S))]
    zero = [i for i, c in enumerate(S) if c % M.q == 0]
    p = sum(w[i] for i in zero) / sum(w)
    pl = sum(y[i] for i in zero) / sum(y)
    resid = M.q * rho - M.fill - (M.D - 1) * (M.q * p - 1) * (-1) ** (M.D - 1)
    mp.dps = 80
    return p, pl, resid


def verb_bases():
    start = time.time()
    report = {}
    for q, dims in (
        (7, list(range(2, 27, 2)) + [172, 174]),
        (9, range(2, 57, 2)),
        (11, range(2, 75, 2)),
    ):
        report[q] = {M.D: kmin(M, 4) for M in machines(q, dims)}
    seven = report[7]
    assert all(
        seven[D] == (0 if D == 2 else 1 if D < 26 else 2) for D in range(2, 27, 2)
    )
    assert seven[172] == 2 and seven[174] == 3
    assert max(report[9].values()) <= 2 and max(report[11].values()) <= 2
    for q in (7, 9, 11):
        steps = [
            D
            for D in sorted(report[q])
            if D - 2 in report[q] and report[q][D] != report[q][D - 2]
        ]
        print(
            f"base {q}: exact certificates at every even D in range, K_min <= {max(report[q].values())}, depth steps at D = {steps}"
            + (", K_min(172) = 2, K_min(174) = 3" if q == 7 else "")
        )
    for q, dims in ((7, range(2, 41, 2)), (9, range(2, 57, 2)), (11, range(2, 75, 2))):
        for M in machines(q, dims):
            assert min(tower(M, 25)) > 0
    for q in (5, 7, 9, 11, 13):
        for D in (8, 12):
            assert min(tower(Machine(q, D), 25)) > 0
    print(
        "V(L) > 0 at every L <= 25: base 7 even D = 2..40, base 9 even D <= 56, base 11 even D <= 74, bases 5..13 at D = 8, 12"
    )
    boxes = {
        3: range(2, 16),
        5: range(2, 27),
        7: range(2, 19),
        9: range(2, 13),
        11: range(2, 13),
    }
    cases = 0
    for q, dims in boxes.items():
        for M in machines(q, dims):
            assert (above(M) > 0) == (M.D % 2 == 1), (q, M.D)
            cases += 1
    for q, D in ((21, 3), (31, 5), (51, 5), (101, 3)):
        M = Machine(q, D)
        assert above(M) > 0
        cases += 1
    print(
        f"sign law sgn(rho - fill/q) = (-1)^(D+1) by exact Sturm counts in {cases} cases: base 3 D = 2..15, base 5 D = 2..26, base 7 D = 2..18, bases 9, 11 D = 2..12, (q, D) = (21, 3), (31, 5), (51, 5), (101, 3)"
    )
    c = (mp.sqrt(7) - 1) / 6
    g7 = 2 * (4 * c**3 + 2 * c**2 - 2 * c - 1)
    assert abs(g7 + (34 + 14 * mp.sqrt(7)) / 27) < mp.mpf("1e-60")
    grid = min(
        2 * (mp.cos(t) + mp.cos(2 * t) + mp.cos(3 * t))
        for t in mp.linspace(0, mp.pi, 4001)
    )
    assert grid >= g7
    print(f"min g_7 = {mp.nstr(g7, 10)} = -(34 + 14 sqrt 7)/27; min g_5 = -9/4")
    print(f"bases done in {time.time() - start:.1f} s")


def verb_odd():
    start = time.time()
    p = (1 << 61) - 1
    for M in machines(3, range(2, 41)):
        delta = bareiss(pencil(M))
        assert delta != 0 and delta % p == detmod(pencil(M), p)
        assert (delta > 0) == (M.D % 2 == 0)
    print(
        "sgn Delta_D = (-1)^D at base 3, D = 2..40, Bareiss with a mod-(2^61 - 1) cross-check"
    )
    for q, dims in ((3, range(7, 81, 6)), (5, range(11, 81, 10))):
        rows = []
        for M in machines(q, dims):
            delta = bareiss(pencil(M))
            assert delta != 0 and delta % p == detmod(pencil(M), p) != 0
            rows.append(M.D)
        print(
            f"Delta_D != 0 at base {q}, odd D = 1 mod {q}, D = {rows[0]}..{rows[-1]} ({len(rows)} dims)"
        )
    rows5 = {}
    for M in machines(5, range(2, 157, 2)):
        e = v2det(M.even)
        assert e < 2 * (M.D - 1) + v2(M.D + 4)
        rows5[M.D] = e
        if M.D in (120, 130):
            assert v2(abs(bareiss(M.even))) == e
    rows3 = []
    for M in machines(3, range(4, 119, 6)):
        e = v2det(M.even)
        assert e < M.D - 1 + v2(M.D + 2)
        rows3.append(e)
    for M in machines(3, range(2, 31)):
        d = bareiss(M.even)
        assert v2det(M.even) == v2(abs(d))
    peak = max(rows5.values())
    high = {D: e for D, e in rows5.items() if e > 26}
    print(
        f"v_2(det M_even) at base 5, even D = 2..156: max {peak}, above 26 at {high} (Bareiss confirms D = 120, 130), D = 4..12 reading {[rows5[D] for D in range(4, 13, 2)]}, below 2(D-1) + v_2(D+4) everywhere"
    )
    print(
        f"v_2(det M_even) at base 3, D = 4 mod 6, D = 4..118: max {max(rows3)}, first five {rows3[:5]}, below D - 1 + v_2(D+2) everywhere; 2-adic elimination matches Bareiss at D = 2..30"
    )
    print(f"odd done in {time.time() - start:.1f} s")


def verb_anchors():
    start = time.time()
    M = Machine(3, 3)
    assert M.even == [[6, 6], [1, 3]] and charpoly(M.even) == [1, -9, 12]
    rho = (9 + mp.sqrt(33)) / 2
    print(
        f"D = 3: M_even = {M.even}, chi = x^2 - 9x + 12, rho = {mp.nstr(rho, 7)}, log_3 rho = {mp.nstr(mp.log(rho, 3), 7)}, log_3 20 - 1 = {mp.nstr(mp.log(20, 3) - 1, 7)}"
    )
    four = Machine(3, 4)
    chi4 = charpoly(four.even)
    assert chi4 == [1, -11, -66] and 16**2 - 11 * 16 - 66 == 14 and four.fill == 48
    print(
        f"D = 4: chi = x^2 - 11x - 66, rho = {mp.nstr((11 + mp.sqrt(385)) / 2, 8)} against fill/3 = 16, chi(16) = 14"
    )
    fails = 0
    for M in machines(3, range(3, 51)):
        B, D, n = M.P, M.D, M.r + 1

        def b(i):
            return B[i] if 0 <= i < len(B) else 0

        rows = [
            b(D - 3 * i) + sum(b(D + j - 3 * i) + b(D - j - 3 * i) for j in range(1, n))
            for i in range(n)
        ]
        sigma = [sum(M.entry(new, c) for new in range(-M.r, M.r + 1)) for c in range(n)]
        formula = [
            (M.fill + 2 * (-1) ** (D - 1) * (D - 1)) // 3
            if c % 3 == 0
            else (M.fill - (-1) ** (D - 1) * (D - 1)) // 3
            for c in range(n)
        ]
        assert sigma == formula
        fails += rows != formula
        if D == 3:
            assert rows == [12, 4] and formula == [8, 6]
    assert fails == 48
    print(
        f"carry-orientation sums equal fill/3 + (2/3)(-1)^(D-1)(D-1)cos(2 pi c/3) at D = 3..50; the transposed even-basis row sums differ at all {fails} dims, D = 3 reading [12, 4] against [8, 6]"
    )
    match = 0
    for M in machines(3, range(2, 51)):
        D = M.D
        top = M.fill + 2 * (D - 1)
        assert top < 3 * 3 ** (D - 1) or M.even == [[2]]
        match += (-1) ** (D + 1) == -1
    print(
        f"rho_D < 3^(D-1) at every D = 2..50 by the pinning rho_D <= (fill + 2(D-1))/3 (rho_2 = 2 < 3), so sgn(rho_D - 3^(D-1)) matches (-1)^(D+1) in {match} of 49 dims; rho_3 - 9 = {mp.nstr(rho - 9, 16)}"
    )
    M5 = Machine(5, 3)
    c5 = charpoly(M5.even)
    r5 = max(mp.polyroots(c5, maxsteps=200, extraprec=200), key=mp.re)
    print(f"base 5, D = 3: rho_3 - 25 = {mp.nstr(mp.re(r5) - 25, 11)}")
    four = tower(Machine(3, 4), 40)
    assert min(four) > 0 and four[1] == 4
    six = tower(Machine(3, 6), 3)
    assert six[3] == 135092 and six[1] < 0
    print(
        f"V(L) > 0 at D = 4 for L <= 40 with V(1) = {four[1]}; D = 6: V(1) = {six[1]}, V(3) = {six[3]}"
    )
    breaks = {}
    for M in machines(3, range(6, 41, 2)):
        V = tower(M, 40)
        assert V[1] < 0
        W = [-(1 if V[k - 1] > 0 else -1) for k in range(1, 41)]
        breaks[M.D] = next((k + 1 for k in range(1, 40) if W[k] == W[k - 1]), None)
    print(
        "first break of the W_k sign alternation, even D = 6..40, k <= 40:",
        " ".join(f"{D}:{k}" for D, k in breaks.items()),
    )
    mp.dps = 60
    worst = mp.mpf(0)
    for D in range(2, 13):
        M = Machine(3, D)
        beta = [1] * (M.r + 1)
        for L in range(5):
            prod = mp.mpf(0)
            for m in range(3**L):
                term = mp.mpf(1)
                for j in range(L):
                    psi = 2 * mp.pi * m * 3**j / 3**L
                    term *= (2 * mp.cos(psi)) ** (D - 1) * (D + 2 * mp.cos(psi))
                prod += term
            worst = max(worst, abs(prod / 3**L - beta[0]) / beta[0])
            beta = M.step(beta)
    assert worst < mp.mpf("1e-50")
    mp.dps = 80
    print(
        f"product formula of Theorem 5.1 against the carry census b_D(L) at D = 2..12, L <= 4: worst relative gap {mp.nstr(worst, 3)}"
    )
    print(f"anchors done in {time.time() - start:.1f} s")


def main():
    verbs = {
        "base3": verb_base3,
        "base5": verb_base5,
        "bases": verb_bases,
        "odd": verb_odd,
        "anchors": verb_anchors,
    }
    chosen = sys.argv[1:] or list(verbs)
    for name in chosen:
        verbs[name]()


if __name__ == "__main__":
    main()
