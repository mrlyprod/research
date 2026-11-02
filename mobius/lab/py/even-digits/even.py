import math
import random
import sys
import time
from fractions import Fraction

import numpy as np

GAMMA = 0.5772156649015329
C_P = math.sqrt(2) - 4 / math.pi

# CLOSED FORMS


def lam(N):
    x0 = (N / math.pi) * (
        math.log(N + 3)
        + GAMMA
        + math.log(math.tan(3 * math.pi / 8 + math.pi / (4 * N)))
    ) + C_P * (N + 1) ** 2 / (8 * N)
    return (N + 1) / 2 + x0 + 0.5 / math.sin(math.pi / (2 * N))


def half(N):
    return list(range((N + 1) // 2))


def even(N):
    return list(range(0, N, 2))


def classes(N):
    return f"{N % 4} mod 4, {N % 8} mod 8, {N % 3} mod 3"


# TRANSFER


def factor_direct(digits, th):
    z = np.zeros(th.shape, dtype=complex)
    for v in digits:
        z += np.exp(2j * np.pi * v * th)
    return np.abs(z)


def factor_ap(n, c, th):
    s = np.sin(np.pi * c * th)
    num = np.sin(np.pi * n * c * th)
    out = np.empty(th.shape)
    small = np.abs(s) < 1e-12
    out[~small] = np.abs(num[~small] / s[~small])
    out[small] = n
    return out


def grid_sum(N, digits, step, level, s, direct):
    y = N**level
    a = np.arange(y, dtype=np.int64)
    out = np.ones(y)
    for j in range(level):
        shift = float((Fraction(N) ** j * s) % 1)
        th = (((pow(N, j, y) * a) % y) / y + shift) % 1.0
        out *= factor_direct(digits, th) if direct else factor_ap(len(digits), step, th)
    return float(out.sum())


SHIFTS = [
    Fraction(0),
    Fraction(1, 2),
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(1, 6),
    Fraction(3819660113, 10**10),
]


def transfer():
    t0 = time.time()
    print(
        "THE EVEN DIGITS AT THE SHIFTED GRID: E = {0, 2, .., base-1} = 2F, F = {0..fill-1}, odd base = 2 fill - 1"
    )
    print(
        "Sigma^E_i(s) = sum over a < base^i of |hat E_i(s + a/base^i)|, read from E's digits at bases 3..41 and through |sin(2 pi fill t)/sin(2 pi t)| above; Sigma^F_i(2s) from F's"
    )
    print(
        "each row: base, its classes, the levels read, max |Sigma^E_i(s)/Sigma^F_i(2s) - 1|, max Sigma^E_i(s)/((3/2) lambda^i), min Sigma^E_i(s)/base^i, Sigma^F_i(1/2)/Sigma^F_i(0) at the top level"
    )
    bases = list(range(3, 42, 2)) + [99, 101, 103, 105, 1001, 1003, 94939, 94941]
    worst_id, worst_up, worst_lo = 0.0, 0.0, math.inf
    seen = set()
    for N in bases:
        fill = (N + 1) // 2
        direct = N <= 41
        cap = 3 * 10**7 if direct else 2 * 10**6
        levels = [
            i for i in range(1, 40) if N**i * (fill if direct else 4) * i <= cap
        ] or [1]
        L = lam(N)
        rid, rup, rlo = 0.0, 0.0, math.inf
        top = {}
        for i in levels:
            for s in SHIFTS:
                se = grid_sum(N, even(N), 2, i, s, direct)
                sf = grid_sum(N, half(N), 1, i, (2 * s) % 1, direct)
                rid = max(rid, abs(se / sf - 1))
                rup = max(rup, se / (1.5 * L**i))
                rlo = min(rlo, se / N**i)
                if i == levels[-1] and s in (Fraction(0), Fraction(1, 4)):
                    top[s] = sf
        assert rid < 1e-9 and rup < 1 and rlo >= 1 - 1e-9
        worst_id, worst_up, worst_lo = (
            max(worst_id, rid),
            max(worst_up, rup),
            min(worst_lo, rlo),
        )
        seen.add((N % 4, N % 3))
        print(
            f"  {N:6d}  {classes(N)}  levels 1..{levels[-1]}  {rid:.1e}  {rup + 5e-7:.6f}  {rlo - 5e-7:.6f}  {top[Fraction(1, 4)] / top[Fraction(0)]:.6f}"
        )
    print(f"classes (base mod 4, base mod 3) sampled: {sorted(seen)}")
    print(
        f"identity Sigma^E_i(s) = Sigma^F_i(2s) to {worst_id:.1e}; Sigma^E_i(s) <= {worst_up + 5e-7:.6f} (3/2) lambda^i; Sigma^E_i(s) >= {worst_lo - 5e-7:.6f} base^i"
    )
    print(f"transfer in {time.time() - t0:.1f}s")


# REGIONS


def mu_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    mu[0] = 0
    sieve = np.ones(X + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(X**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    primes = np.nonzero(sieve)[0]
    for p in primes:
        mu[p::p] = -mu[p::p]
    for p in primes:
        if p * p > X:
            break
        mu[p * p :: p * p] = 0
    return mu, sieve


def strings(N, digits, level):
    out = np.zeros(1, dtype=np.int64)
    d = np.array(digits, dtype=np.int64)
    for i in range(level):
        out = (out[:, None] + d[None, :] * (N**i)).ravel()
    out.sort()
    return out


def convergent(num, den, cap):
    n = num.astype(np.int64).copy()
    d = np.full_like(n, den)
    h1, h2 = np.ones_like(n), np.zeros_like(n)
    k1, k2 = np.zeros_like(n), np.ones_like(n)
    c = n // d
    h1, h2 = c * h1 + h2, h1
    k1, k2 = c * k1 + k2, k1
    n, d = d, n - c * d
    live = d > 0
    while live.any():
        c = np.zeros_like(n)
        c[live] = n[live] // d[live]
        hn, kn = c * h1 + h2, c * k1 + k2
        go = live & (kn <= cap)
        h2, k2 = np.where(go, h1, h2), np.where(go, k1, k2)
        h1, k1 = np.where(go, hn, h1), np.where(go, kn, k1)
        n, d = np.where(go, d, n), np.where(go, n - c * d, 0)
        live = go & (d > 0)
    return h1, k1


def smooth(d, N):
    while d > 1:
        g = math.gcd(d, N)
        if g == 1:
            return False
        while d % g == 0:
            d //= g
    return d == 1


def outside(d, N):
    while True:
        g = math.gcd(d, N)
        if g == 1:
            return d
        while d % g == 0:
            d //= g


def regions_one(N, k, Z, mu):
    fill = (N + 1) // 2
    y = N**k
    DE = np.zeros(y)
    DE[strings(N, even(N), k)] = 1.0
    DF = np.zeros(y)
    DF[strings(N, half(N), k)] = 1.0
    m = mu[:y].astype(np.float64)
    hatE = np.conj(np.fft.fft(DE))
    hatF = np.conj(np.fft.fft(DF))
    Sneg = np.fft.fft(m)
    exact = float(m[DE > 0].sum())
    a = np.arange(y, dtype=np.int64)
    Q = int(y**0.6)
    l, d = convergent(a, y, Q)
    h = np.abs(a * d - l * y)
    assert (d <= Q).all() and (h * Q <= y).all() and (np.gcd(l, d) == 1).all()
    sm = np.array([v >= 1 and smooth(v, N) for v in range(Q + 1)])
    sm2 = sm | np.array([v % 4 == 2 and smooth(v // 2, N) for v in range(Q + 1)])
    A = d >= y**0.4
    C = (~A) & (d < Z) & (h < Z)
    B = (~A) & (~C)
    C2 = C & sm2[d]
    C1 = C & (~sm2[d])
    T = C & sm2[d] & (~sm[d])
    term = hatE * Sneg / y
    parts = {
        nm: term[msk].sum().real
        for nm, msk in (("A", A), ("B", B), ("C1'", C1), ("C2'", C2))
    }
    tot = sum(parts.values())
    assert abs(tot - exact) < 1e-6 * y
    nZ = sum(1 for v in range(1, Z) if sm2[v])
    assert C1.sum() <= 3 * Z * Z and C2.sum() <= 3 * Z * nZ
    for ai in np.nonzero(C1)[0]:
        dd, ll = int(d[ai]), int(l[ai])
        g = math.gcd(2 * ll, dd)
        assert outside(dd // g, N) > 1
    eF = np.abs(hatF) / fill**k
    eE = np.abs(hatE) / fill**k
    print(
        f"  base {N} ({classes(N)}), level {k}, y = {y}, Z = {Z}: exact {int(exact)}, the four regions sum to {tot:.6f}; points A {A.sum()}, B {B.sum()}, C1' {C1.sum()} <= 3 Z^2, C2' {C2.sum()} <= 3 Z N'_Z = {3 * Z * nZ}"
    )
    print(
        f"    the {T.sum()} points with d = 2e: l^1 mass {eE[T].sum():.6f} fill^k under E against {eF[T].sum() + 5e-7:.6f} fill^k under F, largest |hat E_k|/fill^k {eE[T].max():.6f} against {eF[T].max() + 5e-7:.6f}; on C1' largest {eE[C1].max() + 5e-7:.6f}"
    )
    pairs = sorted(set(zip(l[T].tolist(), d[T].tolist())))
    n = np.arange(y)
    worst = 0.0
    for ll, dd in pairs:
        w = m * np.exp(2j * np.pi * ((n * ll) % dd) / dd)
        pm = np.abs(np.cumsum(w)).max()
        sel = T & (l == ll) & (d == dd)
        lhs = np.abs(Sneg[sel])
        rhs = (1 + 2 * np.pi * h[sel] / dd) * pm
        worst = max(worst, float((lhs / rhs).max()))
    assert worst <= 1 + 1e-9
    print(
        f"    partial summation at the {len(pairs)} fractions with d = 2e: |S(a/y)| <= (1 + 2 pi h/d) max_u |sum_(n <= u) mu(n) e(n l/d)|, largest ratio {worst + 5e-7:.6f}"
    )


def progression_identity(N, X, mu):
    checked = 0
    for e in (1, N, N * N):
        inv2 = pow(2, -1, e) if e > 1 else 0

        def M(v, q, r):
            return int(mu[: v + 1][r % q :: q].sum()) if v >= 0 else 0

        def O(v, r):
            s, j = 0, 0
            while v >> j:
                s += M(v >> j, e, (r * pow(inv2, j, e)) % e if e > 1 else 0)
                j += 1
            return s

        for u in (X, X // 3 + 7, X // 10 + 1):
            for r in range(2 * e):
                lhs = M(u, 2 * e, r)
                rhs = O(u, r) if r % 2 else -O(u // 2, (r * inv2) % e if e > 1 else 0)
                assert lhs == rhs
                checked += 1
    return checked


def lemma_draws(N, k, Z, draws, seed):
    fill = (N + 1) // 2
    E = even(N)
    y = N**k
    rng = random.Random(seed)
    cp = 1 - (2 / fill) * (1 - math.cos(math.pi / (4 * N)))

    def logratio(a):
        s = 0.0
        for i in range(k):
            ph = (pow(N, i, y) * a % y) / y
            z = sum(
                complex(math.cos(2 * math.pi * v * ph), math.sin(2 * math.pi * v * ph))
                for v in E
            )
            s += math.log(abs(z) / fill)
        return s

    def mloc(dd):
        return max(1, math.floor(math.log(dd / 2) / math.log(N) + 1e-12) + 1)

    worst, n, kinds = -1e9, 0, {"odd": 0, "2 mod 4": 0, "0 mod 4": 0}
    while n < draws:
        dd = rng.randrange(2, Z)
        if smooth(dd, N) or (dd % 4 == 2 and smooth(dd // 2, N)):
            continue
        ll = rng.randrange(dd)
        if math.gcd(ll, dd) != 1:
            continue
        a = (ll * y) // dd + rng.randrange(-Z, Z)
        hh = abs(a * dd - ll * y)
        if hh >= Z:
            continue
        assert 2 * hh / (dd * y) < y ** (-2 / 3) / (4 * N * (N - 1))
        d2 = dd // math.gcd(2, dd)
        worst = max(
            worst, logratio(a) - math.floor(2 * k / (3 * mloc(d2))) * math.log(cp)
        )
        kinds["odd" if dd % 2 else ("2 mod 4" if dd % 4 == 2 else "0 mod 4")] += 1
        n += 1
    assert worst <= 1e-9
    print(
        f"  base {N} ({classes(N)}), level {k}, d < {Z}, h < {Z}: {draws} seeded draws off d = 2e, kinds {kinds}; floor(2k/(3 m_d')) log c' minus log |hat E_k(a/y)|/fill^k is at least {math.floor(-worst * 10**4) / 10**4:.4f}"
    )


def witness_at(N, k):
    fill = (N + 1) // 2
    y = N**k
    a = (y - 1) // 2
    s = 0.0
    for i in range(k):
        ph = (pow(N, i, y) * a % y) / y
        z = sum(
            complex(math.cos(2 * math.pi * v * ph), math.sin(2 * math.pi * v * ph))
            for v in even(N)
        )
        s += math.log(abs(z) / fill)
    return math.exp(s)


def witness_limit(N):
    fill = (N + 1) // 2
    p, m = 1.0, 1
    while True:
        u = N ** (-m)
        f = abs(math.sin(math.pi * fill * u) / math.sin(math.pi * u)) / fill
        p *= f
        if 1 - f < 1e-16:
            break
        m += 1
    return p


def regions():
    t0 = time.time()
    print(
        "THE DISSECTION ON THE EVEN DIGITS at small base and level, every grid point a mod y"
    )
    print(
        "A: d >= y^(2/5); B: d < y^(2/5), max(d, h) >= Z; C2': d < Z, h < Z, d = e or 2e with e dividing a power of the base; C1': the rest of d < Z, h < Z"
    )
    mu, _ = mu_upto(2 * 10**6)
    for N, k, Z in ((3, 12, 16), (5, 9, 12), (7, 7, 16), (9, 6, 16), (11, 6, 16)):
        regions_one(N, k, Z, mu)
    print(
        "THE ARITHMETIC AT 2 base: M(u; 2e, r) = sum_j M(u/2^j; e, r 2^(-j)) at odd r and minus that sum at u/2, r/2 at even r, exact integers"
    )
    for N in (3, 5, 7, 9, 11):
        print(
            f"  base {N} ({classes(N)}): e in 1, base, base^2, every r mod 2e, three u up to 2 * 10^6: {progression_identity(N, 2 * 10**6, mu)} identities hold"
        )
    print("LEMMA A' THROUGH THE HALVING, and its witness at d = 2")
    for N, k, Z, seed in (
        (5, 45, 40, 1009),
        (7, 36, 40, 1013),
        (9, 30, 40, 1019),
        (11, 28, 40, 1021),
    ):
        lemma_draws(N, k, Z, 1500, seed)
    print(
        "  the witness a = (y-1)/2, d = 2, h = 1, at every level: |hat E_k(a/y)|/fill^k = prod_(m <= k) |hat F(base^(-m))|/fill, falling in k to a positive limit;"
    )
    print(
        "  k_0 is the least level where the pair's bound c'^floor(2k/3) at d = 2 falls under that limit, the ratio read from E's digits at k_0 where k_0 <= 400"
    )
    for N in (3, 5, 7, 9, 11, 101, 94939, 94941):
        fill = (N + 1) // 2
        lim = witness_limit(N)
        cp = 1 - (2 / fill) * (1 - math.cos(math.pi / (4 * N)))
        k0 = math.ceil(1.5 * (math.floor(math.log(lim) / math.log(cp)) + 1))
        tail = ""
        if k0 <= 400:
            r = witness_at(N, k0)
            assert r >= lim * (1 - 1e-9) and r > cp ** math.floor(2 * k0 / 3)
            tail = f", ratio at k_0 {r - 5e-7:.6f} against the bound {cp ** math.floor(2 * k0 / 3) + 5e-7:.6f}"
        print(
            f"    base {N} ({classes(N)}): limit {lim - 5e-7:.6f}, k_0 {k0:.4g}{tail}"
        )
    print(f"regions in {time.time() - t0:.1f}s")


# METER


def members(X, N, ok):
    n = np.arange(X + 1, dtype=np.int64)
    inside = np.ones(X + 1, dtype=bool)
    while n.any():
        inside &= ok[n % N] | (n == 0)
        n //= N
    inside[0] = False
    return inside


def meter():
    t0 = time.time()
    X = 10**7
    mu, isp = mu_upto(X)
    primes = np.nonzero(isp)[0]
    sign = np.where(np.arange(X + 1) % 2 == 0, 1, -1)
    print(
        "THE METER ON THE EVEN DIGITS below 10^7: 2 M_E(2X) = M^tw_F(X) - M_F(X) and M_E(2X + 1) = M_E(2X) at every X, exact integers;"
    )
    print(
        "then readings, evidence for nothing: max |M_E|/A_E and max |M_E|/A_E^(1/2) over A_E >= 100, max |M^tw_F|/A_F^(1/2); and the prime count psi_E(10^7) against kappa_E A_E(10^7)"
    )
    for N in (5, 7, 101, 103, 1009, 10007):
        okF = np.arange(N) <= (N - 1) // 2
        okE = np.arange(N) % 2 == 0
        inF = members(X, N, okF)
        inE = members(X, N, okE)
        assert (inE[0::2][: X // 2 + 1] == inF[: X // 2 + 1]).all() and not inE[
            1::2
        ].any()
        ME = np.cumsum(np.where(inE, mu, 0).astype(np.int64))
        MF = np.cumsum(np.where(inF, mu, 0).astype(np.int64))
        MT = np.cumsum(np.where(inF, mu * sign, 0).astype(np.int64))
        h = X // 2
        assert (2 * ME[0 : 2 * h + 1 : 2] == MT[: h + 1] - MF[: h + 1]).all()
        assert (ME[1 : 2 * h : 2] == ME[0 : 2 * h - 1 : 2]).all()
        AE = np.cumsum(inE)
        AF = np.cumsum(inF)
        se, sf = AE >= 100, AF >= 100
        psi = 0.0
        for p in primes:
            q = int(p)
            while q <= X:
                if inE[q]:
                    psi += math.log(p)
                q *= p
        cop = sum(1 for f in range(N) if okE[f] and math.gcd(f, N) == 1)
        phi = sum(1 for f in range(1, N) if math.gcd(f, N) == 1)
        kap = Fraction(N * cop, phi * ((N + 1) // 2))
        assert kap == Fraction(N, N + 1)
        print(
            f"  base {N} ({classes(N)}): A_E {AE[-1]}, M_E {ME[-1]}, max |M_E|/A_E {np.max(np.abs(ME[se]) / AE[se]):.4f}, max |M_E|/A_E^(1/2) {np.max(np.abs(ME[se]) / np.sqrt(AE[se])):.3f}, max |M^tw_F|/A_F^(1/2) {np.max(np.abs(MT[sf]) / np.sqrt(AF[sf])):.3f}; psi_E {psi:.4f} = {psi / math.log(2):.0f} log 2 against kappa_E A_E = {float(kap) * AE[-1]:.1f}, kappa_E = {kap}"
        )
    print(f"meter in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "transfer"
    {"transfer": transfer, "regions": regions, "meter": meter}[verb]()
