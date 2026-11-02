import bisect
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, getcontext
from fractions import Fraction
from itertools import pairwise
from math import comb, e, gcd, isqrt, log, pi

getcontext().prec = 80

# THE OBJECT

Q = 3
D = 2
CELLS = ((0, 0), (2, 0), (0, 2))
WEIGHT = (Fraction(3, 8), Fraction(3, 8), Fraction(1, 4))
LEVEL = 10
BINS = 24
FOLD_LO, FOLD_HI = 4, 8
MASS_A = (20, 40, 60, 80, 100)
MASS_MAX = 130.0
SPECTRUM_LEVELS = (6, 8, 10)
GRIPENBERG_LEN = 8
NORM_LEN = 14

PROBS = (
    ("1/3,1/3,1/3", (Fraction(1, 3), Fraction(1, 3), Fraction(1, 3))),
    ("2/5,2/5,1/5", (Fraction(2, 5), Fraction(2, 5), Fraction(1, 5))),
    ("3/8,3/8,1/4", (Fraction(3, 8), Fraction(3, 8), Fraction(1, 4))),
    ("1/2,1/4,1/4", (Fraction(1, 2), Fraction(1, 4), Fraction(1, 4))),
)

RENEWAL = (
    ("1/2,1/2,1/2", (Fraction(1, 2), Fraction(1, 2), Fraction(1, 2))),
    ("1/2,1/2,1/4", (Fraction(1, 2), Fraction(1, 2), Fraction(1, 4))),
    ("1/2,1/2,1/3", (Fraction(1, 2), Fraction(1, 2), Fraction(1, 3))),
    ("1/3,1/3,1/3", (Fraction(1, 3), Fraction(1, 3), Fraction(1, 3))),
    ("1/2,1/4,1/4", (Fraction(1, 2), Fraction(1, 4), Fraction(1, 4))),
    ("2/5,2/5,1/5", (Fraction(2, 5), Fraction(2, 5), Fraction(1, 5))),
    ("3/8,3/8,1/4", (Fraction(3, 8), Fraction(3, 8), Fraction(1, 4))),
)

CONTROLS = ("1/3,1/3,1/3", "1/2,1/4,1/4", "3/8,3/8,1/4")

NAPKIN_LENGTH = {
    "1/3,1/3,1/3": (
        "0.652441",
        "0.022937",
        "0.000000",
        "0.0208  0.0094  0.0047  0.0021",
    ),
    "2/5,2/5,1/5": (
        "0.536956",
        "0.019130",
        "0.066973",
        "0.0224  0.0091  0.0043  0.0022",
    ),
    "3/8,3/8,1/4": (
        "0.580858",
        "0.020478",
        "0.044604",
        "0.0205  0.0092  0.0041  0.0022",
    ),
    "1/2,1/4,1/4": (
        "0.427372",
        "0.014471",
        "0.153880",
        "0.0299  0.0153  0.0077  0.0039",
    ),
}

NAPKIN_MASS = {
    "1/2,1/2,1/2": ("1.5849625007", "1.09812  1.09852  1.09812  1.09812  1.09803"),
    "1/2,1/2,1/4": ("1.2715533032", "0.88098  0.88130  0.88098  0.88098  0.88091"),
    "1/2,1/2,1/3": ("1.3646005647", "0.20259  0.14175  0.11707  0.10183  0.09030"),
    "1/3,1/3,1/3": ("1.0000000000", "1.09825  1.09825  1.09825  1.09825  1.09825"),
    "1/2,1/4,1/4": ("1.0000000000", "0.69284  0.69309  0.69284  0.69284  0.69278"),
    "2/5,2/5,1/5": ("1.0000000000", "0.23297  0.19792  0.18623  0.17623  0.16776"),
    "3/8,3/8,1/4": ("1.0000000000", "0.23133  0.16019  0.13671  0.11972  0.10734"),
}

NAPKIN_RUNG0 = {
    "1/3,1/3,1/3": ("[3, 3, 3]", "3", "1.000000000", "1.000000000", "-1.000000000"),
    "2/5,2/5,1/5": (
        "[18/5, 18/5, 9/5]",
        "18/5",
        "0.834043767",
        "1.464973521",
        "-1.165956233",
    ),
    "3/8,3/8,1/4": (
        "[27/8, 27/8, 9/4]",
        "27/8",
        "0.892789261",
        "1.261859507",
        "-1.107210739",
    ),
    "1/2,1/4,1/4": (
        "[9/2, 9/4, 9/4]",
        "9/2",
        "0.630929754",
        "1.261859507",
        "-1.369070246",
    ),
}

NAPKIN_RUNG1 = {
    "hat (1/2, 1, 1/2)": ("0.5000000", "0.5000000", "1.0000000", "1.0000000"),
    "D4": ("0.6830127", "0.7105812", "0.4929285", "0.5500157"),
}

FAILURES = []


def safe_log2(x, digits, mode):
    v = -(Decimal(x.numerator) / Decimal(x.denominator)).ln() / Decimal(2).ln()
    return str(v.quantize(Decimal(1).scaleb(-digits), rounding=mode))


def check(tag, got, want):
    if got != want:
        FAILURES.append("%s reads %s wants %s" % (tag, got, want))
    return got


# THE ARITHMETIC CLASS, EXACT


def factor(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def exponent_rows(w):
    primes = set()
    for x in w:
        primes |= set(factor(x.numerator)) | set(factor(x.denominator))
    primes = sorted(primes)
    rows = []
    for x in w:
        up, dn = factor(x.numerator), factor(x.denominator)
        rows.append(tuple(up.get(p, 0) - dn.get(p, 0) for p in primes))
    return primes, rows


def rank(rows):
    m = [[Fraction(v) for v in r] for r in rows]
    top = 0
    for col in range(len(m[0])):
        p = next((k for k in range(top, len(m)) if m[k][col] != 0), None)
        if p is None:
            continue
        m[top], m[p] = m[p], m[top]
        d = m[top][col]
        m[top] = [v / d for v in m[top]]
        for k in range(len(m)):
            if k != top and m[k][col] != 0:
                f = m[k][col]
                m[k] = [a - f * b for a, b in zip(m[k], m[top])]
        top += 1
    return top


def arithmetic_class(w):
    primes, rows = exponent_rows(w)
    r = rank(rows)
    if r != 1:
        return "nonlattice", None
    seed = next(v for v in rows if any(v))
    g = 0
    for v in seed:
        g = gcd(g, abs(v))
    unit = tuple(v // g for v in seed)
    pos = next(i for i, b in enumerate(unit) if b)
    mult = 0
    for v in rows:
        mult = gcd(mult, abs(v[pos] // unit[pos]))
    span = Fraction(1)
    for p, ex in zip(primes, unit):
        span *= Fraction(p) ** (ex * mult)
    return "lattice", (1 / span if span < 1 else span)


# THE MASS SIDE


def classes(w):
    out = {}
    for x in w:
        out[x] = out.get(x, 0) + 1
    return sorted(out.items(), key=lambda t: (-t[1], -t[0]))


def log_weights(w):
    return [
        ((Decimal(x.denominator) / Decimal(x.numerator)).ln(), g) for x, g in classes(w)
    ]


def dirichlet_root(lw):
    lo, hi = 1e-9, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sum(g * e ** (-mid * float(l)) for l, g in lw) > 1:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def stopping_table(lw, amax):
    rows = []
    caps = [int(amax / float(l)) + 2 for l, _ in lw]

    def walk(i, acc, ks):
        if i == len(lw):
            total, rest = 1, sum(ks)
            for k, (_, g) in zip(ks, lw):
                total *= comb(rest, k) * g**k
                rest -= k
            rows.append((float(acc), total))
            return
        for k in range(caps[i] + 1):
            nxt = acc + k * lw[i][0]
            if float(nxt) > amax:
                break
            walk(i + 1, nxt, ks + [k])

    walk(0, Decimal(0), [])
    rows.sort()
    ws, cum, acc = [], [], 0
    for x, c in rows:
        acc += c
        if ws and ws[-1] == x:
            cum[-1] = acc
        else:
            ws.append(x)
            cum.append(acc)
    return ws, cum


def oscillation(ws, cum, delta, window, points=3000):
    out = []
    for a in MASS_A:
        vals = []
        for i in range(points):
            x = a + i * window / (points - 1)
            k = bisect.bisect_right(ws, x)
            if k:
                vals.append(log(cum[k - 1]) - delta * x)
        out.append(max(vals) - min(vals))
    return out


# THE LENGTH SIDE


def integer_weights(w):
    d = 1
    for x in w:
        d = d * x.denominator // gcd(d, x.denominator)
    return [int(x * d) for x in w], d


def mass_within(level, nums):
    pts = [(0, 0, 1)]
    for _ in range(level):
        nxt = []
        for i, j, m in pts:
            for k, (a, b) in enumerate(CELLS):
                nxt.append((Q * i + a, Q * j + b, m * nums[k]))
        pts = nxt
    side = Q**level
    top = isqrt(2 * side * side) + 2
    hist = [0] * (top + 2)
    for i, j, m in pts:
        s = i * i + j * j
        r = isqrt(s)
        if r * r < s:
            r += 1
        hist[r] += m
    out, acc = [0] * (top + 1), 0
    for r in range(top + 1):
        acc += hist[r]
        out[r] = acc
    return out


def fold(rs, g, bins=BINS):
    b = [[] for _ in range(bins)]
    for r, v in zip(rs, g):
        u = log(r, Q)
        b[int((u - int(u)) * bins) % bins].append(v)
    return [sum(x) / len(x) if x else None for x in b]


def bar(rs, g, cut):
    f1, f2 = fold(rs[:cut], g[:cut]), fold(rs[cut:], g[cut:])
    shared = [abs(a - b) for a, b in zip(f1, f2) if a is not None and b is not None]
    return max(shared), len(shared)


def ripple(w):
    nums, _den = integer_weights(w)
    M = mass_within(LEVEL, nums)
    prev = mass_within(LEVEL - 1, nums)
    lim = len(prev) - 1
    exact = all(
        M[r] == nums[0] * prev[min(r, lim)] for r in range(1, 2 * Q ** (LEVEL - 1))
    )
    alpha = -log(float(w[0]), Q)
    rs = list(range(Q**FOLD_LO, Q**FOLD_HI + 1))
    g = [log(M[r]) - alpha * log(r) for r in rs]
    f = fold(rs, g)
    f = [x - sum(f) / len(f) for x in f]
    drift, shared = bar(rs, g, len(rs) // 2)
    period, _ = bar(rs, g, Q ** ((FOLD_LO + FOLD_HI) // 2) - Q**FOLD_LO)
    decade = []
    for k in range(FOLD_LO, FOLD_HI):
        decade.append(
            max(
                abs((log(M[Q * r]) - alpha * log(Q * r)) - (log(M[r]) - alpha * log(r)))
                for r in range(Q**k, Q ** (k + 1))
            )
        )
    return alpha, f, max(f) - min(f), drift, shared, period, decade, exact


# SAFE ROOTS


def floor_root(num, den, n, digits=7):
    scale = 10**digits
    lo = int(e ** ((log(num) - log(den)) / n) * scale) - 4
    while (lo + 1) ** n * den <= num * scale**n:
        lo += 1
    while lo**n * den > num * scale**n:
        lo -= 1
    return "%d.%0*d" % (lo // scale, digits, lo % scale)


def ceil_root(num, den, n, digits=7):
    scale = 10**digits
    hi = int(e ** ((log(num) - log(den)) / n) * scale) + 4
    while (hi - 1) ** n * den >= num * scale**n:
        hi -= 1
    while hi**n * den < num * scale**n:
        hi += 1
    return "%d.%0*d" % (hi // scale, digits, hi % scale)


# EXACT ARITHMETIC IN Q(sqrt ROOT)

SQRT_PREC = 10**20
SQRT_CACHE = {}
ZERO = (Fraction(0), Fraction(0))


def rat_sqrt_down(x):
    if x <= 0:
        return Fraction(0)
    n, d = x.numerator, x.denominator
    return Fraction(isqrt(n * d * SQRT_PREC * SQRT_PREC), d * SQRT_PREC)


def rat_sqrt_up(x):
    if x <= 0:
        return Fraction(0)
    n, d = x.numerator, x.denominator
    return Fraction(isqrt(n * d * SQRT_PREC * SQRT_PREC) + 1, d * SQRT_PREC)


def root_bounds(root):
    if root not in SQRT_CACHE:
        SQRT_CACHE[root] = (rat_sqrt_down(Fraction(root)), rat_sqrt_up(Fraction(root)))
    return SQRT_CACHE[root]


def sadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def ssub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def smul(x, y, root):
    return (x[0] * y[0] + root * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def sbounds(x, root):
    a, b = x
    if b == 0:
        return a, a
    lo, hi = root_bounds(root)
    if b > 0:
        return a + b * lo, a + b * hi
    return a + b * hi, a + b * lo


def mat_mul(A, B, root):
    n = len(A)
    out = []
    for i in range(n):
        row = []
        for j in range(n):
            acc = ZERO
            for k in range(n):
                acc = sadd(acc, smul(A[i][k], B[k][j], root))
            row.append(acc)
        out.append(tuple(row))
    return tuple(out)


def rho_lower(P, root):
    if len(P) == 1:
        lo, hi = sbounds(P[0][0], root)
        return lo if lo >= 0 else (-hi if hi <= 0 else Fraction(0))
    t = sadd(P[0][0], P[1][1])
    dd = ssub(smul(P[0][0], P[1][1], root), smul(P[0][1], P[1][0], root))
    disc = ssub(smul(t, t, root), (4 * dd[0], 4 * dd[1]))
    tlo, thi = sbounds(t, root)
    tabs = tlo if tlo >= 0 else (-thi if thi <= 0 else Fraction(0))
    dlo, dhi = sbounds(disc, root)
    if dlo >= 0:
        return (tabs + rat_sqrt_down(dlo)) / 2
    if dhi <= 0:
        return rat_sqrt_down(max(Fraction(0), sbounds(dd, root)[0]))
    return tabs / 2


def norm_upper(P, root):
    if len(P) == 1:
        lo, hi = sbounds(P[0][0], root)
        return max(abs(lo), abs(hi))
    G = [[ZERO, ZERO], [ZERO, ZERO]]
    for i in range(2):
        for j in range(2):
            acc = ZERO
            for k in range(2):
                acc = sadd(acc, smul(P[k][i], P[k][j], root))
            G[i][j] = acc
    t = sadd(G[0][0], G[1][1])
    dd = ssub(smul(G[0][0], G[1][1], root), smul(G[0][1], G[1][0], root))
    thi = sbounds(t, root)[1]
    dlo = max(Fraction(0), sbounds(dd, root)[0])
    disc = max(Fraction(0), thi * thi - 4 * dlo)
    return rat_sqrt_up((thi + rat_sqrt_up(disc)) / 2)


def gripenberg(mats, root, maxlen):
    words = [((k,), m) for k, m in enumerate(mats)]
    best, arg = "0.0000000", ()
    for n in range(1, maxlen + 1):
        for word, P in words:
            r = rho_lower(P, root)
            if r <= 0:
                continue
            s = floor_root(r.numerator, r.denominator, n)
            if float(s) > float(best):
                best, arg = s, word
        if n == maxlen:
            break
        words = [
            (word + (k,), mat_mul(m, P, root))
            for word, P in words
            for k, m in enumerate(mats)
        ]
    return best, arg


def norm_bracket(mats, root, maxlen):
    prods, best, arg = list(mats), None, 0
    for n in range(1, maxlen + 1):
        u = max(norm_upper(P, root) for P in prods)
        s = ceil_root(u.numerator, u.denominator, n)
        if best is None or float(s) < float(best):
            best, arg = s, n
        if n == maxlen:
            break
        prods = [mat_mul(m, P, root) for P in prods for m in mats]
    return best, arg


def dl_matrices(c):
    N = len(c) - 1
    cc = lambda k: c[k] if 0 <= k <= N else ZERO
    T0 = tuple(
        tuple(cc(2 * i - j - 1) for j in range(1, N + 1)) for i in range(1, N + 1)
    )
    T1 = tuple(tuple(cc(2 * i - j) for j in range(1, N + 1)) for i in range(1, N + 1))
    return T0, T1


def restrict(T, root):
    N = len(T)
    TB = [[ssub(T[i][j], T[i][j + 1]) for j in range(N - 1)] for i in range(N)]
    C = []
    acc = [ZERO] * (N - 1)
    for i in range(N - 1):
        acc = [sadd(acc[j], TB[i][j]) for j in range(N - 1)]
        C.append(tuple(acc))
    ok = all(sadd(acc[j], TB[N - 1][j]) == ZERO for j in range(N - 1))
    return tuple(C), ok


# THE MULTIFRACTAL SPECTRUM


def partition(w, s):
    return sum(float(x) ** s for x in w)


def tau(w, s):
    return log(partition(w, s), Q)


def alpha_of(w, s):
    z = partition(w, s)
    return -sum(float(x) ** s * log(float(x)) for x in w) / (z * log(Q))


def legendre(w, a):
    lo, hi = -300.0, 300.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if alpha_of(w, mid) > a:
            lo = mid
        else:
            hi = mid
    s = 0.5 * (lo + hi)
    return tau(w, s) + s * a


def coarse(w, level):
    cls = classes(w)
    out = {}

    def walk(i, rest, ks):
        if i == len(cls) - 1:
            ks = ks + [rest]
            mass, cnt, left = Fraction(1), 1, level
            for k, (x, g) in zip(ks, cls):
                mass *= x**k
                cnt *= comb(left, k) * g**k
                left -= k
            hit = out.setdefault(mass, [0, []])
            hit[0] += cnt
            hit[1].append(tuple(ks))
            return
        for k in range(rest + 1):
            walk(i + 1, rest - k, ks + [k])

    walk(0, level, [])
    return cls, out


def lnfac(n):
    if n < 1:
        return 0.0
    return 0.5 * log(2 * pi * n) + n * log(n) - n + 1.0 / (12 * n) - 1.0 / (360 * n**3)


def stirling_count(level, ks, cls):
    v = lnfac(level)
    for k, (_, g) in zip(ks, cls):
        v += -lnfac(k) + k * log(g)
    return v


# THE RUN

print("WEIGHTED DESIGNS")
print()
print("== THE OBJECT ==")
print(
    "base q %d, dimension D %d, cells %s, level %d"
    % (Q, D, " ".join(str(c) for c in CELLS), LEVEL)
)
print(
    "weights %s, a probability vector over a common denominator"
    % ", ".join(str(x) for x in WEIGHT)
)
nums, den = integer_weights(WEIGHT)
print(
    "integer masses %s over %d, so every level-%d cell mass is an integer over %d^%d"
    % (nums, den, LEVEL, den, LEVEL)
)
check("weights sum to one", sum(WEIGHT), Fraction(1))
print()

print("== THE ARITHMETIC CLASS ==")
print(
    "positive rational weights only: a weight has a prime-exponent vector only there, and the logs of the primes are Q-independent"
)
print(
    "under that hypothesis the group generated by the log-weights is cyclic iff the prime-exponent matrix has rank 1"
)
CLASSOF = {}
for name, w in RENEWAL:
    check(
        "%s is a vector of rationals" % name,
        all(isinstance(x, Fraction) for x in w),
        True,
    )
    primes, rows = exponent_rows(w)
    kind, span = arithmetic_class(w)
    CLASSOF[name] = kind
    print(
        "  %-12s primes %-10s rows %-28s rank %d  %s%s"
        % (
            name,
            primes,
            rows,
            rank(rows),
            kind,
            "" if span is None else ", span ln %s" % span,
        )
    )
print()

print("== THE MASS SIDE ==")
print("delta is the root of sum_f w_f^s = 1; N(a) counts words of mass at least e^(-a)")
OSC = {}
for name, w in RENEWAL:
    lw = log_weights(w)
    delta = dirichlet_root(lw)
    ws, cum = stopping_table(lw, MASS_MAX)
    osc = oscillation(ws, cum, delta, log(Q))
    OSC[name] = osc
    spread = (max(osc) - min(osc)) / max(osc)
    falling = all(x > y for x, y in pairwise(osc)) and osc[-1] / osc[0] < 0.8
    trend = "flat" if spread < 0.01 else ("decaying" if falling else "mixed")
    want = "flat" if CLASSOF[name] == "lattice" else "decaying"
    check("trend %s" % name, trend, want)
    check("delta %s" % name, "%.10f" % delta, NAPKIN_MASS[name][0])
    check(
        "oscillation %s" % name,
        "  ".join("%.5f" % v for v in osc),
        NAPKIN_MASS[name][1],
    )
    print(
        "  %-12s delta %.10f  %-11s osc of ln N(a) - delta a over a window of ln %d at a = %s: %s  %s"
        % (
            name,
            delta,
            CLASSOF[name],
            Q,
            ",".join(str(a) for a in MASS_A),
            "  ".join("%.5f" % v for v in osc),
            trend,
        )
    )
print(
    "for a probability vector s = 1 solves the Dirichlet equation, so delta = 1 identically;"
)
print("what weights move on the mass side is the arithmetic class alone")
for name, w in PROBS:
    check(
        "delta one %s" % name, "%.10f" % dirichlet_root(log_weights(w)), "1.0000000000"
    )
print()

print("== THE LENGTH SIDE ==")
print(
    "M(r) is the mass of cells within radius r of the corner fixed point, detrended by"
)
print(
    "alpha = -log_q w_0 and folded into %d bins of log_%d r over whole periods of R in %d^%d..%d^%d"
    % (BINS, Q, Q, FOLD_LO, Q, FOLD_HI)
)
BASE, RIP = None, {}
for name, w in PROBS:
    alpha, f, swing, drift, shared, period, decade, exact = ripple(w)
    if BASE is None:
        BASE, gap = f, 0.0
    else:
        gap = max(abs(a - b) for a, b in zip(f, BASE))
    RIP[name] = (swing, drift, gap)
    check("swing %s" % name, "%.6f" % swing, NAPKIN_LENGTH[name][0])
    check("drift %s" % name, "%.6f" % drift, NAPKIN_LENGTH[name][1])
    check("0/1 gap %s" % name, "%.6f" % gap, NAPKIN_LENGTH[name][2])
    check(
        "decade %s" % name,
        "  ".join("%.4f" % v for v in decade),
        NAPKIN_LENGTH[name][3],
    )
    check("ripple kept %s" % name, swing > 10 * drift, True)
    check("ripple kept on whole periods %s" % name, swing > 10 * period, True)
    check("cross level %s" % name, exact, True)
    print(
        "  %-12s alpha %.9f  swing %.6f on drift bar %.6f  gap against equal weights %.6f"
        % (name, alpha, swing, drift, gap)
    )
    print(
        "  %-12s log_%d periodicity residual by decade: %s   cross-level identity exact %s"
        % ("", Q, "  ".join("%.4f" % v for v in decade), exact)
    )
    print(
        "  %-12s drift bar on %d shared bins, %.6f on a whole-period split of the window"
        % ("", shared, period)
    )
print(
    "equal weights give every cell mass 1, so M(r) is the cell count and the fold is the 0/1 ripple"
)
print()

print("== THE LOCAL DIMENSIONS ==")
for name, w in PROBS:
    amin = -log(float(max(w)), Q)
    amax = -log(float(min(w)), Q)
    ainf = -sum(float(x) * log(float(x)) for x in w) / log(Q)
    print(
        "  %-12s alpha_min %.9f  alpha_max %.9f  alpha at s = 1 %.9f"
        % (name, amin, amax, ainf)
    )
print()

print("== THE JSR LADDER ==")
print(
    "rung 0, the designs: supp phi sits in the unit cell, so the Daubechies-Lagarias matrices"
)
print(
    "are 1 x 1, T_f = [q^D w_f], and the joint spectral radius is q^D max_f w_f in closed form"
)
for name, w in PROBS:
    T = [Q**D * x for x in w]
    jsr = max(T)
    mats = [(((x, Fraction(0)),),) for x in T]
    lo, word = gripenberg(mats, 0, GRIPENBERG_LEN)
    hi, at = norm_bracket(mats, 0, GRIPENBERG_LEN)
    holder = -D - log(float(max(w)), Q)
    print(
        "  %-12s T_f = %-28s JSR %-6s solver [%s, %s]  alpha_Holder(phi) %.9f"
        % (name, "[" + ", ".join(str(x) for x in T) + "]", str(jsr), lo, hi, holder)
    )
    check(
        "rung 0 matrices %s" % name,
        "[" + ", ".join(str(x) for x in T) + "]",
        NAPKIN_RUNG0[name][0],
    )
    check("rung 0 jsr %s" % name, str(jsr), NAPKIN_RUNG0[name][1])
    check(
        "rung 0 alpha_min %s" % name,
        "%.9f" % -log(float(max(w)), Q),
        NAPKIN_RUNG0[name][2],
    )
    check(
        "rung 0 alpha_max %s" % name,
        "%.9f" % -log(float(min(w)), Q),
        NAPKIN_RUNG0[name][3],
    )
    check("rung 0 Holder %s" % name, "%.9f" % holder, NAPKIN_RUNG0[name][4])
    check("rung 0 lower %s" % name, float(lo) <= float(jsr), True)
    check("rung 0 upper %s" % name, float(hi) >= float(jsr), True)
    check("rung 0 collapse %s" % name, float(hi) - float(lo) < 1e-6, True)
print()
print(
    "rung 1, the first overlap at q = 2: the same solver on T_0 = (c_(2i-j-1)), T_1 = (c_(2i-j)),"
)
print(
    "restricted to sum v_i = 0, Gripenberg words to length %d and norm upper bounds to length %d"
    % (GRIPENBERG_LEN, NORM_LEN)
)
MASKS = (
    (
        "hat (1/2, 1, 1/2)",
        0,
        (
            (Fraction(1, 2), Fraction(0)),
            (Fraction(1), Fraction(0)),
            (Fraction(1, 2), Fraction(0)),
        ),
    ),
    (
        "D4",
        3,
        (
            (Fraction(1, 4), Fraction(1, 4)),
            (Fraction(3, 4), Fraction(1, 4)),
            (Fraction(3, 4), Fraction(-1, 4)),
            (Fraction(1, 4), Fraction(-1, 4)),
        ),
    ),
)
SQ3 = (Fraction(1, 4), Fraction(1, 4))
D4_JSR = "0.6830127"
D4_ALPHA = "0.5500157"
print(
    "the hat mask carries alpha = 1 exactly; D4's lower end is the closed form (1 + sqrt 3)/4 and"
)
print(
    "its Holder exponent 2 - log_2(1 + sqrt 3) = %s must sit inside the printed bracket"
    % D4_ALPHA
)
for name, root, c in MASKS:
    even = ZERO
    odd = ZERO
    for i, x in enumerate(c):
        if i % 2 == 0:
            even = sadd(even, x)
        else:
            odd = sadd(odd, x)
    check("mask %s even" % name, even, (Fraction(1), Fraction(0)))
    check("mask %s odd" % name, odd, (Fraction(1), Fraction(0)))
    T0, T1 = dl_matrices(c)
    V0, ok0 = restrict(T0, root)
    V1, ok1 = restrict(T1, root)
    check("invariant subspace %s" % name, ok0 and ok1, True)
    lo, word = gripenberg([V0, V1], root, GRIPENBERG_LEN)
    hi, at = norm_bracket([V0, V1], root, NORM_LEN)
    alo = safe_log2(Fraction(hi), 7, ROUND_FLOOR)
    ahi = safe_log2(Fraction(lo), 7, ROUND_CEILING)
    check("rung 1 lower %s" % name, lo, NAPKIN_RUNG1[name][0])
    check("rung 1 upper %s" % name, hi, NAPKIN_RUNG1[name][1])
    check("rung 1 alpha lower %s" % name, alo, NAPKIN_RUNG1[name][2])
    check("rung 1 alpha upper %s" % name, ahi, NAPKIN_RUNG1[name][3])
    check("rung 1 bracket order %s" % name, float(lo) <= float(hi), True)
    print(
        "  %-18s dim %d  JSR in [%s, %s], the lower at word %s, the upper at length %d"
        % (name, len(V0), lo, hi, "".join(str(k) for k in word), at)
    )
    print("  %-18s alpha = -log_2 JSR in [%s, %s]" % ("", alo, ahi))
    if name == "D4":
        check(
            "D4 lower is (1 + sqrt 3)/4",
            lo,
            floor_root(sbounds(SQ3, 3)[0].numerator, sbounds(SQ3, 3)[0].denominator, 1),
        )
        check("D4 lower reads", lo, D4_JSR)
        closed = safe_log2(sbounds(SQ3, 3)[1], 7, ROUND_CEILING)
        check("D4 closed form", closed, D4_ALPHA)
        check(
            "D4 closed form inside the bracket",
            float(alo) <= float(closed) <= float(ahi),
            True,
        )
print()

print("== THE CONTROLS ==")
print(
    "pre-registered: equal weights reproduce the 0/1 design exactly; log-commensurable weights"
)
print("keep every ripple; only an irrational log ratio may move a mass observable")
print(
    "  %-12s %-11s %-9s %-9s %-8s %-9s %-9s %s"
    % (
        "weights",
        "class",
        "swing",
        "drift",
        "0/1 gap",
        "osc first",
        "osc last",
        "verdict",
    )
)
moved = 0
for name in CONTROLS:
    swing, drift, gap = RIP[name]
    osc = OSC[name]
    kept = swing > 10 * drift
    flat = osc[-1] / osc[0] > 0.9
    if not flat:
        moved += 1
    verdict = "length kept, mass %s" % ("flat" if flat else "moved")
    print(
        "  %-12s %-11s %-9.6f %-9.6f %-8.6f %-9.5f %-9.5f %s"
        % (name, CLASSOF[name], swing, drift, gap, osc[0], osc[-1], verdict)
    )
    check("control ripple %s" % name, kept, True)
check(
    "equal weights reproduce the 0/1 fold", "%.6f" % RIP["1/3,1/3,1/3"][2], "0.000000"
)
check("exactly one control moves the mass side", moved, 1)
print(
    "length observable moves under none of the three; mass observable moves under one of the three"
)
print()

print("== THE MULTIFRACTAL SPECTRUM ==")
print("Moran self-similar measure, ratios all 1/q, open set condition: tau(s) solves")
print(
    "sum_f w_f^s (1/q)^tau(s) = 1, hence tau(s) = log_q sum_f w_f^s, and f(alpha) = inf_s (alpha s + tau(s))"
)
print("  s   sum_f w_f^s        tau(s)        alpha(s)      f(alpha(s))")
for s in (-2, -1, 0, 1, 2, 3):
    zr = sum(x**s for x in WEIGHT)
    a = alpha_of(WEIGHT, s)
    print(
        "  %2d  %-18s %-13.9f %-13.9f %.9f"
        % (s, str(zr), tau(WEIGHT, s), a, tau(WEIGHT, s) + s * a)
    )
check(
    "tau at 0 is the box dimension",
    "%.12f" % tau(WEIGHT, 0),
    "%.12f" % log(len(WEIGHT), Q),
)
check("tau at 1 vanishes", "%.12f" % tau(WEIGHT, 1), "%.12f" % 0.0)
a1 = alpha_of(WEIGHT, 1)
check(
    "f is tangent to the diagonal at s = 1", "%.9f" % legendre(WEIGHT, a1), "%.9f" % a1
)
amin = -log(float(max(WEIGHT)), Q)
amax = -log(float(min(WEIGHT)), Q)
print(
    "alpha ranges over [%.9f, %.9f]; f(alpha_min) = log_q #argmax w = %.9f, f(alpha_max) = %.9f"
    % (
        amin,
        amax,
        log(sum(1 for x in WEIGHT if x == max(WEIGHT)), Q),
        log(sum(1 for x in WEIGHT if x == min(WEIGHT)), Q),
    )
)
print()
print(
    "the level-L box partition carries the moments exactly: sum_i mu_i^s = (sum_f w_f^s)^L"
)
for level in SPECTRUM_LEVELS:
    cls, groups = coarse(WEIGHT, level)
    for s in (-2, -1, 0, 1, 2, 3):
        got = sum(cnt * mass**s for mass, (cnt, _) in groups.items())
        check(
            "partition identity L %d s %d" % (level, s),
            got,
            sum(x**s for x in WEIGHT) ** level,
        )
    print(
        "  level %2d  %d distinct box masses, %d boxes, moments exact at s = -2..3"
        % (level, len(groups), sum(c for c, _ in groups.values()))
    )
print()
print(
    "the coarse-grained spectrum f_L(alpha) = log_q N(alpha) / L against the Legendre transform"
)
print(
    "f_L <= f at every level and every achievable alpha, since N_i mu_i^s <= q^(L tau(s)); the gap"
)
print("is the Stirling volume term of the multinomial count and falls like log L / L")
print(
    "  level  alpha         f_L(alpha)    f(alpha)      deficit    max deficit  max |f_L - Stirling|"
)
DEF = {}
for level in SPECTRUM_LEVELS:
    cls, groups = coarse(WEIGHT, level)
    worst, stir, mid = 0.0, 0.0, None
    for mass, (cnt, ks) in sorted(groups.items()):
        a = -log(float(mass)) / (level * log(Q))
        fl = log(cnt) / (level * log(Q))
        fx = legendre(WEIGHT, a)
        check("f_L under the transform L %d" % level, fl <= fx + 1e-12, True)
        worst = max(worst, fx - fl)
        if len(ks) == 1:
            stir = max(
                stir, abs(fl - stirling_count(level, ks[0], cls) / (level * log(Q)))
            )
        if ks[0][0] * 2 == level:
            mid = (a, fl, fx)
        if len(ks) == 1 and 0 in ks[0]:
            check("endpoint exact L %d" % level, "%.9f" % abs(fx - fl), "0.000000000")
    DEF[level] = worst
    print(
        "  %5d  %-13.9f %-13.9f %-13.9f %-10.6f %-12.6f %.6f"
        % (level, mid[0], mid[1], mid[2], mid[2] - mid[1], worst, stir)
    )
    check("Stirling agreement L %d" % level, stir < 1e-3, True)
check("the deficit falls with the level", DEF[6] > DEF[8] > DEF[10], True)
print()
print(
    "the whole coarse-grained band at level %d, printed row by row"
    % SPECTRUM_LEVELS[-1]
)
print("  k   box mass                  boxes    alpha         f_L(alpha)    f(alpha)")
cls, groups = coarse(WEIGHT, SPECTRUM_LEVELS[-1])
level = SPECTRUM_LEVELS[-1]
for mass, (cnt, ks) in sorted(groups.items(), reverse=True):
    a = -log(float(mass)) / (level * log(Q))
    print(
        "  %2d  %-25s %-8d %-13.9f %-13.9f %.9f"
        % (
            ks[0][0],
            str(mass),
            cnt,
            a,
            log(cnt) / (level * log(Q)),
            legendre(WEIGHT, a),
        )
    )
top = max(log(c) / (level * log(Q)) for c, _ in groups.values())
print(
    "the band's top %.9f rises to tau(0) = log_q |F| = %.9f as the level grows"
    % (top, tau(WEIGHT, 0))
)
check("the band's top sits under the box dimension", top < tau(WEIGHT, 0), True)
print()

print("== THE VERDICT ==")
if FAILURES:
    for line in FAILURES:
        print("FAIL " + line)
    raise SystemExit(1)
print("every assertion holds")
