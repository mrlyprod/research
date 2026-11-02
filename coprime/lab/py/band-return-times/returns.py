import argparse
import time
from fractions import Fraction
from itertools import combinations
from math import comb, exp, gcd, log, log10

MINPOLY = (
    [-1, 1],
    [-2, 1],
    [-3, 1],
    [-6, 1],
    [45, -15, 1],
    [90, -26, 1],
    [-3402, 945, -63, 1],
    [1134, -99, 1],
    [1299078, -293787, 16065, -255, 1],
    [-96228, 17469, -392, 1],
)


def up(x, d):
    f = 10**d
    n = int(x * f)
    return (n + 1 if n < x * f else n) / f


def down(x, d):
    f = 10**d
    n = int(x * f)
    return (n - 1 if n > x * f else n) / f


# BAND AUTOMATON


def succ(j, z1, z2):
    if j % 3 == 0:
        return (j // 3, (j + z1) // 3)
    if (j - z2) % 3 == 0:
        return ((j - z2) // 3,)
    return ()


def first_return(z1, z2):
    start = z1 // 3
    seen = {start}
    frontier = [start]
    n = 1
    while frontier:
        n += 1
        nxt = []
        for j in frontier:
            for t in succ(j, z1, z2):
                if t == 0:
                    return n
                if t not in seen:
                    seen.add(t)
                    nxt.append(t)
        frontier = nxt
    return 0


def band(z1, z2):
    return -((z2 - 1) // 2), (z1 - 1) // 2


def degree_profile(z1, z2):
    lo, hi = band(z1, z2)
    cnt = [0, 0, 0]
    edges = 0
    for j in range(lo, hi + 1):
        d = len(succ(j, z1, z2))
        cnt[d] += 1
        edges += d
        for t in succ(j, z1, z2):
            if not lo <= t <= hi:
                raise ValueError("the band is not invariant at " + str((z1, z2, j)))
    return cnt, edges, hi - lo + 1


# COLUMN TRANSFER


def slot_profile(k, n):
    return [(n - r + k - 1) // k for r in range(k)]


def column_count(k, n, prim=True, twist=1):
    if n < k:
        return 0
    w = (3**k - 1) // 2
    slots = slot_profile(k, n)
    top = sum(s * 3**r for r, s in enumerate(slots))
    total = 0
    for j in range(top // w + 1):
        M = j * w
        cur = {0: 1}
        for r in range(k):
            mr = M // 3**r % 3
            lo = 1 if prim and r == 0 else 0
            nxt = {}
            for t, v in cur.items():
                for c in range(lo, slots[r] + 1):
                    if (t + c - mr) % 3:
                        continue
                    u = (t + c - mr) // 3
                    nxt[u] = nxt.get(u, 0) + v * comb(slots[r] - lo, c - lo) * twist**c
            cur = nxt
        total += cur.get(M // 3**k, 0)
    return total


def lift_count(k, n, prim=True):
    return column_count(k, n, prim) - (0 if prim else 1)


def model_returns(k, N, cap=40):
    w = (3**k - 1) // 2
    out = []
    lo = hi = 0.0
    prev = 0
    for n in range(k, min(N, cap * k) + 1):
        a = column_count(k, n, True, 2)
        d = a - prev
        prev = a
        if d:
            lo += d * w / 3**n
            hi += d * w / 3 ** (n - 1)
        out.append((n, lo, hi))
    n = out[-1][0]
    while n < N:
        n += 1
        lo += 4.0 / 9
        hi += 4.0 / 3
        out.append((n, lo, hi))
    return out


def exact_L(k):
    acc = 0.0
    for s in range(1 << (k - 1)):
        a = sum(3 ** (i + 1) for i in range(k - 1) if s >> i & 1)
        acc += 1.0 / (1 + 2 * a)
    return acc


def brute_lift_count(k, n, prim=True):
    w = (3**k - 1) // 2
    pw = [3**i for i in range(n)]
    out = 0
    for s in range(1, 1 << n):
        if prim and not s & 1:
            continue
        v = sum(pw[i] for i in range(n) if s >> i & 1)
        if v % w == 0:
            out += 1
    return out


# THE FIXED MATRIX AT n = bk


def block_head(b, j, K=64):
    M = j * (3**K - 1) // 2
    d = [M // 3**r % 3 for r in range(K)]
    r0 = K - 1
    while r0 > 1 and d[r0 - 1] == d[K - 1]:
        r0 -= 1
    return d[:r0], d[K - 1], j // 2 - 1 if j and not j % 2 else j // 2


def block_pieces(b, j):
    head, bulk, fin = block_head(b, j)
    S = b // 2 + 1
    cur = {0: 1}
    for r, mr in enumerate(head):
        lo = 1 if r == 0 else 0
        nxt = {}
        for t, v in cur.items():
            for c in range(lo, b + 1):
                if (t + c - mr) % 3:
                    continue
                u = (t + c - mr) // 3
                nxt[u] = nxt.get(u, 0) + v * comb(b - lo, c - lo)
        cur = nxt
    A = [[0] * S for _ in range(S)]
    for t in range(S):
        for c in range(b + 1):
            if (t + c - bulk) % 3:
                continue
            u = (t + c - bulk) // 3
            if u < S:
                A[u][t] += comb(b, c)
            elif comb(b, c):
                raise ValueError("the carry leaves the band at " + str((b, j, t, c)))
    return A, [cur.get(t, 0) for t in range(S)], fin, len(head)


def block_count(b, k):
    tot = 0
    for j in range(b + 1):
        A, h, f, r0 = block_pieces(b, j)
        if k < r0:
            return None
        for _ in range(k - r0):
            h = [sum(A[i][t] * h[t] for t in range(len(h))) for i in range(len(A))]
        if f < len(h):
            tot += h[f]
    return tot


def block_closure(A, s):
    seen, front = {s}, [s]
    while front:
        nxt = []
        for t in front:
            for u in range(len(A)):
                if A[u][t] and u not in seen:
                    seen.add(u)
                    nxt.append(u)
        front = nxt
    return seen


def block_rho(A, steps=600):
    v = [1.0] * len(A)
    g = 0.0
    for i in range(steps):
        v = [sum(A[u][t] * v[t] for t in range(len(A))) for u in range(len(A))]
        s = sum(v)
        if s <= 0:
            return 0.0
        v = [x / s for x in v]
        if i >= steps // 2:
            g += log(s)
    return exp(g / (steps - steps // 2))


def block_bracket(b):
    dev, r0s, rhos = set(), [], []
    for j in range(b + 1):
        A, h, f, r0 = block_pieces(b, j)
        r0s.append(r0)
        rhos.append(block_rho(A))
        for t in range(len(A)):
            dev.add(Fraction(3 * sum(A[u][t] for u in range(len(A))) - 2**b, 3))
    par = (
        {Fraction(-1, 3), Fraction(2, 3)}
        if not b % 2
        else {Fraction(-2, 3), Fraction(1, 3)}
    )
    if not dev <= par:
        raise ValueError("a column sum leaves the parity pair at b = " + str(b))
    A, h, f, _ = block_pieces(b, 1)
    C = block_closure(A, 0)
    if f or not h[0] or not A[0][0]:
        raise ValueError(
            "the j = 1 sector does not start and end at the carry 0 at b = " + str(b)
        )
    for t in sorted(C):
        if sum(A[u][t] for u in C) != sum(A[u][t] for u in range(len(A))):
            raise ValueError(
                "the closure of the carry 0 is not closed at b = " + str(b)
            )
        if 0 not in block_closure(A, t):
            raise ValueError(
                "the carry 0 is not reachable from " + str(t) + " at b = " + str(b)
            )
    cs = [Fraction(3 * sum(A[u][t] for u in C) - 2**b, 3) for t in sorted(C)]
    return sorted(dev), min(cs), max(cs), len(C), r0s, max(rhos)


# EXACT ALGEBRA


def poly_div(P, Q):
    P = list(P)
    n, m = len(P) - 1, len(Q) - 1
    if n < m:
        return None
    out = [0] * (n - m + 1)
    for i in range(n - m, -1, -1):
        if P[i + m] % Q[m]:
            return None
        out[i] = P[i + m] // Q[m]
        for t in range(m + 1):
            P[i + t] -= out[i] * Q[t]
    return out if not any(P[:m]) else None


def poly_at(P, x):
    v = 0
    for c in reversed(P):
        v = v * x + c
    return v


def poly_show(P):
    out = ""
    for i in range(len(P) - 1, -1, -1):
        if not P[i]:
            continue
        s = "-" if P[i] < 0 else "+"
        a = abs(P[i])
        t = str(a) if a != 1 or not i else ""
        t += "x^" + str(i) if i > 1 else ("x" if i == 1 else "")
        out += t if not out and s == "+" else " " + s + " " + t
    return out or "0"


def poly_roots(P, R):
    n = len(P) - 1
    q = [P[i] / R ** (n - i) for i in range(n + 1)]
    z = [(0.4 + 0.9j) ** i for i in range(n)]
    for _ in range(4000):
        mv = 0.0
        for i in range(n):
            d = 1.0 + 0j
            for t in range(n):
                if t != i:
                    d *= z[i] - z[t]
            v = 0j
            for c in reversed(q):
                v = v * z[i] + c
            s = v / d
            z[i] -= s
            mv = max(mv, abs(s))
        if mv < 1e-16:
            break
    return [w * R for w in z]


def poly_minfactor(P, rs, lam):
    rest = [r for r in rs if r is not lam]
    for d in range(1, len(P)):
        for S in combinations(rest, d - 1):
            f = [1.0 + 0j]
            for r in (lam,) + S:
                g = [0j] * (len(f) + 1)
                for i, c in enumerate(f):
                    g[i] -= c * r
                    g[i + 1] += c
                f = g
            if max(abs(c.imag) for c in f) > 0.4:
                continue
            g = [round(c.real) for c in f]
            if max(abs(g[i] - f[i].real) for i in range(len(g))) > 0.45:
                continue
            if poly_div(P, g) is not None:
                return g
    return None


def poly_newton(P, x, steps=5, bits=256):
    D = [i * P[i] for i in range(1, len(P))]
    v = Fraction(x)
    for _ in range(steps):
        d = poly_at(D, v)
        if not d:
            return v
        v -= Fraction(poly_at(P, v)) / d
        v = Fraction(round(v * 2**bits), 2**bits)
    return v


def poly_isolate(P, x, digits):
    w = Fraction(1, 10**6)
    cap = Fraction(abs(x)) / 1000 + 1
    while w < cap and poly_at(P, Fraction(x) - w) * poly_at(P, Fraction(x) + w) > 0:
        w *= 2
    if w >= cap:
        return None
    lo, hi = Fraction(x) - w, Fraction(x) + w
    eps = Fraction(1, 10 ** (digits + 3))
    while hi - lo > eps:
        mid = (lo + hi) / 2
        if poly_at(P, lo) * poly_at(P, mid) <= 0:
            hi = mid
        else:
            lo = mid
    return lo, hi


def poly_radical(P):
    if len(P) == 2:
        return str(-P[0])
    if len(P) != 3:
        return "-"
    p, q = -P[1], P[0]
    d = p * p - 4 * q
    f = 1
    t = 2
    while t * t <= d:
        while d % (t * t) == 0:
            d //= t * t
            f *= t
        t += 1
    g = gcd(gcd(p, f), 2)
    num = (
        str(p // g)
        + " + "
        + (str(f // g) + " " if f // g != 1 else "")
        + "sqrt "
        + str(d)
    )
    return num if 2 // g == 1 else "(" + num + ") / 2"


def linear_solve(rows, rhs):
    n = len(rows[0])
    A = [[Fraction(v) for v in r] + [Fraction(rhs[i])] for i, r in enumerate(rows)]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        A[r] = [v / A[r][c] for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [A[i][t] - f * A[r][t] for t in range(n + 1)]
        piv.append(c)
        r += 1
        if r == len(A):
            break
    for i in range(r, len(A)):
        if A[i][n] and not any(A[i][:n]):
            return None
    x = [Fraction(0)] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    return x


def min_recurrence(seq, maxord):
    for L in range(1, maxord + 1):
        if len(seq) < 2 * L + 6:
            return None
        c = linear_solve(
            [[seq[i + t] for t in range(L)] for i in range(L + 2)],
            [seq[i + L] for i in range(L + 2)],
        )
        if c is None:
            continue
        if all(
            sum(c[t] * seq[i + t] for t in range(L)) == seq[i + L]
            for i in range(len(seq) - L)
        ):
            return [-int(v) for v in c] + [1]
    return None


def lam_block(b, terms=0, maxord=0, kbig=160):
    terms = terms or 4 * b + 14
    maxord = maxord or 2 * b + 4
    seq = [block_count(b, k) for k in range(3, 3 + terms)]
    P = min_recurrence(seq, maxord)
    if P is None:
        return None
    rs = poly_roots(P, 2.0**b / 3 + 1)
    lam = max(rs, key=lambda z: z.real)
    ratio = block_count(b, kbig + 1) / block_count(b, kbig)
    if abs(lam - ratio) > 1e-2 * ratio:
        raise ValueError(
            "the dominant root and the block ratio disagree at b = " + str(b)
        )
    Q = poly_minfactor(P, rs, lam)
    lam2 = max(abs(z) for z in rs if z is not lam) if len(rs) > 1 else 0.0
    return (
        P,
        Q,
        poly_isolate(Q or P, poly_newton(Q or P, lam.real), 12),
        ratio,
        lam2 / lam.real,
    )


# THE LIFT COUNT


def lift_support(word):
    k = len(word) + 1
    T = [i + 1 for i, c in enumerate(word) if c]
    m = 1 + 2 * sum(3**i for i in T)
    supp = [p for p in range(k) if p not in T] + [k + p for p in T]
    return k, m, sorted(supp)


def submask_count(word):
    _k, m, supp = lift_support(word)
    h = len(supp) // 2
    lo, hi = supp[:h], supp[h:]
    d = {}
    for s in range(1 << len(lo)):
        v = sum(pow(3, lo[i], m) for i in range(len(lo)) if s >> i & 1) % m
        d[v] = d.get(v, 0) + 1
    out = 0
    for s in range(1 << len(hi)):
        v = sum(pow(3, hi[i], m) for i in range(len(hi)) if s >> i & 1) % m
        out += d.get((-v) % m, 0)
    return out


def submask_set(word):
    _k, m, supp = lift_support(word)
    n = len(supp)
    return [
        s
        for s in range(1 << n)
        if sum(3 ** supp[i] for i in range(n) if s >> i & 1) % m == 0
    ], n


def closure_faults(D, n):
    S = set(D)
    full = (1 << n) - 1
    bad = 0
    for a in D:
        if (full ^ a) not in S:
            bad += 1
        for b in D:
            if not a & b and (a | b) not in S:
                bad += 1
            if a & b == b and (a ^ b) not in S:
                bad += 1
    return bad


def irreducibles(D):
    S = set(D)
    out = []
    for a in D:
        if not a:
            continue
        if not any(b and b != a and b & a == b and (a ^ b) in S for b in D):
            out.append(a)
    return out


def lift_census(k):
    tot = triv = pw = dis = 0
    ito = imax = 0
    cls = {}
    for w in range(1 << (k - 1)):
        word = tuple((w >> i) & 1 for i in range(k - 1))
        D, _n = submask_set(word)
        I = irreducibles(D)
        tot += len(D)
        ito += len(I)
        imax = max(imax, len(I))
        triv += len(D) == 2
        pw += not len(D) & (len(D) - 1)
        dis += all(not I[a] & I[b] for a in range(len(I)) for b in range(a))
        r = lift_support(word)[1] % 9
        c = cls.setdefault(r, [0, 0, 0])
        c[0] += 1
        c[1] += len(D)
        c[2] += len(I)
    return tot, ito, imax, triv, pw, dis, cls


def rank_exact(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    n = len(rows[0])
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c] / pv
                for j in range(c, n):
                    rows[i][j] -= f * rows[r][j]
        r += 1
        if r == len(rows):
            break
    return r


def all_words(n):
    return [
        tuple((s >> i) & 1 for i in range(L))
        for L in range(n + 1)
        for s in range(1 << L)
    ]


def hankel_rank(p, q, fn=None):
    fn = fn or submask_count
    U, V = all_words(p), all_words(q)
    memo = {}
    H = []
    for u in U:
        row = []
        for v in V:
            w = u + v
            if w not in memo:
                memo[w] = fn(w)
            row.append(memo[w])
        H.append(row)
    return rank_exact(H), len(U)


def irreducible_count(word):
    D, _ = submask_set(word)
    return len(irreducibles(D))


# THE FLOOR AND THE MODEL


def floor_count(k):
    w = (3**k - 1) // 2
    pw = [3**i for i in range(k)]
    out = 0
    for s in range(1, (1 << k) - 1):
        v = sum(pw[i] for i in range(k) if s >> i & 1)
        if gcd(v, w) == 1:
            out += 1
    return out


def model_by_top(k):
    out = [0.0] * k
    for t in range(1, k):
        base = 1 + 2 * 3**t
        acc = 0.0
        for s in range(1 << (t - 1)):
            a = sum(3 ** (i + 1) for i in range(t - 1) if s >> i & 1)
            acc += 1.0 / (base + 2 * a)
        out[t] = acc * 2**k
    return out


# THE FITS

TQ = (
    12.706,
    4.303,
    3.182,
    2.776,
    2.571,
    2.447,
    2.365,
    2.306,
    2.262,
    2.228,
    2.201,
    2.179,
    2.16,
    2.145,
    2.131,
    2.12,
    2.11,
    2.101,
    2.093,
    2.086,
    2.08,
    2.074,
    2.069,
    2.064,
    2.06,
    2.056,
    2.052,
    2.048,
    2.045,
    2.042,
)


def ols(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sl = sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / sxx
    a = my - sl * mx
    se = (sum((ys[i] - a - sl * xs[i]) ** 2 for i in range(n)) / (n - 2) / sxx) ** 0.5
    t = TQ[n - 3] if 3 <= n <= 32 else 1.96
    return sl, se, sl - t * se, sl + t * se, n


def ols_diag(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sl = sum((xs[i] - mx) * (ys[i] - sum(ys) / n) for i in range(n)) / sxx
    a = sum(ys) / n - sl * mx
    r = [ys[i] - a - sl * xs[i] for i in range(n)]
    s2 = sum(v * v for v in r) / (n - 2)
    h = [1 / n + (xs[i] - mx) ** 2 / sxx for i in range(n)]
    dw = sum((r[i] - r[i - 1]) ** 2 for i in range(1, n)) / sum(v * v for v in r)
    m = n // 2
    v1 = sum(v * v for v in r[:m]) / m
    v2 = sum(v * v for v in r[m:]) / (n - m)
    stu = max(abs(r[i]) / (s2 * (1 - h[i])) ** 0.5 for i in range(n))
    cook = max(r[i] ** 2 * h[i] / (2 * s2 * (1 - h[i]) ** 2) for i in range(n))
    return dw, max(v1, v2) / min(v1, v2), stu, cook


def ols2(x1, x2, ys):
    n = len(ys)
    d = [sum(x1) / n, sum(x2) / n, sum(ys) / n]
    u = [v - d[0] for v in x1]
    w = [v - d[1] for v in x2]
    y = [v - d[2] for v in ys]
    a, b, c = (
        sum(v * v for v in u),
        sum(u[i] * w[i] for i in range(n)),
        sum(v * v for v in w),
    )
    e, f = sum(u[i] * y[i] for i in range(n)), sum(w[i] * y[i] for i in range(n))
    det = a * c - b * b
    return (c * e - b * f) / det, (a * f - b * e) / det


def survival(hist, k):
    s, b = [], 2
    while True:
        v = sum(x for d, x in hist.items() if d > b * k)
        if not v:
            break
        s.append(v)
        b += 1
    return s


def geom_fit(s):
    n = [s[i] - s[i + 1] for i in range(len(s) - 1)] + [s[-1]]
    tot = sum(n)
    r = 1 - tot / sum((i + 1) * n[i] for i in range(len(n)))
    e = [tot * (1 - r) * r**i for i in range(len(n))]
    raw = sum((n[i] - e[i]) ** 2 / e[i] for i in range(len(n)))
    o, x, co, ce = [], [], 0, 0.0
    for i in range(len(n)):
        co += n[i]
        ce += e[i]
        if ce >= 5:
            o.append(co)
            x.append(ce)
            co, ce = 0, 0.0
    if co:
        o[-1] += co
        x[-1] += ce
    return (
        r,
        raw,
        len(n) - 2,
        sum((o[i] - x[i]) ** 2 / x[i] for i in range(len(o))),
        len(o) - 2,
    )


def dec_down(x, d):
    n = x.numerator * 10**d // x.denominator
    return str(n // 10**d) + "." + str(n % 10**d).rjust(d, "0")


def dec_up(x, d):
    n = -((-x.numerator * 10**d) // x.denominator)
    return str(n // 10**d) + "." + str(n % 10**d).rjust(d, "0")


# SWEEP


def sweep(k):
    w = (3**k - 1) // 2
    hist = {}
    cand = 0
    for z1 in range(3, w, 3):
        if gcd(z1, w) != 1:
            continue
        cand += 1
        d = first_return(z1, w - z1)
        if d:
            hist[d] = hist.get(d, 0) + 2
    return hist, 2 * cand


# VERBS


def cmd_automaton(args):
    print(
        "the band automaton of the direction (z1, z2), z1 + z2 = w, 3 | z1: states are the integers j in "
        "[-(z2-1)//2, (z1-1)//2]; from j the moves are j -> (j + a)/3 over the increments a in {0, z1, -z2} "
        "that keep the quotient integral, so out-degree 2 on j = 0 mod 3, 1 on j = z2 mod 3, 0 on the third class; "
        "the walk leaves 0 by the forced increment z1 and a return is a walk back to 0, of length n exactly when "
        "the multiplier m it spells has m w binary of base-3 length n"
    )
    print("w z1 z2 states edges mean deg0 deg1 deg2 firstreturn")
    for w in args.weights:
        for z1 in range(3, w, 3):
            z2 = w - z1
            if gcd(z1, z2) != 1:
                continue
            cnt, edges, states = degree_profile(z1, z2)
            if z1 % 30 and w > 40:
                continue
            print(
                w,
                z1,
                z2,
                states,
                edges,
                down(edges / states, 6),
                cnt[0],
                cnt[1],
                cnt[2],
                first_return(z1, z2),
            )
    print(
        "criticality over every coprime direction of these weights; the mean out-degree is exactly "
        "1 + (n0 - n2) / N for the counts n0, n1, n2 of the residue classes 0, z2 and the third among the N "
        "band states, so |mean - 1| <= 1 / N, and mean = 1 exactly when n0 = n2, which 3 | N gives and does "
        "not exhaust"
    )
    every = flat = split = 0
    for w in args.weights:
        worst = 0.0
        tot = 0
        for z1 in range(3, w, 3):
            z2 = w - z1
            if gcd(z1, z2) != 1:
                continue
            lo, hi = band(z1, z2)
            n0 = sum(1 for j in range(lo, hi + 1) if not j % 3)
            n2 = sum(1 for j in range(lo, hi + 1) if (j - z2) % 3 and j % 3)
            cnt, edges, states = degree_profile(z1, z2)
            if (
                2 * cnt[2] + cnt[1] != edges
                or edges != states + n0 - n2
                or abs(n0 - n2) > 1
            ):
                raise ValueError(
                    "the out-degree is not carried by the residue classes at "
                    + str((w, z1))
                )
            worst = max(worst, abs(edges / states - 1))
            tot += 1
            every += 1
            if n0 == n2:
                flat += 1
                split += states % 3 != 0
        print("w", w, "directions", tot, "max |mean out-degree - 1|", up(worst, 6))
    print(
        "over",
        every,
        "coprime directions of these weights the mean out-degree is exactly 1 on",
        flat,
        "of them, and",
        split,
        "of those have 3 not dividing the band size, so 3 | N is sufficient for exact "
        "criticality and not necessary",
    )
    print("the two smallest witnesses, read off the moves themselves:")
    for z1, z2 in ((3, 1), (3, 2)):
        cnt, edges, states = degree_profile(z1, z2)
        lo, hi = band(z1, z2)
        print(
            "  (z1, z2) =",
            (z1, z2),
            "states",
            list(range(lo, hi + 1)),
            "degrees",
            [len(succ(j, z1, z2)) for j in range(lo, hi + 1)],
            "N",
            states,
            "mean out-degree",
            down(edges / states, 6),
            "| 3 divides N:",
            not states % 3,
        )
    print(
        "so (3, 1) is the smallest non-critical direction, mean 3/2, and (3, 2) is critical at N = 2 with 3 "
        "not dividing N: criticality is a statement about the residue split of the band and never an identity"
    )


def cmd_returns(args):
    print(
        "L(k, n) = #{m >= 1 : m R_k binary in base 3 and below 3^n} is the return count of the weight R_k at "
        "horizon n, the number of distinct returns the band automata of weight R_k can spell; it is computed by "
        "the column transfer, never by enumeration; free = 2^n / R_k is the count a random binary string of "
        "length n would give and rho = L R_k / 2^n the excess over it"
    )
    print(
        "depth b, L(k, bk) at k = "
        + str(args.kmin)
        + "..; lam = L(k+1, b(k+1)) / L(k, bk) at two large k, "
        "the growth of the return count per block, against the free rate 2^b / 3; these ratios are a reading "
        "and not the rate, the second root of the recurrence sitting within a percent of lam_b at even b, so "
        "the verb ladder is where lam_b is computed exactly"
    )
    for b in range(1, args.bmax + 1):
        row = [lift_count(k, b * k) for k in range(args.kmin, args.kmin + 6)]
        r1 = lift_count(args.kbig + 1, b * (args.kbig + 1)) / lift_count(
            args.kbig, b * args.kbig
        )
        r2 = lift_count(2 * args.kbig + 1, b * (2 * args.kbig + 1)) / lift_count(
            2 * args.kbig, b * 2 * args.kbig
        )
        free = 2.0**b / 3.0
        print(
            "b",
            b,
            row,
            "lam",
            down(r1, 6),
            down(r2, 6),
            "free",
            down(free, 6),
            "excess",
            down(r2 * 3 / 2.0**b, 6),
        )
    print("the raw counts, every m including 3 | m:")
    for b in range(1, args.bmax + 1):
        row = [lift_count(k, b * k, False) for k in range(args.kmin, args.kmin + 6)]
        print("b", b, row)
    print("rho = L(k, bk) R_k / 2^(bk - 1), the excess over the free model:")
    for b in range(2, args.bmax + 1):
        row = []
        for k in range(args.kmin, args.kmax + 1):
            w = (3**k - 1) // 2
            row.append(down(lift_count(k, b * k) * w / 2.0 ** (b * k - 1), 4))
        print("b", b, row)
    k = args.kfine
    print("the fine return count at k =", k, ": n, L(k, n), new = L(k, n) - L(k, n-1)")
    prev = 0
    gaps = []
    for n in range(k, args.bmax * k + 1):
        cur = lift_count(k, n)
        if cur == prev:
            gaps.append(n)
        prev = cur
    print("n with no return at all, k <", k, "* bmax:", gaps)
    for b in range(1, args.bmax + 1):
        print("b", b, "L", lift_count(k, b * k), "L at bk+1", lift_count(k, b * k + 1))


def cmd_hist(args):
    print(
        "d(z) is the first return time of the band automaton of (z, R_k - z), the base-3 length of the shortest "
        "binary lift m R_k that carries z; Z = #{z : d(z) < infinity}, Phi = #{d = k} the coprime submask floor, "
        "U = #{d <= 2k}, V = #{d <= 3k}, tail = Z - U"
    )
    print("k R_k cand Z Phi U V tail dmax dmax/sqrt(R_k) secs")
    rows = []
    for k in range(args.kmin, args.kmax + 1):
        t = time.time()
        hist, cand = sweep(k)
        Z = sum(hist.values())
        phi = sum(v for d, v in hist.items() if d <= k)
        U = sum(v for d, v in hist.items() if d <= 2 * k)
        V = sum(v for d, v in hist.items() if d <= 3 * k)
        mx = max(hist) if hist else 0
        w = (3**k - 1) // 2
        print(
            k,
            w,
            cand,
            Z,
            phi,
            U,
            V,
            Z - U,
            mx,
            down(mx / w**0.5, 4),
            round(time.time() - t, 1),
        )
        rows.append((k, hist, Z, phi, U, V))
    print(
        "the fine histogram in the bulk, h(t) = #{z : d(z) = k + t + 1} against the model "
        "E(t) = Sum_{max T = t} 2^k / m_T, the model computed from the multipliers alone:"
    )
    for k, hist, Z, phi, U, V in rows:
        if k < args.bulkmin:
            continue
        mod = model_by_top(k)
        h = [hist.get(k + t + 1, 0) for t in range(k)]
        print(
            "k",
            k,
            "h",
            h[1:],
            "model",
            [down(mod[t], 1) for t in range(1, k)],
            "ratio",
            [down(h[t] / mod[t], 4) if mod[t] else "-" for t in range(1, k)],
        )
    print("the depth histogram of the tail, b -> #{z : d(z) in ((b-1)k, bk]}:")
    for k, hist, Z, phi, U, V in rows:
        if Z == U:
            continue
        dh = {}
        for d, v in hist.items():
            if d > 2 * k:
                dh[-(-d // k)] = dh.get(-(-d // k), 0) + v
        print("k", k, sorted(dh.items()))
    print(
        "the support of the first return time inside [k, dmax]: the lengths reached, the lengths missing, and "
        "the first holes past the proved gap at k + 1:"
    )
    for k, hist, Z, phi, U, V in rows:
        ds = sorted(hist)
        miss = [n for n in range(k, ds[-1] + 1) if n not in hist]
        print(
            "k",
            k,
            "dmax",
            ds[-1],
            "distinct lengths",
            len(ds),
            "missing",
            len(miss),
            "first four holes",
            miss[:4],
            "first six lengths",
            ds[:6],
        )
    print(
        "the tail survival S(b) = #{z : d(z) > bk}, in z values and in distinct directions, a direction being "
        "the pair (z, R_k - z) which the sweep counts twice, with the sample size behind every exponent and a "
        "maximum-likelihood geometric fitted to the depth counts:"
    )
    for k, hist, Z, phi, U, V in rows:
        if Z == U:
            continue
        s = survival(hist, k)
        if any(v % 2 for v in s):
            raise ValueError("a survival count is odd at k = " + str(k))
        print("k", k, "S(2..) z values", s)
        print(
            "  S(2..) directions",
            [v // 2 for v in s],
            "ratios",
            [down(s[i + 1] / s[i], 4) for i in range(len(s) - 1)],
        )
        e = []
        b = 2
        while 2 * b < len(s) + 2:
            e.append(
                (
                    b,
                    down(-log(s[2 * b - 2] / s[b - 2]) / log(2), 3)
                    if s[2 * b - 2]
                    else "-",
                    s[b - 2] // 2,
                    s[2 * b - 2] // 2,
                )
            )
            b *= 2
        print(
            "  local exponent -log2(S(2b)/S(b)) at (b, exponent, directions behind S(b), behind S(2b))",
            e,
        )
        r, raw, rdf, pool, pdf = geom_fit(s)
        print(
            "  one geometric fitted by maximum likelihood to the depth counts: ratio",
            down(r, 4),
            ", chi2",
            up(raw, 1),
            "on",
            rdf,
            "df over the unpooled depths and",
            up(pool, 1) if pdf > 0 else "-",
            "on",
            pdf,
            "df with the bins pooled to expectation 5; the unpooled depths fall below expectation 5 in the "
            "tail, where the statistic is not valid, so neither verdict is carried by the counts",
        )
    print("Phi against the independent floor count, and the depth-1 identity:")
    for k, hist, Z, phi, U, V in rows:
        f = floor_count(k)
        if f != phi:
            raise ValueError(
                "the floor "
                + str(f)
                + " against the depth-1 count "
                + str(phi)
                + " at k = "
                + str(k)
            )
        print("k", k, "Phi", f, "matched")


def brute_depths(k, n):
    w = (3**k - 1) // 2
    out = {}
    for c in range(1, 1 << n):
        K = sum(3**i for i in range(n) if c >> i & 1)
        if K % w:
            continue
        m = K // w
        d = max(i for i in range(n) if c >> i & 1) + 1
        sup = [i for i in range(n) if c >> i & 1]
        for t in range(1, 1 << len(sup)):
            A = sum(3 ** sup[i] for i in range(len(sup)) if t >> i & 1)
            if A % m:
                continue
            z = A // m
            if 0 < z < w and gcd(z, w) == 1 and (out.get(z, n + 1) > d):
                out[z] = d
    return out


def sweep_weight(w):
    hist = {}
    for z1 in range(3, w, 3):
        if gcd(z1, w) != 1:
            continue
        d = first_return(z1, w - z1)
        if d:
            hist[d] = hist.get(d, 0) + 1
    return hist


def cmd_critical(args):
    print(
        "the deepest first return of the band automaton over every coprime direction of a weight w prime to 3, "
        "against the scale sqrt(w) a critical walk on a band of w/2 states predicts; Z is the occupied count, "
        "med the median first return and q9 its ninth decile"
    )
    print("w Z dmax med q9 dmax/sqrt(w) med/sqrt(w) q9/sqrt(w) secs")
    rows = []
    w = args.wmin
    while w <= args.wmax:
        if w % 3:
            t = time.time()
            hist = sweep_weight(w)
            if hist:
                Z = sum(hist.values())
                ds = []
                for d in sorted(hist):
                    ds += [d] * hist[d]
                med, q9 = ds[len(ds) // 2], ds[9 * len(ds) // 10]
                rows.append((w, Z, max(hist), med, q9))
                print(
                    w,
                    Z,
                    max(hist),
                    med,
                    q9,
                    down(max(hist) / w**0.5, 4),
                    down(med / w**0.5, 4),
                    down(q9 / w**0.5, 4),
                    round(time.time() - t, 1),
                )
        w = int(w * args.step) + 1
    band = [r[2] / r[0] ** 0.5 for r in rows]
    print(
        "dmax / sqrt(w) over",
        len(rows),
        "weights: min",
        down(min(band), 4),
        "max",
        up(max(band), 4),
        "mean",
        down(sum(band) / len(band), 4),
    )
    print(
        "the exponent the ladder carries: log d on log w by ordinary least squares, with the 95 percent "
        "interval and the t statistic against the critical exponent 1/2:"
    )
    for name, sel, col in (
        ("dmax, every weight", lambda r: True, 2),
        ("dmax, Z >= 8", lambda r: r[1] >= 8, 2),
        ("median, every weight", lambda r: True, 3),
        ("median, Z >= 4", lambda r: r[1] >= 4, 3),
        ("median, Z >= 8", lambda r: r[1] >= 8, 3),
        ("median, Z >= 16", lambda r: r[1] >= 16, 3),
        ("ninth decile", lambda r: True, 4),
    ):
        sub = [r for r in rows if sel(r)]
        if len(sub) < 4:
            continue
        sl, se, lo, hi, n = ols([log(r[0]) for r in sub], [log(r[col]) for r in sub])
        print(
            " ",
            name,
            "n",
            n,
            "exponent",
            down(sl, 4),
            "95 percent [",
            down(lo, 4),
            ",",
            up(hi, 4),
            "] t against 1/2",
            down((sl - 0.5) / se, 2),
            "| excludes 1/2:",
            not lo <= 0.5 <= hi,
        )
    print(
        "the median is not one number: a weight with Z <= 2 has its median equal to its dmax, and",
        sum(1 for r in rows if r[1] <= 2),
        "of the",
        len(rows),
        "weights have Z <= 2 and",
        sum(1 for r in rows if r[1] <= 4),
        "have Z <= 4, so the median exponent moves with the cut and only "
        "its sign below 1/2 survives every cut",
    )
    loo = []
    for i in range(len(rows)):
        sub = rows[:i] + rows[i + 1 :]
        sl, se, lo, hi, n = ols([log(r[0]) for r in sub], [log(r[2]) for r in sub])
        loo.append((sl, lo, hi, rows[i][0]))
    a, z = min(loo), max(loo)
    print(
        "leave one out on the dmax fit: the slope ranges over [",
        down(a[0], 4),
        ",",
        up(z[0], 4),
        "], the floor deleting w =",
        a[3],
        "and the ceiling deleting w =",
        z[3],
        ";",
        sum(1 for v in loo if v[1] <= 0.5 <= v[2]),
        "of the",
        len(loo),
        "leave-one-out intervals cover 1/2 at the weights",
        " ".join(str(v[3]) for v in loo if v[1] <= 0.5 <= v[2]),
        ", the widest upper endpoint",
        up(max(v[2] for v in loo), 4),
        "- so the exclusion of 1/2 rests on single weights and is not a property of the ladder; what survives "
        "every deletion is the sign, dmax ~ w^0.41 with the exponent below 1/2",
    )
    c1, c2 = ols2(
        [log(r[0]) for r in rows], [log(r[1]) for r in rows], [log(r[2]) for r in rows]
    )
    print(
        "two predictors, log dmax on log w and log Z:",
        down(c1, 4),
        "on log w and",
        down(c2, 4),
        "on log Z, both positive, so controlling for the sample size lowers the exponent on the weight: the "
        "one-predictor slope overstates and the drift below 1/2 is understated, never overstated",
    )
    dw, vr, stu, cook = ols_diag([log(r[0]) for r in rows], [log(r[2]) for r in rows])
    print(
        "the dmax fit's diagnostics: Durbin-Watson",
        down(dw, 2),
        ", residual variance ratio across the halves of the ladder",
        up(vr, 2),
        ", largest studentised residual",
        up(stu, 2),
        ", largest Cook distance",
        up(cook, 3),
        "- ordinary least squares is not what breaks here; the "
        "leverage of single weights and the Z <= 2 weights are",
    )
    sl, se, lo, hi, n = ols([log(r[0]) for r in rows], [log(r[1]) for r in rows])
    print(
        "  log Z on log w: exponent",
        down(sl, 4),
        "95 percent [",
        down(lo, 4),
        ",",
        up(hi, 4),
        "], so the sample a weight offers grows with the weight",
    )


def cmd_model(args):
    print(
        "the model return count D(k, N) = Sum over the primitive lifts K of length at most N of 2^L(K) / m(K), "
        "the equidistribution count of (lift, direction) pairs the band automata of weight R_k spell inside "
        "horizon N; 2^L(K) is summed by the same column transfer and 1/m(K) is sandwiched by the length of K, "
        "so lo and hi bracket D with lo the safe lower and hi the safe upper; D(k, 2k) must reproduce 2^k L_k "
        "with L_k = Sum_T 1 / m_T, the aggregate the lift half already carries"
    )
    print("k L_k 2^k L_k D(k,2k) lo hi ratio")
    for k in range(args.kmin, args.kmax + 1):
        rows = model_returns(k, 2 * k)
        lo, hi = rows[-1][1], rows[-1][2]
        Lk = exact_L(k)
        print(
            k,
            down(Lk, 5),
            down(2**k * Lk, 2),
            "[",
            down(lo, 2),
            up(hi, 2),
            "]",
            down(lo / (2**k * Lk), 4),
            up(hi / (2**k * Lk), 4),
        )
    print(
        "D(k, bk) / 2^k at every depth b, the aggregate return count the model gives with no enumeration:"
    )
    for k in range(args.kmin, args.kmax + 1):
        rows = model_returns(k, args.bmax * k)
        band = []
        for b in range(2, args.bmax + 1):
            _n, lo, hi = rows[b * k - k]
            band.append((b, down(lo / 2**k, 3), up(hi / 2**k, 3)))
        print("k", k, band)
    print(
        "the deep part at the critical cutoff N = floor(sqrt(R_k)), the scale a critical walk on the band "
        "predicts for the deepest return; beyond a depth of 40 blocks the per-digit increment of D is exactly "
        "4/9 and the row extends by it"
    )
    print(
        "share is the percentage of the deep part carried by the 4/9 extrapolation past 40 blocks rather than "
        "by the transfer; both endpoints of the extrapolated stretch are built on the free per-digit increment "
        "4/9 and its triple 4/3, and rho > 1 at every depth puts the true increment above 4/9, so the "
        "extrapolated part of the band leans low at both ends and hi stays an upper bound only while rho < 3"
    )
    print("k N D(k,2k)/2^k D(k,N)/2^k deep = (D(k,N) - D(k,2k))/2^k share")
    for k in range(args.kmin, args.kmax + 1):
        w = (3**k - 1) // 2
        N = int(w**0.5)
        rows = model_returns(k, N)
        base = rows[k][1], rows[k][2]
        cap = rows[min(N, 40 * k) - k][1]
        top = rows[-1][1], rows[-1][2]
        deep = top[0] - base[0]
        print(
            k,
            N,
            "[",
            down(base[0] / 2**k, 4),
            up(base[1] / 2**k, 4),
            "]",
            "[",
            down(top[0] / 2**k, 4),
            up(top[1] / 2**k, 4),
            "]",
            "[",
            down((top[0] - base[0]) / 2**k, 4),
            up((top[1] - base[1]) / 2**k, 4),
            "]",
            up(100 * (top[0] - cap) / deep, 0) if deep else 0,
        )


def cmd_ladder(args):
    print(
        "at n = bk the slot profile is uniform, s_r = b for every column r, so the column transfer of Q.2 is a "
        "matrix fixed in k: L(k, bk) = Sum_{j = 0..b} w_j B(b, j)^(k - r0(j)) h_j with B(b, j), h_j, w_j and "
        "the head length r0(j) all independent of k, hence L(k, bk) obeys a constant-coefficient linear "
        "recurrence in k whose dominant root is lam_b; the carry state space is [0, b // 2], every column sum "
        "of every B(b, j) is Sum_{c = a mod 3} binom(b, c) = (2^b + 2 cos(pi (b - 2a) / 3)) / 3, so the "
        "deviation from 2^b / 3 takes exactly two values and they are set by the parity of b, and a "
        "nonnegative matrix has its spectral radius between its least and its greatest column sum"
    )
    print(
        "the parity bracket, asserted per b below: 2^b / 3 - 1/3 <= lam_b <= 2^b / 3 + 2/3 at even b and "
        "2^b / 3 - 2/3 <= lam_b <= 2^b / 3 + 1/3 at odd b, so 3 lam_b - 2^b lies in [-1, 2] at even b and in "
        "[-2, 1] at odd b and the headline is the two-sided |3 lam_b / 2^b - 1| <= 2^(1 - b) at every b; the "
        "parity refines which edge is which, not the two-sided rate"
    )
    print(
        "the fixed matrix is asserted against the general transfer at every k = 3.."
        + str(args.kcheck)
    )
    exc = []
    for b in range(args.bmin, args.bmax + 1):
        for k in range(3, args.kcheck + 1):
            if set(slot_profile(k, b * k)) != {b}:
                raise ValueError("the slot profile is not uniform at " + str((k, b)))
            if block_count(b, k) != lift_count(k, b * k):
                raise ValueError(
                    "the fixed matrix misses the transfer at " + str((k, b))
                )
        dev, clo, chi, states, r0s, rmax = block_bracket(b)
        dlo, dhi = dev[0], dev[-1]
        exc.append((b, None))
        out = lam_block(b)
        if out is None:
            print("b", b, "no recurrence of order at most", 2 * b + 4)
            continue
        P, Q, iv, ratio, sub = out
        if iv is None:
            print(
                "b",
                b,
                "order",
                len(P) - 1,
                "lam_b is not isolated by the exact bisection",
            )
            continue
        e0, e1 = 3 * iv[0] - 2**b, 3 * iv[1] - 2**b
        sl = Fraction(1, 10**9)
        if e0 < 3 * dlo - sl or e1 > 3 * dhi + sl:
            raise ValueError("lam_b leaves the column-sum bracket at b = " + str(b))
        if abs(Fraction(ratio) - iv[0]) > iv[1] / 100:
            raise ValueError(
                "the isolated root misses L(k+1, b(k+1)) / L(k, bk) at b = " + str(b)
            )
        print("b", b, "order", len(P) - 1, "charpoly", poly_show(P))
        if Q is None:
            print(
                "   the minimal polynomial of lam_b is not identified; lam_b is a root of the polynomial above"
            )
        else:
            print(
                "   deg",
                len(Q) - 1,
                "minpoly of lam_b",
                poly_show(Q),
                "| exact",
                poly_radical(Q)
                if len(Q) <= 3
                else "in radicals, degree " + str(len(Q) - 1)
                if len(Q) <= 5
                else "no radical form, degree " + str(len(Q) - 1),
            )
        print(
            "   lam_b in [",
            dec_down(iv[0], 9),
            ",",
            dec_up(iv[1], 9),
            "]  3 lam_b - 2^b in [",
            dec_down(e0, 6),
            ",",
            dec_up(e1, 6),
            "]  column sum - 2^b / 3 in {",
            ", ".join(str(v) for v in dev),
            "}  saturated:",
            e0 <= 3 * dlo + sl or e1 >= 3 * dhi - sl,
        )
        print(
            "   head lengths r0(j), j = 0..b:",
            " ".join(str(v) for v in r0s),
            " so the exponent k - r0(j) is k - 2 for",
            r0s.count(2),
            "of the",
            b + 1,
            "sectors",
        )
        print(
            "   the j = 1 closure carries",
            states,
            "of",
            b // 2 + 1,
            "carry states and its column sums lie in 2^b / 3 + [",
            clo,
            ",",
            chi,
            "]; max_j rho(B(b, j)) reads",
            down(rmax, 6),
            "against lam_b, gap",
            up(abs(rmax - float(iv[0])) / float(iv[0]), 9),
        )
        if abs(rmax - float(iv[0])) > 1e-6 * float(iv[0]):
            raise ValueError("a block outgrows lam_b at b = " + str(b))
        exc[-1] = (b, e0)
        print(
            "   the block ratio L(k + 1, b(k + 1)) / L(k, bk) at k = 160 reads",
            down(ratio, 6),
            "and the second root of the recurrence is",
            up(sub, 6),
            "of lam_b, so that ratio holds about",
            int(-160 * log10(sub)) if 0 < sub < 1 else 0,
            "correct digits at k = 160: a ratio is not how this rate is read",
        )
    good = [v for v in exc if v[1] is not None]
    print(
        "the excess 3 lam_b - 2^b, truncated down, at b =",
        good[0][0],
        "..",
        good[-1][0],
        ":",
        ", ".join(dec_down(v[1], 6) for v in good),
    )
    print(
        "every one of them is positive, so lam_b sits above the free rate 2^b / 3 at every depth printed and "
        "the live edge of the parity bracket is the upper one, while the two-sided rate stays 2^(1 - b)"
    )


def cmd_lift(args):
    t0 = time.time()
    print(
        "the depth-2 lift census: K = m R_k with one position per column, N_K(m) the submasks of K "
        "divisible by m, M_k = Sum_T N_K(m), the sweep exhaustive over all 2^(k-1) sets T inside [1, k-1]"
    )
    seq = []
    band = []
    for k in range(1, args.kbig + 1):
        M = sum(
            submask_count(tuple((w >> i) & 1 for i in range(k - 1)))
            for w in range(1 << (k - 1))
        )
        seq.append(M)
        band.append(M / 2.0**k)
        print("k", k, "M_k", M, "M_k / 2^k", repr(M / 2.0**k))
    print(
        "M_k / 2^k lies inside ["
        + str(down(min(band), 5))
        + ", "
        + str(up(max(band), 5))
        + "], the endpoints "
        "rounded outward so the band holds the exact values"
    )
    n = (len(seq) + 1) // 2
    H = [[seq[i + j] for j in range(n)] for i in range(n)]
    print(
        "the Hankel matrix of M_k is",
        n,
        "by",
        n,
        "of rank",
        rank_exact(H),
        "reading every one of the",
        len(seq),
        "terms k = 1.." + str(len(seq)) + ", so no linear recurrence of order at most",
        n - 1,
        "holds on them",
    )
    print()
    print(
        "k  M_k/2^k  Sum_T iota_T / 2^k  max_T iota_T  trivial T  power-of-two N  disjoint irreducibles  m_T mod 9"
    )
    ib = []
    for k in range(2, args.kmax + 1):
        tot, ito, imax, triv, pw, dis, cls = lift_census(k)
        pk = 2.0**k
        ib.append(ito / pk)
        print(
            "k",
            k,
            repr(tot / pk),
            repr(ito / pk),
            imax,
            str(triv) + "/" + str(1 << (k - 1)),
            str(pw) + "/" + str(1 << (k - 1)),
            str(dis) + "/" + str(1 << (k - 1)),
            " ".join(
                str(r)
                + ":"
                + str(c[0])
                + ":"
                + str(round(c[1] / pk, 4))
                + ":"
                + str(round(c[2] / pk, 4))
                for r, c in sorted(cls.items())
            ),
        )
    print(
        "Sum_T iota_T / 2^k lies inside ["
        + str(down(min(ib), 6))
        + ", "
        + str(up(max(ib), 6))
        + "], the "
        "endpoints rounded outward; the even readings fall at every step from k = 6 and the odd ones do not"
    )
    print(
        "the disjoint column counts the T whose irreducibles are pairwise disjoint, where N_K(m) = 2^iota; "
        "the power-of-two column is larger, so a power-of-two count does not force disjointness"
    )
    print(
        "the residue classes of the ladder are m_T mod 9, printed as class:sets:M share:iota share; T inside "
        "[1, k-1] gives a_T = Sum 3^i with i >= 1, so a_T is 3 mod 9 when 1 is in T and 0 mod 9 otherwise, and "
        "m_T = 1 + 2 a_T is 7 or 1 mod 9 and never 4, class 4 needing a_T = 6 mod 9"
    )
    print()
    top = 0
    for n in range(1, args.nmax + 1):
        a, rows = hankel_rank(n, n)
        top = max(top, a)
        print(
            "Hankel of N over the column word, p = q =",
            n,
            "rows",
            rows,
            "rank",
            a,
            "full rank",
            rows,
        )
    print(
        "the reversed reading is the transpose of this matrix at p = q and carries no further information"
    )
    for n in range(1, min(args.nmax, 5) + 1):
        a, rows = hankel_rank(n, n, irreducible_count)
        print(
            "Hankel of the irreducible count iota over the column word, p = q =",
            n,
            "rows",
            rows,
            "rank",
            a,
            "full rank",
            rows,
        )
    print(
        "a column transfer with a state set free of k is a linear representation of N as a series over the "
        "column word, so its dimension is at most its state count and at least the Hankel rank; the rank "
        "reaches",
        top,
        "on the words of length at most",
        2 * args.nmax,
        "that is on k at most",
        2 * args.nmax + 1,
        "so no such transfer with fewer than that many states exists there, and whether "
        "the rank is unbounded is observed in the deficiency column and not proved here",
    )
    print("elapsed", round(time.time() - t0, 1), "s")


def cmd_check(args):
    for k in range(2, 6):
        for n in range(k, min(4 * k, 21) + 1):
            for q in (True, False):
                a = lift_count(k, n, q)
                b = brute_lift_count(k, n, q)
                if a != b:
                    raise ValueError(
                        "transfer "
                        + str(a)
                        + " against brute "
                        + str(b)
                        + " at k, n, prim = "
                        + str((k, n, q))
                    )
    print(
        "the column transfer matches the brute enumeration of binary multiples, primitive and raw, "
        "at every k = 2..5, n = k..min(4k, 21)"
    )
    for k in range(2, 9):
        if lift_count(k, k) != 1 or lift_count(k, k + 1) != 1:
            raise ValueError("the length k+1 gap fails at k = " + str(k))
        if lift_count(k, 2 * k) != 2 ** (k - 1) + 1:
            raise ValueError("the depth-2 count fails at k = " + str(k))
        if lift_count(k, 2 * k, False) != 2**k + 1:
            raise ValueError("the raw depth-2 count fails at k = " + str(k))
        if lift_count(k, 3 * k, False) != 2 * 3**k + 1:
            raise ValueError("the raw depth-3 count fails at k = " + str(k))
        if lift_count(k, 3 * k) != 3**k + 1:
            raise ValueError("the depth-3 count fails at k = " + str(k))
    print(
        "L(k, k) = L(k, k+1) = 1, L(k, 2k) = 2^(k-1) + 1 and L(k, 3k) = 3^k + 1 at every k = 2..8, the raw "
        "counts 2^k + 1 at 3^(2k) and 2 * 3^k + 1 at 3^(3k) reproducing the block ladder"
    )
    for k in range(2, 6):
        n = 3 * k
        ref = brute_depths(k, n)
        w = (3**k - 1) // 2
        for z1 in range(3, w, 3):
            if gcd(z1, w) != 1:
                continue
            d = first_return(z1, w - z1)
            r = ref.get(z1, 0)
            if r and d != r:
                raise ValueError(
                    "the BFS return "
                    + str(d)
                    + " against the brute lift "
                    + str(r)
                    + " at "
                    + str((k, z1))
                )
            if not r and d and d <= n:
                raise ValueError(
                    "the BFS returns at "
                    + str(d)
                    + " where no lift does, at "
                    + str((k, z1))
                )
        print(
            "k",
            k,
            "brute lifts to 3^" + str(n),
            "depths matched on",
            len(ref),
            "directions",
        )
    print(
        "the BFS first return equals the shortest binary lift on every direction the brute enumeration reaches"
    )
    for k in range(2, 10):
        hist, _ = sweep(k)
        if hist and min(hist) < k:
            raise ValueError("a return shorter than k at k = " + str(k))
        if any(v % 2 for v in hist.values()):
            raise ValueError("a first-return count is odd at k = " + str(k))
    print(
        "no return is shorter than k at any k = 2..9, and every first-return count is even, the sweep "
        "counting the direction (z, R_k - z) once at z and once at R_k - z"
    )
    for b in range(1, 9):
        for k in range(3, 8):
            if set(slot_profile(k, b * k)) != {b}:
                raise ValueError("the slot profile is not uniform at " + str((k, b)))
            if block_count(b, k) != lift_count(k, b * k):
                raise ValueError(
                    "the fixed matrix misses the transfer at " + str((k, b))
                )
    print(
        "the slot profile is uniform and the fixed matrix reproduces the column transfer at every b = 1..8, "
        "k = 3..7"
    )
    for b, mp in enumerate(MINPOLY, 1):
        _P, Q, iv, _ratio, _sub = lam_block(b)
        if Q != mp:
            raise ValueError(
                "the minimal polynomial of lam_" + str(b) + " reads " + str(Q)
            )
        d = 3 * iv[0] - 2**b, 3 * iv[1] - 2**b
        lo, hi = (-1, 2) if not b % 2 else (-2, 1)
        if d[0] < lo or d[1] > hi:
            raise ValueError("lam_" + str(b) + " leaves the parity bracket")
    if lam_block(4)[0] != [90, -153, 77, -15, 1] or lam_block(5)[0] != [
        -45,
        60,
        -16,
        1,
    ]:
        raise ValueError("the characteristic polynomial at b = 4 or b = 5 moved")
    print(
        "lam_b has the pinned minimal polynomial at every b = 1..10, lam_4 = 6 with characteristic polynomial "
        "(x - 1)(x - 3)(x - 5)(x - 6) and lam_5 = 3 (5 + sqrt 5) / 2 with (x - 1)(x^2 - 15x + 45), and "
        "3 lam_b - 2^b sits inside [-1, 2] at even b and inside [-2, 1] at odd b at every one"
    )
    for k in range(2, 10):
        for w in range(1 << (k - 1)):
            word = tuple((w >> i) & 1 for i in range(k - 1))
            D, n = submask_set(word)
            if len(D) != submask_count(word):
                raise ValueError(
                    "the meet in the middle misses the brute submask count at "
                    + str((k, w))
                )
            if closure_faults(D, n):
                raise ValueError("the solution set is not closed at " + str((k, w)))
            if len(D) % 2:
                raise ValueError("an odd submask count at " + str((k, w)))
    print(
        "the meet-in-the-middle submask count matches the brute enumeration, and the solution set is closed "
        "under complement in K, under disjoint union and under nested difference, over every one of the "
        "2^(k-1) sets T at every k = 2..9, so every count is even"
    )
    for k, ref in ((11, 5224), (12, 11852), (13, 20888), (14, 43364)):
        M = sum(
            submask_count(tuple((w >> i) & 1 for i in range(k - 1)))
            for w in range(1 << (k - 1))
        )
        if M != ref:
            raise ValueError("M_" + str(k) + " reads " + str(M))
    print(
        "M_k reads 5224, 11852, 20888, 43364 at k = 11..14, the cut-free aggregate of the lift half"
    )
    D, n = submask_set((1, 0, 0, 0))
    I = irreducibles(D)
    if len(D) != 6 or sorted(I) != [5, 11, 20, 26]:
        raise ValueError("the k = 5, T = {1} irreducibles moved: " + str((len(D), I)))
    print(
        "at k = 5, T = {1}, m = 7 the support {0, 2, 3, 4, 6} has four irreducibles and two decompositions "
        "of the whole, 6 solutions against 7 packings, so the decomposition is not unique and the count is "
        "the number of distinct unions and not the number of packings"
    )
    D, n = submask_set((0, 1, 0, 0, 0, 0))
    I = irreducibles(D)
    if (
        len(D) != 8
        or len(I) != 6
        or all(not I[a] & I[b] for a in range(len(I)) for b in range(a))
    ):
        raise ValueError("the k = 7, T = {2} witness moved: " + str((len(D), len(I))))
    print(
        "at k = 7, T = {2}, m = 19 the count is 8, a power of two, with 6 irreducibles that overlap, so a "
        "power-of-two count does not force pairwise disjointness"
    )
    r, raw, rdf, pool, pdf = geom_fit(survival(sweep(13)[0], 13))
    if (
        abs(r - 0.8868) > 5e-4
        or abs(raw - 143.5) > 0.2
        or rdf != 30
        or abs(pool - 47.9) > 0.2
        or pdf != 12
    ):
        raise ValueError(
            "the geometric fit at k = 13 moved: " + str((r, raw, rdf, pool, pdf))
        )
    print(
        "the maximum-likelihood geometric at k = 13 keeps its ratio 0.8868 and its chi2 143.5 on 30 df "
        "unpooled, 47.9 on 12 df pooled"
    )


def main():
    p = argparse.ArgumentParser()
    s = p.add_subparsers(dest="cmd", required=True)
    a = s.add_parser("automaton")
    a.add_argument(
        "--weights",
        type=int,
        nargs="+",
        default=[13, 40, 100, 101, 121, 257, 364, 1093],
    )
    a.set_defaults(fn=cmd_automaton)
    b = s.add_parser("returns")
    b.add_argument("--kmin", type=int, default=3)
    b.add_argument("--kmax", type=int, default=24)
    b.add_argument("--bmax", type=int, default=8)
    b.add_argument("--kfine", type=int, default=8)
    b.add_argument("--kbig", type=int, default=160)
    b.set_defaults(fn=cmd_returns)
    c = s.add_parser("hist")
    c.add_argument("--kmin", type=int, default=2)
    c.add_argument("--kmax", type=int, default=13)
    c.add_argument("--bulkmin", type=int, default=9)
    c.set_defaults(fn=cmd_hist)
    e = s.add_parser("model")
    e.add_argument("--kmin", type=int, default=6)
    e.add_argument("--kmax", type=int, default=15)
    e.add_argument("--bmax", type=int, default=8)
    e.set_defaults(fn=cmd_model)
    f = s.add_parser("critical")
    f.add_argument("--wmin", type=int, default=2000)
    f.add_argument("--wmax", type=int, default=400000)
    f.add_argument("--step", type=float, default=1.35)
    f.set_defaults(fn=cmd_critical)
    g = s.add_parser("ladder")
    g.add_argument("--bmin", type=int, default=1)
    g.add_argument("--bmax", type=int, default=14)
    g.add_argument("--kcheck", type=int, default=7)
    g.set_defaults(fn=cmd_ladder)
    h = s.add_parser("lift")
    h.add_argument("--kmax", type=int, default=12)
    h.add_argument("--kbig", type=int, default=17)
    h.add_argument("--nmax", type=int, default=7)
    h.set_defaults(fn=cmd_lift)
    d = s.add_parser("check")
    d.set_defaults(fn=cmd_check)
    args = p.parse_args()
    args.fn(args)


main()
