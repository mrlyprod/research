import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction

import numpy as np
from mpmath import iv

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", ".."))
DATA_DIR = os.path.join(REPO, "data", os.path.relpath(HERE, REPO))
sys.path.insert(
    0,
    os.path.join(
        HERE, "..", "..", "..", "..", "mobius", "lab", "py", "mobius-dissection"
    ),
)
import dissection as MD

sys.path.insert(0, os.path.join(HERE, "..", "digit-uniform-bound"))
import ubound as UB

C1 = 2.0 / math.pi
G0 = 0.9625229
C0 = 0.97

# THE CHAIN AGAINST 1/5


def root_z(q, m=1):
    L = math.log(q)
    f = lambda z: (
        (z - m) * (z - 1) ** 2
        - m * (C1 * L * z + G0 * (z - 1) + C1 * (z - 1) ** 2 / (q * z - 1))
    )
    lo, hi = float(m), float(m) + 10.0
    while f(hi) < 0.0:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0.0 else (lo, mid)
    return hi


def clears(q, m=1, e=0.2):
    return root_z(q, m) < q**e * (1 - m / q)


def margin(q, m=1):
    c1 = 2 / iv.pi
    w = iv.mpf(q) ** (iv.mpf(1) / 5) * (1 - iv.mpf(m) / q)
    return (w - m) * (w - 1) ** 2 - iv.mpf(m) * (
        c1 * iv.log(q) * w + iv.mpf(G0) * (w - 1) + c1 * (w - 1) ** 2 / (q * w - 1)
    )


def cap_gap(q):
    c1 = 2 / iv.pi
    return (
        iv.mpf(q) ** (iv.mpf(1) / 5) * (1 - iv.mpf(1) / q)
        - 1
        - iv.sqrt(2 * c1 * iv.log(q) + iv.mpf(C0))
    )


def slope_gap(q):
    c1 = 2 / iv.pi
    return (
        iv.mpf(q) ** (iv.mpf(1) / 5)
        * (1 - iv.mpf(1) / q)
        * iv.sqrt(2 * c1 * iv.log(q) + iv.mpf(C0))
        - 5 * c1
    )


def maynard_alpha(q):
    return math.log((q / (q - 1)) * math.log(q) + 3 * q / (q - 1)) / math.log(q)


def maynard_s(q):
    s = 0
    while (
        math.log((1 + (2 + s + 1) / math.log(q)) * (q / (q - s - 1)) * math.log(q))
        / math.log(q)
        < 0.2
    ):
        s += 1
    return s


def w_budget(q):
    m = 0
    while MD.pb_step3(q, m + 1) < (q - m - 1) * q**-0.8:
        m += 1
    return m


def chain_budget(q):
    m = 0
    while clears(q, m + 1):
        m += 1
    return m


def wall():
    t0 = time.time()
    iv.prec = 120
    gam = (2 / iv.pi) * (iv.euler + iv.log(8 / iv.pi))
    assert float(gam.b) <= G0
    print(
        "THE DIGIT-UNIFORM CHAIN AGAINST 1/5 at one excluded digit: z < base^(1/5)(1 - 1/base) gives alpha_1 < 1/5 at every digit"
    )
    qu = next(q for q in range(3, 10**5) if all(clears(r) for r in range(q, q + 3000)))
    qt = next(q for q in range(qu, 10**5) if float(cap_gap(q).a) > 0.0)
    worst = min((float(margin(q).a), q) for q in range(qu, qt + 1))
    assert worst[0] > 0.0 and float(margin(qu - 1).b) < 0.0
    assert float(cap_gap(qt).a) > 0.0 and float(slope_gap(100).a) > 0.0
    z = root_z(qu)
    a1 = math.ceil(math.log(z * qu / (qu - 1)) / math.log(qu) * 1e6) / 1e6
    print(
        "  certified at 120 bits on [%d, %d]: tightest margin %.4e at base %d; the chain fails at base %d, margin %.4e"
        % (qu, qt, worst[0], worst[1], qu - 1, float(margin(qu - 1).b))
    )
    print(
        "  from base %d the cap 1 + sqrt(2 (2/pi) log base + 0.97) sits below base^(1/5)(1 - 1/base), gap %.4e there,"
        " and the gap grows from base 100 on (slope test %.4f > 0)"
        % (qt, float(cap_gap(qt).a), float(slope_gap(100).a))
    )
    gap = 0.2 - math.log(z * qu / (qu - 1)) / math.log(qu)
    e = math.floor(math.log10(gap))
    gap = math.floor(gap / 10**e * 100) / 100
    print(
        "  so the wall is %d: z = %.6f against %.6f and alpha_1 < %.6f at base %d, 1/5 - alpha_1 >= %.2fe%03d"
        % (qu, z, qu**0.2 * (1 - 1 / qu), a1, qu, gap, e)
    )
    print()
    print(
        "MAYNARD'S WRITTEN CONSTANT, alpha_q <= log((q/(q-1)) log q + 3q/(q-1))/log q, his Section 8"
    )
    qm = next(
        q
        for q in range(2, 10**7)
        if maynard_alpha(q) < 0.2
        and all(maynard_alpha(r) < 0.2 for r in (q + 1, 2 * q, 10 * q))
    )
    print(
        "  least q with alpha_q < 1/5: %d (alpha_q %.9f there, %.9f at %d); at q = 2000001 alpha_q = %.6f < 0.198"
        % (qm, maynard_alpha(qm), maynard_alpha(qm - 1), qm - 1, maynard_alpha(2000001))
    )
    assert maynard_alpha(2000001) < 0.198
    print()
    print(
        "THE MISSING-DIGIT BUDGET at base 10^7 and 10^8: Maynard's C_(q,s) = 1 + (2+s)/log q, the wall condition (W), the chain"
    )
    for q in (10**7, 10**8):
        print(
            "  base %d: Maynard s <= %d, (W) m <= %d, the chain m <= %d"
            % (q, maynard_s(q), w_budget(q), chain_budget(q))
        )
    assert maynard_s(10**8) >= 10
    print(
        "  the chain's own wall at two and three excluded digits: %d and %d"
        % tuple(
            next(
                q
                for q in range(3, 10**5)
                if all(clears(r, m) for r in range(q, q + 3000))
            )
            for m in (2, 3)
        )
    )
    print("runtime %.1f s" % (time.time() - t0))


def window():
    t0 = time.time()
    print(
        "THE DIGIT-UNIFORM WINDOW AGAINST 1/5 at two window digits, one certificate per base for every excluded digit at once"
    )
    rows = {}
    q = 583
    while True:
        rows[q] = UB.uniform_alpha(q, 2)[0]
        if rows[q] >= 0.2:
            break
        q -= 1
    lo = q + 1
    assert all(rows[r] < 0.2 for r in range(lo, 584))
    print(
        "  alpha_1 < 1/5 certified at every base %d..583, the largest bound %.6f at base %d; base %d reads %.6f and fails"
        % (
            lo,
            max(rows[r] for r in range(lo, 584)),
            max(range(lo, 584), key=lambda r: rows[r]),
            q,
            rows[q],
        )
    )
    print(
        "  with the chain from 584 the one-missing-digit wall of the dissection is %d"
        % lo
    )
    print("runtime %.1f s" % (time.time() - t0))


# THE WALL AT EVERY NUMBER OF EXCLUDED DIGITS


def cap_gap_m(q, m):
    c1 = 2 / iv.pi
    return (
        q ** (iv.mpf(1) / 5) * (1 - iv.mpf(m) / q)
        - m
        - iv.sqrt(m * (c1 * iv.log(q) + iv.mpf(C0)))
    )


def slope_m(q, m):
    c1 = 2 / iv.pi
    return (
        q ** (iv.mpf(1) / 5) * iv.sqrt(c1 * iv.log(q) + iv.mpf(C0))
        - iv.mpf(5) / 2 * iv.sqrt(iv.mpf(m)) * c1
    )


def w_exact(q, m):
    n = -(-(q - 2) // 2)
    x = iv.mpf(q)
    phq = (
        4 / iv.pi
        + (2 / iv.pi) * (iv.log(n) + iv.euler + iv.mpf(1) / (2 * n))
        + (1 - 2 / iv.pi) * (x - 2) / x
        + iv.mpf("0.727") / x
    )
    return (x - m) * x ** (-iv.mpf(4) / 5) - iv.sqrt(iv.mpf(m)) - phq


def w_smooth(q, m):
    x = iv.mpf(q)
    s = (
        4 / iv.pi
        + (2 / iv.pi) * (iv.log(x / 2) + iv.euler + 1 / (x - 2))
        + (1 - 2 / iv.pi)
        + iv.mpf("0.727") / x
    )
    return (x - m) * x ** (-iv.mpf(4) / 5) - iv.sqrt(iv.mpf(m)) - s


def certify(fn, lo, hi):
    stack, low = [(lo, hi)], None
    while stack:
        a, b = stack.pop()
        if b - a <= 4:
            for q in range(a, b + 1):
                v = float(fn(iv.mpf(q)).a)
                if v <= 0.0:
                    return None
                if low is None or v < low[0]:
                    low = (v, q)
            continue
        if float(fn(iv.mpf([a, b])).a) > 0.0:
            continue
        mid = (a + b) // 2
        stack += [(a, mid), (mid + 1, b)]
    return low


def least(pred, lo):
    hi = lo
    while not pred(hi):
        hi *= 2
    while not all(pred(q) for q in range(hi, hi + 50)):
        hi *= 2
    lo = hi // 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if all(pred(q) for q in range(mid, mid + 50)) else (mid, hi)
    return hi


def chain_wall(m):
    q0 = least(lambda q: clears(q, m), 86)
    while float(margin(q0 - 1, m).a) > 0.0:
        q0 -= 1
    assert float(margin(q0 - 1, m).b) < 0.0
    if m == 1:
        qt = next(q for q in range(q0, 10**5) if float(cap_gap(q).a) > 0.0)
        assert float(slope_gap(100).a) > 0.0
    else:
        qt = least(lambda q: float(cap_gap_m(iv.mpf(q), m).a) > 0.0, q0)
        assert (
            float(cap_gap_m(iv.mpf(qt), m).a) > 0.0
            and float(slope_m(iv.mpf(qt), m).a) > 0.0
        )
    low = certify(lambda q: margin(q, m), q0, qt)
    assert low is not None
    z = root_z(q0, m)
    gap = 0.2 - math.log(z * q0 / (q0 - m)) / math.log(q0)
    below = (q0 - 1) ** 0.2 * (1 - m / (q0 - 1)) - root_z(q0 - 1, m)
    return q0, qt, low, gap, below


def w_wall(m):
    qs = least(lambda q: float(w_smooth(q, m).a) > 0.0, 327)
    assert qs >= 327 and float(w_smooth(qs, m).a) > 0.0
    q = qs - 1
    while float(w_exact(q, m).a) > 0.0:
        q -= 1
    assert float(w_exact(q, m).b) < 0.0
    at = w_exact(q + 1, m)
    wq = (q + 1 - m) * (q + 1) ** -0.8
    return (
        q + 1,
        qs,
        float(at.a),
        math.log(wq / (wq - float(at.a))) / math.log(q + 1),
        float(w_exact(q, m).b),
    )


def walls():
    t0 = time.time()
    iv.prec = 120
    assert float(((2 / iv.pi) * (iv.euler + iv.log(8 / iv.pi))).b) <= G0
    print(
        "THE WALL PER NUMBER m OF EXCLUDED DIGITS: the least base from which each certificate proves alpha_1 < 1/5 at every choice of m digits"
    )
    print(
        "chain: (z - m)(z - 1)^2 = m((2/pi) log(base) z + gamma'(z - 1) + (2/pi)(z - 1)^2/(base z - 1)) against z < base^(1/5)(1 - m/base)"
    )
    print(
        "(W): sqrt(m) + Phi_base/base < base^(1/5)(1 - m/base), Phi_base the step 3 constant of mobius"
    )
    print(
        " m   chain wall   certified to   root-equation margin   1/5 - alpha_1   w - z one below     (W) wall   smooth from   w - PB at wall   1/5 - alpha_1   w - PB one below   better"
    )
    rows, thin, miss = [], [], []
    for m in range(1, 14):
        c, qt, low, gap, below = chain_wall(m)
        w, qs, wat, wgap, wbelow = w_wall(m)
        e = math.floor(math.log10(gap))
        g = math.floor(gap / 10**e * 100) / 100
        best = "chain" if c < w else "(W)"
        rows.append((m, c, w))
        ew = math.floor(math.log10(wgap))
        gw = math.floor(wgap / 10**ew * 100) / 100
        thin.append((gap, "chain", m, c))
        thin.append((wgap, "(W)", m, w))
        miss.append((-below, "chain", m, c - 1))
        miss.append((-wbelow, "(W)", m, w - 1))
        print(
            "%2d %12d %14d %22.4e %12.2fe%03d %17.3e %12d %13d %16.4e %12.2fe%03d %18.3e   %s"
            % (m, c, qt, low[0], g, e, below, w, qs, wat, gw, ew, wbelow, best)
        )
    cross = next(m for m, c, w in rows if w < c)
    assert all(c < w for m, c, w in rows if m < cross) and all(
        w < c for m, c, w in rows if m >= cross
    )
    print(
        "the chain is the better certificate at m <= %d and (W) at %d <= m <= 13"
        % (cross - 1, cross)
    )
    t, u = min(thin), min(miss)
    et = math.floor(math.log10(t[0]))
    print(
        "the thinnest pass over all 26 walls: 1/5 - alpha_1 >= %.2fe%03d, %s at m = %d, base %d; the closest failure one below: %.3e, %s at m = %d, base %d"
        % (
            math.floor(t[0] / 10**et * 100) / 100,
            et,
            t[1],
            t[2],
            t[3],
            -u[0],
            u[1],
            u[2],
            u[3],
        )
    )
    h = w_smooth(iv.mpf(14) ** 5, 14)
    assert float(h.a) > 0.0
    print(
        "from m = 14 on: the smooth (W) gap at base m^5 is at least %.4f at m = 14 and grows in m, so the (W) wall is below m^5,"
        " while the chain needs z < base^(1/5) with z > m, so its wall is above m^5"
        % (math.floor(float(h.a) * 1e4) / 1e4)
    )
    print("runtime %.1f s" % (time.time() - t0))


# THE SINGULAR SERIES


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            r -= r // p
            while m % p == 0:
                m //= p
        p += 1
    return r - r // m if m > 1 else r


def mob(n):
    k, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            k = -k
        p += 1
    return -k if n > 1 else k


def ramanujan(d, f):
    return sum(mob(d // e) * e for e in range(1, d + 1) if d % e == 0 and f % e == 0)


def kappa(q, F):
    return Fraction(q, phi(q)) * Fraction(
        sum(1 for f in F if math.gcd(f, q) == 1), len(F)
    )


def principal(q, F):
    return sum(
        Fraction(mob(d) * ramanujan(d, f), phi(d))
        for d in range(1, q + 1)
        if q % d == 0
        for f in F
    ) / len(F)


def series():
    t0 = time.time()
    print(
        "THE MAIN TERM of C2: sum over squarefree d | base of mu(d)/phi(d) sum_((l,d)=1) hat F_k(l/d), divided by fill^k,"
    )
    print(
        "against kappa_F = (base/phi(base)) #{f in F : (f, base) = 1}/fill, in exact rationals"
    )
    n = 0
    for q in range(3, 31):
        for m in (1, 2):
            for E in itertools.combinations(range(q), m):
                F = [v for v in range(q) if v not in E]
                assert principal(q, F) == kappa(q, F)
                n += 1
    print("  equal at all %d sets missing one or two digits of every base 3..30" % n)
    for q, E in (
        (10, (5,)),
        (10, (1,)),
        (10, (0, 5)),
        (10, (1, 3)),
        (12, (0, 6)),
        (30, (0, 15, 29)),
    ):
        F = [v for v in range(q) if v not in E]
        s1 = sum(1 for b in E if math.gcd(b, q) == 1)
        v1 = Fraction(q * (phi(q) - s1), (q - 1) * phi(q))
        print(
            "  base %2d missing %-12s kappa_F = %-8s = %.6f; the printed q(phi(q) - s')/((q-1) phi(q)) = %.6f"
            % (q, str(E), str(kappa(q, F)), float(kappa(q, F)), float(v1))
        )
    print("runtime %.1f s" % (time.time() - t0))


# THE REGIONS FOR LAMBDA


def vonmangoldt(N):
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    for p in range(2, math.isqrt(N) + 1):
        if s[p]:
            s[p * p :: p] = False
    ps = np.nonzero(s)[0]
    lam = np.zeros(N + 1)
    lam[ps] = np.log(ps)
    for p in ps[ps <= math.isqrt(N)]:
        pk = int(p) * int(p)
        while pk <= N:
            lam[pk] = math.log(p)
            pk *= int(p)
    return lam


def regions_one(q, e0, k, Z):
    F = [v for v in range(q) if v != e0]
    fill = len(F)
    y = q**k
    D = np.zeros(y)
    D[MD.strings(q, F, k)] = 1.0
    lam = vonmangoldt(y)[:y]
    hatF = np.conj(np.fft.fft(D))
    Sneg = np.fft.fft(lam)
    exact = float(lam[D > 0].sum())
    a = np.arange(y, dtype=np.int64)
    Q = int(y**0.6)
    l, d = MD.convergent(a, y, Q)
    h = np.abs(a * d - l * y)
    sm = np.array([MD.smooth(int(v), q) for v in range(Q + 1)])
    A = d >= y**0.4
    C = (~A) & (d < Z) & (h < Z)
    B = (~A) & (~C)
    C2 = C & sm[d]
    C1 = C & (~sm[d])
    term = hatF * Sneg / y
    parts = {
        nm: term[msk].sum().real
        for nm, msk in (("A", A), ("B", B), ("C1", C1), ("C2", C2))
    }
    tot = sum(parts.values())
    assert abs(tot - exact) < 1e-6 * y
    kap = float(kappa(q, F))
    h0 = C2 & (h == 0)
    mus = np.array([mob(int(v)) for v in d[h0]], dtype=float)
    phs = np.array([phi(int(v)) for v in d[h0]], dtype=float)
    pred = (hatF[h0] * mus / phs).sum().real
    assert abs(pred - kap * fill**k) < 1e-6 * fill**k
    print(
        "base %d, excluded %d, level %d, y = %d, Z = %d, kappa_F = %.6f"
        % (q, e0, k, y, Z, kap)
    )
    print(
        "  exact sum of Lambda over the strings %.4f = %.6f kappa_F fill^k; the four regions sum to it, difference %.2e"
        % (exact, exact / (kap * fill**k), abs(tot - exact))
    )
    print(
        "  the principal characters at the C2 points with h = 0, y mu(d)/phi(d) in place of S: %.6f kappa_F fill^k"
        % (pred / (kap * fill**k))
    )
    print(
        "  C2 read on the grid at h = 0: %.6f kappa_F fill^k"
        % (term[h0].sum().real / (kap * fill**k))
    )
    for nm in ("A", "B", "C1", "C2"):
        msk = {"A": A, "B": B, "C1": C1, "C2": C2}[nm]
        print(
            "  region %-2s points %8d  contribution %+.6f kappa_F fill^k"
            % (nm, msk.sum(), parts[nm] / (kap * fill**k))
        )


def regions():
    t0 = time.time()
    print(
        "THE DISSECTION FOR LAMBDA at small base and level, every grid point a mod y, the regions of mobius"
    )
    print()
    regions_one(10, 5, 6, 16)
    regions_one(10, 1, 6, 16)
    regions_one(5, 2, 9, 12)
    print("runtime %.1f s" % (time.time() - t0))


# THE COUNT

DESIGNS = [
    (10, (5,)),
    (10, (0,)),
    (10, (1,)),
    (10, (9,)),
    (10, (0, 5)),
    (7, (3,)),
    (3, (1,)),
    (5, (1, 3)),
]


def count_upto(x, q, F):
    ds = []
    v = x
    while v:
        ds.append(v % q)
        v //= q
    ds = ds[::-1]
    L, fill, lead = len(ds), len(F), sum(1 for f in F if f > 0)
    tot = sum(lead * fill ** (j - 1) for j in range(1, L))
    for i, g in enumerate(ds):
        tot += sum(1 for f in F if f < g and (i > 0 or f > 0)) * fill ** (L - 1 - i)
        if g not in F:
            return tot
    return tot + 1


def member(v, q, E):
    v = np.asarray(v, dtype=np.int64).copy()
    ok = np.ones(v.shape, dtype=bool)
    while (v > 0).any():
        ok &= ~((v > 0) & np.isin(v % q, E))
        v //= q
    return ok


def census(N):
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    for p in range(2, math.isqrt(N) + 1):
        if s[p]:
            s[p * p :: p] = False
    ps = np.nonzero(s)[0].astype(np.int64)
    del s
    pw, pl = [], []
    for p in ps[ps <= math.isqrt(N)]:
        pk = int(p) * int(p)
        while pk <= N:
            pw.append(pk)
            pl.append(math.log(p))
            pk *= int(p)
    pw = np.array(pw, dtype=np.int64)
    pl = np.array(pl)
    o = np.argsort(pw)
    pw, pl = pw[o], pl[o]
    rng = random.Random(2027)
    xs = sorted(rng.randrange(N // 100, N) for _ in range(12))
    out = []
    for q, E in DESIGNS:
        F = [v for v in range(q) if v not in E]
        mp = member(ps, q, E)
        P = ps[mp]
        cp = np.concatenate([[0.0], np.cumsum(np.log(P))])
        mw = member(pw, q, E)
        W = pw[mw]
        cw = np.concatenate([[0.0], np.cumsum(pl[mw])])
        pts = [q**j for j in range(1, 40) if q**j <= N] + xs
        rows = []
        for x in pts:
            psi = (
                cp[np.searchsorted(P, x, side="right")]
                + cw[np.searchsorted(W, x, side="right")]
            )
            rows.append([x, float(psi), count_upto(x, q, F)])
        out.append({"base": q, "excluded": list(E), "rows": rows})
    return out


def count():
    t0 = time.time()
    N = 10**8
    path = os.path.join(DATA_DIR, "count-%d.json" % N)
    if os.path.exists(path):
        with open(path) as fh:
            data = json.load(fh)
        src = "cached"
    else:
        data = census(N)
        os.makedirs(DATA_DIR, exist_ok=True)
        with open(path, "w") as fh:
            json.dump(data, fh)
        src = "computed"
    for q, E in DESIGNS:
        F = [v for v in range(q) if v not in E]
        cum = np.cumsum(member(np.arange(1, 10**5 + 1), q, E))
        assert all(
            count_upto(x, q, F) == cum[x - 1]
            for x in (1, 9, 10, 11, 99, 100, 12345, 99999, 10**5)
        )
    print(
        "PRIMES ON A DESIGN against the principal-character main term, sum_(n <= x, n in S_F) Lambda(n) / (kappa_F A_F(x)), x <= %d (%s)"
        % (N, src)
    )
    top, seeded = [], []
    for rec in data:
        q, E = rec["base"], tuple(rec["excluded"])
        F = [v for v in range(q) if v not in E]
        kap = float(kappa(q, F))
        pw = [r for r in rec["rows"] if round(math.log(r[0], q), 9).is_integer()]
        rd = [r for r in rec["rows"] if r not in pw]
        tail = (
            ", ".join("%.4f" % (r[1] / (kap * r[2])) for r in pw[-4:])
            if kap
            else "kappa_F = 0"
        )
        rr = [r[1] / (kap * r[2]) for r in rd] if kap else [0.0]
        print(
            "  base %2d missing %-7s kappa_F %.6f  at the last four powers of the base: %s;  at 12 seeded x: %.4f..%.4f"
            % (q, str(E), kap, tail, min(rr), max(rr))
        )
        if kap and any(f + 1 in F for f in F):
            top.append(pw[-1][1] / (kap * pw[-1][2]))
            seeded.extend(rr)
        if not any(f + 1 in F for f in F):
            print(
                "           outside (E): psi_F(%d) = %.4f against kappa_F A_F = %.1f"
                % (pw[-1][0], pw[-1][1], kap * pw[-1][2])
            )
    print(
        "  the sets with two consecutive digits: at the largest power of the base %.4f..%.4f, at the seeded x %.4f..%.4f"
        % (min(top), max(top), min(seeded), max(seeded))
    )
    print("runtime %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "wall"
    {
        "wall": wall,
        "walls": walls,
        "window": window,
        "series": series,
        "regions": regions,
        "count": count,
    }[verb]()
