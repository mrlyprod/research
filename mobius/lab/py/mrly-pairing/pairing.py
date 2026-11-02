import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "..",
        "..",
        "zeta",
        "lab",
        "py",
        "mrly-euler",
    ),
)

import euler

# FAMILIES


def full(q):
    return (1 << q) - 1


def miss(q, d):
    return full(q) & ~(1 << d)


FAMILIES = [(3, 0b011, 14), (4, 0b0011, 11), (5, 0b00011, 9), (10, miss(10, 9), 6)]

CONTROLS = [(2, 0b11, 20), (3, 0b111, 13)]

THETAMAX = {
    (3, 0b011): 0.4869,
    (4, 0b0011): 0.4819,
    (5, 0b00011): 0.4633,
    (10, miss(10, 9)): 0.4871,
}


def digs(q, F):
    return [d for d in range(q) if (F >> d) & 1]


def alpha_of(q, F):
    from math import log

    return log(len(digs(q, F))) / log(q)


def strings(q, F, L):
    import numpy

    out = numpy.zeros(1, dtype=numpy.int64)
    d = numpy.array(digs(q, F), dtype=numpy.int64)
    for i in range(L):
        out = (out[:, None] + d[None, :] * (q**i)).ravel()
    out.sort()
    return out


def indicator(q, F, L):
    import numpy

    v = numpy.zeros(q**L, dtype=numpy.float64)
    v[strings(q, F, L)] = 1.0
    return v


# MOBIUS


def mu_upto(N):
    import numpy

    mu = numpy.ones(N + 1, dtype=numpy.int8)
    mu[0] = 0
    sieve = numpy.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    primes = numpy.nonzero(sieve)[0]
    for p in primes:
        mu[p::p] = -mu[p::p]
    for p in primes:
        if p * p > N:
            break
        mu[p * p :: p * p] = 0
    return mu


def mu_of(vals):
    import numpy

    rem = numpy.array(vals, dtype=numpy.int64)
    mu = numpy.ones(rem.shape, dtype=numpy.int8)
    top = int(rem.max())
    lim = int(top**0.5) + 1
    sieve = numpy.ones(lim + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(lim**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    for p in numpy.nonzero(sieve)[0]:
        d = rem % p == 0
        if not d.any():
            continue
        rem[d] //= p
        mu[d] = -mu[d]
        again = d & (rem % p == 0)
        if again.any():
            mu[again] = 0
            while again.any():
                rem[again] //= p
                again = again & (rem % p == 0)
    mu[rem > 1] = -mu[rem > 1]
    return mu


VERBS = {}

# SPLIT


def levels(q, L):
    import numpy

    lev = numpy.full(q**L, L, dtype=numpy.int16)
    for j in range(L - 1, -1, -1):
        lev[:: q ** (L - j)] = j
    return lev


def split_one(q, F, L):
    from math import log

    import numpy

    N = q**L
    k = len(digs(q, F))
    a = alpha_of(q, F)
    ind = indicator(q, F, L)
    mu = mu_upto(N - 1).astype(numpy.float64)[:N]
    direct = float(numpy.dot(ind, mu))
    Gh = numpy.fft.fft(ind)
    del ind
    Sh = numpy.fft.fft(mu)
    del mu
    paired = float((numpy.conj(Gh) * Sh).sum().real) / N
    zero = float((numpy.conj(Gh[0]) * Sh[0]).real) / N
    e2S = float((numpy.abs(Sh) ** 2).sum()) / N
    del Sh
    absG = numpy.abs(Gh)
    del Gh
    e2G = float((absG**2).sum()) / N
    lev = levels(q, L)
    w = numpy.bincount(lev, weights=absG, minlength=L + 1)
    del absG, lev
    C = float(w.sum())
    lq = log(q)
    ell1 = log(C / N) / (L * lq)
    cand = [
        ((log(w[j]) / lq - L) / L + min(0.5 + j / (2.0 * L), 0.75), j)
        for j in range(L + 1)
        if w[j] > 0
    ]
    best, arg = max(cand)
    return {
        "q": q,
        "F": F,
        "L": L,
        "k": k,
        "alpha": a,
        "direct": direct,
        "paired": paired,
        "e2G": e2G,
        "e2S": e2S,
        "zero": zero,
        "w": w,
        "C": C,
        "ell1": ell1,
        "unif": ell1 + 0.75,
        "best": best,
        "arg": arg,
        "top": float(w[L]) / C,
        "hi": float(w[(L + 1) // 2 :].sum()) / C,
    }


def split():
    import numpy

    print("THE POSITION PAIRING ON THE GRID, EXACT")
    print("M_F(q^L) = q^(-L) sum_a G_L(a/q^L) S_L(a/q^L), both sides finite")
    print(
        "q F        L  M_F(q^L)   pairing      int|G_L|^2 k^L      int|S_L|^2  a=0 term"
    )
    rows = []
    for q, F, L in FAMILIES + CONTROLS:
        r = split_one(q, F, L)
        rows.append(r)
        print(
            "%2d %-8s %2d %10.1f %12.4f %10.0f %8d %11.0f %9.5f"
            % (
                q,
                euler.show(F, q),
                L,
                r["direct"],
                r["paired"],
                r["e2G"],
                r["k"] ** L,
                r["e2S"],
                r["zero"],
            )
        )
    print()
    print("THE LEVEL PROFILE of the l^1 mass, C_L = sum_j k^(L-j) c_j")
    print(
        "share of C_L at level j = L (primitive denominator q^L), at levels j >= L/2,"
    )
    print("and the proved floor m/q for the top level")
    print("q F        L  C_L/q^L      top share  floor m/q  share j >= L/2")
    for r in rows:
        m = r["q"] - r["k"]
        print(
            "%2d %-8s %2d %12.4f %10.6f %10.6f %14.6f"
            % (
                r["q"],
                euler.show(r["F"], r["q"]),
                r["L"],
                r["C"] / r["q"] ** r["L"],
                r["top"],
                m / r["q"],
                r["hi"],
            )
        )
    print()
    print("THE COST, exponents in log_q against the trivial alpha")
    print(
        "cs = Cauchy-Schwarz (alpha+1)/2, unif = l^1 times the uniform GRH max x^(3/4),"
    )
    print("den = the per-denominator split, GRH max x^(1/2) q^(j/2) at level j")
    print(
        "mob = alpha times the thetamax of mobius.md, the measured exponent of M_F itself"
    )
    print(
        "q F        L  alpha    alpha/2  mob      cs       unif     den      C_L/w_j  q/m"
    )
    for r in rows:
        m = r["q"] - r["k"]
        tm = THETAMAX.get((r["q"], r["F"]))
        print(
            "%2d %-8s %2d %8.6f %8.6f %8s %8.6f %8.6f %8.6f %8.4f %8.4f"
            % (
                r["q"],
                euler.show(r["F"], r["q"]),
                r["L"],
                r["alpha"],
                r["alpha"] / 2,
                ("%.6f" % (r["alpha"] * tm)) if tm else "-",
                (r["alpha"] + 1) / 2,
                r["unif"],
                r["best"],
                r["C"] / r["w"][r["arg"]],
                r["q"] / m if m else float("inf"),
            )
        )
    print()
    print("THE LADDER IN L at base 3 F = {0,1}: the per-denominator gain is a constant")
    print("unif - den = log_q(C_L/w_L)/L and C_L/w_L <= q/m = 1.5 at every L")
    print(" L  unif     den      unif-den  L(unif-den) C_L/w_L  top share")
    for L in range(6, 15):
        r = split_one(3, 0b011, L)
        d = r["unif"] - r["best"]
        print(
            "%2d %8.6f %8.6f %9.6f %11.6f %8.4f %10.6f"
            % (L, r["unif"], r["best"], d, L * d, r["C"] / r["w"][r["arg"]], r["top"])
        )
    print()
    print("THE ONE-STEP CONSTANT read off the l^1 norms themselves")
    print(
        "C_L = sum_(a mod q^L) |G_L(a/q^L)|, ratio C_L/C_(L-1) against the grid sup B_q(F)"
    )
    print(
        "q F        L  C_L            C_L/C_(L-1)     sup_t sum_r |g((t+r)/q)|  top share"
    )
    for q, F, L in FAMILIES:
        cs = []
        for l in range(1, min(L, 9) + 1):
            r = split_one(q, F, l)
            cs.append(r["C"])
        d = digs(q, F)
        t = numpy.linspace(0.0, 1.0, 200001)[:-1]
        h = numpy.zeros_like(t)
        for rr in range(q):
            z = numpy.zeros_like(t, dtype=numpy.complex128)
            for dd in d:
                z += numpy.exp(2j * numpy.pi * dd * (t + rr) / q)
            h += numpy.abs(z)
        for l in range(len(cs) - 2, len(cs)):
            print(
                "%2d %-8s %2d %14.6f %15.9f %25.9f %10.6f"
                % (
                    q,
                    euler.show(F, q),
                    l + 1,
                    cs[l],
                    cs[l] / cs[l - 1],
                    h.max(),
                    split_one(q, F, l + 1)["top"],
                )
            )
    print()
    print("THE PRINCIPAL FIBRE carries almost none of the design meter")
    print("the a = 0 term is q^(-L) k^L M(q^L), exponent alpha - 1/2 under RH,")
    print(
        "against the conjectured alpha/2 for M_F itself: alpha - 1/2 < alpha/2 iff alpha < 1"
    )
    print("q F        L  M_F(q^L)   a=0 term   share      alpha-1/2  alpha/2")
    for r in rows:
        sh = r["zero"] / r["direct"] if r["direct"] else float("nan")
        print(
            "%2d %-8s %2d %10.1f %10.5f %10.6f %10.6f %8.6f"
            % (
                r["q"],
                euler.show(r["F"], r["q"]),
                r["L"],
                r["direct"],
                r["zero"],
                sh,
                r["alpha"] - 0.5,
                r["alpha"] / 2,
            )
        )


VERBS["split"] = split

# GLUE

GLUE = [(3, 0b011, 18), (4, 0b0011, 14), (5, 0b00011, 12), (10, miss(10, 9), 7)]

GLUE_CONTROLS = [(2, 0b11, 20), (3, 0b111, 12)]


def coeffs(q, F, L):
    import numpy

    N = q**L
    sf = strings(q, F, L)
    sf = sf[sf >= 1]
    mu = mu_of(sf).astype(numpy.int64)
    c = numpy.zeros(N + 1, dtype=numpy.int32)
    for e, m in zip(sf, mu):
        if m == 0:
            continue
        d = sf[: numpy.searchsorted(sf, N // e, "right")]
        c[d * e] += m
    return c


PHASES = [0.0, 0.25, 0.5, 0.75]


def partials(q, F, L):
    import numpy

    sf = strings(q, F, L)
    sf = sf[sf >= 1]
    mu = mu_of(sf).astype(numpy.int64)
    out = []
    for l in range(1, L):
        row = []
        for c in PHASES:
            x = int(float(q) ** (l + c))
            cut = numpy.searchsorted(sf, x, "right")
            a = numpy.searchsorted(sf, x // sf[:cut], "right")
            row.append((x, int((mu[:cut] * a).sum())))
        out.append(row)
    return out, sf, mu


def glue():

    import numpy

    print("THE COEFFICIENTS OF ZETA_F M_F, c_F(n) = sum_(d e = n, d, e in S_F) mu(e)")
    print("q F        L      n<=q^L    nonzero  first n>1  c_F(n)  max|c_F|  at n")
    for q, F, L in [
        (3, 0b011, 12),
        (4, 0b0011, 9),
        (5, 0b00011, 8),
        (10, miss(10, 9), 6),
    ]:
        c = coeffs(q, F, L)
        nz = numpy.nonzero(c[2:])[0] + 2
        am = int(numpy.abs(c).argmax())
        print(
            "%2d %-8s %2d %10d %10d %10d %7d %9d %5d"
            % (
                q,
                euler.show(F, q),
                L,
                q**L,
                len(nz) + (1 if c[1] else 0),
                int(nz[0]),
                int(c[nz[0]]),
                int(abs(c[am])),
                am,
            )
        )
    print()
    print("THE PARTIAL SUMS P(x) = sum_(n<=x) c_F(n) = sum_(e in S_F) mu(e) A_F(x/e)")
    print("the abscissa of D_F = zeta_F M_F - 1 is read off log|P|/log x against alpha")
    for q, F, L in GLUE + GLUE_CONTROLS:
        a = alpha_of(q, F)
        ps, sf, mu = partials(q, F, L)
        print(
            "q = %d  F = %s  alpha = %.6f  alpha/2 = %.6f"
            % (q, euler.show(F, q), a, a / 2)
        )
        print(
            "  L   P(x) at log_q x = L, L+1/4, L+1/2, L+3/4      P(x)/x^alpha at the same four"
        )
        for l in range(1, L):
            row = ps[l - 1]
            print(
                "  %2d %10d %10d %10d %10d   %9.6f %9.6f %9.6f %9.6f"
                % (
                    (l,)
                    + tuple(v for _, v in row)
                    + tuple(v / float(x) ** a for x, v in row)
                )
            )
    print()
    print(
        "THE LIMIT TEST: sigma_c(D_F) < alpha forces M_F(sigma) -> 0 as sigma -> alpha+"
    )
    print("M_F(sigma) = sum_(n in S_F) mu(n) n^(-sigma), summed to n <= q^L")
    print("q F        L  sigma-alpha  M_F(sigma)   tail bound q^(-L alpha/2)")
    for q, F, L in GLUE:
        a = alpha_of(q, F)
        sf = strings(q, F, L)
        sf = sf[sf >= 1]
        mu = mu_of(sf).astype(numpy.float64)
        x = sf.astype(numpy.float64)
        for eps in [0.2, 0.1, 0.05, 0.02, 0.0]:
            v = float((mu * x ** (-(a + eps))).sum())
            print(
                "%2d %-8s %2d %12.4f %12.6f %18.2e"
                % (q, euler.show(F, q), L, eps, v, float(q) ** (-L * a / 2))
            )


VERBS["glue"] = glue

# INVERSE

SEEDS = {
    (3, 0b011): ("0.720788", "28.60568"),
    (10, miss(10, 9)): ("1.001589", "2.7392"),
}


def dirichlet_inverse(q, F, L):
    import numpy

    N = q**L
    if F == full(q):
        sf = numpy.arange(2, N + 1, dtype=numpy.int64)
    else:
        sf = strings(q, F, L)
        sf = sf[sf >= 2]
        if (F >> 1) & 1 and (L == 0 or F & 1):
            sf = numpy.append(sf, numpy.int64(N))
    nu = numpy.zeros(N + 1, dtype=numpy.int32)
    nu[1] = 1
    lo = 1
    while lo <= N:
        hi = min(2 * lo, N + 1)
        for d in sf[: numpy.searchsorted(sf, hi - 1, "right")]:
            m0 = (lo + d - 1) // d
            m1 = (hi - 1) // d + 1
            if m1 <= m0:
                continue
            nu[m0 * d : (m1 - 1) * d + 1 : d] -= nu[m0:m1]
        lo = hi
    return nu


def refine(q, F):
    from mpmath import mp

    sys.path.insert(
        0,
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "..",
            "..",
            "zeta",
            "lab",
            "py",
            "design-zeta",
        ),
    )
    import design_zeta

    d = design_zeta.Design(q, tuple(digs(q, F)))
    re0, im0 = SEEDS[(q, F)]
    s = mp.mpc(re0, im0)
    for _ in range(8):
        v, _ = d.zeta(s)
        h = mp.mpf(10) ** -8
        v2, _ = d.zeta(s + h)
        s = s - v * h / (v2 - v)
    v, e = d.zeta(s)
    return s, v, e, d.alpha


def series_at(nu, q, L, sigmas):
    import numpy

    acc = [0.0] * len(sigmas)
    out = [[] for _ in sigmas]
    N = q**L
    lo = 1
    marks = [q**l for l in range(1, L + 1)]
    m = 0
    while lo <= N:
        hi = min(lo + (1 << 20), N + 1)
        if m < len(marks) and marks[m] < hi:
            hi = marks[m] + 1
        n = numpy.arange(lo, hi, dtype=numpy.float64)
        v = nu[lo:hi].astype(numpy.float64)
        for i, sg in enumerate(sigmas):
            acc[i] += float((v * n ** (-sg)).sum())
        if m < len(marks) and hi == marks[m] + 1:
            for i in range(len(sigmas)):
                out[i].append(acc[i])
            m += 1
        lo = hi
    return out


def abscissa(q, F, L):
    from mpmath import mp

    sys.path.insert(
        0,
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "..",
            "..",
            "zeta",
            "lab",
            "py",
            "design-zeta",
        ),
    )
    import design_zeta

    d = design_zeta.Design(q, tuple(digs(q, F)))
    nu = dirichlet_inverse(q, F, L)
    a = float(d.alpha)
    r = float(refine(q, F)[0].real)
    sigmas = [round(r + 0.08, 4), round(r + 0.20, 4)]
    out = series_at(nu, q, L, sigmas)
    tgt = [float(1 / d.zeta(mp.mpf(str(sg)))[0]) for sg in sigmas]
    return sigmas, out, tgt, a, r


def ladder_stats(nu, q, L):
    import numpy

    N = q**L
    marks = [q**l for l in range(1, L + 1)]
    sig, mx = [], []
    run = 0
    best = 0
    lo, m = 1, 0
    while lo <= N:
        hi = min(lo + (1 << 20), N + 1)
        if m < len(marks) and marks[m] < hi:
            hi = marks[m] + 1
        c = numpy.cumsum(nu[lo:hi].astype(numpy.int64))
        c += run
        run = int(c[-1])
        numpy.abs(c, out=c)
        best = max(best, int(c.max()))
        if m < len(marks) and hi == marks[m] + 1:
            sig.append(run)
            mx.append(best)
            m += 1
        lo = hi
    return sig, mx


def inverse():
    from math import log

    print("THE DIRICHLET INVERSE OF THE DESIGN INDICATOR, nu_F = 1_(S_F)^(-1)")
    print(
        "zeta_F(s) N_F(s) = 1 exactly, the identity that replaces zeta M = 1 on a design"
    )
    print("control: the full digit set gives nu_F = mu term for term")
    for q in [2, 3]:
        nu = dirichlet_inverse(q, full(q), 17 if q == 2 else 11)
        mu = mu_upto(len(nu) - 1)
        print(
            "  q = %d full, n <= %d, nu_F = mu at every n: %s"
            % (q, len(nu) - 1, bool((nu[1:] == mu[1:]).all()))
        )
    print()
    print(
        "THE ZERO THAT FORCES THE ABSCISSA, refined on the sibling engine lab/py/design-zeta"
    )
    print("sigma_c(N_F) >= Re rho for every zero rho of zeta_F with Re rho > alpha")
    print(
        "q F        zero rho                                   |zeta_F(rho)| bound     alpha"
    )
    from mpmath import mp

    for q, F in [(3, 0b011), (10, miss(10, 9))]:
        r, v, e, a = refine(q, F)
        print(
            "%2d %-8s %-42s %-9s %-9s %s"
            % (
                q,
                euler.show(F, q),
                mp.nstr(r, 16),
                mp.nstr(abs(v), 4),
                mp.nstr(e, 4),
                mp.nstr(a, 10),
            )
        )
    print()
    print("THE DESIGN MERTENS OF nu_F against the design's own mass A_F(x) and against")
    print(
        "the rightmost zero of zeta_F, whose real part is a lower bound for the abscissa"
    )
    for q, F, L in [
        (3, 0b011, 16),
        (4, 0b0011, 12),
        (5, 0b00011, 10),
        (10, miss(10, 9), 7),
    ]:
        a = alpha_of(q, F)
        k = len(digs(q, F))
        nu = dirichlet_inverse(q, F, L)
        raw, mxs = ladder_stats(nu, q, L)
        rz = float(refine(q, F)[0].real) if (q, F) in SEEDS else None
        print(
            "q = %d  F = %s  alpha = %.6f  rightmost censused zero Re = %s"
            % (q, euler.show(F, q), a, ("%.6f" % rz) if rz else "not censused")
        )
        print(
            "  the level ratio of the running maximum against q^alpha = k = %d and" % k
        )
        print("  q^(Re rho) = %s" % ("%.6f" % (q**rz) if rz else "not censused"))
        print("  L   sum nu_F(n)   max|sum|    level ratio  exponent   max/A_F(q^L)")
        prev = 0
        for l in range(1, L + 1):
            x = q**l
            mx = mxs[l - 1]
            ex = log(mx) / log(x) if mx > 1 else 0.0
            print(
                "  %2d %13d %12d %12s %10.6f %14.4f"
                % (
                    l,
                    raw[l - 1],
                    mx,
                    ("%.4f" % (mx / prev)) if prev else "-",
                    ex,
                    mx / float(k**l),
                )
            )
            prev = mx
        del nu

    print()
    print(
        "THE ABSCISSA TEST: partial sums of N_F(sigma) = sum nu_F(n) n^(-sigma) against"
    )
    print("1/zeta_F(sigma) from lab/py/design-zeta, at sigma above and below Re rho;")
    print(
        "the transport theorem, not this table, places sigma_c(N_F) at or above Re rho"
    )
    for q, F, L in [(3, 0b011, 16), (10, miss(10, 9), 7)]:
        sigmas, out, tgt, a, r = abscissa(q, F, L)
        print(
            "q = %d  F = %s  alpha = %.7f  Re rho = %.7f" % (q, euler.show(F, q), a, r)
        )
        print(
            "  sigma    1/zeta_F     partial at q^(L-2)  q^(L-1)      q^L        error"
        )
        for i, sg in enumerate(sigmas):
            p3, p2, p1 = out[i][-3], out[i][-2], out[i][-1]
            print(
                "  %-8.4f %12.6f %14.6f %12.6f %12.6f %11.2e"
                % (sg, tgt[i], p3, p2, p1, abs(p1 - tgt[i]))
            )


VERBS["inverse"] = inverse

# BOX

BOXES = [
    (10, miss(10, 9), "1.00150", "1.00168", "2.73915", "2.73925", 1),
    (3, 0b011, "0.72074", "0.72084", "28.60563", "28.60573", 1),
    (10, miss(10, 9), "0.99900", "1.00050", "2.73810", "2.74030", 0),
]


def box():
    from mpmath import mp

    sys.path.insert(
        0,
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "..",
            "..",
            "zeta",
            "lab",
            "py",
            "design-zeta",
        ),
    )
    import design_zeta

    print(
        "WINDING BOXES FOR THE ZEROS OF zeta_F, argument principle on lab/py/design-zeta"
    )
    print(
        "winding 1 on a rectangle resolves one zero inside it, sampled and not certified,"
    )
    print("pinning Re rho to the box edges; winding 0 resolves no zero there")
    for q, F, x0, x1, y0, y1, want in BOXES:
        d = design_zeta.Design(q, tuple(digs(q, F)))
        f = design_zeta.Cache(d)
        w, mx = design_zeta.box_phase(
            f, mp.mpf(x0), mp.mpf(x1), mp.mpf(y0), mp.mpf(y1), 16
        )
        lo = min(abs(v) for v in f.m.values())
        print(
            "q = %d  F = %s  alpha = %s" % (q, euler.show(F, q), mp.nstr(d.alpha, 10))
        )
        print("  Re in [%s, %s]  Im in [%s, %s]" % (x0, x1, y0, y1))
        print(
            "  winding %s  expected %d  max phase step %s  evaluations %d"
            % (mp.nstr(w, 8), want, mp.nstr(mx, 4), f.n)
        )
        print(
            "  min |zeta_F| on the contour %s  engine bound %s  ratio %s"
            % (mp.nstr(lo, 4), mp.nstr(f.emax, 4), mp.nstr(lo / f.emax, 4))
        )


VERBS["box"] = box


# ONESTEP

GAMMA = 0.5772156649015329

ONESTEP_FAMILIES = [
    (3, 2),
    (10, 9),
    (11, 0),
    (13, 0),
    (100, 0),
    (100, 49),
    (1000, 0),
    (1000, 499),
    (2234, 0),
    (2234, 1116),
    (3690, 0),
    (3690, 1844),
]


def dn(x, d=6):
    from math import floor

    return floor(x * 10.0**d) / 10.0**d


def up(x, d=6):
    from math import ceil

    return ceil(x * 10.0**d) / 10.0**d


def psi_q(q):
    from math import log, pi

    n = -(-(q - 2) // 2)
    return (2 * q / pi) * (log(n) + GAMMA + 1 / (2 * n)) + (1 - 2 / pi) * q


def phi_q(q):
    from math import log, pi

    n = -(-(q - 2) // 2)
    return (
        (4 / pi) * q
        + (2 * q / pi) * (log(n) + GAMMA + 1 / (2 * n))
        + (1 - 2 / pi) * (q - 2)
        + 0.727
    )


def sec_e(q, e0):
    from math import cos, pi

    return 1.0 / cos(pi * (e0 - (q - 1) / 2.0) / q)


def bound_e(q, e0):
    from math import pi

    return (4 / pi) * q + psi_q(q) + q / 2.0 - sec_e(q, e0) / 2.0


def harmonic(n):
    return sum(1.0 / j for j in range(1, int(n) + 1))


def psi2_q(q, exact=False):
    from math import log, pi

    p = q // 2
    if exact:
        h = harmonic(p - 1)
    else:
        h = 0.0 if p < 3 else log(p - 1) + GAMMA + 1 / (2 * (p - 1))
    if q % 2 == 0:
        return (q / pi) * (2 * h - 1 + 1 / p) + (1 - 2 / pi) * q / 2
    return (q / pi) * (2 * h - 1 + 2 / p) + (1 - 2 / pi) * (q / 2 + 1 / (2 * q))


def bound2_e(q, e0):
    from math import pi

    return (4 / pi) * q + psi2_q(q) + q / 2.0 - sec_e(q, e0) / 2.0


def mono_floor(hi=400, exact=False):
    from math import pi

    ok = [psi2_q(q, exact) >= (1 + pi) * q / 2 for q in range(4, hi)]
    for i in range(len(ok)):
        if all(ok[i:]):
            return 4 + i
    return None


def h_argmax_floor(hi=80, n=2001):
    from math import cos, pi

    import numpy

    tau = numpy.linspace(0.0, 0.5, n)
    c = numpy.cos(pi * tau)
    last = 3
    for q in range(4, hi):
        for e0 in range(q):
            b = pi * (e0 - (q - 1) / 2.0) / q
            h = c * (psi2_q(q) - q / 2.0 - numpy.cos(2 * b * tau) / (2 * cos(b)))
            if int(h.argmax()) != 0:
                last = q
    return last + 1, hi - 1


def weight_cost(q):
    from math import pi

    import numpy

    r = numpy.arange(q, dtype=numpy.float64)
    a = 1.0 / numpy.sin(pi * (r + 0.5) / q)
    return 2.0 * float(numpy.sum(1.0 / (a + 1.0)))


def terms(q, e0, t):
    from math import pi, sin

    import numpy

    r = numpy.arange(q, dtype=numpy.float64)
    ph = (t + r) / q
    a = ((-1.0) ** r) * sin(pi * t) / numpy.sin(pi * ph)
    c = e0 - (q - 1) / 2.0
    return numpy.sqrt(a * a + 1.0 - 2.0 * a * numpy.cos(2 * pi * c * ph))


def phase_left(q, e0, t):
    from math import pi, sin

    import numpy

    r = numpy.arange(q, dtype=numpy.float64)
    ph = (t + r) / q
    sgn = numpy.sign(((-1.0) ** r) * sin(pi * t) / numpy.sin(pi * ph))
    c = e0 - (q - 1) / 2.0
    return float((1.0 + sgn * numpy.cos(2 * pi * c * ph)).sum())


def phase_right(q, e0, t):
    from math import cos, pi

    c = e0 - (q - 1) / 2.0
    return q + cos(2 * pi * c * (t - 0.5) / q) / cos(pi * c / q)


def sigma_t(q, e0, t):
    return float(terms(q, e0, t).sum())


def kernel_t(q, t):
    from math import pi, sin

    import numpy

    r = numpy.arange(q, dtype=numpy.float64)
    return float(sin(pi * t) * (1.0 / numpy.sin(pi * (t + r) / q)).sum())


def direct_t(q, e0, t):
    from math import pi

    import numpy

    f = numpy.array([d for d in range(q) if d != e0], dtype=numpy.float64)
    tot = 0.0
    for r in range(q):
        tot += abs(complex(numpy.exp(2j * pi * f * ((t + r) / q)).sum()))
    return tot


ONESTEP_CUT = 2001


def seat(q, e0, n=ONESTEP_CUT):
    import numpy

    best, at = -1.0, 0.0
    for t in numpy.linspace(1e-12, 0.5, n):
        v = sigma_t(q, e0, t)
        if v > best:
            best, at = v, t
    return best, at


def ker_seat(q, n=ONESTEP_CUT):
    import numpy

    best, at = -1.0, 0.0
    for t in numpy.linspace(1e-12, 0.5, n):
        v = kernel_t(q, t)
        if v > best:
            best, at = v, t
    return best, at


def kernel_slack(q, n=ONESTEP_CUT):
    from math import pi, sin

    import numpy

    worst, at = 1e18, 0.0
    for t in numpy.linspace(1e-9, 0.5, n):
        v = (4 / pi) * q + sin(pi * t) * psi2_q(q) - kernel_t(q, t)
        if v < worst:
            worst, at = v, t
    return worst, at


BLOCK = 1 << 20


def psi_grid(q, chord):
    from math import pi

    import numpy

    if not chord:
        n = numpy.ceil((q - 2) / 2)
        return (2 * q / pi) * (numpy.log(n) + GAMMA + 1 / (2 * n)) + (1 - 2 / pi) * q
    p = numpy.floor(q / 2)
    h = numpy.log(p - 1) + GAMMA + 1 / (2 * (p - 1))
    ex = numpy.where(q % 2 == 0, 1.0 / p, 2.0 / p)
    fl = numpy.where(q % 2 == 0, q / 2, q / 2 + 1 / (2 * q))
    return (q / pi) * (2 * h - 1 + ex) + (1 - 2 / pi) * fl


def rhs_grid(q, worst, chord):
    from math import pi

    import numpy

    c = numpy.where(q % 2 == 1, 0.0, 0.5) if worst else (q - 1) / 2
    sec = 1.0 / numpy.cos(pi * c / q)
    return (4 / pi) + psi_grid(q, chord) / q + 0.5 - sec / (2 * q)


def proved_wall(b, worst, hi, chord=False, lo=17):
    import numpy

    first, held = None, 0
    for a in range(lo, hi + 1, BLOCK):
        q = numpy.arange(a, min(a + BLOCK, hi + 1), dtype=numpy.float64)
        ok = numpy.nonzero((q - 1) * q ** (-b) - rhs_grid(q, worst, chord) > 0)[0]
        held += int(ok.size)
        if first is None and ok.size:
            first = int(q[ok[0]])
    return (first, held)


def cost_cross(worst, chord, lo, hi, b=0.75):
    import numpy

    q = numpy.arange(lo, hi, dtype=numpy.float64)
    alpha = numpy.log(q - 1) / numpy.log(q)
    cq = numpy.log(rhs_grid(q, worst, chord)) / numpy.log(q)
    win = (alpha - b - cq) / alpha > 1.0 / (2 * (q - 1) * numpy.log(q))
    i = numpy.nonzero(win)[0]
    if not i.size:
        return None, None
    return int(q[i[0]]), bool(win[i[0] :].all())


def half_max(q):
    from math import pi

    import numpy

    r = numpy.arange(q, dtype=numpy.float64)
    ph = (0.5 + r) / q
    a = ((-1.0) ** r) / numpy.sin(pi * ph)
    c = numpy.arange((q - 1) // 2 + 1, dtype=numpy.float64) - (q - 1) / 2.0
    s = numpy.sqrt(
        a * a + 1.0 - 2.0 * a * numpy.cos(2 * pi * c[:, None] * ph[None, :])
    ).sum(axis=1)
    i = int(numpy.argmax(s))
    return float(s[i]), i


def ceiling_wall(b, worst, hi):
    q = hi
    while q >= 17:
        v, e0 = half_max(q) if worst else (sigma_t(q, 0, 0.5), 0)
        if not ((q - 1) * q ** (-b) > v / q):
            return q + 1, q, e0
        q -= 1
    return 17, None, None


CEIL_HI = 2000

DRAWS = 4000

DRAW_SEED = 1009

DRAW_BASES = [17, 23, 60, 101, 333, 1000, 3690]

MONO2 = 36

RUNGS = [
    ("1/2", 0.75, 3690, 4000000),
    ("13/25", 1417.0 / 1850, 8578, 8000000),
    ("11/20", 913.0 / 1160, 33547, 40000000),
]


def onestep():
    from math import pi

    import numpy

    print("THE EXACT ONE-STEP CONSTANT AT ONE EXCLUDED DIGIT")
    print("B_q(F) = sup_t sum_(r mod q) |g_F((t+r)/q)| with F = {0..q-1} less {e_0}")
    print(
        "reduction: |g_F((t+r)/q)| = |A_r - e(c (t+r)/q)| with A_r = (-1)^r sin(pi t)/sin(pi (t+r)/q)"
    )
    print("and c = e_0 - (q-1)/2, so the sum is exact in O(q) at every t")
    print(" q   e_0    t     reduction        direct sum over F")
    for q, e0, t in [(7, 0, 0.13), (7, 3, 0.5), (11, 5, 0.77), (12, 4, 0.31)]:
        print(
            "%2d %5d %6.2f %16.12f %16.12f"
            % (q, e0, t, sigma_t(q, e0, t), direct_t(q, e0, t))
        )
    print()
    print("THE PHASE IDENTITY the sharpening runs on, exact at every q, e_0 and t")
    print(
        "sum_r (1 + sign(A_r) cos(2 pi c (t+r)/q)) = q + cos(2 pi c (t - 1/2)/q)/cos(pi c/q) >= q + 1"
    )
    print(" q   e_0    t      left             right            q + 1")
    for q, e0, t in [(7, 0, 0.13), (11, 5, 0.77), (100, 0, 0.5), (100, 37, 0.29)]:
        print(
            "%4d %5d %6.2f %16.9f %16.9f %10d"
            % (q, e0, t, phase_left(q, e0, t), phase_right(q, e0, t), q + 1)
        )
    print()
    print(
        "THE SEAT: the grid scan over t in [0, 1/2], cut 1/%d, Sigma symmetric about t = 1/2"
        % (2 * (ONESTEP_CUT - 1))
    )
    print(
        "the t -> 0 value is 2(q-1) exactly, so B_q(F) >= 2(q-1) at every q and every e_0"
    )
    print(
        "the seat is t = 1/2 at every family printed here except q = 11 and q = 13, where it is interior,"
    )
    print("so t = 1/2 is where the constant is read and not where it is proved to sit")
    print(" q   e_0   sup on the cut     at t      Sigma(1/2)         2(q-1)  seat")
    for q, e0 in ONESTEP_FAMILIES:
        v, at = seat(q, e0)
        h = sigma_t(q, e0, 0.5)
        tag = (
            "t = 1/2"
            if abs(at - 0.5) < 1e-9
            else ("t -> 0" if at < 0.4 else "interior")
        )
        print(
            "%5d %5d %18.10f %9.6f %18.10f %8d  %s"
            % (q, e0, dn(v, 10), at, dn(h, 10), 2 * (q - 1), tag)
        )
    print()
    print("THE SPLIT DEFECT: what the triangle |g_F| <= |D_q| + |g_E| throws away")
    print(
        "K_q = sup_t sum_r |D_q((t+r)/q)| is the exact kernel sup and sum_r |g_E| = q exactly at m = 1,"
    )
    print(
        "so the triangle bound is K_q + q and the defect is (K_q + q) - B_q(F); Phi_q is the proved kernel bound"
    )
    print(
        " q   e_0   Phi_q/q     K_q/q      B_q(F)/q   defect/q   Phi_q slack/q  proved bound/q"
    )
    for q, e0 in ONESTEP_FAMILIES:
        if q < 50:
            continue
        k, _ = ker_seat(q)
        b, _ = seat(q, e0)
        print(
            "%5d %5d %10.6f %10.6f %11.6f %10.6f %13.6f %14.6f"
            % (
                q,
                e0,
                up(phi_q(q) / q),
                dn(k / q),
                dn(b / q),
                up((k + q - b) / q),
                up((phi_q(q) - k) / q),
                up(bound_e(q, e0) / q),
            )
        )
    print()
    print("THE CHORD KERNEL BOUND, what replaces Phi_q inside step 3")
    print(
        "1/sin x <= 1/x + (2/pi)(1 - 2/pi) x on (0, pi/2] is the chord of the convex csc x - 1/x,"
    )
    print(
        "and pairing r with q-1-r puts every shifted-grid argument inside (0, pi/2] at t in (0, 1/2],"
    )
    print(
        "so K(t) = sin(pi t) sum_r 1/sin(pi (t+r)/q) <= (4/pi) q + sin(pi t) Psi'_q with"
    )
    print(
        "Psi'_q = (q/pi)(2 H(P-1) - 1 + 1/P) + (1 - 2/pi) q/2 at even q, P = floor(q/2),"
    )
    print("and (q/pi)(2 H(P-1) - 1 + 2/P) + (1 - 2/pi)(q/2 + 1/(2q)) at odd q")
    print(
        " q     Psi_q/q    Psi'_q/q   (K_q - (4/pi)q)/q   old gap/q   new gap/q   lemma slack/q  at t"
    )
    for q, e0 in ONESTEP_FAMILIES:
        if e0 != 0 or q < 50:
            continue
        k, _ = ker_seat(q)
        v, at = kernel_slack(q)
        e = (k - (4 / pi) * q) / q
        print(
            "%5d %10.6f %11.6f %18.6f %11.6f %11.6f %14.6f %5.3f"
            % (
                q,
                up(psi_q(q) / q),
                up(psi2_q(q) / q),
                dn(e),
                up(psi_q(q) / q - e),
                up(psi2_q(q) / q - e),
                dn(v / q),
                at,
            )
        )
    print(
        "H is the harmonic upper bound ln n + gamma + 1/(2n) of grh-power-saving.md at every Psi and Psi' here"
    )
    print(
        "the monotone step of the sharpening needs Psi'_q >= (1 + pi) q/2, first true at q = %d with that H"
        % mono_floor()
    )
    print(
        "and at q = %d with the harmonic number itself, the two floors the convention separates"
        % mono_floor(exact=True)
    )
    f, fhi = h_argmax_floor()
    print(
        "the hypothesis is sufficient and not necessary: the max of h(tau) sits at tau = 0 at every e_0"
    )
    print(
        "from q = %d up on the exhaustive scan 4..%d, so either floor has room above the true one"
        % (f, fhi)
    )
    print()
    print("THE WEIGHT ROUTE, priced: what dropping the singular terms would cost")
    print(
        "the kept weight w_r = |A_r|/(|A_r| + 1) >= s/(1 + s) >= s/2 is worth q/2 at the seat,"
    )
    print(
        "and dropping it by (1 + sign(A_r) cos) <= 2 costs 2 sum_r 1/(|A_r| + 1) = (2 - 4/pi) q at t = 1/2"
    )
    for q in (1000, 3690, 20000):
        print(
            "  q = %6d   cost/q %.6f   limit 2(1 - 2/pi) = %.6f   net q/2 - cost %.6f q"
            % (
                q,
                up(weight_cost(q) / q),
                up(2 * (1 - 2 / pi)),
                dn(0.5 - weight_cost(q) / q),
            )
        )
    print(
        "  the net is negative at every q printed, so the route loses more than it wins"
    )
    print()
    print("THE PROVED SHARPENING, q >= 17 and m = 1")
    print(
        "B_q(F) <= (4/pi) q + Psi_q + q/2 - sec(pi (e_0 - (q-1)/2)/q)/2 with Psi_q = Phi_q - (4/pi) q - 0.000239,"
    )
    print(
        "against the step 3 bound q PB_q(1) = q + Phi_q: a saving of q/2 at every e_0 and of q/2 + q/pi at e_0 in {0, q-1}"
    )
    print(
        "the wall q_0(a) is the least q with (q-1) q^(-b(a)) > B_q(F)/q, scanned exhaustively from q = 17"
    )
    print(
        "rung a  b(a)       q_0 at step 3   sharpened, every e_0   sharpened, e_0 in {0, q-1}   scan"
    )
    print(
        "held is the count of q in the scan meeting the condition, so held = hi - w + 1 is an up-set"
    )
    for name, b, old, hi in RUNGS:
        w1, n1 = proved_wall(b, True, hi)
        w2, n2 = proved_wall(b, False, hi)
        print(
            "%-7s %.8f %14d %22s %28s   %d..%d  held %d %d"
            % (name, b, old, w1, w2, 17, hi, n1, n2)
        )
    print()
    print("THE SAME WALLS WITH THE CHORD KERNEL BOUND, q >= %d and m = 1" % MONO2)
    print(
        "B_q(F) <= (4/pi) q + Psi'_q + q/2 - sec(pi (e_0 - (q-1)/2)/q)/2, Psi'_q = Psi_q - q/2 + 2/pi at even q"
    )
    print(
        "rung a  b(a)       q_0 at step 3   chord, every e_0   chord, e_0 in {0, q-1}   scan"
    )
    for name, b, old, hi in RUNGS:
        w1, n1 = proved_wall(b, True, hi, True, MONO2)
        w2, n2 = proved_wall(b, False, hi, True, MONO2)
        print(
            "%-7s %.8f %14d %18s %24s   %d..%d  held %d %d"
            % (name, b, old, w1, w2, MONO2, hi, n1, n2)
        )
    print(
        "the cost-out crossing of each chord wall is a kept number of lab/rs/mertens-numerology and is not printed here"
    )
    print()
    print(
        "THE CEILING OF THE LEVER: the same wall computed from the measured B_q(F) itself"
    )
    print("this is a measurement of Sigma(1/2) and not an upper bound certificate,")
    print(
        "so it names what an exact constant could ever buy and never enters a statement"
    )
    print(
        "every e_0 is scanned, not the middle one, and the scan runs downward from the top of"
    )
    print(
        "its range, so the printed wall is the last failure plus one and the tail above it is clear"
    )
    print(
        "rung a  b(a)       ceiling, every e_0   last failure   ceiling, e_0 in {0, q-1}   last failure   scan"
    )
    for name, b, old, hi in RUNGS[:1]:
        w1, f1, d1 = ceiling_wall(b, True, CEIL_HI)
        w2, f2, _ = ceiling_wall(b, False, 4 * CEIL_HI)
        print(
            "%-7s %.8f %20s %8s at e_0 = %d %14s %13s          %d..%d and %d..%d"
            % (name, b, w1, f1, d1, w2, f2, 17, CEIL_HI, 17, 4 * CEIL_HI)
        )
    print()
    print("FALSIFICATION: a family with B_q(F) above the proved bound kills the lever")
    print(
        "the ratio B_q(F)/bound over every e_0 at q = 17..60 on the cut and at the seats above"
    )
    worst = None
    for q in range(17, 61):
        for e0 in range(q):
            v = seat(q, e0, 401)[0]
            if worst is None or v / bound_e(q, e0) > worst[0]:
                worst = (v / bound_e(q, e0), q, e0)
    big = max(
        (seat(q, e0)[0] / bound_e(q, e0), q, e0)
        for q, e0 in ONESTEP_FAMILIES
        if q >= 100
    )
    rng = numpy.random.default_rng(DRAW_SEED)
    draws = None
    for _ in range(DRAWS):
        q = int(rng.choice(DRAW_BASES))
        e0 = int(rng.integers(0, q))
        t = float(rng.uniform(1e-9, 1.0 - 1e-9))
        v = sigma_t(q, e0, t) / bound_e(q, e0)
        if draws is None or v > draws[0]:
            draws = (v, q, e0, t)
    print("  q = 17..60, every e_0:  worst ratio %.6f at q = %d, e_0 = %d" % worst)
    print("  the seats above:        worst ratio %.6f at q = %d, e_0 = %d" % big)
    print(
        "  %d draws, seed %d:    worst ratio %.6f at q = %d, e_0 = %d, t = %.6f"
        % (DRAWS, DRAW_SEED, draws[0], draws[1], draws[2], draws[3])
    )
    print("  no family reaches 1, so nothing measured contradicts the sharpening")
    print()
    print(
        "FALSIFICATION OF THE CHORD BOUND, the same three sweeps against (4/pi) q + Psi'_q + q/2 - sec/2"
    )
    worst = None
    for q in range(MONO2, 61):
        for e0 in range(q):
            v = seat(q, e0, 401)[0]
            if worst is None or v / bound2_e(q, e0) > worst[0]:
                worst = (v / bound2_e(q, e0), q, e0)
    big = max(
        (seat(q, e0)[0] / bound2_e(q, e0), q, e0)
        for q, e0 in ONESTEP_FAMILIES
        if q >= 100
    )
    rng = numpy.random.default_rng(DRAW_SEED)
    draws = None
    for _ in range(DRAWS):
        q = int(rng.choice(DRAW_BASES))
        e0 = int(rng.integers(0, q))
        t = float(rng.uniform(1e-9, 1.0 - 1e-9))
        v = sigma_t(q, e0, t) / bound2_e(q, e0)
        if draws is None or v > draws[0]:
            draws = (v, q, e0, t)
    print(
        "  q = %d..60, every e_0: worst ratio %.6f at q = %d, e_0 = %d"
        % ((MONO2,) + worst)
    )
    print("  the seats above:        worst ratio %.6f at q = %d, e_0 = %d" % big)
    print(
        "  %d draws, seed %d:    worst ratio %.6f at q = %d, e_0 = %d, t = %.6f"
        % (DRAWS, DRAW_SEED, draws[0], draws[1], draws[2], draws[3])
    )
    print(
        "  no family reaches 1, so nothing measured contradicts the chord bound either"
    )


VERBS["onestep"] = onestep

# PERDEN

PERDEN_FAMILIES = [(3, [2], 12), (10, [9], 6), (101, [0], 3), (1499, [749], 2)]

PERDEN_POINTWISE = [(10, [9], 5), (10, [9], 6)]

PERDEN_HONEST = [(3, [2], [6, 8, 10, 12]), (10, [9], [4, 5, 6])]

PERDEN_BRUTE = [(3, [2], 4), (3, [2], 5), (10, [9], 2)]

PERDEN_RUNGS = [
    ((1, 2), (3, 4), "both"),
    ((13, 25), (1417, 1850), "Zhang"),
    ((11, 20), (913, 1160), "Zhang"),
    ((4, 7), (4, 5), "both"),
    ((3, 5), (4, 5), "BH"),
    ((2, 3), (5, 6), "BH"),
]


def hat_abs(q, E, j):
    import numpy

    n = q**j
    v = numpy.arange(n, dtype=numpy.int64)
    out = numpy.ones(n, dtype=numpy.complex128)
    for i in range(j):
        u = (v * pow(q, i, n)) % n
        t = u.astype(numpy.float64) / n
        z = numpy.exp(2j * numpy.pi * t)
        den = numpy.where(u == 0, 1.0 + 0j, 1.0 - z)
        g = numpy.where(
            u == 0, complex(q, 0), (1.0 - numpy.exp(2j * numpy.pi * q * t)) / den
        )
        for e in E:
            g = g - z**e
        out *= g
    return numpy.abs(out)


def levmass(q, E, J):
    import numpy

    C, c = [], []
    for j in range(J + 1):
        h = hat_abs(q, E, j)
        C.append(float(h.sum()))
        c.append(1.0 if j == 0 else float(h[numpy.arange(q**j) % q != 0].sum()))
        del h
    return C, c


def smalldens(q, E, L):
    import numpy

    n = q**L
    w = hat_abs(q, E, L)
    a = numpy.arange(n, dtype=numpy.int64)
    d = n // numpy.gcd(a, n)
    bad = ((a % q) != 0) & (d.astype(numpy.float64) <= float(n) ** 0.5)
    return int(bad.sum()), float(w[bad].sum()) / float(w.sum())


def nu(N):
    from math import sqrt

    import numpy

    a = numpy.arange(N, dtype=numpy.int64)
    best = numpy.full(N, float(N), dtype=numpy.float64)
    for Q in range(1, int(2.0 * sqrt(N)) + 2):
        u = (a * Q) % N
        numpy.minimum(best, numpy.minimum(u, N - u).astype(numpy.float64) + Q, out=best)
    return best


def nu_slow(N):
    from math import gcd

    out = []
    for x in range(N):
        bst = float(N)
        for Q in range(1, N + 1):
            for r in range(Q + 1):
                if gcd(r, Q) == 1 or (Q == 1 and r == 0):
                    bst = min(bst, Q + abs(x * Q - r * N))
        out.append(bst)
    return out


def honest(q, E, L, a, b):
    import numpy

    N = q**L
    w = hat_abs(q, E, L)
    con = numpy.minimum(float(N) ** b, float(N) ** a * numpy.sqrt(nu(N)))
    return float((w * con).sum()) / (float(w.sum()) * float(N) ** b)


def brute(q, E, j):
    from cmath import exp
    from math import pi

    n, k = q**j, q - len(E)
    D = [d for d in range(q) if d not in E]
    vals = [0]
    for i in range(j):
        vals = [v + d * q**i for v in vals for d in D]
    assert len(vals) == k**j
    tot, top = 0.0, 0.0
    for a in range(n):
        z = abs(sum(exp(2j * pi * a * v / n) for v in vals))
        tot += z
        if a % q:
            top += z
    return tot, top


def shares(q, E, C, c, L):
    from math import exp, log

    k = q - len(E)
    return [exp((L - j) * log(k) + log(c[j]) - log(C[L])) for j in range(L + 1)]


def perden_exp(q, sh, L, a, b):
    from math import log

    lq = log(q)
    gain = [
        log(sh[j]) / (L * lq) - max(0.0, b - a - j / (2.0 * L)) for j in range(L + 1)
    ]
    ratio = sum(
        sh[j] * q ** (-L * max(0.0, b - a - j / (2.0 * L))) for j in range(L + 1)
    )
    jm = max(range(L + 1), key=lambda j: gain[j])
    return jm, gain[jm], ratio


def perden():
    from fractions import Fraction as Fr
    from math import log

    import numpy

    print("THE PER-DENOMINATOR MOBIUS INPUT IN STEP 5 (T-M10)")
    print(
        "x = q^L, S_mu(theta) = sum_(n < x) mu(n) e(n theta), M_F(q^L) = q^(-L) sum_(a mod q^L) hatF_L(a/q^L) S_mu(-a/q^L)"
    )
    print(
        "group a by its q-power level j: a = a' q^(L-j) with q not dividing a', so a/q^L = a'/q^j"
    )
    print(
        "hatF_L(a'/q^j) = k^(L-j) hatF_j(a'/q^j) because g_F(integer) = k, so the level-j mass is k^(L-j) c_j"
    )
    print(
        "with c_j = sum over a' mod q^j, q not dividing a', of |hatF_j(a'/q^j)|, c_0 = 1, C_j = sum_(i <= j) k^(j-i) c_i"
    )
    print(
        "Baker-Harman PROPOSITION p.194 eq. 6 under hypothesis (4), L(s,chi) zero-free in sigma > a for EVERY"
    )
    print(
        "Dirichlet character: S_mu(theta) << x^(a+eps) Q^(1/2) (1 + x|theta - r/Q|)^(1/2) at (r,Q) = 1;"
    )
    print(
        "the reduced denominator of a'/q^j divides q^j, so the level-j constant is AT MOST"
    )
    print(
        "min(x^(b(a)), x^a q^(j/2)), taking (r,Q) = the frequency itself, where the second factor is 1,"
    )
    print(
        "and the first the uniform THEOREM p.193; that choice of (r,Q) is the exact-frequency COROLLARY,"
    )
    print(
        "and the bracket [m/q, 1] below is proved for it and for it alone. The full PROPOSITION lets ANY"
    )
    print(
        "reduced r/Q serve any frequency: per-frequency constant x^a nu(a)^(1/2) with nu(a) = min_Q (Q +"
    )
    print(
        "||aQ||_(q^L)), since Q(1 + x|a/q^L - r/Q|) = Q + |aQ - r q^L|. That form is measured, not proved,"
    )
    print("in the honest block below, and it buys no exponent either")
    print()
    print(
        "|M_F(q^L)| <<_(q,eps) x^eps q^(-L) sum_(j = 0..L) k^(L-j) c_j min(x^(b(a)), x^a q^(j/2))"
    )
    print()
    print("THE TOP LEVEL PAYS THE UNIFORM PRICE AT EVERY RUNG")
    print(
        "the level-j charge is a + j/(2L), so the level beats b(a) only below j/L = 2(b - a);"
    )
    print(
        "at j = L the charge reads a + 1/2, and a + 1/2 - b(a) > 0 at every rung, exact rationals"
    )
    print(
        "   a     b(a)      source  a + 1/2 - b(a)   crossing 2(b - a)   beats uniform at j = L"
    )
    for (an, ad), (bn, bd), src in PERDEN_RUNGS:
        a, b = Fr(an, ad), Fr(bn, bd)
        assert a + Fr(1, 2) > b and 2 * (b - a) <= Fr(1, 2)
        print(
            "%6s %9s %8s %15s %19s   %s"
            % (
                a,
                b,
                src,
                a + Fr(1, 2) - b,
                2 * (b - a),
                "yes" if a + Fr(1, 2) < b else "no",
            )
        )
    print()
    print(
        "the charge is an upper bound, not the pointwise truth. At composite q the reduced denominator of"
    )
    print(
        "a'/q^L with q not dividing a' can sit far below q^L: at q = 10, L = 6, a' = 5^6 = 15625 gives"
    )
    print(
        "15625/10^6 = 1/64, denominator 2^6 = x^0.301 and charge x^(a + 0.1505), beneath x^(3/4) at a = 1/2."
    )
    print(
        "The displayed inequality is unharmed, being an upper bound, and the mass such frequencies carry is"
    )
    print(
        " q     E    L   top-level a' with reduced denominator <= x^(1/2)    their share of C_L"
    )
    for q, E, L in PERDEN_POINTWISE:
        n, share = smalldens(q, E, L)
        print("%5d %5d %3d %40d %24.6e" % (q, E[0], L, n, share))
    print()
    print(
        "THE LEVEL PROFILE: the l^1 mass is NOT at the low denominators; it decays geometrically downward"
    )
    print(
        "proved: sum_(s mod q) |g_F((t+s)/q)|^2 = qk exactly and |g_F| <= k, so sum_s |g_F| >= q,"
    )
    print(
        "hence C_j >= q C_(j-1) at every j, hence c_j = C_j - k C_(j-1) >= (m/q) C_j,"
    )
    print("and sum_(i <= J) k^(j-i) c_i = k^(j-J) C_J <= (k/q)^(j-J) C_j")
    print(
        "the share column runs over j <= floor(L/2), the levels whose charge is at most x^(3/4), the tie at"
    )
    print(
        "j = L/2 included, so it over-reports the levels the PROPOSITION strictly beats the uniform at;"
    )
    print(
        "the q = 101 and q = 1499 families stop at j = 3 and j = 2, where C_j/C_(j-1) is still moving,"
    )
    print(
        "244.658399 then 234.507307 at q = 101, so their top shares are short rows and not constants"
    )
    print(
        " q     E    j   C_j          C_j/C_(j-1)  q     c_j/C_j    floor m/q   den<=x^(1/2)  cap (k/q)^(j-J)"
    )
    prof = {}
    for q, E, J in PERDEN_FAMILIES:
        C, c = levmass(q, E, J)
        prof[(q, E[0])] = (C, c, J)
        k, m = q - len(E), len(E)
        for j in range(1, J + 1):
            sh = shares(q, E, C, c, j)
            Jc = j // 2
            low = sum(sh[: Jc + 1])
            dec = sum(k ** (j - i) * c[i] for i in range(j + 1))
            assert abs(dec - C[j]) <= 1e-9 * C[j]
            assert C[j] >= q * C[j - 1] * (1 - 1e-12)
            assert c[j] / C[j] >= m / q - 1e-12
            assert low <= (k / q) ** (j - Jc) + 1e-12
            print(
                "%5d %5d %3d %13.4f %11.6f %6d %10.6f %11.6f %13.6f %14.6f"
                % (
                    q,
                    E[0],
                    j,
                    C[j],
                    C[j] / C[j - 1],
                    q,
                    c[j] / C[j],
                    m / q,
                    low,
                    (k / q) ** (j - Jc),
                )
            )
    print()
    print("THE EXPONENT IT BUYS, against the uniform step 5")
    print(
        "unif = log_q(C_L/q^L)/L + b(a), den = max_j [log_q(k^(L-j) c_j/q^L)/L + min(b(a), a + j/(2L))]"
    )
    print(
        "ratio = the per-denominator sum divided by C_L x^(b(a)), proved to lie in [m/q, 1]"
    )
    print(
        "proved: every level gain is at most 0 and the top level's is log_q(c_L/C_L)/L, so the saving is at"
    )
    print(
        "MOST -log_q(top share)/L, itself at most log_q(q/m)/L; that the largest term is the top level is"
    )
    print(
        "printed as argmax j and not proved. Read on the whole bound rather than on its largest term the"
    )
    print(
        "saving is log_q(ratio)/L, smaller again, the two agreeing only through the x^eps that absorbs L+1"
    )
    print("terms. Either reading is a single factor at most q/m and never an exponent")
    print(
        "q     E    L   a      b(a)     unif        den         den - unif  L(den-unif)  -log_q top  cap log_q(q/m)  argmax j  ratio     log_q(ratio)/L"
    )
    for q, E, J in PERDEN_FAMILIES:
        C, c, _ = prof[(q, E[0])]
        k, m = q - len(E), len(E)
        a, b = 0.5, 0.75
        for L in [J] if J < 6 else [J - 2, J]:
            sh = shares(q, E, C, c, L)
            jm, g, ratio = perden_exp(q, sh, L, a, b)
            unif = log(C[L]) / (L * log(q)) - 1 + b
            sp = -log(sh[L]) / log(q)
            assert -log(q / m) / log(q) / L - 1e-12 <= -sp / L <= g <= 1e-12
            assert m / q - 1e-12 <= ratio <= 1 + 1e-12
            print(
                "%5d %5d %3d %6.3f %8.5f %11.6f %11.6f %11.6f %12.6f %11.6f %15.6f %9d %9.6f %11.6f"
                % (
                    q,
                    E[0],
                    L,
                    a,
                    b,
                    unif,
                    unif + g,
                    g,
                    L * g,
                    sp,
                    log(q / m) / log(q),
                    jm,
                    ratio,
                    log(ratio) / (L * log(q)),
                )
            )
    print()
    print(
        "THE FULL MINOR-ARC FORM, the same tool at full strength (measured, not proved)"
    )
    print(
        "any reduced r/Q may serve any frequency, so the per-frequency constant is min(x^(b(a)), x^a nu^(1/2))"
    )
    print(
        "with nu(a) = min_Q (Q + ||aQ||_(q^L)); dropping the coprimality changes nothing, a smaller Q being"
    )
    print(
        "never worse, and the scan Q <= 2 q^(L/2) is exact because Dirichlet gives nu(a) <= 2 q^(L/2)"
    )
    assert float(numpy.abs(nu(81) - numpy.array(nu_slow(81))).max()) == 0.0
    print(
        "nu checked against a full search over every reduced r/Q at q^L = 81: 0 mismatches"
    )
    print(
        "ratio_nu = that weighted sum over C_L x^(b(a)); the exact-frequency ratio is the column beside it"
    )
    print(
        "q     E    L   ratio_nu   log_q(ratio_nu)/L   L times that   exact-frequency ratio"
    )
    for q, E, Ls in PERDEN_HONEST:
        C, c, _ = prof[(q, E[0])]
        for L in Ls:
            rn = honest(q, E, L, 0.5, 0.75)
            gn = log(rn) / (L * log(q))
            _, _, ratio = perden_exp(q, shares(q, E, C, c, L), L, 0.5, 0.75)
            assert 0.0 < rn <= ratio + 1e-12
            print(
                "%5d %5d %3d %10.6f %19.6f %14.6f %22.6f"
                % (q, E[0], L, rn, gn, L * gn, ratio)
            )
    print(
        "ratio_nu rises with L while L times the gain falls, so the full form buys a bounded factor too and"
    )
    print(
        "no exponent; that ratio_nu is bounded below in L is measured over these rows and is not proved"
    )
    print()
    print("THE TRANSFORM, CHECKED WITHOUT ITSELF: brute force over the digit strings")
    print(
        "hatF_j(a/q^j) summed over the k^j allowed n directly, against the product form hat_abs"
    )
    print("q     E    j   C_j brute     C_j product   c_j/C_j brute  c_j/C_j product")
    for q, E, j in PERDEN_BRUTE:
        tot, top = brute(q, E, j)
        C, c = levmass(q, E, j)
        assert abs(tot - C[j]) <= 1e-9 * C[j] and abs(top - c[j]) <= 1e-9 * C[j]
        print(
            "%5d %5d %3d %13.6f %13.6f %14.6f %16.6f"
            % (q, E[0], j, tot, C[j], top / tot, c[j] / C[j])
        )
    print()
    print(
        "WHERE THE PROPOSITION DOES BELONG: the d-form minor arc, not the q-power grid"
    )
    print(
        "on a Dirichlet arc |theta - l/d| <= 1/d^2 with (l,d) = 1 the PROPOSITION reads"
    )
    print(
        "S_mu(theta) << x^(a+eps) d^(1/2) (1 + x/d^2)^(1/2) = x^(a+eps) (d + x/d)^(1/2) <= x^(a+eps) (d^(1/2) + x^(1/2) d^(-1/2))"
    )
    print(
        "and at a = 1/2 that is x^eps ((x d)^(1/2) + x d^(-1/2)), the first and third terms of the d-form"
    )
    print(
        "minor-arc bound Theorem L5 owes, whose middle term x^(4/5) is the unconditional Vaughan piece;"
    )
    print(
        "this row applies eq. 6 at an arbitrary arc denominator d, so it holds only within the printed"
    )
    print(
        "range on Q that the source carries and is unread here, unlike the blocks above, where the"
    )
    print("min caps the PROPOSITION by the uniform THEOREM at the level that decides")
    print()
    print("THE WALL IT MOVES: none, and no wall row is printed here")
    print(
        "the exponent is b(a) + c_q with the SAME c_q, so the certificate (q-1) q^(-b) - PB_q(1, e_0) > 0"
    )
    print(
        "takes no per-denominator quantity at all and the GRH walls stay exactly where verb onestep's chord"
    )
    print("put them; re-evaluating that certificate here would discriminate nothing")


VERBS["perden"] = perden

if __name__ == "__main__":
    VERBS[sys.argv[1] if len(sys.argv) > 1 else "split"]()
