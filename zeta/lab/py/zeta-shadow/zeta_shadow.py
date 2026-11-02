import os
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "design-zeta"))

from design_zeta import Design

# THE LADDER

LADDER = [
    (5, (0, 1)),
    (4, (0, 1)),
    (3, (0, 1)),
    (4, (0, 1, 2)),
    (5, (0, 1, 2, 3)),
    (10, tuple(range(8))),
    (10, tuple(range(9))),
    (20, tuple(range(19))),
    (50, tuple(range(49))),
    (2, (0, 1)),
]

NEW = [(20, tuple(range(19))), (50, tuple(range(49)))]

OLD = [d for d in LADDER if d not in NEW]

SETS = {"ladder": LADDER, "new": NEW, "old": OLD}

NZ = 12

CAPZ = {(20, tuple(range(19))): 6, (50, tuple(range(49))): 6}

DIFF = mp.mpf("1e-6")
TRUST = mp.mpf("0.6")


def line(*a):
    print(" ".join(str(x) for x in a), flush=True)


def tag(q, F):
    miss = tuple(d for d in range(q) if d not in F)
    if not miss:
        return "base " + str(q) + " full"
    if len(miss) <= len(F):
        return "base " + str(q) + " missing {" + ",".join(str(d) for d in miss) + "}"
    return "base " + str(q) + " {" + ",".join(str(d) for d in F) + "}"


def med(xs):
    v = sorted(xs)
    n = len(v)
    if n % 2:
        return v[n // 2]
    return (v[n // 2 - 1] + v[n // 2]) / 2


def words(q, F, L):
    cur = [0]
    for _ in range(L):
        cur = [q * m + a for m in cur for a in F]
    return cur


# THE PRINCIPAL FIBRE


def fibre(q, F, L):
    return (mp.mpf(len(F)) / q) ** L


def arc(q, F, L):
    N = q**L
    tot = mp.mpf(1) / N
    for n in words(q, F, L):
        if n:
            tot += mp.sin(mp.pi * mp.mpf(n) / N) / (mp.pi * n)
    return tot


def check(q, F, L):
    mp.mp.dps = 40
    N = q**L
    D = set(words(q, F, L))
    G = [mp.mpf(0)] * N
    for a in range(N):
        G[a] = mp.fsum([mp.e ** (2j * mp.pi * m * a / N) for m in D])
    worst = mp.mpf(0)
    for n in range(N):
        v = mp.fsum([G[a] * mp.e ** (-2j * mp.pi * n * a / N) for a in range(N)]) / N
        worst = max(worst, abs(v - (1 if n in D else 0)))
    return worst, G[0] / N


def mass(which):
    line(
        "MASS the principal fibre of the position identity. zeta_(F,L)(s) = q^(-L) sum_(a mod q^L)"
    )
    line(
        "  G_L(a/q^L) S_L(s, a/q^L) splits the LEVEL-L polynomial against the partial sum of zeta to"
    )
    line(
        "  q^L, and the a = 0 fibre carries G_L(0)/q^L = (k/q)^L exactly. The level weight falls to 0,"
    )
    line(
        "  so no level is forced: the L = 1 reading c = k/q is a choice and the sweep tests it."
    )
    line(
        "  ARC is the continuous form, int over abs(t) < 1/(2 q^L) of G_L, exact as a sinc sum;"
    )
    line(
        "  the full set's own readings run toward Si(pi)/pi =",
        mp.nstr(mp.si(mp.pi) / mp.pi, 10),
    )
    for q, F in which:
        k = len(F)
        c = mp.mpf(k) / q
        if k**2 <= 4096 and q**2 <= 4096:
            w, f = check(q, F, 2)
            line(
                "   IDENTITY",
                tag(q, F),
                "L 2 largest abs 1_(D_L) - inverse transform",
                mp.nstr(w, 4),
                "a = 0 fibre",
                mp.nstr(f, 12),
                "against (k/q)^2",
                mp.nstr(c**2, 12),
            )
        row = [
            tag(q, F),
            "k/q",
            mp.nstr(c, 10),
            "alpha",
            mp.nstr(mp.log(k) / mp.log(q), 10),
        ]
        for L in (1, 2, 3):
            if k**L > 200000:
                break
            a = arc(q, F, L)
            row += [
                "| L",
                L,
                "fibre",
                mp.nstr(fibre(q, F, L), 8),
                "arc",
                mp.nstr(a, 8),
                "arc/fibre",
                mp.nstr(a / fibre(q, F, L), 8),
            ]
        line("  ", *row)


# THE FIRST-ORDER SHADOW


def deriv(d, s):
    return (d.zeta(s + DIFF)[0] - d.zeta(s - DIFF)[0]) / (2 * DIFF)


def pole_gap(d, s):
    per = 2 * mp.pi / mp.log(d.q)
    best = mp.inf
    for m in range(3):
        x = d.alpha - m
        j = mp.nint(mp.im(s) * mp.log(d.q) / (2 * mp.pi))
        for jj in (j - 1, j, j + 1):
            best = min(best, abs(s - mp.mpc(x, jj * per)))
    return best


def newton(d, s0, trust=TRUST, cap=14):
    s = s0
    for _ in range(cap):
        v, _ = d.zeta(s)
        g = deriv(d, s)
        if g == 0:
            return None
        step = v / g
        if abs(step) > trust:
            step = step * trust / abs(step)
        s = s - step
        if abs(step) < mp.mpf("1e-28"):
            break
    v, _ = d.zeta(s)
    if abs(v) > mp.mpf("1e-16"):
        return None
    if abs(s - s0) > mp.mpf("1.5") or pole_gap(d, s) < mp.mpf("0.02"):
        return None
    return s


def predict(which, nz):
    line(
        "PREDICT at a zeta zero rho_0 one has zeta(rho_0) = 0, so for ANY constant c the split"
    )
    line(
        "  zeta_F = c zeta + E_F gives E_F(rho_0) = zeta_F(rho_0), and a zero of zeta_F near rho_0"
    )
    line(
        "  sits at rho_0 - zeta_F(rho_0)/(c zeta'(rho_0)) to first order. STEP is the constant-free"
    )
    line(
        "  reading c zeta' -> zeta_F', that is Newton's own first step -zeta_F(rho_0)/zeta_F'(rho_0),"
    )
    line(
        "  and it is the primary column. PRED is the same law at the L = 1 fibre reading c = k/q."
    )
    line(
        "  Both see no design zero. FOUND is the zero Newton reaches from rho_0, accepted only at"
    )
    line(
        "  abs(zeta_F) < 1e-16, within 1.5 of rho_0 and 0.02 clear of the pole lattice; a miss is a"
    )
    line(
        "  zero the trust region does not reach and its rung's medians are conditioned on that."
    )
    out = {}
    for q, F in which:
        d = Design(q, F)
        c = mp.mpf(len(F)) / q
        n_use = CAPZ.get((q, F), nz)
        line(
            " ",
            tag(q, F),
            "alpha",
            mp.nstr(d.alpha, 10),
            "k/q",
            mp.nstr(c, 8),
            "zeros",
            n_use,
            "top gamma",
            mp.nstr(mp.im(mp.zetazero(n_use)), 12),
        )
        rows = []
        for n in range(1, n_use + 1):
            g = mp.im(mp.zetazero(n))
            rho = mp.mpc(mp.mpf("0.5"), g)
            v, eb = d.zeta(rho)
            if abs(v) < mp.mpf("1e-20"):
                line(
                    "     n",
                    n,
                    "gamma",
                    mp.nstr(g, 12),
                    "abs zeta_F",
                    mp.nstr(abs(v), 4),
                    "CONTROL E_F = 0 and both predicted offsets are 0",
                )
                continue
            zp = mp.zeta(rho, derivative=1)
            fp = deriv(d, rho)
            cpl = fp / zp
            step = -v / fp
            pred = -v / (c * zp)
            s = newton(d, rho)
            if s is None:
                line(
                    "     n",
                    n,
                    "gamma",
                    mp.nstr(g, 12),
                    "absE",
                    mp.nstr(abs(v), 8),
                    "coupling",
                    mp.nstr(cpl, 8),
                    "step",
                    mp.nstr(abs(step), 8),
                    "pred",
                    mp.nstr(abs(pred), 8),
                    "FOUND none",
                )
                continue
            off = s - rho
            rows.append((step, pred, off, cpl, abs(v), eb, s))
            line(
                "     n",
                n,
                "gamma",
                mp.nstr(g, 12),
                "absE",
                mp.nstr(abs(v), 8),
                "coupling",
                mp.nstr(cpl, 8),
                "off",
                mp.nstr(off, 8),
                "step ratio",
                mp.nstr(abs(step) / abs(off), 8),
                "pred ratio",
                mp.nstr(abs(pred) / abs(off), 8),
                "abs Im off / abs Re off",
                mp.nstr(abs(mp.im(off)) / abs(mp.re(off)), 8),
            )
        out[(q, F)] = rows
        if rows:
            r0 = sorted(abs(st) / abs(o) for st, p, o, cl, av, eb, s in rows)
            r1 = sorted(abs(p) / abs(o) for st, p, o, cl, av, eb, s in rows)
            iz = sorted(
                abs(mp.im(o)) / abs(mp.re(o)) for st, p, o, cl, av, eb, s in rows
            )
            line(
                "     SUMMARY located",
                len(rows),
                "of",
                n_use,
                "| STEP median ratio",
                mp.nstr(med(r0), 8),
                "band [",
                mp.nstr(r0[0], 6),
                ",",
                mp.nstr(r0[-1], 6),
                "] largest abs ratio - 1",
                mp.nstr(max(abs(x - 1) for x in r0), 6),
                "| PRED median ratio",
                mp.nstr(med(r1), 8),
                "band [",
                mp.nstr(r1[0], 6),
                ",",
                mp.nstr(r1[-1], 6),
                "] largest abs ratio - 1",
                mp.nstr(max(abs(x - 1) for x in r1), 6),
                "| median abs coupling - k/q",
                mp.nstr(med([abs(cl - c) for st, p, o, cl, av, eb, s in rows]), 8),
                "median abs coupling - 1",
                mp.nstr(med([abs(cl - 1) for st, p, o, cl, av, eb, s in rows]), 8),
                "| median abs off",
                mp.nstr(med([abs(o) for st, p, o, cl, av, eb, s in rows]), 8),
                "median absE",
                mp.nstr(med([av for st, p, o, cl, av, eb, s in rows]), 8),
                "largest ladder bound",
                mp.nstr(max(eb for st, p, o, cl, av, eb, s in rows), 4),
                "| abs Im off / abs Re off band [",
                mp.nstr(iz[0], 6),
                ",",
                mp.nstr(iz[-1], 6),
                "]",
            )
    return out


# THE RUNGS

ZZ = None


def ordinates():
    global ZZ
    if ZZ is None:
        ZZ = [mp.im(mp.zetazero(n)) for n in range(1, 40)]
    return ZZ


def shadow_stat(rows):
    zz = ordinates()
    ims = sorted(mp.im(s) for st, p, o, cl, av, eb, s in rows)
    band = [t for t in zz if ims[0] - 2 < t < ims[-1] + 2]
    if len(band) < 2 or len(ims) < 3:
        return None, None
    null = (band[-1] - band[0]) / (len(band) - 1) / 4
    ds = [min(abs(t - y) for y in zz) for t in ims]
    return sum(ds) / len(ds), null


def loglog(xs, ys):
    lx = [mp.log(x) for x in xs]
    ly = [mp.log(y) for y in ys]
    n = len(lx)
    mx = sum(lx) / n
    my = sum(ly) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    sxx = sum((a - mx) ** 2 for a in lx)
    b = sxy / sxx
    a = my - b * mx
    res = sum((y - (a + b * x)) ** 2 for x, y in zip(lx, ly))
    tot = sum((y - my) ** 2 for y in ly)
    return b, 1 - res / tot


def rungs(which, nz):
    got = predict(which, nz)
    line(
        "RATE the median paired offset against the missing-digit density m/q = 1 - k/q and against"
    )
    line(
        "  1 - alpha. The two normalisations differ only by log q, so the ladder discriminates them"
    )
    line(
        "  only over the range of (1 - alpha)/(m/q), which is printed. Every median is over the"
    )
    line("  zeros the trust region reaches, and the rungs do not share one height.")
    xs, ys, zs = [], [], []
    for q, F in which:
        rows = got.get((q, F))
        if not rows:
            continue
        a = mp.log(len(F)) / mp.log(q)
        mq = 1 - mp.mpf(len(F)) / q
        if mq == 0:
            continue
        m = med([abs(o) for st, p, o, cl, av, eb, s in rows])
        xs.append(mq)
        ys.append(m)
        zs.append(1 - a)
        line(
            "  ",
            tag(q, F),
            "alpha",
            mp.nstr(a, 10),
            "n",
            len(rows),
            "1-alpha",
            mp.nstr(1 - a, 8),
            "m/q",
            mp.nstr(mq, 8),
            "(1-alpha)/(m/q)",
            mp.nstr((1 - a) / mq, 8),
            "median absE",
            mp.nstr(med([av for st, p, o, cl, av, eb, s in rows]), 8),
            "median abs off",
            mp.nstr(m, 8),
            "/(m/q)",
            mp.nstr(m / mq, 8),
            "/(1-alpha)",
            mp.nstr(m / (1 - a), 8),
            "| median abs Im off",
            mp.nstr(med([abs(mp.im(o)) for st, p, o, cl, av, eb, s in rows]), 8),
            "median abs Re s - 1/2",
            mp.nstr(
                med(
                    [abs(mp.re(s) - mp.mpf("0.5")) for st, p, o, cl, av, eb, s in rows]
                ),
                8,
            ),
        )
    if len(xs) > 2:
        b1, r1 = loglog(xs, ys)
        b2, r2 = loglog(zs, ys)
        c1 = sorted(y / x for x, y in zip(xs, ys))
        c2 = sorted(y / x for x, y in zip(zs, ys))
        g = sorted(z / x for x, z in zip(xs, zs))
        line(
            "   FIT a least squares in the logs, a fit and not a theorem: median abs off scales as"
        )
        line(
            "    (m/q)^",
            mp.nstr(b1, 6),
            "with R2",
            mp.nstr(r1, 6),
            "and as (1-alpha)^",
            mp.nstr(b2, 6),
            "with R2",
            mp.nstr(r2, 6),
        )
        line(
            "    the m/q column spans",
            mp.nstr(c1[-1] / c1[0], 6),
            "and the 1-alpha column",
            mp.nstr(c2[-1] / c2[0], 6),
            "while (1-alpha)/(m/q) itself spans",
            mp.nstr(g[-1] / g[0], 6),
            "over the ladder, so the discrimination is",
            mp.nstr((c2[-1] / c2[0]) / (c1[-1] / c1[0]), 6),
            "inside that gap",
        )
    line(
        "LADDER the shadow statistic of the family row on the zeros this study pairs: mean distance"
    )
    line(
        "  from a located design ordinate to the nearest zeta ordinate over a quarter of the mean"
    )
    line(
        "  gap between consecutive zeta ordinates in the range. The pairing is zeta-zero-first where"
    )
    line(
        "  the family row's is design-zero-first, so this is a parallel ladder and not that row."
    )
    for q, F in which:
        rows = got.get((q, F))
        if not rows:
            continue
        a = mp.log(len(F)) / mp.log(q)
        md, null = shadow_stat(rows)
        if md is None:
            continue
        line(
            "  ",
            tag(q, F),
            "alpha",
            mp.nstr(a, 10),
            "n",
            len(rows),
            "mean dist",
            mp.nstr(md, 8),
            "null gap/4",
            mp.nstr(null, 8),
            "ratio to null",
            mp.nstr(md / null, 8),
        )
    line(
        "SHARPNESS the constant-free step against the k/q reading, pooled by the size of the offset"
    )
    pool = []
    for q, F in which:
        for st, p, o, cl, av, eb, s in got.get((q, F)) or []:
            pool.append((abs(o), abs(st) / abs(o), abs(p) / abs(o)))
    for lo, hi in (
        (0, mp.mpf("0.05")),
        (mp.mpf("0.05"), mp.mpf("0.1")),
        (mp.mpf("0.1"), mp.mpf("0.2")),
        (mp.mpf("0.2"), mp.mpf("0.4")),
        (mp.mpf("0.4"), mp.inf),
    ):
        b = [(r0, r1) for u, r0, r1 in pool if lo <= u < hi]
        if not b:
            continue
        line(
            "   abs off in [",
            mp.nstr(lo, 3),
            ",",
            mp.nstr(hi, 3),
            ") n",
            len(b),
            "STEP median",
            mp.nstr(med([x[0] for x in b]), 8),
            "largest abs ratio - 1",
            mp.nstr(max(abs(x[0] - 1) for x in b), 6),
            "| PRED median",
            mp.nstr(med([x[1] for x in b]), 8),
            "largest abs ratio - 1",
            mp.nstr(max(abs(x[1] - 1) for x in b), 6),
        )


# STUDY


def main():
    argv = sys.argv[1:]
    verb = argv[0] if argv else "mass"
    which = SETS[argv[1]] if len(argv) > 1 else SETS["ladder"]
    t0 = time.time()
    if verb == "mass":
        mass(which)
    if verb == "predict":
        predict(which, NZ)
    if verb == "rungs":
        rungs(which, NZ)
    line("seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
