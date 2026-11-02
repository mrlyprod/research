import os
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "design-zeta"))
sys.path.insert(0, os.path.join(HERE, "..", "zeta-locus"))

from design_zeta import strip
from zeta_locus import (
    ZCache,
    alpha_of,
    build,
    comb,
    disc,
    line,
    polish,
    seek,
    tag,
    teeth,
)

RHO = mp.mpf("0.45")
RHO2 = mp.mpf("0.6")
LO = mp.mpf("-0.92")
HI = mp.mpf("3.02")
TOLT = mp.mpf("0.05")
TOLS = mp.mpf("0.05")

Q3 = [(3, (1,)), (3, (0, 1)), (3, (1, 2)), (3, (0, 1, 2))]
Q4 = [
    (4, (1,)),
    (4, (0, 1)),
    (4, (1, 2)),
    (4, (1, 3)),
    (4, (2, 3)),
    (4, (0, 1, 2)),
    (4, (0, 1, 3)),
    (4, (0, 2, 3)),
    (4, (1, 2, 3)),
    (4, (0, 1, 2, 3)),
]
Q5 = [(5, (0, 1)), (5, (1, 2))]
HALF = [(9, (0, 1, 2)), (16, (0, 1, 2, 3))]
WIDE = [(10, tuple(range(9)))]
CONTROL = [(2, (0, 1))]
RUNGS = [(5, (0, 1, 2, 3)), (10, tuple(range(8)))]
NEAR = [
    (3, (0, 1)),
    (4, (0, 1, 2)),
    (5, (0, 1, 2, 3)),
    (10, tuple(range(8))),
    (10, tuple(range(9))),
    (2, (0, 1)),
    (3, (0, 1, 2)),
    (4, (0, 1, 2, 3)),
]
FIX = [
    (3, (0, 1)),
    (4, (0, 1)),
    (4, (1, 2)),
    (4, (1, 3)),
    (4, (2, 3)),
    (5, (0, 1)),
    (4, (0, 1, 2)),
    (5, (0, 1, 2, 3)),
    (10, tuple(range(8))),
    (2, (0, 1)),
]
TWO = [(3, (0, 1)), (4, (0, 1)), (4, (1, 2)), (4, (1, 3)), (4, (2, 3))]
LADDER = [
    (3, (0, 1)),
    (4, (0, 1, 2)),
    (5, (0, 1, 2, 3)),
    (10, tuple(range(8))),
    (2, (0, 1)),
]
SETS = {
    "q3": Q3,
    "q4": Q4,
    "q5": Q5,
    "half": HALF,
    "wide": WIDE,
    "ctl": CONTROL,
    "near": NEAR,
    "fix": FIX,
    "two": TWO,
    "ladder": LADDER,
    "rungs": RUNGS,
    "q34": Q3 + Q4,
    "all": Q3 + Q4 + Q5 + HALF + CONTROL + RUNGS + WIDE,
}
CAP = {(10, tuple(range(8))): 25.0}

# THE LATTICE


def period(q):
    return 2 * mp.pi / mp.log(q)


def height(q, ymax):
    per = period(q)
    y = mp.mpf(ymax)
    j0 = int(mp.floor(y / per))
    for j in (j0, j0 + 1):
        if abs(y - j * per) < mp.mpf("0.5"):
            y = (j + mp.mpf("0.5")) * per
    return y


def lattice(q, F, c):
    k = len(F)
    g1 = mp.mpf(sum(F))
    out = []
    for rec in c:
        out.append(
            {
                "lvl": 0,
                "j": rec["j"],
                "s0": rec["s0"],
                "r": rec["r"],
                "null": rec["null"],
            }
        )
        s1 = rec["s0"] - 1
        r1 = s1 * g1 * rec["r"] / (k * (q - 1))
        out.append(
            {
                "lvl": 1,
                "j": rec["j"],
                "s0": s1,
                "r": r1,
                "null": rec["null"] or abs(s1) < mp.mpf("1e-15") or g1 == 0,
            }
        )
    return out


def assign(z, lat, rho=RHO):
    for p in lat:
        if abs(z - p["s0"]) < rho:
            if p["null"]:
                if abs(z - p["s0"]) < mp.mpf("1e-7"):
                    return "cofactor", p
                continue
            return "tooth", p
    return "second", None


def classify(zs, lat, rho):
    kinds = [(z,) + assign(z, lat, rho) for z in sorted(zs, key=lambda w: mp.im(w))]
    sec = [z for z, kd, p in kinds if kd == "second"]
    tth = [z for z, kd, p in kinds if kd == "tooth"]
    cof = [z for z, kd, p in kinds if kd == "cofactor"]
    lv1 = [z for z, kd, p in kinds if kd == "tooth" and p["lvl"] == 1]
    return sec, tth, cof, lv1


# THE HUNT


def hunt(d, dp, f, x0, x1, res, seeds, fine=0):
    nx, ny = [(11, 7), (17, 11), (25, 15)][fine]
    extra = [3, 8, 14][fine]
    zs = list(seeds)
    short = []
    for a, b, w, mxx in res:
        want = int(mp.nint(abs(w)))
        got = [z for z in zs if a <= mp.im(z) < b]
        if len(got) >= want:
            continue
        for z0 in seek(f, x0, x1, a, b, want + extra, nx, ny):
            if len(got) >= want:
                break
            z, _ = polish(d, z0, (b - a) / 40, x1 - x0, dp)
            if z is None or mp.re(z) < x0 or mp.re(z) > x1:
                continue
            if mp.im(z) < a - mp.mpf("0.02") or mp.im(z) > b + mp.mpf("0.02"):
                continue
            if any(abs(z - y) < mp.mpf("1e-8") for y in zs):
                continue
            zs.append(z)
            if a <= mp.im(z) < b:
                got.append(z)
        if len(got) < want:
            short.append((a, b, want, len(got)))
    return zs, short


def sweep(q, F, ymax, show=True):
    t0 = time.time()
    y = height(q, ymax)
    a = alpha_of(q, F)
    per = period(q)
    x0, x1 = a + LO, a + HI
    c = comb(q, F, y, verbose=False)
    lat = lattice(q, F, c)
    d = build(q, F, dps=25, tol=mp.mpf(10) ** -10)
    dp = build(q, F, dps=32, tol=mp.mpf(10) ** -18)
    f = ZCache(d)
    res = strip(f, x0, x1, mp.mpf("0.02"), y, max(4, int(y / 2)), 20)
    tot = int(mp.nint(sum(r[2] for r in res)))
    mx = max(r[3] for r in res)
    seeds, discs = [], 0
    for rec in c:
        w = int(mp.nint(rec["w"]))
        if mp.mpf("0.02") + RHO < mp.im(rec["s0"]) < y - RHO:
            discs += w
        if w <= 0:
            continue
        have = (
            [rec["z"]]
            if rec["z"] is not None and abs(rec["z"] - rec["s0"]) < RHO
            else []
        )
        if rec["null"]:
            have.append(rec["s0"])
        if len(have) < w:
            have = teeth(f, dp, rec["s0"], RHO, w, False)
        for z in have:
            if (
                x0 < mp.re(z) < x1
                and mp.mpf("0.02") < mp.im(z) < y
                and not any(abs(z - w2) < mp.mpf("1e-8") for w2 in seeds)
            ):
                seeds.append(z)
    zs, short = hunt(d, dp, f, x0, x1, res, seeds)
    for step in (1, 2):
        if short or len(zs) < tot:
            zs, short = hunt(d, dp, f, x0, x1, res, zs, fine=step)
    sec, tth, cof, lv1 = classify(zs, lat, RHO)
    sec6 = classify(zs, lat, RHO2)[0]
    live = [p for p in lat if not p["null"]]
    us = (
        sorted(min(abs(z - p["s0"]) for p in live) for z in sec) if live and sec else []
    )
    out = {
        "q": q,
        "F": F,
        "y": y,
        "alpha": a,
        "per": per,
        "tot": tot,
        "mx": mx,
        "sec": sec,
        "tth": tth,
        "cof": cof,
        "lv1": lv1,
        "short": short,
        "us": us,
        "lat": lat,
        "sec_n": len(sec),
        "kq": mp.mpf(len(F)) / q,
        "discs": discs,
        "sec6": sec6,
        "y0": mp.mpf(ymax),
        "below": [z for z in sec if mp.im(z) < mp.mpf(ymax)],
    }
    if not show:
        return out
    line(
        "FAMILY",
        tag(q, F),
        "alpha",
        mp.nstr(a, 12),
        "alpha/2",
        mp.nstr(a / 2, 12),
        "k/q",
        mp.nstr(out["kq"], 8),
        "period",
        mp.nstr(per, 10),
        "T",
        mp.nstr(y, 10),
        "strip Re",
        mp.nstr(x0, 8),
        mp.nstr(x1, 8),
    )
    line(
        "  SPLIT winding",
        tot,
        "located",
        len(zs),
        "teeth",
        len(tth),
        "level one teeth",
        len(lv1),
        "cofactor only",
        len(cof),
        "N_2 at rho 0.45",
        len(sec),
        "N_2 at rho 0.6",
        len(sec6),
        "N_2 at rho 0.45 below Im",
        mp.nstr(mp.mpf(ymax), 6),
        len(out["below"]),
        "pole disc zeros",
        discs,
        "winding less discs",
        tot - discs,
        "max phase step",
        mp.nstr(mx, 4),
        "short boxes",
        short,
        "evals",
        f.n,
        "sec",
        round(time.time() - t0, 1),
    )
    line(
        "  RADIUS second family distance to the nearest live pole: min",
        "-" if not us else mp.nstr(us[0], 8),
        "in [0.45,0.6)",
        sum(1 for u in us if u < RHO2),
        "in [0.45,0.9)",
        sum(1 for u in us if u < mp.mpf("0.9")),
        "of",
        len(us),
        "first five",
        " ".join(mp.nstr(u, 8) for u in us[:5]),
    )
    r10 = [p for p in lat if p["lvl"] == 1 and p["j"] <= 1]
    line(
        "  LEVEL ONE residues derived from r_(0,j)",
        " ".join(
            "j "
            + str(p["j"])
            + " abs r "
            + mp.nstr(abs(p["r"]), 6)
            + " null "
            + str(p["null"])
            for p in r10
        ),
    )
    for z in sec:
        t = mp.im(z) / per
        line(
            "  SEC",
            tag(q, F),
            "Re",
            mp.nstr(mp.re(z), 12),
            "Im",
            mp.nstr(mp.im(z), 12),
            "Im/per",
            mp.nstr(t, 10),
            "frac",
            mp.nstr(t - mp.floor(t), 8),
            "Re-alpha/2",
            mp.nstr(mp.re(z) - a / 2, 10),
            "Re-alpha",
            mp.nstr(mp.re(z) - a, 10),
        )
    if sec:
        rm = max(sec, key=lambda z: mp.re(z))
        ra = max(zs, key=lambda z: mp.re(z))
        line(
            "  RIGHTMOST",
            tag(q, F),
            "second family Re",
            mp.nstr(mp.re(rm), 12),
            "at Im",
            mp.nstr(mp.im(rm), 12),
            "right of alpha",
            mp.re(rm) > a,
            "Mertens exponent of nu_F at least",
            mp.nstr(mp.re(rm), 10) if mp.re(rm) > a else "-",
            "| any zero Re",
            mp.nstr(mp.re(ra), 12),
            "at Im",
            mp.nstr(mp.im(ra), 12),
            "right of alpha",
            mp.re(ra) > a,
            "exponent at least",
            mp.nstr(mp.re(ra), 10) if mp.re(ra) > a else "-",
        )
    return out


# SYMMETRY


def symmetry(got):
    line(
        "SYMMETRY expectation on the full set c = 1/2 = alpha/2 and every zero self paired,",
        "on a design no functional equation and so failure at the third zero",
    )
    tally = [0, 0, 0, 0, 0, mp.mpf(0)]
    for key in got:
        o = got[key]
        sec = sorted(o["sec"], key=lambda z: mp.im(z))
        if len(sec) < 3:
            line("  ", tag(*key), "second family", len(sec), "too few to test")
            continue
        cF = (mp.re(sec[0]) + mp.re(sec[1])) / 2
        ims = sorted(mp.im(z) for z in sec)
        mgap = min(ims[i + 1] - ims[i] for i in range(len(ims) - 1))
        res = [mp.re(z) for z in sec]
        blo, bhi = min(res), max(res)
        ov = max(mp.mpf(0), min(bhi, cF + TOLS) - max(blo, cF - TOLS))
        pself = ov / (bhi - blo) if bhi > blo else mp.mpf(1)
        ok, bad, worst, first, selfp = 0, 0, mp.mpf(0), None, 0
        for z in sec[2:]:
            want = 2 * cF - mp.re(z)
            if abs(mp.re(z) - cF) < TOLS:
                ok += 1
                selfp += 1
                continue
            best = mp.inf
            for w in sec:
                if abs(w - z) < mp.mpf("1e-9"):
                    continue
                best = min(best, max(abs(mp.im(w) - mp.im(z)), abs(mp.re(w) - want)))
            if best < max(TOLT, TOLS):
                ok += 1
            else:
                bad += 1
                worst = max(worst, mp.mpf(0) if best == mp.inf else best)
                if first is None:
                    first = (z, want, best)
        a = o["alpha"]
        tally[0] += 1
        tally[1] += ok
        tally[2] += ok + bad
        tally[4] += selfp
        tally[5] += pself * (ok + bad)
        if first is not None and abs(first[0] - sec[2]) < mp.mpf("1e-9"):
            tally[3] += 1
        line(
            "  ",
            tag(*key),
            "c_F",
            mp.nstr(cF, 10),
            "from Im",
            mp.nstr(mp.im(sec[0]), 8),
            mp.nstr(mp.im(sec[1]), 8),
            "c_F - alpha/2",
            mp.nstr(cF - a / 2, 8),
            "c_F - alpha",
            mp.nstr(cF - a, 8),
            "c_F - 1/2",
            mp.nstr(cF - mp.mpf("0.5"), 8),
            "paired",
            str(ok) + "/" + str(ok + bad),
            "worst residual",
            "-" if first is None else mp.nstr(worst, 8),
            "first failure",
            "-"
            if first is None
            else mp.nstr(mp.re(first[0]), 8)
            + "+"
            + mp.nstr(mp.im(first[0]), 8)
            + "i wants Re "
            + mp.nstr(first[1], 8)
            + " nearest "
            + mp.nstr(first[2], 6),
            "self pairs",
            selfp,
            "reflection partners",
            ok - selfp,
            "chance self pair rate",
            mp.nstr(pself, 6),
            "expected",
            mp.nstr(pself * (ok + bad), 6),
            "min ordinate gap",
            mp.nstr(mgap, 8),
            "test tolerance",
            mp.nstr(TOLT, 3),
        )
    line(
        "  CHANCE every pair recorded is a self pair, a real part landing within",
        mp.nstr(TOLS, 3),
        "of c_F; the reflection branch needs two zeros within",
        mp.nstr(TOLT, 3),
        "in Im and the smallest ordinate gap in any design is far above it",
    )
    line(
        "  TALLY self pairs",
        tally[4],
        "reflection partners",
        tally[1] - tally[4],
        "expected self pairs by chance",
        mp.nstr(tally[5], 6),
        "of",
        tally[2],
        "observed",
        tally[1],
        "rate",
        mp.nstr(mp.mpf(tally[1]) / max(1, tally[2]), 6),
        "against chance",
        mp.nstr(tally[5] / max(1, tally[2]), 6),
    )
    line(
        "  TALLY designs tested",
        tally[0],
        "zeros paired",
        str(tally[1]) + "/" + str(tally[2]),
        "designs failing at the third zero",
        tally[3],
    )


def shadow(got):
    line(
        "SHADOW expectation derived from the full set, where the second family IS the zeros of",
        "zeta: if the family is continuous in F its ordinates converge as alpha -> 1",
    )
    line(
        "  the null for an unrelated ordinate set of the same density is a quarter of the mean",
        "gap between consecutive zeta ordinates in the range, and the full set reads exactly 0",
    )
    zz = [mp.im(mp.zetazero(n)) for n in range(1, 40)]
    rows = []
    for key in got:
        o = got[key]
        sec = sorted(o["sec"], key=lambda z: mp.im(z))
        if len(sec) < 3:
            continue
        band = [t for t in zz if mp.im(sec[0]) - 2 < t < mp.im(sec[-1]) + 2]
        if len(band) < 2:
            continue
        gap = (band[-1] - band[0]) / (len(band) - 1)
        near = [min(zz, key=lambda t: abs(mp.im(z) - t)) for z in sec]
        ds = [abs(mp.im(z) - t) for z, t in zip(sec, near)]
        used = len({mp.nstr(t, 12) for t in near})
        dr = [abs(mp.re(z) - mp.mpf("0.5")) for z in sec]
        rows.append(
            (
                o["alpha"],
                tag(*key),
                len(sec),
                sum(ds) / len(ds),
                max(ds),
                gap / 4,
                sum(dr) / len(dr),
                max(dr),
                used,
            )
        )
    for a, t, n, md, xd, null, mr, xr, used in sorted(rows):
        line(
            "  ",
            t,
            "alpha",
            mp.nstr(a, 10),
            "n",
            n,
            "mean dist to a zeta ordinate",
            mp.nstr(md, 8),
            "max",
            mp.nstr(xd, 8),
            "null gap/4",
            mp.nstr(null, 8),
            "ratio to null",
            mp.nstr(md / null, 8),
            "distinct ordinates",
            used,
            "| mean abs(Re - 1/2)",
            mp.nstr(mr, 8),
            "max",
            mp.nstr(xr, 8),
        )


# THE TRANSPORT READING


def mertens(got):
    line(
        "MERTENS a zero at Re rho > alpha forces the Mertens exponent of the design's own",
        "Mobius nu_F to be at least Re rho",
    )
    rows = []
    for key in got:
        o = got[key]
        zs = o["sec"] + o["tth"]
        if not zs:
            continue
        ra = max(zs, key=lambda z: mp.re(z))
        rs = max(o["sec"], key=lambda z: mp.re(z)) if o["sec"] else None
        rows.append((tag(*key), o["alpha"], ra, rs))
    for t, a, ra, rs in rows:
        line(
            "  ",
            t,
            "alpha",
            mp.nstr(a, 10),
            "rightmost zero Re",
            mp.nstr(mp.re(ra), 12),
            "at Im",
            mp.nstr(mp.im(ra), 12),
            "exponent bound",
            mp.nstr(mp.re(ra), 10) if mp.re(ra) > a else "none",
            "rightmost second family Re",
            "-" if rs is None else mp.nstr(mp.re(rs), 12),
            "second family bound",
            "none" if rs is None or mp.re(rs) <= a else mp.nstr(mp.re(rs), 10),
        )
    line(
        "  BOUNDS designs with a proved exponent bound",
        sum(1 for t, a, ra, rs in rows if mp.re(ra) > a),
        "of",
        len(rows),
        "and from the second family alone",
        sum(1 for t, a, ra, rs in rows if rs is not None and mp.re(rs) > a),
    )


# COUNT


def wcount(q, F, ymax):
    t0 = time.time()
    y = height(q, ymax)
    a = alpha_of(q, F)
    per = period(q)
    d = build(q, F, dps=25, tol=mp.mpf(10) ** -10)
    f = ZCache(d)
    res = strip(f, a + LO, a + HI, mp.mpf("0.02"), y, max(4, int(y / 10)), 20)
    tot = int(mp.nint(sum(r[2] for r in res)))
    mx = max(r[3] for r in res)
    tt, t6, j = 0, 0, 1
    while j * per + RHO2 < y:
        s0 = a + 2 * mp.pi * 1j * j / mp.log(q)
        tt += int(mp.nint(disc(f, s0, RHO)[0]))
        t6 += int(mp.nint(disc(f, s0, RHO2)[0]))
        j += 1
    w0 = int(mp.nint(disc(f, mp.mpc(a, 0), RHO)[0]))
    line(
        "  WCOUNT",
        tag(q, F),
        "T",
        mp.nstr(y, 10),
        "winding",
        tot,
        "teeth at rho 0.45",
        tt,
        "N_2 at rho 0.45",
        tot - tt,
        "teeth at rho 0.6",
        t6,
        "N_2 at rho 0.6",
        tot - t6,
        "j=0 disc",
        w0,
        "max phase step",
        mp.nstr(mx, 4),
        "evals",
        f.n,
        "sec",
        round(time.time() - t0, 1),
    )
    return y, tot - tt, tot - t6


def count(which, ymax):
    line(
        "COUNT expectation from the classical N(T) = (T/2 pi) log(T/2 pi e) + O(log T):",
        "c = 1/(2 pi) =",
        mp.nstr(1 / (2 * mp.pi), 10),
        "d = -(1 + log 2 pi)/(2 pi) =",
        mp.nstr(-(1 + mp.log(2 * mp.pi)) / (2 * mp.pi), 10),
    )
    line(
        "  N_2(T) is the winding of Z on the box less the zeros of Z in every pole disc,",
        "both by the argument principle; no zero is located",
    )
    rows = []
    for q, F in which:
        t1, n1, m1 = wcount(q, F, ymax)
        t2, n2, m2 = wcount(q, F, 2 * ymax)
        cf = (mp.mpf(n2) / t2 - mp.mpf(n1) / t1) / (mp.log(t2) - mp.log(t1))
        df = mp.mpf(n1) / t1 - cf * mp.log(t1)
        rows.append((q, F, t1, n1, t2, n2, cf, df, m1, m2))
    line("  READ c_F and d_F from the two heights, N_2(T) = c_F T log T + d_F T")
    for q, F, t1, n1, t2, n2, cf, df, m1, m2 in rows:
        line(
            "   ",
            tag(q, F),
            "alpha",
            mp.nstr(alpha_of(q, F), 10),
            "k",
            len(F),
            "T1",
            mp.nstr(t1, 8),
            "N_2 at rho 0.45",
            n1,
            "at rho 0.6",
            m1,
            "T2",
            mp.nstr(t2, 8),
            "N_2 at rho 0.45",
            n2,
            "at rho 0.6",
            m2,
            "c_F",
            mp.nstr(cf, 10),
            "d_F",
            mp.nstr(df, 10),
            "c_F 2 pi/log q",
            mp.nstr(cf * 2 * mp.pi / mp.log(q), 10),
        )
    groups = {}
    for q, F, t1, n1, t2, n2, cf, df, m1, m2 in rows:
        groups.setdefault(mp.nstr(alpha_of(q, F), 8), []).append((tag(q, F), cf))
    for a in sorted(groups):
        fam = groups[a]
        if len(fam) < 2:
            continue
        line(
            "   EQUAL ALPHA",
            a,
            "c_F spread",
            mp.nstr(max(x for t, x in fam) - min(x for t, x in fam), 8),
            " ".join(t + " " + mp.nstr(x, 8) for t, x in fam),
        )
    ks = {}
    for q, F, t1, n1, t2, n2, cf, df, m1, m2 in rows:
        ks.setdefault(len(F), []).append((tag(q, F), cf))
    for k in sorted(ks):
        fam = ks[k]
        if len(fam) < 2:
            continue
        line(
            "   EQUAL k",
            k,
            "c_F spread",
            mp.nstr(max(x for t, x in fam) - min(x for t, x in fam), 8),
            " ".join(t + " " + mp.nstr(x, 8) for t, x in fam),
        )


# THE FULL SET LIMIT


def limit(got):
    line(
        "LIMIT expectation: if RH is the alpha -> 1 face of MrlyMath the second family's",
        "Re s - alpha/2 contracts to 0 as k/q -> 1, spread 0 on the full set",
    )
    rows = []
    for key in got:
        o = got[key]
        sec = o["sec"]
        if not sec:
            continue
        a = o["alpha"]
        v = [mp.re(z) - a / 2 for z in sec]
        m = sum(v) / len(v)
        sd = mp.sqrt(sum((x - m) ** 2 for x in v) / len(v))
        rows.append(
            (o["kq"], key, a, len(v), m, sd, max(v) - min(v), max(abs(x) for x in v))
        )
    for kq, key, a, n, m, sd, spread, mab in sorted(rows):
        line(
            "  ",
            tag(*key),
            "k/q",
            mp.nstr(kq, 8),
            "alpha",
            mp.nstr(a, 10),
            "1-alpha",
            mp.nstr(1 - a, 10),
            "n",
            n,
            "mean Re-alpha/2",
            mp.nstr(m, 10),
            "sd",
            mp.nstr(sd, 10),
            "spread",
            mp.nstr(spread, 10),
            "max abs",
            mp.nstr(mab, 10),
            "max abs/(1-alpha)",
            "-" if a >= 1 else mp.nstr(mab / (1 - a), 8),
        )
    rat = [
        (mab / (1 - a), sd / (1 - a), abs(m) / (1 - a), tag(*key))
        for kq, key, a, n, m, sd, spread, mab in rows
        if a < 1
    ]
    if rat:
        line(
            "  BAND over the designs with alpha < 1: max abs/(1-alpha) in [",
            mp.nstr(min(x[0] for x in rat), 8),
            ",",
            mp.nstr(max(x[0] for x in rat), 8),
            "] low",
            min(rat)[3],
            "high",
            max(rat)[3],
            "; sd/(1-alpha) in [",
            mp.nstr(min(x[1] for x in rat), 8),
            ",",
            mp.nstr(max(x[1] for x in rat), 8),
            "]; abs mean/(1-alpha) in [",
            mp.nstr(min(x[2] for x in rat), 8),
            ",",
            mp.nstr(max(x[2] for x in rat), 8),
            "]",
        )
    full = [x for x in rows if x[2] >= 1]
    if full:
        line(
            "  FULL SET controls max abs",
            " ".join(tag(*x[1]) + " " + mp.nstr(x[7], 6) for x in full),
        )


# STUDY


def main():
    argv = sys.argv[1:]
    verb = argv[0] if argv else "family"
    ymax = float(argv[1]) if len(argv) > 1 else 40.0
    which = SETS[argv[2]] if len(argv) > 2 else SETS["all"]
    t0 = time.time()
    if verb == "count":
        count(which, ymax)
    else:
        line(
            "FAMILY the zeros of zeta_F left when every comb of the pole lattice is stripped;",
            "a zero within abs(u) <",
            mp.nstr(RHO, 3),
            "of a level zero or level one pole with nonvanishing residue is a tooth",
        )
        line(
            "BLIND the ladder does not reach tolerance in alpha - 1 < Re s < alpha - 0.92,",
            "so every count is a count on alpha - 0.92 < Re s < alpha + 3.02",
        )
        got = {}
        for q, F in which:
            got[(q, F)] = sweep(q, F, min(ymax, CAP.get((q, F), ymax)))
        if verb in ("symmetry", "tests"):
            symmetry(got)
        if verb in ("limit", "tests"):
            limit(got)
            shadow(got)
        if verb in ("family", "tests"):
            mertens(got)
    line("seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
