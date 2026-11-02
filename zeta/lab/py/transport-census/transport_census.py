import os
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "design-zeta"))
sys.path.insert(0, os.path.join(HERE, "..", "zeta-locus"))

import design_zeta as dz
import zeta_locus as zl

# DESIGNS

LOCUS = zl.Q3 + zl.Q4 + zl.Q5 + zl.HALF + zl.CONTROL + zl.WIDE
SHADOW = [(20, tuple(range(19))), (50, tuple(range(49)))]
RUNGS = [(5, (0, 1, 2, 3)), (10, tuple(range(8)))]
FULL = [(2, (0, 1)), (3, (0, 1, 2)), (4, (0, 1, 2, 3))]
SETS = {
    "locus": LOCUS,
    "shadow": SHADOW,
    "ctl": FULL,
    "rungs": RUNGS,
    "b20": SHADOW[:1],
    "b50": SHADOW[1:],
    "all": LOCUS + RUNGS,
}

HALF = mp.mpf("5e-5")
EDGE = mp.mpf(10) ** -6
YMAX = mp.mpf(40)


def line(*a):
    print(" ".join(str(x) for x in a))


def down(x, n=7):
    return mp.nstr(mp.floor(x * mp.mpf(10) ** n) / mp.mpf(10) ** n, n + 2)


def build(q, F):
    return zl.build(q, F, dps=25, tol=mp.mpf(10) ** -10)


# THE ZERO FREE EDGE


def edge(d):
    amin = min(a for a in d.F if a != 0)
    s = d.alpha + mp.mpf("0.05")
    while s < mp.mpf(14):
        v, e = d.zeta(s)
        if mp.power(amin, s) * (v + e) < 2:
            return s, v, e
        s += mp.mpf("0.05")
    return None, None, None


# THE POLES


def poledist(d, x0, x1, y0, y1):
    per = 2 * mp.pi / mp.log(d.q)
    best = mp.inf
    for m in range(3):
        xp = d.alpha - m
        for j in range(int(mp.floor(y0 / per)) - 1, int(mp.ceil(y1 / per)) + 2):
            yp = j * per
            dx = max(x0 - xp, xp - x1, mp.mpf(0))
            dy = max(y0 - yp, yp - y1, mp.mpf(0))
            best = min(best, mp.sqrt(dx * dx + dy * dy))
    return best


# THE SWEEP


def sweep(q, F, ymax, nsub):
    d = build(q, F)
    dp = zl.build(q, F, dps=32, tol=mp.mpf(10) ** -18)
    s1, v1, e1 = edge(d)
    x0, x1 = d.alpha + EDGE, s1
    f = dz.Cache(d, cof=True)
    t0 = time.time()
    res = dz.strip(f, x0, x1, mp.mpf("0.02"), ymax, nsub, 20)
    tot = sum(r[2] for r in res)
    mx = max(r[3] for r in res)
    zs = []
    for a, b, w, m in res:
        want = int(mp.nint(abs(w)))
        if want == 0:
            continue
        hit = 0
        for z0 in zl.seek(f, x0, x1, a, b, 2 * want + 8, 15, 11):
            if hit >= want:
                break
            z, _ = zl.polish(d, z0, (b - a) / 40, x1 - x0, dp)
            if z is None or mp.re(z) < x0 or mp.re(z) > x1:
                continue
            if mp.im(z) < a - mp.mpf("0.5") or mp.im(z) > b + mp.mpf("0.5"):
                continue
            if any(abs(z - y) < mp.mpf("1e-9") for y in zs):
                continue
            hit += 1
            zs.append(z)
    return {
        "d": d,
        "s1": s1,
        "v1": v1,
        "e1": e1,
        "x0": x0,
        "x1": x1,
        "res": res,
        "n": int(mp.nint(tot)),
        "wind": tot,
        "mx": mx,
        "zs": zs,
        "evals": f.n,
        "sec": time.time() - t0,
    }


# THE BOX


def certify(q, F, z, half=HALF, n=16):
    d = build(q, F)
    f = dz.Cache(d)
    x0, x1 = mp.re(z) - half, mp.re(z) + half
    y0, y1 = mp.im(z) - half, mp.im(z) + half
    w, mx = dz.box_phase(f, x0, x1, y0, y1, n)
    lo = min(abs(v) for v in f.m.values())
    return {
        "x0": x0,
        "x1": x1,
        "y0": y0,
        "y1": y1,
        "w": w,
        "mx": mx,
        "lo": lo,
        "e": f.emax,
        "pole": poledist(d, x0, x1, y0, y1),
        "alpha": d.alpha,
    }


VERBS = {}

# CENSUS


def census(which=None, ymax=YMAX):
    which = which or LOCUS
    line("TRANSPORT CENSUS: zeros of zeta_F right of alpha, and the Mertens exponent")
    line("of the design's own Mobius nu_F they force through the transport theorem")
    line(
        "sigma1 is a proved zero free edge: a_min^sigma zeta_F(sigma) < 2 there, so the"
    )
    line("census over alpha < Re s < sigma1 is complete in Re and cut only in height;")
    line(
        "the count is the winding and is exact, the rightmost is the rightmost LOCATED"
    )
    got = {}
    for q, F in which:
        r = sweep(q, F, ymax, int(ymax / 2))
        d = r["d"]
        zs = sorted(r["zs"], key=lambda z: -mp.re(z))
        line(
            " ",
            zl.tag(q, F),
            "alpha",
            mp.nstr(d.alpha, 10),
            "k/q",
            mp.nstr(mp.mpf(d.k) / q, 8),
            "sigma1",
            mp.nstr(r["s1"], 6),
            "a_min^s zeta_F",
            mp.nstr(mp.power(min(a for a in d.F if a != 0), r["s1"]) * r["v1"], 8),
            "bound",
            mp.nstr(r["e1"], 3),
        )
        line(
            "    strip Re",
            mp.nstr(r["x0"], 10),
            mp.nstr(r["x1"], 8),
            "Im to",
            mp.nstr(ymax, 6),
            "Z winding",
            mp.nstr(r["wind"], 8),
            "zeros",
            r["n"],
            "located",
            len(zs),
            "complete",
            len(zs) == r["n"],
            "max phase step",
            mp.nstr(r["mx"], 4),
            "evals",
            r["evals"],
            "sec",
            round(r["sec"], 1),
        )
        if not zs:
            line(
                "    NO ZERO RIGHT OF alpha to height",
                mp.nstr(ymax, 6),
                "- no transport bound on theta(nu_F) from this census",
            )
            got[(q, F)] = None
            continue
        for z in zs:
            line(
                "    ZERO",
                zl.tag(q, F),
                "Re",
                mp.nstr(mp.re(z), 12),
                "Im",
                mp.nstr(mp.im(z), 12),
                "Re-alpha",
                mp.nstr(mp.re(z) - d.alpha, 10),
            )
        z = zs[0]
        b = certify(q, F, z)
        line(
            "    BOX Re [",
            mp.nstr(b["x0"], 12),
            ",",
            mp.nstr(b["x1"], 12),
            "] Im [",
            mp.nstr(b["y0"], 12),
            ",",
            mp.nstr(b["y1"], 12),
            "] winding",
            mp.nstr(b["w"], 8),
            "step",
            mp.nstr(b["mx"], 4),
            "min |zeta_F|",
            mp.nstr(b["lo"], 4),
            "engine bound",
            mp.nstr(b["e"], 4),
            "pole distance",
            mp.nstr(b["pole"], 6),
        )
        ok = abs(b["w"] - 1) < mp.mpf("0.05") and b["lo"] > b["e"]
        line(
            "    THETA",
            zl.tag(q, F),
            "theta(nu_F) >=",
            down(b["x0"]),
            "alpha",
            mp.nstr(d.alpha, 10),
            "gain",
            down(b["x0"] - d.alpha),
            "above 1",
            b["x0"] > 1,
            "certified",
            ok,
        )
        got[(q, F)] = (z, b, len(zs), r["n"])
    return got


VERBS["census"] = census

# LAW


def nullbox(q, F, j, half=HALF):
    d = build(q, F)
    s0 = d.alpha + 2 * mp.pi * 1j * j / mp.log(q)
    out = []
    for cof in (True, False):
        f = dz.Cache(d, cof=cof)
        w, mx = dz.box_phase(
            f,
            mp.re(s0) - half,
            mp.re(s0) + half,
            mp.im(s0) - half,
            mp.im(s0) + half,
            16,
        )
        out.append((w, mx, min(abs(v) for v in f.m.values()), f.emax))
    return s0, out


LANDED = [
    (3, (0, 1), "0.72074", "0.72084", "28.60563", "28.60573", 1),
    (10, tuple(range(9)), "1.00150", "1.00168", "2.73915", "2.73925", 1),
    (10, tuple(range(9)), "0.99900", "1.00050", "2.73810", "2.74030", 0),
]


def repro():
    line("REPRODUCTION of the two landed boxes of lab/py/mrly-pairing and its control,")
    line("same edges, same argument principle, evaluated here")
    for q, F, x0, x1, y0, y1, want in LANDED:
        d = build(q, F)
        f = dz.Cache(d)
        w, mx = dz.box_phase(f, mp.mpf(x0), mp.mpf(x1), mp.mpf(y0), mp.mpf(y1), 16)
        lo = min(abs(v) for v in f.m.values())
        line(
            " ",
            zl.tag(q, F),
            "Re [",
            x0,
            ",",
            x1,
            "] Im [",
            y0,
            ",",
            y1,
            "] winding",
            mp.nstr(w, 8),
            "expected",
            want,
            "step",
            mp.nstr(mx, 4),
            "min |zeta_F|",
            mp.nstr(lo, 4),
            "engine bound",
            mp.nstr(f.emax, 4),
            "pole distance",
            mp.nstr(poledist(d, mp.mpf(x0), mp.mpf(x1), mp.mpf(y0), mp.mpf(y1)), 6),
            "reproduced",
            abs(w - want) < mp.mpf("0.05"),
        )


def law(which=None, ymax=YMAX):
    which = which or LOCUS
    got = census(which, ymax)
    line("LAW: the columns, no fit")
    line(
        "  DESIGN alpha k/q sigma1 zeros-right-of-alpha rightmost-Re gain=Re-alpha above-1 1-in-F"
    )
    rows = []
    for q, F in which:
        g = got[(q, F)]
        d = build(q, F)
        one = 1 in F
        if g is None:
            rows.append((q, F, d.alpha, mp.mpf(d.k) / q, 0, None, None, one))
            continue
        z, b, _, ntot = g
        rows.append(
            (q, F, d.alpha, mp.mpf(d.k) / q, ntot, mp.re(z), mp.re(z) - d.alpha, one)
        )
    for r in sorted(rows, key=lambda t: (t[3], t[2])):
        q, F, a, kq, n, re, gain, one = r
        line(
            "  ROW",
            zl.tag(q, F),
            "alpha",
            mp.nstr(a, 10),
            "k/q",
            mp.nstr(kq, 8),
            "zeros",
            n,
            "rightmost located",
            "-" if re is None else mp.nstr(re, 12),
            "gain",
            "-" if gain is None else mp.nstr(gain, 10),
            "above 1",
            "-" if re is None else (re > 1),
            "1 in F",
            one,
        )
    line("  BY k/q, then by alpha: the gain column above is printed in k/q order and")
    line("  again in alpha order below; equal-key pairs with unequal gain kill a law")
    for r in sorted(rows, key=lambda t: (t[2], t[3])):
        q, F, a, kq, n, re, gain, one = r
        line(
            "  ALPHAORDER",
            zl.tag(q, F),
            "alpha",
            mp.nstr(a, 10),
            "k/q",
            mp.nstr(kq, 8),
            "gain",
            "-" if gain is None else mp.nstr(gain, 10),
        )
    keys = {}
    for r in rows:
        keys.setdefault((mp.nstr(r[2], 10), mp.nstr(r[3], 8)), []).append(r)
    for key, fam in sorted(keys.items()):
        if len(fam) < 2:
            continue
        gs = [r[6] for r in fam if r[6] is not None]
        line(
            "  TIE alpha",
            key[0],
            "k/q",
            key[1],
            "designs",
            " ".join(zl.tag(r[0], r[1]) for r in fam),
            "gains",
            " ".join("-" if r[6] is None else mp.nstr(r[6], 8) for r in fam),
            "spread",
            "-" if len(gs) < 2 else mp.nstr(max(gs) - min(gs), 8),
        )
    ab = [r for r in rows if r[5] is not None and r[5] > 1]
    line("  ABOVE ONE: designs whose rightmost censused zero has Re rho > 1, so the")
    line(
        "  Mertens of nu_F outruns x itself:",
        " ".join(zl.tag(r[0], r[1]) for r in ab),
        "| their k/q",
        " ".join(mp.nstr(r[3], 8) for r in ab),
        "| their alpha",
        " ".join(mp.nstr(r[2], 10) for r in ab),
    )
    line("  COROLLARY: 1 in F and a certified zero rho with Re rho > alpha give")
    line("  sigma_c(N_F) >= Re rho > alpha >= alpha/2, so sum_(n <= x) nu_F(n) is not")
    line(
        "  O(x^(alpha/2 + eps)) and not O(x^(alpha - eps)): no square-root-shaped bound"
    )
    for r in sorted(rows, key=lambda t: t[2]):
        q, F, a, kq, n, re, gain, one = r
        if re is None:
            line(
                "  NOBOUND",
                zl.tag(q, F),
                "no zero right of alpha to height",
                mp.nstr(ymax, 6),
            )
        elif not one:
            line(
                "  NOINVERSE",
                zl.tag(q, F),
                "1 not in F, so 1 not in S_F and the",
                "Dirichlet inverse nu_F does not exist; the zero stands, the",
                "transport does not",
            )
        else:
            b = got[(q, F)][1]
            line(
                "  BOUND",
                zl.tag(q, F),
                "theta(nu_F) >=",
                down(b["x0"]),
                ">",
                "alpha",
                mp.nstr(a, 10),
                ">= alpha/2",
                mp.nstr(a / 2, 10),
                "square-root shape refuted for nu_F",
                True,
            )
    line("  CONTROLS: the cofactor only teeth of the full sets at Re s = alpha = 1")
    for q, F in FULL:
        s0, out = nullbox(q, F, 1)
        (wz, mz, lz, _), (wf, mf, lf, _) = out
        line(
            "  NULLBOX",
            zl.tag(q, F),
            "pole",
            mp.nstr(s0, 12),
            "winding of Z",
            mp.nstr(wz, 8),
            "winding of zeta_F",
            mp.nstr(wf, 8),
            "min |Z|",
            mp.nstr(lz, 4),
            "min |zeta_F|",
            mp.nstr(lf, 4),
            "steps",
            mp.nstr(mz, 4),
            mp.nstr(mf, 4),
        )
    repro()


VERBS["law"] = law

# THE CROSSING


def root(which=None, ymax=YMAX):
    line(
        "THE CROSSING: the root of a_min^sigma zeta_F(sigma) = 2, bisected against the ladder bound, beside the grid edge sigma_1"
    )
    for q, F in which or LOCUS:
        d = build(q, F)
        amin = min(a for a in d.F if a != 0)
        s1, v1, e1 = edge(d)
        lo, hi = s1 - mp.mpf("0.05"), s1
        v, e = d.zeta(lo)
        assert mp.power(amin, lo) * (v - e) > 2, (q, F, lo)
        while hi - lo > mp.mpf(10) ** -12:
            mid = (lo + hi) / 2
            v, e = d.zeta(mid)
            g = mp.power(amin, mid)
            if g * (v - e) > 2:
                lo = mid
            elif g * (v + e) < 2:
                hi = mid
            else:
                break
        up = mp.ceil(hi * 10**5) / 10**5
        line(
            " ",
            zl.tag(q, F),
            "alpha",
            mp.nstr(d.alpha, 10),
            "a_min",
            amin,
            "grid sigma_1",
            mp.nstr(s1, 12),
            "a_min^s zeta_F there",
            mp.nstr(mp.power(amin, s1) * v1, 10),
            "bound",
            mp.nstr(e1, 3),
            "root in",
            mp.nstr(lo, 12),
            mp.nstr(hi, 12),
            "rounded up",
            mp.nstr(up, 7),
        )


VERBS["root"] = root


def main():
    argv = sys.argv[1:]
    verb = argv[0] if argv else "census"
    ymax = mp.mpf(argv[1]) if len(argv) > 1 else YMAX
    which = SETS[argv[2]] if len(argv) > 2 else LOCUS
    t0 = time.time()
    VERBS[verb](which, ymax)
    line("seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
