import os
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "design-zeta"))

from design_zeta import Design, phase, strip

# COFACTOR


def build(q, F, dps=45, tol=None, P=None):
    if tol is None:
        tol = mp.mpf(10) ** -22
    if P is None and len(F) == 1:
        P = 12
    return Design(q, F, P=P, dps=dps, tol=tol)


def ztune(d, s):
    key = int(abs(mp.im(s)) / 8)
    W, L = getattr(d, "wl", {}).get(key, (mp.mpf(6), 8))
    v, e = d.ladder(s, W, L, False)
    while e >= d.tol and W <= 54:
        W += 6
        L += 6
        v, e = d.ladder(s, W, L, False)
    if e < d.tol:
        if not hasattr(d, "wl"):
            d.wl = {}
        d.wl[key] = (max(mp.mpf(6), W - 6), max(8, L - 6))
    return v, e


def zwork(d, s):
    a, e = ztune(d, s)
    return (1 - d.k * mp.power(d.q, -s)) * d.poly_lo(s) + a, e


def cofactor(d, s):
    old = mp.mp.dps
    mp.mp.dps = d.dps + 30 + int(abs(mp.im(s)) / 4)
    v, e = zwork(d, mp.mpmathify(s))
    mp.mp.dps = old
    return +v, +e


class ZCache:
    def __init__(self, d):
        self.d = d
        self.m = {}
        self.n = 0

    def __call__(self, z):
        key = (mp.nstr(mp.re(z), 22), mp.nstr(mp.im(z), 22))
        if key not in self.m:
            self.m[key] = cofactor(self.d, z)[0]
            self.n += 1
        return self.m[key]


def laurent(d, j):
    old = mp.mp.dps
    mp.mp.dps = d.dps + 40 + int(abs(j) * 2)
    L = mp.log(d.q)
    s0 = mp.log(d.k) / L + 2 * mp.pi * 1j * j / L
    h = mp.mpf(10) ** -5
    z0, e0 = zwork(d, s0)
    zp = zwork(d, s0 + h)[0]
    zm = zwork(d, s0 - h)[0]
    z1 = (zp - zm) / (2 * h)
    z2 = (zp + zm - 2 * z0) / (2 * h * h)
    r = z0 / L
    c0 = z1 / L + z0 / 2
    c1 = z2 / L + z1 / 2 + z0 * L / 12
    mp.mp.dps = old
    return +s0, +r, +c0, +c1, +e0


def predict(r, c0, c1):
    u1 = -r / c0
    if c1 == 0:
        return u1, u1
    disc = mp.sqrt(c0 * c0 - 4 * c1 * r)
    roots = [(-c0 + disc) / (2 * c1), (-c0 - disc) / (2 * c1)]
    u2 = min(roots, key=lambda z: abs(z - u1))
    return u1, u2


def disc(f, s0, rho, n=32):
    t, mx = phase(f, lambda u: s0 + rho * mp.exp(2 * mp.pi * 1j * u), n, budget=280)
    return t / (2 * mp.pi), mx


def ring(f, s0, rho, nr=5, na=12):
    pts = [(abs(f(s0)), s0)]
    for i in range(1, nr + 1):
        rr = rho * mp.mpf(i) / (nr + 1)
        for m in range(na):
            pts.append(
                (abs(f(s0 + rr * mp.exp(2 * mp.pi * 1j * mp.mpf(m) / na))), i, m)
            )
    pts.sort(key=lambda t: t[0])
    out = []
    for t in pts:
        if len(t) == 2:
            out.append(t[1])
        else:
            rr = rho * mp.mpf(t[1]) / (nr + 1)
            out.append(s0 + rr * mp.exp(2 * mp.pi * 1j * mp.mpf(t[2]) / na))
    return out


def teeth(f, dp, s0, rho, want, skip):
    out = []
    for z0 in ring(f, s0, rho)[: want + 6]:
        if len(out) >= want:
            break
        z, _ = polish(f.d, z0, rho / 40, rho, dp)
        if z is None or abs(z - s0) > rho:
            continue
        if skip and abs(z - s0) < mp.mpf("1e-9"):
            continue
        if any(abs(z - y) < mp.mpf("1e-9") for y in out):
            continue
        out.append(z)
    return out


def seek(f, x0, x1, y0, y1, take, nx=11, ny=7):
    pts = [
        (
            abs(
                f(mp.mpc(x0 + (x1 - x0) * i / (nx - 1), y0 + (y1 - y0) * m / (ny - 1)))
            ),
            i,
            m,
        )
        for i in range(nx)
        for m in range(ny)
    ]
    pts.sort()
    return [
        mp.mpc(x0 + (x1 - x0) * i / (nx - 1), y0 + (y1 - y0) * m / (ny - 1))
        for v, i, m in pts[:take]
    ]


def muller(d, z0, w, span, tol, steps):
    def f(z):
        if abs(z - z0) > span:
            raise ValueError
        return cofactor(d, z)[0]

    try:
        r = mp.findroot(
            f, [z0 - w, z0, z0 + w * 1j], solver="muller", tol=tol, maxsteps=steps
        )
    except Exception:
        return None
    return r if abs(r - z0) <= span else None


def polish(d, z0, w, span, dp=None):
    r = muller(d, z0, w, span, mp.mpf(10) ** -18, 40)
    if r is None:
        return None, None
    if dp is not None:
        r2 = muller(dp, r, mp.mpf(10) ** -9, mp.mpf("1e-4"), mp.mpf(10) ** -20, 15)
        if r2 is not None:
            r = r2
            return r, abs(cofactor(dp, r)[0])
    return r, abs(cofactor(d, r)[0])


# SWEEP


def line(*a):
    print(" ".join(str(x) for x in a))


def tag(q, F):
    return "q" + str(q) + "F" + "".join(str(a) for a in F)


def alpha_of(q, F):
    return mp.log(len(F)) / mp.log(q)


def row(q, F, alpha, per, z, kind, pred):
    t = mp.im(z) / per
    line(
        "  ROW",
        tag(q, F),
        "Re",
        mp.nstr(mp.re(z), 12),
        "Im",
        mp.nstr(mp.im(z), 12),
        "Im/per",
        mp.nstr(t, 10),
        "frac",
        mp.nstr(t - mp.floor(t), 8),
        "alpha",
        mp.nstr(alpha, 10),
        "k/q",
        mp.nstr(mp.mpf(len(F)) / q, 8),
        "Re-alpha",
        mp.nstr(mp.re(z) - alpha, 10),
        kind,
        pred,
    )


# STUDY

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
COMB_TOL = mp.mpf(10) ** -24


def comb(q, F, ymax, dps=40, tol=COMB_TOL, verbose=True):
    d = build(q, F, dps=dps, tol=tol)
    dp = build(q, F, dps=32, tol=mp.mpf(10) ** -18)
    fw = ZCache(build(q, F, dps=25, tol=mp.mpf(10) ** -10))
    rho = mp.mpf("0.45")
    L = mp.log(q)
    per = 2 * mp.pi / L
    alpha = alpha_of(q, F)
    if verbose:
        line(
            " ",
            tag(q, F),
            "alpha",
            mp.nstr(alpha, 12),
            "period",
            mp.nstr(per, 12),
            "k/q",
            mp.nstr(mp.mpf(len(F)) / q, 8),
            "disc radius",
            mp.nstr(rho, 4),
        )
    out = []
    for j in range(int(mp.ceil(mp.mpf(ymax) / per)) + 1):
        s0, r, c0, c1, e = laurent(d, j)
        u1, u2 = predict(r, c0, c1)
        w, wmx = disc(fw, s0, rho)
        null = abs(r) < 1000 * e
        nz = int(mp.nint(w)) - (1 if null else 0)
        zl = teeth(fw, dp, s0, rho, nz, null) if nz > 0 else []
        z = min(zl, key=lambda y: abs(y - s0)) if zl else None
        rec = {
            "j": j,
            "s0": s0,
            "r": r,
            "R": c0,
            "u1": u1,
            "u2": u2,
            "z": z,
            "w": w,
            "wmx": wmx,
            "null": null,
            "nz": nz,
            "all": zl,
            "v": None if z is None else abs(fw(z)),
        }
        out.append(rec)
        if not verbose:
            continue
        msg = [
            "    j",
            j,
            "r",
            mp.nstr(r, 12),
            "R",
            mp.nstr(c0, 12),
            "u1",
            mp.nstr(u1, 10),
            "u2-u1",
            mp.nstr(abs(u2 - u1), 8),
            "ladder bound",
            mp.nstr(e, 3),
            "Z zeros in disc",
            mp.nstr(w, 8),
            "step",
            mp.nstr(wmx, 4),
            "residue null",
            null,
            "zeta_F zeros in disc",
            nz,
            "located",
            len(zl),
        ]
        if z is None:
            msg += ["tooth", "none"]
        else:
            m1 = abs(z - s0 - u1)
            m2 = abs(z - s0 - u2)
            msg += [
                "tooth",
                mp.nstr(z, 14),
                "|u|",
                mp.nstr(abs(z - s0), 8),
                "|Z|",
                mp.nstr(rec["v"], 3),
                "miss1",
                mp.nstr(m1, 6),
                "miss2",
                mp.nstr(m2, 6),
                "miss2/miss1",
                mp.nstr(m2 / max(m1, mp.mpf(10) ** -40), 5),
            ]
        line(*msg)
    return out


def census(q, F, ymax, nsub, lo=-0.92, hi=3.02, n=20, cmb=None, deep=False):
    d = build(q, F, dps=25, tol=mp.mpf(10) ** -10)
    dp = build(q, F, dps=32, tol=mp.mpf(10) ** -18)
    L = mp.log(q)
    per = 2 * mp.pi / L
    alpha = alpha_of(q, F)
    x0, x1 = alpha + mp.mpf(lo), alpha + mp.mpf(hi)
    f = ZCache(d)
    t0 = time.time()
    res = strip(f, x0, x1, mp.mpf("0.02"), mp.mpf(ymax), nsub, n)
    tot = sum(r[2] for r in res)
    mx = max(r[3] for r in res)
    bx = [r for r in res if abs(r[2]) > mp.mpf("0.2")]
    nulls = [
        mp.im(c["s0"])
        for c in (cmb or [])
        if c["null"] and c["j"] > 0 and mp.im(c["s0"]) <= ymax
    ]
    line(
        " ",
        tag(q, F),
        "alpha",
        mp.nstr(alpha, 10),
        "strip Re",
        mp.nstr(x0, 8),
        mp.nstr(x1, 8),
        "Im to",
        ymax,
        "Z winding",
        mp.nstr(tot, 8),
        "null teeth",
        len(nulls),
        "zeta_F zeros",
        mp.nstr(tot - len(nulls), 8),
        "boxes",
        len(bx),
        "max phase step",
        mp.nstr(mx, 4),
        "evals",
        f.n,
        "sec",
        round(time.time() - t0, 1),
    )
    cum, marks = mp.mpf(0), []
    for a, b, w, mxx in res:
        cum += w
        marks.append((b, cum - len([y for y in nulls if y <= b])))
    line(
        "    N_F(T)",
        " ".join(
            "T " + mp.nstr(T, 6) + " N " + mp.nstr(c, 6)
            for T, c in marks[nsub // 4 - 1 :: max(1, nsub // 4)]
        ),
        "per period",
        mp.nstr((tot - len(nulls)) * per / mp.mpf(ymax), 8),
    )
    zs = [
        (c["z"], "comb j " + str(c["j"]))
        for c in (cmb or [])
        if c["z"] is not None and mp.im(c["z"]) <= ymax + 1
    ]
    for c in cmb or []:
        for y in c["all"]:
            if y is c["z"] or mp.im(y) > ymax + 1:
                continue
            if not any(abs(y - t[0]) < mp.mpf("1e-9") for t in zs):
                zs.append((y, "comb j " + str(c["j"]) + " second tooth"))
    for a, b, w, mxx in bx if deep else []:
        want = int(mp.nint(abs(w)))
        hit = 0
        for z0 in seek(f, x0, x1, a, b, want + 3):
            if hit >= want:
                break
            z, _ = polish(d, z0, (b - a) / 40, x1 - x0, dp)
            if z is None or mp.re(z) < x0 or mp.re(z) > x1:
                continue
            if mp.im(z) < a - 1 or mp.im(z) > b + 1:
                continue
            hit += 1
            if any(abs(z - y[0]) < mp.mpf("1e-9") for y in zs):
                continue
            if any(
                c["null"] and abs(z - c["s0"]) < mp.mpf("1e-9") for c in (cmb or [])
            ):
                continue
            zs.append((z, "second"))
    line(
        "    located",
        len(zs),
        "of zeta_F zeros",
        mp.nstr(tot - len(nulls), 8),
        "complete",
        len(zs) == int(mp.nint(tot)) - len(nulls),
        "deep",
        deep,
        "teeth",
        sum(1 for z, k in zs if k.startswith("comb")),
    )
    for z, kind in sorted(zs, key=lambda y: mp.im(y[0])):
        pred = "-"
        if cmb and kind.startswith("comb"):
            best = min(cmb, key=lambda c: abs(z - c["s0"]))
            pred = (
                mp.nstr(best["s0"] + best["u1"], 12)
                + " miss "
                + mp.nstr(abs(z - best["s0"] - best["u1"]), 6)
            )
        row(q, F, alpha, per, z, kind, pred)
    return zs, tot, res


def laws(got):
    line("LAWS")
    groups = {}
    for q, F in got:
        groups.setdefault(mp.nstr(alpha_of(q, F), 10), []).append((q, F))
    for a in sorted(groups):
        fam = groups[a]
        if len(fam) < 2:
            continue
        line("  EQUAL ALPHA", a, "designs", " | ".join(tag(q, F) for q, F in fam))
        for q, F in fam:
            line(
                "    ",
                tag(q, F),
                "k/q",
                mp.nstr(mp.mpf(len(F)) / q, 8),
                "zeros",
                len(got[(q, F)]),
                "complete",
                got[(q, F)][0][2] if got[(q, F)] else "-",
            )
        for i in range(len(fam)):
            for m in range(i + 1, len(fam)):
                A = [z for z, k, c in got[fam[i]]]
                B = [z for z, k, c in got[fam[m]]]
                if not A or not B:
                    continue
                pairs = [
                    (abs(mp.im(z) - mp.im(y)), abs(mp.re(z) - mp.re(y)), z, y)
                    for z in A
                    for y in B
                ]
                shared = min(abs(z - y) for z in A for y in B)
                pairs.sort()
                di, dr, z, _ = pairs[0]
                near = [p for p in pairs if p[0] < mp.mpf("0.02")]
                worst = max(near, key=lambda p: p[1]) if near else None
                line(
                    "    ",
                    tag(*fam[i]),
                    "vs",
                    tag(*fam[m]),
                    "closest in Im: dIm",
                    mp.nstr(di, 6),
                    "dRe",
                    mp.nstr(dr, 6),
                    "at Im",
                    mp.nstr(mp.im(z), 10),
                    "| shared zero",
                    mp.nstr(shared, 6),
                    "| worst dRe at dIm<0.02",
                    "none"
                    if worst is None
                    else mp.nstr(worst[1], 6)
                    + " at Im "
                    + mp.nstr(mp.im(worst[2]), 10),
                )
    line("  SECOND FAMILY, the zeros left when every tooth is stripped")
    for q, F in got:
        sec = [mp.re(z) for z, k, c in got[(q, F)] if k == "second"]
        a2 = alpha_of(q, F) / 2
        line(
            "    ",
            tag(q, F),
            "alpha/2",
            mp.nstr(a2, 10),
            "n",
            len(sec),
            "Re range",
            "-" if not sec else mp.nstr(min(sec), 8),
            "-" if not sec else mp.nstr(max(sec), 8),
            "complete",
            got[(q, F)][0][2] if got[(q, F)] else "-",
        )
    line("  COMB in frac(Im s log q/2 pi): worst Re gap between zeros of equal frac")
    for q, F in got:
        per = 2 * mp.pi / mp.log(q)
        fr = sorted(((mp.im(z) / per) % 1, mp.re(z)) for z, k, c in got[(q, F)])
        if len(fr) < 2:
            continue
        pairs = [
            abs(fr[i][1] - fr[i + 1][1])
            for i in range(len(fr) - 1)
            if fr[i + 1][0] - fr[i][0] < mp.mpf("0.03")
        ]
        line(
            "    ",
            tag(q, F),
            "n",
            len(fr),
            "Re range",
            mp.nstr(min(x for f, x in fr), 8),
            mp.nstr(max(x for f, x in fr), 8),
            "worst gap at equal frac",
            "none" if not pairs else mp.nstr(max(pairs), 8),
        )


# ROUCHE


def explicit(d, s0, M=60):
    L = mp.log(d.q)
    dm = [mp.mpf(0)] * (M + 1)
    em = [mp.mpf(0)] * (M + 1)
    for arr, tgt in ((d.loglo, dm), (d.logmid, em)):
        for t in arr:
            v = mp.exp(-s0 * t)
            for m in range(M + 1):
                tgt[m] += v
                v = v * (-t) / (m + 1)
    c = [mp.mpf(0)] + [-((-L) ** m) / mp.factorial(m) for m in range(1, M + 1)]
    P = [
        em[m] + mp.fsum([c[i] * dm[m - i] for i in range(1, m + 1)])
        for m in range(M + 1)
    ]
    return P, c


def ptail(d, s0, rho, P, M=60):
    sig = mp.re(s0)
    L = mp.log(d.q)
    h = M // 2 + 1
    sd = mp.fsum([mp.exp((rho - sig) * t) for t in d.loglo])
    cd = mp.fsum(
        [mp.exp((rho - sig) * t) * (rho * t) ** h / mp.factorial(h) for t in d.loglo]
    )
    ce = mp.fsum(
        [mp.exp((rho - sig) * t) * (rho * t) ** h / mp.factorial(h) for t in d.logmid]
    )
    C = mp.expm1(L * rho)
    Ch = (L * rho) ** h * mp.exp(L * rho) / mp.factorial(h)
    rem = ce + Ch * sd + C * cd
    return mp.fsum([abs(P[m]) * rho**m for m in range(2, M + 1)]) + rem


def gbound(d, sig, extra=2):
    if not hasattr(d, "lv"):
        cur = list(d.mid)
        d.lv = []
        for i in range(extra):
            d.lv.append([mp.log(n) for n in cur])
            cur = sorted(d.q * m + a for m in cur for a in d.F)
    t = d.tail(d.P + extra, sig)
    if t == mp.inf:
        return mp.inf
    for lg in d.lv:
        t += mp.fsum([mp.exp(-sig * x) for x in lg])
    return t


def tbound(d, s0, R2):
    sig = mp.re(s0) - R2
    aw = abs(s0) + R2
    rat = d.amax * mp.power(d.q, -d.P)
    tot = mp.mpf(0)
    for l in range(1, 600):
        t = gbound(d, sig + l)
        if t == mp.inf:
            return mp.inf
        term = mp.binomial(aw + l - 1, l) * mp.power(d.q, -sig - l) * d.gam[l] * t
        tot += term
        rl = (aw + l) / (l + 1) * rat
        if rl < 1:
            return tot + term * rl / (1 - rl)
    return mp.inf


def certify(d, s0, Z0c, e0, samp, R, N, rhos, r2s):
    P, _ = explicit(d, s0)
    p1 = mp.log(d.q) * d.poly_lo(s0) - mp.fsum([mp.exp(-s0 * t) * t for t in d.logmid])
    maxe = max(e for v, e in samp)
    w = mp.exp(2 * mp.pi * 1j / N)
    T = []
    for m in range(N):
        u = R * w**m
        Pu = (1 - mp.exp(-mp.log(d.q) * u)) * d.poly_lo(s0 + u) + d.poly_mid(s0 + u)
        T.append(samp[m][0] - Pu)
    c1 = mp.fsum([T[m] * w ** (-m) for m in range(N)]) / (N * R)
    U0 = abs(Z0c) + e0
    best = None
    for R2 in r2s:
        if R2 <= R:
            continue
        B = tbound(d, s0, mp.mpf(R2))
        if B == mp.inf:
            continue
        al = (mp.mpf(R) / R2) ** N
        e1 = maxe / R + (B / R2) * al / (1 - al)
        Z1 = p1 + c1
        L1 = abs(Z1) - e1
        if L1 <= 0:
            continue
        for rho in rhos:
            rho = mp.mpf(rho)
            if rho >= R2:
                continue
            minM = L1 * rho - U0
            if minM <= 0:
                continue
            tau = rho / R2
            marg = minM - ptail(d, s0, rho, P) - B * tau * tau / (1 - tau)
            if best is None or marg > best[0]:
                best = (marg, rho, mp.mpf(R2), minM, B, L1, U0)
    return best, U0


def rouche(q, F, ymax, R=0.2, N=24, verbose=True, pool=1200, secs=None):
    base = build(q, F, dps=40, tol=mp.mpf(10) ** -24)
    fw = ZCache(build(q, F, dps=25, tol=mp.mpf(10) ** -10))
    L = mp.log(q)
    per = 2 * mp.pi / L
    depths, P = [], base.P
    while len(F) ** P <= pool and len(depths) < 3:
        depths.append(P)
        P += 1
    ds = {base.P: base}
    rhos = [mp.mpf(x) / 400 for x in range(6, 181)]
    r2s = [0.22, 0.25, 0.28, 0.31, 0.34, 0.37, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9]
    if verbose:
        line(
            " ",
            tag(q, F),
            "alpha",
            mp.nstr(alpha_of(q, F), 10),
            "sample circle R",
            R,
            "samples",
            N,
            "peel depths tried",
            depths,
        )
    t0 = time.time()
    out = []
    for j in range(int(mp.ceil(mp.mpf(ymax) / per)) + 1):
        old = mp.mp.dps
        mp.mp.dps = base.dps + 40 + int(abs(j) * 2)
        s0 = mp.log(base.k) / L + 2 * mp.pi * 1j * j / L
        Z0c, e0 = zwork(base, s0)
        mp.mp.dps = old
        if abs(Z0c) < 1000 * e0 or abs(Z0c) < mp.mpf(10) ** -15:
            if verbose:
                line("    j", j, "residue null, no tooth to certify")
            out.append({"j": j, "null": True, "cert": False, "skip": False})
            continue
        if secs is not None and time.time() - t0 > secs:
            if verbose:
                line("    j", j, "SKIPPED, design time budget spent")
            out.append({"j": j, "null": False, "cert": False, "skip": True})
            continue
        best, used = None, None
        for P in depths:
            if P not in ds:
                ds[P] = build(q, F, dps=40, tol=mp.mpf(10) ** -24, P=P)
            d = ds[P]
            mp.mp.dps = d.dps + 40 + int(abs(j) * 2)
            z0, ee = zwork(d, s0)
            samp = [
                zwork(d, s0 + mp.mpf(R) * mp.exp(2 * mp.pi * 1j * mp.mpf(m) / N))
                for m in range(N)
            ]
            b, _ = certify(d, s0, z0, ee, samp, mp.mpf(R), N, rhos, r2s)
            mp.mp.dps = old
            if b is not None and (best is None or b[0] > best[0]):
                best, used = b, P
            if b is not None and b[0] > 0:
                break
        n45 = disc(fw, s0, mp.mpf("0.45"))[0]
        rec = {
            "j": j,
            "null": False,
            "skip": False,
            "n45": n45,
            "best": best,
            "P": used,
            "cert": best is not None and best[0] > 0,
        }
        if rec["cert"]:
            rec["nrho"] = disc(fw, s0, best[1])[0]
        out.append(rec)
        if not verbose:
            continue
        if best is None:
            line(
                "    j",
                j,
                "no admissible radius at any depth",
                depths,
                "zeros in 0.45",
                mp.nstr(n45, 6),
            )
        else:
            marg, rho, R2, minM, B, L1, U0b = best
            line(
                "    j",
                j,
                "CERTIFIED" if marg > 0 else "failed",
                "peel P",
                used,
                "margin",
                mp.nstr(marg, 8),
                "rho",
                mp.nstr(rho, 5),
                "R2",
                mp.nstr(R2, 4),
                "min|model|",
                mp.nstr(minM, 6),
                "B_T",
                mp.nstr(B, 6),
                "|Z_1| low",
                mp.nstr(L1, 6),
                "|Z_0| up",
                mp.nstr(U0b, 6),
                "|Z_0/Z_1|",
                mp.nstr(U0b / L1, 6),
                "samples",
                N,
                "zeros in rho",
                "-" if marg <= 0 else mp.nstr(rec["nrho"], 6),
                "zeros in 0.45",
                mp.nstr(n45, 6),
            )
    return out


def main():
    argv = sys.argv[1:]
    verb = argv[0] if argv else "shadow"
    ymax = float(argv[1]) if len(argv) > 1 else 30.0
    sets = {
        "q3": Q3,
        "q4": Q4,
        "q5": Q5,
        "half": HALF,
        "wide": WIDE,
        "ctl": CONTROL,
        "all": Q3 + Q4 + Q5 + HALF + CONTROL + WIDE,
        "rou": Q3 + Q5 + HALF + WIDE,
    }
    which = sets[argv[2]] if len(argv) > 2 else sets["all"]
    t0 = time.time()
    if verb == "rouche":
        line(
            "ROUCHE certificate: exactly one zero of Z, hence of zeta_F, in abs(u) < rho"
        )
        tot = {"cert": 0, "fail": 0, "one": 0, "not one": 0, "null": 0, "skip": 0}
        for q, F in which:
            for rec in rouche(q, F, ymax, secs=300 if q == 10 else 105):
                if rec.get("null"):
                    tot["null"] += 1
                    continue
                if rec.get("skip"):
                    tot["skip"] += 1
                    continue
                tot["cert" if rec["cert"] else "fail"] += 1
                if rec["cert"]:
                    tot[
                        "one" if abs(rec["nrho"] - 1) < mp.mpf("0.2") else "not one"
                    ] += 1
        line(
            "CERTIFIED",
            tot["cert"],
            "FAILED",
            tot["fail"],
            "null poles",
            tot["null"],
            "skipped for budget",
            tot["skip"],
            "certified discs whose argument principle count is one",
            tot["one"],
            "certified discs disagreeing",
            tot["not one"],
        )
    elif verb == "shadow":
        line("SHADOW derived from the residue and the regular part, no fit")
        for q, F in which:
            comb(q, F, ymax)
    elif verb == "census":
        line("CENSUS zeros of the cofactor Z, one strip of width 4, every design deep")
        got = {}
        for q, F in which:
            c = comb(q, F, ymax, verbose=True)
            zs, tot, _ = census(q, F, ymax, int(ymax / 2), cmb=c, deep=True)
            nul = sum(
                1 for x in c if x["null"] and x["j"] > 0 and mp.im(x["s0"]) <= ymax
            )
            done = len(zs) == int(mp.nint(tot)) - nul
            got[(q, F)] = [(z, k, done) for z, k in zs]
        laws(got)
    line("seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
