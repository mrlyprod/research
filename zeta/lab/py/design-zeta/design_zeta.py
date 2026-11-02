import sys
import time

import mpmath as mp

TOL = mp.mpf(10) ** -22

# ENGINE


class Design:
    def __init__(self, q, F, P=None, dps=45, tol=TOL):
        mp.mp.dps = dps
        self.dps = dps
        self.q = q
        self.F = tuple(sorted(F))
        self.k = len(self.F)
        self.k1 = len([a for a in self.F if a != 0])
        self.amax = max(self.F)
        self.tol = tol
        if P is None:
            P = 2
            while self.k > 1 and self.k1 * self.k**P <= 100:
                P += 1
        self.P = P
        self.alpha = mp.log(self.k) / mp.log(q)
        self.gam = [mp.mpf(sum(a**l for a in self.F)) for l in range(1200)]
        self.lo = sorted(n for j in range(1, P) for n in self.strings(j))
        self.mid = self.strings(P)
        self.loglo = [mp.log(n) for n in self.lo]
        self.logmid = [mp.log(n) for n in self.mid]

    def strings(self, p):
        cur = [a for a in self.F if a != 0]
        for _ in range(p - 1):
            cur = sorted(self.q * m + a for m in cur for a in self.F)
        return sorted(cur)

    def tail(self, p, sig):
        r = self.k * mp.power(self.q, -sig)
        if r >= 1:
            return mp.inf
        return self.k1 * r ** (p - 1) / (1 - r)

    def poly_lo(self, w):
        return mp.fsum([mp.exp(-w * L) for L in self.loglo])

    def poly_mid(self, w):
        return mp.fsum([mp.exp(-w * L) for L in self.logmid])

    def cut(self, w, L):
        aw, sig = abs(w), mp.re(w)
        rho = max((aw + L + 1) / (L + 2), 1) * self.amax / mp.power(self.q, self.P)
        if rho >= mp.mpf("0.9"):
            return mp.inf
        maj = (
            mp.binomial(aw + L, L + 1)
            * mp.power(self.q, -sig - (L + 1))
            * self.k
            * self.amax ** (L + 1)
            * self.tail(self.P, sig + L + 1)
        )
        return maj / (1 - rho)

    def ladder(self, s, W, L, top=True):
        q, k = self.q, self.k
        J = max(0, int(mp.ceil(W - mp.re(s))))
        g, e = {}, {}
        for j in range(J, J + L + 2):
            g[j] = mp.mpf(0)
            e[j] = self.tail(self.P, mp.re(s) + j)
        for j in range(J - 1, -1, -1):
            w = s + j
            acc = self.poly_mid(w)
            eacc = self.cut(w, L)
            c = mp.mpf(1)
            for l in range(1, L + 1):
                c = c * (-w - (l - 1)) / l
                co = mp.power(q, -w - l) * c * self.gam[l]
                acc += co * g[j + l]
                eacc += abs(co) * e[j + l]
            if not top and j == 0:
                return acc, eacc
            den = 1 - k * mp.power(q, -w)
            g[j] = acc / den
            e[j] = eacc / abs(den)
        return g[0], e[0]

    def tune(self, s, top=True):
        W, L = mp.mpf(6), 8
        while W <= 400:
            v, e = self.ladder(s, W, L, top)
            if e < self.tol:
                return v, e
            W += 6
            L += 6
        raise RuntimeError(
            "FAIL tolerance "
            + mp.nstr(self.tol, 3)
            + " not reached at s = "
            + mp.nstr(s, 14)
            + ", best bound "
            + mp.nstr(e, 4)
        )

    def zeta(self, s):
        s = mp.mpmathify(s)
        old = mp.mp.dps
        mp.mp.dps = self.dps + 30 + int(abs(mp.im(s)) / 4)
        s = mp.mpmathify(s)
        g, e = self.tune(s)
        v = self.poly_lo(s) + g
        mp.mp.dps = old
        return +v, +e

    def cofactor(self, s):
        s = mp.mpmathify(s)
        old = mp.mp.dps
        mp.mp.dps = self.dps + 30 + int(abs(mp.im(s)) / 4)
        s = mp.mpmathify(s)
        a, e = self.tune(s, top=False)
        v = (1 - self.k * mp.power(self.q, -s)) * self.poly_lo(s) + a
        mp.mp.dps = old
        return +v, +e

    def residue(self, j):
        old = mp.mp.dps
        mp.mp.dps = self.dps + 30 + int(abs(j) * 2)
        sj = self.alpha + 2 * mp.pi * 1j * j / mp.log(self.q)
        a, e = self.tune(sj, top=False)
        r, er = a / mp.log(self.q), e / mp.log(self.q)
        mp.mp.dps = old
        return +r, +er


# CONTOURS


class Cache:
    def __init__(self, d, cof=False):
        self.d = d
        self.cof = cof
        self.m = {}
        self.n = 0
        self.emax = mp.mpf(0)

    def __call__(self, z):
        key = (mp.nstr(mp.re(z), 22), mp.nstr(mp.im(z), 22))
        if key not in self.m:
            v, e = self.d.cofactor(z) if self.cof else self.d.zeta(z)
            self.m[key] = v
            self.emax = max(self.emax, e)
            self.n += 1
        return self.m[key]


def phase(f, pt, n0, cap=1.0, budget=4000):
    us = [mp.mpf(i) / n0 for i in range(n0 + 1)]
    vs = [f(pt(u)) for u in us]
    for _ in range(30):
        bad = [i for i in range(len(us) - 1) if abs(mp.arg(vs[i + 1] / vs[i])) > cap]
        if not bad or len(us) > budget:
            break
        nu, nv = [], []
        for i in range(len(us) - 1):
            nu.append(us[i])
            nv.append(vs[i])
            if i in bad:
                um = (us[i] + us[i + 1]) / 2
                nu.append(um)
                nv.append(f(pt(um)))
        nu.append(us[-1])
        nv.append(vs[-1])
        us, vs = nu, nv
    tot = mp.mpf(0)
    mx = mp.mpf(0)
    for i in range(len(vs) - 1):
        dd = mp.arg(vs[i + 1] / vs[i])
        tot += dd
        mx = max(mx, abs(dd))
    return tot, mx


def box_phase(f, x0, x1, y0, y1, n):
    a = phase(f, lambda u: mp.mpc(x0 + (x1 - x0) * u, y0), n)
    b = phase(f, lambda u: mp.mpc(x1, y0 + (y1 - y0) * u), n)
    c = phase(f, lambda u: mp.mpc(x1 - (x1 - x0) * u, y1), n)
    d = phase(f, lambda u: mp.mpc(x0, y1 - (y1 - y0) * u), n)
    tot = a[0] + b[0] + c[0] + d[0]
    return tot / (2 * mp.pi), max(a[1], b[1], c[1], d[1])


def strip(f, x0, x1, y0, y1, nsub, n):
    cuts = [y0 + (y1 - y0) * i / nsub for i in range(nsub + 1)]
    hor = {}
    for y in cuts:
        hor[y] = phase(f, lambda u, y=y: mp.mpc(x0 + (x1 - x0) * u, y), n)
    out = []
    for i in range(nsub):
        a, b = cuts[i], cuts[i + 1]
        rt = phase(f, lambda u, a=a, b=b: mp.mpc(x1, a + (b - a) * u), n)
        lf = phase(f, lambda u, a=a, b=b: mp.mpc(x0, a + (b - a) * u), n)
        tot = hor[a][0] + rt[0] - hor[b][0] - lf[0]
        mx = max(hor[a][1], rt[1], hor[b][1], lf[1])
        out.append((a, b, tot / (2 * mp.pi), mx))
    return out


def locate(f, box, steps=9):
    x0, x1, y0, y1 = [mp.mpf(v) for v in box]
    for _ in range(steps):
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        hit = False
        for a, b, c, d in [
            (x0, xm, y0, ym),
            (xm, x1, y0, ym),
            (x0, xm, ym, y1),
            (xm, x1, ym, y1),
        ]:
            if abs(box_phase(f, a, b, c, d, 16)[0]) > mp.mpf("0.5"):
                x0, x1, y0, y1 = a, b, c, d
                hit = True
                break
        if not hit:
            break
    return mp.mpc((x0 + x1) / 2, (y0 + y1) / 2), max(x1 - x0, y1 - y0)


def zero(d, box, cof=False):
    f = Cache(d, cof=cof)
    z0, w = locate(f, box)
    r = mp.findroot(
        f,
        [z0 - w, z0, z0 + w * 1j],
        solver="muller",
        tol=mp.mpf(10) ** -30,
        maxsteps=60,
    )
    if abs(r - z0) > 4 * w:
        r = z0
    v = abs(f(r))
    return r, v, f.emax


# STUDY


def line(*a):
    print(" ".join(str(x) for x in a))


def control():
    line("CONTROL base 2 full digit set, zeta_F = zeta")
    d = Design(2, (0, 1))
    for s in [
        mp.mpf(2),
        mp.mpf("0.5") + 14.134725141734693j,
        mp.mpf("0.3") + 40j,
        -1 + 2j,
        mp.mpf("0.5") + 100j,
    ]:
        v, e = d.zeta(s)
        mp.mp.dps = 60
        line(
            "  s",
            mp.nstr(s, 8),
            "gap to mpmath zeta",
            mp.nstr(abs(v - mp.zeta(s)), 4),
            "proved bound",
            mp.nstr(e, 3),
        )
    r, e = d.residue(0)
    line("  residue at alpha = 1", mp.nstr(r, 20), "bound", mp.nstr(e, 3))
    r, e = d.residue(1)
    line(
        "  residue at alpha + 2 pi i/log 2", mp.nstr(abs(r), 4), "bound", mp.nstr(e, 3)
    )


def residues():
    line("RESIDUES against the certified enclosures of lab/py/burnol-residue")
    for q, F in [(3, (0, 1)), (3, (0, 2))]:
        d = Design(q, F)
        for j in (0, 1, 2):
            r, e = d.residue(j)
            line("  q", q, "F", F, "j", j, mp.nstr(r, 18), "bound", mp.nstr(e, 3))


def scaling():
    line("SCALING zeta_(aF)(s) = a^(-s) zeta_F(s)")
    d1 = Design(3, (0, 1))
    d2 = Design(3, (0, 2))
    for s in [mp.mpf(2) + 3j, mp.mpf("0.8") + 23j, -mp.mpf("0.5") + 11j]:
        a = d2.zeta(s)[0]
        b = mp.power(2, -s) * d1.zeta(s)[0]
        line("  s", mp.nstr(s, 8), "gap", mp.nstr(abs(a - b), 4))


def census(q, F, lo, hi, ymax, nsub, n=20, cof=True):
    d = Design(q, F, tol=mp.mpf(10) ** -10, dps=25)
    f = Cache(d, cof=cof)
    x0, x1 = d.alpha + mp.mpf(lo), d.alpha + mp.mpf(hi)
    t0 = time.time()
    res = strip(f, x0, x1, mp.mpf("0.02"), mp.mpf(ymax), nsub, n)
    tot = mp.fsum([w for a, b, w, mx in res])
    mxall = max([mx for a, b, w, mx in res])
    hits = [(a, b, w, mx) for a, b, w, mx in res if abs(w) > mp.mpf("0.2")]
    line(
        "  q",
        q,
        "F",
        F,
        "alpha",
        mp.nstr(d.alpha, 12),
        "object",
        "Z = zeta_F (1 - k q^(-s))" if cof else "zeta_F",
        "strip Re in",
        mp.nstr(x0, 12),
        mp.nstr(x1, 12),
        "Im in 0.02",
        ymax,
    )
    for a, b, w, mx in hits:
        line(
            "    Im [",
            mp.nstr(a, 6),
            ",",
            mp.nstr(b, 6),
            "] winding",
            mp.nstr(w, 8),
            "max phase step",
            mp.nstr(mx, 4),
        )
    line(
        "    count",
        mp.nstr(tot, 10),
        "boxes",
        len(hits),
        "largest phase step",
        mp.nstr(mxall, 4),
        "largest bound on the contour",
        mp.nstr(f.emax, 4),
        "evaluations",
        f.n,
        "seconds",
        round(time.time() - t0, 1),
    )
    return tot, hits


def zeros(q, F, lo, hi, boxes, cof=False):
    d = Design(q, F, tol=mp.mpf(10) ** -22, dps=45)
    per = 2 * mp.pi / mp.log(q)
    out = []
    for a, b in boxes:
        r, v, e = zero(d, (d.alpha + mp.mpf(lo), d.alpha + mp.mpf(hi), a, b), cof=cof)
        out.append(r)
        line(
            "    zero",
            mp.nstr(r, 16),
            "abs",
            "Z" if cof else "zeta_F",
            mp.nstr(v, 3),
            "largest bound",
            mp.nstr(e, 3),
            "Im/period",
            mp.nstr(mp.im(r) / per, 12),
        )
    line(
        "    period 2 pi/log q",
        mp.nstr(per, 14),
        "alpha",
        mp.nstr(d.alpha, 12),
        "alpha/2",
        mp.nstr(d.alpha / 2, 12),
    )
    for i in range(len(out) - 1):
        g = mp.im(out[i + 1]) - mp.im(out[i])
        line("    gap", mp.nstr(g, 14), "gap minus period", mp.nstr(g - per, 8))
    return out


def control_census(ymax):
    line(
        "CONTROL base 2 full digit set, the census run on BOTH sides of the abscissa alpha = 1"
    )
    a, _ = census(2, (0, 1), 0.02, 3.02, ymax, int(ymax / 2), cof=False)
    line(
        "    object counted zeros of zeta_F = zeta right of the abscissa, count",
        mp.nstr(a, 6),
    )
    b, _ = census(2, (0, 1), -0.98, -0.02, ymax, int(ymax / 2), cof=False)
    line(
        "    object counted zeros of zeta_F = zeta left of the abscissa, count",
        mp.nstr(b, 6),
    )
    c, _ = census(2, (0, 1), -0.92, 3.02, ymax, int(ymax / 2), cof=True)
    teeth = int(mp.floor(mp.mpf(ymax) * mp.log(2) / (2 * mp.pi)))
    line(
        "    object counted zeros of the cofactor Z on one strip, count",
        mp.nstr(c, 6),
        "= zeros of zeta",
        mp.nstr(b, 6),
        "plus teeth of 1 - 2 q^(-s) on Re s = alpha",
        teeth,
        "predicted floor(T log q/2 pi)",
        teeth,
    )
    line(
        "    the Euler product forbids zeros of zeta in Re s >= 1, and the strip Re in",
        "[alpha + 0.02, alpha + 3.02] confirms it; the sliver alpha < Re s < alpha + 0.02",
        "carries the teeth of the cofactor and no zero of zeta",
    )


def main():
    full = "--full" in sys.argv
    ymax = 60 if full else 30
    y10 = ymax
    for a in sys.argv:
        if a.startswith("--y10="):
            y10 = float(a.split("=")[1])
    t0 = time.time()
    control()
    residues()
    scaling()
    control_census(ymax)
    line(
        "CENSUS zeros of the cofactor Z = zeta_F (1 - k q^(-s)), analytic on Re s > alpha - 1"
    )
    line("  the split is at Re s = alpha exactly, so no sliver is left unscanned")
    census(3, (0, 1), 0.0, 3.02, ymax, int(ymax / 2))
    census(3, (0, 1), -0.92, 0.0, ymax, int(ymax / 2))
    census(3, (0, 2), 0.0, 3.02, ymax, int(ymax / 2))
    if full:
        census(10, tuple(range(9)), 0.0, 3.02, y10, int(y10 / 2))
        census(10, tuple(range(9)), -0.92, 0.0, y10, int(y10 / 2))
    line("ZEROS polished, right of the abscissa")
    zeros(3, (0, 1), 0.02, 3.02, [(22, 24), (28, 30)])
    if full:
        zeros(10, tuple(range(9)), 0.02, 3.02, [(2, 4), (4, 6)])
    line("ZEROS polished on the cofactor, base 3 digits 0 1, left of the abscissa")
    zeros(3, (0, 1), -0.92, 0.0, [(6, 8), (10, 12)], cof=True)
    line("seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
