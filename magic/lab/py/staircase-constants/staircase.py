import math
import time
from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 50
HALF = mp.mpf(3) / 2


def nine(x):
    x = mp.re(x)
    return mp.nstr(
        x, 9 + max(0, int(mp.floor(mp.log10(abs(x)))) + 1), strip_zeros=False
    )


# CORNERS


def corners(code, dim):
    return [i for i in range(2**dim) if code >> i & 1]


def profile(code, dim):
    a = [0] * (dim + 1)
    for i in corners(code, dim):
        a[i.bit_count()] += 1
    return tuple(a)


def profiles(dim):
    seen = {}
    for code in range(1, 2 ** (2**dim)):
        seen.setdefault(profile(code, dim), []).append(code)
    return seen


def coefficients(a):
    dim = len(a) - 1
    c = [0] * (dim + 1)
    for j, aj in enumerate(a):
        for k in range(j + 1):
            c[k] += aj * math.comb(j, k) * (-1) ** k
    return c


def fill(a, n):
    dim = len(a) - 1
    return sum(aj * n ** (dim - j) * (n - 1) ** j for j, aj in enumerate(a))


def shape(a):
    dim = len(a) - 1
    w = sum(a)
    mean = Fraction(sum(j * aj for j, aj in enumerate(a)), w)
    var = Fraction(sum(j * j * aj for j, aj in enumerate(a)), w) - mean**2
    drift = Fraction(dim, 2) - mean
    return w, drift, Fraction(dim, 8) + drift - var / 2


# RENDER


def tile(code, dim, side):
    axes = np.indices((side,) * dim)
    index = sum((axes[j] % 2) << (dim - 1 - j) for j in range(dim))
    return ((code >> index) & 1).astype(np.int64)


def section_render():
    bad = 0
    checked = 0
    for dim in (1, 2, 3):
        for code in range(1, 2 ** (2**dim)):
            a = profile(code, dim)
            w = sum(a)
            for side in (3, 5, 7, 9):
                bad += int(tile(code, dim, side).sum()) != fill(a, (side + 1) // 2)
                checked += 1
            for side in (2, 4, 6, 8):
                bad += int(tile(code, dim, side).sum()) * 2**dim != w * side**dim
                checked += 1
            levels = {1: 5, 2: 3, 3: 2}[dim]
            word = tile(code, dim, 3)
            for level in range(2, levels + 1):
                word = np.kron(word, tile(code, dim, 2 * level + 1))
                want = math.prod(fill(a, k + 1) for k in range(1, level + 1))
                bad += int(word.sum()) != want
                checked += 1
    print(
        f"render: {checked} counts over every code in dim 1, 2, 3 (odd sides 3..9, even sides 2..8, row words), mismatches {bad}"
    )
    named = {
        (2, 7): (1, 2, 0),
        (2, 14): (0, 2, 1),
        (2, 6): (0, 2, 0),
        (2, 9): (1, 0, 1),
        (3, 23): (1, 3, 0, 0),
        (3, 232): (0, 0, 3, 1),
        (3, 150): (0, 3, 0, 1),
        (3, 105): (1, 0, 3, 0),
    }
    wrong = [k for k, v in named.items() if profile(k[1], k[0]) != v]
    print(
        f"render: crate corner order (corner i = binary digits of i, most significant first), named profiles wrong {len(wrong)}"
    )


# ROOTS


def roots(a):
    n = sp.Symbol("n")
    p = sp.Poly(coefficients(a), n)
    out = []
    for factor, mult in sp.factor_list(p.as_expr())[1]:
        f = sp.Poly(factor, n)
        if f.degree() == 1:
            r = -f.all_coeffs()[1] / f.all_coeffs()[0]
            out += [mp.mpf(sp.Rational(r).p) / sp.Rational(r).q] * mult
        else:
            cs = [int(x) for x in f.all_coeffs()]
            out += list(mp.polyroots(cs, maxsteps=400, extraprec=400)) * mult
    return out


def constant(rs, dim):
    value = HALF_GAMMA**dim / mp.fprod(mp.gamma(2 - r) for r in rs)
    assert abs(mp.im(value)) < mp.mpf(10) ** -40 and mp.re(value) > 0
    return mp.re(value)


HALF_GAMMA = mp.gamma(HALF)


def log_exact_formula(a, rs, level):
    dim = len(a) - 1
    w = sum(a)
    s = level * mp.log(mp.mpf(w) / 2**dim)
    for r in rs:
        s += mp.loggamma(level + 2 - r) - mp.loggamma(2 - r)
    s -= dim * (mp.loggamma(level + HALF) - mp.loggamma(HALF))
    return mp.re(s)


# EXACT


def section_exact(table):
    worst = mp.mpf(0)
    for a, rs in table:
        dim = len(a) - 1
        top, bottom = 1, 1
        marks = {1, 2, 3, 10, 100, 1000}
        for k in range(1, 1001):
            top *= fill(a, k + 1)
            bottom *= (2 * k + 1) ** dim
            if k in marks:
                gap = abs(mp.log(top) - mp.log(bottom) - log_exact_formula(a, rs, k))
                worst = max(worst, gap)
    print(
        f"exact: R_L = (w/2^dim)^L prod Gamma(L+2-r)/Gamma(2-r) / (Gamma(L+3/2)/Gamma(3/2))^dim on {len(table)} profiles at L = 1, 2, 3, 10, 100, 1000, worst log gap {mp.nstr(worst, 3)}"
    )


# ASYMPTOTIC


def direct_logs(a, top):
    dim = len(a) - 1
    w = sum(a)
    big = sp.Symbol("N")
    poly = sp.Poly(
        sum(aj * (big + 1) ** (dim - j) * (big - 1) ** j for j, aj in enumerate(a))
        - w * big**dim,
        big,
    )
    d = [float(x) for x in poly.all_coeffs()] if not poly.is_zero else [0.0]
    side = 2.0 * np.arange(1, top + 1) + 1.0
    return np.log1p(np.polyval(d, side) / (w * side**dim))


def section_asymptotic(table, top=10**6):
    marks = (10**4, 10**5, top)
    worst = {m: 0.0 for m in marks}
    for a, rs in table:
        dim = len(a) - 1
        _w, drift, c1 = shape(a)
        c1_roots = mp.re(sum((mp.mpf(1) / 2 - r) * (mp.mpf(5) / 2 - r) for r in rs) / 2)
        assert abs(c1_roots - mp.mpf(c1.numerator) / c1.denominator) < mp.mpf(10) ** -30
        logs = direct_logs(a, top)
        log_c = float(mp.log(constant(rs, dim)))
        for m in marks:
            got = float(np.sum(logs[:m]))
            want = float(drift) * math.log(m) + log_c + float(c1) / m
            worst[m] = max(worst[m], abs(got - want))
    print(
        f"asymptotic: direct log sum against drift log L + log C + c1/L, c1 = dim/8 + drift - var/2, on {len(table)} profiles"
    )
    for m in marks:
        print(
            f"  L = {m:>8}: worst residual {worst[m]:.2e}, times L^2 {worst[m] * m * m:.2f}"
        )


# CLOSED FORMS


def section_named(byprofile):
    pi = mp.pi
    sq = mp.sqrt
    forms = [
        (1, 1, "sqrt(pi)/2", sq(pi) / 2),
        (1, 2, "sqrt(pi)/2", sq(pi) / 2),
        (2, 7, "3 pi/(4 Gamma(1/3))", 3 * pi / (4 * mp.gamma(mp.mpf(1) / 3))),
        (2, 14, "3 pi/(8 Gamma(2/3))", 3 * pi / (8 * mp.gamma(mp.mpf(2) / 3))),
        (2, 6, "pi/4", pi / 4),
        (2, 9, "cosh(pi/2)/2", mp.cosh(pi / 2) / 2),
        (2, 11, "3 cosh(pi/(2 sqrt 3))/4", 3 * mp.cosh(pi / (2 * sq(3))) / 4),
        (2, 1, "pi/4", pi / 4),
        (2, 3, "sqrt(pi)/2", sq(pi) / 2),
        (3, 23, "pi^(3/2)/(2 Gamma(1/4))", pi**HALF / (2 * mp.gamma(mp.mpf(1) / 4))),
        (3, 232, "pi^(3/2)/(6 Gamma(3/4))", pi**HALF / (6 * mp.gamma(mp.mpf(3) / 4))),
        (
            3,
            150,
            "pi^(3/2)/(8 |Gamma(7/4 - i sqrt3/4)|^2)",
            pi**HALF / (8 * abs(mp.gamma(mp.mpc(7, -sq(3)) / 4)) ** 2),
        ),
        (
            3,
            105,
            "pi^(3/2)/(8 |Gamma(5/4 - i sqrt3/4)|^2)",
            pi**HALF / (8 * abs(mp.gamma(mp.mpc(5, -sq(3)) / 4)) ** 2),
        ),
        (3, 129, "cosh(pi sqrt3/2)/4", mp.cosh(pi * sq(3) / 2) / 4),
        (3, 126, "pi/4", pi / 4),
        (3, 24, "pi/4", pi / 4),
    ]
    print(
        "named: dim, code, profile, w/2^dim, drift, C (Gamma at roots), closed form, gap"
    )
    worst = mp.mpf(0)
    cs = {}
    for dim, code, text, value in forms:
        a = profile(code, dim)
        rs = byprofile[a]
        w, drift, _ = shape(a)
        c = constant(rs, dim)
        cs[(dim, code)] = c
        worst = max(worst, abs(c - value))
        print(
            f"  dim {dim} code {code:>3} {a!s:<14} {Fraction(w, 2**dim)!s:>4} {drift!s:>5} {mp.nstr(c, 12):<15} {nine(c):<12} {text}  {mp.nstr(abs(c - value), 3)}"
        )
    print(f"named: worst closed-form gap {mp.nstr(worst, 3)}")
    mirrors = [
        ((2, 7), (2, 14), "9 sqrt3 pi/64", 9 * sq(3) * pi / 64),
        ((3, 23), (3, 232), "sqrt2 pi^2/24", sq(2) * pi**2 / 24),
        ((1, 1), (1, 2), "pi/4", pi / 4),
    ]
    for f, g, text, value in mirrors:
        print(
            f"  mirror dim {f[0]} codes {f[1]} x {g[1]}: {mp.nstr(cs[f] * cs[g], 15)} against {text} gap {mp.nstr(abs(cs[f] * cs[g] - value), 3)}"
        )


# PAIRING


def reflect(rs):
    out = mp.mpf(1)
    for r in rs:
        if abs(r) < mp.mpf(10) ** -40 or abs(r - 1) < mp.mpf(10) ** -40:
            out *= mp.pi / 4
        else:
            out *= mp.sin(mp.pi * r) / (4 * r * (1 - r))
    return mp.re(out)


def section_pairing(byprofile):
    n = sp.Symbol("n")
    agree = 0
    total = 0
    for dim in (1, 2, 3, 4):
        for a in profiles(dim) if dim < 4 else four_profiles():
            p = sp.Poly(coefficients(a), n)
            q = sp.Poly(p.as_expr().subs(n, 1 - n) * (-1) ** dim, n)
            flipped = p == q
            palindrome = tuple(a) == tuple(reversed(a))
            agree += flipped == palindrome
            total += 1
    print(
        f"pairing: P(1-n) = (-1)^dim P(n) exactly when the profile is a palindrome, {agree} of {total} profiles at dim 1..4"
    )
    worst = mp.mpf(0)
    for a, rs in byprofile.items():
        dim = len(a) - 1
        mirror = tuple(reversed(a))
        worst = max(
            worst,
            abs(constant(rs, dim) * constant(byprofile[mirror], dim) - reflect(rs)),
        )
    print(
        f"pairing: C_F C_mirror = prod sin(pi r)/(4 r (1-r)) on all {len(byprofile)} profiles at dim 1..3, worst gap {mp.nstr(worst, 3)}"
    )
    flips = 0
    palins = 0
    for dim in (2, 3):
        for code in range(1, 2 ** (2**dim)):
            cs = set(corners(code, dim))
            closed = all((2**dim - 1 - i) in cs for i in cs)
            a = profile(code, dim)
            palin = a == tuple(reversed(a))
            flips += closed
            palins += palin
            assert palin or not closed
    print(
        f"pairing: codes closed under the flip of every coordinate {flips}, codes with a palindromic profile {palins}, at dim 2 and 3 together"
    )
    a = profile(11, 2)
    print(
        f"pairing: dim 2 code 11 corners {corners(11, 2)} profile {a}, not flip-closed, roots {[mp.nstr(r, 12) for r in byprofile[a]]}"
    )
    tiny = mp.mpf(10) ** -30

    def core(rs):
        return [r for r in rs if min(abs(r), abs(r - mp.mpf(1) / 2), abs(r - 1)) > tiny]

    def paired(rs, shifts):
        left = list(rs)
        while left:
            r = left.pop()
            hit = next(
                (
                    i
                    for i, s in enumerate(left)
                    if any(abs(s - (1 - r + k)) < tiny for k in shifts)
                ),
                None,
            )
            if hit is None:
                return False
            left.pop(hit)
        return True

    palin = [a for a in byprofile if a == tuple(reversed(a))]
    exact = [a for a, rs in byprofile.items() if paired(core(rs), (0,))]
    shifted = [a for a, rs in byprofile.items() if paired(core(rs), range(-4, 5))]
    extra = [a for a in exact if a not in palin]
    print(
        f"pairing: of {len(byprofile)} profiles at dim 1..3, palindromic {len(palin)}, core roots (outside 0, 1/2, 1) closed under r -> 1-r {len(exact)}, closed up to an integer shift {len(shifted)}"
    )
    print(
        f"pairing: reducible without a palindrome {len(extra)}: "
        + ", ".join(f"{a} code {byprofile_code(a)}" for a in extra)
    )
    rest = [a for a in byprofile if a not in shifted]
    print(
        f"pairing: unreduced {len(rest)}, rational core roots among them {sum(all(abs(mp.im(r)) < tiny and abs(mp.re(r) * 840 - round(mp.re(r) * 840)) < tiny for r in core(byprofile[a])) for a in rest)}"
    )


def byprofile_code(a):
    dim = len(a) - 1
    return next(c for c in range(1, 2 ** (2**dim)) if profile(c, dim) == a)


def four_profiles():
    out = []
    for a0 in range(2):
        for a1 in range(5):
            for a2 in range(7):
                for a3 in range(5):
                    for a4 in range(2):
                        if a0 + a1 + a2 + a3 + a4:
                            out.append((a0, a1, a2, a3, a4))
    return out


# PARITY


def parity_constant(dim, sign):
    us = [
        mp.exp(2j * mp.pi * (k + (0 if sign < 0 else mp.mpf(1) / 2)) / dim)
        for k in range(dim)
    ]
    return mp.re(HALF_GAMMA**dim / mp.fprod(mp.gamma((3 - u) / 2) for u in us))


def zeta_series(dim, sign):
    total = mp.mpf(0)
    for m in range(1, 400):
        tail = (1 - mp.mpf(2) ** (-dim * m)) * mp.zeta(dim * m) - 1
        term = tail / m * (-1 if sign < 0 else (-1) ** (m + 1))
        total += term
        if abs(tail) < mp.mpf(10) ** -60:
            break
    return mp.exp(total)


def section_parity():
    pi = mp.pi
    forms = {
        (2, -1): ("pi/4", pi / 4),
        (2, 1): ("cosh(pi/2)/2", mp.cosh(pi / 2) / 2),
        (4, -1): ("pi cosh(pi/2)/8", pi * mp.cosh(pi / 2) / 8),
        (4, 1): (
            "(cosh(pi/sqrt2) + cos(pi/sqrt2))/4",
            (mp.cosh(pi / mp.sqrt(2)) + mp.cos(pi / mp.sqrt(2))) / 4,
        ),
        (6, -1): ("pi cosh(pi sqrt3/2)/24", pi * mp.cosh(pi * mp.sqrt(3) / 2) / 24),
        (6, 1): (
            "cosh(pi/2)(cosh(pi/2) + cos(pi sqrt3/2))/4",
            mp.cosh(pi / 2) * (mp.cosh(pi / 2) + mp.cos(pi * mp.sqrt(3) / 2)) / 4,
        ),
    }
    print(
        "parity: dim, parity, code, Gamma at roots of unity, gap to the zeta log series exp(-+sum_m (+-1)^m (lambda(dim m) - 1)/m), direct log sum to L = 10^6 plus c1/L, closed form"
    )
    for dim in (2, 3, 4, 5, 6):
        for sign, name in ((-1, "odd"), (1, "even")):
            want = 1 if name == "odd" else 0
            code = sum(1 << i for i in range(2**dim) if i.bit_count() % 2 == want)
            a = profile(code, dim)
            g = parity_constant(dim, sign)
            _, drift, c1 = shape(a)
            assert drift == 0
            logs = np.log1p(sign * (2.0 * np.arange(1, 10**6 + 1) + 1.0) ** (-dim))
            direct = math.exp(float(np.sum(logs)) - float(c1) / 10**6)
            series = zeta_series(dim, sign)
            text, value = forms.get((dim, sign), ("", None))
            gap = (
                "" if value is None else f"closed-form gap {mp.nstr(abs(g - value), 3)}"
            )
            print(
                f"  dim {dim} {name:<4} code {code:<6} {nine(g):<12} series gap {mp.nstr(abs(g - series), 3):<9} direct {direct:.12f}  {text} {gap}"
            )


# FAMILY


def section_family():
    print(
        "family: at most one odd coordinate, dim, code, w/2^dim, drift, C at the roots of the fill polynomial, gap to (sqrt(pi)/2)^dim (dim+1)/Gamma(1/(dim+1)), direct at L = 10^6"
    )
    for dim in range(1, 9):
        a = tuple([1, dim] + [0] * (dim - 1))
        code = sum(1 << i for i in range(2**dim) if i.bit_count() <= 1)
        rs = roots(a)
        w, drift, c1 = shape(a)
        c = constant(rs, dim)
        form = (mp.sqrt(mp.pi) / 2) ** dim * (dim + 1) / mp.gamma(mp.mpf(1) / (dim + 1))
        logs = direct_logs(a, 10**6)
        direct = math.exp(
            float(np.sum(logs)) - float(drift) * math.log(10**6) - float(c1) / 10**6
        )
        print(
            f"  dim {dim} code {code:<20} {Fraction(w, 2**dim)!s:>6} {drift!s:>7} {mp.nstr(c, 12):<15} {nine(c):<12} gap {mp.nstr(abs(c - form), 3):<8} direct {direct:.10f}"
        )


# FENCE


def section_fence(byprofile):
    print(
        "fence: renormalised limit prod_k P((N_k+1)/2)/(w (N_k/2)^dim) along N_k = 3^k, k = 1..60"
    )
    for dim, code in ((1, 1), (2, 7), (3, 23), (2, 6), (2, 9)):
        a = profile(code, dim)
        w = sum(a)
        out = mp.mpf(1)
        for k in range(1, 61):
            side = 3**k
            out *= mp.mpf(fill(a, (side + 1) // 2) * 2**dim) / (w * mp.mpf(side) ** dim)
        print(
            f"  dim {dim} code {code}: {mp.nstr(out, 12)}, row-word constant {mp.nstr(constant(byprofile[a], dim), 12)}"
        )
    a = profile(7, 2)
    for top in (10**3, 10**5):
        primes = [p for p in sp.primerange(3, top)]
        s = mp.mpf(1)
        for p in primes:
            s *= mp.mpf(fill(a, (p + 1) // 2) * 4) / (3 * mp.mpf(p) ** 2)
        print(
            f"  dim 2 code 7 along the odd primes below {top}: renormalised product {mp.nstr(s, 8)}, sum 1/p {mp.nstr(sum(mp.mpf(1) / p for p in primes), 6)}"
        )


# STAIRCASE


def stair_exact(a, rs, n):
    dim = len(a) - 1
    w = sum(a)
    g = mp.barnesg
    top = n * (n + 1) / 2 * mp.log(mp.mpf(w) / 2**dim)
    for r in rs:
        top += mp.log(g(n + 3 - r)) - mp.log(g(3 - r)) - n * mp.loggamma(2 - r)
    top -= dim * (
        mp.log(g(n + mp.mpf(5) / 2)) - mp.log(g(mp.mpf(5) / 2)) - n * mp.loggamma(HALF)
    )
    side = (
        n * (n + 1) / 2 * mp.log(2)
        + mp.log(g(n + mp.mpf(5) / 2))
        - mp.log(g(mp.mpf(5) / 2))
        - n * mp.loggamma(HALF)
    )
    return dim + mp.re(top) / side


def stair_direct(a, n):
    top = sum((n - j + 1) * mp.log(fill(a, j + 1)) for j in range(1, n + 1))
    side = sum((n - j + 1) * mp.log(2 * j + 1) for j in range(1, n + 1))
    return top / side


def section_staircase(byprofile):
    a = profile(7, 2)
    print(
        "staircase: dim 2 code 7, dimension(n) by Barnes G at the roots and by the direct fill sum"
    )
    print(
        "  "
        + ", ".join(
            f"{mp.nstr(stair_exact(a, byprofile[a], n), 10)}" for n in range(1, 6)
        )
    )
    worst = mp.mpf(0)
    for b, rs in byprofile.items():
        for n in (1, 2, 3, 6, 20):
            worst = max(worst, abs(stair_exact(b, rs, n) - stair_direct(b, n)))
    print(
        f"staircase: Barnes G form against the direct sum on all {len(byprofile)} profiles at dim 1..3, n = 1, 2, 3, 6, 20, worst gap {mp.nstr(worst, 3)}"
    )
    top = 2 * 10**5
    marks = (10**2, 10**3, 10**4, 10**5, top)
    print(
        "staircase: gap = dim - dimension, against log(2^dim/w)/(log 2n - 3/2) (staircase) and /(log 2L - 1) (row word); n x error"
    )
    for dim, code in ((1, 1), (2, 7), (2, 9), (3, 23)):
        b = profile(code, dim)
        w = sum(b)
        lead = math.log(2**dim / w)
        logs = direct_logs(b, top)
        row_fill = np.cumsum(logs) + np.arange(1, top + 1) * math.log(w / 2**dim)
        row_side = np.cumsum(np.log(2.0 * np.arange(1, top + 1) + 1.0))
        stair_fill = np.cumsum(row_fill)
        stair_side = np.cumsum(row_side)
        cells = []
        for m in marks:
            stair_gap = -stair_fill[m - 1] / stair_side[m - 1]
            row_gap = -row_fill[m - 1] / row_side[m - 1]
            stair_err = m * (stair_gap - lead / (math.log(2 * m) - 1.5))
            row_err = m * (row_gap - lead / (math.log(2 * m) - 1))
            cells.append(
                f"n={m}: {stair_gap:.6f} ({stair_err:+.2f}), row {row_gap:.6f} ({row_err:+.2f})"
            )
        print(f"  dim {dim} code {code}:")
        for cell in cells:
            print(f"    {cell}")


# MAIN


def main():
    clock = time.time()

    def lap(name):
        nonlocal clock
        print(f"[{name} {time.time() - clock:.1f}s]")
        clock = time.time()

    section_render()
    lap("render")
    byprofile = {}
    for dim in (1, 2, 3):
        for a in profiles(dim):
            byprofile[a] = roots(a)
    table = list(byprofile.items())
    print(
        f"roots: {len(table)} profiles (3 at dim 1, 11 at dim 2, 63 at dim 3) cover all {3 + 15 + 255} nonempty codes"
    )
    lap("roots")
    section_exact(table)
    lap("exact")
    section_asymptotic(table)
    lap("asymptotic")
    section_named(byprofile)
    lap("named")
    section_pairing(byprofile)
    lap("pairing")
    section_parity()
    lap("parity")
    section_family()
    lap("family")
    section_fence(byprofile)
    lap("fence")
    section_staircase(byprofile)
    lap("staircase")


if __name__ == "__main__":
    main()
