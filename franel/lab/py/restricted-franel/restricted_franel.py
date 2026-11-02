import subprocess
import sys
import time
from fractions import Fraction
from math import gcd, log, pi

import numpy as np

BASE3 = (3, frozenset({0, 1}), "base 3, digits {0,1}")
BASE10 = (10, frozenset(range(9)), "base 10, digit 9 missing")

TWO_PI_I = 2j * pi

# THE DESIGN


def holds(base, digits, n):
    if n == 0:
        return True
    while n > 0:
        if n % base not in digits:
            return False
        n //= base
    return True


def flags(base, digits, q):
    keep = np.zeros(q + 1, dtype=bool)
    for n in range(q + 1):
        keep[n] = holds(base, digits, n)
    return keep


def mobius_sieve(n):
    mu = np.ones(n + 1, dtype=np.int8)
    rest = np.arange(n + 1, dtype=np.int64 if n >= 2**31 else np.int32)
    root = int(n**0.5) + 1
    small = np.ones(root + 1, dtype=bool)
    small[:2] = False
    for p in range(2, root + 1):
        if not small[p]:
            continue
        small[p * p :: p] = False
        mu[p::p] *= -1
        mu[p * p :: p * p] = 0
        rest[p::p] //= p
    mu[rest > 1] *= -1
    mu[0] = 0
    return mu


def totients(n):
    phi = np.arange(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if phi[p] == p:
            phi[p::p] -= phi[p::p] // p
    return phi


# THE DILATED MERTENS SUMS


def mertens_dilated(keep, mu, q, d):
    c = np.arange(1, q // d + 1, dtype=np.int64)
    return int(mu[c][keep[d * c]].sum())


def mertens_design(keep, mu, q):
    return mertens_dilated(keep, mu, q, 1)


# THE LITERAL FOURIER SIDE, DENOMINATOR SET


def literal_denominator(keep, q, freqs):
    out = {m: 0j for m in freqs}
    for b in range(1, q + 1):
        if not keep[b]:
            continue
        a = np.arange(1, b + 1, dtype=np.int64)
        a = a[np.gcd(a, b) == 1]
        for m in freqs:
            r = (m * a) % b
            out[m] += complex(np.exp(TWO_PI_I * (r / b)).sum())
    return out


def card_denominator(keep, phi, q):
    return int(phi[1 : q + 1][keep[1 : q + 1]].sum())


def ladder_check(keep, mu, q):
    total = 0j
    worst = 0.0
    at = 0
    bad = 0
    for b in range(1, q + 1):
        if not keep[b]:
            continue
        a = np.arange(1, b + 1, dtype=np.int64)
        a = a[np.gcd(a, b) == 1]
        v = complex(np.exp(TWO_PI_I * (a % b) / b).sum())
        total += v
        e = abs(v - int(mu[b]))
        if e > worst:
            worst, at = e, b
        if round(v.real) != int(mu[b]) or abs(v.imag) > 0.5:
            bad += 1
    return total, worst, at, bad


# THE STRICT SET


def design_list(keep, q):
    return np.nonzero(keep[1 : q + 1])[0] + 1


def strict_literal(keep, q):
    sf = design_list(keep, q)
    total = 0j
    card = 0
    for i, b in enumerate(sf):
        b = int(b)
        a = sf[: i + 1]
        a = a[np.gcd(a, b) == 1]
        card += a.size
        total += complex(np.exp(TWO_PI_I * (a % b) / b).sum())
    return total, card


def strict_ramanujan_literal(keep, b):
    a = np.arange(1, b + 1, dtype=np.int64)
    a = a[keep[a] & (np.gcd(a, b) == 1)]
    return complex(np.exp(TWO_PI_I * (a % b) / b).sum())


def strict_ramanujan_divisor(keep, mu, b):
    total = 0j
    for d in range(1, b + 1):
        if b % d or mu[d] == 0:
            continue
        n = b // d
        ap = np.arange(1, n + 1, dtype=np.int64)
        ap = ap[keep[d * ap]]
        total += int(mu[d]) * complex(np.exp(TWO_PI_I * (ap % n) / n).sum())
    return total


def strict_weights(keep, mu, phi, q):
    sf = design_list(keep, q)
    plain = 0
    counted = 0
    ratio = 0.0
    for i, b in enumerate(sf):
        b = int(b)
        a = sf[: i + 1]
        pf = int((np.gcd(a, b) == 1).sum())
        plain += int(mu[b])
        counted += int(mu[b]) * pf
        ratio += int(mu[b]) * pf / int(phi[b])
    return plain, counted, ratio


# THE GCD KERNEL


def kernel_sum(keep, mu, q):
    ms = {}
    for d in range(1, q + 1):
        v = mertens_dilated(keep, mu, q, d)
        if v:
            ms[d] = v
    total = Fraction(0)
    ks = sorted(ms)
    for d in ks:
        for e in ks:
            total += Fraction(gcd(d, e) ** 2 * ms[d] * ms[e], d * e)
    return total, ms


def fourier_side(keep, q, cut):
    nodes = []
    for b in range(1, q + 1):
        if not keep[b]:
            continue
        for a in range(1, b + 1):
            if gcd(a, b) == 1:
                nodes.append((a, b))
    m = len(nodes)
    num = np.array([a for a, _ in nodes], dtype=np.int64)
    den = np.array([b for _, b in nodes], dtype=np.int64)
    total = 0.0
    step = max(1, 2_000_000 // max(m, 1))
    k = 1
    while k <= cut:
        hi = min(cut, k + step - 1)
        ks = np.arange(k, hi + 1, dtype=np.int64)
        r = (ks[:, None] * num[None, :]) % den[None, :]
        s = np.exp(TWO_PI_I * (r / den[None, :])).sum(axis=1)
        total += float((np.abs(s) ** 2 / ks.astype(float) ** 2).sum())
        k = hi + 1
    return 2.0 * total, m


# THE BACKWARD SCAN


def kernel_float(keep, mu, w, q):
    ms = np.zeros(q + 1)
    for d in range(1, q + 1):
        ms[d] = mertens_dilated(keep, mu, q, d)
    v = ms[1 : q + 1]
    return float(v @ w[:q, :q] @ v)


def gcd_weights(qmax):
    d = np.arange(1, qmax + 1, dtype=np.int64)
    g = np.gcd.outer(d, d).astype(np.float64)
    return g * g / np.outer(d.astype(np.float64), d.astype(np.float64))


def scan_backward(design, qmax):
    keep = keep_of(design, qmax)
    mu = mobius_sieve(qmax)
    w = gcd_weights(qmax)
    best = 0.0
    at = 0
    tail = 0.0
    tail_at = 0
    seen = 0
    bad = 0
    for q in range(1, qmax + 1):
        if not keep[q]:
            continue
        seen += 1
        g = kernel_float(keep, mu, w, q)
        mf = mertens_design(keep, mu, q)
        ident = g * pi * pi / 3.0
        r = 2.0 * mf * mf / ident
        if r > 1.0 + 1e-9:
            bad += 1
        if r > best:
            best, at = r, q
        if q >= SCAN_FLOOR and r > tail:
            tail, tail_at = r, q
    return seen, bad, best, at, tail, tail_at


# THE CONTROL


def farey_delta_square(keep, q):
    nodes = []
    for b in range(1, q + 1):
        if not keep[b]:
            continue
        for a in range(1, b + 1):
            if gcd(a, b) == 1:
                nodes.append(Fraction(a, b))
    nodes.sort()
    m = len(nodes)
    s2 = sum((r - Fraction(j + 1, m)) ** 2 for j, r in enumerate(nodes))
    return m, s2


def reflection_closed(keep, q):
    for b in range(1, q + 1):
        if not keep[b]:
            continue
        for a in range(1, b):
            if gcd(a, b) == 1 and keep[a] != keep[b - a]:
                return False
    return True


# THE VERBS

FREQS = [1, 2, 3, 4, 5, 6, 12]

DEN_SMALL = [(BASE3, 2187), (BASE10, 1000), ((0, None, "full set, control"), 300)]

DEN_LADDER = [
    (BASE3, [243, 2187, 19683, 100000]),
    (BASE10, [100, 1000, 10000]),
    ((0, None, "full set, control"), [100, 1000, 10000]),
]


def keep_of(design, q):
    base, digits, _ = design
    if digits is None:
        k = np.ones(q + 1, dtype=bool)
        return k
    return flags(base, digits, q)


def verb_denominator():
    print("DENOMINATOR SET, THE FREQUENCY-m IDENTITY")
    for design, q in DEN_SMALL:
        keep = keep_of(design, q)
        mu = mobius_sieve(q)
        lit = literal_denominator(keep, q, FREQS)
        print(f"\n{design[2]}, Q = {q}")
        print(
            "    m | sum_(d|m) d M_F(Q/d; d) |    Re literal |    Im literal |       err"
        )
        for m in FREQS:
            exact = sum(
                d * mertens_dilated(keep, mu, q, d)
                for d in range(1, m + 1)
                if m % d == 0
            )
            v = lit[m]
            err = abs(v - exact)
            print(f"{m:5d} | {exact:23d} | {v.real:13.9f} | {v.imag:13.9f} | {err:.3e}")

    print("\nDENOMINATOR SET, FREQUENCY 1 AGAINST THE DESIGN MERTENS SUM")
    print(
        "   set                        |      Q |  M_F(Q) |    Re literal |    Im literal |  agg err |"
        " per-b err |  at b | bad |    s"
    )
    for design, ladder in DEN_LADDER:
        for q in ladder:
            t0 = time.time()
            keep = keep_of(design, q)
            mu = mobius_sieve(q)
            exact = mertens_design(keep, mu, q)
            v, worst, at, bad = ladder_check(keep, mu, q)
            dt = time.time() - t0
            print(
                f"   {design[2]:26s} | {q:6d} | {exact:7d} | {v.real:13.9f} | "
                f"{v.imag:13.9f} | {abs(v - exact):.2e} | {worst:9.2e} | {at:5d} | {bad:3d} | {dt:4.1f}"
            )


SCAN_FLOOR = 100

SCAN = [(BASE3, 2187), (BASE10, 400), ((0, None, "full set, control"), 400)]

ID_CASES = [
    (BASE3, 81, 200000),
    (BASE3, 243, 200000),
    (BASE10, 40, 200000),
    ((0, None, "full set, control"), 40, 200000),
]


def verb_identity():
    print(
        "\nTHE DIGIT-RESTRICTED FRANEL IDENTITY, GCD-KERNEL SIDE AGAINST THE FOURIER SIDE"
    )
    print(
        "   set                        |    Q |    m |      G_F(Q) | (pi^2/3) G_F |  truncated |      gap |  tail bound"
    )
    for design, q, cut in ID_CASES:
        keep = keep_of(design, q)
        mu = mobius_sieve(q)
        g, _ = kernel_sum(keep, mu, q)
        ident = float(g) * pi * pi / 3.0
        lhs, m = fourier_side(keep, q, cut)
        tail = 2.0 * m * m / cut
        print(
            f"   {design[2]:26s} | {q:4d} | {m:4d} | {float(g):11.6f} | {ident:12.6f} | "
            f"{lhs:10.6f} | {abs(ident - lhs):.2e} | {tail:.4f}"
        )

    print("\n   THE RANK FORM, G_F(Q) - 1 = 12 m sum delta^2, exact rationals")
    print(
        "   set                        |    Q |    m |   sum delta^2 |  G_F - 1 == 12 m sum delta^2"
    )
    for design, q in [
        (BASE3, 81),
        (BASE3, 243),
        (BASE10, 40),
        ((0, None, "full set, control"), 40),
    ]:
        keep = keep_of(design, q)
        mu = mobius_sieve(q)
        g, _ = kernel_sum(keep, mu, q)
        m, s2 = farey_delta_square(keep, q)
        print(
            f"   {design[2]:26s} | {q:4d} | {m:4d} | {float(s2):13.10f} | {g - 1 == 12 * m * s2}"
        )

    print("\n   the backward inequality 2 M_F(Q)^2 <= (pi^2/3) G_F(Q), sampled")
    print("   set                        |      Q |  M_F(Q) |   (pi^2/3) G_F |  ratio")
    for design, ladder in [
        (BASE3, [243, 2187]),
        (BASE10, [100, 1000]),
        ((0, None, "full set, control"), [100, 1000]),
    ]:
        for q in ladder:
            keep = keep_of(design, q)
            mu = mobius_sieve(q)
            g, _ = kernel_sum(keep, mu, q)
            mf = mertens_design(keep, mu, q)
            ident = float(g) * pi * pi / 3.0
            print(
                f"   {design[2]:26s} | {q:6d} | {mf:7d} | {ident:14.6f} | {2.0 * mf * mf / ident:.6f}"
            )

    print("\n   the same inequality asserted at EVERY integer Q up to the bound")
    print("   both sides step only at Q in S_F, so scanning S_F covers every integer Q")
    print(
        "   set                        |   Qmax | jumps | violations |  max ratio |  at Q |  max over Q >= 100 |  at Q"
    )
    for design, qmax in SCAN:
        seen, bad, best, at, tail, tail_at = scan_backward(design, qmax)
        assert bad == 0
        print(
            f"   {design[2]:26s} | {qmax:6d} | {seen:5d} | {bad:10d} | {best:10.6f} | {at:5d} | "
            f"{tail:18.6f} | {tail_at:5d}"
        )


STRICT_LADDER = [(BASE3, [243, 2187, 19683, 177147]), (BASE10, [100, 1000, 10000])]


def verb_strict():
    print("\nTHE STRICT SET, THE DIVISOR-SUM IDENTITY")
    print("   set                        |    Q | denominators |  max abs diff")
    for design, q in [(BASE3, 729), (BASE10, 200)]:
        keep = keep_of(design, q)
        mu = mobius_sieve(q)
        worst = 0.0
        n = 0
        for b in range(1, q + 1):
            if not keep[b]:
                continue
            n += 1
            worst = max(
                worst,
                abs(
                    strict_ramanujan_literal(keep, b)
                    - strict_ramanujan_divisor(keep, mu, b)
                ),
            )
        print(f"   {design[2]:26s} | {q:4d} | {n:12d} | {worst:.3e}")

    print("\nTHE STRICT SET, THE FIRST Q WITH A NONZERO IMAGINARY PART")
    for design, q in [(BASE3, 100), (BASE10, 100)]:
        keep = keep_of(design, q)
        run = 0j
        hit = None
        for b in range(1, q + 1):
            if not keep[b]:
                continue
            run += strict_ramanujan_literal(keep, b)
            if hit is None and abs(run.imag) > 1e-9:
                hit = (b, run)
        b, run = hit
        print(f"   {design[2]:26s}: Q = {b}, sum = {run.real:.9f} + {run.imag:.9f} i")

    print("\nTHE STRICT SET, FREQUENCY 1 AGAINST EVERY MERTENS-TYPE CANDIDATE")
    print(
        "   set                        |      Q |     card |     Re T |     Im T |     abs T | abs T/card |"
        "  M_F(Q) |  sum mu phi_F |  sum mu phi_F/phi"
    )
    for design, ladder in STRICT_LADDER:
        for q in ladder:
            keep = keep_of(design, q)
            mu = mobius_sieve(q)
            phi = totients(q)
            t, card = strict_literal(keep, q)
            plain, counted, ratio = strict_weights(keep, mu, phi, q)
            print(
                f"   {design[2]:26s} | {q:6d} | {card:8d} | {t.real:8.3f} | {t.imag:8.3f} | "
                f"{abs(t):9.3f} | {abs(t) / card:10.6f} | {plain:7d} | {counted:13d} | {ratio:17.6f}"
            )


# THE DILATE AUTOMATON


def dilate_matrix(base, digits, d):
    t = [[0] * d for _ in range(d)]
    for r in range(d):
        for e in range(base):
            v = d * e + r
            if v % base in digits:
                t[r][v // base] += 1
    return t


def dilate_accept(base, digits, d):
    return [1 if (r == 0 or holds(base, digits, r)) else 0 for r in range(d)]


def qfree(d, base):
    while True:
        g = gcd(d, base)
        if g == 1:
            return d
        d //= g


def automaton_count(base, digits, d, power):
    t = dilate_matrix(base, digits, d)
    v = [0] * d
    v[0] = 1
    for _ in range(power):
        w = [0] * d
        for i, x in enumerate(v):
            if x:
                for j, y in enumerate(t[i]):
                    if y:
                        w[j] += x * y
        v = w
    return sum(x * y for x, y in zip(v, dilate_accept(base, digits, d)))


def dilate_mask(base, digits, d, x):
    lut = np.zeros(base, dtype=bool)
    for f in digits:
        lut[f] = True
    t = d * np.arange(1, x + 1, dtype=np.int64)
    ok = np.ones(x, dtype=bool)
    while t.any():
        ok &= lut[t % base] | (t == 0)
        t //= base
    return ok


def divisor_counts(base, digits, cut):
    m = np.nonzero(flags(base, digits, cut))[0]
    m = m[m > 0]
    n = np.zeros(cut + 1, dtype=np.int64)
    for d in range(1, cut + 1):
        n[d] = int((m % d == 0).sum())
    return n


def smith_bilinear(n, cut, block=400):
    u = np.sqrt(n[1 : cut + 1].astype(np.float64))
    idx = np.arange(1, cut + 1, dtype=np.int64)
    total = 0.0
    for lo in range(0, cut, block):
        hi = min(lo + block, cut)
        g = np.gcd.outer(idx[lo:hi], idx).astype(np.float64)
        w = g * g / (idx[lo:hi, None] * idx[None, :])
        total += float((w * u[lo:hi, None] * u[None, :]).sum())
    return total


DILATE_D3 = [1, 2, 3, 4, 5, 7, 8, 11, 13, 16, 22, 31]
DILATE_D10 = [1, 2, 3, 4, 5, 7, 11, 13, 121, 243, 729]


def verb_dilate():
    base, digits, name = BASE3
    print(f"\nTHE DILATE AUTOMATON, {name}, d = 2, states are the carries")
    t = dilate_matrix(base, digits, 2)
    acc = dilate_accept(base, digits, 2)
    for r, row in enumerate(t):
        print(f"   carry {r} -> {row}   accept {acc[r]}")
    print(
        f"   row sums {[sum(r) for r in t]}, column sums {[sum(c) for c in zip(*t)]}, |F| = {len(digits)}"
    )
    print("\n    L |    3^L | transfer matrix | brute force |    2^L - 1")
    for lvl in range(1, 13):
        x = base**lvl
        auto = automaton_count(base, digits, 2, lvl) - 1
        brute = int(dilate_mask(base, digits, 2, x - 1).sum())
        print(
            f"   {lvl:2d} | {x:6d} | {auto:15d} | {brute:11d} | {2**lvl - 1:10d}  {'ok' if auto == brute else 'MISMATCH'}"
        )

    print("\nTHE COLUMN-SUM LAW, every dilate d <= 64, both designs")
    for bb, dd, nm in (BASE3, BASE10):
        bad, off = 0, []
        for d in range(1, 65):
            m = dilate_matrix(bb, dd, d)
            if any(sum(c) != len(dd) for c in zip(*m)):
                bad += 1
            if any(sum(r) != len(dd) for r in m):
                off.append(d)
        print(
            f"   {nm:26s} columns off |F| at {bad} of 64; rows off |F| at {len(off)} of 64, every one with gcd(d,q) > 1: {all(gcd(d, bb) > 1 for d in off)}"
        )

    for bb, dd, nm, ds in (
        (BASE3[0], BASE3[1], BASE3[2], DILATE_D3),
        (BASE10[0], BASE10[1], BASE10[2], DILATE_D10),
    ):
        al = np.log(len(dd)) / np.log(bb)
        print(f"\nTHE ACCEPT MASS, {nm}, K_d = A_d(q^24)/|F|^24 against |Acc_d|/d")
        print("      d | gcd(d,q) |     K_d | |Acc_d|/d | d^(alpha-1)")
        for d in ds + ([8, 16, 32] if bb == 10 else []):
            k = automaton_count(bb, dd, d, 24) / len(dd) ** 24
            a = sum(dilate_accept(bb, dd, d))
            print(
                f"   {d:6d} | {gcd(d, bb):8d} | {k:7.4f} | {a / d:9.6f} | {d ** (al - 1):11.6f}"
            )

    for bb, dd, nm, top, ds in (
        (3, BASE3[1], BASE3[2], 12, DILATE_D3),
        (10, BASE10[1], BASE10[2], 7, DILATE_D10[:7]),
    ):
        x = bb**top
        mu = mobius_sieve(x)
        al = np.log(len(dd)) / np.log(bb)
        root = x ** (al / 2)
        print(f"\nTHE DILATED METER, {nm}, x = {bb}^{top} = {x}")
        print(
            "      d |  A_d(x) |  M_F(x;d) | max |M| | log max/log x | same at x/q | (U) ratio | (U\u0027) ratio |  local exponents"
        )
        for d in ds:
            ok = dilate_mask(bb, dd, d, x)
            run = np.cumsum(np.where(ok, mu[1 : x + 1], 0))
            peaks = [
                int(np.abs(run[: bb**lvl]).max()) for lvl in range(top - 4, top + 1)
            ]
            loc = [
                round(float(np.log(peaks[i + 1] / peaks[i]) / np.log(bb)), 3)
                for i in range(len(peaks) - 1)
            ]
            mass, peak = int(ok.sum()), peaks[-1]
            ru = peak / (d ** ((al - 1) / 2) * root)
            rv = peak / (qfree(d, bb) ** ((al - 1) / 2) * root)
            under = np.log(peaks[-2]) / np.log(x / bb)
            print(
                f"   {d:6d} | {mass:7d} | {int(run[-1]):9d} | {peak:7d} | {np.log(peak) / np.log(x):13.6f} | "
                f"{under:11.6f} | {ru:9.4f} | {rv:10.4f} | {loc}"
            )

    print(f"\nTHE LADDER, {BASE3[2]}, the meter against its dilate at d = 2")
    x = 3**12
    mu = mobius_sieve(x)
    r1 = np.cumsum(np.where(dilate_mask(3, BASE3[1], 1, x), mu[1 : x + 1], 0))
    r2 = np.cumsum(np.where(dilate_mask(3, BASE3[1], 2, x), mu[1 : x + 1], 0))
    print("    L |  M_F(3^L) |  M_F(3^L;2) | max |M_F| | max |M_F(;2)|")
    for lvl in range(1, 13):
        n = 3**lvl
        print(
            f"   {lvl:2d} | {int(r1[n - 1]):9d} | {int(r2[n - 1]):11d} | {int(np.abs(r1[:n]).max()):9d} | "
            f"{int(np.abs(r2[:n]).max()):13d}"
        )

    al = np.log(2) / np.log(3)
    print(
        f"\nTHE SMITH REDUCTION, B(Q) = sum gcd(d,e)^2/(d e) sqrt(N_F(Q;d) N_F(Q;e)), {BASE3[2]}"
    )
    print(
        "        Q |  A_F(Q) |        B(Q) | B/Q^alpha | B/(Q^alpha ln Q) | B/(Q^alpha ln^2 Q) | local exponent"
    )
    prev = None
    for lvl in range(4, 9):
        cut = 3**lvl
        b = smith_bilinear(divisor_counts(3, BASE3[1], cut), cut)
        qa = cut**al
        loc = "" if prev is None else f"{np.log(b / prev) / np.log(3):14.3f}"
        print(
            f"   {cut:8d} | {2**lvl:7d} | {b:11.3f} | {b / qa:9.4f} | {b / (qa * np.log(cut)):16.4f} | "
            f"{b / (qa * np.log(cut) ** 2):18.4f} | {loc}"
        )
        prev = b


# THE CONVERSE


def digit_gap(digits):
    f = sorted(digits)
    g = 0
    for x in f[1:]:
        g = gcd(g, x - f[0])
    return g


def brute_dilate(base, digits, d, power):
    return sum(1 for c in range(1, base**power) if holds(base, digits, d * c))


def positive_count(base, digits, d, power):
    return automaton_count(base, digits, d, power) - (1 if 0 in digits else 0)


def subsets(base):
    for m in range(1, 1 << base):
        yield {e for e in range(base) if m >> e & 1}


SCOPE_BASES = (3, 4, 5)
GAP_BASES = range(2, 8)
GAP_POWER = 400


def verb_converse():
    print(
        "\nTHE SCOPE OF D1, D2, D5, automaton against brute force, q = 3, 4, 5, every F, d <= 6, L <= 5"
    )
    tally = {True: [0, 0], False: [0, 0]}
    for base in SCOPE_BASES:
        for digits in subsets(base):
            for d in range(1, 7):
                for lvl in range(1, 6):
                    row = tally[0 in digits]
                    row[0] += 1
                    row[1] += positive_count(base, digits, d, lvl) != brute_dilate(
                        base, digits, d, lvl
                    )
    for key, label in ((True, "0 in F"), (False, "0 not in F")):
        n, bad = tally[key]
        print(f"   {label:12s} {bad:4d} mismatches of {n:4d}")
    base, digits, d, lvl = 3, {1, 2}, 1, 3
    print(
        f"   witness q = 3, F = {{1,2}}, d = 1, L = 3: automaton {positive_count(base, digits, d, lvl)}, "
        f"true {brute_dilate(base, digits, d, lvl)}, |F|^L {len(digits) ** lvl}"
    )

    print(
        f"\nTHE MASS CONSTANT, K_d against |Acc_d|/d at L = {GAP_POWER}, every F with 0 in F, gcd(d,q) = 1, 2 <= d <= 24"
    )
    split = {True: [0, 0], False: [0, 0]}
    odd = []
    for base in GAP_BASES:
        for digits in subsets(base):
            if 0 not in digits or len(digits) < 2:
                continue
            gap = digit_gap(digits)
            for d in range(2, 25):
                if gcd(d, base) != 1:
                    continue
                k = (
                    automaton_count(base, digits, d, GAP_POWER)
                    / len(digits) ** GAP_POWER
                )
                row = split[gcd(d, gap) == 1]
                row[0] += 1
                good = abs(k - sum(dilate_accept(base, digits, d)) / d) < 1e-9
                row[1] += good
                if good == (gcd(d, gap) > 1):
                    odd.append(
                        (
                            base,
                            sorted(digits),
                            d,
                            round(k, 6),
                            sum(dilate_accept(base, digits, d)) / d,
                        )
                    )
    for key, label in ((True, "gcd(d, gap) = 1"), (False, "gcd(d, gap) > 1")):
        n, ok = split[key]
        print(f"   {label:16s} {ok:5d} agree, {n - ok:5d} fail, of {n:5d}")
    print(f"   off the split: {odd}")
    base, digits, d = 3, {0, 2}, 2
    print(
        f"   witness q = 3, F = {{0,2}}, gap 2, d = 2: counts "
        f"{[brute_dilate(base, digits, d, lvl) + 1 for lvl in range(1, 9)]}, |F|^L, so K_2 = 1 against "
        f"|Acc_2|/2 = {sum(dilate_accept(base, digits, d)) / d}"
    )

    base, digits, name = BASE3
    al = np.log(len(digits)) / np.log(base)
    top = 12
    x = base**top
    mu = mobius_sieve(x)
    run = np.cumsum(np.where(dilate_mask(base, digits, 1, x), mu[1 : x + 1], 0))
    peak, tail = int(np.abs(run).max()), int(run[-1])
    print(
        f"\n(U) IS UNSATISFIABLE, {name}, x = {base}^{top}, M_F(x; {base}^j) = M_F(x) = {tail} at every j by D4"
    )
    print(
        "      j |        d = q^j | d^((alpha-1)/2) x^(alpha/2) | |M_F| / bound | max |M| / bound"
    )
    for j in range(top + 1):
        d = base**j
        bound = d ** ((al - 1) / 2) * x ** (al / 2)
        print(
            f"   {j:4d} | {d:14d} | {bound:27.4f} | {abs(tail) / bound:13.3f} | {peak / bound:15.3f}"
        )

    base, digits, name = BASE10
    print(
        f"\nTHE BASE-SMOOTH MASS, {name}, K_d at L = 24 over the base-smooth d <= 1000"
    )
    smooth = sorted(d for d in range(1, 1001) if qfree(d, base) == 1)
    masses = {
        d: automaton_count(base, digits, d, 24) / len(digits) ** 24 for d in smooth
    }
    lo = min(masses, key=masses.get)
    hi = max(masses, key=masses.get)
    print(
        f"   {len(smooth)} such d; K_d in [{masses[lo]:.4f} at d = {lo}, {masses[hi]:.4f} at d = {hi}]"
    )
    print(
        f"   sample {[(d, round(masses[d], 4)) for d in (2, 4, 5, 8, 16, 32, 64, 128, 256, 512, 625)]}"
    )
    ratio = {}
    for d in range(1, 201):
        f = qfree(d, base)
        k = automaton_count(base, digits, d, 24) / len(digits) ** 24
        kf = masses.get(f) or automaton_count(base, digits, f, 24) / len(digits) ** 24
        masses[f] = kf
        ratio[d] = k / kf
    lo = min(ratio, key=ratio.get)
    hi = max(ratio, key=ratio.get)
    print(
        f"   K_d / K_(q-free part) over d <= 200 in [{ratio[lo]:.4f} at d = {lo}, {ratio[hi]:.4f} at d = {hi}]"
    )


# THE SANDWICH


def jordan_two(n):
    j = np.arange(n + 1, dtype=np.int64) ** 2
    for p in range(2, n + 1):
        if j[p] == p * p:
            j[p::p] -= j[p::p] // (p * p)
    return j


def sigma_ratio_max(n):
    s = np.zeros(n + 1)
    for k in range(1, n + 1):
        s[k::k] += 1.0 / k
    return float(s[1:].max()), int(s[1:].argmax()) + 1


def dilate_vector(keep, mu, q):
    x = np.zeros(q + 1, dtype=np.int64)
    for d in range(1, q + 1):
        x[d] = mertens_dilated(keep, mu, q, d)
    return x


def harmonic_exact(x, q):
    y = [Fraction(0)] * (q + 1)
    for f in range(1, q + 1):
        y[f] = sum(Fraction(int(x[f * m]), m) for m in range(1, q // f + 1))
    return y


def jordan_form(x, q):
    j2 = jordan_two(q)
    y = harmonic_exact(x, q)
    return sum(Fraction(int(j2[f]), f * f) * y[f] * y[f] for f in range(1, q + 1)), y


def kernel_double(x, q, block=1024):
    idx = np.arange(1, q + 1, dtype=np.int64)
    v = x[1 : q + 1].astype(np.float64)
    total = 0.0
    for lo in range(0, q, block):
        hi = min(lo + block, q)
        g = np.gcd.outer(idx[lo:hi], idx).astype(np.float64)
        w = g * g / (idx[lo:hi, None] * idx[None, :])
        total += float(v[lo:hi] @ w @ v)
    return total


def divisor_maps(q, mu):
    from scipy.sparse import csr_matrix

    rows, cols, harm, mob = [], [], [], []
    for f in range(1, q + 1):
        m = np.arange(1, q // f + 1, dtype=np.int64)
        rows.append(np.full(m.size, f - 1, dtype=np.int64))
        cols.append(f * m - 1)
        harm.append(1.0 / m)
        mob.append(mu[m] / m)
    r, c = np.concatenate(rows), np.concatenate(cols)
    a = csr_matrix((np.concatenate(harm), (r, c)), shape=(q, q))
    b = csr_matrix((np.concatenate(mob), (r, c)), shape=(q, q))
    return a, b


def top_singular(a, steps=3000, seed=1):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(a.shape[1])
    v /= np.linalg.norm(v)
    at = a.T.tocsr()
    for _ in range(steps):
        w = at @ (a @ v)
        n = np.linalg.norm(w)
        if n == 0:
            return 0.0
        v = w / n
    return float(np.linalg.norm(a @ v))


SANDWICH = [
    (BASE3, [81, 243, 729, 2187, 6561]),
    ((0, None, "full set, control"), [40, 81, 243]),
]

SANDWICH_EXACT = 729

RANK_CASES = [
    (3, {0, 1}, 40),
    (10, set(range(9)), 40),
    (3, {0, 2}, 40),
    (4, {0, 2, 3}, 40),
    (5, {0, 2, 4}, 40),
    (10, {0, 2, 5, 7}, 40),
]


def verb_sandwich():
    print(
        "\nTHE SANDWICH, G_F(Q) = sum_f (J_2(f)/f^2) y_f^2 with y_f = sum_m x_(fm)/m, x_d = M_F(Q/d; d)"
    )
    print(
        "   6/pi^2 = %.6f; exact means the Jordan form and the gcd double sum agree as rationals"
        % (6 / pi**2)
    )
    print(
        "   set                        |     Q |      G_F(Q) | sum x^2 | sum y^2 | G/sum y^2 | G/sum x^2 |"
        " x_1^2/sum x^2 |   exact | double gap"
    )
    rows = []
    for design, ladder in SANDWICH:
        for q in ladder:
            t0 = time.time()
            keep = keep_of(design, q)
            mu = mobius_sieve(q)
            x = dilate_vector(keep, mu, q)
            g, y = jordan_form(x, q)
            sx = int((x[1:] ** 2).sum())
            sy = float(sum(v * v for v in y[1:]))
            exact = "-"
            if q <= SANDWICH_EXACT:
                gd, _ = kernel_sum(keep, mu, q)
                exact = str(gd == g)
            gap = abs(float(g) - kernel_double(x, q))
            rows.append((design, q, mu, x, float(g), sx, sy))
            print(
                f"   {design[2]:26s} | {q:5d} | {float(g):11.6f} | {sx:7d} | {sy:7.3f} | {float(g) / sy:9.6f} | "
                f"{float(g) / sx:9.6f} | {int(x[1]) ** 2 / sx:13.6f} | {exact:>7s} | {gap:.1e}   {time.time() - t0:4.1f} s"
            )

    print(
        "\n   THE TWO MAPS, l^2 operator norms on vectors indexed by d <= Q, against the Schur bound N_Q^(1/2)"
    )
    print(
        "   N_Q = H_Q max_(d<=Q) sigma(d)/d; the power readings are Rayleigh quotients, lower bounds on the norms"
    )
    print(
        "   set                        |     Q |    H_Q | max sigma/d | at d |    N_Q | (6/pi^2)/N_Q |  N_Q^(1/2) | 1 + ln Q |"
        " norm x->y | norm y->x | sqrt(sum y^2/sum x^2)"
    )
    for design, q, mu, x, g, sx, sy in rows:
        a, b = divisor_maps(q, mu)
        na, nb = top_singular(a), top_singular(b)
        hq = float((1.0 / np.arange(1, q + 1)).sum())
        sm, at = sigma_ratio_max(q)
        nq = hq * sm
        print(
            f"   {design[2]:26s} | {q:5d} | {hq:6.4f} | {sm:11.6f} | {at:4d} | {nq:6.2f} | {6 / pi**2 / nq:12.4f} | "
            f"{np.sqrt(nq):10.6f} | {1 + np.log(q):8.4f} | {na:9.6f} | {nb:9.6f} | {np.sqrt(sy / sx):9.6f}"
        )

    print(
        "\n   THE RANK FORM'S CONSTANT, G_F(Q) - 12 m sum delta^2 as an exact rational at every Q in S_F up to the bound"
    )
    print(
        "   base | digits          | 1 in F |  Qmax | jumps | constants | G_F(26) | 12 m S2 at 26"
    )
    for base, digits, qmax in RANK_CASES:
        keep = flags(base, frozenset(digits), qmax)
        mu = mobius_sieve(qmax)
        consts = set()
        at26 = ""
        for q in range(1, qmax + 1):
            if not keep[q]:
                continue
            g, _ = kernel_sum(keep, mu, q)
            m, s2 = farey_delta_square(keep, q)
            consts.add(g - 12 * m * s2)
            if q == 26:
                at26 = f"{float(g):.6f} | {float(12 * m * s2):.6f}"
        print(
            f"   {base:4d} | {sorted(digits)!s:15s} | {1 in digits!s:6s} | {qmax:5d} | {len(consts):5d} | {sorted(consts)} | {at26}"
        )


# THE ACCEPTING LAW


def padded_strings(base, digits, level):
    s = np.zeros(1, dtype=np.int64)
    for j in range(level):
        s = np.concatenate([s + f * base**j for f in sorted(digits)])
    return np.sort(s)


def design_below(base, digits, level):
    s = padded_strings(base, digits, level)
    s = s[s > 0]
    if 1 in digits:
        s = np.append(s, base**level)
    return s


def divisor_table(members, cut):
    n = np.zeros(cut + 1, dtype=np.int64)
    for m in members.tolist():
        r = int(m**0.5)
        while r * r > m:
            r -= 1
        while (r + 1) * (r + 1) <= m:
            r += 1
        for d in range(1, r + 1):
            if m % d == 0:
                n[d] += 1
                if m // d != d:
                    n[m // d] += 1
    return n


def count_multiples(members, cut):
    n = np.zeros(cut + 1, dtype=np.int64)
    for lo in range(1, cut + 1, 64):
        ds = np.arange(lo, min(lo + 64, cut + 1), dtype=np.int64)
        n[lo : lo + len(ds)] = (members[None, :] % ds[:, None] == 0).sum(axis=1)
    return n


def smith_jordan(n, cut):
    x = np.sqrt(n[: cut + 1].astype(np.float64))
    j2 = jordan_two(cut).astype(np.float64)
    total = 0.0
    for f in np.nonzero(n[1 : cut + 1])[0] + 1:
        f = int(f)
        y = float((x[f::f] / np.arange(1, cut // f + 1, dtype=np.float64)).sum())
        total += j2[f] / (f * f) * y * y
    return total


def sigma_over_d(cut):
    s = np.zeros(cut + 1)
    for k in range(1, cut + 1):
        s[k::k] += 1.0 / k
    return s


def roots_l1(base, digits, level, d):
    a = np.arange(1, d, dtype=np.int64)
    total = np.ones(d - 1)
    for j in range(level):
        t = (a * pow(base, j, d) % d) / d
        g = np.zeros(d - 1, dtype=np.complex128)
        for f in digits:
            g += np.exp(TWO_PI_I * f * t)
        total *= np.abs(g)
    return float(total.sum())


def mult_order(base, d):
    if gcd(base, d) != 1:
        return 0
    k, v = 1, base % d
    while v != 1:
        v = v * base % d
        k += 1
    return k


def pair_carry_matrix(base, digits, r):
    dist = {0: 1}
    for _ in range(r):
        nxt = {}
        for v, c in dist.items():
            for f in digits:
                nxt[v + f] = nxt.get(v + f, 0) + c
        dist = nxt
    top = max(dist)
    bound = top // (base - 1) + 1
    states = list(range(-bound, bound + 1))
    m = {c: {} for c in states}
    for c in states:
        for u, cu in dist.items():
            for v, cv in dist.items():
                w = u - v + c
                if w % base == 0 and w // base in m:
                    m[c][w // base] = m[c].get(w // base, 0) + cu * cv
    reach = {0}
    frontier = [0]
    while frontier:
        c = frontier.pop()
        for c2 in m[c]:
            if c2 not in reach:
                reach.add(c2)
                frontier.append(c2)
    back = {0}
    frontier = [0]
    while frontier:
        c = frontier.pop()
        for c0 in states:
            if c in m[c0] and c0 not in back:
                back.add(c0)
                frontier.append(c0)
    cls = sorted(reach & back)
    mat = np.array([[m[c].get(c2, 0) for c2 in cls] for c in cls], dtype=np.float64)
    return cls, mat


def pair_count(base, digits, r, t):
    cls, mat = pair_carry_matrix(base, digits, r)
    i = cls.index(0)
    v = np.zeros(len(cls))
    v[i] = 1
    exact = [[Fraction(int(x)) for x in row] for row in mat]
    vec = [Fraction(0)] * len(cls)
    vec[i] = Fraction(1)
    for _ in range(t):
        vec = [
            sum(vec[a] * exact[a][b] for a in range(len(cls))) for b in range(len(cls))
        ]
    return int(vec[i])


def perron_certificate(mat, floor, steps=2000):
    v = np.ones(len(mat))
    for _ in range(steps):
        v = mat @ v
        v /= v.max()
    rho = float((mat @ v).max() / v.max())
    vq = [Fraction(round(x * 10**9), 10**9) for x in v]
    exact = [[Fraction(int(x)) for x in row] for row in mat]
    quot = [
        sum(exact[a][b] * vq[b] for b in range(len(mat))) / vq[a]
        for a in range(len(mat))
    ]
    return rho, min(quot), min(quot) > floor


def verb_accepting():
    base, digits, name = BASE3
    fill = len(digits)
    al = np.log(fill) / np.log(base)
    print(
        f"\nTHE SURROGATE, {name}: sum tau(m) <= B(Q) <= (1 + ln Q) sum_d (sigma(d)/d) N_F(Q;d) <= (1 + ln Q)^2 sum tau(m)"
    )
    print(
        "    L |      Q | A_F(Q) | sum tau |        B(Q) | B/sum tau | amgm bound | bound/B | (1+lnQ)^2 sum tau | B/(Q^a ln^2 Q) | double gap"
    )
    for lvl in range(4, 13):
        t0 = time.time()
        cut = base**lvl
        members = design_below(base, digits, lvl)
        n = divisor_table(members, cut)
        tau = int(n.sum())
        b = smith_jordan(n, cut)
        s = sigma_over_d(cut)
        amgm = (1 + np.log(cut)) * float((s[1:] * n[1:]).sum())
        top = (1 + np.log(cut)) ** 2 * tau
        gap = ""
        if lvl <= 8:
            gap = f"{abs(b - smith_bilinear(n, cut)):.1e}"
        print(
            f"   {lvl:2d} | {cut:6d} | {len(members):6d} | {tau:7d} | {b:11.3f} | {b / tau:9.4f} | {amgm:10.1f} | "
            f"{amgm / b:7.3f} | {top:17.1f} | {b / (cut**al * np.log(cut) ** 2):14.4f} | {gap:>10s}   {time.time() - t0:.1f} s"
        )

    print(
        f"\nTHE LAW METER, {name}: R(Q, d) = N_F(Q; d) d_co / A_F(Q) over d <= Q^(1/2), Q = 3^L"
    )
    print(
        "    L |  d <= | max R | argmax | max R at coprime d | argmax | cells within 0.02 of the max | R(Q,4) | 1 + 2^(1-L/2) - 2^(2-L) | N at coprime argmax | Mobius peak there | peak/N^(1/2)"
    )
    for lvl in range(8, 17):
        cut = base**lvl
        members = design_below(base, digits, lvl)
        top = int(np.sqrt(cut))
        while top * top > cut:
            top -= 1
        n = count_multiples(members, top)
        a_f = len(members)
        ds = np.arange(1, top + 1)
        co = np.array([qfree(int(d), base) for d in ds])
        r = n[1:] * co / a_f
        best = float(r.max())
        arg = int(ds[r.argmax()])
        near = ds[r >= best - 0.02].tolist()
        copr = np.gcd(ds, base) == 1
        bestc = float(r[copr].max())
        argc = int(ds[copr][r[copr].argmax()])
        r4 = float(r[3])
        pred = 1 + 2 ** (1 - lvl / 2) - 2 ** (2 - lvl) if lvl % 2 == 0 else float("nan")
        near_s = " ".join(str(v) for v in near[:12]) + (
            " ..." if len(near) > 12 else ""
        )
        cs = np.sort(members[members % argc == 0] // argc)
        mu = mobius_sieve(cut // argc)
        peak = int(np.abs(np.cumsum(mu[cs])).max())
        print(
            f"   {lvl:2d} | {top:5d} | {best:.4f} | {arg:6d} | {bestc:18.4f} | {argc:6d} | {near_s:28s} | {r4:.4f} | {pred:.4f} | "
            f"{cs.size:19d} | {peak:17d} | {peak / cs.size**0.5:12.4f}"
        )

    print(
        f"\nTHE RIPPLE, {name}: A_4(3^L/4) / A_F(3^L/4) against #Acc_4/4 = 3/4, and N_F(3^L; 4) against 2^(L-2) + 2^(L/2-1) - 1"
    )
    print(
        "    L | N_F(3^L;4) = A_4(3^L/4) | A_F(floor(3^L/4)) |  ratio | 2^(L-2)+2^(L/2-1)-1 | A_4(3^L)/2^L"
    )
    for lvl in (8, 10, 12, 14):
        cut = base**lvl
        members = design_below(base, digits, lvl)
        n4 = int((members % 4 == 0).sum())
        af4 = int((members <= cut // 4).sum())
        wide = design_below(base, digits, lvl + 2)
        a4 = (
            int(((wide <= 4 * cut) & (wide % 4 == 0)).sum())
            + 1
            - int(holds(base, digits, 4 * cut))
        )
        print(
            f"   {lvl:2d} | {n4:23d} | {af4:17d} | {n4 / af4:.4f} | {2 ** (lvl - 2) + 2 ** (lvl // 2 - 1) - 1:19d} | {a4 / 2**lvl:.4f}"
        )

    print(
        f"\nTHE REPUNITS, {name}: N_F(3^(kt); R_t) at R_t = (3^t - 1)/2 against 2^t + 1, 2 3^t + 1 and 6^t + 5^t + 3^t + 1"
    )
    print(
        "    t |  R_t | k = 2 | 2^t+1 | k = 3 | 2 3^t+1 | k = 4 | 6^t+5^t+3^t+1 | R at k = 2 | x = 2 3^t + 2 | (U') yard | peak | ratio | peak/N^(1/2)"
    )
    for t in range(2, 9):
        rep = (base**t - 1) // 2
        cells = []
        for k in (2, 3, 4):
            if k * t <= 16:
                cells.append(int((design_below(base, digits, k * t) % rep == 0).sum()))
            else:
                cells.append(None)
        q = base ** (2 * t)
        x = q // rep
        mu = mobius_sieve(x)
        m = design_below(base, digits, 2 * t)
        cs = m[m % rep == 0] // rep
        cs = cs[cs <= x]
        run = np.cumsum(np.where(np.isin(np.arange(x + 1), cs), mu[: x + 1], 0))
        peak = int(np.abs(run).max())
        yard = rep ** ((al - 1) / 2) * x ** (al / 2)
        f = lambda v: f"{v:5d}" if v is not None else "    -"
        print(
            f"   {t:2d} | {rep:4d} | {f(cells[0])} | {2**t + 1:5d} | {f(cells[1])} | {2 * 3**t + 1:7d} | {f(cells[2])} | "
            f"{6**t + 5**t + 3**t + 1:13d} | {cells[0] * rep / 2 ** (2 * t):10.4f} | {x:13d} | {yard:9.3f} | {peak:4d} | {peak / yard:5.3f} | {peak / len(cs) ** 0.5:12.3f}"
        )

    print(
        f"\nTHE ENERGIES, {name}: K_r(t) = #(u_1..u_r, v_1..v_r) in D_t^(2r) with u_1+..+u_r = v_1+..+v_r, carry class of 0, Perron root rho_r, Lambda(2r) = 3 rho_r"
    )
    print(
        "    r | carries | K_r(1) K_r(2) K_r(3) | brute K_r(2) |      rho_r |  Lambda(2r) | Lambda/4^r | certificate min (Mv)_c/v_c > 4^r/3 | char poly"
    )
    for r in range(1, 6):
        cls, mat = pair_carry_matrix(base, digits, r)
        counts = [pair_count(base, digits, r, t) for t in (1, 2, 3)]
        strings = design_below(base, digits, 2)
        strings = np.append(strings[strings < base**2], 0)
        sums = {}
        for tup in np.array(np.meshgrid(*([strings] * r))).reshape(r, -1).T:
            v = int(tup.sum())
            sums[v] = sums.get(v, 0) + 1
        brute = sum(c * c for c in sums.values())
        rho, low, ok = perron_certificate(mat, Fraction(4**r, 3))
        poly = np.rint(np.poly(mat)).astype(np.int64)
        print(
            f"   {r:2d} | {cls!s:>16s} | {counts[0]:6d} {counts[1]:6d} {counts[2]:6d} | {brute:12d} | {rho:10.6f} | {3 * rho:11.6f} | "
            f"{3 * rho / 4**r:10.6f} | {float(low):12.6f} {ok!s:>5s} | {' '.join(str(c) for c in poly)}"
        )

    print(
        f"\nTHE ROOTS OF UNITY, {name}: E_d = sum over a != 0 mod d of abs(hat F_L(a/d)) / 2^L over d <= 3^(L/2)"
    )
    print(
        "    L |  d <= | max E_d | argmax | ord_d(3) | E at d = 3^(L/2) - 1 | (3^t - 1)/2^t - 1 | N_F(3^L; 3^(L/2)-1) | R there"
    )
    for lvl in range(8, 17, 2):
        t0 = time.time()
        cut = base**lvl
        top = base ** (lvl // 2)
        members = design_below(base, digits, lvl)
        e = np.array(
            [
                roots_l1(base, digits, lvl, d) / fill**lvl if d > 1 else 0.0
                for d in range(1, top + 1)
            ]
        )
        arg = int(e.argmax()) + 1
        pinned = top - 1
        t = lvl // 2
        exact = (base**t - 1) / fill**t - 1
        npin = int((members % pinned == 0).sum())
        print(
            f"   {lvl:2d} | {top:5d} | {e.max():7.4f} | {arg:6d} | {mult_order(base, arg):8d} | {e[pinned - 1]:20.4f} | "
            f"{exact:17.4f} | {npin:19d} | {npin * pinned / len(members):.4f}   {time.time() - t0:.1f} s"
        )


# THE REPUNIT DILATE

REPUNIT_DESIGNS = (
    (3, frozenset({0, 1}), "base 3, digits {0,1}", 17, 24, 10),
    (4, frozenset({0, 1}), "base 4, digits {0,1}", 13, 18, 10),
)


def repunit(base, t):
    return (base**t - 1) // (base - 1)


def marked_block(base, digits, t):
    parts = [np.array([base**t], dtype=np.int64)]
    for k in range(t):
        parts.append(
            base**k * (base * (base - 1) * padded_strings(base, digits, t - 1 - k) + 1)
        )
    return np.sort(np.concatenate(parts))


def shifted_sum(base, digits, mu, s):
    return int(mu[base * (base - 1) * padded_strings(base, digits, s) + 1].sum())


def repunit_walk(mu, affine, last):
    run = np.cumsum(mu[affine])
    end = int(run[-1])
    peak = int(np.abs(run).max())
    at = int(affine[int(np.abs(run).argmax())])
    total = end + int(mu[last])
    if abs(total) > peak:
        peak, at = abs(total), last
    return end, total, peak, at


def repunit_pari(base, ts):
    lines = []
    for t in ts:
        lines.append(
            f"tt=getabstime(); n=2^{t}; s=0; mx=0; "
            f"for(i=0,n-1, c={base - 1}*fromdigits(binary(i),{base})+1; s+=moebius(c); mx=max(mx,abs(s))); "
            f'e=s; s+=moebius({base}^{t}+1); mx=max(mx,abs(s)); print({t}," ",e," ",s," ",mx," ",getabstime()-tt);'
        )
    out = subprocess.run(
        ["gp", "-q"],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    rows = {}
    for line in out.split("\n"):
        if line.strip():
            t, e, s, mx, ms = line.split()
            rows[int(t)] = (int(e), int(s), int(mx), int(ms) / 1000)
    return rows


def verb_repunit():
    for base, digits, name, top_sieve, top_pari, top_brute in REPUNIT_DESIGNS:
        t0 = time.time()
        mu = mobius_sieve(base**top_sieve + 1)
        sieve_time = time.time() - t0
        pari = repunit_pari(base, range(2, top_pari + 1))
        alpha = log(len(digits)) / log(base)
        print(
            f"\nTHE REPUNIT DILATE, {name}: R_t = (b^t - 1)/(b - 1), x_t = floor(b^(2t)/R_t) = (b - 1)(b^t + 1), "
            f"R_t^(-1) S_F below x_t = {{(b - 1) m + 1 : m in B_t}} union {{b^t + 1}}, N = 2^t + 1, T_s = sum over B_s of mu(b(b - 1) w + 1), Y_t = R_t^((alpha - 1)/2) x_t^(alpha/2) the (U') yardstick"
        )
        print(
            "    t |        R_t |          x_t |        N | T_(t-1) | T_(t-2) | mu(b^t+1) | M_F(x_t;R_t) | max abs M_F |         at y | N^(1/2) |  N^(0.6) | peak/N^(1/2) | peak/N^(0.6) |   peak/Y_t | checks | pari s"
        )
        for t in range(2, top_pari + 1):
            rep = repunit(base, t)
            n = 2**t + 1
            x = (base - 1) * (base**t + 1)
            assert x == base ** (2 * t) // rep
            ends, totals, peaks = pari[t][:3]
            yard = rep ** ((alpha - 1) / 2) * x ** (alpha / 2)
            if t <= top_sieve:
                b = padded_strings(base, digits, t)
                aff = (base - 1) * b + 1
                assert (np.sort(base**t - (base - 1) * b) == aff).all()
                assert (marked_block(base, digits, t) == aff).all()
                ts = [shifted_sum(base, digits, mu, s) for s in (t - 1, t - 2)]
                end, total, peak, at = repunit_walk(mu, aff, base**t + 1)
                assert end == ts[0] + int(mu[base]) * ts[1]
                assert (end, total, peak) == (ends, totals, peaks)
                checks = "set marked shifted pari"
                if t <= top_brute:
                    m = design_below(base, digits, 2 * t)
                    cs = np.sort(m[m % rep == 0] // rep)
                    assert cs.size == n and (cs == np.append(aff, base**t + 1)).all()
                    checks += " brute"
                tcol = f"{ts[0]:7d} | {ts[1]:7d}"
                ycol = f"{at:12d}"
            else:
                end, total, peak = ends, totals, peaks
                checks = "pari"
                tcol = "      - |       -"
                ycol = "           -"
            print(
                f"   {t:2d} | {rep:10d} | {x:12d} | {n:8d} | {tcol} | {total - end:9d} | {total:12d} | {peak:11d} | {ycol} | "
                f"{n**0.5:7.1f} | {n**0.6:8.1f} | {peak / n**0.5:12.4f} | {peak / n**0.6:12.4f} | {peak / yard:10.3f} | {checks:23s} | {pari[t][3]:.1f}"
            )
        print(
            f"   sieve to {base}^{top_sieve} + 1 in {sieve_time:.1f} s, design total {time.time() - t0:.1f} s"
        )

    base, digits, name = 3, frozenset({0, 2}), "base 3, digits {0,2}"
    print(
        f"\nTHE SCALED DESIGN, {name}: S_F = 2 S_(0,1), so d^(-1) S_F = (d/2)^(-1) S_(0,1) at even d and 2 (d^(-1) S_(0,1)) at odd d"
    )
    print(
        "    t | 3^t-1 | dilate at 3^t - 1 is the {0,1} repunit dilate | M_F(x_t; 3^t - 1) | {0,1} reading | R_t | R_t dilate obeys the law"
    )
    mu = mobius_sieve(3**8 + 2)
    ones = frozenset({0, 1})
    for t in range(2, 7):
        rep = repunit(3, t)
        m = design_below(base, digits, 2 * t)
        m1 = design_below(3, ones, 2 * t)
        aff = np.append(2 * padded_strings(3, ones, t) + 1, 3**t + 1)
        full = np.sort(m[m % (2 * rep) == 0] // (2 * rep))
        half = np.sort(m[m % rep == 0] // rep)
        same = bool(full.size == aff.size and (full == aff).all())
        if rep % 2:
            law = 2 * aff
        else:
            law = np.sort(m1[m1 % (rep // 2) == 0] // (rep // 2))
            law = law[law <= 2 * (3**t + 1)]
        obeys = bool(half.size == law.size and (half == law).all())
        print(
            f"   {t:2d} | {3**t - 1:5d} | {same!s:>43s} | {int(mu[full].sum()):17d} | {int(mu[aff].sum()):13d} | {rep:3d} | {obeys!s:>24s}"
        )


def main():
    verbs = {
        "denominator": verb_denominator,
        "identity": verb_identity,
        "strict": verb_strict,
        "dilate": verb_dilate,
        "converse": verb_converse,
        "sandwich": verb_sandwich,
        "accepting": verb_accepting,
        "repunit": verb_repunit,
    }
    want = sys.argv[1:] or list(verbs)
    t0 = time.time()
    for v in want:
        verbs[v]()
    print(f"\ntotal runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
