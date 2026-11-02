from fractions import Fraction
from math import gcd

import numpy as np

# GENERAL LAW


def shadow_mask(base, keep, level):
    span = base**level
    idx = np.arange(span)
    good = np.ones(span, dtype=bool)
    rest = idx.copy()
    for _ in range(level):
        good &= np.isin(rest % base, list(keep))
        rest //= base
    return good


def cell_index(total, span, scale):
    step = total // (span * scale)
    return (np.arange(total) // step) % span


def cross_mean(mask_f, mask_g, m, n):
    span_f = len(mask_f)
    span_g = len(mask_g)
    a = span_f * m
    b = span_g * n
    total = a * b // gcd(a, b)
    hit = mask_f[cell_index(total, span_f, m)] & mask_g[cell_index(total, span_g, n)]
    return Fraction(int(hit.sum()), total)


def cross_cov(mask_f, mask_g, m, n):
    mu_f = Fraction(int(mask_f.sum()), len(mask_f))
    mu_g = Fraction(int(mask_g.sum()), len(mask_g))
    return cross_mean(mask_f, mask_g, m, n) - mu_f * mu_g


# PARITY CHECK


def parity_integral(m, n):
    total = m * n // gcd(m, n)
    i = np.arange(total)
    sm = 1 - 2 * (((m * i) // total) % 2)
    sn = 1 - 2 * (((n * i) // total) % 2)
    return Fraction(int((sm * sn).sum()), total)


def parity_mean(m):
    i = np.arange(m)
    return Fraction(int((1 - 2 * (((m * i) // m) % 2)).sum()), m)


def parity_master(m, n):
    g = gcd(m, n)
    if (m // g) % 2 == 1 and (n // g) % 2 == 1:
        return Fraction(g * g, m * n)
    return Fraction(0)


def check_parity():
    bad = [
        (m, n)
        for m in range(1, 41)
        for n in range(m, 41)
        if parity_integral(m, n) != parity_master(m, n)
    ]
    assert bad == [], "master integral mismatches %s" % bad[:4]
    for m, n in [(3, 5), (3, 9), (5, 15)]:
        g = gcd(m, n)
        want = Fraction(g * g, m * n)
        got = parity_integral(m, n)
        assert got == want, "parity (%d,%d): got %s want %s" % (m, n, got, want)
    named = ", ".join(
        "(%d,%d) %s" % (m, n, parity_integral(m, n))
        for m, n in [(3, 5), (3, 9), (5, 15)]
    )
    print("parity master integral, all pairs to 40: 0 mismatches; %s" % named)
    odd = [n for n in range(1, 41) if n % 2 == 1]
    tree = []
    for m in odd:
        for n in odd:
            if n < m:
                continue
            g = gcd(m, n)
            cov = (parity_integral(m, n) - parity_mean(m) * parity_mean(n)) / 4
            assert cov == Fraction(g * g - 1, 4 * m * n), "tree cov (%d,%d)" % (m, n)
            if g == 1:
                tree.append(cov)
    assert set(tree) == {Fraction(0)}, "tree coprime covariance not zero"
    strip = shadow_mask(2, {1}, 1)
    live = 0
    coprime = 0
    for m in range(1, 41):
        for n in range(m, 41):
            g = gcd(m, n)
            got = cross_cov(strip, strip, m, n)
            want = (
                Fraction(g * g, 4 * m * n)
                if (m // g) % 2 and (n // g) % 2
                else Fraction(0)
            )
            assert got == want, "odd-strip cov (%d,%d): got %s want %s" % (
                m,
                n,
                got,
                want,
            )
            if g == 1:
                coprime += 1
                live += got != 0
    print(
        "tree parity layers: covariance (g^2-1)/(4mn), zero at all %d coprime odd pairs to 40"
        % len(tree)
    )
    print(
        "1-periodic odd strip: covariance g^2/(4mn), %s at (3,5), nonzero at %d of the %d coprime pairs to 40"
        % (cross_cov(strip, strip, 3, 5), live, coprime)
    )


# CARPET STACK


def carpet_mask(level):
    return shadow_mask(3, {0, 2}, level)


def chi3(k):
    r = k % 3
    return 0 if r == 0 else (1 if r == 1 else -1)


def carpet_law_level1(m, n):
    g = gcd(m, n)
    a, b = m // g, n // g
    return Fraction(2 * chi3(a) * chi3(b), 9 * a * b)


def carpet_kernel(level, lo=1, hi=41):
    q = 3**level
    mask = carpet_mask(level)
    table = {}
    hits = 0
    reduce_fail = 0
    zero_fail = 0
    live = 0
    for m in range(lo, hi):
        for n in range(m, hi):
            g = gcd(m, n)
            a, b = m // g, n // g
            cov = cross_cov(mask, mask, m, n)
            if cov != cross_cov(mask, mask, a, b):
                reduce_fail += 1
            if (cov == 0) != (a % q == 0 or b % q == 0):
                zero_fail += 1
            if g == 1 and cov != 0:
                live += 1
            key = (a % q, b % q)
            val = cov * a * b
            if key in table:
                hits += 1
                if table[key] != val:
                    raise AssertionError("residue kernel breaks at %s" % (key,))
            table[key] = val
    return table, hits, reduce_fail, zero_fail, live


def carpet_2d(level, m, n):
    mask = carpet_mask(level)
    span = len(mask)
    a = span * m
    b = span * n
    total = a * b // gcd(a, b)
    fm = mask[cell_index(total, span, m)]
    fn = mask[cell_index(total, span, n)]
    joint = int((np.outer(fm, fm) & np.outer(fn, fn)).sum())
    return (
        Fraction(joint, total * total)
        - Fraction(int(fm.sum()) * int(fn.sum()), total * total) ** 2
    )


def check_carpet():
    mask = carpet_mask(1)
    for m in range(1, 41):
        for n in range(1, 41):
            got = cross_cov(mask, mask, m, n)
            want = carpet_law_level1(m, n)
            assert got == want, "carpet L=1 (%d,%d): got %s want %s" % (m, n, got, want)
    print(
        "carpet L=1: covariance = (2/9) chi(m') chi(n')/(m'n') on all 1600 ordered pairs to 40; "
        "(1,2) %s, (2,5) %s, (1,3) %s"
        % (
            cross_cov(mask, mask, 1, 2),
            cross_cov(mask, mask, 2, 5),
            cross_cov(mask, mask, 1, 3),
        )
    )
    tables = {}
    for level in (1, 2, 3):
        table, hits, reduce_fail, zero_fail, live = carpet_kernel(level)
        assert reduce_fail == 0, "reduction fails at L=%d" % level
        assert zero_fail == 0, "zero law fails at L=%d" % level
        tables[level] = table
        print(
            "carpet L=%d: 820 pairs to 40, reduction Cov(m,n)=Cov(m/g,n/g) exact, "
            "kernel mod %d agrees on %d repeated residue keys, zero set is 3^%d dividing m' or n', "
            "%d of the 490 coprime pairs carry nonzero covariance"
            % (level, 3**level, hits, level, live)
        )
    units = [r for r in range(9) if r % 3]
    rows = " ; ".join(
        "%d: %s" % (u, " ".join(str(tables[2][(min(u, v), max(u, v))]) for v in units))
        for u in units
    )
    print("carpet L=2 kernel G_2(a,b) = Cov * a b on units mod 9, rows a = %s" % rows)
    sep = tables[2][(1, 1)] * tables[2][(4, 4)] - tables[2][(1, 4)] ** 2
    print(
        "carpet L=2 kernel is not separable: G(1,1)G(4,4) - G(1,4)^2 = %s, so no theta(a)theta(b) form"
        % sep
    )
    row3 = " ".join(str(tables[3][(1, v)]) for v in range(1, 27) if v % 3)
    print("carpet L=3 kernel row G_3(1,b) over units b mod 27: %s" % row3)
    for level in (1, 2):
        for m, n in [(1, 2), (2, 5), (1, 3), (4, 7)]:
            got = carpet_2d(level, m, n)
            one = cross_mean(carpet_mask(level), carpet_mask(level), m, n)
            mu = Fraction(2, 3) ** (2 * level)
            assert got == (one - mu) * (one + mu), "2d carpet (%d,%d) L=%d" % (
                m,
                n,
                level,
            )
    print(
        "carpet 2D field: cell-counted covariance equals Cov_1D times (mean + (2/3)^(2L)) at "
        "(1,2),(2,5),(1,3),(4,7) and L=1,2, so the 2D and 1D zero sets coincide"
    )


# LEVELS


def level_product_mask(level):
    top = carpet_mask(1)
    span = 3**level
    out = np.ones(span, dtype=bool)
    for i in range(level):
        out &= top[(np.arange(span) // 3 ** (level - 1 - i)) % 3]
    return out


def check_levels():
    for level in (1, 2, 3, 4):
        assert (level_product_mask(level) == carpet_mask(level)).all(), (
            "level split at %d" % level
        )
    nested = all(
        bool(carpet_mask(3)[i]) <= bool(carpet_mask(2)[i % 9]) for i in range(27)
    )
    assert nested, "f_3 not contained in f_2(3x)"
    print(
        "levels: f_L(x) = product of f_1(3^i x) over i < L exact as masks at L = 1,2,3,4; "
        "f_L(nx) <= f_(L-1)(3nx) pointwise, the gap carrying measure (2/3)^(L-1)/3"
    )
    masks = {level: carpet_mask(level) for level in (1, 2, 3)}
    breaches = 0
    unpredicted = 0
    total = 0
    for level in (1, 2, 3):
        for other in (1, 2, 3):
            for m in range(1, 41):
                for n in range(1, 41):
                    g = gcd(m, n)
                    a, b = m // g, n // g
                    cov = cross_cov(masks[level], masks[other], m, n)
                    dead = b % 3**level == 0 or a % 3**other == 0
                    total += 1
                    if dead and cov != 0:
                        breaches += 1
                    if not dead and cov == 0:
                        unpredicted += 1
    assert breaches == 0 and unpredicted == 0, "cross-level zero law breached"
    print(
        "levels: Cov(f_L(mx), f_L'(nx)) = 0 exactly when 3^L divides n/g or 3^L' divides m/g, "
        "on all %d ordered triples (L,L' in 1..3, m,n to 40); the two levels enter asymmetrically"
        % total
    )
    print(
        "levels: cross-level examples Cov(f_1(x), f_2(3x)) = %s and Cov(f_2(x), f_1(3x)) = %s, "
        "so level L' at scale 3n is not redundant against level L at scale n"
        % (cross_cov(masks[1], masks[2], 1, 3), cross_cov(masks[2], masks[1], 1, 3))
    )


# DESIGNS


def design_mask(base, removed):
    return shadow_mask(base, set(range(base)) - {removed}, 1)


def design_law_33(removed_f, removed_g, m, n):
    g = gcd(m, n)
    a, b = m // g, n // g
    if a % 3 == 0 or b % 3 == 0:
        return Fraction(0)
    root = complex(-0.5, 3**0.5 / 2)
    w = root ** ((a * removed_g - b * removed_f) % 3)
    w *= (1 - root ** ((-b) % 3)) * (1 - root ** (a % 3))
    return Fraction(round(2 * w.real), 27 * a * b)


def design_law_23(removed_g, m, n):
    g = gcd(m, n)
    a, b = m // g, n // g
    if b % 2 == 0 or a % 3 == 0:
        return Fraction(0)
    root = complex(-0.5, 3**0.5 / 2)
    u = root ** ((a * removed_g) % 3) * (1 - root ** (a % 3))
    return Fraction(round(2 * u.real), 18 * a * b)


def check_designs():
    strip = shadow_mask(2, {1}, 1)
    base3 = {removed: design_mask(3, removed) for removed in (0, 1, 2)}
    seen = set()
    for rf in (0, 1, 2):
        for rg in (0, 1, 2):
            for m in range(1, 41):
                for n in range(1, 41):
                    got = cross_cov(base3[rf], base3[rg], m, n)
                    want = design_law_33(rf, rg, m, n)
                    assert got == want, "base-3 pair (%d,%d) at (%d,%d)" % (
                        rf,
                        rg,
                        m,
                        n,
                    )
                    g = gcd(m, n)
                    if got != 0:
                        seen.add(got * 27 * (m // g) * (n // g))
    print(
        "designs, base 3 level 1: Cov = c/(27 m'n') with c in %s, zero exactly when 3 divides m'n', "
        "on all 9 ordered design pairs and 1600 scale pairs to 40"
        % sorted(int(v) for v in seen)
    )
    zero = 0
    for rg in (0, 1, 2):
        for m in range(1, 41):
            for n in range(1, 41):
                got = cross_cov(strip, base3[rg], m, n)
                assert got == design_law_23(rg, m, n), (
                    "strip vs base-3 %d at (%d,%d)" % (rg, m, n)
                )
                if rg == 1 and got == 0:
                    zero += 1
    print(
        "designs, base 2 against base 3: Cov(p(mx), h_r(nx)) = c/(18 m'n') with c in -3,0,3, "
        "zero when n/g is even or 3 divides m/g"
    )
    print(
        "designs: the odd strip is orthogonal to the Sierpinski shadow at every scale pair, "
        "%d of %d pairs exactly zero, because the centred shadow is even and the centred strip odd "
        "under x -> -x" % (zero, 1600)
    )


# GRAM


def jordan2(k):
    out = k * k
    d = 2
    r = k
    while d * d <= r:
        if r % d == 0:
            out = out // (d * d) * (d * d - 1)
            while r % d == 0:
                r //= d
        d += 1
    if r > 1:
        out = out // (r * r) * (r * r - 1)
    return out


def exact_det(rows):
    size = len(rows)
    work = [list(row) for row in rows]
    det = Fraction(1)
    for col in range(size):
        pivot = next((r for r in range(col, size) if work[r][col] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            det = -det
        det *= work[col][col]
        inv = Fraction(1) / work[col][col]
        for r in range(col + 1, size):
            factor = work[r][col] * inv
            if factor:
                for c in range(col, size):
                    work[r][c] -= factor * work[col][c]
    return det


def check_gram():
    for cap in range(1, 13):
        odds = list(range(1, 2 * cap + 2, 2))
        rows = [[Fraction(gcd(a, b) ** 2, a * b) for b in odds] for a in odds]
        got = exact_det(rows)
        want = Fraction(1)
        for k in odds:
            want *= Fraction(jordan2(k), k * k)
        assert got == want, "gram det at K=%d: got %s want %s" % (cap, got, want)
        assert got > 0, "gram det not positive at K=%d" % cap
    odds = list(range(1, 26, 2))
    final = Fraction(1)
    for k in odds:
        final *= Fraction(jordan2(k), k * k)
    print(
        "gram: det[gcd(m,n)^2/(mn)] over odd m,n <= 2K+1 equals prod J_2(k)/k^2 over the same odds, "
        "exact at K = 1..12, so Smith's factor-closed determinant applies to the odd index set"
    )
    print(
        "gram: the value is prod over odd k <= 2K+1 of prod over p | k of (1 - p^-2); at K = 12 it is "
        "%s, positive, so the odd parity layers are linearly independent in L^2" % final
    )


def main():
    check_parity()
    check_carpet()
    check_levels()
    check_designs()
    check_gram()


if __name__ == "__main__":
    main()
