import math
import sys
import time
from fractions import Fraction
from functools import cache

import numpy as np

# CORNERS


def corners(code, dim):
    return [i for i in range(2**dim) if code >> i & 1]


def weight(i):
    return i.bit_count()


def profile(code, dim):
    a = [0] * (dim + 1)
    for i in corners(code, dim):
        a[weight(i)] += 1
    return a


def code_of(cs):
    return sum(1 << c for c in cs)


def mirror(code, dim):
    full = 2**dim - 1
    return code_of(full ^ c for c in corners(code, dim))


def complement(code, dim):
    return (1 << 2**dim) - 1 - code


def flip(code, dim, c):
    return code_of(i ^ c for i in corners(code, dim))


def bits(dim):
    codes = np.arange(2 ** (2**dim), dtype=np.int64)
    return (codes[:, None] >> np.arange(2**dim)[None, :]) & 1


# WALSH


def hadamard(dim):
    s = np.arange(2**dim)
    pc = np.array([weight(int(v)) for v in (s[:, None] & s[None, :]).ravel()]).reshape(
        2**dim, 2**dim
    )
    return (-1) ** pc


def spectra(dim):
    return bits(dim) @ hadamard(dim)


def levels(h, dim):
    lv = np.array([weight(s) for s in range(2**dim)])
    return np.stack([h[..., lv == j].sum(axis=-1) for j in range(dim + 1)], axis=-1)


def square_levels(h, dim):
    return levels(h * h, dim)


def evaluate(lev, dim, n):
    return sum(int(lev[j]) * n ** (dim - j) for j in range(dim + 1))


# RENDER


@cache
def parity_index(n, dim):
    grid = np.indices((n,) * dim).reshape(dim, -1) % 2
    return sum(grid[k] << (dim - 1 - k) for k in range(dim))


def histogram(n, dim):
    return np.bincount(parity_index(n, dim), minlength=2**dim).astype(np.int64)


def literal(code, n, dim):
    table = np.array([code >> i & 1 for i in range(2**dim)], dtype=bool)
    return int(np.count_nonzero(table[parity_index(n, dim)]))


def section_expansion():
    bad = [0, 0]
    checks = [0, 0]
    for dim in (1, 2, 3):
        lev = levels(spectra(dim), dim)
        for code in range(1, 2 ** (2**dim)):
            w = code.bit_count()
            for n in range(1, 10):
                f = literal(code, n, dim)
                want = evaluate(lev[code], dim, n) if n % 2 else w * n**dim
                checks[n % 2] += 1
                bad[n % 2] += 2**dim * f != want
    rng = np.random.default_rng(7919)
    sample = rng.choice(np.arange(1, 2**16), size=2000, replace=False)
    lev4 = levels(spectra(4), 4)
    for code in sample:
        code = int(code)
        w = code.bit_count()
        for n in range(1, 10):
            f = literal(code, n, 4)
            want = evaluate(lev4[code], 4, n) if n % 2 else w * n**4
            checks[n % 2] += 1
            bad[n % 2] += 16 * f != want
    b4 = bits(4)
    for n in range(1, 10):
        f = b4 @ histogram(n, 4)
        if n % 2:
            want = sum(lev4[:, j] * n ** (4 - j) for j in range(5))
        else:
            want = b4.sum(axis=1) * n**4
        checks[n % 2] += len(f)
        bad[n % 2] += int(np.count_nonzero(16 * f != want))
    print(
        f"expansion: 2^dim fill(N) = sum_j 2^dim W_j N^(dim-j) at odd N 1..9, {checks[1]} checks, {bad[1]} mismatches; w N^dim at even N 2..8, {checks[0]} checks, {bad[0]} mismatches"
    )
    print(
        "  literal grids: all 273 nonempty codes at dim 1..3 and 2000 seeded codes at dim 4, sides 1..9; literal corner histogram x all 65536 codes at dim 4"
    )
    for dim, code in ((1, 1), (2, 7), (2, 11), (3, 23), (3, 232)):
        lev = levels(spectra(dim), dim)[code]
        ws = ", ".join(str(Fraction(int(x), 2**dim)) for x in lev)
        print(
            f"  dim {dim} code {code}: W = ({ws}), 2^dim fill = {poly_text(lev, dim)}, fill(3) = {literal(code, 3, dim)}"
        )


def poly_text(lev, dim):
    terms = []
    for j in range(dim + 1):
        c = int(lev[j])
        if c == 0:
            continue
        p = dim - j
        mono = "" if p == 0 else ("N" if p == 1 else f"N^{p}")
        coef = str(abs(c)) if (abs(c) != 1 or p == 0) else ""
        body = coef + ("*" if coef and mono else "") + mono
        terms.append(("- " if c < 0 else "+ ") + body)
    s = " ".join(terms)
    return s[2:] if s.startswith("+ ") else "-" + s[2:]


# PROFILE


def krawtchouk(j, w, dim):
    return sum(
        (-1) ** i * math.comb(w, i) * math.comb(dim - w, j - i) for i in range(j + 1)
    )


def section_profile():
    bad = total = 0
    for dim in (1, 2, 3, 4):
        lev = levels(spectra(dim), dim)
        b = bits(dim)
        wt = np.array([weight(i) for i in range(2**dim)])
        a = np.stack([b[:, wt == w].sum(axis=1) for w in range(dim + 1)], axis=1)
        K = np.array(
            [[krawtchouk(j, w, dim) for w in range(dim + 1)] for j in range(dim + 1)]
        )
        bad += int(np.count_nonzero(a @ K.T != lev))
        total += lev.size
        assert np.array_equal(K @ K, 2**dim * np.eye(dim + 1, dtype=int))
    print(
        f"profile: 2^dim W_j = sum_w a_w K_j(w) on every code at dim 1..4, {total} level sums, {bad} mismatches; K^2 = 2^dim I at dim 1..4"
    )


# MIRROR


def section_mirror():
    for dim in (1, 2, 3, 4):
        lev = levels(spectra(dim), dim)
        sign = np.array([(-1) ** (dim - j) for j in range(dim + 1)])
        mbad = 0
        swap, crit, selfdual, other = set(), set(), set(), set()
        for code in range(1, 2 ** (2**dim)):
            neg = sign * lev[code]
            mbad += not np.array_equal(neg, (-1) ** dim * lev[mirror(code, dim)])
            void = -lev[code].copy()
            void[0] += 2**dim
            if np.array_equal(neg, (-1) ** dim * void):
                swap.add(code)
            if np.array_equal(neg, -((-1) ** dim) * void):
                other.add(code)
            a = profile(code, dim)
            if all(a[j] + a[dim - j] == math.comb(dim, j) for j in range(dim + 1)):
                crit.add(code)
            if complement(code, dim) == mirror(code, dim):
                selfdual.add(code)
        formula = 1
        for j in range(dim + 1):
            m = math.comb(dim, j)
            if 2 * j < dim:
                formula *= math.comb(2 * m, m)
            elif 2 * j == dim:
                formula *= math.comb(m, m // 2) if m % 2 == 0 else 0
        nonself = sorted(swap - selfdual)
        print(
            f"mirror dim {dim}: fill_F(-N) = (-1)^dim fill_F'(N) fails on {mbad} codes; swap {len(swap)}, criterion {len(crit)}, equal {swap == crit}, formula {formula}, self-dual {len(selfdual)} all inside {selfdual <= swap}, other sign {len(other)}"
        )
        if dim <= 3:
            print(f"  swapping, not self-dual: {len(nonself)}, least {nonself[:6]}")
        if dim == 2:
            print(
                f"  dim 2 swapping codes {sorted(swap)}; code 7 swaps: {7 in swap}; code 3 swaps: {3 in swap}"
            )
        if dim == 3:
            print(
                f"  dim 3 code 1 swaps: {1 in swap}; code 23 swaps: {23 in swap}; code 27 swaps: {27 in swap}, self-dual: {27 in selfdual}, corners {[format(c, '03b') for c in corners(27, 3)]}, profile {profile(27, 3)}"
            )


# ROOTS


def k_coefficients(code, dim):
    a = profile(code, dim)
    c = [0] * (dim + 1)
    for j, aj in enumerate(a):
        poly = np.poly1d([1.0])
        for _ in range(dim - j):
            poly = poly * np.poly1d([1.0, 0.0])
        for _ in range(j):
            poly = poly * np.poly1d([1.0, -1.0])
        c = np.polyadd(c, aj * poly.coeffs)
    return np.poly1d(c)


def section_roots():
    worst = 0.0
    seen = set()
    count = 0
    for dim in (1, 2, 3):
        lev = levels(spectra(dim), dim)
        for code in range(1, 2 ** (2**dim)):
            key = (dim, tuple(profile(code, dim)))
            if key in seen:
                continue
            seen.add(key)
            count += 1
            P = k_coefficients(code, dim)
            R = np.poly1d([float(lev[code][j]) for j in range(dim, -1, -1)])
            scale = float(np.abs(lev[code]).sum())
            for r in P.roots:
                if abs(r - 0.5) < 1e-9:
                    continue
                mu = 1 / (2 * r - 1)
                worst = max(worst, abs(R(mu)) / (scale * max(1.0, abs(mu)) ** dim))
            grid = np.linspace(-0.999, 0.999, 2001)
            assert np.all(R(grid) > 0)
    lev = levels(spectra(2), 2)[7]
    R = np.poly1d([float(lev[j]) for j in range(2, -1, -1)])
    print(
        f"roots: {count} profiles at dim 1..3, every root r != 1/2 of P_F(n) sends mu = 1/(2r-1) to a zero of R(mu) = sum W_j mu^j, worst scaled residual {worst:.1e}; R > 0 on (-1, 1) at every profile"
    )
    print(
        f"  dim 2 code 7: R roots {sorted(np.round(R.roots, 12))}, P_F roots {sorted(np.round(k_coefficients(7, 2).roots, 12))}"
    )


# DRIFT


def downset(code, dim):
    cs = set(corners(code, dim))
    return all((c & ~(1 << b)) in cs for c in cs for b in range(dim) if c >> b & 1)


def bichromatic(code, dim):
    cs = set(corners(code, dim))
    return sum(
        1
        for c in range(2**dim)
        for b in range(dim)
        if not c >> b & 1 and ((c in cs) != ((c | 1 << b) in cs))
    )


def section_drift():
    bad = total = 0
    dbad = off = above = 0
    dcount = []
    for dim in (1, 2, 3, 4):
        lev = levels(spectra(dim), dim)
        dcount.append(0)
        for code in range(1, 2 ** (2**dim)):
            a = profile(code, dim)
            w = sum(a)
            mean = Fraction(sum(j * aj for j, aj in enumerate(a)), w)
            drift = Fraction(dim, 2) - mean
            total += 1
            bad += drift != Fraction(int(lev[code][1]), 2 * int(lev[code][0]))
            u = bichromatic(code, dim)
            above += int(lev[code][1]) > u
            if downset(code, dim):
                dcount[-1] += 1
                dbad += int(lev[code][1]) != u
            else:
                off += int(lev[code][1]) != u
    print(
        f"drift: dim/2 - mean = W_1/(2 W_0) on all {total} nonempty codes at dim 1..4, {bad} mismatches"
    )
    print(
        f"  I = 2^-dim sum_x s(f, x) = U/2^(dim-1): 2^dim W_1 > U on {above} codes; 2^dim W_1 = U on the nonempty down-sets {dcount} at dim 1..4, {dbad} mismatches; off down-sets it fails on {off} codes"
    )


# STABILITY


def section_stability():
    bad = checks = 0
    for dim in (1, 2, 3):
        h = spectra(dim)
        sq = square_levels(h, dim)
        for n in range(1, 10, 2):
            hist = histogram(n, dim)
            for code in range(1, 2 ** (2**dim)):
                fills = [
                    int(sum(hist[i] for i in corners(flip(code, dim, c), dim)))
                    for c in range(2**dim)
                ]
                inside = sum(fills[c] for c in corners(code, dim))
                checks += 2
                bad += 2**dim * inside != evaluate(sq[code], dim, n)
                bad += sum(fills) != code.bit_count() * n**dim
    polys = sorted(
        {poly_text(levels(spectra(3), 3)[flip(23, 3, c)], 3) for c in range(8)}
    )
    print(
        f"stability: 2^dim sum_(c in F) fill_(F+c)(N) = sum_j 2^(2 dim) W^j N^(dim-j) and sum over all c = w N^dim, every code at dim 1..3, odd sides 1..9; {checks} checks, {bad} mismatches"
    )
    print(
        f"  dim 3 code 23 flip orbit, 8 fills of {len(polys)} polynomials in N (times 8): {polys}"
    )


# ENTROPY


def section_entropy():
    for dim, code in ((1, 1), (1, 2), (2, 7), (3, 23)):
        lev = levels(spectra(dim), dim)[code]
        w = code.bit_count()
        inf = math.log2(2**dim / w)
        row = []
        for n in (2, 3, 9, 27, 81, 243, 729):
            ratio = Fraction(
                evaluate(lev, dim, n) if n % 2 else w * n**dim, 2**dim * n**dim
            )
            row.append(f"{n}: {-math.log2(ratio):.9f}")
        slope = (int(lev[1]) / int(lev[0])) / math.log(2)
        n = 729
        ratio = Fraction(evaluate(lev, dim, n), 2**dim * n**dim)
        print(
            f"entropy dim {dim} code {code}: bits per level {', '.join(row)}, infinity {inf:.9f}; N (limit - cost) at 729 {n * (inf + math.log2(ratio)):.6f} against (W_1/W_0)/ln 2 = {slope:.6f}"
        )


# THRESHOLD


def partial(n, t):
    return sum(math.comb(n, j) for j in range(t + 1))


def level_code(dim, allowed):
    return code_of(i for i in range(2**dim) if weight(i) in allowed)


def section_threshold(top=2048):
    bad = 0
    for dim in range(1, 9):
        H = hadamard(dim)
        for t in range(dim + 1):
            code = level_code(dim, set(range(t + 1)))
            vec = np.array([code >> i & 1 for i in range(2**dim)], dtype=np.int64)
            l = levels(vec @ H, dim)
            bad += int(l[0]) != partial(dim, t)
            bad += int(l[1]) != (t + 1) * math.comb(dim, t + 1)
    print(
        f"threshold: at most t odd has 2^dim W_0 = S(dim, t), 2^dim W_1 = (t+1) C(dim, t+1) at dim 1..8, every t, {bad} mismatches"
    )
    maj = []
    for dim in range(1, 9):
        r = Fraction(partial(dim, dim // 2), 2**dim)
        want = (
            Fraction(1, 2)
            if dim % 2
            else Fraction(1, 2) + Fraction(math.comb(dim, dim // 2), 2 ** (dim + 1))
        )
        maj.append(f"{r}{'' if r == want else ' MISMATCH'}")
    print(f"  majority, at most dim/2 odd, ratio at dim 1..8: {', '.join(maj)}")
    lev4 = levels(spectra(4), 4)
    one = level_code(4, {0, 1})
    two = level_code(4, {2, 3, 4})
    print(
        f"  dim 4: at most one odd code {one} ratio {Fraction(int(lev4[one][0]), 16)}, 2^4 fill = {poly_text(lev4[one], 4)}; at least two odd code {two} ratio {Fraction(int(lev4[two][0]), 16)}, is the mirror of at most two odd: {two == mirror(level_code(4, {0, 1, 2}), 4)}"
    )
    for dim, t in ((2, 1), (3, 1), (4, 2)):
        code = level_code(dim, set(range(t + 1)))
        print(
            f"  dim {dim} at most {t} odd, code {code}: 2^dim fill = {poly_text(levels(spectra(dim), dim)[code], dim)}"
        )
    groups = {}
    row = [1]
    for n in range(1, top + 1):
        row = [1] + [row[i] + row[i + 1] for i in range(len(row) - 1)] + [1]
        s = 0
        for t in range(n):
            s += row[t]
            v = (s & -s).bit_length() - 1
            groups.setdefault((s >> v, n - v), []).append((n, t))
    assert partial(274, 52) == 8 * partial(271, 51)
    half = groups.pop((1, 1))
    odd_majority = half == [(n, (n - 1) // 2) for n in range(1, top + 1, 2)]
    shared = {k: m for k, m in groups.items() if len(m) > 1}
    sizes = sorted({len(m) for m in shared.values()})
    print(
        f"  coincidences of W_0 = S(dim, t)/2^dim among 0 <= t < dim <= {top}: 1/2 is held by exactly the odd majorities: {odd_majority}; {len(shared)} further shared values, group sizes {sizes}:"
    )
    for k, m in sorted(shared.items(), key=lambda kv: kv[1][0]):
        odd, e = k
        if odd == 1:
            label = f"2^-{e}"
        elif odd == (1 << e) - 1:
            label = f"1 - 2^-{e}"
        else:
            label = f"odd/2^{e}, odd of {odd.bit_length()} bits"
        print(f"    {label}: {m}")


# FLAT


def section_flat():
    counts = []
    for dim in (1, 2, 3, 4):
        lev = levels(spectra(dim), dim)
        free = (
            set(np.nonzero(np.all(lev[:, 1:dim] == 0, axis=1))[0].tolist())
            if dim > 1
            else set(range(2 ** (2**dim)))
        )
        crit = set()
        formula_bad = 0
        for code in range(2 ** (2**dim)):
            a = profile(code, dim)
            ev = {Fraction(a[w], math.comb(dim, w)) for w in range(0, dim + 1, 2)}
            od = {Fraction(a[w], math.comb(dim, w)) for w in range(1, dim + 1, 2)}
            if len(ev) <= 1 and len(od) <= 1:
                crit.add(code)
                alpha = ev.pop() if ev else Fraction(0)
                beta = od.pop() if od else Fraction(0)
                l = lev[code]
                formula_bad += Fraction(int(l[0]), 2**dim) != (alpha + beta) / 2
                formula_bad += (
                    Fraction(int(l[dim]), 2**dim) != (alpha - beta) / 2
                    if dim > 0
                    else 0
                )
        counts.append(len(free))
        print(
            f"flat dim {dim}: W_j = 0 for 0 < j < dim on {len(free)} of {2 ** (2**dim)} codes, level-parity criterion {len(crit)}, equal {free == crit}, W_0 = (alpha+beta)/2 and W_dim = (alpha-beta)/2 fails {formula_bad}"
        )
        if dim == 2:
            print(
                f"  dim 2 lower-order-free codes {sorted(free)}; code 11 fill: 2^2 fill = {poly_text(lev[11], 2)}"
            )
        if dim == 4:
            h = spectra(4)
            g = -2 * h
            g[:, 0] += 16
            bent = set(np.nonzero(np.all(np.abs(g) == 4, axis=1))[0].tolist())
            x = code_of(
                c
                for c in range(16)
                if ((c >> 3 & 1) & (c >> 2 & 1)) ^ ((c >> 1 & 1) & (c & 1))
            )
            print(
                f"  dim 4: {len(bent)} bent codes, {len(bent & free)} of them lower-order-free; x1 x2 + x3 x4 is code {x}, bent {x in bent}, profile {profile(x, 4)}, 2^4 fill = {poly_text(lev[x], 4)}"
            )


# MAIN

SECTIONS = {
    "expansion": section_expansion,
    "profile": section_profile,
    "mirror": section_mirror,
    "roots": section_roots,
    "drift": section_drift,
    "stability": section_stability,
    "entropy": section_entropy,
    "threshold": section_threshold,
    "flat": section_flat,
}

if __name__ == "__main__":
    verbs = sys.argv[1:] or list(SECTIONS)
    for verb in verbs:
        start = time.perf_counter()
        SECTIONS[verb]()
        print(f"  [{verb} {time.perf_counter() - start:.1f}s]")
