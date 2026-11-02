from fractions import Fraction
from itertools import product
from math import pi, prod

import numpy as np
from mpmath import zeta

CORNERS = [tuple(v) for v in product((0, 1), repeat=3)]
Z3 = float(zeta(3))
Z2 = float(zeta(2))


def code_set(mask):
    return frozenset(v for i, v in enumerate(CORNERS) if mask >> i & 1)


def code_mask(P):
    return sum(1 << i for i, v in enumerate(CORNERS) if v in P)


def origin_filled(P):
    return 1 if (0, 0, 0) in P else 0


def prime_factors(n):
    ps, d = [], 2
    while d * d <= n:
        if n % d == 0:
            ps.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        ps.append(n)
    return ps


def radical(n):
    r = 1
    for p in prime_factors(n):
        r *= p
    return r


def mobius(n):
    ps = prime_factors(n)
    r = 1
    for p in ps:
        r *= p
    return 0 if r != n else (-1) ** len(ps)


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def digits_with(q, step, parity):
    return sum(1 for d in range(0, q, step) if d % 2 == parity)


def corner_count(P, q, step=1):
    return sum(prod(digits_with(q, step, c) for c in v) for v in P)


def bracket(P, q):
    k = corner_count(P, q)
    return sum(
        Fraction(mobius(e) * corner_count(P, q, e), k) for e in divisors(radical(q))
    )


def rational_part(P, q):
    d = len(next(iter(P)))
    r = bracket(P, q)
    for p in prime_factors(q):
        r *= Fraction(p**d, p**d - 1)
    return r


def bits(P):
    return [v[0] * 4 + v[1] * 2 + v[2] for v in sorted(P)]


def difference_span(P):
    L = bits(P)
    basis = []
    for w in sorted({a ^ b for a in L for b in L}, reverse=True):
        v = w
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    H = {0}
    for b in basis:
        H |= {h ^ b for h in H}
    return len(basis), H


def classify(P):
    s2, H = difference_span(P)
    if s2 == 3:
        return "span", s2
    if all(x in H for x in bits(P)):
        return "inside", s2
    return "offset", s2


def weight(P):
    return tuple(sum(1 for v in P if sum(v) == j) for j in range(4))


def is_subgroup(P):
    L = set(bits(P))
    return 0 in L and all((a ^ b) in L for a in L for b in L)


def pinned_odd(P):
    return any(all(v[i] == 1 for v in P) for i in range(3))


def even_band(P):
    return Fraction(8, 7) * (1 - Fraction(origin_filled(P), len(P)))


def odd_limits(P, q):
    kind, s2 = classify(P)
    naive = rational_part(P, q)
    if kind == "span":
        return [naive]
    corrected = Fraction(8, 7) * naive * (1 - Fraction(1, 2**s2))
    if kind == "inside":
        return [corrected]
    return [Fraction(8, 7) * naive, corrected]


def odd_limit_over_bases(P):
    kind, s2 = classify(P)
    if kind == "span":
        return Fraction(1)
    if kind == "inside":
        return Fraction(8, 7) * (1 - Fraction(1, 2**s2))
    return None


def even_point_fraction(P, q, n):
    e, o = Fraction(q + 1, 2), Fraction(q - 1, 2)
    w = {v: e ** (3 - sum(v)) * o ** sum(v) for v in P}
    k = sum(w.values())
    total = Fraction(0)
    for t in CORNERS:
        lam = (
            sum(w[v] * (-1) ** (t[0] * v[0] + t[1] * v[1] + t[2] * v[2]) for v in P) / k
        )
        total += lam**n
    return total / 8


def predicted_level(P, q, n):
    if q % 2 == 0:
        return rational_part(P, q)
    return Fraction(8, 7) * rational_part(P, q) * (1 - even_point_fraction(P, q, n))


def design(P, q):
    return np.array(
        [v for v in product(range(q), repeat=3) if (v[0] % 2, v[1] % 2, v[2] % 2) in P],
        dtype=np.int64,
    )


def triangular(rows):
    M = [list(r) for r in rows]
    piv = 0
    for col in range(3):
        r = next((i for i in range(piv, len(M)) if M[i][col]), None)
        if r is None:
            continue
        M[piv], M[r] = M[r], M[piv]
        for i in range(piv + 1, len(M)):
            while M[i][col]:
                f = M[piv][col] // M[i][col]
                M[piv] = [a - f * b for a, b in zip(M[piv], M[i])]
                M[piv], M[i] = M[i], M[piv]
        piv += 1
    return piv, M


def lattice(P, q):
    F = design(P, q).tolist()
    rows = [[a - b for a, b in zip(v, F[0])] for v in F[1:]]
    if not rows:
        return 0, 0
    rank, M = triangular(rows)
    if rank < 3:
        return rank, 0
    return 3, abs(M[0][0] * M[1][1] * M[2][2])


def stack(F, q, n):
    pts = np.zeros((1, 3), dtype=np.int64)
    for _ in range(n):
        pts = (pts[:, None, :] * q + F[None, :, :]).reshape(-1, 3)
    return pts


def visible(pts):
    g = np.gcd(np.gcd(pts[:, 0], pts[:, 1]), pts[:, 2])
    return int(np.count_nonzero(g == 1))


def visible_fraction(P, q, n):
    F = design(P, q)
    if n == 1:
        return visible(F) / len(F)
    base = stack(F, q, n - 1)
    hits = 0
    for v in F:
        hits += visible(base * q + v)
    return hits / len(F) ** n


CARPET = frozenset([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
NET = frozenset([(1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)])
TREE = frozenset([(0, 0, 0), (0, 0, 1)])
VOID = frozenset([(0, 0, 0), (1, 1, 1)])
AXES = frozenset([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
TOP = frozenset([(1, 1, 1)])
PINNED = frozenset([(1, 1, 0), (1, 0, 0)])
TWIN_SUB = frozenset([(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)])
TWIN_SPAN = frozenset([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 1, 1)])
FLAT_CARPET = frozenset([(0, 0), (1, 0), (0, 1)])

FAMILIES = [("carpet", CARPET), ("net", NET), ("tree", TREE), ("void", VOID)]

CASES = [
    ("carpet", CARPET, 2, 12),
    ("carpet", CARPET, 4, 5),
    ("net", NET, 4, 5),
    ("tree", TREE, 3, 7),
    ("tree", TREE, 5, 5),
    ("void", VOID, 3, 8),
    ("void", VOID, 5, 5),
    ("twin subgroup", TWIN_SUB, 3, 6),
    ("twin spanning", TWIN_SPAN, 3, 6),
    ("twin subgroup", TWIN_SUB, 4, 4),
    ("twin spanning", TWIN_SPAN, 4, 4),
]


def band_report(codes):
    bad = [
        (code_mask(P), q)
        for P in codes
        for q in range(2, 41, 2)
        if rational_part(P, q) != even_band(P)
    ]
    print("even bases: 255 codes x even q<=40, exceptions to (8/7)(1-t/|P|):", len(bad))
    band = sorted({even_band(P) for P in codes})
    print(
        "even band values of delta*zeta(3):", len(band), " ".join(str(x) for x in band)
    )
    flat = {rational_part(FLAT_CARPET, q) for q in range(2, 41, 2)}
    v = float(min(flat)) / Z2
    print(
        "D=2 parity carpet at even q<=40: delta*zeta(2) =",
        " ".join(str(x) for x in flat),
        "delta =",
        repr(v),
        "16/(3 pi^2) =",
        repr(16 / (3 * pi**2)),
    )


def similarity_report(codes):
    bad = 0
    for P in codes:
        for q in range(3, 76, 2):
            for e in divisors(radical(q)):
                if corner_count(P, q, e) != corner_count(P, q // e):
                    bad += 1
    print("odd bases: 255 codes x odd q<=75, exceptions to k_e(q)=k(q/e):", bad)


def class_report(codes):
    kinds = {}
    for P in codes:
        kind, s2 = classify(P)
        kinds.setdefault(kind, []).append(s2)
    print(
        "class counts: span",
        len(kinds["span"]),
        "inside",
        len(kinds["inside"]),
        "offset",
        len(kinds["offset"]),
        "total",
        sum(len(v) for v in kinds.values()),
    )
    print(
        "inside by s2:",
        {s: kinds["inside"].count(s) for s in sorted(set(kinds["inside"]))},
    )
    classes = {}
    for P in codes:
        classes.setdefault(weight(P), set()).add(classify(P)[0])
    mixed = [k for k in classes.values() if len(k) > 1]
    print(
        "weight classes:",
        len(classes),
        "mixing regimes:",
        len(mixed),
        "spanning against not:",
        sum(1 for k in mixed if "span" in k),
        "corrected against no-limit:",
        sum(1 for k in mixed if k == {"inside", "offset"}),
    )
    subgroups = sorted(code_mask(P) for P in codes if is_subgroup(P))
    print("subgroup codes:", len(subgroups), " ".join(str(m) for m in subgroups))
    stable = sorted(
        code_mask(P) for P in codes if odd_limit_over_bases(P) == even_band(P)
    )
    print(
        "parity-stable codes:",
        len(stable),
        "equal to subgroup codes:",
        stable == subgroups,
    )
    print("codes with two subsequential limits:", 255 - len(stable))


def lattice_report(codes):
    full = [P for P in codes if lattice(P, 2)[0] == 3]
    print(
        "q=2: full-rank codes:",
        len(full),
        "smallest |P|:",
        min(len(P) for P in full),
        "indices:",
        sorted({lattice(P, 2)[1] for P in full}),
    )
    print(
        "even q=4,6,8,10: lattice indices over all codes:",
        sorted({lattice(P, q)[1] for P in codes for q in (4, 6, 8, 10)}),
    )
    print(
        "q=3: codes of lattice rank <=1:",
        sum(1 for P in codes if lattice(P, 3)[0] <= 1),
        "codes with a coordinate pinned odd:",
        sum(1 for P in codes if pinned_odd(P)),
    )


def family_report():
    print("family table, delta*zeta(3):", " ".join(name for name, _ in FAMILIES))
    print("  every even:", " ".join(str(even_band(P)) for _, P in FAMILIES))
    for q in (3, 5, 7, 9, 11):
        print(
            "  q=%-2d       :" % q,
            " ".join(str(odd_limits(P, q)[-1]) for _, P in FAMILIES),
        )
    print("  odd limit  :", " ".join(str(odd_limit_over_bases(P)) for _, P in FAMILIES))
    print(
        "carpet rational parts at odd q=3,5,7,9,11:",
        " ".join("%.5f" % float(odd_limits(CARPET, q)[-1]) for q in (3, 5, 7, 9, 11)),
        "k(q):",
        " ".join(str(corner_count(CARPET, q)) for q in (3, 5, 7, 9, 11)),
        "dips 1/k(7) = 1/%d, k(3)/k(9) = %d/%d"
        % (corner_count(CARPET, 7), corner_count(CARPET, 3), corner_count(CARPET, 9)),
    )
    print(
        "net: odd q<=81 with rational part <= 1:",
        [q for q in range(3, 82, 2) if rational_part(NET, q) <= 1],
        "B(9) =",
        bracket(NET, 9),
        "= 1 - %d/%d" % (corner_count(NET, 3), corner_count(NET, 9)),
        "k(1) =",
        corner_count(NET, 1),
    )


def exact_report():
    codes = [code_set(m) for m in range(1, 256)]
    band_report(codes)
    similarity_report(codes)
    class_report(codes)
    lattice_report(codes)
    family_report()


def measured_report():
    for name, P, q, n in CASES:
        k = len(design(P, q))
        got = visible_fraction(P, q, n)
        level = float(predicted_level(P, q, n)) / Z3
        lims = (
            [float(x) / Z3 for x in odd_limits(P, q)]
            if q % 2
            else [float(rational_part(P, q)) / Z3]
        )
        print(
            "%-14s q=%d n=%2d k=%3d measured=%.6f level=%.6f limits=%s"
            % (name, q, n, k, got, level, " ".join("%.6f" % x for x in lims))
        )
    print(
        "pinned {110,100} q=3 odd levels:",
        " ".join("%.6f" % visible_fraction(PINNED, 3, n) for n in (3, 5, 7)),
        "outside the trichotomy, whose value is %.6f"
        % (float(odd_limits(PINNED, 3)[0]) / Z3),
    )


def oscillation_report(name, P, q, levels):
    lims = [float(x) / Z3 for x in odd_limits(P, q)]
    print("%s q=%d limits: odd %.6f even %.6f" % (name, q, lims[0], lims[1]))
    for n in range(1, levels + 1):
        print(
            "  n=%d measured=%.6f level=%.6f"
            % (n, visible_fraction(P, q, n), float(predicted_level(P, q, n)) / Z3)
        )


def main():
    exact_report()
    print()
    measured_report()
    print()
    oscillation_report("top corner", TOP, 5, 8)
    print()
    oscillation_report("axes", AXES, 3, 7)


if __name__ == "__main__":
    main()
