import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import groupby, pairwise
from math import comb, isqrt

from sympy import Poly, discriminant, resultant, symbols

T, N = symbols("t n")
WORKERS = min(8, os.cpu_count() or 1)
CHUNKS = 256
DIMS = (2, 3, 4, 5, 6)
VERBS = {}

# BOX


def sig_bounds(dim, cap=None):
    return [
        min(comb(dim, j), cap if cap else comb(dim, j)) + 1 for j in range(1, dim + 1)
    ]


def box_size(bounds):
    out = 1
    for b in bounds:
        out *= b
    return out


def unrank(index, bounds):
    out = [1]
    for b in bounds:
        out.append(index % b)
        index //= b
    return tuple(out)


def multiplicity(dim, sig):
    out = 1
    for j, s in enumerate(sig):
        out *= comb(comb(dim, j), s)
    return out


def weight(sig):
    return list(sig)


# POLYNOMIALS


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            out[i + j] += u * v
    return out


def polypow(a, k):
    out = [1]
    for _ in range(k):
        out = polymul(out, a)
    return out


def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def poly_of(coeffs, var=T):
    return Poly(list(reversed(coeffs)), var, domain="ZZ")


def factors(coeffs):
    body = trim(list(coeffs))
    if body == [0]:
        return 0, []
    lead, parts = poly_of(body).factor_list()
    out = []
    for g, mult in parts:
        out.append((tuple(int(c) for c in reversed(g.all_coeffs())), mult))
    return int(lead), sorted(out)


def star(g):
    d = len(g) - 1
    out = [0] * (d + 1)
    for j, c in enumerate(g):
        if c:
            term = polymul([0] * j + [c], polypow([1, 1], d - j))
            for i, v in enumerate(term):
                out[i] += v
    return out


def fill_coeffs(dim, sig):
    out = [0] * (dim + 1)
    for j, s in enumerate(sig):
        if s:
            term = polymul([0] * j + [s], polypow([1, 1], dim - j))
            for i, v in enumerate(term):
                out[i] += v
    return out


def corners(dim, code):
    return [
        tuple((c >> i) & 1 for i in range(dim))
        for c in range(1 << dim)
        if (code >> c) & 1
    ]


def signature_of(dim, code):
    sig = [0] * (dim + 1)
    for c in corners(dim, code):
        sig[sum(c)] += 1
    return tuple(sig)


def text(coeffs, var="n"):
    body = []
    for i, c in reversed(list(enumerate(coeffs))):
        if not c:
            continue
        head = "" if c == 1 and i else str(c)
        body.append(head + ("" if i == 0 else var if i == 1 else f"{var}^{i}"))
    return " + ".join(body) if body else "0"


# NORM


def norm_report(dim):
    seen = {}
    bad = []
    lifted = 0
    linear = Counter()
    stray = 0
    for code in range(1 << (1 << dim)):
        sig = signature_of(dim, code)
        if sig not in seen:
            body = trim(list(sig))
            fill = fill_coeffs(dim, sig)
            m = len(body) - 1
            lead, parts = factors(sig)
            build = polypow([1, 1], dim - m)
            for g, mult in parts:
                for _ in range(mult):
                    build = polymul(build, star(g))
            build = [lead * c for c in build]
            res = Poly(
                resultant(poly_of(body, T).as_expr(), N - T * (N + 1), T),
                N,
                domain="ZZ",
            )
            res = polymul(
                [int(c) for c in reversed(res.all_coeffs())], polypow([1, 1], dim - m)
            )
            units = Counter(
                g[1] + 1
                for g, mult in parts
                for _ in range(mult)
                if len(g) == 2 and g[0] == 1
            )
            wild = sum(1 for g, mult in parts if len(g) == 2 and g[0] != 1)
            seen[sig] = (build == fill, res == fill, m < dim, units, wild)
        row = seen[sig]
        if not (row[0] and row[1]):
            bad.append((code, sig))
        lifted += row[2]
        if sig[0] == 1:
            linear += row[3]
            stray += row[4]
    return len(seen), bad, lifted, linear, stray


def disc_pairs(dim):
    ok = 0
    off = 0
    for index in range(box_size(sig_bounds(dim))):
        sig = unrank(index, sig_bounds(dim))
        for g, mult in factors(weight(sig))[1]:
            if len(g) >= 3:
                a = int(discriminant(poly_of(list(g), T)))
                b = int(discriminant(poly_of(star(list(g)), N)))
                ok += a == b
                off += a != b
    return ok, off


def verb_norm():
    print("NORM FORM")
    print(
        "  P(n) = (n+1)^(D-deg W) cont(W) prod over irreducible factors g of W of g*(n), g*(n) = (n+1)^deg g g(n/(n+1))"
    )
    print(
        "  cont(W) is the content, which is 1 on the origin-filled box; the constant is never the leading coefficient"
    )
    print("  the resultant form is P(n) = (n+1)^(D-deg W) Res_t(W(t), n - t(n+1))")
    for dim in (2, 3, 4):
        sigs, bad, lifted, linear, stray = norm_report(dim)
        print(
            f"  D = {dim}: {1 << (1 << dim)} designs, {sigs} signatures, factored form exact on {(1 << (1 << dim)) - len(bad)}, resultant form exact on {(1 << (1 << dim)) - len(bad)}"
        )
        print(
            f"    designs needing the (n+1)^(D-deg W) lift, s_D = 0: {lifted} of {1 << (1 << dim)}"
        )
        print(
            f"    origin-filled designs: linear factors are (a n + 1) with a counted {dict(sorted(linear.items()))}"
        )
        print(f"    origin-filled linear factors not of the form (a n + 1): {stray}")
        ok, off = disc_pairs(dim)
        print(
            f"    disc(g) = disc(g*) over the {box_size(sig_bounds(dim))} origin-filled signatures, {ok} factor slots of degree at least 2 with g(1) != 0, mismatches {off}"
        )
    sig = (0, 2, 1)
    _, parts = factors(sig)
    print(
        "  origin empty breaks the law: signature {} has W = {} and fill {} = {}".format(
            sig,
            text(list(sig), "t"),
            text(fill_coeffs(2, sig)),
            " ".join(
                f"({text(star(list(g)))})" for g, mult in parts for _ in range(mult)
            ),
        )
    )


VERBS["norm"] = verb_norm

# SWEEP


def label_of(top, signs):
    if top == 1:
        return "Q"
    if top == 2:
        return "quadratic " + (
            "mixed" if len(signs) > 1 else ("imaginary" if -1 in signs else "real")
        )
    return f"degree {top}"


def quad_field(g):
    return fundamental_part(g[1] * g[1] - 4 * g[2] * g[0])


def sweep_chunk(job):
    dim, cap, lo, hi = job
    bounds = sig_bounds(dim, cap)
    sigs = Counter()
    codes = Counter()
    degrees = Counter()
    keep = defaultdict(set)
    carried = Counter()
    for index in range(lo, hi):
        sig = unrank(index, bounds)
        parts = factors(sig)[1]
        top = 1
        signs = set()
        fields = set()
        for g, mult in parts:
            d = len(g) - 1
            if d < 2:
                continue
            top = max(top, d)
            degrees[d] += mult
            if d == 2:
                signs.add(-1 if g[1] * g[1] - 4 * g[2] * g[0] < 0 else 1)
                fields.add(quad_field(g))
            if 2 <= d <= 5:
                keep[d].add(g)
        tag = label_of(top, signs)
        sigs[tag] += 1
        weight_ = multiplicity(dim, sig)
        codes[tag] += weight_
        pure = len(trim(list(sig))) == 3
        for field in fields:
            carried[(field, "sigs")] += 1
            carried[(field, "designs")] += weight_
            if pure:
                carried[(field, "pure sigs")] += 1
                carried[(field, "pure designs")] += weight_
            if top == 2:
                carried[(field, "class sigs")] += 1
                carried[(field, "class designs")] += weight_
    return sigs, codes, degrees, {d: s for d, s in keep.items()}, carried


SWEEPS = {}


def sweep(dim, cap=None):
    if (dim, cap) in SWEEPS:
        return SWEEPS[(dim, cap)]
    total = box_size(sig_bounds(dim, cap))
    edges = [total * i // CHUNKS for i in range(CHUNKS + 1)]
    jobs = [
        (dim, cap, edges[i], edges[i + 1])
        for i in range(CHUNKS)
        if edges[i] < edges[i + 1]
    ]
    sigs, codes, degrees, carried = Counter(), Counter(), Counter(), Counter()
    keep = defaultdict(set)
    if len(jobs) == 1:
        results = [sweep_chunk(jobs[0])]
    else:
        with ProcessPoolExecutor(max_workers=WORKERS) as pool:
            results = list(pool.map(sweep_chunk, jobs))
    for a, b, c, d, e in results:
        sigs += a
        codes += b
        degrees += c
        carried += e
        for key, value in d.items():
            keep[key] |= value
    SWEEPS[(dim, cap)] = (total, sigs, codes, degrees, keep, carried)
    return SWEEPS[(dim, cap)]


def order(tag):
    return (
        ("Q", "quadratic imaginary", "quadratic real", "quadratic mixed").index(tag)
        if tag.startswith(("Q", "quad"))
        else 4 + int(tag.split()[-1])
    )


def verb_ladder():
    print("THE LADDER BY DEGREE")
    print(
        "  box: s_0 = 1, 0 <= s_j <= C(D, j); a signature carries C(D, j) choose s_j oriented designs"
    )
    for dim in DIMS:
        total, sigs, codes, degrees, _, _ = sweep(dim)
        print(
            f"  D = {dim}: {total} signatures, {sum(codes.values())} oriented designs"
        )
        for tag in sorted(sigs, key=order):
            print(f"    {tag:20s} signatures {sigs[tag]:8d}  designs {codes[tag]}")
        print(
            f"    irreducible factor slots by degree: {dict(sorted(degrees.items()))}"
        )


VERBS["ladder"] = verb_ladder

# TABLES

TABLES = {
    (3, (1, 1)): [
        23,
        31,
        44,
        59,
        76,
        83,
        87,
        104,
        107,
        108,
        116,
        135,
        139,
        140,
        152,
        172,
        175,
        199,
        200,
        204,
        211,
        212,
        216,
        231,
        239,
        243,
        244,
        247,
        255,
        268,
        283,
        300,
        307,
        324,
        327,
        331,
        335,
        339,
        351,
        356,
        364,
        367,
        379,
        411,
        419,
        424,
        431,
        436,
        439,
        440,
        451,
        459,
        460,
        472,
        484,
        491,
        492,
        499,
        503,
        515,
        516,
        519,
        524,
        527,
        543,
        547,
        563,
        567,
        588,
        620,
        628,
        643,
        648,
        652,
        655,
        671,
        675,
        676,
        679,
        680,
        687,
        695,
        696,
        707,
        716,
        728,
        731,
        743,
        744,
        748,
        751,
        755,
        756,
        759,
        771,
        780,
        804,
        808,
        812,
        815,
    ],
    (3, (3, 0)): [
        49,
        81,
        148,
        169,
        229,
        257,
        316,
        321,
        361,
        404,
        469,
        473,
        564,
        568,
        621,
        697,
        733,
        756,
        761,
        785,
        788,
        837,
        892,
        940,
        961,
        985,
        993,
        1016,
        1076,
        1101,
        1129,
        1229,
        1257,
        1300,
        1304,
        1345,
        1369,
        1373,
        1384,
        1396,
        1425,
        1436,
        1489,
        1492,
        1509,
        1524,
        1556,
        1573,
        1593,
        1620,
        1708,
        1765,
        1772,
        1825,
        1849,
        1901,
        1929,
        1937,
        1940,
        1944,
        1957,
        2021,
        2024,
        2057,
        2089,
        2101,
        2177,
        2213,
        2228,
        2233,
        2241,
        2292,
        2296,
        2300,
        2349,
        2429,
        2505,
        2557,
        2589,
        2597,
        2636,
        2673,
        2677,
        2700,
        2708,
        2713,
        2777,
        2804,
        2808,
        2836,
        2857,
        2917,
        2920,
        2941,
        2981,
        2993,
        3021,
        3028,
        3124,
        3132,
    ],
    (4, (0, 2)): [
        117,
        125,
        144,
        189,
        225,
        229,
        256,
        257,
        272,
        320,
        333,
        392,
        400,
        432,
        441,
        512,
        513,
        549,
        576,
        576,
        592,
        605,
        656,
        657,
        697,
        761,
        784,
        788,
        832,
        837,
        873,
        892,
        981,
        985,
        1008,
        1008,
        1016,
        1025,
        1040,
        1040,
        1076,
        1088,
        1088,
        1089,
        1129,
        1161,
        1168,
        1197,
        1197,
        1225,
        1229,
        1257,
        1264,
        1280,
        1372,
        1384,
        1396,
        1413,
        1421,
        1424,
        1436,
        1489,
        1492,
        1509,
        1521,
        1525,
        1552,
        1556,
        1568,
        1593,
        1600,
        1616,
        1629,
        1728,
        1737,
        1765,
        1805,
        1808,
        1809,
        1813,
        1825,
        1856,
        1872,
        1929,
        1936,
        1937,
        1940,
        1953,
        1953,
        2021,
        2048,
        2048,
        2057,
        2061,
        2089,
        2112,
        2112,
        2133,
        2156,
        2156,
    ],
    (4, (2, 1)): [
        275,
        283,
        331,
        400,
        448,
        475,
        491,
        507,
        563,
        643,
        688,
        731,
        751,
        775,
        848,
        976,
        1024,
        1099,
        1107,
        1156,
        1192,
        1255,
        1323,
        1328,
        1371,
        1375,
        1399,
        1423,
        1424,
        1456,
        1472,
        1472,
        1475,
        1588,
        1600,
        1728,
        1732,
        1775,
        1791,
        1792,
        1823,
        1856,
        1879,
        1927,
        1931,
        1963,
        1968,
        1975,
        1984,
        1984,
        2000,
        2048,
        2051,
        2068,
        2092,
        2096,
        2116,
        2151,
        2183,
        2191,
        2219,
        2243,
        2284,
        2312,
        2319,
        2327,
        2375,
        2412,
        2443,
        2475,
        2480,
        2488,
        2563,
        2608,
        2619,
        2687,
        2696,
        2704,
        2736,
        2763,
        2764,
        2767,
        2787,
        2816,
        2824,
        2843,
        2859,
        2911,
        2943,
        3008,
        3052,
        3119,
        3163,
        3175,
        3188,
        3216,
        3223,
        3267,
        3271,
        3275,
    ],
    (4, (4, 0)): [
        725,
        1125,
        1600,
        1957,
        2000,
        2048,
        2225,
        2304,
        2525,
        2624,
        2777,
        3600,
        3981,
        4205,
        4225,
        4352,
        4400,
        4525,
        4752,
        4913,
    ],
    (5, (1, 2)): [
        1609,
        1649,
        1777,
        2209,
        2297,
        2617,
        2665,
        2869,
        3017,
        3089,
        3233,
        3369,
        3857,
        3889,
        4169,
        4261,
        4409,
        4417,
        4429,
        4432,
        4477,
        4549,
        4597,
        4757,
        4817,
        4897,
        5025,
        5164,
        5437,
        5501,
        5584,
        5653,
        5753,
        5864,
        5913,
        6241,
        6449,
        6581,
        6757,
        6793,
        7096,
        7177,
        7265,
        7333,
        7373,
        7376,
        7672,
        7684,
        7717,
        7909,
        8073,
        8105,
        8249,
        8329,
        8357,
        8529,
        8705,
        8752,
        8945,
        8968,
        9065,
        9137,
        9412,
        9437,
        9489,
        9552,
        9584,
        9664,
        9701,
        9808,
        9829,
        10229,
        10277,
        10329,
        10381,
        10449,
        10492,
        10532,
        10589,
        10609,
        10729,
        10825,
        10832,
        10933,
        11317,
        11332,
        11469,
        11693,
        11809,
        11876,
        11993,
        12184,
        12205,
        12349,
        12389,
        12440,
        12481,
        12517,
        12533,
        12752,
    ],
    (5, (3, 1)): [
        4511,
        4903,
        5519,
        5783,
        7031,
        7367,
        7463,
        8519,
        8647,
        9439,
        9759,
        10407,
        11119,
        11243,
        11551,
        12447,
        13219,
        13523,
        13799,
        13883,
    ],
    (5, (5, 0)): [
        14641,
        24217,
        36497,
        38569,
        65657,
        70601,
        81509,
        81589,
        89417,
        101833,
        106069,
        117688,
        122821,
        124817,
        126032,
        135076,
        138136,
        138917,
        144209,
        147109,
    ],
}

QUADRATIC_REACH = 200


def squarefree(m):
    m = abs(m)
    k = 2
    while k * k <= m:
        if m % (k * k) == 0:
            return False
        k += 1
    return True


def fundamental(d):
    if d % 4 == 1:
        return squarefree(d)
    if d % 4 == 0:
        m = d // 4
        return m % 4 in (2, 3) and squarefree(m)
    return False


def sign_of(sig):
    return (-1) ** sig[1]


def table(degree, sig):
    if degree == 2:
        return [d for d in range(2, QUADRATIC_REACH) if fundamental(sign_of(sig) * d)]
    return TABLES.get((degree, sig))


# FIELDS


def square(value):
    if value <= 0:
        return 0
    root = isqrt(value)
    return root if root * root == value else 0


def gp_text(coeffs):
    return "+".join(f"{c}*x^{i}" for i, c in enumerate(coeffs) if c)


def gp_run(lines, per):
    if not lines:
        return []
    handle, name = tempfile.mkstemp(suffix=".gp")
    with os.fdopen(handle, "w") as fh:
        fh.write("\n".join(lines) + "\nquit\n")
    done = subprocess.run(
        ["gp", "-q", "-s", "500000000", name],
        capture_output=True,
        text=True,
        check=False,
    )
    os.unlink(name)
    values = [int(v) for v in done.stdout.split()]
    if len(values) != per * len(lines):
        raise RuntimeError(f"gp returned {len(values)} values for {len(lines)} lines")
    return [values[i : i + per] for i in range(0, len(values), per)]


def guarded(disc, raw):
    if not disc or raw % disc or disc % 4 not in (0, 1):
        return 0
    return square(raw // disc)


def pari_fields(keys):
    lines = [
        f"P={gp_text(list(reversed(g)))};print(nfdisc(P));print(polsturm(P));print(poldisc(P))"
        for g in keys
    ]
    out = {}
    for g, (disc, r1, raw) in zip(keys, gp_run(lines, 3)):
        d = len(g) - 1
        index = guarded(disc, raw)
        out[g] = (
            d,
            (r1, (d - r1) // 2),
            disc if index else None,
            index,
            raw,
            "value" if index else "guard",
        )
    return out


def pari_classes(carriers, need):
    bodies = [gp_text(list(reversed(g))) for g in carriers]
    pairs = [(i, j) for i in range(len(bodies)) for j in range(i)]
    rows = gp_run([f"print(nfisisom({bodies[i]},{bodies[j]})!=0)" for i, j in pairs], 1)
    same = dict(zip(pairs, [row[0] for row in rows]))
    reps = []
    for i in range(len(bodies)):
        if all(not same[(i, j)] for j in reps):
            reps.append(i)
        if len(reps) >= need:
            break
    return len(reps)


def field_data(keys):
    keys = sorted(keys)
    if not shutil.which("gp"):
        raise SystemExit("verb fields needs PARI, and gp is not on PATH")
    out = {}
    step = 20000
    for start in range(0, len(keys), step):
        out.update(pari_fields(keys[start : start + step]))
    return out


def candidates(delta):
    out = set()
    f = 1
    while f * f <= abs(delta):
        if delta % (f * f) == 0 and (delta // (f * f)) % 4 in (0, 1):
            out.add(delta // (f * f))
        f += 1
    return out


CACHE = {}


def collect():
    if CACHE:
        return CACHE["data"], CACHE["reach"]
    per = {}
    for dim in DIMS:
        per[dim] = sweep(dim)[4]
    keys = set()
    for dim in DIMS:
        for d in per[dim]:
            keys |= per[dim][d]
    data = field_data(keys)
    reach = {}
    for dim in DIMS:
        classes = defaultdict(set)
        open_slots = defaultdict(list)
        members = defaultdict(list)
        for d in per[dim]:
            for g in per[dim][d]:
                degree, sig, disc, _, _, _ = data[g]
                if disc is None:
                    open_slots[(degree, sig)].append(g)
                else:
                    classes[(degree, sig)].add(disc)
                    members[(degree, sig)].append(g)
        reach[dim] = (classes, open_slots, members)
    CACHE["data"], CACHE["reach"] = data, reach
    return data, reach


def risky(key, value, open_slots, data):
    return any(
        value in {abs(v) for v in candidates(data[g][4])}
        for g in open_slots.get(key, ())
    )


def run_and_gap(reached, listed, key=None, open_slots=None, data=None):
    run = []
    gap = None
    owed = []
    for value in listed:
        if value in reached:
            run.append(value)
            continue
        if data is not None and risky(key, value, open_slots, data):
            owed.append(value)
        gap = value
        break
    return run, gap, owed


def distinct_fields(carriers, need, cap=24):
    if not carriers:
        return 0, False
    return pari_classes(carriers[:cap], need), len(carriers) > cap


def verb_fields():
    print("FIELD DISCRIMINANTS BY DEGREE AND SIGNATURE")
    print(
        "  every distinct irreducible factor of W of degree 2..5 over the origin-filled box at D = 2..6"
    )
    print(
        "  field discriminant by PARI nfdisc on the reversed monic model x^d g(1/x), signature by polsturm, guarded against poldisc"
    )
    print(
        "  a value is accepted only if it is 0 or 1 mod 4 and divides the polynomial discriminant with a square quotient, the quotient being the square of the index"
    )
    print(
        "  the polynomial discriminant is used only for that guard, for the index and for the candidate bound on an unresolved factor"
    )
    print(
        "  tables of smallest field discriminants: degree 2 generated here, degrees 3, 4, 5 read from the LMFDB"
    )
    data, reach = collect()
    print(f"  {len(data)} distinct irreducible factors of degree 2..5")
    indexes = Counter(v[3] for v in data.values() if v[2] is not None and v[3])
    print(f"  index [O_K : Z[1/theta]]: {dict(sorted(indexes.items())[:8])}")
    modes = Counter(v[5] for v in data.values())
    bad = sum(n for m, n in modes.items() if m != "value")
    print(f"  nfdisc values failing the guard: {bad} of {len(data)}")
    for dim in DIMS:
        classes, open_slots, _ = reach[dim]
        print(f"  D = {dim}, box {sig_bounds(dim)}")
        for key in sorted(set(classes) | set(open_slots)):
            degree, sig = key
            listed = table(degree, sig)
            reached = {abs(v) for v in classes.get(key, ())}
            run, gap, owed = run_and_gap(reached, listed or [], key, open_slots, data)
            print(
                f"    degree {degree} signature {sig}: {len(classes.get(key, ()))} distinct field discriminants, {len(open_slots.get(key, ()))} unresolved"
            )
            if listed:
                print(
                    "      run {}".format(
                        " ".join(str(v) for v in run) if run else "empty"
                    )
                )
                print(
                    "      first gap {}{}".format(
                        gap if gap is not None else "none inside the table",
                        ""
                        if not owed
                        else f" (owed: an unresolved factor could carry {owed}, no probe)",
                    )
                )
        layer_report(dim, classes)


VERBS["fields"] = verb_fields

# LAYER


def fundamental_part(d):
    f = 1
    best = d
    while f * f <= abs(d):
        if d % (f * f) == 0 and fundamental(d // (f * f)):
            best = d // (f * f)
        f += 1
    return best


def pure_layer(dim):
    out = defaultdict(set)
    for b in range(dim + 1):
        for c in range(1, comb(dim, 2) + 1):
            d = b * b - 4 * c
            root = int(abs(d) ** 0.5)
            if d >= 0 and root * root == d:
                continue
            out["imaginary" if d < 0 else "real"].add(fundamental_part(d))
    return out


def layer_report(dim, classes):
    reach = 4 * comb(dim, 2)
    pure = pure_layer(dim)
    told = {
        "imaginary": {d for d in classes.get((2, (0, 1)), set())},
        "real": {d for d in classes.get((2, (2, 0)), set())},
    }
    window = {d for d in range(-reach, 0) if fundamental(d)}
    print(
        f"    reach 4 C(D, 2) = {reach}, fundamental discriminants inside it {len(window)}"
    )
    for kind in ("imaginary", "real"):
        extra = sorted(told[kind] - pure[kind])
        short = sorted(pure[kind] - told[kind])
        print(
            "      {}: pure layer {} fields, census {} fields, census only {}, pure only {}".format(
                kind,
                len(pure[kind]),
                len(told[kind]),
                extra if extra else "none",
                short if short else "none",
            )
        )
    print(
        "      imaginary pure layer equals the window: {}".format(
            pure["imaginary"] == window
        )
    )


# SWAP


def splits(coeffs):
    return all(len(g) <= 2 for g, mult in factors(coeffs)[1])


def swap_chunk(job):
    dim, lo, hi = job
    bounds = sig_bounds(dim)
    codes = Counter()
    sigs = Counter()
    for index in range(lo, hi):
        sig = unrank(index, bounds)
        void = [comb(dim, j) - sig[j] for j in range(1, dim + 1)]
        key = (splits(list(sig)), splits(void))
        codes[key] += multiplicity(dim, sig)
        sigs[key] += 1
    return codes, sigs


def verb_swap():
    print("THE FILL AND VOID TABLE")
    print("  P(n) is the fill, V(n) = (2n+1)^D - P(n) the void, V(0) = 0 so V = n Q(n)")
    print(
        "  P splits means P is a product of linear factors over Q, V splits means Q is"
    )
    print("  the empty void, W = (1+t)^D, is counted as splitting")
    for dim in (3, 4):
        total = box_size(sig_bounds(dim))
        edges = [total * i // CHUNKS for i in range(CHUNKS + 1)]
        jobs = [
            (dim, edges[i], edges[i + 1])
            for i in range(CHUNKS)
            if edges[i] < edges[i + 1]
        ]
        codes, sigs = Counter(), Counter()
        with ProcessPoolExecutor(max_workers=WORKERS) as pool:
            for a, b in pool.map(swap_chunk, jobs):
                codes += a
                sigs += b
        designs = sum(codes.values())
        print(
            f"  D = {dim}: {designs} oriented designs with the origin filled, {total} signatures"
        )
        for key in ((True, True), (True, False), (False, True), (False, False)):
            tag = "P{} V{}".format("+" if key[0] else "-", "+" if key[1] else "-")
            print(
                f"    {tag} designs {codes[key]:6d} share {codes[key] / designs:.4f} signatures {sigs[key]}"
            )
        full = fill_coeffs(dim, [comb(dim, j) for j in range(dim + 1)])
        solid = polypow([1, 2], dim)
        print(f"    fill of the full box equals (2n+1)^D: {full == solid}")


VERBS["swap"] = verb_swap

# HUNTER


def merged(degree):
    rows = []
    limit = None
    sigs = (
        sorted({k[1] for k in TABLES if k[0] == degree})
        if degree > 2
        else [(0, 1), (2, 0)]
    )
    for sig in sigs:
        listed = table(degree, sig)
        seen = Counter()
        for value in listed:
            seen[value] += 1
            rows.append((value, sig, seen[value]))
        limit = listed[-1] if limit is None else min(limit, listed[-1])
    return sorted(rows), limit


def witnessed(key, value, copy, dim, reach, data, fields):
    classes, open_slots, members = reach[dim]
    reached = {abs(v) for v in classes.get(key, ())}
    if value not in reached:
        return False, value if risky(key, value, open_slots, data) else None
    if copy == 1 or not fields:
        return True, None
    carriers = [g for g in members.get(key, ()) if abs(data[g][2]) == value]
    found, _ = distinct_fields(carriers, copy)
    return found >= copy, None


def verb_hunter():
    print("THE HUNTER BOUND")
    print(
        "  B is the largest bound with every field of degree d and abs discriminant at most B reached by the box"
    )
    print(
        "  the merge is signature aware: one discriminant in two signatures is two fields, and a discriminant listed twice in one signature is two fields resolved by field isomorphism"
    )
    print(
        "  the merged table is valid only below the smallest per-signature table maximum, printed per degree"
    )
    print("  the box is s_0 = 1, 0 <= s_j <= C(D, j); its height is max_j C(D, j)")
    data, reach = collect()
    for degree in (2, 3, 4, 5):
        rows, limit = merged(degree)
        print(
            f"  degree {degree}: {len(rows)} table entries, merge valid below {limit}"
        )
        for dim in (3, 4, 5, 6):
            marks = {}
            for fields in (True, False):
                bound = 0
                miss = None
                owed = None
                for value, group in groupby(rows, key=lambda row: row[0]):
                    if value > limit or miss:
                        break
                    for value, sig, copy in group:
                        ok, short = witnessed(
                            (degree, sig), value, copy, dim, reach, data, fields
                        )
                        if not ok:
                            owed = short
                            miss = (value, sig, copy)
                            break
                    if not miss:
                        bound = value
                marks[fields] = (bound, miss, owed)
            field_bound, miss, owed = marks[True]
            disc_bound = marks[False][0]
            tail = (
                "none inside the table"
                if miss is None
                else "{} at signature {}{}".format(
                    miss[0], miss[1], " listed twice" if miss[2] > 1 else ""
                )
            )
            print(
                "    D = {} box {} height {}: B by field {}, by discriminant {}, first miss {}{}".format(
                    dim,
                    sig_bounds(dim),
                    max(comb(dim, j) for j in range(dim + 1)),
                    field_bound or "none",
                    disc_bound or "none",
                    tail,
                    ""
                    if owed is None
                    else f" (owed: an unresolved factor could carry {owed}, no probe)",
                )
            )


VERBS["hunter"] = verb_hunter

# SIGN


def field_table(carried):
    out = defaultdict(dict)
    for (field, what), value in carried.items():
        out[field][what] = value
    return out


def pure_window(dim):
    return 4 * comb(dim, 2), dim * dim - 4


def quaddisc_check(discs):
    if not shutil.which("gp"):
        return None
    discs = sorted(discs)
    rows = gp_run([f"print(quaddisc({d}))" for d in discs], 1)
    return sum(1 for d, row in zip(discs, rows) if row[0] != fundamental_part(d))


def verb_sign():
    print("THE QUADRATIC FIELDS BY SIGN")
    print(
        "  a signature carries the quadratic field K when some irreducible quadratic factor 1 + b t + c t^2 of W has field discriminant d_K"
    )
    print(
        "  d_K is the fundamental part of b^2 - 4c by exact integer arithmetic, checked against PARI quaddisc on every distinct b^2 - 4c"
    )
    print(
        "  pure means deg W = 2, the layer (1, b, c); class means the top degree of W is 2, the quadratic column of the ladder"
    )
    print(
        "  the pure window is abs(d_K) <= 4 C(D, 2) on the imaginary side and d_K <= D^2 - 4 on the real side; outside the pure layer means no pure signature carries the field"
    )
    tables = {}
    discs = set()
    for dim in DIMS:
        _, sigs, _, _, keep, carried = sweep(dim)
        tables[dim] = field_table(carried)
        discs |= {g[1] * g[1] - 4 * g[2] * g[0] for g in keep.get(2, ())}
    mismatch = quaddisc_check(discs)
    print(
        "  {} distinct order discriminants b^2 - 4c over D = 2..6, quaddisc mismatches {}".format(
            len(discs), "not run, gp missing" if mismatch is None else mismatch
        )
    )
    first = {}
    for dim in DIMS:
        for field in tables[dim]:
            first.setdefault(field, dim)
    for dim in DIMS:
        table = tables[dim]
        imaginary = sorted((f for f in table if f < 0), reverse=True)
        real = sorted(f for f in table if f > 0)
        neg, pos = pure_window(dim)
        by_class = Counter()
        for tag in ("quadratic imaginary", "quadratic real", "quadratic mixed"):
            by_class[tag] = sweep(dim)[2].get(tag, 0)
        print(
            "  D = {}: {} imaginary fields, {} real fields; quadratic class designs {} + {} + {} = {}".format(
                dim,
                len(imaginary),
                len(real),
                by_class["quadratic imaginary"],
                by_class["quadratic real"],
                by_class["quadratic mixed"],
                sum(by_class.values()),
            )
        )
        for kind, fields, window in (
            ("imaginary", imaginary, neg),
            ("real", real, pos),
        ):
            new = [f for f in fields if first[f] == dim]
            beyond = [f for f in fields if not table[f].get("pure sigs", 0)]
            print(
                "    {}: window {}, new at this D {}, outside the pure layer {}".format(
                    kind,
                    window,
                    " ".join(str(f) for f in new) if new else "none",
                    " ".join(str(f) for f in beyond) if beyond else "none",
                )
            )
            print(
                "      {:>6s} {:>8s} {:>22s} {:>8s} {:>22s} {:>8s} {:>22s}".format(
                    "d_K",
                    "sigs",
                    "designs",
                    "pure",
                    "pure designs",
                    "class",
                    "class designs",
                )
            )
            for f in fields:
                row = table[f]
                print(
                    "      {:>6d} {:>8d} {:>22d} {:>8d} {:>22d} {:>8d} {:>22d}".format(
                        f,
                        row.get("sigs", 0),
                        row.get("designs", 0),
                        row.get("pure sigs", 0),
                        row.get("pure designs", 0),
                        row.get("class sigs", 0),
                        row.get("class designs", 0),
                    )
                )
    print(
        "  least D per field, imaginary: {}".format(
            " ".join(
                f"{f}:{first[f]}"
                for f in sorted((f for f in first if f < 0), reverse=True)
            )
        )
    )
    print(
        "  least D per field, real: {}".format(
            " ".join(f"{f}:{first[f]}" for f in sorted(f for f in first if f > 0))
        )
    )
    for dim in DIMS:
        table = tables[dim]
        for kind, fields in (
            ("imaginary", sorted((f for f in table if f < 0), reverse=True)),
            ("real", sorted(f for f in table if f > 0)),
        ):
            designs = [table[f]["designs"] for f in fields]
            sigs = [table[f]["sigs"] for f in fields]
            print(
                "  D = {} {}: designs per field fall with abs(d_K) {}, signatures per field fall with abs(d_K) {}".format(
                    dim,
                    kind,
                    "strictly"
                    if all(a > b for a, b in pairwise(designs))
                    else (
                        "weakly" if all(a >= b for a, b in pairwise(designs)) else "no"
                    ),
                    "strictly"
                    if all(a > b for a, b in pairwise(sigs))
                    else ("weakly" if all(a >= b for a, b in pairwise(sigs)) else "no"),
                )
            )
    for dim in DIMS:
        neg, pos = pure_window(dim)
        pure_imag = {
            f for f in tables[dim] if f < 0 and tables[dim][f].get("pure sigs", 0)
        }
        window = {d for d in range(-neg, 0) if fundamental(d)}
        print(
            "  D = {}: the pure imaginary fields equal the window: {}; fields with one signature each: {}".format(
                dim,
                pure_imag == window,
                " ".join(
                    str(f) for f in sorted(tables[dim]) if tables[dim][f]["sigs"] == 1
                )
                or "none",
            )
        )


VERBS["sign"] = verb_sign

# BEYOND


def root_bound(dim):
    r = 2 ** (1 / dim) - 1
    return int(1 / (r * r)), int(2 / r)


NODES = [0]


def root_tails(dim, b, c):
    disc = b * b - 4 * c
    if disc < 0:
        roots = [complex(-b, (-disc) ** 0.5) / (2 * c)]
    else:
        roots = [
            complex((-b + disc**0.5) / (2 * c)),
            complex((-b - disc**0.5) / (2 * c)),
        ]
    powers = [[theta**j for j in range(dim + 1)] for theta in roots]
    tails = [
        [
            sum(comb(dim, i) * abs(theta) ** i for i in range(j + 1, dim + 1))
            * (1 + 1e-9)
            + 1e-9
            for j in range(dim + 1)
        ]
        for theta in roots
    ]
    return powers, tails


def cofactors(dim, b, c, top, prune=True):
    bounds = [comb(dim, j) for j in range(dim + 1)]
    powers, tails = root_tails(dim, b, c)
    found = []

    def rec(h, sums):
        NODES[0] += 1
        k = len(h) - 1
        if k + 2 <= dim and h[k] >= 1:
            s3 = b * h[k] + (c * h[k - 1] if k >= 1 else 0)
            s4 = c * h[k]
            if 0 <= s3 <= bounds[k + 1] and 1 <= s4 <= bounds[k + 2]:
                w = [0] * (k + 3)
                for j in range(k + 3):
                    w[j] = (
                        (h[j] if j <= k else 0)
                        + (b * h[j - 1] if 1 <= j <= k + 1 else 0)
                        + (c * h[j - 2] if j >= 2 else 0)
                    )
                found.append(tuple(w) + (0,) * (dim - k - 2))
        if k + 1 > top or k + 3 > dim:
            return
        j = k + 1
        base = b * h[k] + (c * h[k - 1] if k >= 1 else 0)
        for wj in range(bounds[j] + 1):
            nxt = [sums[r] + wj * powers[r][j] for r in range(len(sums))]
            if prune and any(abs(nxt[r]) > tails[r][j] for r in range(len(sums))):
                continue
            rec(h + [wj - base], nxt)

    rec([1], [complex(1)] * len(powers))
    return found


def divide(w, b, c):
    h = []
    for j in range(len(w) - 2):
        h.append(
            w[j] - (b * h[j - 1] if j >= 1 else 0) - (c * h[j - 2] if j >= 2 else 0)
        )
    return h


def beyond_scan(dim, top, prune=True):
    cmax, bmax = root_bound(dim)
    out = defaultdict(list)
    pairs = 0
    NODES[0] = 0
    for c in range(1, cmax + 1):
        for b in range(-bmax, bmax + 1):
            disc = b * b - 4 * c
            if disc == 0 or (disc > 0 and (b < 1 or square(disc))):
                continue
            pairs += 1
            field = fundamental_part(disc)
            for w in cofactors(dim, b, c, top, prune):
                out[field].append((w, b, c))
    return out, pairs


def pure_dim(field):
    dim = 2
    while field not in pure_layer(dim)["imaginary" if field < 0 else "real"]:
        dim += 1
    return dim


def verb_beyond():
    print("BEYOND THE PURE LAYER")
    print(
        "  root bound: a root theta of a box W at D has abs(theta) >= 2^(1/D) - 1 = R, since 1 = abs(sum_(j>=1) s_j theta^j) <= (1 + abs(theta))^D - 1"
    )
    print(
        "  a quadratic factor 1 + b t + c t^2 has root product 1/c, so c <= 1/R^2; a real one has small root 2/(b + sqrt(d)), so b + sqrt(d) <= 2/R"
    )
    print(
        "  walk: every W = (1 + b t + c t^2) h with h in Z[t] of degree at most H and W in the box, over all (b, c) inside the root bound; exhaustive when H = D - 2, a lower bound on the census otherwise"
    )
    print(
        "  the pure layer is the set of fields carried by (1, b, c) alone, computed from the definition; pure window means its bounding box"
    )
    print(
        "  pruning: with S_j = sum_(i <= j) W_i theta^i at each root theta of the factor, W(theta) = 0 gives abs(S_j) <= sum_(i > j) C(D, i) abs(theta)^i, so a prefix breaking it is dropped"
    )
    print(
        "  the pruning is validated by an unpruned control at D = 6, H = 4 and D = 7, H = 3, which must return the same W"
    )
    print("  at D <= 6 the cut is checked against the full census of the sign verb")
    for dim, top in ((6, 4), (7, 3)):
        started = time.time()
        pruned, _ = beyond_scan(dim, top, True)
        nodes = NODES[0]
        plain, _ = beyond_scan(dim, top, False)
        same = {k: sorted(v) for k, v in pruned.items()} == {
            k: sorted(v) for k, v in plain.items()
        }
        print(
            f"  control D = {dim} H = {top}: pruned {sum(len(v) for v in pruned.values())} W over {nodes} nodes, unpruned {sum(len(v) for v in plain.values())} W over {NODES[0]} nodes, same W {same}; {time.time() - started:.1f} s"
        )
    full = {}
    for dim in DIMS:
        full[dim] = set(field_table(sweep(dim)[5]))
    print(
        "  family (1 - t + c t^2)(1 + t) = 1 + (c - 1) t^2 + c t^3 with c = C(D, 2) + 1, in the box iff C(D, 2) + 1 <= C(D, 3), which is D >= 6: order discriminant -(2 D (D - 1) + 3)"
    )
    print(
        "    "
        + " ".join(
            "D={}:{}{}".format(
                d,
                -(2 * d * (d - 1) + 3),
                "" if fundamental(-(2 * d * (d - 1) + 3)) else "(order)",
            )
            for d in range(6, 21)
        )
    )
    outside = Counter()
    counts = []
    for dim, top in ((3, 1), (4, 2), (5, 3), (6, 4), (7, 5), (8, 4), (9, 2), (10, 2)):
        started = time.time()
        cmax, bmax = root_bound(dim)
        neg, pos = pure_window(dim)
        found, pairs = beyond_scan(dim, top)
        fields = set(found)
        imaginary = sorted((f for f in fields if f < 0), reverse=True)
        real = sorted(f for f in fields if f > 0)
        counts.append((dim, len(imaginary), len(real)))
        print(
            "  D = {}: root bound c <= {}, b + sqrt(d) <= {}; {} pairs (b, c); H = {} ({}); {} W over {} nodes; pure window {} and {}; {:.1f} s".format(
                dim,
                cmax,
                bmax,
                pairs,
                top,
                "exhaustive" if top >= dim - 2 else "cut",
                sum(len(v) for v in found.values()),
                NODES[0],
                neg,
                pos,
                time.time() - started,
            )
        )
        if dim in full:
            missed = sorted(full[dim] - fields)
            print(
                "    cut against the census: {} of {} fields reached, extra {}, missed {}".format(
                    len(fields & full[dim]),
                    len(full[dim]),
                    sorted(fields - full[dim]) or "none",
                    missed or "none",
                )
            )
        pure = pure_layer(dim)
        for kind, listed, window in (
            ("imaginary", imaginary, neg),
            ("real", real, pos),
        ):
            beyond = [f for f in listed if f not in pure[kind]]
            print(
                "    {}: {} fields, pure layer {} fields, outside the pure layer {}".format(
                    kind,
                    len(listed),
                    len(pure[kind]),
                    " ".join(str(f) for f in beyond) if beyond else "none",
                )
            )
            for f in beyond:
                w, b, c = min(found[f])
                print(
                    "      {} carried by signature {} = ({})({})".format(
                        f, w, text([1, b, c], "t"), text(divide(w, b, c), "t")
                    )
                )
            if beyond:
                print(
                    "      pure D of each field outside the pure layer: {}".format(
                        " ".join(f"{f}:{pure_dim(f)}" for f in beyond)
                    )
                )
                for f in beyond:
                    outside[pure_dim(f) - dim] += 1
            run = []
            for d in range(2, 400):
                if fundamental(-d if kind == "imaginary" else d):
                    if (-d if kind == "imaginary" else d) in fields:
                        run.append(d)
                    else:
                        break
            print(
                "      gapless run to {}, first gap {}".format(
                    run[-1] if run else "none", d
                )
            )

    print(
        f"  fields found outside the pure layer over D = 3..10, by pure D minus D: {dict(sorted(outside.items()))}"
    )
    print(
        "  fields per D over D = 3..10, imaginary then real: {}".format(
            " ".join(f"{d}:{a}+{b}" for d, a, b in counts)
        )
    )


VERBS["beyond"] = verb_beyond

# POLYA

MISSED = {
    "2.2.37.1": (37, [-9, -1, 1]),
    "3.3.316.1": (316, [2, -4, -1, 1]),
    "4.4.1125.1": (1125, [1, 4, -4, -1, 1]),
    "5.5.14641.1": (14641, [-1, 3, 3, -4, -1, 1]),
    "2.0.87.1": (-87, [22, 1, 1]),
}


def shifted(m, shift):
    d = len(m) - 1
    out = [0] * (d + 1)
    for j, a in enumerate(m):
        term = [a * ((-1) ** j)]
        term = polymul(term, polypow([shift, 1], j))
        for i, v in enumerate(term):
            out[i] += v
    if out[d] < 0:
        out = [-v for v in out]
    return out


def real_ceiling(m):
    rows = gp_run(
        [
            f"P={gp_text(m)};print(polsturm(P));print(if(polsturm(P),floor(vecmax(polrootsreal(P)))+1,0))"
        ],
        2,
    )
    return rows[0][0], rows[0][1]


def polya_exponent(g):
    w = list(g)
    m = 0
    while min(w) < 0:
        w = polymul(w, [1, 1])
        m += 1
    return m, w


def room(g):
    d = len(g) - 1
    e = 0
    while any(g[k] > comb(d + e, k) for k in range(d + 1)):
        e += 1
    return e


def least_dim(w):
    dim = len(w) - 1
    while any(w[j] > comb(dim, j) for j in range(len(w))):
        dim += 1
    return dim


def powers_reznick(g):
    d = len(g) - 1
    big = max(abs(g[j]) / comb(d, j) for j in range(d + 1))
    script = f"P={gp_text(g)};r=polrootsreal(deriv(P)*(1+x)-{d}*P);m=min(1,pollead(P));for(i=1,#r,if(r[i]>=0,m=min(m,subst(P,x,r[i])/(1+r[i])^{d})));print(floor(m*10^9))"
    lam = gp_run([script], 1)[0][0] / 10**9
    return big, lam, (d * d - d) * big / (2 * lam) - d


def verb_polya():
    print("POLYA FOR THE FIELDS THE BOX MISSES AT D = 6")
    print(
        "  a field K with generator gamma gives theta = -1/(gamma + N), N above every real conjugate, whose primitive minimal polynomial g has g(0) = 1, positive leading coefficient and no root in [0, oo)"
    )
    print(
        "  Polya: (1 + t)^m g(t) has positive coefficients for m large, with m > (d^2 - d) L(g)/(2 lambda(g)) - d, L = max_j abs(g_j)/C(d, j), lambda = inf_(t >= 0) g(t)/(1+t)^d"
    )
    print(
        "  W = (1 + t)^m g then sits in the box at every D with W_j <= C(D, j); the least such D over m is printed, an upper bound on the D where K first appears"
    )
    print(
        "  the generator is the LMFDB polynomial of the field, its discriminant recomputed by nfdisc; lambda is the minimum over t = 0, the critical points of g(t)/(1+t)^d on [0, oo) and the limit lc(g), rounded down at nine decimals; N runs over eight translates and the least D is kept"
    )
    if not shutil.which("gp"):
        raise SystemExit("verb polya needs PARI, and gp is not on PATH")
    for label, (disc, m) in MISSED.items():
        check = gp_run([f"print(nfdisc({gp_text(m)}))"], 1)[0][0]
        r1, least = real_ceiling(m)
        rows = []
        for shift in range(least, least + 8):
            g = list(reversed(shifted(m, shift)))
            assert g[0] == 1 and g[-1] > 0
            exponent, _ = polya_exponent(g)
            big, lam, bound = powers_reznick(g)
            best = min(
                least_dim(polymul(g, polypow([1, 1], exponent + extra)))
                for extra in range(40)
            )
            rows.append((best, shift, g, exponent, bound, big, lam, room(g)))
        best, shift, g, exponent, bound, big, lam, e = min(rows)
        print(
            "  {} disc {} (nfdisc {}): r1 = {}, N from {}; best translate N = {}, g = {}".format(
                label, disc, check, r1, least, shift, text(g, "t")
            )
        )
        print(
            f"    Polya exponent m = {exponent} (least m with (1 + t)^m g nonnegative), Powers-Reznick m > {bound:.2f} (L = {big:.3f}, lambda = {lam:.6f}); room e = {e}; least D over m is {best}"
        )


VERBS["polya"] = verb_polya

# MAIN


def main():
    asked = sys.argv[1:] or ["all"]
    if asked == ["all"]:
        asked = [
            "norm",
            "ladder",
            "sign",
            "fields",
            "swap",
            "hunter",
            "beyond",
            "polya",
        ]
    started = time.time()
    for name in asked:
        VERBS[name]()
    print()
    print(f"wall time {time.time() - started:.1f} s on {WORKERS} workers")


if __name__ == "__main__":
    main()
