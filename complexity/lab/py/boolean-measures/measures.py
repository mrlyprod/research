import csv
import os
import random
import time
from fractions import Fraction
from itertools import permutations
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DIMS = (3, 4)
MEASURES = ("s", "bs", "C", "dt", "deg", "dnf", "cnf")
KEYS = (
    ("genus",),
    ("gf2deg",),
    ("pop",),
    ("gf2deg", "pop"),
    ("fingerprint",),
    ("genus", "fingerprint"),
    ("genus", "gf2deg", "pop", "fingerprint"),
)
LETTERS = "xyzw"
SEED = 20260725
TRIALS = 50000
SAMPLE = 20000
LEVELS = (1, 2, 3, 4, 6, 8, 12, 16)
NAMED = (
    "bang dim 4, code 27",
    "bang dim 4, code 281",
    "bang dim 4, code 855",
    "bang dim 4, code 1911",
    "bang dim 4, code 7128",
)
POP16 = np.array([i.bit_count() for i in range(1 << 16)], dtype=np.uint8)


def corner_maps(d):
    n = 1 << d
    maps = []
    for pi in permutations(range(d)):
        for flip in range(n):
            m = [0] * n
            for c in range(n):
                image = 0
                for a in range(d):
                    b = ((c >> (d - 1 - pi[a])) & 1) ^ ((flip >> (d - 1 - a)) & 1)
                    image |= b << (d - 1 - a)
                m[c] = image
            maps.append([1 << target for target in m])
    return maps


def catalog(d):
    n = 1 << d
    maps = corner_maps(d)
    seen = bytearray(1 << n)
    out = []
    for code in range(1 << n):
        if seen[code]:
            continue
        cells = [i for i in range(n) if (code >> i) & 1]
        orbit = set()
        for m in maps:
            moved = 0
            for i in cells:
                moved |= m[i]
            orbit.add(moved)
        for member in orbit:
            seen[member] = 1
        out.append((code, sorted(orbit)))
    return out


def bit(f, x):
    return (f >> x) & 1


def differing(f, x, n):
    return f ^ ((1 << n) - 1) if bit(f, x) else f


def cubes(d):
    n = 1 << d
    out = {}
    for fixed in range(n):
        for vals in range(n):
            if vals & fixed == vals:
                out[(fixed, vals)] = sum(1 << y for y in range(n) if y & fixed == vals)
    return out


def constant(f, mask):
    return f & mask in (0, mask)


def sensitivity(f, d):
    n = 1 << d
    best = 0
    for x in range(n):
        m = differing(f, x, n)
        best = max(best, sum(bit(m, x ^ (1 << i)) for i in range(d)))
    return best


def pack(blocks):
    memo = {}

    def go(used):
        if used in memo:
            return memo[used]
        top = 0
        for b in blocks:
            if b & used == 0:
                top = max(top, 1 + go(used | b))
        memo[used] = top
        return top

    return go(0)


def block_sensitivity(f, d):
    n = 1 << d
    best = 0
    for x in range(n):
        m = differing(f, x, n)
        best = max(best, pack([x ^ y for y in range(n) if bit(m, y)]))
    return best


def certificate(f, d, cube):
    n = 1 << d
    order = sorted(range(n), key=int.bit_count)
    best = 0
    for x in range(n):
        for fixed in order:
            if constant(f, cube[(fixed, x & fixed)]):
                best = max(best, fixed.bit_count())
                break
    return best


def dt_depth(f, d, cube):
    memo = {}

    def go(fixed, vals):
        key = (fixed, vals)
        if key in memo:
            return memo[key]
        if constant(f, cube[key]):
            out = 0
        else:
            out = d
            for i in range(d):
                if not (fixed >> i) & 1:
                    low = go(fixed | 1 << i, vals)
                    high = go(fixed | 1 << i, vals | 1 << i)
                    out = min(out, 1 + max(low, high))
        memo[key] = out
        return out

    return go(0, 0)


def mobius(f, d, gf2):
    n = 1 << d
    coef = [bit(f, x) for x in range(n)]
    for i in range(d):
        step = 1 << i
        for s in range(n):
            if s & step:
                coef[s] = coef[s] ^ coef[s ^ step] if gf2 else coef[s] - coef[s ^ step]
    return coef


def top_weight(coef):
    return max((s.bit_count() for s, c in enumerate(coef) if c), default=-1)


def anf(f, d):
    words = []
    for s, c in enumerate(mobius(f, d, True)):
        if c:
            word = "".join(LETTERS[a] for a in range(d) if (s >> (d - 1 - a)) & 1)
            words.append((len(word), word or "1"))
    return " + ".join(word for _, word in sorted(words)) or "0"


def min_cover(primes, need):
    best = [len(primes)]

    def go(left, count):
        if left == 0:
            best[0] = min(best[0], count)
            return
        if count + 1 >= best[0]:
            return
        cells = [y for y in range(left.bit_length()) if bit(left, y)]
        y = min(cells, key=lambda c: sum(1 for p in primes if bit(p, c)))
        for p in primes:
            if bit(p, y):
                go(left & ~p, count + 1)

    go(need, 0)
    return best[0]


def cover_size(f, d, cube, target):
    n = 1 << d
    need = f if target else f ^ ((1 << n) - 1)
    if need == 0:
        return 0
    good = {key: mask for key, mask in cube.items() if mask & need == mask}
    primes = []
    for (fixed, vals), mask in good.items():
        freed = (
            (fixed & ~(1 << i), vals & ~(1 << i)) for i in range(d) if (fixed >> i) & 1
        )
        if all(key not in good for key in freed):
            primes.append(mask)
    return min_cover(primes, need)


def all_measures(f, d, cube):
    return {
        "s": sensitivity(f, d),
        "bs": block_sensitivity(f, d),
        "C": certificate(f, d, cube),
        "dt": dt_depth(f, d, cube),
        "deg": top_weight(mobius(f, d, False)),
        "dnf": cover_size(f, d, cube, 1),
        "cnf": cover_size(f, d, cube, 0),
    }


def fingerprint(code, d):
    coef = [0] * (d + 1)
    for i in range(1 << d):
        if bit(code, i):
            ones = i.bit_count()
            for j in range(ones + 1):
                coef[d - ones + j] += comb(ones, j) * (-1) ** (ones - j)
    return tuple(reversed(coef))


def poly_text(coef):
    d = len(coef) - 1
    out = ""
    for j, c in enumerate(coef):
        if c == 0:
            continue
        power = d - j
        size = "" if abs(c) == 1 and power else str(abs(c))
        var = "" if power == 0 else "k" if power == 1 else f"k^{power}"
        sign = "-" if c < 0 else "+"
        out += f" {sign} {size}{var}" if out else f"{'-' if c < 0 else ''}{size}{var}"
    return out or "0"


def render_index(d, side):
    grid = np.indices((side,) * d).reshape(d, -1) & 1
    weights = np.array([1 << (d - 1 - a) for a in range(d)])
    return (grid * weights[:, None]).sum(axis=0)


def fitted(code, d, renders):
    m = d + 1
    table = np.array([bit(code, i) for i in range(1 << d)], dtype=np.int64)
    rows = []
    for k in range(1, m + 1):
        fill = int(table[renders[k]].sum())
        rows.append([Fraction(k ** (d - j)) for j in range(m)] + [Fraction(fill)])
    for col in range(m):
        pivot = next(r for r in range(col, m) if rows[r][col] != 0)
        rows[col], rows[pivot] = rows[pivot], rows[col]
        lead = rows[col][col]
        rows[col] = [v / lead for v in rows[col]]
        for r in range(m):
            if r != col and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[col])]
    return tuple(int(row[m]) for row in rows)


def is_levelset(code, d):
    value = [-1] * (d + 1)
    for i in range(1 << d):
        w = i.bit_count()
        if value[w] not in (-1, bit(code, i)):
            return False
        value[w] = bit(code, i)
    return True


def is_pin(code, d):
    cells = [i for i in range(1 << d) if bit(code, i)]
    if not cells:
        return False
    fixed = sum(1 for i in range(d) if len({bit(c, i) for c in cells}) == 1)
    return len(cells) == 1 << (d - fixed)


def genus(orbit, d):
    if any(is_levelset(c, d) for c in orbit):
        return "iso"
    if any(is_pin(c, d) for c in orbit):
        return "axis"
    return "comp"


def rows_of(d):
    cube = cubes(d)
    out = []
    for code, orbit in catalog(d):
        row = {
            "name": f"bang dim {d}, code {code}",
            "code": code,
            "orbit": len(orbit),
            "genus": genus(orbit, d),
            "gf2deg": top_weight(mobius(code, d, True)),
            "pop": code.bit_count(),
            "fingerprint": fingerprint(code, d),
        }
        row.update(all_measures(code, d, cube))
        out.append(row)
    return out


def csv_cells(row):
    fp = " ".join(str(c) for c in row["fingerprint"])
    head = [
        row["name"],
        str(row["code"]),
        str(row["orbit"]),
        row["genus"],
        str(row["gf2deg"]),
        str(row["pop"]),
        fp,
    ]
    return head + [str(row[m]) for m in MEASURES]


def csv_diff(rows, d):
    with open(os.path.join(HERE, f"measures_d{d}.csv"), newline="") as handle:
        stored = list(csv.reader(handle))
    compared = 0
    bad = 0
    for row, kept in zip(rows, stored[1:]):
        for mine, theirs in zip(csv_cells(row)[1:], kept[1:]):
            compared += 1
            bad += mine != theirs
    bad += abs(len(rows) - len(stored) + 1) * 13
    return compared, bad


def key_of(row, fields):
    return tuple(row[field] for field in fields)


def groups(rows, fields):
    out = {}
    for row in rows:
        out.setdefault(key_of(row, fields), []).append(row)
    return out


def splits(group, measure):
    return len({row[measure] for row in group}) > 1


def determination(rows):
    out = {}
    for measure in MEASURES:
        out[measure] = "none"
        for fields in KEYS:
            if not any(
                splits(group, measure) for group in groups(rows, fields).values()
            ):
                out[measure] = " + ".join(fields)
                break
    return out


def pin_code(d, r):
    return sum(1 << i for i in range(1 << d) if i >> (d - r) == 0)


def popcount(a):
    return POP16[a & 0xFFFF].astype(np.int8) + POP16[a >> 16].astype(np.int8)


def families(d):
    out = []

    def go(remaining, chosen):
        if remaining == 0:
            out.append(tuple(chosen))
            return
        low = remaining & -remaining
        rest = remaining ^ low
        go(rest, chosen)
        sub = rest
        while True:
            go(rest ^ sub, chosen + [low | sub])
            if sub == 0:
                break
            sub = (sub - 1) & rest

    go((1 << d) - 1, [])
    return out


def profile(d, codes):
    n = 1 << d
    full = np.uint32((1 << n) - 1)
    cube = cubes(d)
    by_size = {}
    for family in families(d):
        by_size.setdefault(len(family), []).append(family)
    by_codim = {}
    for fixed in range(n):
        by_codim.setdefault(fixed.bit_count(), []).append(fixed)
    s = np.zeros(len(codes), np.int8)
    bs = np.zeros(len(codes), np.int8)
    C = np.zeros(len(codes), np.int8)
    for x in range(n):
        m = codes ^ (((codes >> x) & 1) * full)
        nb = np.uint32(sum(1 << (x ^ (1 << i)) for i in range(d)))
        s = np.maximum(s, popcount(m & nb))
        bx = np.zeros(len(codes), np.int8)
        for k in range(1, d + 1):
            hit = np.zeros(len(codes), bool)
            for family in by_size[k]:
                fm = np.uint32(sum(1 << (x ^ b) for b in family))
                hit |= (m & fm) == fm
            if not hit.any():
                break
            bx[hit] = k
        bs = np.maximum(bs, bx)
        cx = np.full(len(codes), d, np.int8)
        for c in range(d - 1, -1, -1):
            hit = np.zeros(len(codes), bool)
            for fixed in by_codim[c]:
                mask = np.uint32(cube[(fixed, x & fixed)])
                v = codes & mask
                hit |= (v == 0) | (v == mask)
            if not hit.any():
                break
            cx[hit] = c
        C = np.maximum(C, cx)
    return s, bs, C


def histogram(s, bs, C):
    keys, counts = np.unique(np.stack([s, bs, C], axis=1), axis=0, return_counts=True)
    return [
        (tuple(int(v) for v in key), int(count)) for key, count in zip(keys, counts)
    ]


def print_profile(label, s, bs, C):
    print(
        f"{label}: designs {len(s)}, C != bs {int((C != bs).sum())}, bs - s > 1 {int((bs - s > 1).sum())}"
    )
    for (a, b, c), count in histogram(s, bs, C):
        print(f"  (s, bs, C) = ({a}, {b}, {c}): {count}")


def table_line(row, d):
    cells = [row["name"], row["genus"], str(row["gf2deg"]), str(row["pop"])]
    if d == 4:
        cells.append(poly_text(row["fingerprint"]))
    return " | ".join(cells + [str(row[m]) for m in MEASURES])


def counts(rows, field):
    out = {}
    for row in rows:
        out[row[field]] = out.get(row[field], 0) + 1
    return ", ".join(f"{key}: {out[key]}" for key in sorted(out))


def influence(d):
    n = 1 << d
    codes = np.arange(1 << n, dtype=np.uint32)
    full = np.uint32((1 << n) - 1)
    total = np.zeros(len(codes), np.int64)
    top = np.zeros(len(codes), np.int8)
    for x in range(n):
        m = codes ^ (((codes >> x) & 1) * full)
        nb = np.uint32(sum(1 << (x ^ (1 << i)) for i in range(d)))
        sx = popcount(m & nb)
        total += sx
        top = np.maximum(top, sx)
    count = len(codes)
    mean = Fraction(int(total.sum()), count * n)
    var = Fraction(int((total * total).sum()), count * n * n) - mean * mean
    full_share = Fraction(int((top == d).sum()), count)
    one_edge = int((total == 2).sum())
    edges = sorted(set((total // 2).tolist()))
    return mean, var, full_share, one_edge, edges


def sample_bits():
    rng = random.Random(SEED)
    codes = []
    for _ in range(TRIALS):
        f = 0
        for i in range(32):
            f |= rng.getrandbits(1) << i
        codes.append(f)
    return np.array(codes, dtype=np.uint32)


def sample_words():
    rng = random.Random(SEED)
    uniform = [rng.getrandbits(32) for _ in range(SAMPLE)]
    thinned = []
    for _ in range(SAMPLE):
        p = LEVELS[rng.randint(0, len(LEVELS) - 1)]
        code = 0
        for i in range(32):
            if rng.randint(0, 31) < p:
                code |= 1 << i
        thinned.append(code)
    return np.array(uniform, dtype=np.uint32), np.array(thinned, dtype=np.uint32)


def main():
    start = time.time()
    store = {}
    for d in DIMS:
        rows = rows_of(d)
        store[d] = rows
        print(f"D={d} classes {len(rows)}, designs {sum(row['orbit'] for row in rows)}")
        compared, bad = csv_diff(rows, d)
        print(f"D={d} csv cells compared {compared}, mismatches {bad}")
        renders = {k: render_index(d, 2 * k - 1) for k in range(1, d + 2)}
        bad = sum(fitted(row["code"], d, renders) != row["fingerprint"] for row in rows)
        print(
            f"D={d} fill polynomial fitted from rendered grids against closed form, mismatches {bad} of {len(rows)}"
        )
    print("D=3 table: name | genus | gf2deg | pop | " + " | ".join(MEASURES))
    for row in store[3]:
        print("  " + table_line(row, 3))
    print(f"D=3 genus split {counts(store[3], 'genus')}")
    print(f"D=3 gf2deg histogram {counts(store[3], 'gf2deg')}")
    deeper = [
        f"{row['name']} (bs {row['bs']}, dt {row['dt']})"
        for row in store[3]
        if row["dt"] > row["bs"]
    ]
    print(f"D=3 dt > bs classes {len(deeper)}: {deeper}")
    for d in DIMS:
        cube = cubes(d)
        for r in range(d + 1):
            got = all_measures(pin_code(d, r), d, cube)
            got["gf2deg"] = top_weight(mobius(pin_code(d, r), d, True))
            six = {got[m] for m in ("s", "bs", "C", "dt", "deg", "gf2deg")}
            print(
                f"D={d} pin of {r} axes: six measures {sorted(six)}, dnf {got['dnf']}, cnf {got['cnf']}"
            )
    d3 = store[3]
    d4 = store[4]
    by_name = {row["name"]: row for row in d4}
    fp3 = groups(d3, ("fingerprint",))
    shared3 = [group for group in fp3.values() if len(group) > 1]
    print(
        f"D=3 distinct fill polynomials {len(fp3)} of {len(d3)}, shared by two or more {len(shared3)}"
    )
    for group in shared3:
        names = ", ".join(
            f"{row['name']} ({row['genus']}, gf2deg {row['gf2deg']})" for row in group
        )
        print(f"D=3 collision {poly_text(group[0]['fingerprint'])}: {names}")
    report = determination(d3)
    print(
        "D=3 coarsest key determining each measure: "
        + ", ".join(f"{m} {report[m]}" for m in MEASURES)
    )
    fp4 = groups(d4, ("fingerprint",))
    shared4 = sum(1 for group in fp4.values() if len(group) > 1)
    print(
        f"D=4 distinct fill polynomials {len(fp4)} of {len(d4)}, shared by two or more {shared4}"
    )
    report = determination(d4)
    print(
        "D=4 coarsest key determining each measure: "
        + ", ".join(f"{m} {report[m]}" for m in MEASURES)
    )
    big = [group for group in groups(d4, KEYS[-1]).values() if len(group) > 1]
    parted = sum(1 for group in big if any(splits(group, m) for m in MEASURES))
    total = sum(1 for group in big for m in MEASURES if splits(group, m))
    print(
        f"D=4 full-key groups of two or more {len(big)}, splitting some measure {parted}, measure-splits {total}"
    )
    pairs = [group for group in big if len(group) == 2 and splits(group, "bs")]
    pairs.sort(key=lambda group: min(row["code"] for row in group))
    first = pairs[0]
    names = " and ".join(row["name"] for row in first)
    split_names = [m for m in MEASURES if splits(first, m)]
    agree = [m for m in MEASURES if not splits(first, m)]
    print(
        f"D=4 size-two full-key groups splitting bs {len(pairs)}, smallest code {first[0]['code']}: {names}"
    )
    print(
        f"  shared key genus {first[0]['genus']}, gf2deg {first[0]['gf2deg']}, pop {first[0]['pop']}, fill {poly_text(first[0]['fingerprint'])}"
    )
    print(f"  split {len(split_names)} of 7 measures {split_names}, agree {agree}")
    print(
        "D=4 witness table: name | genus | gf2deg | pop | fill | "
        + " | ".join(MEASURES)
    )
    for row in first:
        print("  " + table_line(row, 4))
    for d, rows in store.items():
        gaps = [row["name"] for row in rows if row["C"] != row["bs"]]
        print(f"D={d} C != bs classes: {gaps}")
    for d, rows in store.items():
        seps = [
            f"{row['name']} (s {row['s']}, bs {row['bs']}, orbit {row['orbit']})"
            for row in rows
            if row["s"] != row["bs"]
        ]
        print(f"D={d} s != bs classes: {seps}")
    both = d3 + d4
    print(
        f"deg <= s^2 violations across both catalogs: {sum(1 for row in both if row['deg'] > row['s'] ** 2)} of {len(both)}"
    )
    tight = [
        f"{row['name']} (deg {row['deg']}, s {row['s']}, orbit {row['orbit']})"
        for row in d4
        if row["s"] >= 2 and row["deg"] == row["s"] ** 2
    ]
    print(f"D=4 deg = s^2 with s >= 2: {tight}")
    for name in NAMED:
        print(f"{name} anf {anf(by_name[name]['code'], 4)}")
    for d in range(1, 5):
        mean, var, full_share, one_edge, edges = influence(d)
        print(
            f"D={d} influence mean {mean}, variance {var}, P[s(f) = D] {full_share}, designs with one bichromatic edge {one_edge}"
        )
        print(f"D={d} realised bichromatic edge counts {edges}")
    codes = np.arange(1 << 16, dtype=np.uint32)
    print_profile(
        "D=4 exhaustive (s, bs, C) profile by the subcube route", *profile(4, codes)
    )
    print_profile(
        f"D=5 spot check, {TRIALS} uniform designs, seed {SEED}, one bit per draw",
        *profile(5, sample_bits()),
    )
    uniform, thinned = sample_words()
    print_profile(
        f"D=5 spot check, {SAMPLE} uniform designs, seed {SEED}, one word per draw",
        *profile(5, uniform),
    )
    print_profile(
        f"D=5 spot check, {SAMPLE} thinned designs, densities {LEVELS} of 32",
        *profile(5, thinned),
    )
    print(
        f"domain: every design at D = 1..4, {TRIALS} + {SAMPLE} + {SAMPLE} sampled designs at D = 5, {time.time() - start:.1f} s"
    )


if __name__ == "__main__":
    main()
