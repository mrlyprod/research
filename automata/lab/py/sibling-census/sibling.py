import sys
import time
from itertools import combinations, product
from math import gcd

import numpy as np

SEED = 1729
OUT = []


def say(line):
    print(line)
    OUT.append(line)


def check(cond, what):
    if not cond:
        print("FAIL " + what)
        sys.exit(1)


def corner_filled(code, coords):
    idx = 0
    for c in coords:
        idx = 2 * idx + (c & 1)
    return (code >> idx) & 1


def tile(code, dim, side):
    grid = np.zeros((side,) * dim, dtype=np.uint8)
    for coords in product(range(side), repeat=dim):
        grid[coords] = corner_filled(code, coords)
    return grid


def fractal(code, dim, side, level):
    base = tile(code, dim, side)
    out = base
    for _ in range(level - 1):
        out = np.kron(out, base)
    return out


def mask_offsets(code, dim, side, level):
    grid = fractal(code, dim, side, level)
    centre = (side**level - 1) // 2
    cells = np.argwhere(grid == 1) - centre
    keep = np.any(cells != 0, axis=1)
    return [tuple(int(v) for v in row) for row in cells[keep]]


# DICTIONARY


def rule_table(number):
    return np.array([(number >> i) & 1 for i in range(8)], dtype=np.uint8)


def outer_totalistic_elementary():
    found = []
    for number in range(256):
        t = rule_table(number)
        ok = True
        for c in range(2):
            for count in range(3):
                vals = {
                    int(t[4 * l + 2 * c + r])
                    for l in range(2)
                    for r in range(2)
                    if l + r == count
                }
                if len(vals) > 1:
                    ok = False
        if ok:
            found.append(number)
    return found


def dictionary():
    ot = outer_totalistic_elementary()
    say(f"dictionary: outer-totalistic elementary rules {len(ot)}")
    check(len(ot) == 64, "64 outer-totalistic elementary rules")
    check(
        110 not in ot and 30 not in ot and 184 not in ot,
        "110, 30, 184 not outer-totalistic",
    )
    check(
        90 in ot and 150 in ot and 204 in ot and 232 in ot,
        "90, 150, 204, 232 outer-totalistic",
    )
    t = rule_table(110)
    say(
        f"dictionary: rule 110 on 001 -> {t[1]} and on 100 -> {t[4]}, so it is not outer-totalistic"
    )
    say("dictionary: the 64 are " + " ".join(str(n) for n in ot))
    m1 = mask_offsets(1, 1, 3, 1)
    check(
        sorted(m1) == [(-1,), (1,)],
        "bang dim 1, code 1 at side 3 is the elementary mask",
    )
    m7 = mask_offsets(7, 2, 3, 1)
    check(
        len(m7) == 8 and (0, 0) not in m7,
        "bang dim 2, code 7 at side 3 popped is Moore",
    )
    say(
        "dictionary: bang dim 1, code 1 at side 3 is [101]; bang dim 2, code 7 at side 3 popped is the 8-cell Moore mask"
    )


# TOWER


def popcount_le_one_code(dim):
    code = 0
    for i in range(2**dim):
        if i.bit_count() <= 1:
            code |= 1 << i
    return code


def not_all_odd_code(dim):
    return (
        (1 << (2**dim)) - 2 + 1
        if dim == 0
        else (1 << (2**dim)) - 1 - (1 << (2**dim - 1))
    )


def tower():
    codes_a = [popcount_le_one_code(d) for d in range(1, 5)]
    codes_b = [not_all_odd_code(d) for d in range(1, 5)]
    say(f"tower: popcount <= 1 codes {codes_a}; not all odd codes {codes_b}")
    check(codes_a == [3, 7, 23, 279], "popcount <= 1 codes")
    check(codes_b == [1, 7, 127, 32767], "not all odd codes")
    for d in range(1, 5):
        ma = mask_offsets(codes_a[d - 1], d, 3, 1)
        mb = mask_offsets(codes_b[d - 1], d, 3, 1)
        fa = int(tile(codes_a[d - 1], d, 3).sum())
        fb = int(tile(codes_b[d - 1], d, 3).sum())
        check(fa == 2 ** (d - 1) * (d + 2), "popcount <= 1 tile fill 2^(D-1)(D+2)")
        check(fb == 3**d - 1 and len(mb) == 3**d - 1, "not all odd tile is Moore")
        check(
            all(max(abs(v) for v in o) == 1 for o in mb),
            "Moore offsets within Chebyshev 1",
        )
        same = sorted(ma) == sorted(mb)
        check(same == (d <= 2), "masks coincide exactly at D <= 2")
        say(
            f"tower: D={d} popcount<=1 tile {fa} mask {len(ma)}, not-all-odd tile {fb} mask {len(mb)}, masks equal {same}"
        )
    menger = mask_offsets(23, 3, 3, 1)
    check(
        len(menger) == 20 and all(sum(1 for v in o if v == 0) <= 1 for o in menger),
        "Menger mask is 20 cells with at most one zero offset",
    )
    cantor3 = mask_offsets(1, 1, 3, 3)
    check(
        sorted(abs(o[0]) for o in cantor3) == [5, 5, 7, 7, 11, 11, 13, 13],
        "Cantor level 3 is +-5 +-7 +-11 +-13",
    )
    say(
        f"tower: Cantor level 3 offsets {sorted(o[0] for o in cantor3)}, Menger level 1 mask {len(menger)} cells"
    )


# DECOUPLING


def hnf2(rows):
    rows = [list(r) for r in rows if r != (0, 0) and r != [0, 0]]
    while True:
        first = [r for r in rows if r[0] != 0]
        if len(first) <= 1:
            break
        first.sort(key=lambda r: abs(r[0]))
        p = first[0]
        new = [p]
        for r in first[1:]:
            q = r[0] // p[0]
            new.append([r[0] - q * p[0], r[1] - q * p[1]])
        rows = new + [r for r in rows if r[0] == 0]
        rows = [r for r in rows if r != [0, 0]]
    pivot = [r for r in rows if r[0] != 0]
    rest = [r for r in rows if r[0] == 0]
    g = 0
    for r in rest:
        g = gcd(g, abs(r[1]))
    return pivot, g


def lattice_index(offsets, dim):
    if dim == 1:
        g = 0
        for (o,) in offsets:
            g = gcd(g, abs(o))
        return g if g else None
    if not offsets:
        return None
    arr = np.array(offsets, dtype=np.int64)
    v1 = arr[0]
    cross = v1[0] * arr[:, 1] - v1[1] * arr[:, 0]
    nz = np.flatnonzero(cross)
    if len(nz) == 0:
        return None
    basis = [tuple(int(x) for x in v1), tuple(int(x) for x in arr[nz[0]])]
    while True:
        pivot, g = hnf2(basis)
        a, b = pivot[0]
        if a < 0:
            a, b = -a, -b
        x = arr[:, 0]
        y = arr[:, 1]
        bad = x % a != 0
        q = x // a
        bad |= (y - q * b) % g != 0
        idx = np.flatnonzero(bad)
        if len(idx) == 0:
            return a * g
        basis = [(a, b), (0, g), tuple(int(x) for x in arr[idx[0]])]


def eca_ring_step(state, offsets, table_bits, kind):
    if kind == "general":
        idx = state.astype(np.int64)
        for k, o in enumerate(offsets):
            idx |= np.roll(state, -o).astype(np.int64) << (k + 1)
        return table_bits[idx]
    count = np.zeros_like(state, dtype=np.int64)
    for o in offsets:
        count += np.roll(state, -o)
    birth, survive = table_bits
    return np.where(state == 1, survive[count], birth[count]).astype(np.uint8)


def decoupling():
    rng = np.random.default_rng(SEED)
    seen = {}
    rows = []
    for dim in (1, 2):
        for code in range(2 ** (2**dim)):
            for side in (3, 5, 7, 9):
                for level in (1, 2, 3):
                    offs = mask_offsets(code, dim, side, level)
                    key = (dim, frozenset(offs))
                    if key in seen:
                        continue
                    idx = lattice_index(offs, dim)
                    seen[key] = idx
                    rows.append((dim, code, side, level, len(offs), idx))
    hist = {}
    for dim, code, side, level, m, idx in rows:
        k = (dim, "rank<D" if idx is None else str(idx))
        hist[k] = hist.get(k, 0) + 1
    say(
        f"decoupling: distinct masks {len(rows)}; index histogram {sorted(hist.items())}"
    )
    for dim in (1, 2):
        for code in range(2 ** (2**dim)):
            items = [
                f"n{side}L{level}:{m}:{'r' if idx is None else idx}"
                for d, c, side, level, m, idx in rows
                if d == dim and c == code and (dim == 1 or idx != 1)
            ]
            if items:
                say(
                    f"decoupling: D={dim} code={code} side.level:cells:index "
                    + " ".join(items)
                )
    one_d = {
        (side, level): lattice_index(mask_offsets(1, 1, side, level), 1)
        for side in (3, 5, 7, 9)
        for level in (1, 2, 3)
    }
    check(
        one_d[(3, 1)] == 1 and one_d[(3, 2)] == 2 and one_d[(3, 3)] == 1,
        "Cantor tower index 1, 2, 1 at levels 1..3",
    )
    check(
        one_d[(5, 1)] == 2 and one_d[(9, 1)] == 2 and one_d[(7, 1)] == 1,
        "parity tile index by side",
    )
    diag = lattice_index(mask_offsets(9, 2, 3, 1), 2)
    check(diag == 2, "diagonal 4-mask has index 2")
    vn = lattice_index(mask_offsets(6, 2, 3, 1), 2)
    check(vn == 1, "von Neumann 4-mask has index 1")
    tested = 0
    for dim, code, side, level, m, idx in rows:
        if dim != 1 or idx is None or idx == 1:
            continue
        offs = [o[0] for o in mask_offsets(code, dim, side, level)]
        small = [o // idx for o in offs]
        n = 12 * idx * 5
        for kind in ("general", "life"):
            if kind == "general" and m > 12:
                continue
            if kind == "general":
                table = rng.integers(0, 2, size=2 ** (m + 1)).astype(np.uint8)
                bits = table
            else:
                bits = (
                    rng.integers(0, 2, size=m + 1).astype(np.uint8),
                    rng.integers(0, 2, size=m + 1).astype(np.uint8),
                )
            x = rng.integers(0, 2, size=n).astype(np.uint8)
            ys = [x[i::idx].copy() for i in range(idx)]
            for _ in range(40):
                x = eca_ring_step(x, offs, bits, kind)
                ys = [eca_ring_step(y, small, bits, kind) for y in ys]
                for i in range(idx):
                    check(
                        np.array_equal(x[i::idx], ys[i]),
                        f"interleaving D=1 code={code} side={side} level={level} kind={kind}",
                    )
            tested += 1
    say(
        f"decoupling: interleaving equality checked on {tested} (mask, kind) pairs at D=1, 40 steps each"
    )


# COMPOSITES

BLOCKS = np.arange(512, dtype=np.int64)
ROW = [(BLOCKS >> (6 - 3 * i)) & 7 for i in range(3)]
CENTRE = (BLOCKS >> 4) & 1
POP9 = np.array([int(b).bit_count() for b in range(512)], dtype=np.int64)
OUTER = POP9 - CENTRE
RULES = np.arange(256, dtype=np.int64)
TABLE = ((RULES[:, None] >> np.arange(8)[None, :]) & 1).astype(np.int64)


def life_like_table(birth, survive):
    b = np.isin(OUTER, list(birth))
    s = np.isin(OUTER, list(survive))
    return np.where(CENTRE == 1, s, b).astype(np.uint8)


def composites():
    idx = 4 * TABLE[:, ROW[0]] + 2 * TABLE[:, ROW[1]] + TABLE[:, ROW[2]]
    h = ((RULES[None, :, None] >> idx[:, None, :]) & 1).astype(np.uint8)
    flat = h.reshape(65536, 512)
    packed = np.packbits(flat, axis=1).view(np.uint64)
    distinct = np.unique(packed, axis=0).shape[0]
    say(f"composites: 65536 ordered pairs give {distinct} distinct 9-input rules")
    xor9 = life_like_table((1, 3, 5, 7), (0, 2, 4, 6, 8))
    check(
        np.array_equal(h[150, 150], xor9),
        "150 then 150 is the nine-cell XOR B1357/S02468",
    )
    check(np.array_equal(h[105, 105], xor9), "105 then 105 is the nine-cell XOR too")
    for f in range(256):
        check(
            np.array_equal(h[f, 204], TABLE[f, ROW[1]].astype(np.uint8)),
            "f then 204 is f on the middle row",
        )
        check(
            np.array_equal(h[f, 170], TABLE[f, ROW[2]].astype(np.uint8)),
            "f then 170 is f on the row below",
        )
        check(
            np.array_equal(h[f, 240], TABLE[f, ROW[0]].astype(np.uint8)),
            "f then 240 is f on the row above",
        )
        check(
            np.array_equal(h[f, 51], 1 - h[f, 204])
            and np.array_equal(h[f, 85], 1 - h[f, 170])
            and np.array_equal(h[f, 15], 1 - h[f, 240]),
            "the complementing three negate",
        )
    say(
        "composites: f then 204 / 170 / 240 is f on the middle / lower / upper row for all 256 f; then 51 / 85 / 15 is its negation"
    )
    life = life_like_table((3,), (2, 3))
    hit = np.flatnonzero(np.all(flat == life[None, :], axis=1))
    check(len(hit) == 0, "Life is not a composite")
    say(f"composites: B3/S23 occurs among the 65536 composites {len(hit)} times")
    ot = np.ones(65536, dtype=bool)
    for c in range(2):
        for n in range(9):
            cols = np.flatnonzero((CENTRE == c) & (OUTER == n))
            sub = flat[:, cols]
            ot &= sub.min(axis=1) == sub.max(axis=1)
    say(f"composites: outer-totalistic composites {int(ot.sum())} ordered pairs")
    named = {}
    for k in np.flatnonzero(ot):
        f, g = divmod(int(k), 256)
        birth = tuple(
            n
            for n in range(9)
            if flat[k, np.flatnonzero((CENTRE == 0) & (OUTER == n))[0]]
        )
        survive = tuple(
            n
            for n in range(9)
            if flat[k, np.flatnonzero((CENTRE == 1) & (OUTER == n))[0]]
        )
        named.setdefault((birth, survive), []).append((f, g))
    say(f"composites: distinct life-like composites {len(named)}")
    for (birth, survive), pairs in sorted(
        named.items(), key=lambda kv: (len(kv[1]), kv[0])
    ):
        bs = "B" + "".join(map(str, birth)) + "/S" + "".join(map(str, survive))
        shown = " ".join(f"{f}.{g}" for f, g in pairs[:12])
        say(
            f"composites: {bs} from {len(pairs)} pairs f.g: {shown}"
            + (" ..." if len(pairs) > 12 else "")
        )
    cell = [(i, j) for i in range(3) for j in range(3)]

    def perm_blocks(mapping):
        out = np.zeros(512, dtype=np.int64)
        for b in range(512):
            v = 0
            for i, j in cell:
                si, sj = mapping(i, j)
                bit = (b >> (8 - (3 * si + sj))) & 1
                v |= bit << (8 - (3 * i + j))
            out[b] = v
        return out

    syms = [
        lambda i, j: (j, i),
        lambda i, j: (i, 2 - j),
        lambda i, j: (2 - i, j),
        lambda i, j: (j, 2 - i),
    ]
    transpose = perm_blocks(syms[0])
    tsym = np.all(flat == flat[:, transpose], axis=1)
    full = tsym.copy()
    for m in syms[1:]:
        pb = perm_blocks(m)
        full &= np.all(flat == flat[:, pb], axis=1)
    say(
        f"composites: transpose-symmetric {int(tsym.sum())} pairs, full dihedral symmetry {int(full.sum())} pairs"
    )
    swapped = h.transpose(1, 0, 2).reshape(65536, 512)[:, transpose]
    same_order = np.all(flat == swapped, axis=1)
    say(
        f"composites: rows-first equals columns-first for {int(same_order.sum())} of 65536 ordered pairs"
    )
    check(
        bool(ot[150 * 256 + 150]) and bool(full[150 * 256 + 150]),
        "150.150 is life-like and dihedral",
    )
    return h


# CANTOR LIFE

CANTOR = (-13, -11, -7, -5, 5, 7, 11, 13)
NAMED = (
    ("B3/S23", (3,), (2, 3)),
    ("B36/S23", (3, 6), (2, 3)),
    ("B2/S", (2,), ()),
    ("B3/S012345678", (3,), tuple(range(9))),
    ("B3678/S34678", (3, 6, 7, 8), (3, 4, 6, 7, 8)),
    ("B1357/S02468", (1, 3, 5, 7), (0, 2, 4, 6, 8)),
)
RING = 1024
SOUP_STEPS = 2000
DENSITIES = (0.1, 0.25, 0.5)
SOUP_SEEDS = 5
WIDTH = 14
HORIZON = 256
REACH = 13


def masks_of(birth, survive, m):
    b = np.zeros(m + 1, dtype=np.uint8)
    s = np.zeros(m + 1, dtype=np.uint8)
    b[list(birth)] = 1
    s[list(survive)] = 1
    return b, s


def count_1d(x, offsets):
    c = np.zeros(x.shape, dtype=np.uint8)
    for o in offsets:
        c += np.roll(x, -o, axis=-1)
    return c


def life_step(x, c, b, s):
    return np.where(x == 1, s[c], b[c]).astype(np.uint8)


def fate_name(period, disp):
    if disp != 0:
        return "mover"
    return "still" if period == 1 else "oscillator"


def cantor_soups():
    rng = np.random.default_rng(SEED)
    for name, birth, survive in NAMED:
        b, s = masks_of(birth, survive, 8)
        for dens in DENSITIES:
            x = (rng.random((SOUP_SEEDS, RING)) < dens).astype(np.uint8)
            seen = [{} for _ in range(SOUP_SEEDS)]
            fates = ["undecided"] * SOUP_SEEDS
            for t in range(SOUP_STEPS + 1):
                for i in range(SOUP_SEEDS):
                    if fates[i] != "undecided":
                        continue
                    if not x[i].any():
                        fates[i] = "dies"
                        continue
                    key = x[i].tobytes()
                    if key in seen[i]:
                        p = t - seen[i][key]
                        fates[i] = "fixed" if p == 1 else f"period {p}"
                    seen[i][key] = t
                if t < SOUP_STEPS:
                    x = life_step(x, count_1d(x, CANTOR), b, s)
            final = [int(r.sum()) for r in x]
            hist = np.bincount(count_1d(x, CANTOR).ravel(), minlength=9)
            say(
                f"cantor {name} density {dens}: fates {fates}; final live {final} of {RING}; count histogram {hist.tolist()}"
            )


def cantor_xor_period():
    kernel = {0} | set(CANTOR)

    def power_two(k):
        out = {}
        for o in kernel:
            r = (o * (1 << k)) % RING
            out[r] = out.get(r, 0) ^ 1
        return {r for r, v in out.items() if v}

    check(
        power_two(8) == {0},
        "Cantor XOR kernel to the 256th power is 1 on the ring 1024",
    )
    check(power_two(7) != {0}, "Cantor XOR kernel to the 128th power is not 1")
    say(
        f"cantor B1357/S02468: kernel^128 has support {sorted(power_two(7))} and kernel^256 = 1 on the ring 1024, so the period divides 256 and is 256 for a generic soup"
    )


def cantor_patterns():
    seeds = [1] + [
        (1 << (w - 1)) | 1 | (rest << 1)
        for w in range(2, WIDTH + 1)
        for rest in range(1 << (w - 2))
    ]
    check(len(seeds) == 8192, "8192 seeds of width at most 14")
    b, s = masks_of((3,), (2, 3), 8)
    narrow = WIDTH + 2 * REACH * 16 + 2
    wide = WIDTH + 2 * REACH * HORIZON + 2
    x = np.zeros((len(seeds), narrow), dtype=np.uint8)
    start = (narrow - WIDTH) // 2
    for i, code in enumerate(seeds):
        for k in range(WIDTH):
            x[i, start + k] = (code >> k) & 1
    ids = np.arange(len(seeds))
    seen = [{} for _ in seeds]
    fates = {}
    witness = {}
    shift = 0
    t = 0
    while t <= HORIZON and len(ids):
        if t == 16:
            pad = (wide - narrow) // 2
            x = np.pad(x, ((0, 0), (pad, pad)))
            shift = pad
        alive = x.any(axis=1)
        left = np.argmax(x, axis=1)
        right = x.shape[1] - 1 - np.argmax(x[:, ::-1], axis=1)
        keep = []
        for row, i in enumerate(ids):
            if not alive[row]:
                fates[i] = ("death", 0, 0, 0)
                continue
            key = x[row, left[row] : right[row] + 1].tobytes()
            store = seen[i]
            if key in store:
                t0, l0 = store[key]
                p = t - t0
                v = int(left[row]) - shift - l0
                fates[i] = (fate_name(p, v), p, v, int(x[row].sum()))
                witness[i] = key
                continue
            store[key] = (t, int(left[row]) - shift)
            keep.append(row)
        if t == HORIZON:
            for row in keep:
                fates[ids[row]] = ("undecided", 0, 0, int(x[row].sum()))
            break
        x = x[keep]
        ids = ids[keep]
        x = life_step(x, count_1d(x, CANTOR), b, s)
        t += 1
    tally = {}
    for i, (kind, p, v, n) in fates.items():
        tally[kind] = tally.get(kind, 0) + 1
    say(
        f"cantor B3/S23 patterns: {sorted(tally.items())} over {len(seeds)} seeds, horizon {HORIZON}"
    )

    def smallest(kind):
        cands = [(n, p, i) for i, (k, p, v, n) in fates.items() if k == kind]
        if not cands:
            return None
        n, p, i = min(cands)
        cells = [k for k, c in enumerate(witness[i]) if c]
        return n, p, i, cells

    for kind in ("still", "oscillator"):
        w = smallest(kind)
        if w:
            n, p, i, cells = w
            say(
                f"cantor B3/S23 smallest {kind}: {n} cells, period {p}, cells {cells}, from seed {seeds[i]:b}"
            )
        else:
            say(f"cantor B3/S23 smallest {kind}: not found")
    movers = [(i, p, v, n) for i, (k, p, v, n) in fates.items() if k == "mover"]
    for i, p, v, n in sorted(movers, key=lambda m: (m[3], m[1]))[:10]:
        cells = [k for k, c in enumerate(witness[i]) if c]
        say(
            f"cantor B3/S23 mover: {n} cells, period {p}, displacement {v}, cells {cells}, from seed {seeds[i]:b}"
        )
    if not movers:
        say("cantor B3/S23 mover: not found over width <= 14")
    periods = sorted({p for k, p, v, n in fates.values() if k == "oscillator"})
    say(f"cantor B3/S23 oscillator periods found {periods}")


# MENGER LIFE

MENGER = tuple(
    o for o in product((-1, 0, 1), repeat=3) if sum(1 for v in o if v == 0) <= 1
)
TORUS = 32
MENGER_STEPS = 200
MENGER_DENSITIES = (0.15, 0.3)
MENGER_SEEDS = 2
FIELD = 24
BOX = 5
MOVER_SEEDS = 200
MOVER_STEPS = 128


def axis_sum(x, axis):
    return x + np.roll(x, 1, axis) + np.roll(x, -1, axis)


def count_menger(x):
    box = axis_sum(axis_sum(axis_sum(x, -1), -2), -3)
    faces = sum(np.roll(x, d, axis) for axis in (-1, -2, -3) for d in (1, -1))
    return box - x - faces


def menger_rules():
    subsets = [c for k in (1, 2) for c in combinations((3, 4, 5, 6), k)]
    return [(b, s) for b in subsets for s in subsets] + [((3,), (2, 3))]


def rule_name(birth, survive):
    return "B" + "".join(map(str, birth)) + "/S" + "".join(map(str, survive))


def menger_soups():
    rng = np.random.default_rng(SEED)
    rules = menger_rules()
    check(
        len(rules) == 101 and len(MENGER) == 20,
        "100 grid rules plus B3/S23 on the 20-cell mask",
    )
    verdicts = {}
    for birth, survive in rules:
        b, s = masks_of(birth, survive, 20)
        x = np.zeros(
            (len(MENGER_DENSITIES) * MENGER_SEEDS, TORUS, TORUS, TORUS), dtype=np.uint8
        )
        for j, dens in enumerate(MENGER_DENSITIES):
            for k in range(MENGER_SEEDS):
                x[j * MENGER_SEEDS + k] = rng.random((TORUS,) * 3) < dens
        seen = [{} for _ in range(x.shape[0])]
        fates = ["active"] * x.shape[0]
        for t in range(MENGER_STEPS + 1):
            for i in range(x.shape[0]):
                if fates[i] != "active":
                    continue
                if not x[i].any():
                    fates[i] = "dies"
                    continue
                key = x[i].tobytes()
                if key in seen[i]:
                    p = t - seen[i][key]
                    fates[i] = "fixed" if p == 1 else f"period {p}"
                seen[i][key] = t
            if t < MENGER_STEPS:
                x = life_step(x, count_menger(x), b, s)
        dens = [round(float(r.mean()), 3) for r in x]
        fates = [
            "explodes" if f == "active" and d > 0.4 else f for f, d in zip(fates, dens)
        ]
        verdicts[(birth, survive)] = (fates, dens)
    tally = {}
    for fates, dens in verdicts.values():
        for f in fates:
            k = f.split()[0]
            tally[k] = tally.get(k, 0) + 1
    say(f"menger soups: 101 rules x 4 runs, run fates {sorted(tally.items())}")
    quiet = []
    active = []
    for (birth, survive), (fates, dens) in verdicts.items():
        is_quiet = all(f != "active" and f != "explodes" for f in fates)
        say(
            f"menger {rule_name(birth, survive)}: {fates} density {dens}"
            + (" QUIET" if is_quiet else "")
        )
        if is_quiet:
            quiet.append((birth, survive))
        elif all(f == "active" for f in fates):
            active.append(max(dens))
    say(
        f"menger soups: {len(quiet)} rules quiet at both densities; {len(active)} rules active in all four runs with final density between {min(active)} and {max(active)}; 0 runs above 0.4"
    )
    return rules, quiet


def bbox(x):
    axes = [
        np.flatnonzero(x.any(axis=tuple(a for a in range(3) if a != k)))
        for k in range(3)
    ]
    return [int(a[0]) for a in axes], [int(a[-1]) for a in axes]


def menger_movers(rules, quiet):
    rng = np.random.default_rng(SEED + 1)
    found = {}
    tallies = {}
    for birth, survive in rules:
        b, s = masks_of(birth, survive, 20)
        tally = {}
        for seed_no in range(MOVER_SEEDS):
            box = rng.integers(0, 2, size=(BOX,) * 3).astype(np.uint8)
            if not box.any():
                tally["death"] = tally.get("death", 0) + 1
                continue
            x = np.zeros((FIELD,) * 3, dtype=np.uint8)
            lo = (FIELD - BOX) // 2
            x[lo : lo + BOX, lo : lo + BOX, lo : lo + BOX] = box
            origin = np.zeros(3, dtype=np.int64)
            seen = {}
            kind = "undecided"
            for t in range(MOVER_STEPS + 1):
                if not x.any():
                    kind = "death"
                    break
                mn, mx = bbox(x)
                extent = max(hi - lo_ + 1 for lo_, hi in zip(mn, mx))
                if extent > FIELD - 2:
                    kind = "growing"
                    break
                shift = [FIELD // 2 - (lo_ + hi + 1) // 2 for lo_, hi in zip(mn, mx)]
                x = np.roll(x, shift, axis=(0, 1, 2))
                origin -= np.array(shift)
                mn = [m + sh for m, sh in zip(mn, shift)]
                mx = [m + sh for m, sh in zip(mx, shift)]
                crop = x[mn[0] : mx[0] + 1, mn[1] : mx[1] + 1, mn[2] : mx[2] + 1]
                key = (crop.shape, crop.tobytes())
                pos = tuple(int(m) for m in np.array(mn) + origin)
                if key in seen:
                    t0, pos0 = seen[key]
                    p = t - t0
                    v = tuple(a - c for a, c in zip(pos, pos0))
                    kind = fate_name(p, 1 if any(v) else 0)
                    if kind == "mover":
                        cells = [
                            tuple(int(c) for c in cell)
                            for cell in np.argwhere(crop == 1)
                        ]
                        entry = (int(crop.sum()), p, v, cells, seed_no)
                        if (birth, survive) not in found or entry < found[
                            (birth, survive)
                        ]:
                            found[(birth, survive)] = entry
                    break
                seen[key] = (t, pos)
                if t < MOVER_STEPS:
                    x = life_step(x, count_menger(x), b, s)
            tally[kind] = tally.get(kind, 0) + 1
        tallies[(birth, survive)] = tally
        if (
            len(tally) > 1
            or (birth, survive) in quiet
            or (birth, survive) == ((3,), (2, 3))
        ):
            say(
                f"menger movers {rule_name(birth, survive)}: {sorted(tally.items())}"
                + (" QUIET" if (birth, survive) in quiet else "")
            )
    kinds = {}
    for tally in tallies.values():
        for k, v in tally.items():
            kinds[k] = kinds.get(k, 0) + v
    say(f"menger movers: seed fates over {len(rules)} rules {sorted(kinds.items())}")
    say(
        f"menger movers: {len(found)} of {len(rules)} rules carry a mover from {MOVER_SEEDS} seeds in a {BOX}^3 box, {sum(1 for r in found if r in quiet)} among the {len(quiet)} quiet rules"
    )
    for (birth, survive), (n, p, v, cells, seed_no) in sorted(
        found.items(), key=lambda kv: kv[1][:2]
    ):
        say(
            f"menger mover {rule_name(birth, survive)}: {n} cells, period {p}, displacement {v}, seed {seed_no}, cells {cells}"
        )


def menger_block():
    b, s = masks_of((3,), (2, 3), 20)
    depth = 32
    x = np.zeros((8, 8, depth), dtype=np.uint8)
    x[3:5, 3:5, depth // 2] = 1
    block = x[:, :, depth // 2].copy()
    row = np.zeros(depth, dtype=np.uint8)
    row[depth // 2] = 1
    for t in range(1, 9):
        x = life_step(x, count_menger(x), b, s)
        row = np.roll(row, 1) ^ np.roll(row, -1)
        profile = x.sum(axis=(0, 1))
        check(
            np.array_equal((profile == 4).astype(np.uint8), row)
            and np.all((profile == 0) | (profile == 4)),
            f"block stack is rule 90 at t={t}",
        )
        for z in np.flatnonzero(row):
            check(np.array_equal(x[:, :, z], block), f"block shape kept at t={t}")
    say(
        "menger B3/S23 block: a 2x2 plane block becomes a rule 90 stack of blocks along the normal, checked to t=8"
    )


# PRODUCTS

TORUS2 = 256
PRODUCT_STEPS = 512
PRODUCT_G = (204, 170, 110, 54, 30, 90)
FIELD2 = 128
SEED_STEPS = 96


def eca_axis(x, table, axis):
    l = np.roll(x, 1, axis)
    r = np.roll(x, -1, axis)
    return table[4 * l.astype(np.int64) + 2 * x + r]


def composite_step(x, f_table, g_table):
    return eca_axis(eca_axis(x, f_table, -1), g_table, -2)


def products_torus():
    rng = np.random.default_rng(SEED + 2)
    f_table = rule_table(110)
    for g in PRODUCT_G:
        g_table = rule_table(g)
        for seed_no in range(2):
            x = (rng.random((TORUS2, TORUS2)) < 0.3).astype(np.uint8)
            dens = [float(x.mean())]
            churn = []
            seen = {x.tobytes(): 0}
            period = 0
            for t in range(1, PRODUCT_STEPS + 1):
                nxt = composite_step(x, f_table, g_table)
                churn.append(float((nxt != x).mean()))
                x = nxt
                dens.append(float(x.mean()))
                key = x.tobytes()
                if key in seen and not period:
                    period = t - seen[key]
                seen[key] = t
            say(
                f"product 110.{g} seed {seed_no}: density start {dens[0]:.3f} min {min(dens):.3f} max {max(dens):.3f} end {dens[-1]:.3f}, churn {np.mean(churn[-64:]):.3f}, periodic {'period ' + str(period) if period else 'no'}"
            )


def products_seeds():
    table = rule_table(110)
    n = 511
    x = np.zeros((n, FIELD2, FIELD2), dtype=np.uint8)
    at = FIELD2 - 16
    for code in range(1, 512):
        for i in range(3):
            for j in range(3):
                x[code - 1, at + i, at + j] = (code >> (8 - (3 * i + j))) & 1

    def corner(x):
        rows = x.any(axis=2)
        cols = x.any(axis=1)
        return (
            np.argmax(rows, axis=1),
            np.argmax(cols, axis=1),
            FIELD2 - 1 - np.argmax(rows[:, ::-1], axis=1),
            FIELD2 - 1 - np.argmax(cols[:, ::-1], axis=1),
        )

    y0, x0, y1, x1 = corner(x)
    pops = [x.sum(axis=(1, 2))]
    for t in range(1, SEED_STEPS + 1):
        x = composite_step(x, table, table)
        ya, xa, yb, xb = corner(x)
        check(
            np.all(ya == y0 - t)
            and np.all(xa == x0 - t)
            and np.all(yb == y1)
            and np.all(xb == x1),
            f"110.110 bounding box law at t={t}",
        )
        pops.append(x.sum(axis=(1, 2)))
    final = pops[-1]
    say(
        f"product 110.110 seeds: all 511 grow, upper-left corner moves (-1,-1) per step and the lower-right corner is fixed, checked to t={SEED_STEPS}"
    )
    say(
        f"product 110.110 seeds: population at t={SEED_STEPS} min {int(final.min())} median {int(np.median(final))} max {int(final.max())}; single cell {int(final[(1 << 4) - 1])}"
    )


def main():
    t0 = time.time()
    say(f"seed {SEED}")
    dictionary()
    tower()
    decoupling()
    composites()
    cantor_soups()
    cantor_xor_period()
    cantor_patterns()
    menger_block()
    rules, quiet = menger_soups()
    menger_movers(rules, quiet)
    products_torus()
    products_seeds()
    say(f"elapsed {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
