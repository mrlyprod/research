import time
from math import comb, log2

import numpy as np

BINOM = np.array([comb(8, k) for k in range(9)], dtype=np.int64)
POP8 = np.array([i.bit_count() for i in range(256)], dtype=np.int64)
POP9 = np.array([i.bit_count() for i in range(512)], dtype=np.int64)
SYM = ((np.arange(512, dtype=np.int64)[:, None] >> POP8[None, :]) & 1).astype(np.uint8)
SUBSET = (np.arange(512, dtype=np.int64)[:, None] >> np.arange(9)[None, :]) & 1
FILLS = SUBSET @ BINOM
BITS = (np.arange(512, dtype=np.int64)[:, None] >> np.arange(9)[None, :]) & 1
EDGE = SUBSET @ np.array([comb(7, k - 1) if k else 0 for k in range(9)], dtype=np.int64)
ORDER = np.argsort(POP8, kind="stable")
STARTS = np.concatenate([[0], np.cumsum(BINOM)[:-1]])
NAMED = (
    ((3,), (2, 3)),
    ((3, 6), (2, 3)),
    ((2,), ()),
    ((3,), (0, 1, 2, 3, 4, 5, 6, 7, 8)),
    ((3, 6, 7, 8), (3, 4, 6, 7, 8)),
    ((3, 5, 6, 7, 8), (5, 6, 7, 8)),
    ((1, 3, 5, 7), (1, 3, 5, 7)),
    ((4, 6, 7, 8), (3, 5, 6, 7, 8)),
    ((3, 6), (1, 2, 5)),
    ((3, 6, 8), (2, 4, 5)),
    ((), ()),
    ((0, 1, 2, 3, 4, 5, 6, 7, 8), (0, 1, 2, 3, 4, 5, 6, 7, 8)),
    ((), (0, 1, 2, 3, 4, 5, 6, 7, 8)),
    ((0, 1, 2, 3, 4, 5, 6, 7, 8), ()),
)
LIFE = ((3,), (2, 3))
HORIZON = 64
POWERS = (1, 2, 4, 8, 16, 32, 64)
A071053 = (1, 3, 3, 5, 3, 9, 5, 11, 3, 9, 9, 15, 5, 15, 11, 21, 3, 9, 9, 15)
A246035 = (
    1,
    9,
    9,
    25,
    9,
    81,
    25,
    121,
    9,
    81,
    81,
    225,
    25,
    225,
    121,
    441,
    9,
    81,
    81,
    225,
)
A160239 = (
    1,
    8,
    8,
    24,
    8,
    64,
    24,
    112,
    8,
    64,
    64,
    192,
    24,
    192,
    112,
    416,
    8,
    64,
    64,
    192,
)


def moore(d):
    mask = np.ones(3**d, dtype=np.uint8)
    center = 0
    for _ in range(d):
        center = center * 3 + 1
    mask[center] = 0
    return mask.reshape([3] * d)


def residue_corners(d, base=2):
    return [
        [(i // base ** (d - 1 - j)) % base for j in range(d)] for i in range(base**d)
    ]


def tile_from_code(code, d, side=3):
    corners = residue_corners(d)
    filled = [(code >> i) & 1 for i in range(len(corners))]
    out = np.zeros([side] * d, dtype=np.uint8)
    for flat in range(side**d):
        v = tuple((flat // side ** (d - 1 - j)) % side for j in range(d))
        corner = 0
        for coord in v:
            corner = corner * 2 + coord % 2
        out[v] = filled[corner]
    return out


def subset_code(counts):
    code = 0
    for k in counts:
        code |= 1 << k
    return code


def rule_table(birth, survive):
    return np.concatenate([SYM[subset_code(birth)], SYM[subset_code(survive)]])


def rule_name(birth, survive):
    return (
        "B"
        + "".join(str(k) for k in sorted(birth))
        + "/S"
        + "".join(str(k) for k in sorted(survive))
    )


def design_code(table):
    code = 0
    for x in np.flatnonzero(table):
        code |= 1 << int(x)
    return code


def block_table(table):
    out = np.zeros(512, dtype=np.uint8)
    for y in range(512):
        bits = [(y >> (8 - p)) & 1 for p in range(9)]
        centre = bits[4]
        count = sum(bits) - centre
        out[y] = table[centre * 256 + (2**count - 1)]
    return out


def anf(tables):
    a = np.atleast_2d(tables).astype(np.uint8).copy()
    for b in range(a.shape[1].bit_length() - 1):
        step = 1 << b
        view = a.reshape(a.shape[0], -1, 2, step)
        view[:, :, 1, :] ^= view[:, :, 0, :]
    return a


def degrees(tables):
    a = anf(tables)
    pops = POP9 if a.shape[1] == 512 else POP8
    return np.where(a.astype(bool), pops[None, :], -1).max(axis=1)


def walsh(table):
    w = table.astype(np.int64).copy()
    for b in range(9):
        step = 1 << b
        view = w.reshape(-1, 2, step)
        low = view[:, 0, :].copy()
        high = view[:, 1, :].copy()
        view[:, 0, :] = low + high
        view[:, 1, :] = low - high
    return w


def level_sums(table):
    sums = np.zeros(10, dtype=np.int64)
    np.add.at(sums, POP9, walsh(table))
    return sums


def level_sums_closed(table):
    weights = np.bincount(POP9[table.astype(bool)], minlength=10)
    poly = np.zeros(10, dtype=np.int64)
    for w, count in enumerate(weights):
        if not count:
            continue
        term = np.array([1], dtype=np.int64)
        for _ in range(9 - w):
            term = np.convolve(term, [1, 1])
        for _ in range(w):
            term = np.convolve(term, [1, -1])
        poly += count * term
    return poly


def is_pin(table):
    support = np.flatnonzero(table)
    if support.size == 0:
        return False
    bits = ((support[:, None] >> np.arange(9)[None, :]) & 1).sum(axis=0)
    fixed = int(((bits == 0) | (bits == support.size)).sum())
    return support.size == 1 << (9 - fixed)


def is_level_set(table):
    for t in range(512):
        classes = POP9[np.arange(512) ^ t]
        lows = np.zeros(10, dtype=np.int64)
        highs = np.zeros(10, dtype=np.int64)
        np.add.at(lows, classes, table)
        np.add.at(highs, classes, 1 - table)
        if not np.any((lows > 0) & (highs > 0)):
            return True
    return False


def pin_brute():
    out = np.zeros((512, 512), dtype=bool)
    for b in range(512):
        tables = np.concatenate(
            [np.repeat(SYM[b][None, :], 512, axis=0), SYM], axis=1
        ).astype(np.int64)
        sizes = tables.sum(axis=1)
        fixed = ((tables @ BITS == 0) | (tables @ BITS == sizes[:, None])).sum(axis=1)
        out[b] = (sizes > 0) & (sizes == (1 << (9 - fixed)))
    return out


def pin_grid():
    fill = FILLS[:, None] + FILLS[None, :]
    outer = EDGE[:, None] + EDGE[None, :]
    centre = np.broadcast_to(FILLS[None, :], fill.shape)
    fixed = 8 * ((outer == 0) | (outer == fill)) + ((centre == 0) | (centre == fill))
    return fill == (1 << (9 - fixed))


def level_set_grid():
    flags = np.zeros((512, 512), dtype=bool)
    xs = np.arange(512, dtype=np.int64)
    classes = POP9[xs[None, :] ^ xs[:, None]]
    for chosen in range(1024):
        levels = ((chosen >> classes) & 1).astype(np.int64)
        sums_b = np.add.reduceat(levels[:, :256][:, ORDER], STARTS, axis=1)
        sums_s = np.add.reduceat(levels[:, 256:][:, ORDER], STARTS, axis=1)
        good_b = ((sums_b == 0) | (sums_b == BINOM[None, :])).all(axis=1)
        good_s = ((sums_s == 0) | (sums_s == BINOM[None, :])).all(axis=1)
        keep = good_b & good_s
        if not keep.any():
            continue
        codes_b = ((sums_b > 0) << np.arange(9)[None, :]).sum(axis=1)
        codes_s = ((sums_s > 0) << np.arange(9)[None, :]).sum(axis=1)
        flags[codes_b[keep], codes_s[keep]] = True
    return flags


def degree_grid():
    grid = np.zeros((512, 512), dtype=np.int8)
    for b in range(512):
        tables = np.concatenate([np.repeat(SYM[b][None, :], 512, axis=0), SYM], axis=1)
        grid[b] = degrees(tables)
    return grid


def step(grid, birth, survive, mask):
    pad = np.pad(grid, 1)
    counts = np.zeros(grid.shape, dtype=np.int64)
    for di in range(3):
        for dj in range(3):
            if mask[di, dj]:
                counts += pad[di : di + grid.shape[0], dj : dj + grid.shape[1]]
    return np.where(grid == 1, np.isin(counts, survive), np.isin(counts, birth)).astype(
        np.uint8
    )


def step150(row):
    pad = np.pad(row, 1)
    return ((pad[:-2] + pad[1:-1] + pad[2:]) % 2).astype(np.uint8)


def t1_moore():
    print("T1 moore = carpet")
    for d in (1, 2, 3):
        code = (1 << (1 << d)) - 1 - (1 << ((1 << d) - 1))
        tile = tile_from_code(code, d)
        assert np.array_equal(tile, moore(d)), d
        print(
            f"  D={d} code {code} level-1 side-3 tile = moore({d}), fill {int(tile.sum())} of {3**d}"
        )
    assert (1 << 4) - 1 - (1 << 3) == 7
    print(
        "  the plane case is bang dim 2, code 7, fill 8 of 9, dimension log(8)/log(3) = 1.892789"
    )


def t2_life():
    print("T2 life as a design")
    table = rule_table(*LIFE)
    fill = int(table.sum())
    code = design_code(table)
    lam = fill / 512
    print(
        f"  {rule_name(*LIFE)} fill {fill} of 512 = {comb(8, 3)} + {comb(8, 2)} + {comb(8, 3)}, lambda {fill}/512"
    )
    print(f"  code {fill} bits, hex {code:0128x}")
    print(
        f"  code under the 3x3 block order, hex {design_code(block_table(table)):0128x}"
    )
    assert code.bit_count() == fill
    deg = int(degrees(table[None, :])[0])
    print(f"  GF(2) degree {deg}, monomials {int(anf(table[None, :]).sum())}")
    sums = level_sums(table)
    closed = level_sums_closed(table)
    assert np.array_equal(sums, closed)
    print(
        "  walsh level sums " + " ".join(f"S{k}={int(v)}" for k, v in enumerate(sums))
    )
    energy = int((walsh(table).astype(np.int64) ** 2).sum())
    assert energy == 512 * fill
    print(
        f"  parseval sum of W(S)^2 over the 512 subsets {energy} = 512 * fill = {512 * fill}"
    )
    print(
        f"  popcount 4 both ways: f(c=1,|n|=3) = {int(table[256 + 7])}, f(c=0,|n|=4) = {int(table[15])}"
    )
    print(f"  level set {is_level_set(table)}, pin {is_pin(table)}, genus compound")
    print(
        f"  dimension log2(fill) = {log2(fill):.6f} = 9 + log2({fill}/512) = {9 + log2(lam):.6f}"
    )
    assert abs(log2(fill) - (9 + log2(lam))) < 1e-12
    print(f"  outer-totalistic rules 2^18 = {1 << 18}, totalistic 2^10 = {1 << 10}")


def t3_census(pins, levels):
    print("T3 the lambda census")
    dist = np.zeros(257, dtype=np.int64)
    dist[0] = 1
    for k in range(9):
        shifted = np.zeros_like(dist)
        shifted[BINOM[k] :] = dist[: 257 - BINOM[k]]
        dist = dist + shifted
    hist = np.convolve(dist, dist)
    brute = np.bincount((FILLS[:, None] + FILLS[None, :]).ravel(), minlength=513)
    assert np.array_equal(hist, brute)
    assert int(hist.sum()) == 1 << 18
    seen = np.flatnonzero(hist)
    missing = [v for v in range(513) if hist[v] == 0]
    print(
        f"  histogram total {int(hist.sum())}, distinct fills {seen.size} of 513, mirror symmetric {np.array_equal(hist, hist[::-1])}"
    )
    print(f"  unreachable fills {missing}")
    print(
        f"  fill 0 and 512 count {int(hist[0])} each, peak fill {int(np.argmax(hist))} count {int(hist.max())}"
    )
    print(f"  rules sharing the fill of {rule_name(*LIFE)}: {int(hist[140])}")
    quarters = [int(hist[v]) for v in (64, 128, 140, 256, 384, 448)]
    print(f"  counts at fills 64 128 140 256 384 448: {quarters}")
    return hist


def t3_rules(grid, pins, levels, hist):
    print("  rule | fill | lambda | dimension | deg | genus | count at that fill")
    for birth, survive in NAMED:
        b, s = subset_code(birth), subset_code(survive)
        fill = int(FILLS[b] + FILLS[s])
        genus = "iso" if levels[b, s] else ("axis" if pins[b, s] else "comp")
        dim = f"{log2(fill):.6f}" if fill else "none"
        print(
            f"  {rule_name(birth, survive)} | {fill} | {fill}/512 | {dim} | {int(grid[b, s])} | {genus} | {int(hist[fill])}"
        )


def t3_genus(pins, levels, grid):
    total = 1 << 18
    iso = int(levels.sum())
    axis = int((pins & ~levels).sum())
    comp = total - iso - axis
    print(f"  genus over the 2^18: iso {iso}, axis only {axis}, compound {comp}")
    totalistic = set()
    for chosen in range(1024):
        b = subset_code(k for k in range(9) if (chosen >> k) & 1)
        s = subset_code(k for k in range(9) if (chosen >> (k + 1)) & 1)
        assert levels[b, s]
        totalistic.add((b, s))
    assert len(totalistic) == 1024
    print(
        f"  totalistic rules {len(totalistic)} of the 2^18, every one a level set of the full popcount"
    )
    hist = np.bincount(grid.ravel() + 1, minlength=11)
    print(
        "  degree histogram "
        + " ".join(f"{k - 1}:{int(v)}" for k, v in enumerate(hist) if v)
    )


def t4_affine(grid):
    print("T4 the affine life-like rules")
    sym_deg = degrees(SYM)
    predicted = np.zeros((512, 512), dtype=np.int8)
    for b in range(512):
        diff = b ^ np.arange(512)
        alt = sym_deg[diff]
        predicted[b] = np.where(diff == 0, sym_deg[b], np.maximum(sym_deg[b], 1 + alt))
    assert np.array_equal(predicted, grid)
    print(
        "  deg(B,S) = deg(B) when B = S, else max(deg(B), 1 + deg(B xor S)): holds on all 2^18"
    )
    flat = np.argwhere(grid <= 1)
    names = []
    for b, s in flat:
        birth = [k for k in range(9) if (b >> k) & 1]
        survive = [k for k in range(9) if (s >> k) & 1]
        names.append(rule_name(birth, survive))
    print(f"  degree <= 1 rules: {len(names)}")
    for name in names:
        print(f"    {name}")
    assert len(names) == 8
    expected = {
        "B/S",
        "B012345678/S012345678",
        "B/S012345678",
        "B012345678/S",
        "B1357/S1357",
        "B1357/S02468",
        "B02468/S1357",
        "B02468/S02468",
    }
    assert set(names) == expected
    print("  the four degenerate rules and the four parity rules, no others")
    sym_hist = np.bincount(sym_deg + 1, minlength=10)
    print(
        "  degree histogram of the 512 count sets "
        + " ".join(f"{k - 1}:{int(v)}" for k, v in enumerate(sym_hist) if v)
    )
    for d in range(9):
        assert int((sym_deg <= d).sum()) == 1 << (d + 1), d
    for d in range(1, 9):
        assert int((grid <= d).sum()) == 2 * 4**d, d
    print(
        "  count sets of degree at most d number 2^(d+1), so rules of degree at most d number 2 * 4^d for d = 1..8"
    )


def t5_fredkin():
    print("T5 fredkin = rule 150 tensor rule 150")
    mask = moore(2)
    side = 2 * HORIZON + 3
    mid = side // 2
    for birth, survive, label, centred in (
        ((1, 3, 5, 7), (0, 2, 4, 6, 8), "B1357/S02468", True),
        ((1, 3, 5, 7), (1, 3, 5, 7), "B1357/S1357", False),
    ):
        grid = np.zeros((side, side), dtype=np.uint8)
        grid[mid, mid] = 1
        row = np.zeros(side, dtype=np.uint8)
        row[mid] = 1
        pops = []
        for t in range(HORIZON + 1):
            outer = np.outer(row, row)
            if centred:
                assert np.array_equal(grid, outer), t
            pops.append(int(grid.sum()))
            if t in POWERS or t == 0:
                copies = np.zeros((side, side), dtype=np.uint8)
                for di in (-t, 0, t) if t else (0,):
                    for dj in (-t, 0, t) if t else (0,):
                        copies[mid + di, mid + dj] = 1
                if not centred and t:
                    copies[mid, mid] = 0
                assert np.array_equal(grid, copies), (label, t)
            if not centred:
                shifted = outer.copy()
                shifted[mid, mid] ^= 1
                if t in POWERS:
                    assert np.array_equal(grid, shifted), (label, t)
            grid = step(grid, birth, survive, mask)
            row = step150(row)
        squares = []
        line = np.zeros(side, dtype=np.uint8)
        line[mid] = 1
        for t in range(HORIZON + 1):
            squares.append(int(line.sum()) ** 2)
            line = step150(line)
        if centred:
            assert pops == squares
            print(
                f"  {label} slices equal the outer product cell for cell to t = {HORIZON}"
            )
            print(
                f"  {label} population = rule 150 population squared, first 16: {pops[:16]}"
            )
            print(f"  {label} nine copies at t = {POWERS}")
            assert tuple(pops[:20]) == A246035
        else:
            for t in POWERS:
                assert pops[t] == squares[t] - 1, (label, t)
            print(
                f"  {label} equals the outer product with the centre copy removed at every t = 2^j"
            )
            print(
                f"  {label} eight copies at t = {POWERS}, first 16 populations: {pops[:16]}"
            )
            assert tuple(pops[:20]) == A160239
    line = np.zeros(side, dtype=np.uint8)
    line[mid] = 1
    row150 = []
    for _ in range(20):
        row150.append(int(line.sum()))
        line = step150(line)
    assert tuple(row150) == A071053
    print(f"  rule 150 population, first 20: {row150}")
    print(
        "  20 terms each against OEIS A071053, its square A246035 for B1357/S02468, A160239 for B1357/S1357"
    )


def main():
    start = time.time()
    t1_moore()
    t2_life()
    pins = pin_grid()
    assert np.array_equal(pins, pin_brute())
    assert bool(pins[subset_code(LIFE[0]), subset_code(LIFE[1])]) == is_pin(
        rule_table(*LIFE)
    )
    levels = level_set_grid()
    grid = degree_grid()
    hist = t3_census(pins, levels)
    t3_rules(grid, pins, levels, hist)
    t3_genus(pins, levels, grid)
    t4_affine(grid)
    t5_fredkin()
    print(
        f"domain: every one of the 2^18 life-like rules, all 2^18 axial checks against the closed form, {time.time() - start:.1f} s"
    )


if __name__ == "__main__":
    main()
