import math
import os
import zlib

import numpy as np
from scipy.signal import correlate

SIDE = 243
LEVEL = 5
SEED_LIFE = 1729
SEED_RANDOM = 8128
SEED_NOISE = 496
SEED_CHECK = 6561
LIFE_SIDE = 16
LIFE_GENS = 6
LIFE_SEEDS = 8
LIFE_DENSITY = 0.3
CAP = 6000000
CARPET2 = 7
CARPET3 = 495
LETTERS = ((7, "carpet"), (14, "net"), (3, "htree"), (5, "vtree"), (9, "void"))
MAGIC_SIDES = ((3, 5), (5, 3), (3, 7), (7, 3))
RULES = (
    ((3,), (2, 3)),
    ((3, 6), (2, 3)),
    ((1, 3, 5, 7), (1, 3, 5, 7)),
    ((1, 3, 5, 7), (0, 2, 4, 6, 8)),
    ((0, 2, 4, 6, 8), (1, 3, 5, 7)),
    ((0, 2, 4, 6, 8), (0, 2, 4, 6, 8)),
    ((2,), ()),
    ((3, 4), (3, 4)),
)
BAYER = np.array([[0, 8, 2, 10], [12, 4, 14, 6], [3, 11, 1, 9], [15, 7, 13, 5]])
CLASSES = ("kronecker", "magic", "life")
HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))

# TILES


def tile(code, base, side):
    i = np.arange(side) % base
    bit = base * i[:, None] + i[None, :]
    return ((code >> bit) & 1).astype(np.uint8)


def kron(a, b):
    return np.kron(a, b).astype(np.uint8)


def power(a, level):
    out = a
    for _ in range(level - 1):
        out = kron(out, a)
    return out


def square_orbit(code):
    g = tile(code, 3, 3)
    seen = set()
    for flip in (g, g[:, ::-1]):
        h = flip
        for _ in range(4):
            h = np.rot90(h)
            seen.add(
                int(
                    sum(int(h[r, c]) << (3 * r + c) for r in range(3) for c in range(3))
                )
            )
    return seen


def base3_reps():
    seen = set()
    reps = []
    for code in range(1, 512):
        if code in seen:
            continue
        orbit = square_orbit(code)
        seen |= orbit
        reps.append(min(orbit))
    return sorted(reps)


def magic_word(codes, sides):
    out = tile(codes[0], 2, sides[0])
    for code, side in zip(codes[1:], sides[1:]):
        out = kron(out, tile(code, 2, side))
    return out


# LIFE


def life_step(grid, birth, survive):
    n = np.zeros(grid.shape, dtype=np.int8)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr or dc:
                n += np.roll(np.roll(grid, dr, axis=0), dc, axis=1)
    born = np.isin(n, birth) & (grid == 0)
    kept = np.isin(n, survive) & (grid == 1)
    return (born | kept).astype(np.uint8)


def rule_name(birth, survive):
    return "rule birth [%s], survive [%s], wrap" % (
        " ".join(str(d) for d in birth),
        " ".join(str(d) for d in survive),
    )


def life_seeds():
    rng = np.random.default_rng(SEED_LIFE)
    out = [np.zeros((LIFE_SIDE, LIFE_SIDE), dtype=np.uint8)]
    out[0][LIFE_SIDE // 2, LIFE_SIDE // 2] = 1
    for _ in range(LIFE_SEEDS - 1):
        out.append((rng.random((LIFE_SIDE, LIFE_SIDE)) < LIFE_DENSITY).astype(np.uint8))
    return out


# THE CATALOG


def build_catalog(depth):
    reps = base3_reps()
    raw = []
    for code in range(16):
        raw.append(
            ("kronecker", "bang dim 2, code %d, level 1" % code, tile(code, 2, 2))
        )
    for code in range(16):
        raw.append(
            (
                "kronecker",
                "bang dim 2, code %d, level 2" % code,
                power(tile(code, 2, 2), 2),
            )
        )
    for level in range(1, depth + 1):
        for code in reps:
            raw.append(
                (
                    "kronecker",
                    "bang dim 2, base 3, code %d, level %d" % (code, level),
                    power(tile(code, 3, 3), level),
                )
            )
    for sides in MAGIC_SIDES:
        for ca, na in LETTERS:
            for cb, nb in LETTERS:
                raw.append(
                    (
                        "magic",
                        "%s(%d), %s(%d)" % (na, sides[0], nb, sides[1]),
                        magic_word((ca, cb), sides),
                    )
                )
    for birth, survive in RULES:
        name = rule_name(birth, survive)
        for si, seed in enumerate(life_seeds()):
            grid = seed
            for gen in range(1, LIFE_GENS + 1):
                grid = life_step(grid, birth, survive)
                raw.append(("life", "%s s%d g%d" % (name, si, gen), grid.copy()))
    atoms = []
    index = {}
    collisions = 0
    for cls, label, arr in raw:
        if not arr.any():
            continue
        key = (arr.shape, arr.tobytes())
        if key in index:
            collisions += 1
            continue
        index[key] = len(atoms)
        atoms.append(
            {
                "cls": cls,
                "label": label,
                "arr": arr,
                "h": int(arr.shape[0]),
                "w": int(arr.shape[1]),
                "ones": int(arr.sum()),
            }
        )
    return atoms, len(raw), collisions, len(reps)


# CORPORA


def corpus_tree():
    return power(tile(CARPET2, 2, 3), LEVEL)


def corpus_text():
    need = (SIDE * SIDE + 7) // 8
    with open(os.path.join(TREE, "README.md"), "rb") as fh:
        data = fh.read()
    assert len(data) >= need
    bits = np.unpackbits(np.frombuffer(data[:need], dtype=np.uint8))
    return bits[: SIDE * SIDE].reshape(SIDE, SIDE).copy()


def corpus_random():
    return (np.random.default_rng(SEED_RANDOM).random((SIDE, SIDE)) < 0.5).astype(
        np.uint8
    )


def corpus_halftone():
    axis = (np.arange(SIDE) - (SIDE - 1) / 2) / ((SIDE - 1) / 2)
    grey = np.clip(1.0 - np.sqrt(axis[:, None] ** 2 + axis[None, :] ** 2), 0.0, 1.0)
    thresh = (
        BAYER[np.arange(SIDE)[:, None] % 4, np.arange(SIDE)[None, :] % 4] + 0.5
    ) / 16.0
    return (grey > thresh).astype(np.uint8)


# THE SCORE

TABLES = {}


def log2c_table(n):
    if n not in TABLES:
        if n <= 512:
            TABLES[n] = np.array([math.log2(math.comb(n, d)) for d in range(n + 1)])
        else:
            i = np.arange(1, n + 1, dtype=np.float64)
            TABLES[n] = np.concatenate(([0.0], np.cumsum(np.log2((n - i + 1) / i))))
    return TABLES[n]


def log2c(n, k):
    if n <= 512:
        return math.log2(math.comb(n, k))
    if k == 0 or k == n:
        return 0.0
    return (
        math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    ) / math.log(2.0)


def window_sums(x, h, w):
    ii = np.zeros((x.shape[0] + 1, x.shape[1] + 1), dtype=np.int64)
    ii[1:, 1:] = np.cumsum(np.cumsum(x.astype(np.int64), axis=0), axis=1)
    return ii[h:, w:] - ii[:-h, w:] - ii[h:, :-w] + ii[:-h, :-w]


def mismatch(x, atom):
    corr = correlate(
        x.astype(np.float64), atom["arr"].astype(np.float64), mode="valid", method="fft"
    )
    assert np.max(np.abs(corr - np.rint(corr))) < 1e-6
    return (
        window_sums(x, atom["h"], atom["w"])
        + atom["ones"]
        - 2 * np.rint(corr).astype(np.int64)
    )


def check_mismatch(x, atoms):
    rng = np.random.default_rng(SEED_CHECK)
    bad = 0
    for _ in range(200):
        atom = atoms[int(rng.integers(len(atoms)))]
        h, w = atom["h"], atom["w"]
        if h > x.shape[0] or w > x.shape[1]:
            continue
        r = int(rng.integers(x.shape[0] - h + 1))
        c = int(rng.integers(x.shape[1] - w + 1))
        direct = int((x[r : r + h, c : c + w] ^ atom["arr"]).sum())
        bad += direct != int(mismatch(x, atom)[r, c])
    return bad


def floor_cells(log2a, positions):
    n = 1
    while n - (log2a + math.log2(positions) + math.log2(n + 1)) <= 0:
        n += 1
    return n


# THE PURSUIT


def candidates(x, atoms, tables):
    height, width = x.shape
    log2a = math.log2(len(atoms))
    sav_all, who_all, pos_all, cell_all = [], [], [], []
    live = 0
    for ai, atom in enumerate(atoms):
        h, w = atom["h"], atom["w"]
        if h > height or w > width:
            continue
        nw = h * w
        cols = width - w + 1
        positions = (height - h + 1) * cols
        base = log2a + math.log2(positions) + math.log2(nw + 1)
        if nw - base <= 0:
            continue
        live += 1
        sav = (nw - base) - tables[nw][mismatch(x, atom)]
        hit = np.nonzero(sav > 0)
        if hit[0].size == 0:
            continue
        sav_all.append(sav[hit].astype(np.float32))
        who_all.append(np.full(hit[0].size, ai, dtype=np.int32))
        pos_all.append((hit[0] * cols + hit[1]).astype(np.int32))
        cell_all.append(np.full(hit[0].size, nw, dtype=np.int32))
    if not sav_all:
        return (
            np.zeros(0, np.float32),
            np.zeros(0, np.int32),
            np.zeros(0, np.int32),
            np.zeros(0, np.int32),
            live,
        )
    return (
        np.concatenate(sav_all),
        np.concatenate(who_all),
        np.concatenate(pos_all),
        np.concatenate(cell_all),
        live,
    )


def sweep(x, atoms, sav, who, pos, cell, key):
    height, width = x.shape
    score = sav if key == "flat" else sav / cell
    order = np.argsort(-score, kind="stable")
    capped = order.size > CAP
    if capped:
        order = order[:CAP]
    covered = np.zeros((height, width), dtype=bool)
    flat = bytearray(height * width)
    out = []
    for k in order:
        ai = int(who[k])
        atom = atoms[ai]
        h, w = atom["h"], atom["w"]
        cols = width - w + 1
        r, c = divmod(int(pos[k]), cols)
        if flat[r * width + c]:
            continue
        if covered[r : r + h, c : c + w].any():
            continue
        covered[r : r + h, c : c + w] = True
        for rr in range(r, r + h):
            flat[rr * width + c : rr * width + c + w] = b"\x01" * w
        out.append((ai, r, c))
    return out, capped


def describe(x, atoms, placements):
    height, width = x.shape
    n = height * width
    log2a = math.log2(len(atoms))
    covered = np.zeros((height, width), dtype=bool)
    rebuilt = np.zeros((height, width), dtype=np.uint8)
    per_class = dict.fromkeys(CLASSES, 0.0)
    used = {c: {} for c in CLASSES}
    cost = 0.0
    for ai, r, c in placements:
        atom = atoms[ai]
        h, w = atom["h"], atom["w"]
        nw = h * w
        mask = x[r : r + h, c : c + w] ^ atom["arr"]
        d = int(mask.sum())
        this = (
            log2a
            + math.log2((height - h + 1) * (width - w + 1))
            + math.log2(nw + 1)
            + log2c(nw, d)
        )
        cost += this
        per_class[atom["cls"]] += nw - this
        row = used[atom["cls"]].setdefault(ai, [0, 0.0])
        row[0] += 1
        row[1] += nw - this
        covered[r : r + h, c : c + w] = True
        rebuilt[r : r + h, c : c + w] = atom["arr"] ^ mask
    rebuilt[~covered] = x[~covered]
    assert np.array_equal(rebuilt, x)
    u = int((~covered).sum())
    u1 = int(x[~covered].sum())
    enum = math.log2(u + 1) + log2c(u, u1) if u else 0.0
    residual = 1.0 + min(float(u), enum)
    header = math.log2(n + 1)
    return {
        "bits": header + cost + residual,
        "header": header,
        "residual": residual,
        "uncovered": u,
        "credit": float(u) - residual,
        "per_class": per_class,
        "used": used,
        "placements": len(placements),
    }


def deflate_bits(x):
    co = zlib.compressobj(9, zlib.DEFLATED, -15)
    return 8 * len(co.compress(np.packbits(x.ravel()).tobytes()) + co.flush())


def run(x, atoms, tables, key="dense"):
    sav, who, pos, cell, live = candidates(x, atoms, tables)
    placements, capped = sweep(x, atoms, sav, who, pos, cell, key)
    full = describe(x, atoms, placements)
    bare = describe(x, atoms, [])
    assert (
        abs(
            sum(full["per_class"].values())
            - full["header"]
            + full["credit"]
            - (x.size - full["bits"])
        )
        < 1e-6
    )
    won = full["bits"] < bare["bits"]
    rep = dict(full if won else bare)
    rep["bits"] += 1.0
    rep.update(
        {
            "mode": "atoms" if won else "bare",
            "pursued": len(placements),
            "pursuit": full["bits"] + 1.0,
            "bare": bare["bits"] + 1.0,
            "live": live,
            "cands": int(sav.size),
            "capped": capped,
            "raw": x.size,
            "zlib": deflate_bits(x),
            "ones": int(x.sum()),
        }
    )
    return rep


# THE REARRANGEMENT


def rearrange(a, p, q):
    m, n = a.shape
    assert m % p == 0 and n % q == 0
    return (
        a.reshape(p, m // p, q, n // q)
        .transpose(0, 2, 1, 3)
        .reshape(p * q, (m // p) * (n // q))
    )


def kron_spectrum(a, p, q):
    return np.linalg.svd(rearrange(a.astype(np.float64), p, q), compute_uv=False)


def nearest_factors(a, p, q):
    u, s, vt = np.linalg.svd(rearrange(a.astype(np.float64), p, q), full_matrices=False)
    outer, inner = s[0] * u[:, 0], s[0] * vt[0]
    if outer[np.argmax(np.abs(outer))] < 0:
        outer, inner = -outer, -inner
    return outer.reshape(p, q), inner.reshape(a.shape[0] // p, a.shape[1] // q), s


def block_means(a, p, q):
    m, n = a.shape
    return a.reshape(p, m // p, q, n // q).mean(axis=(1, 3))


def threshold(m):
    return (m > 0.5 * m.max()).astype(np.uint8)


def code_of(t):
    side = t.shape[0]
    return int(
        sum(int(t[r, c]) << (side * r + c) for r in range(side) for c in range(side))
    )


def peel(a, sides):
    out, rest = [], a
    for side in sides:
        outer, inner, s = nearest_factors(rest, side, side)
        out.append((code_of(threshold(outer)), s[1] / s[0]))
        rest = threshold(inner)
    return out, rest


# THE RUN


def table(rows):
    print(
        "%-12s %8s %8s %8s %9s %9s %9s"
        % ("corpus", "ones", "raw", "zlib", "codebook", "vs raw", "vs zlib")
    )
    for name, r in rows:
        bits = math.ceil(r["bits"])
        print(
            "%-12s %8d %8d %8d %9d %9d %9d"
            % (
                name,
                r["ones"],
                r["raw"],
                r["zlib"],
                bits,
                r["raw"] - bits,
                r["zlib"] - bits,
            )
        )


def classes(rows):
    print(
        "%-12s %6s %5s %5s %5s %10s %10s %10s %9s"
        % (
            "corpus",
            "place",
            "kron",
            "magic",
            "life",
            "bits K",
            "bits M",
            "bits L",
            "residual",
        )
    )
    for name, r in rows:
        u = r["used"]
        print(
            "%-12s %6d %5d %5d %5d %10.0f %10.0f %10.0f %9.0f"
            % (
                name,
                r["placements"],
                len(u["kronecker"]),
                len(u["magic"]),
                len(u["life"]),
                r["per_class"]["kronecker"],
                r["per_class"]["magic"],
                r["per_class"]["life"],
                r["credit"],
            )
        )


def top_atoms(atoms, rep, k):
    rows = []
    for cls in CLASSES:
        for ai, (n, bits) in rep["used"][cls].items():
            rows.append((bits, n, cls, atoms[ai]))
    rows.sort(key=lambda t: -t[0])
    for bits, n, cls, atom in rows[:k]:
        print(
            "    %-8s %-34s %3d x %-3d fill %5d  %4d placements  %8.0f bits"
            % (cls, atom["label"], atom["h"], atom["w"], atom["ones"], n, bits)
        )


def main():
    atoms, built, collisions, reps = build_catalog(2)
    assert reps == 101
    log2a = math.log2(len(atoms))
    sizes = sorted({a["h"] * a["w"] for a in atoms})
    tables = {s: log2c_table(s) for s in sizes}
    counts = {}
    for a in atoms:
        counts[(a["cls"], a["h"])] = counts.get((a["cls"], a["h"]), 0) + 1

    print("== THE CATALOG ==")
    print(
        "%d atoms built, %d duplicates dropped, %d empty frames dropped, %d kept"
        % (built, collisions, built - collisions - len(atoms), len(atoms))
    )
    for key in sorted(counts):
        print("  %-9s %3d x %-3d %4d atoms" % (key[0], key[1], key[1], counts[key]))
    print(
        "the base-3 class is the %d nonempty square-group orbits of the 511 nonempty plane codes"
        % reps
    )
    print(
        "naming one atom costs log2(%d) = %.4f bits, naming a position at most log2(%d) = %.4f"
        % (len(atoms), log2a, SIDE * SIDE, math.log2(SIDE * SIDE))
    )
    fl = floor_cells(log2a, SIDE * SIDE)
    dead = sum(1 for a in atoms if a["h"] * a["w"] < fl)
    print(
        "an atom of n cells pays only when n > log2(%d) + log2(%d) + log2(n+1), so no atom below %d cells can ever be placed"
        % (len(atoms), SIDE * SIDE, fl)
    )
    print(
        "%d of the %d atoms sit below that floor: they can never pay and still charge every other atom their share of the name"
        % (dead, len(atoms))
    )
    print()

    corpora = [
        (
            "tree render",
            "bang dim 2, code 7 at level %d, side %d" % (LEVEL, SIDE),
            corpus_tree(),
        ),
        (
            "text",
            "the first %d bytes of research/README.md, unpacked MSB first"
            % ((SIDE * SIDE + 7) // 8),
            corpus_text(),
        ),
        (
            "random",
            "uniform bits from numpy default_rng(%d)" % SEED_RANDOM,
            corpus_random(),
        ),
        (
            "halftone",
            "the radial gradient 1 - r ordered-dithered by the 4 x 4 Bayer matrix",
            corpus_halftone(),
        ),
    ]
    bad = sum(check_mismatch(x, atoms) for _, _, x in corpora)
    print(
        "the correlation mismatch count agrees with a direct window comparison on 800 sampled placements, %d failures"
        % bad
    )
    assert bad == 0
    print()

    rows = [(name, run(x, atoms, tables)) for name, _, x in corpora]
    print("== THE TABLE ==")
    table(rows)
    print()
    classes(rows)
    print()
    for (name, source, x), (_, r) in zip(corpora, rows):
        print("%-12s %s" % (name, source))
        print(
            "    %d of %d cells covered, %d uncovered, reconstruction exact; %d atoms could be placed, %d candidate placements scored positive, cap hit %s"
            % (
                r["raw"] - r["uncovered"],
                r["raw"],
                r["uncovered"],
                r["live"],
                r["cands"],
                r["capped"],
            )
        )
        print(
            "    the pursuit takes %d placements and describes the corpus in %.0f bits, the placement-free description takes %.0f, the encoder emits the %s mode"
            % (r["pursued"], r["pursuit"], r["bare"], r["mode"])
        )
        print(
            "    saving = %.0f atom bits - %.1f header + %.0f residual credit = %.0f"
            % (
                sum(r["per_class"].values()),
                r["header"],
                r["credit"],
                r["raw"] - r["bits"],
            )
        )
        for cls in CLASSES:
            if r["used"][cls]:
                fills = [
                    atoms[ai]["ones"] / (atoms[ai]["h"] * atoms[ai]["w"])
                    for ai in r["used"][cls]
                ]
                print(
                    "    %-9s %2d distinct atoms placed, their densities running %.4f to %.4f"
                    % (cls, len(fills), min(fills), max(fills))
                )
        top_atoms(atoms, r, 4)
    print()

    print("== THE CONTROL ==")
    rnd = dict(rows)["random"]
    print(
        "on uniform bits the codebook spends %.0f bits against a raw %d and a zlib %d, a saving of %.0f, with %d placements"
        % (
            rnd["bits"],
            rnd["raw"],
            rnd["zlib"],
            rnd["raw"] - rnd["bits"],
            rnd["pursued"],
        )
    )
    print(
        "it reads under zlib there only because deflate expands an incompressible stream by %d bits; neither compresses"
        % (rnd["zlib"] - rnd["raw"])
    )
    assert rnd["pursued"] == 0 and rnd["raw"] - rnd["bits"] <= 0
    txt = dict(rows)["text"]
    print(
        "on text the codebook places nothing either, so its %.0f bit saving is the residual entropy code and not an atom"
        % (txt["raw"] - txt["bits"])
    )
    assert txt["pursued"] == 0
    print()

    print("== THE GREEDY KEY ==")
    x = corpus_tree()
    sav, who, pos, cell, _ = candidates(x, atoms, tables)
    for key in ("dense", "flat"):
        pl, _ = sweep(x, atoms, sav, who, pos, cell, key)
        rep = describe(x, atoms, pl)
        print(
            "  key %-5s: %4d placements, codebook %6.0f bits, saving %6.0f, K/M/L bits %6.0f %6.0f %6.0f"
            % (
                key,
                rep["placements"],
                rep["bits"] + 1.0,
                x.size - rep["bits"] - 1.0,
                rep["per_class"]["kronecker"],
                rep["per_class"]["magic"],
                rep["per_class"]["life"],
            )
        )
    print(
        "the table above uses the dense key, saving per cell; the flat key, saving per placement, spends the plane on large loose atoms"
    )
    print()

    print("== THE DEPTH LADDER ==")
    print(
        "%-7s %7s %7s %9s %9s %9s %7s %10s %6s %8s %6s %8s"
        % (
            "depth",
            "atoms",
            "log2 A",
            "codebook",
            "vs raw",
            "vs zlib",
            "place",
            "kron bits",
            "r pl",
            "r save",
            "t pl",
            "t save",
        )
    )
    for depth in (2, 3, 4, 5):
        cat, _, _, _ = build_catalog(depth)
        tab = {s: log2c_table(s) for s in sorted({a["h"] * a["w"] for a in cat})}
        rep = run(corpus_tree(), cat, tab)
        rnd2 = run(corpus_random(), cat, tab)
        txt2 = run(corpus_text(), cat, tab)
        assert rnd2["pursued"] == 0 and rnd2["raw"] - rnd2["bits"] <= 0
        bits = math.ceil(rep["bits"])
        print(
            "%-7d %7d %7.4f %9d %9d %9d %7d %10.0f %6d %8.0f %6d %8.0f"
            % (
                depth,
                len(cat),
                math.log2(len(cat)),
                bits,
                rep["raw"] - bits,
                rep["zlib"] - bits,
                rep["placements"],
                rep["per_class"]["kronecker"],
                rnd2["pursued"],
                rnd2["raw"] - rnd2["bits"],
                txt2["pursued"],
                txt2["raw"] - txt2["bits"],
            )
        )
    print(
        "the random control takes zero placements at every depth, so a deeper catalog never manufactures a saving where there is none"
    )
    print(
        "text takes a handful from depth 3 on, each worth under two bits and each lowering the corpus total, so the encoder drops them and emits the placement-free mode"
    )
    print(
        "the level-%d tile of the tree corpus is itself an atom at depth %d, so the last row is recognition and not compression"
        % (LEVEL, LEVEL)
    )
    print()

    print("== THE KRONECKER SPECTRUM ==")
    carpet = corpus_tree()
    for j in (1, 2, 3, 4):
        side = 3**j
        outer, inner, s = nearest_factors(carpet, side, side)
        assert np.array_equal(threshold(outer), power(tile(CARPET2, 2, 3), j))
        assert np.array_equal(threshold(inner), power(tile(CARPET2, 2, 3), LEVEL - j))
        print(
            "  split %3d x %-3d sigma_1 = %.4f, sigma_2/sigma_1 = %.3e; the rank-one factors are carpet levels %d and %d exactly"
            % (side, side, s[0], s[1] / s[0], j, LEVEL - j)
        )
    got = code_of(threshold(nearest_factors(carpet, 3, 3)[0]))
    assert got == CARPET3
    print(
        "  at the 3 x 3 split the recovered base-3 code is %d, the carpet's own level-1 code"
        % got
    )
    print()

    word = (7, 9, 14, 3, 5)
    names = dict(LETTERS)
    read, rest = peel(magic_word(word, (3,) * 5), (3,) * 4)
    got = [c for c, _ in read] + [code_of(rest)]
    want = [code_of(tile(c, 2, 3)) for c in word]
    assert got == want
    print("== THE PEEL ==")
    print(
        "  the magic word %s renders at side %d and the rearranged SVD peels it letter by letter"
        % (", ".join("%s(3)" % names[c] for c in word), SIDE)
    )
    print(
        "  recovered base-3 codes %s, the word's own %s, ratios sigma_2/sigma_1 %s"
        % (got, want, ", ".join("%.1e" % r for _, r in read))
    )
    print()

    print("== THE NOISE DIAL ==")
    rng = np.random.default_rng(SEED_NOISE)
    trials = 40
    svd, mean = [], []
    for k in range(51):
        p = k / 100.0
        a = b = 0
        for _ in range(trials):
            noisy = (carpet ^ (rng.random(carpet.shape) < p)).astype(np.uint8)
            a += code_of(threshold(nearest_factors(noisy, 3, 3)[0])) == CARPET3
            b += code_of(threshold(block_means(noisy, 3, 3))) == CARPET3
        svd.append((p, a))
        mean.append((p, b))
    print(
        "  flipping a fraction p of the %d bits, %d seeds per p on the grid 0.00 to 0.50 in steps of 0.01"
        % (SIDE * SIDE, trials)
    )
    for name, curve in (("rank one", svd), ("block mean", mean)):
        clean = max(
            p for p, ok in curve if all(o == trials for q, o in curve if q <= p)
        )
        first = min((p for p, ok in curve if ok < trials), default=None)
        half = min((p for p, ok in curve if ok * 2 < trials), default=None)
        none = min((p for p, ok in curve if ok == 0), default=None)
        print(
            "  %-11s every seed to p = %.2f, first failure p = %.2f, half the seeds by p = %.2f, none from p = %.2f"
            % (name, clean, first, half, none)
        )
        print(
            "    %s"
            % "  ".join("%.2f:%d" % (p, ok) for p, ok in curve if 0.20 <= p <= 0.35)
        )
    fill = (8.0 / 9.0) ** (LEVEL - 1)
    star = fill / (1.0 + 2.0 * fill)
    half_mean = min(p for p, ok in mean if ok * 2 < trials)
    print(
        "  the block-mean detector's expected margin vanishes at f/(1 + 2f) = %.6f with f = (8/9)^%d = %.6f the density of a filled corner block,"
        % (star, LEVEL - 1, fill)
    )
    print(
        "  a closed form read off the fill law and not off this sweep; the sweep keeps more than half its seeds at p = %.2f and loses them by p = %.2f,"
        % (half_mean - 0.01, half_mean)
    )
    print(
        "  one grid step above the closed form because the threshold divides the largest of eight noisy filled corners rather than their mean"
    )
    assert abs(half_mean - star) <= 0.02
    print(
        "  the rank-one factor of the rearrangement holds %.2f further and falls off a cliff two grid steps wide"
        % (min(p for p, ok in svd if ok * 2 < trials) - half_mean)
    )
    for p in (0.00, 0.10, 0.20, 0.30, 0.40, 0.50):
        noisy = (carpet ^ (rng.random(carpet.shape) < p)).astype(np.uint8)
        s = kron_spectrum(noisy, 3, 3)
        print(
            "    p = %.2f: sigma_2/sigma_1 = %.4f, recovered code %d against %d"
            % (
                p,
                s[1] / s[0],
                code_of(threshold(nearest_factors(noisy, 3, 3)[0])),
                CARPET3,
            )
        )
    print(
        "  the spectrum itself degrades smoothly and reads nothing, while the thresholded rank-one factor stays exact far past it"
    )
    print()
    print("every assertion passed")


main()
