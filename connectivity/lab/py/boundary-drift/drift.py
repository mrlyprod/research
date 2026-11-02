import sys
from collections import defaultdict
from fractions import Fraction
from itertools import product

import numpy as np
from scipy import ndimage

DESIGNS = {
    "gasket": (2, 2, [(0, 0), (0, 1), (1, 0)]),
    "diagonal": (2, 2, [(0, 0), (1, 1)]),
    "seven-of-eight": (
        2,
        3,
        [c for c in product(range(2), repeat=3) if c != (1, 1, 1)],
    ),
    "carpet": (3, 2, [c for c in product(range(3), repeat=2) if c != (1, 1)]),
    "sponge": (3, 3, [c for c in product(range(3), repeat=3) if c.count(1) <= 1]),
    "code-79": (3, 2, [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)]),
    "base-4-L": (4, 2, [(0, c) for c in range(4)] + [(r, 0) for r in range(1, 4)]),
}

KRON_LEVELS = {
    "gasket": 8,
    "diagonal": 8,
    "seven-of-eight": 8,
    "carpet": 8,
    "sponge": 5,
    "code-79": 6,
    "base-4-L": 5,
}

DIGIT_LEVELS = {
    "gasket": 8,
    "diagonal": 8,
    "seven-of-eight": 5,
    "carpet": 5,
    "sponge": 4,
    "code-79": 5,
    "base-4-L": 4,
}

NOTE_VALUES = {
    ("gasket", 5): "2.0082",
    ("gasket", 6): "2.0027",
    ("gasket", 7): "2.0009",
    ("gasket", 8): "2.0003",
    ("diagonal", 5): "4.0000",
    ("seven-of-eight", 3): "1.8542",
    ("carpet", 4): "0.8633",
    ("sponge", 3): "2.2560",
}

failures = []


def check(ok, what):
    if not ok:
        failures.append(what)
        print(f"FAIL {what}")


def tile_numbers(base, dim, cells):
    cellset = set(cells)
    pairs = sum(
        1
        for cell in cells
        for axis in range(dim)
        if cell[:axis] + (cell[axis] + 1,) + cell[axis + 1 :] in cellset
    )
    contacts = []
    for axis in range(dim):
        high = {c[:axis] + c[axis + 1 :] for c in cells if c[axis] == base - 1}
        low = {c[:axis] + c[axis + 1 :] for c in cells if c[axis] == 0}
        contacts.append(len(high & low))
    return len(cells), pairs, contacts


def closed_edges(f, pairs, contact, level):
    return Fraction(pairs * (f**level - contact**level), f - contact)


def closed_boundary(dim, f, pairs, contact, level):
    return 2 * dim - 2 * closed_edges(f, pairs, contact, level) / f**level


def kron_grid(base, dim, cells, level):
    tile = np.zeros((base,) * dim, dtype=np.uint8)
    for cell in cells:
        tile[cell] = 1
    out = tile
    for _ in range(level - 1):
        out = np.kron(out, tile)
    return out.astype(bool)


def kron_edges(occupied):
    total = 0
    for axis in range(occupied.ndim):
        low = tuple(
            slice(None, -1) if a == axis else slice(None) for a in range(occupied.ndim)
        )
        high = tuple(
            slice(1, None) if a == axis else slice(None) for a in range(occupied.ndim)
        )
        total += int(np.count_nonzero(occupied[low] & occupied[high]))
    return total


def digit_occupied(point, base, cellset, level):
    for _ in range(level):
        digits = tuple(x % base for x in point)
        if digits not in cellset:
            return False
        point = tuple(x // base for x in point)
    return True


def digit_census(base, dim, cells, level):
    cellset = set(cells)
    side = base**level
    occupied = {
        p
        for p in product(range(side), repeat=dim)
        if digit_occupied(p, base, cellset, level)
    }
    edges, cross = 0, 0
    block = base ** (level - 1)
    for p in occupied:
        for axis in range(dim):
            q = p[:axis] + (p[axis] + 1,) + p[axis + 1 :]
            if q in occupied:
                edges += 1
                if p[axis] // block != q[axis] // block:
                    cross += 1
    return occupied, edges, cross


def designs():
    print("DESIGNS")
    for name, (base, dim, cells) in DESIGNS.items():
        f, pairs, contacts = tile_numbers(base, dim, cells)
        check(len(set(contacts)) == 1, f"{name} contact counts differ by axis")
        contact = contacts[0]
        check(f != contact, f"{name} fill equals contact")
        print(f"- {name}: base {base}, dim {dim}, f {f}, A {pairs}, c {contact}")
        for level in range(1, KRON_LEVELS[name] + 1):
            grid = kron_grid(base, dim, cells, level)
            n, e = int(grid.sum()), kron_edges(grid)
            check(n == f**level, f"{name} {level} cell count")
            want = closed_edges(f, pairs, contact, level)
            check(e == want, f"{name} {level} kron edges {e} against {want}")
            drift = Fraction(2 * dim * n - 2 * e, n)
            check(
                drift == closed_boundary(dim, f, pairs, contact, level),
                f"{name} {level} boundary",
            )
            tag = NOTE_VALUES.get((name, level))
            if tag is not None:
                check(f"{float(drift):.4f}" == tag, f"{name} {level} note value {tag}")
            if level <= DIGIT_LEVELS[name]:
                occupied, de, cross = digit_census(base, dim, cells, level)
                check(len(occupied) == n and de == e, f"{name} {level} digit rule")
                if level >= 2:
                    check(
                        cross == pairs * contact ** (level - 1),
                        f"{name} {level} cross-block pairs {cross}",
                    )
            print(f"  level {level}: N {n} E {e} boundary/cell {float(drift):.6f}")
        limit = Fraction(2 * dim) - Fraction(2 * pairs, f - contact)
        print(f"  limit {limit} = {float(limit):.6f}, ratio {Fraction(contact, f)}")


def gasket_tree():
    print("GASKET")
    base, dim, cells = DESIGNS["gasket"]
    for level in range(1, 9):
        grid = kron_grid(base, dim, cells, level)
        n, e = int(grid.sum()), kron_edges(grid)
        _, pieces = ndimage.label(grid)
        check(pieces == 1 and e == n - 1, f"gasket {level} tree")
        check(4 * n - 2 * e == 2 * n + 2, f"gasket {level} perimeter 2N+2")
        print(f"  level {level}: N {n} E {e} pieces {pieces} perimeter {4 * n - 2 * e}")
    base, dim, cells = DESIGNS["carpet"]
    for level in range(1, 5):
        grid = kron_grid(base, dim, cells, level)
        n, e = int(grid.sum()), kron_edges(grid)
        _, pieces = ndimage.label(grid)
        check(pieces == 1 and e >= n - 1, f"carpet {level} connected")
        check(4 * n - 2 * e < 2 * n + 2, f"carpet {level} below the tree bound")
        print(
            f"  carpet level {level}: N {n} E {e} pieces {pieces} "
            f"perimeter {4 * n - 2 * e}"
        )


def tree_codes(base, levels):
    found = []
    for code in range(1, 2 ** (base * base)):
        cells = [(i // base, i % base) for i in range(base * base) if code >> i & 1]
        f = len(cells)
        ok = f >= 2
        rows = []
        for level in range(1, levels + 1) if ok else ():
            grid = kron_grid(base, 2, cells, level)
            n, e = int(grid.sum()), kron_edges(grid)
            _, pieces = ndimage.label(grid)
            if not (pieces == 1 and e == n - 1):
                ok = False
                break
            rows.append(Fraction(4 * n - 2 * e, n) == 2 + Fraction(2, f**level))
        if ok:
            check(all(rows), f"base {base} code {code} drift")
            found.append((code, f))
    return found


def census():
    print("TREE CODES")
    for base, levels in ((2, 4), (3, 3)):
        found = tree_codes(base, levels)
        fills = sorted({f for _, f in found})
        print(f"- base {base}: {len(found)} tree codes, fills {fills}")
        check(len(found) == {2: 8, 3: 15}[base], f"base {base} tree count")


def code_cells(code):
    return [(i // 2, i % 2) for i in range(4) if code >> i & 1]


def code_axes(code):
    cells = code_cells(code)
    f, _, contacts = tile_numbers(2, 2, cells)
    cellset = set(cells)
    per_axis = [
        sum(1 for c in cells if c[:axis] + (c[axis] + 1,) + c[axis + 1 :] in cellset)
        for axis in range(2)
    ]
    return f, per_axis, contacts


def word_grid(word):
    out = np.ones((1, 1), dtype=bool)
    for code in word:
        tile = np.zeros((2, 2), dtype=bool)
        for cell in code_cells(code):
            tile[cell] = True
        out = np.kron(out, tile)
    return out


def word_formula(word, letters):
    total = 0
    for axis in range(2):
        for i in range(len(word)):
            term = letters[word[i]][1][axis]
            for j in range(i):
                term *= letters[word[j]][0]
            for j in range(i + 1, len(word)):
                term *= letters[word[j]][2][axis]
            total += term
    return total


def words():
    print("WORDS")
    letters = {code: code_axes(code) for code in range(1, 16)}
    seen = 0
    for length in range(1, 5):
        count = 0
        for word in product(range(1, 16), repeat=length):
            count += 1
            e = kron_edges(word_grid(word))
            if e != word_formula(word, letters):
                check(False, f"word {word} edges")
        seen += count
        print(f"  length {length}: {count} words")
    print(f"- {seen} words of length 1 to 4 over the 15 codes")
    library = [c for c in range(1, 16) if c.bit_count() in (2, 3)]
    for alphabet, length, name in (
        (range(1, 16), 2, "15 codes"),
        (library, 2, "10-code library"),
        (library, 3, "10-code library"),
    ):
        groups = defaultdict(set)
        for word in product(alphabet, repeat=length):
            groups[tuple(sorted(word))].add(kron_edges(word_grid(word)))
        many = [k for k in groups if len(set(k)) > 1]
        differ = sum(1 for k in many if len(groups[k]) > 1)
        print(
            f"- length {length}, {name}: {differ} of {len(many)} multisets differ in E"
        )
        if (length, name) == (2, "15 codes"):
            check((differ, len(many)) == (78, 105), "length 2 sensitivity")
        if (length, name) == (3, "10-code library"):
            check((differ, len(many)) == (188, 210), "length 3 sensitivity")


def sections():
    print("DOMAIN")
    print(
        "- designs "
        + ", ".join(DESIGNS)
        + "; kron levels "
        + ", ".join(f"{n} 1..{v}" for n, v in KRON_LEVELS.items())
    )
    print(
        "- digit-rule levels "
        + ", ".join(f"{n} 1..{v}" for n, v in DIGIT_LEVELS.items())
    )
    designs()
    gasket_tree()
    census()
    words()
    print("FAILURES" if failures else "ok: every check passed")


if __name__ == "__main__":
    sections()
    sys.exit(1 if failures else 0)
