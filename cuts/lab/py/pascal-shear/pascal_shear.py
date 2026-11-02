import sys
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
GOULD_URL = "https://oeis.org/A001316/b001316.txt"
GASKET_URL = "https://oeis.org/A047999/b047999.txt"
GASKET_CACHE = "a047999.txt"
CODE = 7
TOP_LEVEL = 9
DIGIT_LEVEL = 14
EXACT = 128
ROWS = 1024
SHEAR_LEVEL = 6
TRUNCATION_LEVEL = 8


def popcount(value):
    return value.bit_count()


def tile_of(code):
    grid = [[0, 0], [0, 0]]
    for index, (x, y) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
        grid[x][y] = (code >> index) & 1
    return grid


def kron(a, b):
    ra, ca = len(a), len(a[0])
    rb, cb = len(b), len(b[0])
    out = [[0] * (ca * cb) for _ in range(ra * rb)]
    for i in range(ra):
        for j in range(ca):
            if a[i][j]:
                for u in range(rb):
                    for v in range(cb):
                        out[i * rb + u][j * cb + v] = b[u][v]
    return out


def fractal(tile, level):
    out = [row[:] for row in tile]
    for _ in range(1, level):
        out = kron(out, tile)
    return out


def pascal_mod2(rows):
    out = [[1]]
    for n in range(1, rows):
        previous = out[n - 1]
        row = [1]
        for k in range(1, n):
            row.append((previous[k - 1] + previous[k]) & 1)
        row.append(1)
        out.append(row)
    return out


def carry_count(i, j):
    count = 0
    carry = 0
    while i or j or carry:
        total = (i & 1) + (j & 1) + carry
        carry = 1 if total >= 2 else 0
        count += carry
        i >>= 1
        j >>= 1
    return count


def binomial(n, k):
    k = min(k, n - k)
    out = 1
    for i in range(k):
        out = out * (n - i) // (i + 1)
    return out


def valuation2(value):
    return (value & -value).bit_length() - 1


def parse_terms(text):
    terms = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        terms.append((int(parts[0]), int(parts[1])))
    return terms


def read_cache(name):
    return parse_terms((HERE / name).read_text())


def read_live(url):
    request = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    try:
        with urllib.request.urlopen(request, timeout=30) as handle:
            return parse_terms(handle.read().decode())
    except (urllib.error.URLError, OSError, TimeoutError):
        return None


def triangular_rows(count):
    rows = 0
    seen = 0
    while seen < count:
        rows += 1
        seen += rows
    return rows


def level_sets():
    print("KRONECKER LEVEL SET OF CODE 7")
    tile = tile_of(CODE)
    print(f"  tile [[{tile[0][0]}, {tile[0][1]}], [{tile[1][0]}, {tile[1][1]}]]")
    for level in range(1, TOP_LEVEL + 1):
        grid = fractal(tile, level)
        side = 1 << level
        cells = {(i, j) for i in range(side) for j in range(side) if grid[i][j]}
        anded = {(i, j) for i in range(side) for j in range(side) if i & j == 0}
        agree = cells == anded and len(cells) == 3**level
        print(
            f"  L={level:2d}  side {side:4d}  cells {len(cells):6d}"
            f"  3^L {3**level:6d}  set == (i AND j == 0) {agree}"
        )
    same = all(
        sum(1 << (level - popcount(i)) for i in range(1 << level)) == 3**level
        for level in range(1, DIGIT_LEVEL + 1)
    )
    print(f"  digit sum gives 3^L for L = 1..{DIGIT_LEVEL}: {same}")


def kummer():
    print("KUMMER ON EXACT BINOMIALS")
    faults = 0
    for i in range(EXACT):
        for j in range(EXACT):
            value = binomial(i + j, i)
            if valuation2(value) != carry_count(i, j):
                faults += 1
            if (value & 1 == 1) != (i & j == 0):
                faults += 1
    print(f"  0 <= i, j < {EXACT}: {EXACT * EXACT} binomials, {faults} faults")


def recurrence(triangle):
    print("PASCAL MOD 2 FROM THE ADDITIVE RECURRENCE")
    faults = 0
    entries = 0
    for n, row in enumerate(triangle):
        for k, value in enumerate(row):
            entries += 1
            if (value == 1) != (k & (n - k) == 0):
                faults += 1
    print(f"  rows 0..{ROWS - 1}: {entries} entries, {faults} mismatched cells")


def shear(triangle):
    print("THE SHEAR (i, j) -> (i, i + j)")
    side = 1 << SHEAR_LEVEL
    cells = {(i, j) for i in range(side) for j in range(side) if i & j == 0}
    image = {(i, i + j) for (i, j) in cells}
    odd = {
        (k, n)
        for n, row in enumerate(triangle[: 2 * side - 1])
        for k, value in enumerate(row)
        if value and k < side and n - k < side
    }
    print(
        f"  L={SHEAR_LEVEL}: {len(cells)} cells map onto {len(image)} points,"
        f" odd entries in range {len(odd)}, bijection {image == odd}"
    )


def gould(triangle):
    print("ANTIDIAGONAL POPULATION")
    faults = sum(1 for n, row in enumerate(triangle) if sum(row) != 1 << popcount(n))
    print(f"  row sums == 2^popcount(n) for n = 0..{ROWS - 1}, {faults} faults")
    side = 1 << TRUNCATION_LEVEL
    counts = [0] * (2 * side - 1)
    for i in range(side):
        for j in range(side):
            if i & j == 0:
                counts[i + j] += 1
    inside = all(counts[n] == 1 << popcount(n) for n in range(side))
    outside = all(counts[n] < 1 << popcount(n) for n in range(side, 2 * side - 1))
    print(
        f"  inside [0, 2^{TRUNCATION_LEVEL})^2 the count is 2^popcount(n)"
        f" for n < 2^L {inside}, strictly smaller above {outside}"
    )


def oeis(triangle):
    print("OEIS")
    flat = [value for row in triangle for value in row]
    sums = [sum(row) for row in triangle]
    terms = read_live(GOULD_URL)
    if terms is None:
        print("  A001316 live read unavailable, its b-file is too large to keep here")
    else:
        wrong = sum(
            1
            for n, (index, value) in enumerate(terms)
            if index != n or value != 1 << popcount(n)
        )
        print(
            f"  A001316 live: {len(terms)} terms, n = 0..{len(terms) - 1},"
            f" {wrong} differences"
        )
        head = [value for _, value in terms[:ROWS]]
        print(f"  A001316 first {ROWS} terms == our row sums: {head == sums}")
    cache = read_cache(GASKET_CACHE)
    gasket = [("cache", cache)]
    live = read_live(GASKET_URL)
    if live is None:
        print("  A047999 live read unavailable, cache only")
    else:
        gasket.append(("live", live))
    for label, terms in gasket:
        wrong = sum(
            1
            for n, (index, value) in enumerate(terms)
            if index != n or value != flat[n]
        )
        rows = triangular_rows(len(terms))
        print(
            f"  A047999 {label}: {len(terms)} terms, rows 0..{rows - 1},"
            f" {wrong} differences"
        )
    if live is not None:
        print(f"  A047999 cache == live term for term: {cache == live}")


def main():
    print(f"DOMAIN levels 1..{TOP_LEVEL}, binomials < {EXACT}, Pascal rows < {ROWS}")
    level_sets()
    kummer()
    triangle = pascal_mod2(ROWS)
    recurrence(triangle)
    shear(triangle)
    gould(triangle)
    oeis(triangle)
    return 0


if __name__ == "__main__":
    sys.exit(main())
