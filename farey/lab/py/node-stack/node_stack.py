import hashlib
from fractions import Fraction
from math import floor, gcd, lcm
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
PARITY_CHECK_N = 15
PARITY_RENDER_N = 31
CARPET_CHECKS = ((12, 1), (6, 2))
CARPET_CORNER_N = 12
FAREY_MAX_N = 30
DISCREPANCY_QS = (30, 60, 90)
TOP = 6
RENDER_R = 1024
SAVE_R = 512
SHEET_R = 256
CARPET_RENDER_L = 2
CARPET_RENDER_N = 12
INK_FLOOR = 70


def parity_inked(i, j, n):
    return 0 <= i < n and 0 <= j < n and i % 2 == 1 and j % 2 == 1


def parity_edges(n):
    v = {
        (k, j)
        for k in range(n + 1)
        for j in range(n)
        if parity_inked(k - 1, j, n) != parity_inked(k, j, n)
    }
    h = {
        (i, m)
        for m in range(n + 1)
        for i in range(n)
        if parity_inked(i, m - 1, n) != parity_inked(i, m, n)
    }
    return v, h


def parity_edges_form(n):
    v = {(k, j) for k in range(1, n) for j in range(1, n, 2)}
    h = {(i, m) for m in range(1, n) for i in range(1, n, 2)}
    return v, h


def parity_corners(n):
    return {
        (Fraction(k, n), Fraction(m, n))
        for i in range(1, n, 2)
        for j in range(1, n, 2)
        for k in (i, i + 1)
        for m in (j, j + 1)
    }


def parity_corners_form(n):
    return {(Fraction(k, n), Fraction(m, n)) for k in range(1, n) for m in range(1, n)}


def parity_corner_stack(n_max):
    out = {}
    for n in range(1, n_max + 1, 2):
        for p in parity_corners(n):
            out[p] = out.get(p, 0) + 1
    return out


def parity_corner_form(b, d, n_max):
    m = lcm(b, d)
    if m % 2 == 0:
        return 0
    return (n_max // m + 1) // 2


def parity_corner_predicted(n_max):
    out = {}
    for b in range(1, n_max + 1, 2):
        for d in range(1, n_max + 1, 2):
            if lcm(b, d) > n_max:
                continue
            bright = parity_corner_form(b, d, n_max)
            for a in range(1, b):
                if gcd(a, b) != 1:
                    continue
                for c in range(1, d):
                    if gcd(c, d) == 1:
                        out[(Fraction(a, b), Fraction(c, d))] = bright
    return out


def carpet_inked(i, j, level):
    for _ in range(level):
        if i % 3 == 1 and j % 3 == 1:
            return False
        i //= 3
        j //= 3
    return True


def carpet_cell(k, j, level, side):
    s = 3**level
    return 0 <= k < side and 0 <= j < side and carpet_inked(k % s, j % s, level)


def carpet_edges(n, level):
    side = 3**level * n
    v = {
        (k, j)
        for k in range(side + 1)
        for j in range(side)
        if carpet_cell(k - 1, j, level, side) != carpet_cell(k, j, level, side)
    }
    h = {
        (i, m)
        for m in range(side + 1)
        for i in range(side)
        if carpet_cell(i, m - 1, level, side) != carpet_cell(i, m, level, side)
    }
    return v, h


def edge_residues(level):
    s = 3**level
    rows = {}
    for j in range(s):
        r = [
            x
            for x in range(s)
            if carpet_inked((j - 1) % s, x, level) != carpet_inked(j, x, level)
        ]
        if r:
            rows[j] = r
    return sorted(rows), rows


def carpet_edge_lines(n, level):
    side = 3**level * n
    v, _ = carpet_edges(n, level)
    return {Fraction(k, side) for (k, j) in v if 0 < k < side}


def carpet_line_lit(a, b, n, level, residues):
    s = 3**level
    if (s * n) % b:
        return False
    return (a * (s * n // b)) % s in residues


def v3(x):
    t = 0
    while x % 3 == 0:
        x //= 3
        t += 1
    return t


def carpet_lines_reach(n_max, level):
    out = set()
    for b in range(3, 3**level * n_max + 1, 3):
        if b // 3 ** min(v3(b), level) > n_max:
            continue
        for a in range(1, b):
            if gcd(a, b) == 1:
                out.add(Fraction(a, b))
    return out


def carpet_edge_segments(n_max, level):
    side_lines = {}
    for n in range(1, n_max + 1):
        side = 3**level * n
        v, _ = carpet_edges(n, level)
        for k, j in v:
            if k == 0 or k == side:
                continue
            side_lines.setdefault(Fraction(k, side), []).append(
                (Fraction(j, side), Fraction(j + 1, side))
            )
    return side_lines


def top_segments(lines, count):
    out = []
    for x, ivs in lines.items():
        pts = sorted({p for iv in ivs for p in iv})
        for i in range(len(pts) - 1):
            lo, hi = pts[i], pts[i + 1]
            b = sum(1 for (u, w) in ivs if u <= lo and hi <= w)
            if b:
                out.append((b, hi - lo, x, lo, hi))
    out.sort(key=lambda t: (-t[0], -t[1], t[2], t[3]))
    return out[:count]


def parity_edge_lines(n):
    v, _ = parity_edges(n)
    return {Fraction(k, n) for (k, j) in v}


def parity_line_brightness(b, n_max):
    if b % 2 == 0:
        return 0
    return (n_max // b + 1) // 2


def carpet_line_brightness(b, n_max, level):
    s = v3(b)
    if s == 0:
        return 0
    return n_max // (b // 3 ** min(s, level)) - n_max // b


def check_line_brightness_parity(n_max):
    count = {}
    for n in range(1, n_max + 1, 2):
        for f in parity_edge_lines(n):
            count[f] = count.get(f, 0) + 1
    bad = sum(
        1 for f, c in count.items() if parity_line_brightness(f.denominator, n_max) != c
    )
    invented = sum(
        1
        for b in range(1, n_max + 1)
        for a in range(1, b)
        if gcd(a, b) == 1
        and parity_line_brightness(b, n_max) != count.get(Fraction(a, b), 0)
    )
    return len(count), bad + invented


def check_line_brightness_carpet(n_max, level):
    count = {}
    for n in range(1, n_max + 1):
        for f in carpet_edge_lines(n, level):
            count[f] = count.get(f, 0) + 1
    bad = sum(
        1
        for f, c in count.items()
        if carpet_line_brightness(f.denominator, n_max, level) != c
    )
    invented = sum(
        1
        for b in range(1, 3**level * n_max + 1)
        for a in range(1, b)
        if gcd(a, b) == 1
        and carpet_line_brightness(b, n_max, level) != count.get(Fraction(a, b), 0)
    )
    return len(count), bad + invented


def parity_segment_check(n_max):
    odds = list(range(1, n_max + 1, 2))
    lines = {}
    for n in odds:
        v, _ = parity_edges(n)
        for k, j in v:
            lines.setdefault(Fraction(k, n), []).append(
                (n, Fraction(j, n), Fraction(j + 1, n))
            )
    pts = sorted({Fraction(j, n) for n in odds for j in range(n + 1)})
    mids = [(pts[i] + pts[i + 1]) / 2 for i in range(len(pts) - 1)]
    bad = 0
    for x, ivs in lines.items():
        b = x.denominator
        for y in mids:
            literal = sum(1 for (n, u, w) in ivs if u <= y <= w)
            form = sum(1 for n in odds if n % b == 0 and floor(n * y) % 2 == 1)
            bad += literal != form
    return len(lines) * len(mids), bad


def carpet_segment_check(n_max, level):
    s = 3**level
    _, rows = edge_residues(level)
    lines = {}
    for n in range(1, n_max + 1):
        side = s * n
        v, _ = carpet_edges(n, level)
        for k, j in v:
            if 0 < k < side:
                lines.setdefault(Fraction(k, side), []).append(
                    (n, Fraction(j, side), Fraction(j + 1, side))
                )
    pts = sorted(
        {Fraction(j, s * n) for n in range(1, n_max + 1) for j in range(s * n + 1)}
    )
    mids = [(pts[i] + pts[i + 1]) / 2 for i in range(len(pts) - 1)]
    bad = 0
    for x, ivs in lines.items():
        a, b = x.numerator, x.denominator
        for y in mids:
            literal = sum(1 for (n, u, w) in ivs if u <= y <= w)
            form = 0
            for n in range(1, n_max + 1):
                side = s * n
                if side % b:
                    continue
                j = (a * (side // b)) % s
                if j in rows and floor(side * y) % s in rows[j]:
                    form += 1
            bad += literal != form
    return len(lines) * len(mids), bad


def carpet_corners(n, level):
    side = 3**level * n
    out = set()
    for k in range(side + 1):
        for m in range(side + 1):
            if any(
                carpet_cell(k - 1 + du, m - 1 + dv, level, side)
                for du in (0, 1)
                for dv in (0, 1)
            ):
                out.add((Fraction(k, side), Fraction(m, side)))
    return out


def carpet_corner_stack(n_max, level):
    out = {}
    for n in range(1, n_max + 1):
        for p in carpet_corners(n, level):
            out[p] = out.get(p, 0) + 1
    return out


def carpet_corner_form(b, d, n_max):
    m = lcm(b, d)
    return n_max // (m // gcd(m, 3))


def totients(q):
    phi = list(range(q + 1))
    for p in range(2, q + 1):
        if phi[p] == p:
            for k in range(p, q + 1, p):
                phi[k] -= phi[k] // p
    return phi


def farey_sizes(q):
    phi = totients(q)
    return sum(phi[1:]), sum(phi[k] for k in range(3, q + 1, 3))


def farey_list(q, step=1):
    out = []
    for b in range(step, q + 1, step):
        for a in range(1, b + 1):
            if gcd(a, b) == 1:
                out.append(Fraction(a, b))
    out.sort()
    return out


def landau(nodes):
    m = len(nodes)
    total = Fraction(0)
    for i, f in enumerate(nodes, start=1):
        total += abs(f - Fraction(i, m))
    return total


def digest(lines):
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def stack_digest(stack):
    return digest(["%s,%s:%d" % (p[0], p[1], b) for p, b in sorted(stack.items())])


def check_parity_edges(n_max):
    bad = 0
    for n in range(1, n_max + 1, 2):
        if parity_edges(n) != parity_edges_form(n):
            bad += 1
    return bad


def check_parity_corners(n_max):
    bad = 0
    for n in range(1, n_max + 1, 2):
        if parity_corners(n) != parity_corners_form(n):
            bad += 1
    return bad


def check_carpet_lines(n_max, level):
    residues, _ = edge_residues(level)
    rs = set(residues)
    bad = 0
    for n in range(1, n_max + 1):
        literal = carpet_edge_lines(n, level)
        form = {
            f
            for f in literal | carpet_lines_reach(n_max, level)
            if carpet_line_lit(f.numerator, f.denominator, n, level, rs)
        }
        if literal != form:
            bad += 1
    return bad


def blockmax(a, side):
    f = a.shape[0] // side
    return a.reshape(side, f, side, f).max(axis=(1, 3))


def grey(a):
    m = int(a.max())
    v = np.zeros(a.shape, dtype=np.float64)
    lit = a > 0
    v[lit] = INK_FLOOR + (255.0 - INK_FLOOR) * (a[lit] - 1) / max(m - 1, 1)
    return (255 - np.round(v)).astype(np.uint8)


def save_png(acc, path, side):
    Image.fromarray(grey(blockmax(acc, side)), mode="L").save(path, optimize=True)
    return path.stat().st_size, int(acc.max())


def contact_sheet(accs, path, side):
    tiles = [grey(blockmax(a, side)) for a in accs]
    sheet = np.vstack([np.hstack(tiles[:2]), np.hstack(tiles[2:])])
    Image.fromarray(sheet, mode="L").save(path, optimize=True)
    return path.stat().st_size


def edge_raster(layers, r):
    acc = np.zeros((r, r), dtype=np.int32)
    for side, v, h in layers:
        for k, j in v:
            c = min(k * r // side, r - 1)
            y0 = j * r // side
            y1 = min(max((j + 1) * r // side, y0 + 1), r)
            acc[y0:y1, c] += 1
        for i, m in h:
            c = min(m * r // side, r - 1)
            x0 = i * r // side
            x1 = min(max((i + 1) * r // side, x0 + 1), r)
            acc[c, x0:x1] += 1
    return acc


def corner_raster(stack, r, scale):
    acc = np.zeros((r, r), dtype=np.int32)
    for (x, y), b in stack.items():
        cx = min(x.numerator * r // x.denominator, r - 1)
        cy = min(y.numerator * r // y.denominator, r - 1)
        rad = 1 + int(scale * b)
        x0, x1 = max(cx - rad, 0), min(cx + rad + 1, r)
        y0, y1 = max(cy - rad, 0), min(cy + rad + 1, r)
        np.maximum(acc[y0:y1, x0:x1], b, out=acc[y0:y1, x0:x1])
    return acc


def figures():
    parity_layers = [(n,) + parity_edges(n) for n in range(1, PARITY_RENDER_N + 1, 2)]
    carpet_layers = [
        (3**CARPET_RENDER_L * n,) + carpet_edges(n, CARPET_RENDER_L)
        for n in range(1, CARPET_RENDER_N + 1)
    ]
    accs = [
        edge_raster(parity_layers, RENDER_R),
        corner_raster(parity_corner_stack(PARITY_RENDER_N), RENDER_R, 1.0),
        edge_raster(carpet_layers, RENDER_R),
        corner_raster(
            carpet_corner_stack(CARPET_RENDER_N, CARPET_RENDER_L), RENDER_R, 0.25
        ),
    ]
    names = [
        "edges-parity.png",
        "corners-parity.png",
        "edges-carpet.png",
        "corners-carpet.png",
    ]
    for acc, name in zip(accs, names):
        size, peak = save_png(acc, HERE / name, SAVE_R)
        print("figure", name, "bytes", size, "peak", peak)
    print(
        "figure nodes-sheet.png bytes",
        contact_sheet(accs, HERE / "nodes-sheet.png", SHEET_R),
    )


def main():
    print(
        "parity edge set literal against form, odd n <= %d, mismatches"
        % PARITY_CHECK_N,
        check_parity_edges(PARITY_CHECK_N),
    )
    print(
        "parity corner set literal against form, odd n <= %d, mismatches"
        % PARITY_CHECK_N,
        check_parity_corners(PARITY_CHECK_N),
    )
    literal = parity_corner_stack(PARITY_CHECK_N)
    predicted = parity_corner_predicted(PARITY_CHECK_N)
    print(
        "parity corner stack N = %d, lit points" % PARITY_CHECK_N,
        len(literal),
        "predicted",
        len(predicted),
        "missed",
        len(set(literal) - set(predicted)),
        "invented",
        len(set(predicted) - set(literal)),
        "value mismatches",
        sum(1 for p in literal if literal[p] != predicted.get(p)),
    )
    print("parity corner stack digest", stack_digest(literal), stack_digest(predicted))
    top = sorted(literal.items(), key=lambda kv: (-kv[1], kv[0]))[:TOP]
    print("parity corner top", ["%s,%s:%d" % (p[0], p[1], b) for p, b in top])
    print(
        "parity edge lines N = %d, lit lines and form breaches" % PARITY_CHECK_N,
        check_line_brightness_parity(PARITY_CHECK_N),
    )
    print(
        "parity edge segments N = %d, tests and form breaches" % PARITY_CHECK_N,
        parity_segment_check(PARITY_CHECK_N),
    )
    for level in (1, 2):
        js, rows = edge_residues(level)
        print(
            "carpet level %d edge residues J_%d" % (level, level),
            js,
            "full nonzero",
            js == list(range(1, 3**level)),
            "rows",
            {j: rows[j] for j in js},
        )
    for n_max, level in CARPET_CHECKS:
        print(
            "carpet edge lines N = %d L = %d, layer mismatches" % (n_max, level),
            check_carpet_lines(n_max, level),
        )
        lit = set()
        for n in range(1, n_max + 1):
            lit |= carpet_edge_lines(n, level)
        reach = carpet_lines_reach(n_max, level)
        print(
            "carpet lit lines N = %d L = %d" % (n_max, level),
            len(lit),
            "reach form",
            len(reach),
            "equal",
            lit == reach,
        )
        print(
            "carpet line brightness N = %d L = %d, lit lines and form breaches"
            % (n_max, level),
            check_line_brightness_carpet(n_max, level),
        )
        segs = top_segments(carpet_edge_segments(n_max, level), TOP)
        print(
            "carpet top segments N = %d L = %d" % (n_max, level),
            ["%d at x = %s on [%s, %s]" % (b, x, lo, hi) for (b, w, x, lo, hi) in segs],
        )
        print(
            "carpet edge segments N = %d L = %d, tests and form breaches"
            % (n_max, level),
            carpet_segment_check(n_max, level),
        )
    stack = carpet_corner_stack(CARPET_CORNER_N, 1)
    bad = sum(
        1
        for p, b in stack.items()
        if carpet_corner_form(p[0].denominator, p[1].denominator, CARPET_CORNER_N) != b
    )
    print(
        "carpet corner stack N = %d L = 1, lit corners" % CARPET_CORNER_N,
        len(stack),
        "form mismatches",
        bad,
        "digest",
        stack_digest(stack),
    )
    top = sorted(stack.items(), key=lambda kv: (-kv[1], kv[0]))[:TOP]
    print("carpet corner top", ["%s,%s:%d" % (p[0], p[1], b) for p, b in top])
    tile = carpet_corners(1, 2)
    grid = {(Fraction(k, 9), Fraction(m, 9)) for k in range(10) for m in range(10)}
    print(
        "carpet level 2 tile vertices",
        len(grid),
        "corners",
        len(tile),
        "interior to a hole",
        sorted("%s,%s" % p for p in grid - tile),
    )
    phi = totients(3 * CARPET_CORNER_N)
    print(
        "plain Farey pair count at Q = %d" % (3 * CARPET_CORNER_N),
        (1 + sum(phi[1:])) ** 2,
    )
    print(
        "restricted Farey sizes, L = 1, N = 1..%d" % FAREY_MAX_N,
        [farey_sizes(3 * n) for n in range(1, FAREY_MAX_N + 1)],
    )
    for q in DISCREPANCY_QS:
        plain = landau(farey_list(q))
        restricted = landau(farey_list(q, 3))
        print(
            "landau Q = %d" % q,
            "plain %.6f" % float(plain),
            "restricted %.6f" % float(restricted),
        )
    figures()


if __name__ == "__main__":
    main()
