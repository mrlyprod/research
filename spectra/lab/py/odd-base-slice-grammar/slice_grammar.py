from fractions import Fraction

import numpy as np
from mpmath import log as mlog
from mpmath import mp, mpf
from mpmath import sqrt as msqrt

mp.dps = 40


def digit_sums(b, rule):
    if rule == "odd":
        marked = {d for d in range(b) if d % 2 == 1}
    else:
        marked = {(b - 1) // 2}
    plain = np.zeros(b, dtype=np.int64)
    special = np.zeros(b, dtype=np.int64)
    for d in range(b):
        if d in marked:
            special[d] = 1
        else:
            plain[d] = 1
    square = np.convolve(plain, plain)
    total = np.convolve(square, plain) + 3 * np.convolve(special, square)
    return [int(v) for v in total]


def at(c, i):
    if i < 0 or i >= len(c):
        return 0
    return c[i]


def transfer(b, c):
    row_mid = (
        at(c, (3 * b - 3) // 2),
        at(c, (3 * b - 1) // 2),
        at(c, (3 * b - 5) // 2),
    )
    row_low = (at(c, (b - 3) // 2), at(c, (b - 1) // 2), at(c, (b - 5) // 2))
    row_high = (
        at(c, (5 * b - 3) // 2),
        at(c, (5 * b - 1) // 2),
        at(c, (5 * b - 5) // 2),
    )
    return row_mid, row_low, row_high


def apply_rows(rows, v):
    return tuple(r[0] * v[0] + r[1] * v[1] + r[2] * v[2] for r in rows)


def tile_counts(b, c, levels):
    rows = transfer(b, c)
    states = [(0, 1, 0), (1, 0, 0), (0, 0, 1)]
    hexes = [states[1][0]]
    tris = [states[0][0] + states[2][0]]
    for _ in range(levels):
        states = [apply_rows(rows, v) for v in states]
        hexes.append(states[1][0])
        tris.append(states[0][0] + states[2][0])
    return hexes, tris


def layer_counts(b, c, level):
    dist = {0: 1}
    nz = [(s, v) for s, v in enumerate(c) if v]
    for k in range(level):
        shift = b**k
        nxt = {}
        for s0, v0 in dist.items():
            for s1, v1 in nz:
                key = s0 + s1 * shift
                nxt[key] = nxt.get(key, 0) + v0 * v1
        dist = nxt
    m = (3 * b**level - 1) // 2
    return dist.get(m - 1, 0), dist.get(m, 0) + dist.get(m - 2, 0)


def substitution(b, c, levels=9):
    hexes, tris = tile_counts(b, c, levels)
    m00, m10 = hexes[1], tris[1]
    denom = tris[1]
    m01 = Fraction(hexes[2] - m00 * hexes[1], denom)
    m11 = Fraction(tris[2] - m10 * hexes[1], denom)
    mat = [[Fraction(m00), m01], [Fraction(m10), m11]]
    ok = all(x.denominator == 1 for row in mat for x in row)
    for n in range(1, levels + 1):
        ph = mat[0][0] * hexes[n - 1] + mat[0][1] * tris[n - 1]
        pt = mat[1][0] * hexes[n - 1] + mat[1][1] * tris[n - 1]
        if ph != hexes[n] or pt != tris[n]:
            ok = False
    ints = [[int(x) for x in row] for row in mat] if ok else None
    return ints, hexes, tris, ok


def spectral(tr, det, b):
    lam = (mpf(tr) + msqrt(mpf(tr) ** 2 - 4 * det)) / 2
    return lam, mlog(lam) / mlog(b)


def side(tr, det, fill, b):
    q = Fraction(fill, b)
    disc = Fraction(tr * tr - 4 * det)
    u = Fraction(tr) - 2 * q
    if u > 0:
        return 1
    if u == 0:
        return 1 if disc > 0 else 0
    val = disc - u * u
    return 1 if val > 0 else (0 if val == 0 else -1)


def closed_form(b):
    if b % 4 == 3:
        return [
            [Fraction(3 * (b + 1) * (3 * b - 1), 16), Fraction((b + 1) * (b + 5), 32)],
            [Fraction(3 * (b + 1) ** 2, 8), Fraction(3 * (b + 1) ** 2, 16)],
        ]
    return [
        [Fraction(3 * b * b + 6 * b + 7, 16), Fraction(3 * (b - 1) * (b + 3), 32)],
        [Fraction(3 * (b - 1) * (3 * b + 5), 8), Fraction((b + 3) ** 2, 16)],
    ]


def rule_text(tr, det):
    return "x%d %s%d" % (tr, "+" if -det >= 0 else "-", abs(det))


def report(b, rule, levels=9):
    c = digit_sums(b, rule)
    fill = sum(c)
    mat, hexes, tris, ok = substitution(b, c, levels)
    tr = mat[0][0] + mat[1][1]
    det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
    lam, dim = spectral(tr, det, b)
    solid = mlog(mpf(fill)) / mlog(b)
    return {
        "b": b,
        "rule": rule,
        "fill": fill,
        "cells": b**3,
        "matrix": mat,
        "closes": ok,
        "trace": tr,
        "det": det,
        "recurrence": rule_text(tr, det),
        "lam": lam,
        "dim": dim,
        "solid": solid,
        "minus_one": solid - 1,
        "hexes": hexes,
        "tris": tris,
        "side": side(tr, det, fill, b),
    }


def main():
    print("SLICE GRAMMAR, ODD BASES, AT MOST ONE ODD COORDINATE")
    bases = [3, 5, 7, 9]
    rows = {b: report(b, "odd") for b in bases}
    for b in bases:
        r = rows[b]
        print(
            "b=%d fill=%d/%d matrix=%s closes=%s rule=%s dim_slice=%.4f d-1=%.4f"
            % (
                b,
                r["fill"],
                r["cells"],
                r["matrix"],
                r["closes"],
                r["recurrence"],
                float(r["dim"]),
                float(r["minus_one"]),
            )
        )
    print("four rules: %s" % " / ".join(rows[b]["recurrence"] for b in bases))
    print("dimensions: %s" % " / ".join("%.4f" % float(rows[b]["dim"]) for b in bases))
    print("d-1: %s" % " / ".join("%.4f" % float(rows[b]["minus_one"]) for b in bases))

    print()
    print("BASE THREE EXTERNAL TARGET")
    print("matrix %s" % rows[3]["matrix"])
    print("hexagons %s" % ", ".join(str(v) for v in rows[3]["hexes"][:7]))
    print("triangles %s" % ", ".join(str(v) for v in rows[3]["tris"][:7]))

    print()
    print("INDEPENDENT LAYER CENSUS, NO TILE REDUCTION")
    for b, top in ((3, 6), (5, 5), (7, 4), (9, 4)):
        c = digit_sums(b, "odd")
        hexes, tris = tile_counts(b, c, top)
        agree = True
        for n in range(top + 1):
            dh, dt = layer_counts(b, c, n)
            if (dh, dt) != (hexes[n], tris[n]):
                agree = False
        print("b=%d levels 0..%d agree=%s" % (b, top, agree))

    print()
    print("CLOSED FORM AGAINST CENSUS")
    bad = []
    for b in range(3, 22, 2):
        r = report(b, "odd")
        cf = closed_form(b)
        got = [[Fraction(x) for x in row] for row in r["matrix"]]
        if cf != got or not r["closes"]:
            bad.append(b)
    print("odd bases 3..21 matched=%d mismatched=%s" % (10 - len(bad), bad))

    print()
    print("MOD FOUR SPLIT")
    wrong = []
    tested = 0
    for b in range(3, 402, 2):
        r = report(b, "odd", levels=5)
        tested += 1
        want = 1 if b % 4 == 3 else -1
        if r["side"] != want or not r["closes"]:
            wrong.append(b)
    print("odd bases 3..401 tested=%d exceptions=%s" % (tested, wrong))

    print()
    print("BASE FIVE, THE TWO DIGIT RULES")
    odd5 = rows[5]
    mid5 = report(5, "middle")
    print(
        "middle-digit rule fills %d of %d dim_slice=%.6f d-1=%.6f excess=%+.3e"
        % (
            mid5["fill"],
            mid5["cells"],
            float(mid5["dim"]),
            float(mid5["minus_one"]),
            float(mid5["dim"] - mid5["minus_one"]),
        )
    )
    print(
        "odd-coordinate rule fills %d of %d dim_slice=%.4f d-1=%.4f excess=%+.3e"
        % (
            odd5["fill"],
            odd5["cells"],
            float(odd5["dim"]),
            float(odd5["minus_one"]),
            float(odd5["dim"] - odd5["minus_one"]),
        )
    )
    print("middle-digit matrix %s rule %s" % (mid5["matrix"], mid5["recurrence"]))


main()
