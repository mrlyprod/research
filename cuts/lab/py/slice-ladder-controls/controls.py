import itertools
import math
from fractions import Fraction

CODES = range(1, 16)


def tile(code):
    return {(j & 1, (j >> 1) & 1) for j in range(4) if code & (1 << j)}


def kron(a, b, side_b):
    return {(ra * side_b + rb, ca * side_b + cb) for ra, ca in a for rb, cb in b}


def word_array(word):
    cells, side = tile(word[0]), 2
    for code in word[1:]:
        cells = kron(cells, tile(code), 2)
        side *= 2
    return cells, side


def profile(cells, side):
    p = [0] * (2 * side - 1)
    for r, c in cells:
        p[r + c] += 1
    return p


def convolve(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                out[i + j] += x * y
    return out


def dilate(p, k):
    out = [0] * ((len(p) - 1) * k + 1)
    for i, x in enumerate(p):
        out[i * k] = x
    return out


def report_profile_identity():
    bad = 0
    for length in (2, 3):
        for word in itertools.product(CODES, repeat=length):
            head, side = word_array(word[:-1])
            tail = tile(word[-1])
            lhs = profile(kron(head, tail, 2), side * 2)
            rhs = convolve(dilate(profile(head, side), 2), profile(tail, 2))
            bad += lhs != rhs
    print(f"profile identity mismatches over words of length 2 and 3: {bad}")


def digit_polynomial(D):
    p = [0] * (2 * D + 1)
    for v in itertools.product(range(3), repeat=D):
        if sum(1 for x in v if x == 1) <= 1:
            p[sum(v)] += 1
    return p


def central_slice(D, level):
    p = digit_polynomial(D)
    acc = [1]
    for j in range(level):
        acc = convolve(acc, dilate(p, 3**j))
    return acc[(len(acc) - 1) // 2]


def cross_section_vertices(D):
    if D % 2 == 0:
        return math.comb(D, D // 2)
    return math.comb(D, (D - 1) // 2) * ((D + 1) // 2)


def report_vertex_identity(top):
    row, bad = [], 0
    for D in range(2, top + 1):
        lhs, rhs = central_slice(D, 1), cross_section_vertices(D)
        bad += lhs != rhs
        row.append(lhs)
    print(
        f"level-1 slice against cross-section vertices, D = 2..{top}: {bad} mismatches"
    )
    print(f"level-1 slice counts D = 2..{top}: {row}")


def fit_order_two(seq):
    s0, s1, s2, s3 = (Fraction(x) for x in seq[:4])
    det = s1 * s1 - s0 * s2
    return (s1 * s2 - s0 * s3) / det, (s1 * s3 - s2 * s2) / det


def report_rung(D, levels):
    seq = [central_slice(D, level) for level in range(1, levels + 1)]
    c1, c2 = fit_order_two(seq)
    holds = all(seq[n] == c1 * seq[n - 1] + c2 * seq[n - 2] for n in range(2, len(seq)))
    root = (float(c1) + math.sqrt(float(c1 * c1 + 4 * c2))) / 2
    print(f"D = {D} census, levels 1..{levels}: {seq}")
    print(
        f"D = {D} recurrence a(n) = {c1}a(n-1) + {c2}a(n-2), "
        f"fitted on four terms, holds on every term: {holds}"
    )
    print(
        f"D = {D} dominant root {root:.9f}, "
        f"slice dimension {math.log(root) / math.log(3):.9f}"
    )


def report_staircase(bases, top):
    for n in range(1, top + 1):
        num = sum(
            (n - j) * math.log(q * q - ((q - 1) // 2) ** 2)
            for j, q in enumerate(bases[:n])
        )
        den = sum((n - j) * math.log(q) for j, q in enumerate(bases[:n]))
        print(f"staircase n = {n}, bases {bases[:n]}, dimension {num / den:.9f}")


def main():
    print("fill assumed for base q is q^2 - ((q-1)/2)^2")
    report_profile_identity()
    report_vertex_identity(14)
    report_rung(4, 6)
    report_staircase([3, 5, 7, 9, 11], 5)


main()
