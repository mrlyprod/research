import time
from itertools import product
from math import gcd

import numpy as np

BASE2_DIMS = (2, 3, 4)
IDENTITY_CASES = ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (5, 2), (7, 2))
COMPOSITE_BASES = (4, 6, 8, 9, 10, 12)
PRIME_BASES = (2, 3, 5, 7, 11, 13)
CLAIMED_TABLE = {2: (2, 2, 3, 4), 3: (4, 4, 6, 6)}


def cells(q, d):
    return list(product(range(q), repeat=d))


def strides(q, d):
    return [q ** (d - 1 - a) for a in range(d)]


def mod_inverse(a, q):
    old_r, r = a % q, q
    old_s, s = 1, 0
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
    return old_s % q


def vandermonde(q):
    return [[pow(i, j, q) for j in range(q)] for i in range(q)]


def vandermonde_det_mod(q):
    out = 1
    for i in range(q):
        for j in range(i + 1, q):
            out = out * (j - i) % q
    return out


def invert_mod(matrix, q):
    n = len(matrix)
    aug = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if gcd(aug[r][col] % q, q) == 1), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        inv = mod_inverse(aug[col][col], q)
        aug[col] = [x * inv % q for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] % q:
                factor = aug[r][col]
                aug[r] = [(x - factor * y) % q for x, y in zip(aug[r], aug[col])]
    return [row[n:] for row in aug]


def axis_transform(values, q, d, rows):
    size = q**d
    coeff = list(values)
    for stride in strides(q, d):
        nxt = [0] * size
        for base in range(size):
            if (base // stride) % q:
                continue
            line = [coeff[base + t * stride] for t in range(q)]
            for e in range(q):
                nxt[base + e * stride] = sum(rows[e][t] * line[t] for t in range(q)) % q
        coeff = nxt
    return coeff


def gfq_eval(coeff, q, cs):
    out = []
    for x in cs:
        acc = 0
        for value, e in zip(coeff, cs):
            if value:
                term = value
                for xi, ei in zip(x, e):
                    term = term * pow(xi, ei, q) % q
                acc = (acc + term) % q
        out.append(acc)
    return out


def int_mobius(values, q, d):
    size = q**d
    coeff = list(values)
    for stride in strides(q, d):
        for level in range(q - 1, 0, -1):
            for index in range(size):
                if (index // stride) % q == level:
                    coeff[index] -= coeff[index - stride]
    return coeff


def int_eval(coeff, cs):
    out = []
    for x in cs:
        acc = 0
        for value, m in zip(coeff, cs):
            if value and all(xi >= mi for xi, mi in zip(x, m)):
                acc += value
        out.append(acc)
    return out


def degree_of(coeff, cs):
    return max((sum(e) for value, e in zip(coeff, cs) if value), default=-1)


def tensor_power(m, d, q):
    out = np.array(m, dtype=np.int64)
    for _ in range(d - 1):
        out = np.kron(out, np.array(m, dtype=np.int64)) % q
    return out


def eval_matrix(q, d):
    cs = cells(q, d)
    return np.array(
        [
            [np.prod([pow(xi, ei, q) for xi, ei in zip(x, e)]) % q for e in cs]
            for x in cs
        ],
        dtype=np.int64,
    )


def int_basis_matrix(q, d):
    cs = cells(q, d)
    return np.array(
        [[int(all(xi >= mi for xi, mi in zip(x, m))) for m in cs] for x in cs],
        dtype=np.int64,
    )


def int_mobius_matrix(q, d):
    n = q**d
    columns = [int_mobius([int(i == j) for i in range(n)], q, d) for j in range(n)]
    return np.array(columns, dtype=np.int64).T


def xor_anf_masks(table, d):
    out = [0] * (1 << d)
    for s in range(1 << d):
        acc = 0
        sub = s
        while True:
            acc ^= table[sub]
            if sub == 0:
                break
            sub = (sub - 1) & s
        out[s] = acc
    return out


def signed_subset_real(table, d):
    out = [0] * (1 << d)
    for s in range(1 << d):
        acc = 0
        sub = s
        while True:
            if table[sub]:
                acc += -1 if (s ^ sub).bit_count() & 1 else 1
            if sub == 0:
                break
            sub = (sub - 1) & s
        out[s] = acc
    return out


def reversed_mask(s, d):
    return sum(((s >> i) & 1) << (d - 1 - i) for i in range(d))


def code_to_table(code, n):
    return [(code >> i) & 1 for i in range(n)]


def named(q, d):
    cs = cells(q, d)
    rules = (
        ("void", lambda c: all(v == c[0] for v in c)),
        ("tree", lambda c: all(a == 0 or c[a] == 0 for a in range(d))),
        ("carpet", lambda c: sum(c) <= 1),
        ("net", lambda c: sum(c) >= d - 1),
    )
    return [(name, [int(rule(c)) for c in cs]) for name, rule in rules]


def fractal_shapes():
    carpet = [int(c != (1, 1)) for c in cells(3, 2)]
    sponge = [int(c.count(1) <= 1) for c in cells(3, 3)]
    return (("sierpinski carpet", 2, carpet), ("menger sponge", 3, sponge))


def histogram(counts):
    return ", ".join(f"{k}: {counts[k]}" for k in sorted(counts))


def verdict(ok):
    return "PASS" if ok else "FAIL"


def block_inverse():
    print("THE INVERSE VANDERMONDE OVER GF(q)")
    for q in PRIME_BASES:
        v = vandermonde(q)
        inv = invert_mod(v, q)
        ok = inv is not None and np.array_equal(
            np.array(inv) @ np.array(v) % q, np.eye(q, dtype=np.int64)
        )
        print(
            f"  q={q}  gcd(det, q) = {gcd(vandermonde_det_mod(q), q)}  invertible {inv is not None}  {verdict(ok)}"
        )
    inv2 = invert_mod(vandermonde(2), 2)
    print(f"  q=2  inverse {inv2}  {verdict(inv2 == [[1, 0], [1, 1]])}")
    for q in COMPOSITE_BASES:
        inv = invert_mod(vandermonde(q), q)
        print(
            f"  q={q}  gcd(det, q) = {gcd(vandermonde_det_mod(q), q)}  invertible {inv is not None}  {verdict(inv is None)}"
        )
    print()


def block_identity():
    print("THE TENSORED TRANSFORM INVERTS EVALUATION, ALL VALUE TABLES AT ONCE")
    for q, d in IDENTITY_CASES:
        n = q**d
        t = tensor_power(invert_mod(vandermonde(q), q), d, q)
        e = eval_matrix(q, d)
        eye = np.eye(n, dtype=np.int64)
        ok = np.array_equal(t @ e % q, eye) and np.array_equal(e @ t % q, eye)
        print(
            f"  q={q} D={d}  T E = E T = I on {n} coordinates  covers all {q}^{n} tables  {verdict(ok)}"
        )
    print("THE INTEGER MOBIUS INVERTS THE DOWNWARD-CLOSED BASIS")
    for q, d in (
        (2, 2),
        (2, 3),
        (2, 4),
        (3, 2),
        (3, 3),
        (4, 2),
        (4, 3),
        (6, 2),
        (9, 2),
    ):
        n = q**d
        ok = np.array_equal(
            int_basis_matrix(q, d) @ int_mobius_matrix(q, d), np.eye(n, dtype=np.int64)
        )
        print(f"  q={q} D={d}  B M = I over Z on {n} coordinates  {verdict(ok)}")
    print()


def block_base2():
    print("BASE 2, EVERY DESIGN, AGAINST THE CLASSICAL XOR ANF")
    rows = invert_mod(vandermonde(2), 2)
    for d in BASE2_DIMS:
        cs = cells(2, d)
        n = 1 << d
        bad_rt = bad_anf = bad_naive = bad_int = 0
        for code in range(1 << n):
            table = code_to_table(code, n)
            coeff = axis_transform(table, 2, d, rows)
            if gfq_eval(coeff, 2, cs) != table:
                bad_rt += 1
            masked = xor_anf_masks([table[reversed_mask(s, d)] for s in range(n)], d)
            support = [int(v != 0) for v in coeff]
            if [masked[reversed_mask(s, d)] for s in range(n)] != support:
                bad_anf += 1
            if masked != support:
                bad_naive += 1
            integer = int_mobius(table, 2, d)
            real = signed_subset_real([table[reversed_mask(s, d)] for s in range(n)], d)
            if (
                int_eval(integer, cs) != table
                or [real[reversed_mask(s, d)] for s in range(n)] != integer
            ):
                bad_int += 1
        ok = bad_rt == 0 and bad_anf == 0 and bad_int == 0
        print(
            f"  D={d}  designs {1 << n:>5}  GF(2) roundtrip fails {bad_rt}  XOR ANF diffs after reversal {bad_anf}  diffs without reversal {bad_naive}  integer fails {bad_int}  {verdict(ok)}"
        )
    print()


def block_base3():
    print("BASE 3, D=2, ALL 512 DESIGNS")
    cs = cells(3, 2)
    rows = invert_mod(vandermonde(3), 3)
    bad_gf = bad_int = 0
    hist_gf, hist_int = {}, {}
    for code in range(512):
        table = code_to_table(code, 9)
        coeff = axis_transform(table, 3, 2, rows)
        if gfq_eval(coeff, 3, cs) != table:
            bad_gf += 1
        integer = int_mobius(table, 3, 2)
        if int_eval(integer, cs) != table:
            bad_int += 1
        dg, di = degree_of(coeff, cs), degree_of(integer, cs)
        hist_gf[dg] = hist_gf.get(dg, 0) + 1
        hist_int[di] = hist_int.get(di, 0) + 1
    print(
        f"  designs 512  GF(3) roundtrip fails {bad_gf}  integer roundtrip fails {bad_int}  {verdict(bad_gf == 0 and bad_int == 0)}"
    )
    print(f"  GF(3) degree histogram    {histogram(hist_gf)}")
    print(f"  integer degree histogram  {histogram(hist_int)}")
    print(
        f"  designs of GF(3) degree 1: {hist_gf.get(1, 0)}  {verdict(hist_gf.get(1, 0) == 0)}"
    )
    print()


def block_table():
    print("THE BASE-3 DEGREE TABLE OF THE FOUR RULES")
    rows = invert_mod(vandermonde(3), 3)
    ok = True
    for d in (2, 3):
        cs = cells(3, d)
        for (name, table), claim in zip(named(3, d), CLAIMED_TABLE[d]):
            dg = degree_of(axis_transform(table, 3, d, rows), cs)
            di = degree_of(int_mobius(table, 3, d), cs)
            ok = ok and dg == claim
            print(
                f"  D={d} {name:<7} cells {sum(table):>2}/{3**d}  GF(3) deg {dg}  page {claim}  integer deg {di}  {verdict(dg == claim)}"
            )
    print(f"  all eight entries  {verdict(ok)}")
    print()
    print("THE FRACTALS THE ROWS ARE NAMED AFTER")
    for name, d, shape in fractal_shapes():
        cs = cells(3, d)
        dg = degree_of(axis_transform(shape, 3, d, rows), cs)
        di = degree_of(int_mobius(shape, 3, d), cs)
        rule = dict(named(3, d))["carpet"]
        rdg = degree_of(axis_transform(rule, 3, d, rows), cs)
        print(
            f"  {name:<18} D={d}  cells {sum(shape):>2}/{3**d}  GF(3) deg {dg}  integer deg {di}  ceiling D(q-1) = {2 * d}"
        )
        print(
            f"  carpet row         D={d}  cells {sum(rule):>2}/{3**d}  GF(3) deg {rdg}  same set as the fractal {rule == shape}"
        )
    print()


def main():
    start = time.time()
    block_inverse()
    block_identity()
    block_base2()
    block_base3()
    block_table()
    print(
        f"domain: base 2 at D = {', '.join(map(str, BASE2_DIMS))} exhaustive, base 3 at D = 2 exhaustive, base 3 at D = 3 by the matrix identity"
    )
    print(f"wall {time.time() - start:.1f}s")


if __name__ == "__main__":
    main()
