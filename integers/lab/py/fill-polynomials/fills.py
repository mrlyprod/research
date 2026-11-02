import math
import operator
import os
import time
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import permutations, product

from sympy import Poly
from sympy.abc import n as VAR

KMAX = 10
DISC_DIMS = (2, 3, 4, 5, 6)
PARALLEL_FROM = 5
CHUNKS = 96
LADDER_DEPTH = 400


def corners(d):
    return [tuple((i >> (d - 1 - a)) & 1 for a in range(d)) for i in range(1 << d)]


def block_tables(values, combine):
    tables = []
    for start in range(0, len(values), 8):
        block = values[start : start + 8]
        table = [0] * (1 << len(block))
        for mask in range(1, 1 << len(block)):
            low = mask & -mask
            table[mask] = combine(table[mask ^ low], block[low.bit_length() - 1])
        tables.append(table)
    return tables


def block_apply(tables, code, combine):
    out = 0
    for table in tables:
        out = combine(out, table[code & (len(table) - 1)])
        code >>= 8
    return out


def sum_tables(values):
    return block_tables(values, operator.add)


def or_tables(mapping):
    return block_tables([1 << target for target in mapping], operator.or_)


def block_sum(tables, code):
    return block_apply(tables, code, operator.add)


def block_or(tables, code):
    return block_apply(tables, code, operator.or_)


def fill_tables(d, k):
    return sum_tables([k ** (d - sum(c)) * (k - 1) ** sum(c) for c in corners(d)])


def polymul(a, b):
    out = [a[0] * b[0] * 0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            out[i + j] = out[i + j] + left * right
    return out


def scaled_basis(xs):
    basis = []
    for i, xi in enumerate(xs):
        term = [Fraction(1)]
        for j, xj in enumerate(xs):
            if i != j:
                term = polymul(term, [Fraction(-xj), Fraction(1)])
                term = [value / Fraction(xi - xj) for value in term]
        basis.append(term)
    scale = 1
    for term in basis:
        for value in term:
            scale = scale * value.denominator // math.gcd(scale, value.denominator)
    return scale, [[int(value * scale) for value in term] for term in basis]


def group_maps(d):
    cs = corners(d)
    index = {c: i for i, c in enumerate(cs)}
    maps = []
    for perm in permutations(range(d)):
        for pattern in range(1 << d):
            flip = tuple((pattern >> (d - 1 - a)) & 1 for a in range(d))
            maps.append(
                tuple(index[tuple(c[perm[a]] ^ flip[a] for a in range(d))] for c in cs)
            )
    return maps


def classify(d, total):
    maps = [or_tables(mapping) for mapping in group_maps(d)]
    rep_of = [-1] * total
    reps = []
    for code in range(total):
        if rep_of[code] >= 0:
            continue
        reps.append(code)
        for tables in maps:
            rep_of[block_or(tables, code)] = code
    return reps, rep_of


def a129824(d):
    out = 1
    for w in range(d + 1):
        out *= 1 + math.comb(d, w)
    return out


def sweep(d):
    size = 1 << d
    xs = list(range(1, d + 2))
    scale, basis = scaled_basis(xs)
    fit = [fill_tables(d, k) for k in xs]
    held = [(k, fill_tables(d, k)) for k in range(d + 2, KMAX + 1)]
    pop = sum_tables([1] * size)
    polys = []
    failures = 0
    for code in range(1 << size):
        ys = [block_sum(tables, code) for tables in fit]
        coef = [sum(ys[i] * basis[i][j] for i in range(d + 1)) for j in range(d + 1)]
        if any(c % scale for c in coef):
            failures += 1
            coef = [0] * (d + 1)
        else:
            coef = [c // scale for c in coef]
            if coef[d] != block_sum(pop, code):
                failures += 1
            for k, tables in held:
                if sum(c * k**j for j, c in enumerate(coef)) != block_sum(tables, code):
                    failures += 1
                    break
        polys.append(tuple(coef))
    return polys, failures


def polynomial_report():
    out = []
    for d in range(1, 5):
        polys, failures = sweep(d)
        reps, rep_of = classify(d, len(polys))
        members = {}
        for code, rep in enumerate(rep_of):
            members.setdefault(rep, []).append(code)
        lower = sum(1 for r in reps if len({polys[c][:-1] for c in members[r]}) > 1)
        lead = sum(1 for r in reps if len({polys[c][-1] for c in members[r]}) > 1)
        out.append(
            (
                d,
                len(polys),
                len(reps),
                len(set(polys)),
                a129824(d),
                lower,
                lead,
                failures,
            )
        )
    return out


def sig_bounds(d):
    return [math.comb(d, w) + 1 for w in range(d + 1)]


def unrank(index, bounds):
    f = [0] * len(bounds)
    for i in range(len(bounds) - 1, -1, -1):
        f[i] = index % bounds[i]
        index //= bounds[i]
    return f


def poly_of_sig(f, d):
    coef = [0] * (d + 1)
    for w, fw in enumerate(f):
        if fw:
            for j in range(d - w + 1):
                coef[w + j] += fw * math.comb(d - w, j)
    return coef


def factor_list(coef):
    parts = []
    for g, mult in Poly(list(reversed(coef)), VAR, domain="ZZ").factor_list()[1]:
        parts.extend([g] * mult)
    return parts


def quadratic_discs(coef):
    if not any(coef):
        return []
    out = []
    for g in factor_list(coef):
        if g.degree() == 2:
            a, b, c = g.all_coeffs()
            out.append(int(b * b - 4 * a * c))
    return out


def evaluate_scaled(coef, q):
    top = len(coef) - 1
    return sum(c * (-1) ** i * q ** (top - i) for i, c in enumerate(coef))


def peel_unit_linear(coef):
    c = list(coef)
    while len(c) > 1:
        lead = c[-1]
        taken = 0
        for q in range(1, abs(lead) + 1):
            if lead % q == 0 and evaluate_scaled(c, q) == 0:
                taken = q
                break
        if not taken:
            break
        out = [0] * (len(c) - 1)
        out[0] = c[0]
        for i in range(1, len(out)):
            out[i] = c[i] - taken * out[i - 1]
        c = out
    return c


def remainder_disc(coef):
    rest = peel_unit_linear(coef)
    if len(rest) == 3:
        return rest[1] * rest[1] - 4 * rest[2] * rest[0]
    return None


def disc_chunk(job):
    d, lo, hi = job
    bounds = sig_bounds(d)
    every, remainder = set(), set()
    for index in range(lo, hi):
        coef = poly_of_sig(unrank(index, bounds), d)
        every.update(quadratic_discs(coef))
        value = remainder_disc(coef)
        if value is not None:
            remainder.add(value)
    return every, remainder


def disc_sets(d, workers):
    total = 1
    for b in sig_bounds(d):
        total *= b
    if d < PARALLEL_FROM:
        every, remainder = disc_chunk((d, 0, total))
        return every, remainder, total
    edges = [total * i // CHUNKS for i in range(CHUNKS + 1)]
    jobs = [(d, edges[i], edges[i + 1]) for i in range(CHUNKS)]
    every, remainder = set(), set()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for part, rest in pool.map(disc_chunk, jobs):
            every |= part
            remainder |= rest
    return every, remainder, total


def gapless_run(found):
    out = []
    for v in range(-3, -LADDER_DEPTH - 1, -1):
        if v % 4 not in (0, 1):
            continue
        if v not in found:
            break
        out.append(v)
    return out


def polytext(coef):
    return " + ".join(
        f"{c}n^{i}" if i > 1 else (f"{c}n" if i == 1 else str(c))
        for i, c in reversed(list(enumerate(coef)))
        if c
    )


def quartic_splits(d):
    hits = []
    for f in product(*[range(b) for b in sig_bounds(d)]):
        rest = peel_unit_linear(poly_of_sig(list(f), d))
        if len(rest) != 5:
            continue
        parts = factor_list(rest)
        if len(parts) == 2 and all(g.degree() == 2 for g in parts):
            hits.append((f, rest, parts))
    return hits


def solid_cells(mask, s):
    return [
        (x, y, z)
        for x in range(s)
        for y in range(s)
        for z in range(s)
        if (mask >> ((x & 1) | ((y & 1) << 1) | ((z & 1) << 2))) & 1
    ]


def census_values(mask, s):
    cells = solid_cells(mask, s)
    solid = set(cells)
    m = s + 2
    shift = m * m * m
    exposed = 0
    vertices, edges, faces = set(), set(), set()
    for x, y, z in cells:
        for p in (
            (x - 1, y, z),
            (x + 1, y, z),
            (x, y - 1, z),
            (x, y + 1, z),
            (x, y, z - 1),
            (x, y, z + 1),
        ):
            if p not in solid:
                exposed += 1
        for a in (0, 1):
            xa = (x + a) * m
            for b in (0, 1):
                yb = (xa + y + b) * m
                for c in (0, 1):
                    vertices.add(yb + z + c)
            faces.add((x * m + y) * m + z + a)
            faces.add(shift + (x * m + y + a) * m + z)
            faces.add(2 * shift + ((x + a) * m + y) * m + z)
        for b in (0, 1):
            for c in (0, 1):
                edges.add((x * m + y + b) * m + z + c)
                edges.add(shift + ((x + b) * m + y) * m + z + c)
                edges.add(2 * shift + ((x + b) * m + y + c) * m + z)
    v, e, f, fills = len(vertices), len(edges), len(faces), len(cells)
    return [fills, s**3 - fills, exposed, v, e, f, v - e + f - fills]


def census_polys(mask, fam, scale, basis):
    side = (lambda i: 2 * i + 1) if fam == "odd" else (lambda i: 2 * i)
    samples = {i: census_values(mask, side(i)) for i in range(1, 7)}
    out = []
    bad = 0
    for obs in range(7):
        ys = [samples[1 + i][obs] for i in range(4)]
        coef = [
            Fraction(sum(ys[i] * basis[i][j] for i in range(4)), scale)
            for j in range(4)
        ]
        for h in (5, 6):
            if sum(c * h**j for j, c in enumerate(coef)) != samples[h][obs]:
                bad += 1
        out.append(coef)
    return out, bad, samples


def closed_fill(mask, k):
    return sum(
        k ** (3 - j.bit_count()) * (k - 1) ** j.bit_count()
        for j in range(8)
        if (mask >> j) & 1
    )


def factor_multisets(value, count, low=1):
    if count == 1:
        return [[value]] if value >= low else []
    out = []
    a = low
    while a**count <= value:
        if value % a == 0:
            for tail in factor_multisets(value // a, count - 1, a):
                out.append([a] + tail)
        a += 1
    return out


def divisor_shape(coef):
    if all(c == 0 for c in coef):
        return False
    k = 0
    while coef[k] == 0:
        k += 1
    q = [c / coef[k] for c in coef[k:]]
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    if len(q) == 1 or any(value.denominator != 1 for value in q):
        return False
    q = [int(value) for value in q]
    if q[-1] <= 0:
        return False
    for multiset in factor_multisets(q[-1], len(q) - 1):
        built = [1]
        for a in multiset:
            built = polymul(built, [1, a])
        if built == q:
            return True
    return False


def census_locked(mask, scale, basis):
    bad = 0
    drift = 0
    unlocked = False
    for fam in ("odd", "even"):
        polys, misses, samples = census_polys(mask, fam, scale, basis)
        bad += misses
        if fam == "odd":
            drift = sum(
                1 for i in range(1, 7) if samples[i][0] != closed_fill(mask, i + 1)
            )
        for coef in polys:
            if divisor_shape(coef):
                unlocked = True
    return not unlocked, bad, drift


def lock_predicate(mask):
    if mask & 1:
        return 0
    size = mask.bit_count()
    inner = sum(
        1
        for c in range(8)
        for i in range(3)
        if (mask >> c) & 1 and c < (c ^ (1 << i)) and (mask >> (c ^ (1 << i))) & 1
    )
    if size == 5 and inner == 4:
        return 1
    if inner == 0 and (mask >> 7) & 1 and size != 2:
        return 2
    return 0


def lock_report():
    scale, basis = scaled_basis([1, 2, 3, 4])
    path, edgeless, locked = [], [], []
    holdout = 0
    agree = 0
    grid = 0
    for mask in range(256):
        clause = lock_predicate(mask)
        if clause == 1:
            path.append(mask)
        elif clause == 2:
            edgeless.append(mask)
        is_locked, bad, drift = census_locked(mask, scale, basis)
        holdout += bad
        grid += drift == 0
        if is_locked:
            locked.append(mask)
        agree += (clause != 0) == is_locked
    return path, edgeless, sorted(path + edgeless), locked, agree, holdout, grid


def main():
    workers = min(8, os.cpu_count() or 1)
    started = time.time()
    print("DOMAIN")
    print("  every design at D = 1..4 for the polynomial sweep, 2^(2^D) of them")
    print("  every signature at D = 2..6 for the discriminants")
    print("  every design at D = 3 for the lock, all 256")
    print()
    print("FILL POLYNOMIALS")
    for d, designs, classes, distinct, closed, lower, lead, bad in polynomial_report():
        print(
            f"  D = {d}: {designs} designs, {classes} classes, {distinct} distinct polynomials, "
            f"A129824 = {closed}"
        )
        print(
            f"    lower coefficients split {lower} of {classes} classes, leading splits {lead} of {classes}, "
            f"failures {bad}"
        )
    print()
    print("QUARTIC REMAINDERS SPLITTING INTO TWO QUADRATICS AT D = 4")
    hits = quartic_splits(4)
    print(f"  count {len(hits)}")
    for f, rest, parts in hits:
        body = " ".join(f"({polytext(list(reversed(g.all_coeffs())))})" for g in parts)
        print(f"  f={tuple(f)}: {polytext(rest)} = {body}")
    print()
    print("QUADRATIC-FACTOR DISCRIMINANTS")
    for d in DISC_DIMS:
        every, remainder, total = disc_sets(d, workers)
        wide = gapless_run(every)
        narrow = gapless_run(remainder)
        print(f"  D = {d}: {total} signatures")
        print(
            f"    all quadratic factors: gapless -3..{wide[-1]}, length {len(wide)}, deepest {min(every)}"
        )
        print(
            f"    peeled remainder only: gapless -3..{narrow[-1]}, length {len(narrow)}, deepest {min(remainder)}"
        )
    print()
    print("THE CENSUS LOCK AT D = 3")
    path, edgeless, predicted, locked, agree, holdout, grid = lock_report()
    print(
        f"  path clause, origin-free with 5 corners and 4 edges: {len(path)} designs {path}"
    )
    print(
        f"  edgeless clause, origin-free with 111 and size not 2: {len(edgeless)} designs {edgeless}"
    )
    print(f"  predicate union: {len(predicted)} designs {predicted}")
    print(f"  census rule locked: {len(locked)} designs {locked}")
    print(f"  predicate and census rule agree on {agree} of 256 designs")
    print(f"  held-out census fits that missed: {holdout}")
    print(f"  grid fills match the closed form at k = 2..7 on {grid} of 256 designs")
    print()
    print(f"wall time {time.time() - started:.1f} s on {workers} workers")


if __name__ == "__main__":
    main()
