from fractions import Fraction
from itertools import combinations, product

import sympy as sp

DIM = 2
BASES = (2, 3)
TOP = 24


def wallis(top):
    even = {2: Fraction(16, 3)}
    odd = {3: Fraction(3, 8)}
    e, o = even[2], odd[3]
    k = 1
    while 2 * k + 3 <= top + 1:
        e *= Fraction(4 * k * (k + 1), (2 * k + 1) * (2 * k + 3))
        o *= Fraction((2 * k + 1) * (2 * k + 3), (2 * k + 2) * (2 * k + 4))
        if 2 * k + 2 <= top:
            even[2 * k + 2] = e
        if 2 * k + 3 <= top:
            odd[2 * k + 3] = o
        k += 1
    return even, odd


def squarefree_divisors(q):
    primes = sp.primefactors(q)
    out = []
    for size in range(len(primes) + 1):
        for pick in combinations(primes, size):
            value = 1
            for p in pick:
                value *= p
            out.append((value, (-1) ** size))
    return out


def bracket(design, q):
    total = len(design)
    acc = Fraction(0)
    for e, mu in squarefree_divisors(q):
        count = sum(1 for cell in design if all(c % e == 0 for c in cell))
        acc += mu * Fraction(count, total)
    return acc


def base_factor(q, dim):
    factor = Fraction(1)
    for p in sp.primefactors(q):
        factor /= 1 - Fraction(1, p**dim)
    return factor


def designs(q, dim):
    cells = list(product(range(q), repeat=dim))
    for mask in range(1, 1 << len(cells)):
        picked = [cells[i] for i in range(len(cells)) if mask >> i & 1]
        if len(picked) >= 2:
            yield tuple(picked)


def main():
    even, odd = wallis(TOP)
    print("MISMATCH: design coprime densities against the Version L family")
    print("  bracket B(F) = Sum_{e | rad(q)} mu(e) k_e / k")
    print("  delta = B(F) (1/zeta(D)) Prod_{p | q} (1 - p^-D)^-1")
    print(
        f"  at D = {DIM} that is delta = numerator/Pi^2 with numerator = 6 B(F) factor"
    )
    print(f"  Version L swept at d = 2..{TOP}")

    numerators = {}
    for q in BASES:
        factor = base_factor(q, DIM)
        seen = {}
        count = 0
        for design in designs(q, DIM):
            count += 1
            num = 6 * bracket(design, q) * factor
            seen.setdefault(num, 0)
            seen[num] += 1
        numerators[q] = seen
        shown = ", ".join(str(n) for n in sorted(seen))
        print(f"  base {q}: {count} designs, {len(seen)} distinct numerators: {shown}")

    print("  even d: numerator against numerator")
    hits = []
    for q in BASES:
        for d, r in sorted(even.items()):
            if r in numerators[q]:
                hits.append((q, d, r, numerators[q][r]))
    for q, d, r, mult in hits:
        print(f"  match: base {q}, d = {d}, numerator {r}, carried by {mult} designs")
    print(f"  total even-d matches over d = 2..{TOP}: {len(hits)}")

    zero_brackets = sum(1 for q in BASES for n in numerators[q] if n == 0)
    print(f"  designs of zero bracket in these bases: {zero_brackets}")
    print("  odd d: Version L is a nonzero rational and every density here is")
    print("  a nonzero rational over Pi^2, so a match would make Pi^2 rational")
    print(f"  odd Version L values swept: {sorted(odd)}")

    print("  the one match in closed form")
    for q, d, r, _ in hits:
        print(
            f"  base {q}, d = {d}: ({r})/Pi^2 = {sp.N(sp.Rational(r.numerator, r.denominator) / sp.pi**2, 16)}"
        )


main()
