from fractions import Fraction
from math import isqrt

import mpmath
import numpy as np
import sympy as sp

N2 = 200000
N_ROT = 6000
DPS = 30
EXPONENTS = (2, 3, 4)


def parity_classes(limit):
    bound = isqrt(limit)
    keys = []
    for a in range(-bound, bound + 1):
        square = a * a
        inner = isqrt(limit - square)
        b = np.arange(-inner, inner + 1, dtype=np.int64)
        norms = square + b * b
        slot = (a & 1) * 2 + (b & 1)
        keys.append(slot * (limit + 1) + norms)
    flat = np.bincount(np.concatenate(keys), minlength=4 * (limit + 1))
    return flat.reshape(4, limit + 1).astype(np.int64)


def jacobi_r2(limit):
    out = np.zeros(limit + 1, dtype=np.int64)
    for d in range(1, limit + 1, 4):
        out[d::d] += 4
    for d in range(3, limit + 1, 4):
        out[d::d] -= 4
    return out


def shifted(source, step, limit):
    out = np.zeros(limit + 1, dtype=np.int64)
    out[::step] = source[: limit // step + 1]
    return out


def doubling_identity(cls, total, limit):
    predicted = shifted(total, 4, limit)
    return int(np.count_nonzero(cls[0][1:] != predicted[1:]))


def rotation_identity(cls, limit):
    mixed = cls[1] + cls[2]
    predicted = shifted(mixed, 2, limit)
    return int(np.count_nonzero(cls[3][1:] != predicted[1:]))


def swap_identity(cls):
    return int(np.count_nonzero(cls[1][1:] != cls[2][1:]))


def odd_odd_points(norm):
    out = []
    bound = isqrt(norm)
    for i in range(-bound, bound + 1):
        if i % 2 == 0:
            continue
        rest = norm - i * i
        if rest < 0:
            continue
        j = isqrt(rest)
        if j * j != rest or j % 2 == 0:
            continue
        out.append((i, j))
        out.append((i, -j))
    return out


def rotation_as_a_map(cls, limit):
    moved = 0
    faults = 0
    for norm in range(2, limit + 1, 2):
        points = odd_odd_points(norm)
        images = set()
        for i, j in points:
            u = (i + j) // 2
            v = (i - j) // 2
            if u * u + v * v != norm // 2 or (u & 1) == (v & 1):
                faults += 1
            images.add((u, v))
        moved += len(points)
        if len(images) != len(points) or len(points) != int(cls[3][norm]):
            faults += 1
        half = norm // 2
        if len(images) != int(cls[1][half]) + int(cls[2][half]):
            faults += 1
    return moved, faults


def pinned_polynomial():
    t, whole = sp.symbols("t S")
    ee, eo, oe, oo = sp.symbols("S_ee S_eo S_oe S_oo")
    solution = sp.solve(
        [
            ee - t**2 * whole,
            oo - t * (eo + oe),
            eo - oe,
            ee + eo + oe + oo - whole,
        ],
        [ee, eo, oe, oo],
        dict=True,
    )[0]
    a_ee, a_eo, a_oe, a_oo = sp.symbols("a_ee a_eo a_oe a_oo")
    filled = (
        a_ee * solution[ee]
        + a_eo * solution[eo]
        + a_oe * solution[oe]
        + a_oo * solution[oo]
    )
    claimed = whole * (a_ee * t**2 + a_oo * t * (1 - t) + (a_eo + a_oe) * (1 - t) / 2)
    return sp.simplify(filled - claimed) == 0


def dirichlet_beta(s):
    return (
        mpmath.zeta(s, mpmath.mpf(1) / 4) - mpmath.zeta(s, mpmath.mpf(3) / 4)
    ) / 4**s


def shell_sums(cls, total, limit):
    live = np.nonzero(total[: limit + 1])[0]
    live = live[live > 0]
    totals = {s: [mpmath.mpf(0)] * 4 for s in EXPONENTS}
    columns = [cls[slot] for slot in range(4)]
    for norm in live.tolist():
        base = mpmath.mpf(norm)
        powers = {2: base * base}
        powers[3] = powers[2] * base
        powers[4] = powers[3] * base
        for slot in range(4):
            weight = int(columns[slot][norm])
            if weight == 0:
                continue
            value = mpmath.mpf(weight)
            for s in EXPONENTS:
                totals[s][slot] += value / powers[s]
    return totals


def indicators(code):
    return [(code >> k) & 1 for k in range(4)]


def q_polynomial(a, s):
    t = Fraction(1, 2**s)
    return a[0] * t * t + a[3] * t * (1 - t) + Fraction(a[1] + a[2], 2) * (1 - t)


def to_mpf(value):
    return mpmath.mpf(value.numerator) / mpmath.mpf(value.denominator)


def closed_form_gaps(totals, limit):
    zeta = {s: mpmath.zeta(s) for s in EXPONENTS}
    beta = {s: dirichlet_beta(s) for s in EXPONENTS}
    worst = {}
    for s in EXPONENTS:
        biggest = mpmath.mpf(0)
        for code in range(1, 16):
            a = indicators(code)
            predicted = 4 * zeta[s] * beta[s] * to_mpf(q_polynomial(a, s))
            measured = mpmath.fsum(totals[s][slot] for slot in range(4) if a[slot])
            tail = mpmath.mpf(0)
            if s == 2:
                tail = sum(a) * mpmath.pi / (4 * limit)
            gap = abs(measured + tail - predicted)
            biggest = max(biggest, gap)
        worst[s] = biggest
    return worst


def sci(value):
    return mpmath.nstr(value, 2, strip_zeros=False)


def main():
    mpmath.mp.dps = DPS
    cls = parity_classes(N2)
    total = cls.sum(axis=0)
    print(
        f"DOMAIN  parity classes of Z^2 to n = {N2}, rotation as a map to n = {N_ROT}"
    )
    print()
    print("JACOBI")
    jacobi = jacobi_r2(N2)
    bad = int(np.count_nonzero(total[1:] != jacobi[1:]))
    print(f"  r2(n) = 4 (d1(n) - d3(n))            mismatches {bad} of {N2}")
    print()
    print("THE THREE SHELL RELATIONS AS INTEGER IDENTITIES ON THE COEFFICIENTS")
    bad_ee = doubling_identity(cls, total, N2)
    bad_oo = rotation_identity(cls, N2)
    bad_eo = swap_identity(cls)
    print(f"  S_ee = t^2 S    r2_ee(n) = r2(n/4)    mismatches {bad_ee} of {N2}")
    print(f"  S_oo = t S_mix  r2_oo(n) = r2_mix(n/2) mismatches {bad_oo} of {N2}")
    print(f"  S_eo = S_oe     r2_eo(n) = r2_oe(n)   mismatches {bad_eo} of {N2}")
    print()
    print("THE 45 DEGREE ROTATION AS A MAP")
    moved, faults = rotation_as_a_map(cls, N_ROT)
    print(
        f"  (i,j) -> ((i+j)/2, (i-j)/2) on every odd-odd point of even norm to n = {N_ROT}"
    )
    print(f"  odd-odd points transported {moved}, faults {faults}")
    print()
    print("THE THREE RELATIONS PIN Q_c")
    print(
        f"  solving the four linear relations returns Q_c exactly: {pinned_polynomial()}"
    )
    print()
    print("SELF TEST OF THE BETA VALUES")
    print(f"  beta(2) against Catalan   {sci(abs(dirichlet_beta(2) - mpmath.catalan))}")
    print(
        f"  beta(3) against pi^3/32   {sci(abs(dirichlet_beta(3) - mpmath.pi**3 / 32))}"
    )
    print()
    print("THE CLOSED FORM AGAINST TRUNCATED LATTICE SUMS")
    totals = shell_sums(cls, total, N2)
    worst = closed_form_gaps(totals, N2)
    for s in EXPONENTS:
        note = " after the leading tail k pi / (4N)" if s == 2 else ""
        print(
            f"  s = {s}  worst gap over the 15 nonempty designs {sci(worst[s])}{note}"
        )


main()
