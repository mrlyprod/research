import json
import subprocess
import sys
import time
from fractions import Fraction
from itertools import product
from math import comb

import sympy

LO = 2
HI = 40
PRIMES = list(sympy.primerange(2, 100))

# DIGIT POLYNOMIAL


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def digit_poly(D):
    base = [0] * (2 * D - 1)
    for k in range(D):
        base[2 * k] = comb(D - 1, k)
    return polymul(base, [1, D, 1])


def digit_poly_census(D):
    out = [0] * (2 * D + 1)
    for v in product((0, 1, 2), repeat=D):
        if sum(1 for x in v if x == 1) <= 1:
            out[sum(v)] += 1
    return out


def full_carry_census(D):
    H = half_width(D)
    M = [[0] * (2 * H + 1) for _ in range(2 * H + 1)]
    for v in product((0, 1, 2), repeat=D):
        if sum(1 for x in v if x == 1) <= 1:
            for c in range(-H, H + 1):
                t = c + D - sum(v)
                if t % 3 == 0:
                    M[c + H][t // 3 + H] += 1
    return M


def coefficient(P, x):
    return P[x] if 0 <= x < len(P) else 0


def half_width(D):
    return (D - 1) // 2


# THREE BUILDS OF THE EVEN BLOCK


def census_block(D):
    P = digit_poly(D)
    n = half_width(D) + 1
    return [
        [
            coefficient(P, c + D - 3 * cp)
            + (coefficient(P, -c + D - 3 * cp) if c else 0)
            for cp in range(n)
        ]
        for c in range(n)
    ]


def value_block(D):
    P = digit_poly(D)
    n = half_width(D) + 1
    return [
        [
            coefficient(P, c + D - 3 * cp)
            + (coefficient(P, c + D + 3 * cp) if cp else 0)
            for cp in range(n)
        ]
        for c in range(n)
    ]


def full_carry(D):
    P = digit_poly(D)
    H = half_width(D)
    return [
        [coefficient(P, c + D - 3 * cp) for cp in range(-H, H + 1)]
        for c in range(-H, H + 1)
    ]


def split_blocks(D):
    H = half_width(D)
    F = full_carry(D)
    at = lambda c, cp: F[c + H][cp + H]
    even = [
        [at(c, cp) + (at(c, -cp) if cp else 0) for cp in range(H + 1)]
        for c in range(H + 1)
    ]
    odd = [[at(c, cp) - at(c, -cp) for cp in range(1, H + 1)] for c in range(1, H + 1)]
    return F, even, odd


# EXACT DETERMINANTS


def bareiss(A):
    n = len(A)
    if n == 0:
        return 1
    M = [row[:] for row in A]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            p = next((r for r in range(k + 1, n) if M[r][k]), None)
            if p is None:
                return 0
            M[k], M[p] = M[p], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def fraction_det(A):
    n = len(A)
    M = [[Fraction(x) for x in row] for row in A]
    det = Fraction(1)
    for k in range(n):
        p = next((r for r in range(k, n) if M[r][k]), None)
        if p is None:
            return 0
        if p != k:
            M[k], M[p] = M[p], M[k]
            det = -det
        det *= M[k][k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            if f:
                for j in range(k, n):
                    M[i][j] -= f * M[k][j]
    assert det.denominator == 1
    return int(det)


def sympy_det(A):
    return int(sympy.Matrix(A).det(method="bareiss"))


def valuation(x, p):
    assert x != 0
    x = abs(x)
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


# SELF TEST


def selftest():
    for D in range(2, 11):
        assert digit_poly(D) == digit_poly_census(D), (
            "digit polynomial disagrees with the digit census at D=%d" % D
        )
    print(
        "digit polynomial (1+t^2)^(D-1)(1+Dt+t^2) against the 3^D digit census, at most one middle digit: D = 2..10, exact"
    )
    for D in range(2, 10):
        assert full_carry(D) == full_carry_census(D), (
            "carry matrix disagrees with the digit-vector transitions at D=%d" % D
        )
    print(
        "carry matrix M[c, c'] = P[c + D - 3c'] against the transitions c -> (c + D - s)/3 over all digit vectors: D = 2..9, exact"
    )
    assert [list(r) for r in zip(*census_block(3))] == [[6, 6], [1, 3]] and bareiss(
        census_block(3)
    ) == 12
    assert (
        bareiss(census_block(4)) == -66
        and census_block(4)[0][0] + census_block(4)[1][1] == 11
    )
    print(
        "controls: D = 3 block is the transpose of [[6, 6], [1, 3]] with trace 9 and det 12, D = 4 trace 11 and det -66 (char poly x^2 - 11x - 66)"
    )
    for D in range(LO, HI + 1):
        A, B = census_block(D), value_block(D)
        F, E, O = split_blocks(D)
        assert E == B, (
            "full-matrix projection disagrees with the column-fold build at D=%d" % D
        )
        d = bareiss(A)
        assert d != 0, "det M_even is zero at D=%d" % D
        assert fraction_det(A) == d and bareiss(B) == d and fraction_det(E) == d, (
            "even determinants disagree at D=%d" % D
        )
        if D <= 30:
            assert sympy_det(A) == d, "sympy disagrees at D=%d" % D
        assert bareiss(F) == d * bareiss(O), (
            "det of the full carry matrix is not det(even) det(odd) at D=%d" % D
        )
    print(
        "det M_even, D = 2..40: census fold, value fold and full-matrix projection agree through Bareiss and Fraction elimination; sympy agrees to D = 30; det(full) = det(even) det(odd) at every D; no zero"
    )
    n_of = lambda D: half_width(D) + 1
    v2 = {D: valuation(bareiss(census_block(D)), 2) for D in range(5, HI + 1, 2)}
    over = [D for D, v in v2.items() if v > n_of(D)]
    equal = [D for D, v in v2.items() if v == n_of(D)]
    assert over == [5, 7] and equal == [9, 15], (
        "v_2 <= n at odd D: violations %s, equality %s" % (over, equal)
    )
    print(
        "reproduces the Smith-cascade row at odd D = 5..39: v_2(det M_even) <= n except D = 5, 7, equality at D = 9, 15 only"
    )
    d7 = bareiss(census_block(7))
    assert valuation(d7, 2) == 7
    print(
        "reproduces the ladder row: v_2(det M_even) = %d at D = 7, det = %d"
        % (valuation(d7, 2), d7)
    )


# TABLE


def determinants():
    return {D: bareiss(census_block(D)) for D in range(LO, HI + 1)}


def verb_det():
    t = time.time()
    dets = determinants()
    print("DET: v_p(det M_even), base 3, middle-digit design, D = %d..%d" % (LO, HI))
    print(
        "%4s %3s %6s %4s %4s %4s %4s" % ("D", "n", "digits", "v_2", "v_3", "v_5", "v_7")
    )
    for D, d in dets.items():
        print(
            "%4d %3d %6d %4d %4d %4d %4d"
            % (
                D,
                half_width(D) + 1,
                len(str(abs(d))),
                *(valuation(d, p) for p in (2, 3, 5, 7)),
            )
        )
    for p in (2, 3, 5, 7):
        vs = {D: valuation(d, p) for D, d in dets.items()}
        pos = [D for D, v in vs.items() if v > 0]
        zero = [D for D, v in vs.items() if v == 0]
        top = max(vs.values())
        print(
            "p = %d: first D with v_p > 0 is %s; positive at %d of %d; max %d at D = %d; zero at D = %s"
            % (
                p,
                pos[0] if pos else None,
                len(pos),
                len(vs),
                top,
                max(D for D, v in vs.items() if v == top),
                zero,
            )
        )
    v3 = {D: valuation(d, 3) for D, d in dets.items()}
    ratios = [v3[D] / D**2 for D in range(20, HI + 1)]
    print(
        "v_3 / D^2 over D = 20..%d lies in [%.4f, %.4f] against 1/9 = %.4f; at D = 40 it is %.4f"
        % (HI, min(ratios), max(ratios), 1 / 9, v3[HI] / HI**2)
    )
    print(
        "v_3 weakly increasing over D = 3..%d: %s"
        % (HI, all(v3[D] <= v3[D + 1] for D in range(3, HI)))
    )
    print(
        "v_3(3k) = v_3(3k + 1) for every k = 1..%d: %s"
        % (
            (HI - 1) // 3,
            all(v3[3 * k] == v3[3 * k + 1] for k in range(1, (HI - 1) // 3 + 1)),
        )
    )
    print(
        "v_3(3k) - k^2 for k = 1..%d: %s"
        % (HI // 3, [v3[3 * k] - k * k for k in range(1, HI // 3 + 1)])
    )
    print("v_3 at D = 3k: %s" % [v3[3 * k] for k in range(1, HI // 3 + 1)])
    print("seconds: %.2f" % (time.time() - t))


def verb_primes():
    dets = determinants()
    print(
        "PRIMES: v_p(det M_even) > 0 for each prime p < 100 over D = %d..%d" % (LO, HI)
    )
    alive = []
    for p in PRIMES:
        vs = {D: valuation(d, p) for D, d in dets.items()}
        pos = [D for D, v in vs.items() if v > 0]
        if not pos:
            alive.append(p)
            print("p = %2d: v_p = 0 at every D" % p)
        else:
            print(
                "p = %2d: v_p > 0 at %2d of %d, first D = %2d (v_p = %d), max v_p = %d"
                % (p, len(pos), len(vs), pos[0], vs[pos[0]], max(vs.values()))
            )
    print(
        "primes p <= 50 keeping v_p = 0 at every D = %d..%d: %s"
        % (LO, HI, [p for p in alive if p <= 50] or "none")
    )
    print(
        "primes p < 100 keeping v_p = 0 at every D = %d..%d: %s"
        % (LO, HI, alive or "none")
    )
    print(
        "random-residue expectation of the count at p = 5, 7, 47: %.1f, %.1f, %.1f of %d"
        % (len(dets) / 5, len(dets) / 7, len(dets) / 47, len(dets))
    )


# FULL FACTORISATION THROUGH PARI


def gp_factor(dets):
    script = "".join(
        'f=factor(%d);print(%d," ",f[,1]~," ",f[,2]~);\n' % (abs(d), D)
        for D, d in dets.items()
    )
    out = subprocess.run(
        ["gp", "-q"],
        input=script + "quit\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    rows = {}
    for line in out.strip().splitlines():
        D, ps, es = (
            line.split(" ", 2)[0],
            *line.split(" ", 1)[1].replace("] [", "]|[").split("|"),
        )
        rows[int(D)] = dict(zip(json.loads(ps), json.loads(es)))
    return rows


def verb_factor():
    dets = determinants()
    rows = gp_factor(dets)
    print(
        "FACTOR: det M_even completely factored by PARI, D = %d..%d; prime factors above 50 listed"
        % (LO, HI)
    )
    for D, d in dets.items():
        f = rows[D]
        assert abs(d) == sympy.prod([p**e for p, e in f.items()]), (
            "factorisation does not multiply back at D=%d" % D
        )
        assert all(valuation(d, p) == f.get(p, 0) for p in PRIMES), (
            "valuation disagrees with the factorisation at D=%d" % D
        )
        big = sorted((p, e) for p, e in f.items() if p > 50)
        print(
            "D = %2d: largest prime %d (%d digits), above 50: %s"
            % (
                D,
                max(f),
                len(str(max(f))),
                " ".join("%d" % p if e == 1 else "%d^%d" % (p, e) for p, e in big)
                or "-",
            )
        )
    top = max((max(f), D) for D, f in rows.items())
    print("largest prime factor over D = %d..%d: %d at D = %d" % (LO, HI, *top))
    count = sum(1 for f in rows.values() if max(f) > 50)
    print("D with a prime factor above 50: %d of %d" % (count, len(rows)))
    seen = set().union(*rows.values())
    below = list(sympy.primerange(2, 400))
    missing = [q for q in below if q not in seen]
    print(
        "primes below 400 dividing some det: %d of %d; the least ten dividing none: %s"
        % (len(below) - len(missing), len(below), missing[:10])
    )
    print(
        "dims with q | det for q = 53, 59, 61, 67, 71: %s"
        % [sum(1 for f in rows.values() if q in f) for q in (53, 59, 61, 67, 71)]
    )


def main():
    verb = sys.argv[1] if len(sys.argv) > 1 else "all"
    verbs = {
        "selftest": selftest,
        "det": verb_det,
        "primes": verb_primes,
        "factor": verb_factor,
    }
    if verb == "all":
        for f in (selftest, verb_det, verb_primes):
            f()
    else:
        verbs[verb]()


if __name__ == "__main__":
    main()
