import sys
import time
from fractions import Fraction
from math import comb

import numpy as np
from mpmath import mp, mpf

LO = 2
HI = 120
TERMS = 140
PRODUCT = "1.3461220067642173"
FIT = "1.34612251727283689"
DETS = [2, 12, -66, 3120, -14310, 3894912]
MARKS = (60, 80, 100, 120)


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


def even_block(D):
    P = digit_poly(D)
    n = (D - 1) // 2 + 1

    def at(new, old):
        i = old + D - 3 * new
        return P[i] if 0 <= i < len(P) else 0

    return [[at(k, c) + (at(k, -c) if c else 0) for c in range(n)] for k in range(n)]


def bareiss(A):
    n = len(A)
    M = [row[:] for row in A]
    sign, prev = 1, 1
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


def trunc(x, digits):
    e = int(mp.floor(mp.log10(x)))
    whole = int(mp.floor(x * mpf(10) ** (digits - 1 - e)))
    s = str(whole)
    return s[: e + 1] + "." + s[e + 1 :] if e >= 0 else "0." + "0" * (-e - 1) + s


def show(x, digits=20):
    return mp.nstr(x, digits)


def cosines(terms):
    return [mp.cos(2 * mp.pi / mpf(3) ** k) for k in range(2, terms + 1)]


def product(terms):
    return mp.fprod(cosines(terms))


def verb_product():
    start = time.time()
    seen = {}
    for dps in (60, 120):
        mp.dps = dps
        p = product(TERMS)
        seen[dps] = (trunc(p, 40), trunc(1 / p, 40))
        print(
            f"dps = {dps}, {TERMS - 1} terms: prod_(k>=2) cos(2 pi/3^k) = {seen[dps][0]}, inverse r_inf = {seen[dps][1]}"
        )
    assert seen[60] == seen[120]
    print("the two precisions agree on all 40 truncated digits of both")
    mp.dps = 120
    tail = mp.pi**2 / (4 * mpf(9) ** TERMS)
    print(
        f"tail: 1 - prod_(k>{TERMS}) cos(2 pi/3^k) <= pi^2/(4 9^{TERMS}) = {show(tail, 3)}, below the working precision"
    )
    p = product(TERMS)
    steps = [(k, abs(product(k) - product(k - 1))) for k in (5, 10, 20, 40)]
    print(
        "geometric convergence, |prod(k) - prod(k-1)|: "
        + ", ".join(f"k = {k}: {show(s, 3)}" for k, s in steps)
    )
    print(
        f"ratio of consecutive steps at k = 40: {show(steps[3][1] / abs(product(39) - product(38)), 6)} against 1/9 = {show(mpf(1) / 9, 6)}"
    )
    r = 1 / p
    print(f"r_inf = {trunc(r, 40)}")
    print(f"2 r_inf = {trunc(2 * r, 20)}")
    print(f"A = 2/(3 prod) = {trunc(2 / (3 * p), 20)}")
    print(
        f"prod - 3/4 = {show(p - mpf(3) / 4, 12)}, r_inf - 4/3 = {show(r - mpf(4) / 3, 12)}, 2 r_inf - 8/3 = {show(2 * r - mpf(8) / 3, 12)}"
    )
    print(
        f"the product candidate {PRODUCT} is r_inf rounded to 17 digits: {show(mpf(PRODUCT) - r, 3)}"
    )
    fit = mpf(FIT)
    agree = next(i for i in range(len(FIT)) if trunc(r, 30)[i] != FIT[i])
    print(
        f"the fit candidate {FIT} differs from r_inf by {show(fit - r, 6)} ({show((fit - r) / r, 3)} relative) and shares its first {agree - 1} digits"
    )
    print(f"product done in {time.time() - start:.1f} s")


def perron(D):
    A = even_block(D)
    n = len(A)
    fill = 2 ** (D - 1) * (D + 2)
    shift = mpf(fill) / 3
    B = mp.matrix(
        [[mpf(A[i][j]) - (shift if i == j else 0) for j in range(n)] for i in range(n)]
    )
    x = mp.matrix([mpf(1)] * n)
    eps = mpf(10) ** (10 - mp.dps)
    for _ in range(1000):
        y = mp.lu_solve(B, x)
        y = y / sum(y[i] for i in range(n))
        image = [sum(mpf(A[i][j]) * y[j] for j in range(n)) for i in range(n)]
        ratios = [image[i] / y[i] for i in range(n)]
        x = y
        if max(ratios) - min(ratios) < eps * max(ratios):
            break
    v = [Fraction(mp.nstr(y[i], mp.dps)) for i in range(n)]
    assert all(c > 0 for c in v)
    exact = [sum(A[i][j] * v[j] for j in range(n)) / v[i] for i in range(n)]
    return min(exact), max(exact), fill, A


def real(f):
    return mpf(f.numerator) / f.denominator


def neville(points):
    xs = [mpf(1) / D for D, _ in points]
    P = [v for _, v in points]
    n = len(points)
    for level in range(1, n):
        for i in range(n - level):
            P[i] = (xs[i + level] * P[i] - xs[i] * P[i + 1]) / (xs[i + level] - xs[i])
    return P[0]


def least_squares(points, order):
    A = mp.matrix([[mpf(1) / mpf(D) ** j for j in range(order + 1)] for D, _ in points])
    b = mp.matrix([v for _, v in points])
    return mp.lu_solve(A.T * A, A.T * b)[0]


def verb_rate():
    start = time.time()
    mp.dps = 100
    dets = [bareiss(even_block(D)) for D in range(2, 8)]
    assert dets == DETS
    print(f"det M_even at dim = 2..7: {dets}, the known cross-check")
    rho, width, fill = {}, {}, {}
    for D in range(LO, HI + 1):
        lo, hi, fill[D], A = perron(D)
        rho[D] = real(lo + hi) / 2
        width[D] = real(hi - lo)
        assert width[D] < mpf(10) ** (-70) * rho[D]
        top = max(np.linalg.eigvals(np.array(A, dtype=float)).real)
        assert abs(top - float(rho[D])) < 1e-8 * float(rho[D])
    print(
        f"rho_dim enclosed by exact rational Collatz-Wielandt bounds at every dim = {LO}..{HI}, widest relative width {show(max(width[D] / rho[D] for D in rho), 3)}, numpy's dominant eigenvalue within 1e-8 at every dim"
    )
    print(
        f"rho_3 = {trunc(rho[3], 30)} against (9 + sqrt 33)/2 = {trunc((9 + mp.sqrt(33)) / 2, 30)}"
    )
    print(f"rho_{HI} = {trunc(rho[HI], 60)}")
    print(f"fill_{HI}/3 = {trunc(mpf(fill[HI]) / 3, 60)}")
    p = product(TERMS)
    cosk = cosines(TERMS)
    rp, rf = 1 / p, mpf(FIT)
    A = 2 / (3 * p)
    delta = {D: rho[D] - mpf(fill[D]) / 3 for D in rho}
    ratio = {D: abs(delta[D] / delta[D - 1]) for D in range(LO + 1, HI + 1)}
    form = {
        D: (2 * mpf(D - 1) / 3)
        * mp.fprod(c ** (D - 1) * (D + 2 * c) / (D + 2) for c in cosk)
        for D in rho
    }
    resid = {D: abs(delta[D]) / form[D] - 1 for D in rho}
    scaled = {r: {D: abs(delta[D]) * r**D / (D - 1) for D in rho} for r in (rp, rf)}
    print(
        f"candidates: product r = {trunc(rp, 20)}, fit r = {FIT}; A = 2/(3 prod) = {trunc(A, 16)}"
    )
    print(
        "dim, delta = rho - fill/3, |delta_dim/delta_(dim-1)|, |delta| r^dim/(dim-1) under the product and under the fit, |delta|/form - 1 for form = (2(dim-1)/3) prod cos^(dim-1) (dim + 2 cos)/(dim+2)"
    )
    for D in range(LO, HI + 1):
        row = f"{D:4d} {show(delta[D], 20):>26} "
        row += f"{show(ratio[D], 20):>22} " if D > LO else " " * 23
        row += f"{show(scaled[rp][D], 16):>19} {show(scaled[rf][D], 16):>19} {show(resid[D], 4):>12}"
        print(row)
    print(
        f"the form residual alternates and falls geometrically: {show(resid[60], 4)} at dim = 60, {show(resid[120], 4)} at dim = 120, consecutive ratio {show(resid[120] / resid[119], 6)}"
    )
    print(
        "Richardson in 1/dim of the one-step ratio, Neville on the last order+1 dims, orders 1..4"
    )
    for D in MARKS:
        cells = [
            show(neville([(d, ratio[d]) for d in range(D - k, D + 1)]), 16)
            for k in range(1, 5)
        ]
        print(f"  dim = {D}: " + ", ".join(cells))
    print("same, on dims of one parity (step 2), orders 1..5")
    for D in (100, 120):
        cells = [
            show(neville([(d, ratio[d]) for d in range(D - 2 * k, D + 1, 2)]), 16)
            for k in range(1, 6)
        ]
        print(f"  dim = {D}: " + ", ".join(cells))
    print(
        "least-squares fit c + a_1/dim + .. + a_order/dim^order of the one-step ratio, orders 1..3"
    )
    for hi in (100, 120):
        cells = [
            show(least_squares([(d, ratio[d]) for d in range(60, hi + 1)], k), 16)
            for k in range(1, 4)
        ]
        print(f"  dim = 60..{hi}: " + ", ".join(cells))
    for parity, name in ((0, "even"), (1, "odd")):
        split = least_squares(
            [(d, ratio[d]) for d in range(60, 121) if d % 2 == parity], 3
        )
        print(f"  order 3 on {name} dims of 60..120: {show(split, 16)}")
    fit3 = least_squares([(d, ratio[d]) for d in range(60, 121)], 3)
    print(
        f"the order-3 fit over dim = 60..120 reads {show(fit3, 18)}, inverse {show(1 / fit3, 18)}: {show(p - fit3, 3)} below the product, the size of the order-3 truncation at dim ~ 100"
    )
    print(
        f"the retired fit {FIT} sits {show(p - 1 / rf, 3)} below the product on the ratio scale"
    )
    print(
        "the ratio deflated by the form's own dim-dependence, (dim-1)/(dim-2) prod (dim + 2 cos)(dim+1)/((dim+2)(dim - 1 + 2 cos)), against the product"
    )
    for D in MARKS:
        F = (
            mpf(D - 1)
            / (D - 2)
            * mp.fprod(
                (D + 2 * c) * (D + 1) / ((D + 2) * (D - 1 + 2 * c)) for c in cosk
            )
        )
        d = ratio[D] / F
        print(
            f"  dim = {D}: {trunc(d, 24)}, relative gap to the product {show(d / p - 1, 3)}"
        )
    print(
        f"Richardson of the scaled excess at dim = {HI}, orders 1..3, against A = {trunc(A, 16)}"
    )
    for name, r in (("product", rp), ("fit", rf)):
        cells = [
            neville([(d, scaled[r][d]) for d in range(HI - k, HI + 1)])
            for k in range(1, 4)
        ]
        print(
            f"  {name}: "
            + ", ".join(show(c, 16) for c in cells)
            + ", gaps "
            + ", ".join(show(c - A, 3) for c in cells)
        )
    for r, name in ((rp, "product"), (rf, "fit")):
        print(
            f"|delta| r^dim under the {name}: {show(abs(delta[60]) * r**60, 12)} at dim = 60, {show(abs(delta[100]) * r**100, 12)} at dim = 100"
        )
    print(
        "verdict: the one-step ratio and the scaled excess converge to the product; the fit's digits past the seventh are the order-3 truncation residue"
    )
    print(f"rate done in {time.time() - start:.1f} s")


def main():
    verbs = {"product": verb_product, "rate": verb_rate}
    chosen = sys.argv[1:] or list(verbs)
    for name in chosen:
        verbs[name]()


if __name__ == "__main__":
    main()
