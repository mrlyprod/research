import math
import sys
import time
from fractions import Fraction
from functools import cache
from itertools import combinations, pairwise

import mpmath as mp
import sympy as sp

CODE = (1, 3, 0, 0)
R = sp.Rational
N = sp.symbols("n")
B = sp.symbols("b", positive=True)
FIVE = {"g": (0, 0), "g+b": (0, -1), "g+1": (1, 0), "g+b-1": (-1, -1), "g+b+1": (1, -1)}

# DIGIT POLYNOMIALS


def weights(code):
    w = [0, 0, 0, 0]
    for p in range(8):
        if code >> p & 1:
            w[p.bit_count()] += 1
    return tuple(w)


def convolve(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                if y:
                    out[i + j] += x * y
    return out


BASIS = {}


def digit_polynomial(b, w):
    if b not in BASIS:
        e = [1 - k % 2 for k in range(b)]
        o = [k % 2 for k in range(b)]
        ee, eo, oo = convolve(e, e), convolve(e, o), convolve(o, o)
        BASIS[b] = [convolve(ee, e), convolve(ee, o), convolve(eo, o), convolve(oo, o)]
    return [sum(w[k] * BASIS[b][k][s] for k in range(4)) for s in range(3 * b - 2)]


def at(p, s):
    return p[s] if 0 <= s < len(p) else 0


def automaton(b, w=CODE):
    p, g = digit_polynomial(b, w), 3 * (b - 1) // 2
    return [[at(p, c + g - b * d) for d in (-1, 0, 1)] for c in (-1, 0, 1)]


def block(b, w=CODE):
    p, g = digit_polynomial(b, w), 3 * (b - 1) // 2
    return [
        [at(p, g), at(p, g + b)],
        [2 * at(p, g + 1), at(p, g + b - 1) + at(p, g + b + 1)],
    ]


def fill(b, w):
    n = (b + 1) // 2
    return sum(w[k] * n ** (3 - k) * (n - 1) ** k for k in range(4))


def spectra_block(b):
    if b % 4 == 3:
        return [
            [3 * (b + 1) * (3 * b - 1) // 16, (b + 1) * (b + 5) // 32],
            [3 * (b + 1) ** 2 // 8, 3 * (b + 1) ** 2 // 16],
        ]
    return [
        [(3 * b * b + 6 * b + 7) // 16, 3 * (b - 1) * (b + 3) // 32],
        [3 * (b - 1) * (3 * b + 5) // 8, (b + 3) ** 2 // 16],
    ]


def matmul(a, c):
    return [
        [sum(a[i][k] * c[k][j] for k in range(len(c))) for j in range(len(c[0]))]
        for i in range(len(a))
    ]


def word_count(word, letter, centre):
    m = None
    for b in reversed(word):
        x = letter(b)
        m = x if m is None else matmul(m, x)
    return m[centre][centre]


def brute(word, code):
    side = math.prod(word)
    target = 3 * (side - 1) // 2
    keep = [code >> p & 1 for p in range(8)]

    def ok(x, y, z):
        for b in reversed(word):
            if not keep[(x % b & 1) | (y % b & 1) << 1 | (z % b & 1) << 2]:
                return False
            x, y, z = x // b, y // b, z // b
        return True

    return sum(
        ok(x, y, target - x - y)
        for x in range(side)
        for y in range(max(0, target - x - side + 1), min(side, target - x + 1))
    )


# THE WORDS


def verb_words():
    print(
        "WORDS: brute-force slice counts against the carry product, finest letter on the left"
    )
    expect = {(3, 5): 60, (5, 3): 72, (3, 5, 7): 2412, (7, 5, 3): 2688, (5, 7): 300}
    bad = []
    for word in [
        (3, 5),
        (5, 3),
        (3, 5, 7),
        (7, 5, 3),
        (5, 7),
        (7, 5),
        (3, 3, 5),
        (5, 3, 3),
        (3, 7, 5),
        (5, 7, 9),
        (9, 7, 5),
        (3, 3, 3),
    ]:
        x = brute(word, 23)
        y = word_count(word, automaton, 1)
        z = word_count(word, block, 0)
        if not x == y == z or expect.get(word, x) != x:
            bad.append(word)
        print("word=%s brute=%d automaton=%d block=%d" % (word, x, y, z))
    for code in (105, 150, 126, 127, 7, 129, 255):
        for word in [(3, 5), (5, 3), (3, 5, 7), (7, 5, 3)]:
            x, y = (
                brute(word, code),
                word_count(word, lambda b: automaton(b, weights(code)), 1),
            )
            if x != y:
                bad.append((code, word))
        print(
            "code=%d words (3,5) (5,3) (3,5,7) (7,5,3): %s"
            % (code, [brute(wd, code) for wd in [(3, 5), (5, 3), (3, 5, 7), (7, 5, 3)]])
        )
    leak = [
        b
        for b in range(3, 42, 2)
        if any(
            at(digit_polynomial(b, (1, 3, 3, 1)), c + 3 * (b - 1) // 2 - b * d)
            for c in (-1, 0, 1)
            for d in (-2, 2)
        )
    ]
    print(
        "no carry leaves abs(c) <= 1 for the full cube at odd b = 3..41: %s"
        % (not leak)
    )
    print("mismatches=%s" % bad)


# THE LETTER, SYMBOLIC


@cache
def symbolic_coefficient(k, c, d, parity):
    if (c + d + parity + 1 - k) % 2:
        return sp.Integer(0), 2
    t = (c + 3 * (N - 1) - (2 * N - 1) * d - k) / 2
    total, start = sp.Integer(0), 2
    for size in range(4):
        for sub in combinations(range(3), size):
            x = sp.expand(t - sum(N - (1 if i < k else 0) for i in sub))
            poly = sp.Poly(x, N)
            alpha, beta = poly.coeff_monomial(N), poly.coeff_monomial(1)
            if alpha > 0:
                total += (-1) ** size * (x + 2) * (x + 1) / 2
                start = max(start, int(sp.ceiling(-beta / alpha)))
            else:
                start = max(start, int(sp.floor(-beta / alpha)) + 1)
    return sp.expand(total), start


def symbolic_block(w, parity):
    coef, start = {}, 2
    for name, (c, d) in FIVE.items():
        total = sp.Integer(0)
        for k in range(4):
            if w[k]:
                x, s0 = symbolic_coefficient(k, c, d, parity)
                total += w[k] * x
                start = max(start, s0)
        coef[name] = sp.expand(total)
    s = sp.Matrix(
        [[coef["g"], coef["g+b"]], [2 * coef["g+1"], coef["g+b-1"] + coef["g+b+1"]]]
    )
    return s, start


def symbolic_fill(w):
    return sp.expand(sum(w[k] * N ** (3 - k) * (N - 1) ** k for k in range(4)))


def in_b(expr):
    return sp.expand(expr.subs(N, (B + 1) / 2))


def parity_masses(w, parity):
    even, odd = R(w[0] + w[2], 4), R(w[1] + w[3], 4)
    return (odd, even) if parity == 0 else (even, odd)


def limit_letter(u, v):
    return sp.Matrix([[3 * u / 4, v / 8], [3 * v / 2, u / 4]])


def perron(m):
    lam = max(m.eigenvals(), key=lambda e: sp.N(e))
    right = (m - lam * sp.eye(2)).nullspace()[0]
    left = (m.T - lam * sp.eye(2)).nullspace()[0]
    return sp.radsimp(lam), right, left


def first_order(a, r):
    lam, right, left = perron(a)
    return sp.radsimp(
        sp.simplify((left.T * r * right)[0] / (lam * (left.T * right)[0]))
    ), lam


CLASS = {0: "b = 3 mod 4", 1: "b = 1 mod 4"}


def expansion(parity):
    s, start = symbolic_block(CODE, parity)
    a = s.applyfunc(lambda e: sp.Poly(in_b(e), B).coeff_monomial(B**2))
    first = s.applyfunc(lambda e: sp.Poly(in_b(e), B).coeff_monomial(B))
    second = s.applyfunc(lambda e: sp.Poly(in_b(e), B).coeff_monomial(1))
    return s, start, a, first, second


def verb_letter():
    print(
        "LETTER: the carry block of one letter over b^2 as an exact quadratic in t = 1/b, and its limit"
    )
    for parity in (0, 1):
        s, start, a, first, second = expansion(parity)
        bad = []
        for b in range(3, 102, 2):
            n = (b + 1) // 2
            if n % 2 != parity:
                continue
            if n >= start and sp.Matrix(s.subs(N, n)) != sp.Matrix(block(b)):
                bad.append(b)
            if sp.Matrix(spectra_block(b)) != sp.Matrix(block(b)):
                bad.append(("spectra", b))
        print(
            "%s: S_b/b^2 = A + B t + B_2 t^2 = %s + t %s + t^2 %s, extraction stable from b >= %d, mismatches at odd b = 3..101: %s"
            % (
                CLASS[parity],
                a.tolist(),
                first.tolist(),
                second.tolist(),
                2 * start - 1,
                bad,
            )
        )
        u, w = parity_masses(CODE, parity)
        q = {0: R(1, 4), 1: R(3, 4)}
        f = {-1: R(1, 8), 0: R(3, 4), 1: R(1, 8)}
        three = all(
            sp.Poly(
                in_b(
                    sum(
                        CODE[k] * symbolic_coefficient(k, c, d, parity)[0]
                        for k in range(4)
                    )
                ),
                B,
            ).coeff_monomial(B**2)
            == f[d] * q[(c + d + 1 - parity) % 2]
            for c in (-1, 0, 1)
            for d in (-1, 0, 1)
        )
        fb = sp.Poly(in_b(symbolic_fill(CODE)), B)
        print(
            "  M_b[c, c']/b^2 -> f(c') q(c + c' + g) on all nine carries, f = 1/8, 3/4, 1/8, q = 1/4 even, 3/4 odd: %s; fill/b^3 = %s (1 + (%s)/b + O(b^-2))"
            % (
                three,
                fb.coeff_monomial(B**3),
                fb.coeff_monomial(B**2) / fb.coeff_monomial(B**3),
            )
        )
        lam, _, _ = perron(a)
        print(
            "  limit = [[3u/4, v/8], [3v/2, u/4]] with u = %s, v = %s: %s; trace %s, det %s, lambda = %s = %.6f"
            % (u, w, a == limit_letter(u, w), a.trace(), a.det(), lam, float(lam))
        )
        t = in_b(s.trace())
        dt = in_b(s.det())
        root = sp.factor(t / 2) + sp.sqrt(sp.factor(t**2 / 4 - dt))
        print(
            "  exact root rho_b = %s, disc/16 = %s"
            % (
                sp.simplify(root).subs(sp.Abs(B - 1), B - 1),
                sp.factor(t**2 - 4 * dt) / 16,
            )
        )
        mu, _ = first_order(a, first)
        print("  rho_b/b^2 = lambda (1 + mu/b + O(b^-2)), mu = %s" % mu)
        print(
            "  (log b)(log_b rho_b - log_b fill + 1) = log(2 lambda) + (mu - 3/2)/b + O(b^-2) = %.6f + (%s)/b, sign %s"
            % (
                math.log(2 * float(lam)),
                sp.radsimp(mu - R(3, 2)),
                "+" if 2 * lam > 1 else "-",
            )
        )
    a3, a1 = expansion(0)[2], expansion(1)[2]
    lam, _, _ = perron(a1 * a3)
    print(
        "pair: lambda(A_1 A_3) = %s, per letter sqrt = %.6f; [A_3, A_1] = %s"
        % (lam, math.sqrt(float(lam)), (a3 * a1 - a1 * a3).tolist())
    )


# THE DRIFT


def ink_logs(sides, marks, rows):
    v, scale, side2, out = [1.0, 0.0], 0.0, 0.0, {}
    for i, b in enumerate(sides, 1):
        s = spectra_block(b)
        x = [[s[0][0] / b**2, s[0][1] / b**2], [s[1][0] / b**2, s[1][1] / b**2]]
        if rows:
            v = [v[0] * x[0][0] + v[1] * x[1][0], v[0] * x[0][1] + v[1] * x[1][1]]
        else:
            v = [x[0][0] * v[0] + x[0][1] * v[1], x[1][0] * v[0] + x[1][1] * v[1]]
        top = max(v)
        v = [y / top for y in v]
        scale += math.log(top)
        side2 += 2 * math.log(b)
        if i in marks:
            out[i] = scale + math.log(v[0]) - math.log(0.75 + 0.25 * math.exp(-side2))
    return out


WORDS = {
    "sides 3, 7, 11, .. coarsest first": (lambda k: 4 * k - 1, False),
    "sides 5, 9, 13, .. coarsest first": (lambda k: 4 * k + 1, False),
    "sides 3, 5, 7, .. coarsest first": (lambda k: 2 * k + 1, False),
    "sides 3, 5, 7, .. finest first": (lambda k: 2 * k + 1, True),
}


def exponents():
    a3, b3 = expansion(0)[2:4]
    a1, b1 = expansion(1)[2:4]
    mu3, l3 = first_order(a3, b3)
    mu1, l1 = first_order(a1, b1)
    nu, lp = first_order(a1 * a3, b1 * a3 + a1 * b3)
    nu2, _ = first_order(a3 * a1, b3 * a1 + a3 * b1)
    return {
        "3": (l3, mu3 / 4),
        "1": (l1, mu1 / 4),
        "pair": (sp.sqrt(lp), nu / 4),
        "pair reversed": (sp.sqrt(lp), nu2 / 4),
    }, (a3, a1, b3, b1)


def verb_drift():
    print(
        "DRIFT: first-order perturbation of the limiting letter, then local fits at two lengths"
    )
    ex, (_a3, _a1, b3, b1) = exponents()
    print(
        "first-order term identical in both classes: %s = (3/8) [[1, 1/2], [2, 1]]: %s"
        % (b3.tolist(), b3 == b1 == R(3, 8) * sp.Matrix([[1, R(1, 2)], [2, 1]]))
    )
    for key, (lam, gamma) in ex.items():
        print(
            "%s: lambda = %s, gamma = %s = %.9f"
            % (key, sp.nsimplify(lam), sp.radsimp(gamma), float(gamma))
        )
    marks = [4000, 8000, 16000, 32000]
    keys = ["3", "1", "pair", "pair"]
    for (name, (side, rows)), key in zip(WORDS.items(), keys):
        lam, gamma = map(float, ex[key])
        out = ink_logs([side(k) for k in range(1, marks[-1] + 1)], set(marks), rows)
        est = [
            (out[2 * m] - out[m] - m * math.log(lam)) / math.log(2) for m in marks[:-1]
        ]
        rich = [2 * y - x for x, y in pairwise(est)]
        print(
            "%s: derived gamma %.6f, local fits L -> 2L at L = 4000, 8000, 16000: %s; Richardson at 8000, 16000: %s"
            % (
                name,
                gamma,
                ", ".join("%.6f" % e for e in est),
                ", ".join("%.6f" % e for e in rich),
            )
        )


# THE CONSTANTS


def neville(xs, ys):
    p = list(ys)
    for k in range(1, len(xs)):
        for i in range(len(xs) - 1, k - 1, -1):
            p[i] = (p[i] * xs[i - k] - p[i - 1] * xs[i]) / (xs[i - k] - xs[i])
    return p[-1]


def constants_along(side, rows, marks, lam, gamma):
    v, scale, side2, out = [mp.mpf(1), mp.mpf(0)], mp.mpf(0), mp.mpf(0), {}
    for i in range(1, max(marks) + 1):
        b = side(i)
        s = spectra_block(b)
        bb = mp.mpf(b) ** 2
        x = [[s[0][0] / bb, s[0][1] / bb], [s[1][0] / bb, s[1][1] / bb]]
        if rows:
            v = [v[0] * x[0][0] + v[1] * x[1][0], v[0] * x[0][1] + v[1] * x[1][1]]
        else:
            v = [x[0][0] * v[0] + x[0][1] * v[1], x[1][0] * v[0] + x[1][1] * v[1]]
        top = max(v)
        v = [y / top for y in v]
        scale += mp.log(top)
        side2 += 2 * mp.log(b)
        if i in marks:
            out[i] = mp.exp(
                scale
                + mp.log(v[0])
                - mp.log(mp.mpf(3) / 4 + mp.exp(-side2) / 4)
                - i * mp.log(lam)
                - gamma * mp.log(i)
            )
    return out


def blink(a3, a1):
    lam, right, _ = perron(a1 * a3)
    coarse = sp.radsimp((a3 * right)[0] / (sp.sqrt(lam) * right[0]))
    lam, _, left = perron(a3 * a1)
    fine = sp.radsimp((left.T * a3)[0] / (sp.sqrt(lam) * left[0]))
    return coarse, fine


def verb_constants():
    print(
        "CONSTANTS: ink / (lambda^L L^gamma) at L = 2^j, Neville-extrapolated in 1/L to two depths"
    )
    mp.mp.dps = 40
    ex, (a3, a1, _, _) = exponents()
    keys = ["3", "1", "pair", "pair"]
    found = {}
    for (name, (side, rows)), key in zip(WORDS.items(), keys):
        lam, gamma = (mp.mpf(sp.N(x, 50)) for x in ex[key])
        lengths = [256 * 2**j for j in range(7)]
        for shift, tag in ((0, "even L"), (1, "odd L")):
            if key != "pair" and shift:
                continue
            marks = [m + shift for m in lengths]
            out = constants_along(side, rows, set(marks), lam, gamma)
            xs = [mp.mpf(1) / m for m in marks]
            ys = [out[m] for m in marks]
            deep, shallow = neville(xs, ys), neville(xs[:-1], ys[:-1])
            digits = (
                int(-mp.log10(abs(deep - shallow) / deep)) if deep != shallow else 40
            )
            found[name, tag] = deep
            print(
                "%s, %s: C(L=%d) = %s; extrapolated %s and %s, agreeing to %d digits"
                % (
                    name,
                    tag,
                    marks[-1],
                    mp.nstr(ys[-1], 12),
                    mp.nstr(deep, 18),
                    mp.nstr(shallow, 18),
                    digits,
                )
            )
    coarse, fine = blink(a3, a1)
    for name, exact in (
        ("sides 3, 5, 7, .. coarsest first", coarse),
        ("sides 3, 5, 7, .. finest first", fine),
    ):
        ratio = found[name, "odd L"] / found[name, "even L"]
        print(
            "%s: blink C_odd/C_even = %s, derived %s = %s"
            % (
                name,
                mp.nstr(ratio, 15),
                sp.nsimplify(exact),
                mp.nstr(mp.mpf(sp.N(exact, 50)), 15),
            )
        )
    print(
        "row word, finest-first over coarsest-first: even L %s, odd L %s"
        % (
            mp.nstr(
                found["sides 3, 5, 7, .. finest first", "even L"]
                / found["sides 3, 5, 7, .. coarsest first", "even L"],
                15,
            ),
            mp.nstr(
                found["sides 3, 5, 7, .. finest first", "odd L"]
                / found["sides 3, 5, 7, .. coarsest first", "odd L"],
                15,
            ),
        )
    )


# THE DESIGNS


def exact_sign(s, q):
    if s[0][1] * s[1][0] == 0:
        if s[0][0] == 0:
            return None
        return (s[0][0] > q) - (s[0][0] < q)
    t = s[0][0] + s[1][1]
    chi = q * q - t * q + s[0][0] * s[1][1] - s[0][1] * s[1][0]
    if 2 * q < t or chi < 0:
        return 1
    return 0 if chi == 0 else -1


def last_root(p):
    p = sp.Poly(sp.expand(p), N)
    if p.is_zero or p.degree() == 0:
        return 0
    ends = [hi for (lo, hi), _ in p.intervals()]
    return int(sp.ceiling(max(ends))) + 1 if ends else 0


def eventual(w, parity):
    s, start = symbolic_block(w, parity)
    f = symbolic_fill(w)
    if s[0, 1] == 0:
        polys, guard = [(2 * N - 1) * s[0, 0] - f], s[0, 0]
    else:
        t, d = s.trace(), s.det()
        polys, guard = (
            [
                t * (2 * N - 1) - 2 * f,
                f**2 - t * f * (2 * N - 1) + d * (2 * N - 1) ** 2,
            ],
            s[0, 1] * s[1, 0],
        )
    roots = max(last_root(p) for p in polys + [guard])
    zeros = [r for r in sp.Poly(sp.expand(guard), N).real_roots()]
    top = max(start, roots)
    n = top + 2 + (top + parity) % 2
    b = 2 * n - 1
    return (
        exact_sign(
            [[int(x) for x in row] for row in s.subs(N, n).tolist()],
            Fraction(int(f.subs(N, n)), b),
        ),
        b,
        start,
        roots,
        max(zeros) if zeros else None,
    )


def verb_designs():
    print(
        "DESIGNS: all 256 dim 3 parity designs at infinite side, and the sign of the census exponent against solid dimension minus 1 at every odd base"
    )
    classes = {}
    for code in range(1, 256):
        classes.setdefault(weights(code), []).append(code)
    tally, mism, late, law, reach, flat, guards, blind = (
        {},
        [],
        [],
        0,
        [0, 0, 0],
        [],
        {},
        set(),
    )
    for w, codes in sorted(classes.items()):
        e, o = w[0] + w[2], w[1] + w[3]
        row = []
        for parity in (0, 1):
            u, v = parity_masses(w, parity)
            lam = (2 * u + sp.sqrt(u**2 + 3 * v**2)) / 4
            s, start = symbolic_block(w, parity)
            for b in range(2 * start - 1, 62, 2):
                if ((b + 1) // 2) % 2 == parity and [
                    [int(x) for x in r] for r in s.subs(N, (b + 1) // 2).tolist()
                ] != block(b, w):
                    mism.append((w, b, "form"))
            if s.applyfunc(
                lambda x: sp.Poly(in_b(x), B).coeff_monomial(B**2)
            ) != limit_letter(u, v):
                mism.append((w, parity, "limit"))
            sign, b_star, start, roots, zero = eventual(w, parity)
            if zero is not None:
                guards[w, parity] = zero
            if s[0, 0] == 0 and s[1, 1] == 0:
                blind.update(codes)
            reach = [max(reach[0], start), max(reach[1], b_star), max(reach[2], roots)]
            predicted = int(sp.sign(u - v))
            if sign != predicted or predicted != int(sp.sign(o - e)) * (1 - 2 * parity):
                mism.append((w, parity, "sign"))
            for b in range(3, b_star + 1, 2):
                if ((b + 1) // 2) % 2 == parity:
                    x = exact_sign(block(b, w), Fraction(fill(b, w), b))
                    if x != sign:
                        late.append((w, b, x, sign))
                    if x != sign and not (x is None and sign == -1):
                        law += len(codes)
            row.append((sign, lam))
        if e == o:
            for b in range(3, 42, 2):
                m = automaton(b, w)
                if any(sum(r) * b != fill(b, w) for r in m):
                    flat.append((w, b))
        key = (row[0][0], row[1][0])
        tally[key] = tally.get(key, 0) + len(codes)
        print(
            "weights=%s e=%d o=%d designs=%d: 8 lambda/(e + o) = %.6f at 3 mod 4, %.6f at 1 mod 4, sign %s"
            % (
                w,
                e,
                o,
                len(codes),
                float(8 * row[0][1] / (e + o)),
                float(8 * row[1][1] / (e + o)),
                key,
            )
        )
    print(
        "tally of signs (3 mod 4, 1 mod 4) over codes 1..255: %s"
        % dict(sorted(tally.items()))
    )
    print(
        "symbolic forms stable from b >= %d; every real root of every sign polynomial below n = %d, b = %d; every odd base up to b = %d checked exactly"
        % (2 * reach[0] - 1, reach[2], 2 * reach[2] - 1, reach[1])
    )
    print(
        "largest real root of the case guard, s00 on a triangular class and s01 s10 otherwise: %s at %s"
        % ((max(guards.values()), max(guards, key=guards.get)) if guards else "none")
    )
    print(
        "codes whose block has a zero diagonal on a class, count zero at every odd level there: %d"
        % len(blind)
    )
    print(
        "odd bases where the exact sign differs from the large-side sign: %s (None is an empty slice)"
        % late
    )
    print(
        "designs violating sign = sgn(o - e) at 3 mod 4 and sgn(e - o) at 1 mod 4, an empty slice counted below: %d"
        % law
    )
    print(
        "e = o designs with an automaton row sum other than fill/b at odd b = 3..41: %s"
        % flat
    )
    print("mismatches=%s" % mism)


# THE CLOSED FORM AT SIDES 3 MOD 4


def verb_closed():
    print(
        "CLOSED: sides 3, 7, 11, .. coarsest first, the generating function and the constant"
    )
    k, z = sp.symbols("k z")
    a = sp.Matrix([[18, 1], [12, 6]])
    c = sp.Matrix([[-6, 1], [0, 0]])
    b = 4 * k - 1
    letter = sp.Matrix(
        [
            [3 * (b + 1) * (3 * b - 1) / 16, (b + 1) * (b + 5) / 32],
            [3 * (b + 1) ** 2 / 8, 3 * (b + 1) ** 2 / 16],
        ]
    )
    print(
        "S_(4k-1) = (k/2)(k K_1 + K_0) with K_1 = %s, K_0 = %s: %s"
        % (
            a.tolist(),
            c.tolist(),
            sp.expand(letter - k / 2 * (k * a + c)) == sp.zeros(2),
        )
    )
    d = 1 - 24 * z + 96 * z**2
    p = sp.Matrix([1 - 6 * z, 12 * z])
    lhs = (sp.eye(2) - z * a) * (d * p.diff(z) - R(3, 4) * d.diff(z) * p)
    rhs = d * (a + c) * p
    print(
        "W = (1 - 24z + 96z^2)^(-3/4) (1 - 6z, 12z) solves (1 - z K_1) W' = (K_1 + K_0) W, W(0) = e_0: %s"
        % (sp.expand(lhs - rhs) == sp.zeros(2, 1))
    )
    y = [Fraction(1), Fraction(18)]
    for m in range(1, 60):
        y.append(((24 * m + 18) * y[m] - (96 * m + 48) * y[m - 1]) / (m + 1))
    v, bad = [1, 0], []
    for L in range(1, 61):
        s = spectra_block(4 * L - 1)
        v = [s[0][0] * v[0] + s[0][1] * v[1], s[1][0] * v[0] + s[1][1] * v[1]]
        series = y[L] - 6 * y[L - 1]
        if Fraction(math.factorial(L) ** 2, 2**L) * series != v[0]:
            bad.append(L)
        if L <= 4:
            print("L=%d count=%d" % (L, v[0]))
    print(
        "count_L = (L!)^2 2^(-L) [z^L] (1 - 6z)(1 - 24z + 96z^2)^(-3/4) at L = 1..60: mismatches=%s"
        % bad
    )
    mp.mp.dps = 40
    exact = (
        mp.gamma(mp.mpf(3) / 4)
        * (1 + mp.sqrt(3))
        / (3 * (mp.sqrt(3) - 1) ** (mp.mpf(3) / 4))
    )
    lam = (3 + mp.sqrt(3)) / 8
    marks = [256 * 2**j for j in range(7)]
    out = constants_along(lambda i: 4 * i - 1, False, set(marks), lam, mp.mpf(1) / 4)
    print(
        "C_3 = Gamma(3/4)(1 + sqrt(3))/(3 (sqrt(3) - 1)^(3/4)) = %s; extrapolated %s"
        % (
            mp.nstr(exact, 25),
            mp.nstr(
                neville([mp.mpf(1) / m for m in marks], [out[m] for m in marks]), 20
            ),
        )
    )


# RUN

VERBS = {
    "words": verb_words,
    "letter": verb_letter,
    "drift": verb_drift,
    "constants": verb_constants,
    "closed": verb_closed,
    "designs": verb_designs,
}

if __name__ == "__main__":
    for name in sys.argv[1:] or list(VERBS):
        start = time.time()
        VERBS[name]()
        print("%s: %.1fs" % (name, time.time() - start))
