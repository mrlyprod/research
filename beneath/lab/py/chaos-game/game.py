import math
import random
import sys
from fractions import Fraction
from itertools import combinations, product

import numpy as np

SEED = 20261030


def corners(dim, base=2):
    return list(product(range(base), repeat=dim))


def step(x, a, base):
    return tuple((xi + ai) / Fraction(base) for xi, ai in zip(x, a))


def play(x, draws, base):
    for a in draws:
        x = step(x, a, base)
    return x


def cell(x, k, base):
    return tuple(math.floor(xi * base**k) for xi in x)


def word_cell(word, base):
    k = len(word)
    return tuple(
        sum(word[j][i] * base ** (k - 1 - j) for j in range(k))
        for i in range(len(word[0]))
    )


def starts(dim, base):
    pts = [
        tuple(Fraction(0) for _ in range(dim)),
        tuple(Fraction(base**3 - 1, base**3) for _ in range(dim)),
        tuple(Fraction(1, 3 + i) for i in range(dim)),
        tuple(Fraction(2 * i + 5, 7 * i + 13) for i in range(dim)),
    ]
    return pts


def tv(p, q):
    keys = set(p) | set(q)
    return sum(abs(p.get(c, 0) - q.get(c, 0)) for c in keys) / 2


def register():
    bad, checked = 0, 0
    for base, dim in [(2, 1), (2, 2), (3, 2), (2, 3)]:
        alphabet = corners(dim, base)
        for x0 in starts(dim, base):
            for n in range(5):
                for draws in product(alphabet, repeat=n):
                    x = play(x0, draws, base)
                    closed = tuple(
                        sum(
                            draws[n - j][i] * Fraction(1, base**j)
                            for j in range(1, n + 1)
                        )
                        + x0[i] / base**n
                        for i in range(dim)
                    )
                    checked += 1
                    bad += x != closed
                    for k in range(1, n + 1):
                        rev = tuple(reversed(draws[n - k :]))
                        checked += 1
                        bad += cell(x, k, base) != word_cell(rev, base)
    print(
        f"register: x_n = sum a_(n+1-j) b^-j + b^-n x_0 and cell_k(x_n) = reversed last k draws, {checked} checks, {bad} mismatches"
    )
    assert bad == 0


def weights(F, uniform):
    if uniform:
        return {a: Fraction(1, len(F)) for a in F}
    total = sum(range(1, len(F) + 1))
    return {a: Fraction(i + 1, total) for i, a in enumerate(F)}


def cell_law(x0, F, p, n, k, base):
    law = {}
    for draws in product(F, repeat=n):
        mass = math.prod((p[a] for a in draws), start=Fraction(1))
        c = cell(play(x0, draws, base), k, base)
        law[c] = law.get(c, 0) + mass
    return law


def target(F, p, k, base):
    law = {}
    for word in product(F, repeat=k):
        c = word_cell(word, base)
        law[c] = law.get(c, 0) + math.prod((p[a] for a in word), start=Fraction(1))
    return law


def designs():
    out = []
    for base, dim, kmax in [(2, 1, 4), (2, 2, 3), (3, 1, 3)]:
        alphabet = corners(dim, base)
        for r in range(1, len(alphabet) + 1):
            for F in combinations(alphabet, r):
                out.append((base, dim, list(F), kmax))
    carpet = [a for a in corners(2, 3) if a != (1, 1)]
    out.append((3, 2, carpet, 3))
    flake = [a for a in corners(3, 2) if 0 < sum(a) < 3]
    out.append((2, 3, flake, 2))
    return out


def forget():
    cases = late = early = early_uniform = 0
    exact_late = positive_early = bound_early = 0
    worst = Fraction(0)
    for base, dim, F, kmax in designs():
        for uniform in (True, False):
            p = weights(F, uniform)
            for k in range(1, kmax + 1):
                want = target(F, p, k, base)
                for x0 in starts(dim, base):
                    for n in range(k + 2):
                        if len(F) ** n > 4096:
                            continue
                        d = tv(cell_law(x0, F, p, n, k, base), want)
                        cases += 1
                        if n >= k:
                            late += 1
                            exact_late += d == 0
                            worst = max(worst, d)
                        elif len(F) >= 2:
                            early += 1
                            positive_early += d > 0
                            if uniform:
                                early_uniform += 1
                                bound_early += d >= 1 - Fraction(
                                    len(F) ** n, len(F) ** k
                                )
    print(
        f"forget: {len(designs())} designs, uniform and weighted, {cases} (design, weights, k, start, n) cases by exact enumeration of every draw word"
    )
    print(
        f"forget: n >= k: {exact_late} of {late} exactly the product law p^k on cells, largest distance {worst}"
    )
    print(
        f"forget: n < k, card F >= 2: {positive_early} of {early} at positive distance; uniform: {bound_early} of {early_uniform} at distance >= 1 - card F^(n-k)"
    )
    assert (
        exact_late == late and positive_early == early and bound_early == early_uniform
    )


def de_bruijn(m, k):
    a = [0] * (m * k)
    seq = []

    def db(t, p):
        if t > k:
            if k % p == 0:
                seq.extend(a[1 : p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, m):
                a[t] = j
                db(t + 1, t)

    db(1, 1)
    return seq


def euler_count(m, k):
    states = m ** (k - 1)
    used = {0}
    total = 0

    def walk(s, left):
        nonlocal total
        if left == 0:
            total += s == 0
            return
        for a in range(m):
            w = s * m + a
            if w not in used:
                used.add(w)
                walk(w % states, left - 1)
                used.discard(w)

    walk(0, m**k - 1)
    return total


def paint():
    rows = []
    for base, dim, F, kmax in designs():
        m = len(F)
        if m < 2:
            continue
        x0 = tuple(Fraction(1, 3 + i) for i in range(dim))
        for k in range(1, kmax + 1):
            seq = [F[i] for i in de_bruijn(m, k)]
            lap = len(seq)
            sched = seq + seq + seq[: k - 1]
            x, seen, ok = x0, [], True
            for n, a in enumerate(sched, 1):
                x = step(x, a, base)
                if n >= k:
                    seen.append(cell(x, k, base))
            first = seen[:lap]
            second = seen[lap : 2 * lap]
            cells = {word_cell(w, base) for w in product(F, repeat=k)}
            ok = (
                len(set(first)) == lap == len(cells)
                and set(first) == cells
                and set(second) == cells
                and len(second) == lap
            )
            rows.append((base, dim, m, k, lap + k - 1, ok))
    good = sum(r[-1] for r in rows)
    print(
        f"paint: {good} of {len(rows)} (design, k) de Bruijn schedules visit every level-k cell exactly once at steps k..m^k+k-1 and again once per lap of m^k on the second lap"
    )
    assert good == len(rows)
    for m, k in [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (4, 1), (4, 2)]:
        best = math.factorial(m) ** (m ** (k - 1)) // m**k
        count = euler_count(m, k)
        print(
            f"paint: m={m} k={k}: perfect cyclic paintings {count}, (m!)^(m^(k-1))/m^k = {best}, schedule length m^k+k-1 = {m**k + k - 1}"
        )
        assert count == best


def gauss(M, rhs):
    n = len(M)
    A = [row[:] + [r] for row, r in zip(M, rhs)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        inv = 1 / A[c][c]
        A[c] = [v * inv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [vr - f * vc for vr, vc in zip(A[r], A[c])]
    return [A[r][n] for r in range(n)]


def series(m, k, allowed):
    states = m ** (k - 1)
    M = [[Fraction(int(r == c)) for c in range(states)] for r in range(states)]
    for w in allowed:
        s, t = w // m, w % states
        M[s][t] -= Fraction(1, m)
    y = gauss(M, [Fraction(1)] * states)
    return (k - 1) + sum(y) / Fraction(m) ** (k - 1)


def cover_ie(m, k):
    words = m**k
    full = (1 << words) - 1
    total = Fraction(0)
    for mask in range(1, full + 1):
        allowed = [w for w in range(words) if not (mask >> w) & 1]
        sign = 1 if mask.bit_count() % 2 else -1
        total += sign * series(m, k, allowed)
    return total


def cover_dp(m, k):
    words, states = m**k, m ** (k - 1)
    full = (1 << words) - 1
    E = {full: [Fraction(0)] * states}
    for S in sorted(range(full), key=lambda s: -s.bit_count()):
        M = [[Fraction(int(r == c)) for c in range(states)] for r in range(states)]
        rhs = [Fraction(1)] * states
        for s in range(states):
            for a in range(m):
                w = s * m + a
                t = w % states
                if (S >> w) & 1:
                    M[s][t] -= Fraction(1, m)
                else:
                    rhs[s] += E[S | (1 << w)][t] / m
        E[S] = gauss(M, rhs)
    return (k - 1) + sum(E[0]) / states


def cover_mc(m, k, runs, rng):
    words = m**k
    states = m ** (k - 1)
    total = 0
    for _ in range(runs):
        seen, s, t, got = bytearray(words), 0, 0, 0
        for _ in range(k - 1):
            s = (s * m + rng.randrange(m)) % states
            t += 1
        while got < words:
            w = s * m + rng.randrange(m)
            t += 1
            if not seen[w]:
                seen[w] = 1
                got += 1
            s = w % states
        total += t
    return total / runs


def cover():
    rng = random.Random(SEED)
    for m, k in [(2, 1), (3, 1), (4, 1), (2, 2), (3, 2), (2, 3), (4, 2), (2, 4)]:
        n = m**k
        coupon = (k - 1) + n * sum(Fraction(1, j) for j in range(1, n + 1))
        ie = cover_ie(m, k)
        dp = cover_dp(m, k)
        mc = cover_mc(m, k, 20000, rng)
        agree = (
            "the subset recursion agrees"
            if dp == ie
            else "THE SUBSET RECURSION DISAGREES"
        )
        print(
            f"cover: m={m} k={k}: E[T] = {float(ie):.6f} exactly {ie if ie.denominator < 10**8 else 'a rational of ' + str(len(str(ie.denominator))) + '-digit denominator'}, coupon (k-1) + m^k H(m^k) = {float(coupon):.6f}, excess {float(ie - coupon):+.6f}, {agree}, sampled {mc:.3f} over 20000 runs"
        )
        assert dp == ie
        if k == 1:
            assert ie == coupon


CORNERS = 4
CYCLE = [0, 1, 3, 2]


def vec(c):
    return (c & 1, c >> 1)


def corner(v):
    return v[0] + 2 * v[1]


def symmetries():
    out = []
    for swap in (False, True):
        for fx in (0, 1):
            for fy in (0, 1):

                def g(c, swap=swap, fx=fx, fy=fy):
                    x, y = vec(c)
                    if swap:
                        x, y = y, x
                    return corner((x ^ fx, y ^ fy))

                out.append([g(c) for c in range(CORNERS)])
    return out


def act(perm, code):
    out = 0
    for a in range(CORNERS):
        for b in range(CORNERS):
            if (code >> (4 * a + b)) & 1:
                out |= 1 << (4 * perm[a] + perm[b])
    return out


def pairs(code):
    return [
        (a, b)
        for a in range(CORNERS)
        for b in range(CORNERS)
        if (code >> (4 * a + b)) & 1
    ]


def code_of(rule):
    return sum(1 << (4 * a + b) for a, b in rule)


def letters(code):
    return sorted({c for p in pairs(code) for c in p})


def reach(code, a):
    seen, todo = {a}, [a]
    while todo:
        u = todo.pop()
        for x, y in pairs(code):
            if x == u and y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def irreducible(code):
    L = letters(code)
    return bool(L) and all(reach(code, a) >= set(L) for a in L)


def is_product(code):
    L = letters(code)
    return bool(L) and set(pairs(code)) == {(a, b) for a in L for b in L}


def adjacency(code):
    A = np.zeros((CORNERS, CORNERS), dtype=np.int64)
    for a, b in pairs(code):
        A[a, b] = 1
    return A


def pdiv(p, q):
    p = [Fraction(c) for c in p]
    out = []
    while len(p) >= len(q):
        f = p[0] / q[0]
        out.append(f)
        p = [a - f * b for a, b in zip(p, list(q) + [0] * (len(p) - len(q)))][1:]
    return out, p


def pgcd(p, q):
    p, q = [Fraction(c) for c in p], [Fraction(c) for c in q]
    while any(q):
        while q and q[0] == 0:
            q = q[1:]
        _, r = pdiv(p, q)
        while r and r[0] == 0:
            r = r[1:]
        p, q = q, r
    return [c / p[0] for c in p]


def rho(code):
    p = list(charpoly(code))
    deg = len(p) - 1
    dp = [c * (deg - i) for i, c in enumerate(p[:-1])]
    g = pgcd(p, dp)
    sq, _ = pdiv(p, g)
    roots = np.roots([float(c) for c in sq]) if len(sq) > 1 else []
    return float(max([abs(r) for r in roots], default=0.0))


def charpoly(code):
    A = adjacency(code)
    n = CORNERS
    M = np.zeros((n, n), dtype=object)
    coeffs = [1]
    I = np.identity(n, dtype=object)
    Ao = A.astype(object)
    for k in range(1, n + 1):
        M = Ao.dot(M) + coeffs[-1] * I if k > 1 else I.copy()
        AM = Ao.dot(M)
        c = -sum(AM[i, i] for i in range(n)) // k
        coeffs.append(c)
    return tuple(coeffs)


def mixing_index(code):
    L = letters(code)
    succ = {a: [b for x, b in pairs(code) if x == a] for a in L}
    P = [[Fraction(int(b in succ[a]), len(succ[a])) for b in L] for a in L]
    Q = [row[:] for row in P]
    for m in range(1, len(L) + 2):
        if all(row == Q[0] for row in Q):
            return m
        Q = [
            [sum(Q[i][t] * P[t][j] for t in range(len(L))) for j in range(len(L))]
            for i in range(len(L))
        ]
    return None


def live(rel, allowed):
    S = set(allowed)
    while True:
        T = {c for c in S if any((c, d) in rel and d in S for d in S)}
        if T == S:
            return S
        S = T


def closed_cells(rel, L, level):
    words = [(c,) for c in L]
    for _ in range(level - 1):
        words = [w + (d,) for w in words for d in L if (w[-1], d) in rel]
    feas = {}
    for last in L:
        for e in product((-1, 0, 1), repeat=2):
            allowed = [
                c
                for c in L
                if all(
                    (e[i] == 0)
                    or (e[i] == -1 and vec(c)[i] == 0)
                    or (e[i] == 1 and vec(c)[i] == 1)
                    for i in range(2)
                )
            ]
            lv = live(rel, allowed)
            feas[last, e] = any((last, t) in rel for t in lv)
    out = set()
    side = 2**level
    for w in words:
        base = word_cell(tuple(vec(c) for c in w), 2)
        for e in product((-1, 0, 1), repeat=2):
            if feas[w[-1], e]:
                c = (base[0] + e[0], base[1] + e[1])
                if 0 <= c[0] < side and 0 <= c[1] < side:
                    out.add(c)
    return frozenset(out)


def point_cells(rel, L, level):
    words = [(c,) for c in L]
    for _ in range(level - 1):
        words = [w + (d,) for w in words for d in L if (w[-1], d) in rel]
    loops = [
        w
        for n in range(1, 5)
        for w in product(L, repeat=n)
        if all((w[i], w[(i + 1) % n]) in rel for i in range(n))
    ]
    side = 2**level
    out = set()
    for v in words:
        for pre in [()] + [(c,) for c in L]:
            for w in loops:
                path = v + pre + w
                if any((path[i], path[i + 1]) not in rel for i in range(len(path) - 1)):
                    continue
                if (w[-1], w[0]) not in rel:
                    continue
                pt = []
                for i in range(2):
                    head = sum(
                        Fraction(vec(c)[i], 2 ** (j + 1)) for j, c in enumerate(v + pre)
                    )
                    per = sum(
                        Fraction(vec(c)[i], 2 ** (j + 1)) for j, c in enumerate(w)
                    ) / (1 - Fraction(1, 2 ** len(w)))
                    pt.append(head + per / 2 ** len(v + pre))
                opts = []
                for i in range(2):
                    s = pt[i] * side
                    f = math.floor(s)
                    opts.append([f - 1, f] if s == f else [f])
                for cx in opts[0]:
                    for cy in opts[1]:
                        if 0 <= cx < side and 0 <= cy < side:
                            out.add((cx, cy))
    return frozenset(out)


def game_set(code):
    L = letters(code)
    rel = {(b, a) for a, b in pairs(code)}
    return rel, L


def named():
    def rule(f):
        return code_of(
            [(a, b) for a in range(CORNERS) for b in range(CORNERS) if f(a, b)]
        )

    pos = {c: i for i, c in enumerate(CYCLE)}

    def ahead(a, s):
        return CYCLE[(pos[a] + s) % 4]

    return [
        ("free game", rule(lambda a, b: True)),
        ("no repeat", rule(lambda a, b: b != a)),
        ("not one place anticlockwise", rule(lambda a, b: b != ahead(a, 1))),
        ("not the opposite corner", rule(lambda a, b: b != ahead(a, 2))),
        ("not a neighbour", rule(lambda a, b: b in (a, ahead(a, 2)))),
        ("only a neighbour", rule(lambda a, b: b in (ahead(a, 1), ahead(a, 3)))),
        ("only anticlockwise", rule(lambda a, b: b == ahead(a, 1))),
        ("the shift rule, y_next = x_prev", rule(lambda a, b: vec(b)[1] == vec(a)[0])),
        ("gasket, corners 0 1 2", rule(lambda a, b: a != 3 and b != 3)),
    ]


def census():
    syms = symmetries()
    assert len({tuple(p) for p in syms}) == 8
    canon = {}
    for code in range(1 << 16):
        canon[code] = min(act(p, code) for p in syms)
    classes = sorted(set(canon.values()))
    burnside = Fraction(
        sum(sum(1 for c in range(1 << 16) if act(p, c) == c) for p in syms), len(syms)
    )
    print(
        f"census: {len(classes)} classes of the 65536 one-step rules on the square under its 8 symmetries, Burnside {burnside}"
    )
    assert burnside == len(classes)
    irr = [c for c in classes if irreducible(c)]
    prods = [c for c in irr if is_product(c)]
    rhos_all = {round(rho(c), 9) for c in classes}
    rhos_irr = {}
    for c in irr:
        rhos_irr.setdefault(round(rho(c), 9), []).append(c)
    polys = {charpoly(c) for c in classes}
    print(
        f"census: {len(irr)} irreducible classes (strongly connected on the corners they use), {len(prods)} of them products L x L, codes {prods}"
    )
    print(
        f"census: distinct rho {len(rhos_all)} over all classes, {len(rhos_irr)} over the irreducible ones; distinct characteristic polynomials {len(polys)}"
    )
    by_size = {}
    for c in irr:
        by_size[len(letters(c))] = by_size.get(len(letters(c)), 0) + 1
    print(
        f"census: irreducible classes by corners used {dict(sorted(by_size.items()))}"
    )
    integer = sorted({r for r in rhos_irr if abs(r - round(r)) < 1e-9})
    print(
        f"census: integer rho over the irreducible classes {integer}; non-products among them: {sum(1 for r in integer for c in rhos_irr[r] if not is_product(c))}"
    )
    exact = {}
    for c in irr:
        m = mixing_index(c)
        exact.setdefault(m, []).append(c)
    print(
        f"census: least m with P^m rank one, P the uniform move matrix: {', '.join(f'{k}: {len(v)}' for k, v in sorted(exact.items(), key=lambda kv: (kv[0] is None, kv[0] or 0)))} classes (None: never)"
    )
    gsets = {}
    for r in range(1, CORNERS + 1):
        for G in combinations(range(CORNERS), r):
            rel = {(a, b) for a in G for b in G}
            gsets[G] = [closed_cells(rel, list(G), lv) for lv in range(1, 7)]
    worst, unseparated, levels = 0, [], {}
    for c in irr:
        if is_product(c):
            G = tuple(letters(c))
            rel, L = game_set(c)
            assert all(
                closed_cells(rel, L, lv) == gsets[G][lv - 1] for lv in range(1, 5)
            )
            continue
        rel, L = game_set(c)
        open_ = set(gsets)
        lv = 0
        while open_ and lv < 6:
            lv += 1
            mine = closed_cells(rel, L, lv)
            open_ = {G for G in open_ if gsets[G][lv - 1] == mine}
        if open_:
            unseparated.append(c)
        else:
            levels[lv] = levels.get(lv, 0) + 1
            worst = max(worst, lv)
    print(
        f"census: irreducible non-products whose set differs from every design set at some level <= 6: {sum(levels.values())} of {len(irr) - len(prods)}, by least separating level {dict(sorted(levels.items()))}, unseparated {unseparated}"
    )
    assert not unseparated
    rng = random.Random(SEED)
    sample = rng.sample(irr, 150) + prods
    checks = mism = 0
    for c in sample:
        rel, L = game_set(c)
        for lv in (1, 2, 3):
            checks += 1
            mism += closed_cells(rel, L, lv) != point_cells(rel, L, lv)
    for G in gsets:
        rel = {(a, b) for a in G for b in G}
        for lv in (1, 2, 3):
            checks += 1
            mism += gsets[G][lv - 1] != point_cells(rel, list(G), lv)
    print(
        f"census: closed cells by the extremal-tail rule against cells of exact eventually periodic points, {checks} (set, level) pairs at levels 1..3, {mism} mismatches"
    )
    assert mism == 0
    restricted_forgetting(irr)
    print(
        "census: named games, rule as (previous, next) corner pairs, corner c = x + 2y, cycle 0 1 3 2"
    )
    for name, code in named():
        r = rho(code)
        print(
            f"  {name}: code {code}, class {canon[code]}, card W {len(pairs(code))}, rho {r:.9f}, dim {math.log2(r) if r > 0.5 else float('nan'):.6f}, kappa {math.log2(len(pairs(code))) / 2 - math.log2(r):.6f}, product {is_product(code)}, irreducible {irreducible(code)}, exact after m = {mixing_index(code) if irreducible(code) else '-'}"
        )


def restricted_law(code, first, n, k):
    succ = {a: [b for x, b in pairs(code) if x == a] for a in letters(code)}
    paths = {(first,): Fraction(1)}
    for _ in range(n - 1):
        nxt = {}
        for path, mass in paths.items():
            for b in succ[path[-1]]:
                nxt[path + (b,)] = nxt.get(path + (b,), 0) + mass / len(succ[path[-1]])
        paths = nxt
    law = {}
    for path, mass in paths.items():
        c = word_cell(tuple(vec(a) for a in reversed(path[-k:])), 2)
        law[c] = law.get(c, 0) + mass
    return law


def stationary_law(code, k):
    L = letters(code)
    m = len(L)
    succ = {a: [b for x, b in pairs(code) if x == a] for a in L}
    P = [[Fraction(int(b in succ[a]), len(succ[a])) for b in L] for a in L]
    M = [[P[j][i] - (1 if i == j else 0) for j in range(m)] for i in range(m)]
    M[-1] = [Fraction(1)] * m
    pi = gauss(M, [Fraction(0)] * (m - 1) + [Fraction(1)])
    paths = {(a,): pi[i] for i, a in enumerate(L)}
    for _ in range(k - 1):
        nxt = {}
        for path, mass in paths.items():
            for b in succ[path[-1]]:
                nxt[path + (b,)] = mass / len(succ[path[-1]])
        paths = nxt
    law = {}
    for path, mass in paths.items():
        c = word_cell(tuple(vec(a) for a in reversed(path)), 2)
        law[c] = law.get(c, 0) + mass
    return law


def restricted_forgetting(irr):
    finite, never, ok_finite, ok_never = 0, 0, 0, 0
    rng = random.Random(SEED)
    pool = [c for c in irr if mixing_index(c) is not None and not is_product(c)]
    pool += rng.sample([c for c in irr if mixing_index(c) is None], 40)
    for code in pool:
        nu = mixing_index(code)
        for k in (1, 2):
            want = stationary_law(code, k)
            dist = {
                n: max(tv(restricted_law(code, a, n, k), want) for a in letters(code))
                for n in range(k, k + 6)
            }
            if nu is None:
                never += 1
                ok_never += all(d > 0 for d in dist.values())
            else:
                finite += 1
                ok_finite += all((d == 0) == (n >= k + nu) for n, d in dist.items())
    print(
        f"census: restricted games, first draw any corner of the rule: level-k cell law exactly stationary from n = k + m on and not before, m the rank-one index, in {ok_finite} of {finite} (class, k) cases with finite m; positive distance at every n = k..k+5 in {ok_never} of {never} sampled cases with no finite m"
    )
    assert ok_finite == finite and ok_never == never


def width_three():
    pos = {c: i for i, c in enumerate(CYCLE)}
    nb = lambda a, b: (pos[b] - pos[a]) % 4 in (1, 3)
    W = [
        (a, b, c)
        for a in range(4)
        for b in range(4)
        for c in range(4)
        if not (a == b and nb(b, c))
    ]
    states = [(a, b) for a in range(4) for b in range(4)]
    idx = {s: i for i, s in enumerate(states)}
    A = np.zeros((16, 16))
    for a, b, c in W:
        A[idx[a, b], idx[b, c]] += 1
    r = float(max(abs(np.linalg.eigvals(A))))
    kappa = math.log2(len(W)) / 3 - math.log2(r)
    lumped = float(max(abs(np.linalg.eigvals(np.array([[1.0, 1.0], [1.0, 3.0]])))))
    print(
        f"census: no neighbour after a repeat, width 3: card W {len(W)}, rho {r:.9f}, lumped [[1, 1], [1, 3]] rho {lumped:.9f} = 2 + sqrt 2, dim {math.log2(r):.6f}, kappa {kappa:.6f}"
    )
    assert abs(r - 2 - math.sqrt(2)) < 1e-9 and abs(lumped - r) < 1e-9


def square():
    census()
    width_three()


GASKET = [(0, 0), (1, 0), (0, 1)]
DIAGONAL = [(0, 0), (1, 1)]


def magic_law(x0, word, n, k):
    law = {}
    for draws in product(*[word[j % len(word)] for j in range(n)]):
        mass = math.prod(
            (Fraction(1, len(word[j % len(word)])) for j in range(n)), start=Fraction(1)
        )
        c = cell(play(x0, draws, 2), k, 2)
        law[c] = law.get(c, 0) + mass
    return law


def magic_target(letters_, k):
    cells = [word_cell(w, 2) for w in product(*letters_)]
    return {c: Fraction(1, len(cells)) for c in cells}


def phase_painting(word, k):
    p = len(word)
    edges = {}
    for r in range(p):
        for u in product(*[word[(r - k + 1 + j) % p] for j in range(k - 1)]):
            edges[(r, u)] = list(word[r])
    indeg = {}
    for (r, u), outs in edges.items():
        for a in outs:
            v = ((r + 1) % p, (u + (a,))[1:])
            indeg[v] = indeg.get(v, 0) + 1
    if not all(len(outs) == indeg.get(v, 0) for v, outs in edges.items()):
        return None
    start = (0, tuple(word[(j - k + 1) % p][0] for j in range(k - 1)))
    left = {v: list(outs) for v, outs in edges.items()}
    stack, circuit = [(start, None)], []
    while stack:
        v, a = stack[-1]
        if left[v]:
            b = left[v].pop()
            stack.append((((v[0] + 1) % p, (v[1] + (b,))[1:]), (v, b)))
        else:
            stack.pop()
            if a is not None:
                circuit.append(a)
    circuit.reverse()
    assert len(circuit) == sum(len(o) for o in edges.values())
    return start, circuit


def magic():
    x0 = (Fraction(1, 3), Fraction(1, 5))
    words = [(GASKET, DIAGONAL), (GASKET, GASKET, DIAGONAL), (DIAGONAL, GASKET)]
    cases = exact = early = early_pos = 0
    for word in words:
        p = len(word)
        for k in range(1, 4):
            for n in range(k + 3):
                law = magic_law(x0, word, n, k)
                rev = [word[(n - 1 - j) % p] for j in range(k)]
                d = tv(law, magic_target(rev, k))
                cases += 1
                if n >= k:
                    exact += d == 0
                else:
                    early += 1
                    early_pos += d > 0
    print(
        f"magic: non-autonomous game, step n draws from the letter s_n of a periodic word; level-k cell law at n >= k equals uniform on the magic design of the reversed window (s_n, ..., s_(n-k+1)) in {exact} of {cases - early} cases, positive distance in {early_pos} of {early} cases n < k"
    )
    assert exact == cases - early and early_pos == early
    for word in [(GASKET, DIAGONAL), (GASKET, GASKET, DIAGONAL)]:
        p = len(word)
        sizes = [len(w) for w in word]
        for k in range(1, 5):
            cond = all(sizes[r] == sizes[(r - k) % p] for r in range(p))
            found = phase_painting(word, k)
            verdict = "no perfect lap"
            if found:
                (_r0, u0), circuit = found
                window, seen = list(u0), {}
                for (r, u), a in circuit:
                    window.append(a)
                    c = word_cell(tuple(reversed(window[-k:])), 2)
                    seen.setdefault(r, []).append(c)
                good = all(
                    sorted(seen.get(r, []))
                    == sorted(magic_target([word[(r - j) % p] for j in range(k)], k))
                    for r in range(p)
                )
                verdict = f"perfect lap of {len(circuit)} steps, every phase visits every cell of its design exactly once: {good}"
                assert good
            print(
                f"magic: letter sizes {sizes} k={k}: sizes balanced {cond}, Euler circuit {found is not None}; {verdict}"
            )
            assert cond == (found is not None)


MOSAIC = {1: {0: 1, 1: 2, 2: 1}, 2: {0: 2, 3: 1}}


def substitute(grid):
    side = grid.shape[0] * 2
    new = np.zeros((side, side), dtype=int)
    for yy in range(grid.shape[0]):
        for xx in range(grid.shape[1]):
            t = grid[yy, xx]
            if t:
                for c, v in MOSAIC[t].items():
                    x, y = vec(c)
                    new[2 * yy + y, 2 * xx + x] = v
    return new


def mosaic():
    types = sorted(MOSAIC)
    S = np.array(
        [[sum(1 for t in MOSAIC[u].values() if t == v) for v in types] for u in types]
    )
    r = float(max(abs(np.linalg.eigvals(S))))
    grid, counts, grids = np.array([[1]]), [], {}
    for level in range(1, 9):
        grid = substitute(grid)
        grids[level] = grid
        direct = int((grid > 0).sum())
        matrix = int(
            (np.linalg.matrix_power(S, level) @ np.ones(len(types), dtype=int))[0]
        )
        assert direct == matrix
        counts.append(direct)
    level = 6
    cells = {
        (xx, yy)
        for yy in range(2**level)
        for xx in range(2**level)
        if grids[level][yy, xx]
    }
    incoming = {
        v: [(u, c) for u in types for c, t in MOSAIC[u].items() if t == v]
        for v in types
    }
    rng = random.Random(SEED)
    state, hist, visited, stray = 1, [], set(), 0
    for _ in range(200000):
        u, c = rng.choice(incoming[state])
        hist.append(vec(c))
        state = u
        if state == 1 and len(hist) >= level:
            cc = word_cell(tuple(reversed(hist[-level:])), 2)
            if cc in cells:
                visited.add(cc)
            else:
                stray += 1
    print(
        f"mosaic: two types, masks {MOSAIC}; substitution matrix {S.tolist()}, rho {r:.9f}, dim {math.log2(r):.6f}; filled cells by substitution equal e_A S^level 1 at levels 1..8: {counts}"
    )
    print(
        f"mosaic: graph game, next map drawn from the edges into the current type: {len(visited)} of {len(cells)} level-6 cells of the type-A design visited in 200000 steps, {stray} visits outside it"
    )
    assert stray == 0 and visited == cells


def slice_graph(F, sigma):
    w = {d: d[0] + d[1] for d in F}
    lo, hi = min(w.values()), max(w.values())
    V, todo, edges = {sigma}, [sigma], []
    while todo:
        h = todo.pop()
        for d in F:
            h2 = 2 * h - w[d]
            if lo <= h2 <= hi:
                edges.append((h, h2, d))
                if h2 not in V:
                    V.add(h2)
                    todo.append(h2)
    live = set(V)
    while True:
        nxt = {h for h in live if any(a == h and b in live for a, b, _ in edges)}
        if nxt == live:
            break
        live = nxt
    edges = [(a, b, d) for a, b, d in edges if a in live and b in live]
    order = sorted(live)
    idx = {h: i for i, h in enumerate(order)}
    A = np.zeros((len(order), len(order)), dtype=np.int64)
    for a, b, _ in edges:
        A[idx[a], idx[b]] += 1
    return order, idx, A, edges, lo, hi


def slices():
    sigma = Fraction(2, 3)
    for name, F in [
        ("code 7, corners 00 10 01, weights 0 1 1", GASKET),
        ("code 11, corners 00 10 11, weights 0 1 2", [(0, 0), (1, 0), (1, 1)]),
    ]:
        order, idx, A, edges, lo, hi = slice_graph(F, sigma)
        r = float(max(abs(np.linalg.eigvals(A.astype(float)))))
        agree = True
        rows = []
        for level in range(1, 13):
            mat = int(
                (
                    np.linalg.matrix_power(A, level)
                    @ np.ones(len(order), dtype=np.int64)
                )[idx[sigma]]
            )
            target = sigma * 2**level
            brute = sum(
                1
                for word in product(F, repeat=level)
                if lo
                <= target
                - sum((d[0] + d[1]) * 2 ** (level - 1 - j) for j, d in enumerate(word))
                <= hi
            )
            agree &= mat == brute
            rows.append(brute)
        print(
            f"slices: {name}, line x + y = 2/3: carry states {[str(h) for h in order]}, matrix {A.tolist()}, rho {r:.9f}, dim {math.log2(r):.6f}; level-L cells meeting the line by integer enumeration and by the matrix agree at L = 1..12: {agree}, counts {rows}"
        )
        assert agree
        rng = random.Random(SEED)
        into = {}
        for a, b, d in edges:
            into.setdefault(b, []).append((a, d))
        state, x, hist, worst, visited, level = sigma, (0.5, 0.25), [], 0.0, set(), 6
        start_gap = abs(x[0] + x[1] - float(sigma))
        cyl = set()
        for word in product(F, repeat=level):
            h = sigma
            okw = True
            for d in word:
                h = 2 * h - (d[0] + d[1])
                if h not in idx:
                    okw = False
                    break
            if okw:
                cyl.add(word_cell(word, 2))
        stray = 0
        for n in range(1, 100001):
            u, d = rng.choice(into[state])
            x = ((x[0] + d[0]) / 2, (x[1] + d[1]) / 2)
            hist.append(d)
            state = u
            if n > 60:
                worst = max(worst, abs(x[0] + x[1] - float(state)))
            if state == sigma and n >= level:
                c = word_cell(tuple(reversed(hist[-level:])), 2)
                if c in cyl:
                    visited.add(c)
                else:
                    stray += 1
        print(
            f"slices: carry-state game: after 60 steps |x + y - state| <= {worst:.1e} (start gap {start_gap}); level-6 cylinders of the slice visited at state 2/3: {len(visited)} of {len(cyl)}, {stray} outside"
        )
        assert stray == 0 and visited == cyl and worst < 1e-12


def families():
    magic()
    mosaic()
    slices()


VERBS = {
    "register": register,
    "forget": forget,
    "paint": paint,
    "cover": cover,
    "census": square,
    "families": families,
}

if __name__ == "__main__":
    names = sys.argv[1:] or list(VERBS)
    for name in names:
        VERBS[name]()
