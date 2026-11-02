from math import gcd

GASKET = ((0, 0), (1, 0), (0, 1))
FIB = [0, 1, 1]
while len(FIB) < 200:
    FIB.append(FIB[-1] + FIB[-2])


def die(name, got, want):
    raise SystemExit("FAIL %s\n  got  %r\n  want %r" % (name, got, want))


def check(name, got, want):
    if got != want:
        die(name, got, want)
    print("ok   %s" % name)


# GASKET POINTS


def points(n):
    pts = [(0, 0)]
    for _ in range(n):
        nxt = []
        for x, y in pts:
            for dx, dy in GASKET:
                nxt.append((3 * x + dx, 3 * y + dy))
        pts = nxt
    return pts


def in_gasket(x, y, n):
    if x >= 3**n or y >= 3**n:
        return False
    while x or y:
        a, b = x % 3, y % 3
        if a > 1 or b > 1 or (a and b):
            return False
        x //= 3
        y //= 3
    return True


def is_shift_pair(s, t):
    if min(s, t) != 1:
        return False
    u = max(s, t)
    while u % 3 == 0:
        u //= 3
    return u == 1


# CARRY AUTOMATON


def carry(s, t):
    index = {(0, 0): 0}
    order = [(0, 0)]
    edges = []
    i = 0
    while i < len(order):
        a, b = order[i]
        out = []
        for d in range(3):
            u = (s * d + a) % 3
            v = (t * d + b) % 3
            if u < 2 and v < 2:
                nxt = ((s * d + a) // 3, (t * d + b) // 3)
                if nxt not in index:
                    index[nxt] = len(order)
                    order.append(nxt)
                out.append((u, v, index[nxt]))
        edges.append(out)
        i += 1
    return edges


def carry_live(edges):
    back = [[] for _ in edges]
    for i, out in enumerate(edges):
        for u, v, j in out:
            back[j].append(i)
    seen = {0}
    stack = [0]
    while stack:
        x = stack.pop()
        for y in back[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    keep = sorted(seen)
    place = {v: i for i, v in enumerate(keep)}
    return [[(u, v, place[j]) for (u, v, j) in edges[x] if j in place] for x in keep]


def ray_mass(z1, z2, n):
    e = carry_live(carry(z1, z2))
    k = len(e)
    adj = [[j for (u, v, j) in e[i] if not (u and v)] for i in range(k)]
    cur = [0] * k
    cur[0] = 1
    out = [1]
    for _ in range(n):
        nxt = [0] * k
        for i in range(k):
            c = cur[i]
            if c:
                for j in adj[i]:
                    nxt[j] += c
        cur = nxt
        out.append(cur[0])
    return [v - 1 for v in out]


# FREE AUTOMATON AS A CONSTRAINED TENSOR SQUARE


def free_automaton(s, t):
    index = {(0, 0, 0, 0): 0}
    order = [(0, 0, 0, 0)]
    edges = []
    i = 0
    while i < len(order):
        a1, a2, b1, b2 = order[i]
        out = []
        for dx in range(3):
            for dy in range(3):
                u = ((s * dx + a1) % 3, (s * dy + a2) % 3)
                v = ((t * dx + b1) % 3, (t * dy + b2) % 3)
                if u in GASKET and v in GASKET:
                    nxt = (
                        (s * dx + a1) // 3,
                        (s * dy + a2) // 3,
                        (t * dx + b1) // 3,
                        (t * dy + b2) // 3,
                    )
                    if nxt not in index:
                        index[nxt] = len(order)
                        order.append(nxt)
                    out.append(index[nxt])
        edges.append(out)
        i += 1
    return edges


def free_returns(edges, n):
    count = [0] * len(edges)
    count[0] = 1
    out = []
    for _ in range(n + 1):
        out.append(count[0])
        nxt = [0] * len(edges)
        for i, outs in enumerate(edges):
            for e in outs:
                nxt[e] += count[i]
        count = nxt
    return out


def tensor_returns(s, t, n):
    e = carry_live(carry(s, t))
    k = len(e)
    S = [[j for (u, v, j) in e[i]] for i in range(k)]
    U = [[j for (u, v, j) in e[i] if u] for i in range(k)]
    V = [[j for (u, v, j) in e[i] if v] for i in range(k)]
    W = [[j for (u, v, j) in e[i] if u and v] for i in range(k)]
    X = [[0] * k for _ in range(k)]
    X[0][0] = 1
    out = [1]
    for _ in range(n):
        half = []
        for A in (S, U, V, W):
            Y = []
            for i in range(k):
                acc = [0] * k
                for p in A[i]:
                    xp = X[p]
                    for q in range(k):
                        if xp[q]:
                            acc[q] += xp[q]
                Y.append(acc)
            half.append(Y)
        Ys, Yu, Yv, Yw = half
        Z = [[0] * k for _ in range(k)]
        for j in range(k):
            for i in range(k):
                a = 0
                for q in S[j]:
                    a += Ys[i][q]
                for q in U[j]:
                    a -= Yu[i][q]
                for q in V[j]:
                    a -= Yv[i][q]
                for q in W[j]:
                    a += Yw[i][q]
                Z[i][j] = a
        X = Z
        out.append(X[0][0])
    return out


# BOX CONSTRUCTION


def box(s, t, n):
    W = (3**n - 1) // (2 * max(s, t))
    c = 0
    for z1 in range(1, W):
        for z2 in range(1, W - z1 + 1):
            if in_gasket(s * z1, s * z2, n) and in_gasket(t * z1, t * z2, n):
                c += 1
    return c


def fibre(s, t, n):
    binary = {sum(((u >> i) & 1) * 3**i for i in range(n)) for u in range(2**n)}
    return sum(
        1 for u in range(1, 3**n // t + 1) if s * u in binary and t * u in binary
    )


# CENSUS


def ray_census(n):
    by = {}
    for x, y in points(n):
        if x and y:
            g = gcd(x, y)
            by.setdefault((x // g, y // g), []).append(g)
    return by


def residual(n):
    by = ray_census(n)
    R = 0
    weight = {}
    pair = {}
    for r, gs in by.items():
        if len(gs) < 2:
            continue
        rs = r[0] + r[1]
        for a in gs:
            for b in gs:
                if a == b:
                    continue
                d = gcd(a, b)
                if is_shift_pair(a // d, b // d):
                    continue
                R += 1
                weight[d * rs] = weight.get(d * rs, 0) + 1
                k = (a // d, b // d)
                if k[0] > k[1]:
                    k = (k[1], k[0])
                pair[k] = pair.get(k, 0) + 1
    return R, weight, pair, by


def no_adjacent(n):
    out = []
    for m in range(1, 3 ** (n - 1)):
        x, prev, ok = m, 0, True
        while x:
            d = x % 3
            if d > 1 or (d and prev):
                ok = False
                break
            prev, x = d, x // 3
        if ok:
            out.append(m)
    return out


def is_three_power_ratio(a, b):
    if a > b:
        a, b = b, a
    if b % a:
        return False
    q = b // a
    while q % 3 == 0:
        q //= 3
    return q == 1


# LAW 1: THE CONSTRAINED TENSOR SQUARE


def law_tensor():
    bad = 0
    n = 9
    tested = 0
    for s in range(1, 40):
        for t in range(s + 1, 40):
            if gcd(s, t) != 1:
                continue
            tested += 1
            if free_returns(free_automaton(s, t), n) != tensor_returns(s, t, n):
                bad += 1
    check(
        "tensor square reproduces B(s,t) return counts, n <= 9", (tested, bad), (473, 0)
    )
    sizes = []
    for s, t in ((365, 1094), (41, 122), (25, 52), (31, 40)):
        sizes.append((len(carry_live(carry(s, t))), len(free_automaton(s, t))))
    check(
        "carry states against reachable B states",
        sizes,
        [(729, 26931), (81, 835), (38, 393), (35, 354)],
    )


# LAW 2: THE BOX CONSTRUCTION


def law_box():
    n = 9
    bad = 0
    tested = 0
    for s in range(1, 30):
        for t in range(s + 1, 60):
            if gcd(s, t) != 1:
                continue
            tested += 1
            want = free_returns(free_automaton(s, t), n)[n] - 1 - 2 * fibre(s, t, n)
            if want != box(s, t, n):
                bad += 1
    check("box construction against B(s,t), n = 9", (tested, bad), (812, 0))
    got = []
    for s, t in ((365, 1094), (41, 122), (122, 123), (1, 2460), (2431, 2458)):
        for n in (9, 12):
            a = tensor_returns(s, t, n)[n] - 1 - 2 * fibre(s, t, n)
            got.append((a, box(s, t, n)))
    check(
        "box against tensor at large multipliers, n = 9 and 12",
        [x for x in got if x[0] != x[1]],
        [],
    )
    check(
        "witness counts at the four extreme pairs, n = 12",
        [g[0] for g in got],
        [2, 180, 50, 172, 50, 172, 2, 12, 2, 8],
    )


# LAW 3: WEIGHT FOUR AND THE GOLDEN CEILING


def law_weights():
    R = {}
    weight = {}
    pair = {}
    for n in range(4, 14):
        R[n], weight[n], pair[n], _ = residual(n)
    check(
        "R(n) for n = 4..13",
        [R[n] for n in range(4, 14)],
        [20, 88, 432, 1624, 5512, 15896, 46064, 124928, 335704, 863848],
    )
    scaled = 0
    bad = 0
    for n in range(5, 14):
        for w, c in weight[n].items():
            if w % 3 == 0:
                scaled += 1
                if weight[n - 1].get(w // 3, 0) != c:
                    bad += 1
    check("R_{3w}(n) = R_w(n-1)", (scaled, bad), (1869, 0))
    check(
        "no witness of weight below four",
        sorted({min(weight[n]) for n in range(4, 14)}),
        [4],
    )
    fours = []
    ok = True
    for n in range(4, 13):
        F = no_adjacent(n)
        if len(F) != FIB[n + 1] - 1:
            ok = False
        c = sum(
            1
            for a in F
            for b in F
            if a != b and gcd(a, b) == 1 and not is_three_power_ratio(a, b)
        )
        fours.append(2 * c)
        if 2 * c != weight[n].get(4, 0):
            ok = False
    check(
        "R_4(n) counted by the no-adjacent-ones set, n = 4..12",
        (ok, fours),
        (True, [12, 36, 108, 336, 988, 2596, 6672, 17480, 45720]),
    )
    check(
        "|F_n| = Fib(n+1) - 1 for n = 4..12",
        [len(no_adjacent(n)) for n in range(4, 13)],
        [FIB[n + 1] - 1 for n in range(4, 13)],
    )
    tops = []
    for n in range(6, 14):
        thr = (3**n - 1) // 10
        big = [v for k, v in pair[n].items() if k[1] > thr]
        tops.append((len(big), sorted(set(big))))
    check(
        "every pair above (3^n-1)/10 contributes exactly 4, n = 6..13",
        tops,
        [
            (18, [4]),
            (57, [4]),
            (163, [4]),
            (402, [4]),
            (1019, [4]),
            (2702, [4]),
            (7060, [4]),
            (18607, [4]),
        ],
    )
    check(
        "heaviest ray mass is Fib(n+1) - 1, n = 1..40",
        ray_mass(1, 3, 40),
        [FIB[n + 1] - 1 for n in range(41)],
    )
    N = 40
    ref = [FIB[n + 1] - 1 for n in range(N + 1)]
    tested = 0
    breaches = 0
    ties = []
    for z1 in range(1, 121):
        for z2 in range(z1, 241):
            if gcd(z1, z2) != 1:
                continue
            tested += 1
            m = ray_mass(z1, z2, N)
            for n in range(1, N + 1):
                if m[n] > ref[n]:
                    breaches += 1
            if m[N] == ref[N]:
                ties.append((z1, z2))
    check(
        "golden ceiling M_n(z) <= Fib(n+1) - 1 over the box, n <= 40",
        (tested, breaches, ties),
        (13158, 0, [(1, 3)]),
    )
    N = 45
    ref = [FIB[m + 1] - 1 for m in range(N + 1)]
    bin3 = [sum(((u >> i) & 1) * 3**i for i in range(8)) for u in range(1, 256)]
    nad = no_adjacent(8)
    families = [
        [(a, b) for a in bin3 for b in bin3 if a < b],
        [(a, b) for a in nad for b in nad if a < b],
        [(1, t) for t in range(2, 3000)],
        [(a, a + 1) for a in range(1, 1500)],
        [(a, 3 * a - 1) for a in range(1, 1200)],
        [(a, 3 * a + 1) for a in range(1, 1200)],
    ]
    sizes = []
    breaches = 0
    for family in families:
        seen = 0
        for a, b in family:
            if gcd(a, b) != 1:
                continue
            seen += 1
            m = ray_mass(a, b, N)
            for j in range(1, N + 1):
                if m[j] > ref[j]:
                    breaches += 1
        sizes.append(seen)
    check(
        "golden ceiling on six adversarial families, n <= 45",
        (sizes, sum(sizes), breaches),
        ([16940, 253, 2998, 1499, 1199, 1199], 24088, 0),
    )


# LAW 4: THE RESIDUAL TO LEVEL 17


def law_deep():
    import numpy as np

    E, R = [], []
    for n in range(13, 18):
        X = np.zeros(1, dtype=np.int32)
        Y = np.zeros(1, dtype=np.int32)
        for _ in range(n):
            X = np.concatenate((3 * X, 3 * X + 1, 3 * X))
            Y = np.concatenate((3 * Y, 3 * Y, 3 * Y + 1))
        keep = (X > 0) & (Y > 0)
        X, Y = X[keep], Y[keep]
        del keep
        g = np.gcd(X, Y)
        X //= g
        Y //= g
        del g
        key = X.astype(np.int64)
        key *= 3**n
        key += Y
        del X, Y
        key.sort()
        idx = np.flatnonzero(np.concatenate(([True], key[1:] != key[:-1])))
        m = np.diff(np.concatenate((idx, [len(key)]))).astype(np.int64)
        e = int((m * m).sum())
        del key, idx, m
        E.append(e)
        R.append(e - (3**n - 2 ** (n + 1) + 1) - (3**n - 4 * 2**n + 2 * n + 3))
    check("E(n) for n = 13..17", E, [4003372, 11679626, 34050692, 99800950, 292848756])
    check("R(n) for n = 13..17", R, [863848, 2211960, 5549452, 14100688, 35354824])
    peak3 = max(range(len(R)), key=lambda i: R[i] / 3 ** (13 + i))
    fall3 = all(R[i + 1] / 3 ** (14 + i) < R[i] / 3 ** (13 + i) for i in range(4))
    p = (3 + 5**0.5) / 2
    fallp = all(R[i + 1] / p ** (14 + i) < R[i] / p ** (13 + i) for i in range(4))
    check(
        "R/3^n and R/phi^2n both fall through n = 17",
        (peak3, fall3, fallp),
        (0, True, True),
    )
    check(
        "R(17)/3^17 and R(17)/phi^34 truncated down",
        (int(R[4] / 3**17 * 10**7), int(R[4] / p**17 * 10**7)),
        (2737709, 27724831),
    )


# LAW 5: THE CEILING AS A THEOREM


def direction_automaton(a, b):
    index = {0: 0}
    order = [0]
    edges = []
    i = 0
    while i < len(order):
        c = order[i]
        out = []
        for eps in (0, b, -a):
            if (c + eps) % 3 == 0:
                nxt = (c + eps) // 3
                if nxt not in index:
                    index[nxt] = len(order)
                    order.append(nxt)
                out.append((eps, index[nxt]))
        edges.append(out)
        i += 1
    return order, edges


def direction_live(a, b):
    order, edges = direction_automaton(a, b)
    back = [[] for _ in edges]
    for i, out in enumerate(edges):
        for eps, j in out:
            back[j].append(i)
    seen = {0}
    stack = [0]
    while stack:
        x = stack.pop()
        for y in back[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    keep = sorted(seen)
    place = {v: i for i, v in enumerate(keep)}
    return (
        [order[x] for x in keep],
        [[(eps, place[j]) for (eps, j) in edges[x] if j in place] for x in keep],
    )


def direction_mass(a, b, n):
    _order, edges = direction_live(a, b)
    cur = [0] * len(edges)
    cur[0] = 1
    out = [1]
    for _ in range(n):
        nxt = [0] * len(edges)
        for i in range(len(edges)):
            if cur[i]:
                for eps, j in edges[i]:
                    nxt[j] += cur[i]
        cur = nxt
        out.append(cur[0])
    return [v - 1 for v in out]


def state_profile(a, b, n):
    order, edges = direction_live(a, b)
    tab = [[1 if order[i] == 0 else 0 for i in range(len(edges))]]
    for _ in range(n):
        prev = tab[-1]
        tab.append([sum(prev[j] for (eps, j) in edges[i]) for i in range(len(edges))])
    return [max(r) for r in tab]


def first_returns(a, b, n):
    _order, edges = direction_live(a, b)
    cur = [0] * len(edges)
    cur[0] = 1
    f = [0] * (n + 1)
    for m in range(1, n + 1):
        nxt = [0] * len(edges)
        for i in range(len(edges)):
            if cur[i]:
                for eps, j in edges[i]:
                    nxt[j] += cur[i]
        f[m] = nxt[0]
        nxt[0] = 0
        cur = nxt
    return f


def valuation3(x):
    k = 0
    while x % 3 == 0:
        x //= 3
        k += 1
    return k


def double_branch(order, edges):
    degs = [len(o) for o in edges]
    return any(
        degs[i] == 2 and all(degs[j] == 2 for (eps, j) in edges[i])
        for i in range(len(edges))
    )


def is_shift_ray(a, b):
    lo, hi = min(a, b), max(a, b)
    if lo != 1:
        return False
    while hi % 3 == 0:
        hi //= 3
    return hi == 1


def certificate_holds(a, b, den, alpha, beta):
    order, edges = direction_live(a, b)
    if len(order) != len(alpha) or len(order) != len(beta):
        return False
    if alpha[0] != den or beta[0] != 0:
        return False
    for i in range(len(edges)):
        if alpha[i] < 0:
            return False
        if sum(alpha[j] for (eps, j) in edges[i]) > alpha[i] + beta[i]:
            return False
        if sum(beta[j] for (eps, j) in edges[i]) > alpha[i]:
            return False
    return True


def binary_multiples(w, n):
    cur = [0] * w
    cur[0] = 1
    for i in range(n):
        p = pow(3, i, w)
        nxt = list(cur)
        for r in range(w):
            if cur[r]:
                nxt[(r + p) % w] += cur[r]
        cur = nxt
    return cur[0] - 1


def shift_product(j, n):
    if n <= j:
        return 1
    total = 1
    for r in range(j):
        m = len(range(r, n - j, j))
        total *= FIB[m + 2]
    return total


CERTIFICATES = {
    (1, 90): (18, [18, 0, 5, 6, 9, 2, 10, 4, 8], [0, 18, 4, -4, 5, 6, 8, 1, 2]),
    (4, 117): (
        40,
        [40, 0, 5, 14, 6, 13, 11, 27, 17, 1],
        [0, 40, 1, -1, 5, 14, 6, 13, 11, 4],
    ),
    (9, 73): (
        381,
        [381, 0, 46, 49, 78, 17, 101, 23, 66, 163, 39, 83, 232, 32, 62, 149],
        [0, 381, 32, -32, 46, 49, 62, 16, 17, 101, 23, 66, 149, 14, 39, 83],
    ),
    (9, 82): (18, [18, 0, 5, 6, 9, 2, 10, 4, 8], [0, 18, 4, -4, 5, 6, 8, 1, 2]),
    (9, 235): (
        2013,
        [
            2013,
            0,
            16,
            66,
            27,
            55,
            43,
            121,
            70,
            176,
            113,
            297,
            183,
            473,
            296,
            770,
            479,
            1243,
            775,
            11,
        ],
        [
            0,
            2013,
            11,
            -11,
            16,
            66,
            27,
            55,
            43,
            121,
            70,
            176,
            113,
            297,
            183,
            473,
            296,
            770,
            479,
            5,
        ],
    ),
    (10, 81): (18, [18, 0, 5, 6, 9, 2, 10, 4, 8], [0, 18, 4, -4, 5, 6, 8, 1, 2]),
    (13, 108): (
        40,
        [40, 0, 5, 14, 6, 13, 11, 27, 17, 1],
        [0, 40, 1, -1, 5, 14, 6, 13, 11, 4],
    ),
    (27, 217): (
        2013,
        [
            2013,
            0,
            16,
            66,
            27,
            55,
            43,
            121,
            70,
            176,
            113,
            297,
            183,
            473,
            296,
            770,
            479,
            1243,
            775,
            11,
        ],
        [
            0,
            2013,
            11,
            -11,
            16,
            66,
            27,
            55,
            43,
            121,
            70,
            176,
            113,
            297,
            183,
            473,
            296,
            770,
            479,
            5,
        ],
    ),
    (27, 226): (
        34,
        [34, 0, 0, 1, 1, 0, 1, 1, 2, 1, 3, 5, 8, 13, 20, 1],
        [0, 34, 1, -1, 0, 1, 1, 0, 1, 1, 2, 3, 5, 8, 14, -1],
    ),
}


def law_ceiling():
    tested = mismatch = 0
    for a in range(1, 40):
        for b in range(1, 40):
            if gcd(a, b) != 1:
                continue
            tested += 1
            if direction_mass(a, b, 20) != ray_mass(a, b, 20):
                mismatch += 1
    check(
        "direction carry automaton against the gasket-digit automaton",
        (tested, mismatch),
        (947, 0),
    )
    box = occupied = bounded = single = k1 = nodouble = 0
    twoplus = anomaly = states = renewal = 0
    exceptions = []
    twos = []
    for z1 in range(1, 121):
        for z2 in range(z1, 241):
            if gcd(z1, z2) != 1:
                continue
            box += 1
            d = (z1 % 3 == 0) + (z2 % 3 == 0) + ((z1 + z2) % 3 == 0)
            if d > 1:
                twoplus += 1
            order, edges = direction_live(z1, z2)
            if len(order) == 1:
                continue
            if d == 0:
                anomaly += 1
            occupied += 1
            states = max(states, len(order))
            if all(-z1 <= 2 * c <= z2 for c in order):
                bounded += 1
            degs = [len(o) for o in edges]
            classes = {order[i] % 3 for i in range(len(order)) if degs[i] == 2}
            if max(degs) <= 2 and len(classes) <= 1:
                single += 1
            q = z1 if z1 % 3 == 0 else (z2 if z2 % 3 == 0 else z1 + z2)
            if valuation3(q) == 1:
                k1 += 1
            if not double_branch(order, edges):
                nodouble += 1
            else:
                exceptions.append((z1, z2))
            f = first_returns(z1, z2, 46)
            if f[1] == 1 and all(
                sum(f[j] * FIB[m + 1 - j] for j in range(2, m + 1)) <= FIB[m - 1]
                for m in range(2, 47)
            ):
                renewal += 1
            if f[2]:
                twos.append((z1, z2))
    check(
        "carry states of a direction lie in [-a/2, b/2]",
        (box, occupied, states, bounded),
        (13158, 218, 37, 218),
    )
    check("out-degree at most two, branch states in one class mod 3", single, 218)
    check(
        "at most one of z1, z2, w is divisible by three, and no occupied "
        "direction has none of them divisible",
        (twoplus, anomaly),
        (0, 0),
    )
    check(
        "directions settled by the branch argument over the box",
        (k1, nodouble, exceptions),
        (
            107,
            206,
            [
                (1, 9),
                (1, 27),
                (1, 81),
                (1, 90),
                (4, 117),
                (9, 73),
                (9, 82),
                (9, 235),
                (10, 81),
                (13, 108),
                (27, 217),
                (27, 226),
            ],
        ),
    )
    check(
        "renewal criterion on every occupied direction of the box, n <= 46",
        (renewal, twos),
        (218, [(1, 3)]),
    )
    good = sorted(
        k
        for k, v in CERTIFICATES.items()
        if certificate_holds(k[0], k[1], v[0], v[1], v[2])
    )
    check(
        "Fibonacci certificates for the nine non-shift exceptions",
        good,
        [
            (1, 90),
            (4, 117),
            (9, 73),
            (9, 82),
            (9, 235),
            (10, 81),
            (13, 108),
            (27, 217),
            (27, 226),
        ],
    )
    identity = all(
        FIB[p + 2] * FIB[q + 2] == FIB[p + q + 3] - FIB[p + 1] * FIB[q + 1]
        for p in range(60)
        for q in range(60)
    )
    cases = over = strict = 0
    for j in range(1, 14):
        for m in range(1, 46):
            cases += 1
            if shift_product(j, m) > FIB[m + 1]:
                over += 1
            if j >= 2 and m >= 2 and shift_product(j, m) >= FIB[m + 1]:
                strict += 1
    check(
        "Fibonacci product identity and the shift-ray ceiling",
        (identity, cases, over, strict),
        (True, 585, 0, 0),
    )
    tested = breaches = 0
    for z1 in range(1, 31):
        for z2 in range(z1, 61):
            if gcd(z1, z2) != 1:
                continue
            tested += 1
            if direction_mass(z1, z2, 12)[12] > binary_multiples(z1 + z2, 12):
                breaches += 1
    check(
        "ray mass is at most the count of binary multiples of the weight",
        (tested, breaches),
        (829, 0),
    )
    check(
        "binary multiples of the weight outgrow the ceiling",
        ([binary_multiples(w, 24) for w in (4, 10, 28, 82)], FIB[25] - 1),
        ([4196351, 1683971, 613817, 228519], 75024),
    )
    profile = state_profile(1, 9, 12)
    fails = 0
    for z1 in range(1, 121):
        for z2 in range(z1, 241):
            if gcd(z1, z2) != 1:
                continue
            h = state_profile(z1, z2, 20)
            if any(h[m] > h[m - 1] + h[m - 2] for m in range(2, 21)):
                fails += 1
    check(
        "the state maximum breaks the Fibonacci recursion",
        (profile[:7], profile[4] > profile[3] + profile[2], fails),
        ([1, 1, 1, 2, 4, 6, 9], True, 8),
    )


# LAW 6: THE PROVED CASES ON THE ADVERSARIAL FAMILIES


def law_families():
    bin3 = [sum(((u >> i) & 1) * 3**i for i in range(8)) for u in range(1, 256)]
    nad = no_adjacent(8)
    families = [
        [(a, b) for a in bin3 for b in bin3 if a < b],
        [(a, b) for a in nad for b in nad if a < b],
        [(1, t) for t in range(2, 3000)],
        [(a, a + 1) for a in range(1, 1500)],
        [(a, 3 * a - 1) for a in range(1, 1200)],
        [(a, 3 * a + 1) for a in range(1, 1200)],
    ]
    seen = set()
    multiplicity = 0
    for family in families:
        for a, b in family:
            if gcd(a, b) != 1:
                continue
            multiplicity += 1
            seen.add((a, b))
    inside = {(a, b) for (a, b) in seen if a <= 120 and b <= 240}
    zero = nodouble = shift = 0
    left = []
    for a, b in sorted(seen - inside):
        order, edges = direction_live(a, b)
        if len(order) == 1:
            zero += 1
        elif not double_branch(order, edges):
            nodouble += 1
        elif is_shift_ray(a, b):
            shift += 1
        else:
            left.append((a, b))
    check(
        "the six families overlap, and their union splits",
        (
            multiplicity,
            len(seen),
            len(inside),
            len(seen) - len(inside),
            zero,
            nodouble,
            shift,
            len(left),
        ),
        (24088, 23435, 717, 22718, 20945, 1693, 3, 77),
    )
    breaches = 0
    wn, wd = 0, 1
    for a, b in left:
        m = direction_mass(a, b, 60)
        for j in range(1, 61):
            top = FIB[j + 1] - 1
            if m[j] > top:
                breaches += 1
            if top and m[j] * wd > wn * top:
                wn, wd = m[j], top
    check(
        "the directions left to the enumeration hold to n = 60",
        (breaches, -(-wn * 10000 // wd)),
        (0, 1516),
    )


# LAW 7: THE GOLDEN POTENTIAL


def qnorm(p, q, d):
    if d < 0:
        p, q, d = -p, -q, -d
    g = gcd(gcd(abs(p), abs(q)), d)
    if g > 1:
        p //= g
        q //= g
        d //= g
    return (p, q, d)


def qadd(x, y):
    return qnorm(x[0] * y[2] + y[0] * x[2], x[1] * y[2] + y[1] * x[2], x[2] * y[2])


def qsub(x, y):
    return qnorm(x[0] * y[2] - y[0] * x[2], x[1] * y[2] - y[1] * x[2], x[2] * y[2])


def qmul(x, y):
    return qnorm(
        x[0] * y[0] + x[1] * y[1], x[0] * y[1] + x[1] * y[0] + x[1] * y[1], x[2] * y[2]
    )


def qinv(x):
    p, q, d = x
    n = p * p + p * q - q * q
    return qnorm(d * (p + q), -d * q, n)


def qsgn(x):
    p, q, _d = x
    hi, lo = 2 * p + q, q
    if hi >= 0 and lo >= 0:
        return 0 if hi == 0 and lo == 0 else 1
    if hi <= 0 and lo <= 0:
        return 0 if hi == 0 and lo == 0 else -1
    s = hi * hi - 5 * lo * lo
    if hi > 0:
        return 1 if s > 0 else (0 if s == 0 else -1)
    return -1 if s > 0 else (0 if s == 0 else 1)


QZERO = (0, 0, 1)
QONE = (1, 0, 1)
QPHI = (0, 1, 1)
QINVPHI = (-1, 1, 1)
QINVPHI2 = (2, -1, 1)


def golden_potential(a, b):
    order, edges = direction_live(a, b)
    n = len(order)
    if n == 1:
        return None
    m = n - 1
    rows = [[QZERO] * (m + 1) for _ in range(m)]
    for r in range(m):
        rows[r][r] = QPHI
        for eps, j in edges[r + 1]:
            if j == 0:
                rows[r][m] = qadd(rows[r][m], QONE)
            else:
                rows[r][j - 1] = qsub(rows[r][j - 1], QONE)
    for col in range(m):
        piv = None
        for r in range(col, m):
            if qsgn(rows[r][col]):
                piv = r
                break
        if piv is None:
            return "singular"
        rows[col], rows[piv] = rows[piv], rows[col]
        scale = qinv(rows[col][col])
        rows[col] = [qmul(v, scale) if qsgn(v) else QZERO for v in rows[col]]
        for r in range(m):
            if r != col and qsgn(rows[r][col]):
                f = rows[r][col]
                rows[r] = [
                    qsub(rows[r][k], qmul(f, rows[col][k]))
                    if qsgn(rows[col][k])
                    else rows[r][k]
                    for k in range(m + 1)
                ]
    u = [QONE] + [rows[r][m] for r in range(m)]
    return u, order, edges


def potential_value(u, edges):
    tot = QZERO
    for eps, j in edges[0]:
        if j:
            tot = qadd(tot, u[j])
    return tot


def potential_valid(u, edges):
    if u[0] != QONE or any(qsgn(v) <= 0 for v in u):
        return False
    for i in range(1, len(edges)):
        s = QZERO
        for eps, j in edges[i]:
            s = qadd(s, u[j])
        if qsgn(qsub(qmul(QPHI, u[i]), s)) < 0:
            return False
    return True


def stress_directions():
    out = set()
    for j in range(1, 9):
        out.add((1, 3**j))
    for i in range(7):
        for j in range(8):
            a, b = 3**i, 3**i + 3**j
            g = gcd(a, b)
            a, b = a // g, b // g
            if b > a:
                out.add((a, b))
    for k in range(1, 7):
        for m in range(1, 40):
            for a, b in (
                (m, 3**k * m + 1),
                (1, 3**k * m),
                (3**k, 3**k + m),
                (3**k, 3**k + 3 * m + 1),
            ):
                if a > b:
                    a, b = b, a
                if a >= 1 and b > a and b <= 6000 and gcd(a, b) == 1:
                    out.add((a, b))
    return out


def law_partition():
    box = []
    for z1 in range(1, 121):
        for z2 in range(z1, 241):
            if gcd(z1, z2) == 1:
                box.append((z1, z2))
    occupied = certified = passing = 0
    over = []
    values = {}
    for a, b in box:
        r = golden_potential(a, b)
        if r is None:
            continue
        occupied += 1
        u, order, edges = r
        if potential_valid(u, edges):
            certified += 1
        tot = potential_value(u, edges)
        values.setdefault(tot, []).append((a, b))
        if qsgn(qsub(QINVPHI2, tot)) >= 0:
            passing += 1
        else:
            over.append((a, b))
    check(
        "the golden potential certifies the box away from the shift rays",
        (occupied, certified, passing, over),
        (218, 218, 214, [(1, 3), (1, 9), (1, 27), (1, 81)]),
    )
    ranked = []
    for _ in range(3):
        top = None
        for k in values:
            if k not in ranked and (top is None or qsgn(qsub(k, top)) > 0):
                top = k
        ranked.append(top)
    check(
        "the top of the golden potential over the box",
        ([(k, sorted(values[k])) for k in ranked], len(values)),
        (
            [
                ((-1, 1, 1), [(1, 3), (1, 9), (1, 27), (1, 81)]),
                ((2, -1, 1), [(1, 12), (3, 10), (4, 9)]),
                ((-14, 9, 2), [(1, 90), (9, 82), (10, 81)]),
            ],
            57,
        ),
    )
    tail = peak = 0
    for a, b in box:
        r = golden_potential(a, b)
        if r is None:
            continue
        u, order, edges = r
        tot = potential_value(u, edges)
        f = first_returns(a, b, 46)
        s = QZERO
        p = QONE
        for j in range(1, 47):
            p = qmul(p, QINVPHI)
            if j >= 2 and f[j]:
                s = qadd(s, qmul((f[j], 0, 1), p))
        if f[1] == 1 and qsgn(qsub(qmul(QINVPHI, tot), s)) >= 0:
            tail += 1
        if qsgn(qsub(QINVPHI, tot)) >= 0:
            peak += 1
    check(
        "the potential dominates the first-return series and no direction "
        "beats the shift rays",
        (tail, peak),
        (218, 218),
    )
    seen = set(box)
    bin3 = [sum(((u >> i) & 1) * 3**i for i in range(8)) for u in range(1, 256)]
    nad = no_adjacent(8)
    families = [
        [(a, b) for a in bin3 for b in bin3 if a < b],
        [(a, b) for a in nad for b in nad if a < b],
        [(1, t) for t in range(2, 3000)],
        [(a, a + 1) for a in range(1, 1500)],
        [(a, 3 * a - 1) for a in range(1, 1200)],
        [(a, 3 * a + 1) for a in range(1, 1200)],
    ]
    for family in families:
        for a, b in family:
            if gcd(a, b) == 1:
                seen.add((a, b))
    seen |= stress_directions()
    total = len(seen)
    occupied = certified = passing = states = 0
    over = []
    gap = []
    best = QZERO
    attain = []
    for a, b in sorted(seen):
        r = golden_potential(a, b)
        if r is None:
            continue
        occupied += 1
        u, order, edges = r
        states = max(states, len(order))
        if potential_valid(u, edges):
            certified += 1
        tot = potential_value(u, edges)
        if qsgn(qsub(QINVPHI2, tot)) >= 0:
            passing += 1
            if qsgn(qsub(tot, best)) > 0:
                best, attain = tot, [(a, b)]
            elif tot == best:
                attain.append((a, b))
        else:
            over.append((a, b, tot))
            if tot != QINVPHI:
                gap.append((a, b, tot))
    check(
        "the golden potential over the box, the six families and the "
        "high-valuation stress list",
        (total, occupied, certified, states, passing),
        (36037, 1995, 1995, 256, 1987),
    )
    check(
        "every direction over the criterion is a shift ray at exactly one over phi",
        (over, gap),
        ([(1, 3**j, QINVPHI) for j in range(1, 9)], []),
    )
    check(
        "the criterion is attained exactly on the supergolden directions",
        (best, attain),
        (QINVPHI2, [(1, 12), (3, 10), (4, 9)]),
    )
    zero = residue = 0
    for z1, z2 in box:
        if z1 % 3 and z2 % 3:
            residue += 1
        if len(direction_live(z1, z2)[0]) == 1:
            zero += 1
    check(
        "occupancy needs three to divide one coordinate",
        (len(box), residue, zero, len(box) - zero),
        (13158, 6566, 12940, 218),
    )


# LAW 8: THE DEGREE POTENTIAL


def split3(a, b):
    q, p = (a, b) if a % 3 == 0 else (b, a)
    k = 0
    q1 = q
    while q1 % 3 == 0:
        q1 //= 3
        k += 1
    return p, k, q1


def phipow(e):
    r = QONE
    for _ in range(abs(e)):
        r = qmul(r, QPHI if e > 0 else QINVPHI)
    return r


def degree_potential(edges):
    return [QONE if len(o) == 2 else QINVPHI for o in edges]


def super_solution(edges, pi):
    for i in range(1, len(edges)):
        s = QZERO
        for eps, j in edges[i]:
            s = qadd(s, pi[j])
        if qsgn(qsub(qmul(QPHI, pi[i]), s)) < 0:
            return False
    return True


def sweep(edges, d):
    pi = degree_potential(edges)
    pi[0] = QONE
    for _ in range(d):
        nxt = [QONE] + [QZERO] * (len(edges) - 1)
        for i in range(1, len(edges)):
            s = QZERO
            for eps, j in edges[i]:
                s = qadd(s, pi[j])
            nxt[i] = qmul(QINVPHI, s)
        pi = nxt
    tot = QZERO
    for eps, j in edges[0]:
        if j:
            tot = qadd(tot, pi[j])
    return tot


def base3(x):
    d = []
    while x:
        d.append(x % 3)
        x //= 3
    return d or [0]


def burst_floor(a, b, k, q1):
    sign = 1 if b % 3 == 0 else -1
    out = []
    for m in range(1, 3**k):
        if m % 3 != 1:
            continue
        x, ok = m, True
        while x:
            if x % 3 > 1:
                ok = False
                break
            x //= 3
        if ok:
            out.append(sign * q1 * m)
    return out


def law_degree():
    seen = set()
    for z1 in range(1, 121):
        for z2 in range(z1, 241):
            if gcd(z1, z2) == 1:
                seen.add((z1, z2))
    bin3 = [sum(((u >> i) & 1) * 3**i for i in range(8)) for u in range(1, 256)]
    nad = no_adjacent(8)
    families = [
        [(a, b) for a in bin3 for b in bin3 if a < b],
        [(a, b) for a in nad for b in nad if a < b],
        [(1, t) for t in range(2, 3000)],
        [(a, a + 1) for a in range(1, 1500)],
        [(a, 3 * a - 1) for a in range(1, 1200)],
        [(a, 3 * a + 1) for a in range(1, 1200)],
    ]
    for family in families:
        for a, b in family:
            if gcd(a, b) == 1:
                seen.add((a, b))
    seen |= stress_directions()
    triple = eligible = occupied = 0
    resbad = burst = burstbad = 0
    valid = nodouble = doublevalid = settled = 0
    depths = {}
    missed = []
    k1 = k1class = 0
    k1bad = []
    quantbad = []
    attain = []
    for a, b in sorted(seen):
        p, k, q1 = split3(a, b)
        if a % 3 == 0 or b % 3 == 0:
            triple += 1
            if (q1 - p) % 3 == 0:
                eligible += 1
        order, edges = direction_live(a, b)
        if len(order) == 1:
            continue
        occupied += 1
        if (q1 - p) % 3:
            resbad += 1
        f = first_returns(a, b, max(k, 2))
        if any(f[j] for j in range(2, k + 1)):
            burst += 1
        if not double_branch(order, edges):
            nodouble += 1
        pi = degree_potential(edges)
        if super_solution(edges, pi):
            valid += 1
            if double_branch(order, edges):
                doublevalid += 1
            hit = None
            for d in range(25):
                if qsgn(qsub(QINVPHI2, sweep(edges, d))) >= 0:
                    hit = d
                    break
            if hit is None:
                missed.append((a, b))
            else:
                settled += 1
                depths[hit] = depths.get(hit, 0) + 1
        else:
            missed.append((a, b))
        u = golden_potential(a, b)[0]
        tot = potential_value(u, edges)
        place = {c: i for i, c in enumerate(order)}
        floor = burst_floor(a, b, k, q1)
        rung = QZERO
        for c in floor:
            if c in place:
                rung = qadd(rung, u[place[c]])
        if len(floor) != 2 ** (k - 1) or qmul(phipow(-(k - 1)), rung) != tot:
            burstbad += 1
        if 2 * p <= 3**k * q1:
            c = p if b % 3 == 0 else -p
            if c not in place or u[place[c]] != QINVPHI:
                burstbad += 1
        if k != 1:
            continue
        k1 += 1
        t = valuation3(q1 - p) if q1 != p else None
        if t is None:
            continue
        bound = qmul(QINVPHI, qsub(QONE, phipow(-max(t, 2))))
        if qsgn(qsub(bound, tot)) < 0:
            quantbad.append((a, b, t))
        if tot == bound:
            attain.append((a, b))
        if t <= 2:
            k1class += 1
            if qsgn(qsub(QINVPHI2, tot)) < 0:
                k1bad.append((a, b, t))
    check(
        "occupancy needs the residue match q1 = p mod 3",
        (triple, eligible, occupied, resbad),
        (20193, 15556, 1995, 0),
    )
    check(
        "no first return has length between two and v3(q)", (occupied, burst), (1995, 0)
    )
    check(
        "the burst identity and the value at the near predecessor",
        (occupied, burstbad),
        (1995, 0),
    )
    check(
        "the degree potential is a super-solution beyond the branch case",
        (occupied, nodouble, valid, doublevalid),
        (1995, 1902, 1968, 66),
    )
    check(
        "the swept degree potential settles all but the shift rays and "
        "twenty-one directions",
        (
            settled,
            sorted(depths.items()),
            [z for z in missed if is_shift_ray(*z)],
            [z for z in missed if not is_shift_ray(*z)],
        ),
        (
            1966,
            [(1, 1804), (3, 101), (4, 44), (5, 11), (6, 6)],
            [(1, 3**j) for j in range(1, 9)],
            [
                (1, 756),
                (1, 2196),
                (1, 2214),
                (1, 2268),
                (1, 2430),
                (9, 2188),
                (10, 2187),
                (13, 1080),
                (13, 3267),
                (27, 730),
                (27, 2188),
                (28, 729),
                (28, 2187),
                (40, 1053),
                (81, 2188),
                (82, 2187),
                (91, 2214),
                (121, 3159),
                (243, 2188),
                (244, 2187),
                (819, 2539),
            ],
        ),
    )
    check(
        "the golden partition bound at v3(q) = 1",
        (k1, k1class, k1bad, quantbad, attain),
        (757, 512, [], [], [(1, 12), (3, 10)]),
    )
    short = []
    for a in range(1, 260):
        for b in range(1, 1100):
            if gcd(a, b) != 1 or (a % 3 and b % 3):
                continue
            if len(direction_live(a, b)[0]) == 1:
                continue
            f = first_returns(a, b, 4)
            if f[2] or f[3]:
                short.append((min(a, b), max(a, b), f[2], f[3]))
    check(
        "the first returns of length two and three are classified",
        sorted(set(short)),
        [(1, 3, 1, 0), (1, 9, 0, 1), (1, 12, 0, 1), (3, 10, 0, 1), (4, 9, 0, 1)],
    )
    tested = bad = 0
    quartic = qmul(qmul(QPHI, QPHI), qmul(QPHI, QPHI))
    for a in range(1, 40):
        for b in range(1, 120):
            if gcd(a, b) != 1 or (a % 3 and b % 3):
                continue
            order, edges = direction_live(a, b)
            if len(order) == 1:
                continue
            tested += 1
            tot = potential_value(golden_potential(a, b)[0], edges)
            den = qsub(QINVPHI2, qmul(QINVPHI, tot))
            closed = qinv(den) if qsgn(den) > 0 else None
            mass = direction_mass(a, b, 46)
            partial = QZERO
            w = QONE
            for n in range(47):
                partial = qadd(partial, qmul((mass[n] + 1, 0, 1), w))
                w = qmul(w, QINVPHI)
            layer = QZERO
            for m in range(1, 3**11 // (a + b)):
                x, y = base3(a * m), base3(b * m)
                if (
                    max(x) < 2
                    and max(y) < 2
                    and not any(dx and dy for dx, dy in zip(x, y))
                ):
                    layer = qadd(layer, phipow(-len(base3((a + b) * m))))
            ok = closed is None or qsgn(qsub(closed, partial)) >= 0
            ok = ok and (closed is not None and qsgn(qsub(quartic, closed)) >= 0) == (
                qsgn(qsub(QINVPHI2, tot)) >= 0
            )
            ok = ok and (qsgn(qsub(QPHI, layer)) >= 0 or qsgn(qsub(tot, QINVPHI2)) > 0)
            bad += not ok
    check(
        "the criterion restated as a golden series and as a multiplier count",
        (tested, bad),
        (111, 0),
    )


def main():
    law_tensor()
    law_box()
    law_weights()
    law_deep()
    law_ceiling()
    law_families()
    law_partition()
    law_degree()
    print("gasket witness weights: every law holds")


main()
