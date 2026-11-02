from itertools import product

import numpy as np

CARPET = frozenset([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
NET = frozenset([(1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)])
TREE = frozenset([(0, 0, 0), (0, 0, 1)])


def corners(P, q):
    return np.array(
        [v for v in product(range(q), repeat=3) if (v[0] % 2, v[1] % 2, v[2] % 2) in P],
        dtype=np.int64,
    )


def level_points(F, q, n):
    pts = np.zeros((1, 3), dtype=np.int64)
    for _ in range(n):
        pts = (pts[:, None, :] * q + F[None, :, :]).reshape(-1, 3)
    return pts


def prime_factors(m):
    ps, d = [], 2
    while d * d <= m:
        if m % d == 0:
            ps.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        ps.append(m)
    return ps


def squarefree_upto(limit):
    out = []
    for d in range(1, limit + 1):
        ps, prod = prime_factors(d), 1
        for p in ps:
            prod *= p
        if prod == d:
            out.append((d, (-1) ** len(ps)))
    return out


def squarefree_divisors(m):
    out = [(1, 1)]
    for p in prime_factors(m):
        out = out + [(d * p, -mu) for d, mu in out]
    return out


def height_counts(F, q, n):
    weights = {}
    for x in F.sum(axis=1).tolist():
        weights[x] = weights.get(x, 0) + 1
    top, arr = max(weights), np.array([1], dtype=np.int64)
    for _ in range(n):
        L = len(arr)
        new = np.zeros(q * (L - 1) + top + 1, dtype=np.int64)
        for x, c in weights.items():
            new[x : x + q * L : q] += c * arr
        arr = new
    return arr


def aggregate_local(F, q, n, p):
    vs, cur = [tuple(v) for v in F.tolist()], {(0, 0, 0): 1}
    for _ in range(n):
        nxt = {}
        for (a, b, c), t in cur.items():
            for v in vs:
                k = ((q * a + v[0]) % p, (q * b + v[1]) % p, (q * c + v[2]) % p)
                nxt[k] = nxt.get(k, 0) + t
        cur = nxt
    tot = sum(t for r, t in cur.items() if sum(r) % p == 0)
    return cur.get((0, 0, 0), 0), tot


def slice_divisor_counts(F, q, n, s, ds):
    m = n // 2
    lo, hi = level_points(F, q, m), level_points(F, q, n - m)
    scale = q**m
    slo, shi = lo.sum(axis=1).tolist(), hi.sum(axis=1).tolist()
    lol, hil, out = lo.tolist(), hi.tolist(), {}
    for d in ds:
        tab = {}
        for j, v in enumerate(lol):
            k = (slo[j], v[0] % d, v[1] % d, v[2] % d)
            tab[k] = tab.get(k, 0) + 1
        tot = 0
        for i, v in enumerate(hil):
            need = s - shi[i] * scale
            if need >= 0:
                tot += tab.get(
                    (
                        need,
                        (-scale * v[0]) % d,
                        (-scale * v[1]) % d,
                        (-scale * v[2]) % d,
                    ),
                    0,
                )
        out[d] = tot
    return out


def mobius_and_prime_slices(P, q, n):
    F = corners(P, q)
    pts = level_points(F, q, n)
    s = pts.sum(axis=1)
    smax = 3 * (q**n - 1)
    N = np.bincount(s, minlength=smax + 1)
    A = np.bincount(s[np.gcd.reduce(pts, axis=1) == 1], minlength=smax + 1)
    total = np.zeros(smax + 1, dtype=np.int64)
    for d, mu in squarefree_upto(q**n - 1):
        total += mu * np.bincount(s[(pts % d == 0).all(axis=1)], minlength=smax + 1)
    hidden = max(
        [
            int(N[h] - A[h])
            for h in range(2, smax + 1)
            if N[h] and prime_factors(h) == [h]
        ]
        + [0]
    )
    return int(np.count_nonzero(total[1:] != A[1:])), hidden


def peel_mismatches(P, q, n):
    F = corners(P, q)
    a, b = level_points(F, q, n), level_points(F, q, n - 1)
    smax = 3 * (q**n - 1)
    lhs = np.bincount(a.sum(axis=1)[(a % q == 0).all(axis=1)], minlength=smax + 1)
    Nb = np.bincount(b.sum(axis=1), minlength=3 * (q ** (n - 1) - 1) + 1)
    origin, bad = (0, 0, 0) in P, 0
    for h in range(smax + 1):
        r = int(Nb[h // q]) if (origin and h % q == 0 and h // q < len(Nb)) else 0
        bad += int(lhs[h]) != r
    return bad


def parity_dichotomy(q, n):
    F = corners(TREE, q)
    base = level_points(F, q, n - 1)
    ev = od = seen = evengcd = 0
    for v in F.tolist():
        ch = base * q + np.array(v, dtype=np.int64)
        s, g = ch.sum(axis=1), np.gcd.reduce(ch, axis=1)
        pos = s > 0
        e, o = pos & (s % 2 == 0), pos & (s % 2 == 1)
        ev += int(np.count_nonzero(e))
        od += int(np.count_nonzero(o))
        seen += int(np.count_nonzero(e & (g == 1)))
        evengcd += int(np.count_nonzero(o & (g % 2 == 0)))
    return ev, od, seen, evengcd


def walk_lambdas(P, q):
    F = corners(P, q).tolist()
    return {
        t: sum((-1) ** (t[0] * v[0] + t[1] * v[1] + t[2] * v[2]) for v in F) / len(F)
        for t in product((0, 1), repeat=3)
    }


def central_line(P, q, n):
    F = corners(P, q)
    sc = 3 * (q**n - 1) // 2
    ds = squarefree_divisors(sc)
    cnt = slice_divisor_counts(F, q, n, sc, [d for d, _ in ds])
    ps = prime_factors(sc)
    loc = {p: cnt[p] / cnt[1] for p in ps}
    ind = 1.0
    for p in ps:
        ind *= 1 - loc[p]
    return sc, ps, cnt[1], sum(mu * cnt[d] for d, mu in ds), loc, ind


def order(p, q):
    o, x = 1, q % p
    while x != 1:
        x, o = x * q % p, o + 1
    return o


print(
    "DOMAIN carpet q=3 to n=6, net q=3 to n=7, tree q=3 to n=7, carpet q=4 to n=3, carpet q=5 to n=4"
)
print()
print("EXACT LAWS AT EVERY HEIGHT")
for name, P, q, nmax in [
    ("carpet3", CARPET, 3, 4),
    ("net3", NET, 3, 4),
    ("carpet4", CARPET, 4, 3),
    ("carpet5", CARPET, 5, 3),
]:
    for n in range(1, nmax + 1):
        bad, hidden = mobius_and_prime_slices(P, q, n)
        peel = peel_mismatches(P, q, n) if q in (3, 5) else "not a prime base"
        print(
            f"{name} n={n} mobius mismatches={bad} peel mismatches={peel} "
            f"max hidden on a prime slice={hidden}"
        )
pts = level_points(corners(NET, 3), 3, 7)
print(
    f"net3 n=7 points={len(pts)} with 3 dividing the gcd="
    f"{int(np.count_nonzero(np.gcd.reduce(pts, axis=1) % 3 == 0))}"
)
print()

print("PARITY DICHOTOMY, TREE q=3")
for n in (6, 7):
    ev, od, seen, evengcd = parity_dichotomy(3, n)
    print(
        f"tree3 n={n} even-height points={ev} visible among them={seen} "
        f"odd-height points={od} even gcds among them={evengcd}"
    )
print()

print("AGGREGATED LOCAL FACTOR OVER THE HEIGHTS DIVISIBLE BY p, CARPET q=3 n=6")
F = corners(CARPET, 3)
for p in (5, 7, 11, 13):
    div, tot = aggregate_local(F, 3, 6, p)
    print(f"p={p}: {div}/{tot} = {div / tot:.6f} vs 1/p^2 = {1 / p**2:.6f}")
e3 = level_points(F, 3, 3)
bad = sum(
    aggregate_local(F, 3, 3, p)
    != (
        int(np.count_nonzero((e3 % p == 0).all(axis=1))),
        int(np.count_nonzero(e3.sum(axis=1) % p == 0)),
    )
    for p in (2, 5, 7, 11, 13)
)
print(f"transfer against enumeration at n=3 over p=2,5,7,11,13: mismatches={bad}")
lam = walk_lambdas(CARPET, 3)
num = sum(v**6 for v in lam.values()) / 8 * 20**6
den = (1 + lam[(1, 1, 1)] ** 6) / 2 * 20**6
div, tot = aggregate_local(F, 3, 6, 2)
print(f"p=2 counted {div}/{tot} = {div / tot:.7f}")
print(f"p=2 walk formula {round(num)}/{round(den)} = {num / den:.7f}")
print()

print("CENTRAL SLICE")
for q, nmax in ((3, 6), (5, 4)):
    for n in range(1, nmax + 1):
        sc, ps, N, A, loc, ind = central_line(CARPET, q, n)
        print(
            f"q={q} n={n} s*={sc}={'*'.join(map(str, ps))} N={N} A={A} delta={A / N:.5f} "
            f"independence-product={ind:.5f} locals "
            + " ".join(f"{p}:{loc[p]:.5f}" for p in ps)
        )
print()

print("CENTRAL COUNT AND ITS PEEL, CARPET q=3")
hc = [height_counts(corners(CARPET, 3), 3, n) for n in range(15)]
cen = [int(hc[n][3 * (3**n - 1) // 2]) for n in range(15)]
pel = [int(hc[n - 1][(3**n - 1) // 2]) for n in range(2, 15)]
print("N(s*_n) n=1..8:", ", ".join(str(cen[n]) for n in range(1, 9)))
print("N^(3)(s*_n) n=2..8:", ", ".join(str(pel[n - 2]) for n in range(2, 9)))
for n in range(2, 8):
    direct = slice_divisor_counts(corners(CARPET, 3), 3, n, 3 * (3**n - 1) // 2, [3])[3]
    print(
        f"n={n} peel {direct} vs level {n - 1} at R_n={(3**n - 1) // 2}: {pel[n - 2]} "
        f"equal={direct == pel[n - 2]}"
    )
off = sum(cen[n] != 9 * cen[n - 1] - 12 * cen[n - 2] for n in range(2, 15))
off += sum(pel[i] != 9 * pel[i - 1] - 12 * pel[i - 2] for i in range(2, len(pel)))
print(f"(9,-12) recurrence exceptions on both sequences to n=14: {off}")
print(
    f"peel ratio n=14 = {pel[12] / cen[14]:.9f} vs (sqrt(33)-5)/8 = {(33**0.5 - 5) / 8:.10f}"
)
print()

print("THE CENTRAL BILL")
bad = 0
for q in (3, 5, 7):
    for n in range(1, 9):
        ps = prime_factors(3 * (q**n - 1) // 2)
        bad += 3 not in ps
        bad += (2 in ps) != (q % 4 == 1 or n % 2 == 0)
        for p in ps:
            if p not in (2, 3) and q % p:
                bad += n % order(p, q) != 0
print(f"bill rule exceptions over q=3,5,7 and n=1..8: {bad}")
print(
    f"R_7(3) = {(3**7 - 1) // 2}, foreign bill at q=3 n=7 = "
    f"{[p for p in prime_factors(3 * (3**7 - 1) // 2) if p != 3]}"
)
print(f"2^1092 mod 1093^2 = {pow(2, 1092, 1093**2)}")
