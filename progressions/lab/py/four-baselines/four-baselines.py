import math
import os
import sys
import time
from functools import cache
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
DATA_DIR = os.path.join(REPO, "data", os.path.relpath(HERE, REPO))
T0 = time.time()

LIMIT = 10**12
BASES = range(2, 41)
BEHREND_DCAP = 40
TYPE_CAP = 40000
REACH = 10**18
BEYOND_LEVELS = range(8, 12)
LITERATURE = [(19197041, 60974), (515749566, 524288)]
BINOM_N = 64
LIFT_MAXIMISERS = {
    2: 1,
    3: 2,
    4: 2,
    5: 4,
    6: 4,
    7: 6,
    8: 6,
    9: 12,
    10: 4,
    11: 20,
    12: 24,
    13: 52,
    14: 32,
    15: 120,
    16: 96,
    17: 80,
    18: 30,
    19: 18,
    20: 170,
    21: 282,
    22: 90,
    23: 528,
    24: 444,
    25: 140,
    26: 112,
    27: 648,
    28: 16,
    29: 112,
    30: 288,
    31: 480,
    32: 256,
    33: 4400,
    34: 912,
    35: 54,
    36: 5232,
    37: 40,
    38: 5184,
    39: 40,
    40: 180,
}
MAXIMISERS = {}
MASKS = {}

# SEQUENCES

A065825 = [
    1,
    2,
    4,
    5,
    9,
    11,
    13,
    14,
    20,
    24,
    26,
    30,
    32,
    36,
    40,
    41,
    51,
    54,
    58,
    63,
    71,
    74,
    82,
    84,
    92,
    95,
    100,
    104,
    111,
    114,
    121,
    122,
    137,
    145,
    150,
    157,
    163,
    165,
    169,
    174,
    194,
    204,
    209,
]
A003002 = [
    0,
    1,
    2,
    2,
    3,
    4,
    4,
    4,
    4,
    5,
    5,
    6,
    6,
    7,
    8,
    8,
    8,
    8,
    8,
    8,
    9,
    9,
    9,
    9,
    10,
    10,
    11,
    11,
    11,
    11,
    12,
    12,
    13,
    13,
    13,
    13,
    14,
    14,
    14,
    14,
    15,
    16,
    16,
    16,
    16,
    16,
    16,
    16,
    16,
    16,
    16,
    17,
    17,
    17,
    18,
    18,
    18,
    18,
    19,
    19,
    19,
    19,
    19,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    21,
    21,
    21,
    22,
    22,
    22,
    22,
]
A304884 = [
    1,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    5,
    5,
    6,
    5,
    6,
    6,
    6,
    6,
    7,
    7,
    8,
    8,
    8,
    8,
    8,
    8,
    8,
    8,
    9,
    8,
    10,
    8,
    10,
    9,
]


def r3(n):
    if n > A065825[-1]:
        return None
    return sum(1 for k in A065825 if k <= n)


# BRUTE


def first3ap(values):
    s = set(values)
    v = sorted(s)
    for i, x in enumerate(v):
        for z in v[i + 1 :]:
            if (x + z) % 2 == 0 and (x + z) // 2 in s:
                return (x, (x + z) // 2, z)
    return None


def has3ap(values):
    return first3ap(values) is not None


def r3_brute(n):
    best = 0
    chosen, mark = [], set()

    def go(start):
        nonlocal best
        best = max(best, len(chosen))
        for x in range(start, n):
            if len(chosen) + (n - x) <= best:
                return
            if not any((x + y) % 2 == 0 and (x + y) // 2 in mark for y in chosen):
                chosen.append(x)
                mark.add(x)
                go(x + 1)
                chosen.pop()
                mark.discard(x)

    go(0)
    return best


def digits(n, m):
    out = []
    while n:
        out.append(n % m)
        n //= m
    return out[::-1]


# DESIGNS


def signed_set(F):
    return {f + g - 2 * h for f in F for g in F for h in F}


def closes(b, F):
    S = signed_set(F)
    if has3ap(F):
        return True
    return (
        (b in S and -1 in S)
        or (b in S and -b - 1 in S and 1 in S)
        or (-b in S and 1 in S)
        or (-b in S and b + 1 in S and -1 in S)
    )


def level_set(b, F, level):
    vals = np.zeros(1, dtype=np.int64)
    for i in range(level):
        vals = (
            vals[:, None] + np.array(sorted(F), dtype=np.int64)[None, :] * b**i
        ).ravel()
    return vals


def brute_free(b, F, level=3):
    vals = np.sort(level_set(b, F, level))
    s = vals[:, None] + vals[None, :]
    ok = (s % 2 == 0) & (np.arange(len(vals))[:, None] < np.arange(len(vals))[None, :])
    mid = s[ok] // 2
    return not bool(np.isin(mid, vals).any())


def cyclic_sets(b, target):
    half = [0] * b
    for x in range(b):
        half[(2 * x) % b] |= 1 << x
    full = (1 << b) - 1

    def add(blocked, F, x):
        blocked |= 1 << x
        for f in F + [x]:
            blocked |= (
                (1 << ((2 * x - f) % b)) | (1 << ((2 * f - x) % b)) | half[(x + f) % b]
            )
        return blocked

    found = []

    def go(F, blocked, last):
        if len(F) == target:
            found.append(list(F))
            return
        avail = full & ~blocked & ~((1 << (last + 1)) - 1)
        if len(F) + avail.bit_count() < target:
            return
        while avail:
            low = avail & -avail
            x = low.bit_length() - 1
            avail ^= low
            if len(F) + 1 + avail.bit_count() < target:
                return
            go(F + [x], add(blocked, F, x), x)

    go([0], add(0, [], 0), 0)
    return found


def build_designs():
    designs = {2: (1, [0])}
    MAXIMISERS[2] = [[0]]
    for b in BASES:
        if b == 2:
            continue
        sets = cyclic_sets(b, A304884[b - 1])
        assert len(sets) == LIFT_MAXIMISERS[b], (b, len(sets))
        for F in sets[:: max(1, len(sets) // 4)] + [sets[-1]]:
            assert not closes(b, F), (b, F)
            assert brute_free(b, F), (b, F)
        MAXIMISERS[b] = sets
        designs[b] = (len(sets[0]), sets[0])
    return designs


def masks(b):
    if b not in MASKS:
        sets = MAXIMISERS[b]
        member = np.zeros((len(sets), b), dtype=bool)
        for i, F in enumerate(sets):
            member[i, F] = True
        cum = np.zeros((len(sets), b), dtype=np.int64)
        cum[:, 1:] = np.cumsum(member, axis=1)[:, :-1]
        MASKS[b] = (member, cum)
    return MASKS[b]


def family_count(b, N):
    member, cum = masks(b)
    fill = len(MAXIMISERS[b][0])
    t = digits(N - 1, b)
    total = np.zeros(len(member), dtype=np.int64)
    alive = np.ones(len(member), dtype=bool)
    for i, c in enumerate(t):
        total += np.where(alive, cum[:, c], 0) * fill ** (len(t) - 1 - i)
        alive &= member[:, c]
    total += alive
    j = int(total.argmax())
    return int(total[j]), j


def best_other(N, skip, designs):
    best = (0, 0, None)
    for c in designs:
        if c != skip:
            cnt, j = family_count(c, N)
            if (cnt, c) > (best[0], best[1]):
                best = (cnt, c, MAXIMISERS[c][j])
    return best


# GREEDY


def greedy(N):
    M = N - 1
    t = digits(M, 3)
    total, k = 0, len(t)
    for i, c in enumerate(t):
        total += min(c, 2) * 2 ** (k - 1 - i)
        if c >= 2:
            return total
    return total + 1


# BEHREND

WAYS = {}


def ways_for(d):
    if d not in WAYS:
        lmax = len(digits(REACH - 1, 2 * d - 1))
        width = lmax * (d - 1) ** 2 + 1
        sq = [v * v for v in range(d)]
        cur = np.zeros(width, dtype=np.int64)
        cur[0] = 1
        rows = [cur]
        for _ in range(lmax - 1):
            nxt = np.zeros(width, dtype=np.int64)
            for s in sq:
                nxt[s:] += cur[: width - s]
            rows.append(nxt)
            cur = nxt
        WAYS[d] = rows
    return WAYS[d]


def behrend(N, d):
    t = digits(N - 1, 2 * d - 1)
    L = len(t)
    rows = ways_for(d)
    width = len(rows[0])
    total = np.zeros(width, dtype=np.int64)
    pref, alive = 0, True
    for i, c in enumerate(t):
        for v in range(min(c, d)):
            s = pref + v * v
            total[s:] += rows[L - 1 - i][: width - s]
        if c >= d:
            alive = False
            break
        pref += c * c
    if alive:
        total[pref] += 1
    q = int(total.argmax())
    return int(total[q]), q


def behrend_best(N):
    best = (1, 1, 0)
    for d in range(2, BEHREND_DCAP + 1):
        if 2 * d - 1 > N:
            break
        size, q = behrend(N, d)
        if size > best[0]:
            best = (size, d, q)
    return best


# INDEPENDENT COUNTS


def greedy_term(n):
    out, p = 0, 1
    while n:
        if n & 1:
            out += p
        n >>= 1
        p *= 3
    return out


def greedy_alt(N):
    lo, hi = 0, 1 << 62
    while lo < hi:
        mid = (lo + hi) // 2
        if greedy_term(mid) >= N:
            hi = mid
        else:
            lo = mid + 1
    return lo


def behrend_alt(N, d, q):
    t = digits(N - 1, 2 * d - 1)

    @cache
    def go(i, tight, rem):
        if rem < 0:
            return 0
        if i == len(t):
            return 1 if rem == 0 else 0
        top = min(t[i], d - 1) if tight else d - 1
        return sum(go(i + 1, tight and v == t[i], rem - v * v) for v in range(top + 1))

    return go(0, True, q)


def design_count(b, F, N):
    t = digits(N - 1, b)
    fill, allowed = len(F), set(F)
    total = 0
    for i, c in enumerate(t):
        total += sum(1 for f in F if f < c) * fill ** (len(t) - 1 - i)
        if c not in allowed:
            return total
    return total + 1


# SALEM-SPENCER

BINOM = np.zeros((BINOM_N + 1, BINOM_N + 1), dtype=np.int64)
for n in range(BINOM_N + 1):
    BINOM[n][0] = 1
    for k in range(1, n + 1):
        BINOM[n][k] = BINOM[n - 1][k - 1] + (BINOM[n - 1][k] if k <= n - 1 else 0)

TYPES = {}


def types_for(d, L):
    key = (d, L)
    if key not in TYPES:
        if math.comb(L + d - 1, d - 1) > TYPE_CAP:
            TYPES[key] = None
        else:
            cuts = np.array(
                list(combinations(range(L + d - 1), d - 1)), dtype=np.int64
            ).reshape(-1, d - 1)
            pad = np.concatenate(
                [np.full((len(cuts), 1), -1), cuts, np.full((len(cuts), 1), L + d - 1)],
                axis=1,
            )
            TYPES[key] = np.diff(pad, axis=1) - 1
    return TYPES[key]


def multinomial(U):
    cs = np.cumsum(U, axis=1)
    out = np.ones(len(U), dtype=np.int64)
    for j in range(U.shape[1]):
        out *= BINOM[cs[:, j], U[:, j]]
    return out


def spencer(N, d):
    t = digits(N - 1, 2 * d - 1)
    L = len(t)
    T = types_for(d, L)
    if T is None:
        return None
    total = np.zeros(len(T), dtype=np.int64)
    pref, alive = np.zeros(d, dtype=np.int64), True
    for i, c in enumerate(t):
        for v in range(min(c, d)):
            p = pref.copy()
            p[v] += 1
            U = T - p
            ok = (U >= 0).all(axis=1)
            if ok.any():
                total[ok] += multinomial(U[ok])
        if c >= d:
            alive = False
            break
        pref[c] += 1
    if alive:
        total[(T == pref).all(axis=1)] += 1
    j = int(total.argmax())
    return int(total[j]), tuple(int(x) for x in T[j])


def spencer_best(N):
    best = (1, 1, ())
    for d in range(2, BEHREND_DCAP + 1):
        if 2 * d - 1 > N:
            break
        got = spencer(N, d)
        if got is not None and got[0] > best[0]:
            best = (got[0], d, got[1])
    return best


# CONTROLS


def enumerate_class(N, d, key):
    m = 2 * d - 1
    L = len(digits(N - 1, m))
    classes = {}
    for n in range(N):
        t = digits(n, m)
        t = [0] * (L - len(t)) + t
        if max(t, default=0) >= d:
            continue
        classes.setdefault(key(t, d), []).append(n)
    return classes


def controls(designs):
    for n in range(1, 24):
        assert r3_brute(n) == A003002[n], n
    assert all(r3(n) == A003002[n] for n in range(len(A003002)))
    print(
        f"ok r3: brute force n=1..23 and the A065825 count against A003002 n=0..{len(A003002) - 1}"
    )
    for N in list(range(2, 400)) + [3**7, 3**8 + 1, 4321, 37**7, 515749566]:
        assert greedy(N) == greedy_alt(N), N
        if N < 5000:
            assert greedy(N) == sum(1 for n in range(N) if 2 not in digits(n, 3)), N
    print(
        "ok greedy: digit DP against the binary-to-ternary inverse at N=2..399, 3^7, 3^8+1, 4321, 37^7, 515749566, and against enumeration below 5000"
    )
    witness = first3ap(enumerate_class(60, 3, lambda t, d: sum(t))[2])
    assert witness == (2, 6, 10), witness
    print(
        f"ok digit-sum class is NOT progression-free: d=3, base 5, digit sum 2 holds {witness}"
    )
    tested = 0
    for d, N in [
        (2, 200),
        (3, 137),
        (3, 625),
        (4, 343),
        (4, 1000),
        (5, 600),
        (5, 6561),
    ]:
        sq = enumerate_class(N, d, lambda t, d: sum(x * x for x in t))
        ty = enumerate_class(N, d, lambda t, d: tuple(t.count(v) for v in range(d)))
        size, q = behrend(N, d)
        assert size == max(len(c) for c in sq.values()) and len(sq[q]) == size, (d, N)
        assert not has3ap(sq[q]), (d, N)
        got = spencer(N, d)
        assert got is not None
        assert got[0] == max(len(c) for c in ty.values()), (d, N)
        assert len(ty[got[1]]) == got[0] and not has3ap(ty[got[1]]), (d, N)
        assert all(len(c) <= size for c in ty.values())
        tested += 1
    for d, N in [(6, 37**7), (9, 10**12), (3, 3**20 + 5), (12, 10**11 + 7)]:
        size, q = behrend(N, d)
        assert behrend_alt(N, d, q) == size, (d, N)
        assert all(
            behrend_alt(N, d, k) <= size for k in range(max(0, q - 15), q + 16)
        ), (d, N)
    print(
        f"ok behrend and digit-type classes: DP against enumeration at {tested} (d, N), maximiser progression-free, every type class inside a square-sum shell, memoised recount at 4 large (d, N)"
    )
    for b, (fill, F) in designs.items():
        assert design_count(b, F, b**3) == fill**3, b
        if b >= 3:
            assert fill == A304884[b - 1], b
            for N in (200, 1000, 4097):
                allowed = set(F)
                assert design_count(b, F, N) == sum(
                    1 for n in range(N) if all(x in allowed for x in digits(n, b))
                ), (b, N)
    for b in designs:
        for N in (1000, 4097):
            assert family_count(b, N)[0] == max(
                design_count(b, F, N) for F in MAXIMISERS[b]
            ), (b, N)
    print(
        "ok designs: fill equals A304884 at bases 3..40, the maximiser counts equal the census lift column at bases 2..40, criterion and level-3 brute force agree (free) on a sample of each base, truncated counts against enumeration over every maximiser"
    )


# TABLE


def fmt_set(F):
    return "{" + ", ".join(map(str, F)) + "}"


def row_data(N, reps, designs):
    b, L = max(reps, key=lambda bl: (designs[bl[0]][0] ** bl[1], -bl[0]))
    fill = designs[b][0]
    design = fill**L
    g = greedy(N)
    ss = spencer_best(N)
    bh = behrend_best(N)
    assert ss[0] <= bh[0], N
    exact = r3(N)
    known = exact if exact is not None else r3(A065825[-1])
    tr = best_other(N, b, designs)
    others = {"greedy": g, "salem-spencer": ss[0], "behrend": bh[0], "r3": known}
    top = max(others.values())
    names = [k for k, v in others.items() if v == max(top, design, tr[0])]
    if tr[0] == max(top, design, tr[0]):
        names.append("other-base")
    if design == max(top, design, tr[0]):
        names.append("design")
    return {
        "N": N,
        "b": b,
        "L": L,
        "fill": fill,
        "design": design,
        "greedy": g,
        "ss": ss,
        "beh": bh,
        "r3": exact,
        "tr": tr,
        "names": names,
        "win": design > top,
        "family": design <= top and tr[0] > top,
        "ratio": design / top,
        "best_other": top,
    }


def line(r):
    r3s = "-" if r["r3"] is None else str(r["r3"])
    mark = "*" if r["win"] else "+" if r["family"] else " "
    return (
        f"{r['N']:>14}{mark}{str(r['b']) + '^' + str(r['L']):>9} {r['fill']:>3} {r['design']:>12} {r['greedy']:>12} "
        f"{r['ss'][1]:>4} {r['ss'][0]:>12} {r['beh'][1]:>3} {r['beh'][2]:>5} {r['beh'][0]:>12} {r['tr'][1]:>3} {r['tr'][0]:>12} {r3s:>4} {r['ratio']:>6.3f}  {'='.join(r['names'])}"
    )


def classical(N):
    return max(greedy(N), behrend_best(N)[0], r3(min(N, A065825[-1])))


def interval(b, F, L):
    size = len(F) ** L
    lo = max(F) * (b**L - 1) // (b - 1) + 1
    if lo < 2 or classical(lo) >= size:
        return None
    hi = lo
    while classical(2 * hi) < size:
        hi *= 2
    top = 2 * hi
    while top - hi > 1:
        mid = (hi + top) // 2
        if classical(mid) < size:
            hi = mid
        else:
            top = mid
    return lo, hi


HEAD = f"{'N':>14} {'design':>9} {'|F|':>3} {'size':>12} {'greedy':>12} {'SS d':>4} {'SS size':>12} {'B d':>3} {'B q':>5} {'B size':>12} {'ob':>3} {'other base':>12} {'r3':>4} {'ratio':>6}  WINNER"


def main():
    designs = build_designs()
    controls(designs)
    print()
    print(
        "DESIGNS: base, digits, dimension, the least set found with 0 (level-3 brute force free)"
    )
    for b in BASES:
        fill, F = designs[b]
        print(
            f"  b={b:>2} |F|={fill:>2} dim={math.log(fill) / math.log(b):.4f} F={fmt_set(F)}"
        )
    reach = {}
    for b in BASES:
        L = 1
        while b**L <= LIMIT:
            reach.setdefault(b**L, []).append((b, L))
            L += 1
    rows = [row_data(N, reach[N], designs) for N in sorted(reach)]
    print()
    print(
        f"TABLE: {len(rows)} values of N = b^L <= {LIMIT:.0e}, b = 2..40; * marks the design |F|^L strictly above greedy, SS, Behrend and r3 (exact to N = {A065825[-1]}, the value at {A065825[-1]} above); + marks a row the design misses while another base's design truncated to N is strictly above them; ob = the best other base truncated to N, ratio = design over the best of the four"
    )
    lines = [HEAD] + [line(r) for r in rows]
    for text in lines:
        print(text)
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(os.path.join(DATA_DIR, "four-baselines.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    wins = [r for r in rows if r["win"]]
    fam = [r for r in rows if r["family"]]
    print()
    print(
        f"WINNERS: {len(wins)} rows with the design |F|^L strictly above all four; {len(fam)} further rows where another base truncated to N is"
    )
    by_base = {}
    for r in fam:
        by_base.setdefault(r["tr"][1], []).append(r)
    for c in sorted(by_base):
        group = by_base[c]
        top = max(group, key=lambda r: r["tr"][0] / r["best_other"])
        print(
            f"  + base {c} truncated: {len(group)} rows, N from {group[0]['N']} to {group[-1]['N']}, largest lead over the four x{top['tr'][0] / top['best_other']:.4f} at N={top['N']} ({top['tr'][0]} against {top['best_other']}) F={fmt_set(top['tr'][2])}"
        )
    for r in wins:
        N, size, bsize, d, q = (
            r["N"],
            r["design"],
            r["beh"][0],
            r["beh"][1],
            r["beh"][2],
        )
        assert greedy_alt(N) == r["greedy"] and behrend_alt(N, d, q) == bsize
        e4, e2 = math.log(N) ** 0.25, math.log2(N) ** 0.25
        print(
            f"  N={N} = {r['b']}^{r['L']} F={fmt_set(designs[r['b']][1])} design={size}"
        )
        print(
            f"    greedy={r['greedy']} (x{size / r['greedy']:.4f}) SS={r['ss'][0]} d={r['ss'][1]} (x{size / r['ss'][0]:.2f}) Behrend={bsize} d={d} q={q} (x{size / bsize:.3f}) other base={r['tr'][0]} b={r['tr'][1]} (x{size / r['tr'][0]:.3f})"
        )
        print(
            f"    Elkin: (ln N)^(1/4)={e4:.3f}, (log2 N)^(1/4)={e2:.3f}; Behrend x (ln N)^(1/4) = {bsize * e4:.0f}, x (log2 N)^(1/4) = {bsize * e2:.0f}; constant needed to reach the design: {size / (bsize * e4):.3f} (ln), {size / (bsize * e2):.3f} (log2); greedy is untouched by it"
        )
    print()
    print(
        "NEIGHBOURS: the rows on each side of every winner in the table order, and the design's own base across levels"
    )
    index = {r["N"]: i for i, r in enumerate(rows)}
    shown = set()
    for r in wins:
        i = index[r["N"]]
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(rows) and j not in shown:
                shown.add(j)
                s = rows[j]
                tag = "WIN" if s["win"] else "family win" if s["family"] else "lose"
                print(
                    f"  N={s['N']} {s['b']}^{s['L']} design={s['design']} best of the four={s['best_other']} ratio={s['ratio']:.4f} other base={s['tr'][0]} (b={s['tr'][1]}) {tag}"
                )
    for b in sorted({r["b"] for r in wins}):
        print(f"  base {b} across levels (design / greedy):")
        for r in rows:
            if r["b"] == b:
                tag = "WIN" if r["win"] else "family win" if r["family"] else "lose"
                print(
                    f"    L={r['L']} N={r['N']} design/greedy={r['design'] / r['greedy']:.4f} design/best of the four={r['ratio']:.4f} {tag}"
                )
        fill, F = designs[b]
        print(
            f"  base {b} beyond the box (design vs greedy, Behrend, digit-type class):"
        )
        for L in BEYOND_LEVELS:
            N = b**L
            g, bh, ss = greedy(N), behrend_best(N), spencer_best(N)
            print(
                f"    L={L} N={N} design={fill**L} greedy={g} (x{fill**L / g:.4f}) Behrend={bh[0]} d={bh[1]} (x{fill**L / bh[0]:.3f}) SS={ss[0]} d={ss[1]}"
            )
    print()
    print(
        "SPANS: a design of level L is a set inside [0, S], S = max F (b^L - 1)/(b - 1), so it sits in [0, N) for every N > S; classical(N) = max(greedy, Behrend, r3 lower bound) is nondecreasing in N, so the design wins on one interval"
    )
    spans = []
    for b in BASES:
        if b == 2:
            continue
        fill, F = designs[b]
        for L in range(1, 100):
            if b**L > LIMIT:
                break
            got = interval(b, F, L)
            if got:
                spans.append((b, L, fill**L, got))
    print(
        f"  {len(spans)} designs (b, L) with b^L <= {LIMIT:.0e} win on a nonempty interval of N"
    )
    for b, L, size, (lo, hi) in spans:
        c_lo, c_hi, c_b = classical(lo), classical(hi), classical(b**L)
        assert hi == greedy_term(size - 1), (b, L)
        print(
            f"  b={b} L={L} size={size} N in [{lo}, {hi}] ({hi - lo + 1} values) holds b^L={b**L}: {lo <= b**L <= hi}; classical at lo={c_lo} at b^L={c_b} at hi={c_hi}; 3^23 inside: {lo <= 3**23 <= hi}"
        )
    print()
    print(
        "LITERATURE POINTS: sets in 1..N of the stated size, from the comment of the A065825 record (construction unread), against the greedy count at the same N"
    )
    for N, size in LITERATURE:
        bh = behrend_best(N)
        fam_best = best_other(N, 0, designs)
        print(
            f"  N={N} set of {size}: greedy={greedy(N)} (set x{size / greedy(N):.4f}) Behrend={bh[0]} d={bh[1]}; best design truncated to N over every maximiser of every base={fam_best[0]} at base {fam_best[1]} F={fmt_set(fam_best[2])} (set x{size / fam_best[0]:.3f})"
        )
    losers = [r for r in rows if r["ratio"] < 1 and r["r3"] is None]
    print()
    print("CLOSEST LOSSES with r3 unknown and ratio below 1 (design / best other)")
    for r in sorted(losers, key=lambda r: -r["ratio"])[:8]:
        print(
            f"  N={r['N']} = {r['b']}^{r['L']} ratio={r['ratio']:.4f} winner={'='.join(r['names'])}"
        )
    ties = [r["N"] for r in rows if r["ratio"] == 1]
    print(
        f"TIES: {len(ties)} values of N at ratio exactly 1, all powers of 3: {all(3 ** round(math.log(N, 3)) == N for N in ties)}"
    )
    caps = [r for r in rows if r["beh"][1] == BEHREND_DCAP]
    print(
        f"CAPS: Behrend optimum at the d={BEHREND_DCAP} cap in {len(caps)} rows; type classes capped at {TYPE_CAP} types per (d, L); digit-type classes searched for d up to {max(r['ss'][1] for r in rows)}"
    )
    print(f"time {time.time() - T0:.1f} s")


if __name__ == "__main__":
    sys.exit(main())
