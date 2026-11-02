import json
import os
import subprocess
import sys
import time
from math import ceil, comb, floor, isqrt, log

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", ".."))
DATA_DIR = os.path.join(REPO, "data", os.path.relpath(HERE, REPO))

DENSITY = 0.7

# DESIGN


def digits(base, code):
    return [
        (i, j) for i in range(base) for j in range(base) if code >> (base * i + j) & 1
    ]


def level(base, code, n):
    F = digits(base, code)
    cells = {(0, 0)}
    for _ in range(n):
        cells = {(base * r + i, base * c + j) for (r, c) in cells for (i, j) in F}
    return cells


def signed(base, code):
    return sum((-1) ** (i + j) for (i, j) in digits(base, code))


def imbalance_law(base, code, n):
    s, fill = signed(base, code), len(digits(base, code))
    return s**n if base % 2 else fill ** (n - 1) * s


def neighbours(cells, cell):
    r, c = cell
    return [q for q in ((r + 1, c), (r - 1, c), (r, c + 1), (r, c - 1)) if q in cells]


# MATCHING


def has_perfect_matching(cells):
    black = [p for p in cells if (p[0] + p[1]) % 2 == 0]
    if 2 * len(black) != len(cells):
        return False
    mate = {}
    for root in black:
        seen = set()
        stack = [(root, iter(neighbours(cells, root)))]
        found = False
        while stack and not found:
            _u, it = stack[-1]
            advanced = False
            for v in it:
                if v in seen:
                    continue
                seen.add(v)
                if v not in mate:
                    w = v
                    for x, _ in reversed(stack):
                        nxt = mate.get(x)
                        mate[w] = x
                        mate[x] = w
                        w = nxt
                    found = True
                    break
                stack.append((mate[v], iter(neighbours(cells, mate[v]))))
                advanced = True
                break
            if not advanced and not found:
                stack.pop()
        if not found:
            return False
    return True


def brute_count(cells, side):
    states = {0: 1}
    for p in range(side * side):
        r, c = divmod(p, side)
        nxt = {}
        for mask, ways in states.items():
            covered = mask & 1
            rest = mask >> 1
            if (r, c) not in cells or covered:
                if (r, c) not in cells and covered:
                    continue
                nxt[rest] = nxt.get(rest, 0) + ways
                continue
            if c + 1 < side and (r, c + 1) in cells and not rest & 1:
                key = rest | 1
                nxt[key] = nxt.get(key, 0) + ways
            if (r + 1, c) in cells:
                key = rest | 1 << (side - 1)
                nxt[key] = nxt.get(key, 0) + ways
        states = nxt
    return states.get(0, 0)


# KASTELEYN


def kasteleyn(cells, side, marked=None, holes=True):
    black = sorted(p for p in cells if (p[0] + p[1]) % 2 == 0)
    white = sorted(p for p in cells if (p[0] + p[1]) % 2 == 1)
    if len(black) != len(white):
        return None
    col = {p: k for k, p in enumerate(white)}
    below = {}
    for c in range(side):
        count = 0
        for r in range(side - 1, -1, -1):
            below[(r, c)] = count
            if (r, c) not in cells and holes:
                count += 1
    rows = []
    for b in black:
        row = {}
        for w in neighbours(cells, b):
            if w[0] == b[0]:
                left = min(b, w, key=lambda q: q[1])
                sign = (-1) ** below[left]
            else:
                sign = (-1) ** b[1]
            row[col[w]] = (sign, bool(marked and marked(b, w)))
        rows.append(row)
    return rows


def gp_matrix(rows, value="x"):
    triples = []
    for i, row in enumerate(rows):
        for k, (sign, mark) in row.items():
            entry = f"{sign}*{value}" if mark else str(sign)
            triples.append(f"[{i + 1},{k + 1},{entry}]")
    n = len(rows)
    return f"E=[{','.join(triples)}];M=matrix({n},{n});for(k=1,#E,M[E[k][1],E[k][2]]=E[k][3]);"


def gp(script):
    out = subprocess.run(
        ["gp", "-q", "-D", "parisizemax=4000000000", "-D", "threadsizemax=4000000000"],
        input=script,
        capture_output=True,
        text=True,
        check=True,
    )
    return out.stdout.split()


def kasteleyn_count(cells, side, holes=True):
    rows = kasteleyn(cells, side, holes=holes)
    if rows is None:
        return 0
    if not rows:
        return 1
    return abs(int(gp(gp_matrix(rows) + "print(matdet(M));")[0]))


# STRUCTURE


def components(cells):
    seen, out = set(), []
    for p in cells:
        if p in seen:
            continue
        comp, k = [p], 0
        seen.add(p)
        while k < len(comp):
            for q in neighbours(cells, comp[k]):
                if q not in seen:
                    seen.add(q)
                    comp.append(q)
            k += 1
        out.append(comp)
    return out


def sealed_cells(F, side):
    wide = any((i, j + 1) in F for (i, j) in F)
    tall = any((i + 1, j) in F for (i, j) in F)
    top = side - 1
    for C in components(F):
        if sum((-1) ** (i + j) for (i, j) in C) == 0:
            continue
        reach = any(
            (wide and j == top and (i, 0) in F)
            or (wide and j == 0 and (i, top) in F)
            or (tall and i == top and (0, j) in F)
            or (tall and i == 0 and (top, j) in F)
            for (i, j) in C
        )
        if not reach:
            return True
    return False


def sealed(base, code, m=1):
    return sealed_cells(level(base, code, m), base**m)


def one_way(base, code):
    F, top = set(digits(base, code)), base - 1
    for flip in (False, True):
        G = {(j, i) for (i, j) in F} if flip else F
        tall = any((i + 1, j) in G for (i, j) in G) and any(
            (0, j) in G and (top, j) in G for j in range(base)
        )
        if not tall:
            return G
    return None


def run_returns(base, code):
    G, top = one_way(base, code), base - 1
    rows = (
        [i for i in range(base) if (i, 0) in G and (i, top) in G]
        if any((i, j + 1) in G for (i, j) in G)
        else []
    )
    subsets = [
        frozenset(r for k, r in enumerate(rows) if mask >> k & 1)
        for mask in range(2 ** len(rows))
    ]
    step = {
        A: [
            B
            for B in subsets
            if has_perfect_matching(G - {(i, 0) for i in A} - {(i, top) for i in B})
        ]
        for A in subsets
    }
    seen, frontier = set(), [frozenset()]
    while frontier:
        nxt = []
        for A in frontier:
            for B in step[A]:
                if not B:
                    return True
                if B not in seen:
                    seen.add(B)
                    nxt.append(B)
        frontier = nxt
    return False


def orbit(base, code):
    F, top, out = digits(base, code), base - 1, set()
    for k in range(8):
        image = []
        for i, j in F:
            for _ in range(k % 4):
                i, j = j, top - i
            image.append((i, top - j) if k >= 4 else (i, j))
        out.add(sum(1 << (base * i + j) for i, j in image))
    return min(out)


def picture(base, code):
    F = set(digits(base, code))
    return " ".join(
        "".join("#" if (i, j) in F else "." for j in range(base)) for i in range(base)
    )


# VERBS

LATE = (
    (5, 15571455),
    (5, 19627890),
    (5, 19757811),
    (5, 19920882),
    (5, 32709486),
    (5, 32715747),
    (5, 32715771),
    (5, 32912238),
    (7, 136308971855667),
)


def first_level(base, code, top, cap):
    for n in range(1, top + 1):
        if len(digits(base, code)) ** n > cap:
            return "cap"
        if has_perfect_matching(level(base, code, n)):
            return n
    return None


def imbalance():
    for base, top in ((2, 6), (3, 4)):
        values, wrong = {}, 0
        for code in range(1, 2 ** (base * base)):
            s = signed(base, code)
            values[s] = values.get(s, 0) + 1
            for n in range(1, top + 1):
                cells = level(base, code, n)
                direct = sum(1 if (r + c) % 2 == 0 else -1 for (r, c) in cells)
                wrong += direct != imbalance_law(base, code, n)
        zero = values.get(0, 0)
        print(
            f"base {base}: {2 ** (base * base) - 1} codes, levels 1..{top}, law against the coordinate count wrong {wrong} times"
        )
        print(
            f"  s = 0 at {zero} codes, C({base * base},{base * base // 2}) - 1 = {comb(base * base, base * base // 2) - 1}"
        )
        print(f"  codes per s: {dict(sorted(values.items()))}")


def census():
    for base, top in ((2, 4), (3, 4), (4, 3)):
        tally, never = {}, []
        for code in range(1, 2 ** (base * base)):
            if signed(base, code) != 0:
                continue
            n0 = first_level(base, code, top, 7000)
            tally[n0] = tally.get(n0, 0) + 1
            if n0 is None:
                never.append(code)
        print(
            f"base {base}: first tileable level over the s = 0 codes, levels 1..{top}: {tally}"
        )
        print(
            f"  untileable through level {top}: {len(never)} codes in {len({orbit(base, c) for c in never})} orbits"
        )
        rest = never
        for m in (1, 2, 3):
            hit = [
                c
                for c in rest
                if len(digits(base, c)) ** m <= 70000 and sealed(base, c, m)
            ]
            rest = [c for c in rest if c not in set(hit)]
            print(
                f"  a sealed unbalanced component at level {m}: {len(hit)} more codes; {len(rest)} codes in {len({orbit(base, c) for c in rest})} orbits left{', e.g. ' + str(rest[:6]) if rest and m == 3 else ''}"
            )
        if base == 3:
            print(f"  base 3 untileable codes: {never}")
        if rest:
            print(
                f"  the {len(rest)} left tile level {top + 1}: {sum(has_perfect_matching(level(base, c, top + 1)) for c in rest)}, carry a sealed unbalanced component at level {top + 1}: {sum(sealed(base, c, top + 1) for c in rest)}"
            )
            ways = [
                c for c in range(1, 2 ** (base * base)) if one_way(base, c) is not None
            ]
            tiling = [c for c in ways if has_perfect_matching(level(base, c, 1))]
            print(
                f"  blocks meeting in one direction only: {len(ways)} codes, the {len(tiling)} that tile level 1 all have a run walk back to the empty set: {all(run_returns(base, c) for c in tiling)}"
            )
            dead = [
                c
                for c in rest
                if one_way(base, c) is not None and not run_returns(base, c)
            ]
            rest = [c for c in rest if c not in set(dead)]
            print(
                f"  one direction only and no run walk back to the empty set: {len(dead)} codes in {len({orbit(base, c) for c in dead})} orbits {sorted({orbit(base, c) for c in dead})}; {len(rest)} codes in {len({orbit(base, c) for c in rest})} orbits left {sorted({orbit(base, c) for c in rest})}"
            )
    print("late codes: no tiling at level 1, a tiling at level 2")
    for base, code in LATE:
        one, two = level(base, code, 1), level(base, code, 2)
        T2 = kasteleyn_count(two, base * base)
        brute = brute_count(two, base * base) if base == 5 else "-"
        print(
            f"  base {base} code {code} [{picture(base, code)}] fill {len(one)} s {signed(base, code)} sealed {sealed(base, code)} T(1) {kasteleyn_count(one, base)} {brute_count(one, base)} match(2) {has_perfect_matching(two)} T(2) {T2} brute {brute} orbit {orbit(base, code)}"
        )


def deficiency(cells):
    mate, size = {}, 0
    for root in (p for p in cells if (p[0] + p[1]) % 2 == 0):
        seen, stack, found = set(), [(root, iter(neighbours(cells, root)))], False
        while stack and not found:
            _u, it = stack[-1]
            advanced = False
            for v in it:
                if v in seen:
                    continue
                seen.add(v)
                if v not in mate:
                    w = v
                    for x, _ in reversed(stack):
                        nxt = mate.get(x)
                        mate[w], mate[x] = x, w
                        w = nxt
                    found = True
                    break
                stack.append((mate[v], iter(neighbours(cells, mate[v]))))
                advanced = True
                break
            if not advanced and not found:
                stack.pop()
        size += found
    return len(cells) - 2 * size


def search(draws=None):
    import random

    plan = draws or ((5, 1500000), (6, 300000), (7, 100000))
    for base, count in plan:
        rng = random.Random(base)
        codes = {
            code
            for code in (
                sum(1 << k for k in range(base * base) if rng.random() < DENSITY)
                for _ in range(count)
            )
            if code and signed(base, code) == 0
        }
        untileable = unsealed = tested3 = third = 0
        late = set()
        for code in sorted(codes):
            one = level(base, code, 1)
            if has_perfect_matching(one):
                continue
            untileable += 1
            if sealed(base, code):
                continue
            unsealed += 1
            d2 = deficiency(level(base, code, 2))
            if d2 == 0:
                late.add(code)
                continue
            if d2 < len(one) * deficiency(one) and len(one) ** 3 <= 16000:
                tested3 += 1
                third += has_perfect_matching(level(base, code, 3))
        orbits = sorted({orbit(base, c) for c in late})
        print(
            f"base {base}: {count} draws, {len(codes)} distinct codes with s = 0, {untileable} untileable at level 1, {unsealed} of those unsealed, {len(late)} tile first at level 2 ({len(orbits)} orbits {orbits}), {tested3} tried at level 3 (fill^3 <= 16000, level-2 deficiency below fill times level-1 deficiency), {third} tile first at level 3"
        )
        if base == 5:
            neighbourhood(
                sorted(set(orbits) | {orbit(b, code) for b, code in LATE if b == 5})
            )


def neighbourhood(seeds, base=5, radius=4):
    from itertools import combinations

    seen, late, tested, third = set(), set(), 0, 0
    for seed in seeds:
        for r in range(1, radius + 1):
            for flip in combinations(range(base * base), r):
                code = seed
                for k in flip:
                    code ^= 1 << k
                key = orbit(base, code)
                if key in seen:
                    continue
                seen.add(key)
                if (
                    signed(base, code) != 0
                    or has_perfect_matching(level(base, code, 1))
                    or sealed(base, code)
                ):
                    continue
                if has_perfect_matching(level(base, code, 2)):
                    late.add(key)
                    continue
                if len(digits(base, code)) ** 3 <= 16000:
                    tested += 1
                    third += has_perfect_matching(level(base, code, 3))
    print(
        f"  neighbourhood of {len(seeds)} late base-{base} orbits to Hamming distance {radius}: {len(seen)} orbits, {len(late)} late orbits, {tested} orbits untileable at levels 1 and 2 with fill^3 <= 16000 tried at level 3, {third} tile there"
    )


CODES = ((2, 15, 6), (3, 63, 4), (3, 495, 4), (3, 255, 4))


def load_store(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def remember(key, T):
    path = os.path.join(DATA_DIR, "counts.json")
    store = load_store(path)
    store[key] = str(T)
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(store, fh)


def cached_count(base, code, n):
    path = os.path.join(DATA_DIR, "counts.json")
    store = load_store(path)
    key = f"{base},{code},{n}"
    if key in store:
        return int(store[key])
    T = kasteleyn_count(level(base, code, n), base**n)
    remember(key, T)
    return T


def rectangle(side):
    script = f"default(realprecision,{side * side // 3 + 50});print(round(prod(j=1,{side // 2},prod(k=1,{side // 2},4*cos(Pi*j/{side + 1})^2+4*cos(Pi*k/{side + 1})^2))));"
    return int(gp(script)[0])


def fibonacci(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def control():
    import random

    rng = random.Random(1)
    agree = textbook = tileable = 0
    for _ in range(300):
        side = rng.randint(2, 7)
        cells = {(r, c) for r in range(side) for c in range(side) if rng.random() < 0.8}
        T = brute_count(cells, side)
        tileable += T > 0
        agree += kasteleyn_count(cells, side) == T and has_perfect_matching(cells) == (
            T > 0
        )
        textbook += kasteleyn_count(cells, side, holes=False) == T
    print(
        f"control: 300 random cell sets, side 2..7, density 0.8: Kasteleyn with the hole signs and the matching test agree with brute force at {agree}, the column signs alone at {textbook}; {tileable} sets tile"
    )
    for base, code, n in ((3, 495, 1), (3, 495, 2)):
        print(
            f"  carpet level {n}: column signs alone give |det K| = {kasteleyn_count(level(base, code, n), base**n, holes=False)}"
        )


def count():
    control()
    for base, code, top in CODES:
        print(f"base {base} code {code} [{picture(base, code)}]")
        for n in range(1, top + 1):
            start = time.time()
            cells = level(base, code, n)
            T = kasteleyn_count(cells, base**n)
            checks = []
            if base**n <= 9:
                checks.append(f"brute {brute_count(cells, base**n) == T}")
            if code == 15:
                checks.append(f"product formula {rectangle(2**n) == T}")
            if code == 63:
                checks.append(
                    f"F(3^n+1)^(2^(n-1)) {fibonacci(3**n + 1) ** (2 ** (n - 1)) == T}"
                )
            root = isqrt(T)
            two = (T & -T).bit_length() - 1
            odd = isqrt(T >> two)
            shown = (
                str(T)
                if len(str(T)) <= 60
                else f"{str(T)[:20]}...{str(T)[-20:]} ({len(str(T))} digits)"
            )
            square = odd * odd == T >> two
            print(
                f"  level {n}: cells {len(cells)}, T = {shown}, square {root * root == T}, T = 2^{two} times {'the square of ' + (str(odd) if len(str(odd)) <= 40 else str(len(str(odd))) + '-digit ' + str(odd)[:12] + '...') if square else 'a non-square'}, {', '.join(checks)}, {time.time() - start:.1f} s"
            )
            remember(f"{base},{code},{n}", T)


def crossing_polynomial(base, code, n):
    cells, unit = level(base, code, n), base ** (n - 1)
    rows = kasteleyn(
        cells,
        base**n,
        lambda p, q: (p[0] // unit, p[1] // unit) != (q[0] // unit, q[1] // unit),
    )
    edges = sum(1 for row in rows for _, mark in row.values() if mark)
    script = (
        gp_matrix(rows, "x")
        + f"v=vector({edges + 1},t,matdet(subst(M,x,t-1)));P=polinterpolate(vector({edges + 1},t,t-1),v);P=P*sign(subst(P,x,1));c=content(P);print(vector(poldegree(P)+1,k,polcoeff(P,k-1)));print(c);print(issquare(P/c));"
    )
    out = gp(script)
    return (
        edges,
        [int(a) for a in "".join(out[:-2]).strip("[]").split(",")],
        int(out[-2]),
        out[-1] == "1",
    )


def cross():
    for base, code, n in (
        (3, 495, 2),
        (3, 495, 3),
        (3, 255, 2),
        (3, 255, 3),
        (2, 15, 2),
        (2, 15, 3),
        (2, 15, 4),
    ):
        start = time.time()
        fill = len(digits(base, code))
        edges, N, content, square = crossing_polynomial(base, code, n)
        T, below, unit = (
            cached_count(base, code, n),
            cached_count(base, code, n - 1),
            cached_count(base, code, 1),
        )
        mean = sum(k * a for k, a in enumerate(N)) / sum(N)
        support = [k for k, a in enumerate(N) if a]
        mode = max(range(len(N)), key=lambda k: N[k])
        print(
            f"base {base} code {code} level {n}: {edges} edges cross the {fill} blocks of level {n - 1}"
        )
        print(
            f"  sum N_k = T({n}) {sum(N) == T}, N_0 = T({n - 1})^{fill} {N[0] == below**fill}, small blocks only T(1)^({fill}^{n - 1}) = {unit**fill ** (n - 1)}"
        )
        print(
            f"  N_k nonzero exactly at the even k from 0 to {support[-1]}: {support == list(range(0, support[-1] + 1, 2))}, mode {mode}, mean {mean:.6f}, mean per crossing edge {mean / edges:.6f}"
        )
        print(f"  the polynomial is its content {content} times a square: {square}")
        if (base, code, n) == (3, 495, 2):
            inner = [0] * 9
            for k in range(5):
                inner[2 * k] += comb(4, k) * 2 ** (4 - k)
            inner[8] += 1
            outer = [
                sum(
                    inner[i] * inner[k - i] for i in range(max(0, k - 8), min(k, 8) + 1)
                )
                for k in range(17)
            ]
            print(f"  P_2(x) = ((x^2 + 2)^4 + x^8)^2: {outer == N}")
        print(
            f"  share of tilings that respect the blocks N_0/T = {N[0] / T:.6e}, T/N_0 = {T / N[0]:.6f}"
        )
        print(
            f"  N_k = {N if len(N) <= 30 else N[:6]} {'' if len(N) <= 30 else '...'}, {time.time() - start:.1f} s"
        )


def exposure(base, code, depth):
    F, top = set(digits(base, code)), base - 1
    step = {}
    for E in range(16):
        out = []
        for i, j in F:
            child = 0
            child |= (
                1
                if (i > 0 and (i - 1, j) not in F)
                or (i == 0 and (E & 1 or (top, j) not in F))
                else 0
            )
            child |= (
                2
                if (i < top and (i + 1, j) not in F)
                or (i == top and (E & 2 or (0, j) not in F))
                else 0
            )
            child |= (
                4
                if (j > 0 and (i, j - 1) not in F)
                or (j == 0 and (E & 4 or (i, top) not in F))
                else 0
            )
            child |= (
                8
                if (j < top and (i, j + 1) not in F)
                or (j == top and (E & 8 or (i, 0) not in F))
                else 0
            )
            out.append(child)
        step[E] = out
    share, history = {15: 1.0}, []
    for _ in range(depth):
        nxt = {}
        for E, p in share.items():
            for child in step[E]:
                nxt[child] = nxt.get(child, 0.0) + p / len(F)
        share = nxt
        history.append(dict(share))
    return history, all(
        step[child][k] == child for E in range(16) for k, child in enumerate(step[E])
    )


def hadamard(share):
    total = 0.0
    for E, p in share.items():
        d = 4 - E.bit_count()
        if d == 0:
            return float("-inf")
        total += p * log(d)
    return total / 4


def contact(base, code):
    F, top = set(digits(base, code)), base - 1
    rows = (
        sum((i, 0) in F and (i, top) in F for i in range(base))
        if any((i, j + 1) in F for i, j in F)
        else 0
    )
    cols = (
        sum((0, j) in F and (top, j) in F for j in range(base))
        if any((i + 1, j) in F for i, j in F)
        else 0
    )
    return max(rows, cols)


def down(x, places=9):
    return f"{floor(x * 10**places) / 10**places:.{places}f}"


def up(x, places=6):
    return f"{ceil(x * 10**places) / 10**places:.{places}f}"


def growth():
    G = 0.915965594177219015054603514932384110774
    phi = (1 + 5**0.5) / 2
    exact = {15: ("G/pi", G / 3.141592653589793), 63: ("log(phi)/2", log(phi) / 2)}
    for base, code, top in CODES:
        fill = len(digits(base, code))
        print(f"base {base} code {code} [{picture(base, code)}] fill {fill}")
        L = [log(cached_count(base, code, n)) / fill**n for n in range(1, top + 1)]
        for n in range(1, top + 1):
            inc = L[n - 1] - L[n - 2] if n > 1 else None
            ratio = (L[n - 1] - L[n - 2]) / (L[n - 2] - L[n - 3]) if n > 2 else None
            print(
                f"  level {n}: log T / cells = {down(L[n - 1])}"
                + (f", increment {inc:.6f}" if inc else "")
                + (f", increment ratio {ratio:.4f}" if ratio else "")
            )
        history, idempotent = exposure(base, code, 400)
        H = [hadamard(history[n - 1]) for n in (1, 2, 4, 10, 20, 40, 60, 400)]
        print(
            f"  each digit's exposure map is idempotent: {idempotent}; Hadamard bound H_n at n = 1, 2, 4, 10, 20, 40, 60, 400: {', '.join(f'{h:.13f}' for h in H)}"
        )
        print(f"  bracket: {down(L[-1], 6)} <= theta <= {up(H[-1] + 1e-9)}")
        q, last = contact(base, code) / fill, (L[-1] - L[-2]) / (L[-2] - L[-3])
        tail = lambda r: (L[-1] - L[-2]) * r / (1 - r)
        print(
            f"  contact ratio q = {contact(base, code)}/{fill}; geometric tail from level {top}: with q {L[-1] + tail(q):.6f}, with the last ratio {L[-1] + tail(last):.6f}"
        )
        if code in exact:
            name, value = exact[code]
            print(f"  closed form {name} = {value:.12f}")
    print(
        f"code 27 [{picture(3, 27)}] has no crossing edge, T(n) = 2^(4^(n-1)): {all(cached_count(3, 27, n) == 2**4 ** (n - 1) for n in (1, 2, 3))}, theta = log(2)/4 = {log(2) / 4:.12f}"
    )


VERBS = {
    "imbalance": imbalance,
    "census": census,
    "search": search,
    "count": count,
    "cross": cross,
    "growth": growth,
}

if __name__ == "__main__":
    for name in sys.argv[1:] or list(VERBS):
        start = time.time()
        print(f"# {name}")
        VERBS[name]()
        print(f"# {name} {time.time() - start:.1f} s")
