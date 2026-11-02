import math
import sys
import time
from fractions import Fraction
from math import gcd, isqrt, log

import numpy as np

A030979_URL = (
    "https://raw.githubusercontent.com/oeis/oeisdata/main/seq/A030/A030979.seq"
)

# MEMBERSHIP


def held(p, q, n):
    m = 2 * q
    a = p % m
    seen = set()
    while a not in seen:
        if a > q:
            return False
        seen.add(a)
        a = n * a % m
    return True


def held_ifs(x, n):
    seen = set()
    while x not in seen:
        if x < 0 or x > 1:
            return False
        seen.add(x)
        y = n * x
        d = math.floor(y)
        if d % 2:
            if y != d:
                return False
            d -= 1
        x = y - d
    return True


def residues(p, q):
    return [r for r in range(1, 2 * q, 2) if held(p, q, r)]


def parity_rule(p, q, r):
    s = p % q
    seen = set()
    while s not in seen:
        if s and (s - p) % 2:
            return False
        seen.add(s)
        s = r * s % q
    return True


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def phi(n):
    return sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


def order(r, q):
    k, a = 1, r % q
    while a != 1:
        a = a * r % q
        k += 1
    return k


def primitive_root(q):
    return next(g for g in range(2, q) if order(g, q) == q - 1)


def prime_share(p, q):
    m = q - 1
    while m % 2 == 0:
        m //= 2
    g = primitive_root(q)
    units = 0
    for d in range(1, m + 1):
        if m % d:
            continue
        h = pow(g, (q - 1) // d, q)
        coset = {p * pow(h, i, q) % q for i in range(d)}
        if all((s - p) % 2 == 0 for s in coset):
            units += phi(d)
    return Fraction(1 + units, q)


def least_period(bits):
    n = len(bits)
    return next(
        t
        for t in range(1, n + 1)
        if n % t == 0 and all(bits[i] == bits[i % t] for i in range(n))
    )


def period(top):
    t0 = time.time()
    checks = law_bad = 0
    for q in range(2, top + 1):
        for p in range(q + 1):
            if gcd(p, q) != 1:
                continue
            res = set(residues(p, q))
            last = 6 * q + 1 if q <= 40 else 2 * q + 1
            for n in range(3, last + 1, 2):
                checks += 1
                law_bad += held_ifs(Fraction(p, q), n) != (n % (2 * q) in res)
    print(
        f"law: orbit rule against the digit walk at every p/q in [0, 1], q <= {top}, odd sides 3..6q+1 at q <= 40 and 3..2q+1 above: {checks} checks, {law_bad} failures"
    )
    count = sym_bad = flip_bad = rule_bad = unit_bad = prime_bad = low_bad = (
        high_bad
    ) = 0
    best = []
    low_eq = []
    short = []
    primes_big = []
    for q in range(2, top + 1):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            count += 1
            res = residues(p, q)
            rs = set(res)
            share = Fraction(len(res), q)
            if res != residues(q - p, q):
                sym_bad += 1
            floor_ = Fraction(1, 2) if q == 2 else Fraction(2, q)
            ceil_ = Fraction(q + 1, 2 * q) if q % 2 else Fraction(1, 2)
            low_bad += share < floor_
            high_bad += share > ceil_
            if share == floor_:
                low_eq.append(q)
            best.append((share, p, q))
            bits = [held(p, q, 2 * i + 1) for i in range(q)]
            t = least_period(bits)
            if t < q:
                short.append((p, q, t))
            if q % 2 == 0:
                flip_bad += any(
                    ((q - r) % (2 * q) in rs) != (r in rs) for r in range(1, 2 * q, 2)
                )
            else:
                rule_bad += any(
                    parity_rule(p, q, r) != (r in rs) for r in range(1, 2 * q, 2)
                )
                units = sum(1 for r in res if gcd(r, q) == 1)
                seen = {}
                for r in range(1, q):
                    if gcd(r, q) != 1:
                        continue
                    g = frozenset(pow(r, i, q) for i in range(order(r, q)))
                    seen[g] = all((p * s % q - p) % 2 == 0 for s in g)
                formula = sum(phi(len(g)) for g, ok in seen.items() if ok)
                unit_bad += formula != units
                if is_prime(q):
                    ps = prime_share(p, q)
                    prime_bad += ps != share
                    if p == 1 and ps > Fraction(2, q):
                        primes_big.append(q)
    best.sort(key=lambda r: (-r[0], r[2], r[1]))
    print(f"fractions p/q in (0,1), q <= {top}: {count}")
    print(
        f"symmetry p -> q-p: {sym_bad} failures; even q, r -> q-r: {flip_bad} failures"
    )
    print(
        f"odd q parity rule mod q: {rule_bad} failures; unit formula over cyclic subgroups: {unit_bad} failures; prime formula: {prime_bad} failures"
    )
    print(
        f"lower bound 2/q (1/2 at q=2): {low_bad} failures, attained at {len(low_eq)} fractions, least q {sorted(set(low_eq))[:12]}"
    )
    print(f"upper bound (q+1)/(2q) odd q, 1/2 even q: {high_bad} failures")
    print("largest shares:", ", ".join(f"{p}/{q} {s}" for s, p, q in best[:8]))
    above = [(p, q, s) for s, p, q in best if q > 3 and s > Fraction(1, 2)]
    print(f"shares above 1/2 with q > 3: {len(above)}")
    odd = max((s, -q, p) for s, p, q in best if q % 2 and q >= 5)
    half = sorted({q for s, p, q in best if q % 2 == 0 and s == Fraction(1, 2)})
    print(
        f"largest share at odd q >= 5: {odd[2]}/{-odd[1]} {odd[0]}; even q with some share 1/2: {len(half)}, first {half[:10]}"
    )
    print(
        f"least period in n = (N-1)/2 below q: {len(short)} fractions, first {short[:6]}"
    )
    print(f"odd primes q <= {top} with share(1/q) > 2/q: {primes_big}")
    print(f"runtime {time.time() - t0:.2f} s")


# EISENSTEIN


def eisenstein(top):
    t0 = time.time()
    tests = bad = full_match = 0
    for q in range(3, top + 1):
        if not is_prime(q):
            continue
        for p in range(1, q):
            c = sum(1 for n in range(3, q, 2) if (p * n // q) % 2)
            a = p if p % 2 else q - p
            leg = pow(a, (q - 1) // 2, q)
            tests += 1
            bad += (c % 2 == 1) != (leg == q - 1)
            full = sum(p * n // q for n in range(1, 2 * q, 2))
            assert full == (2 * p - 1) * (q - 1) // 2 + p
            full_match += (full % 2 == 1) == (leg == q - 1)
    print(
        f"odd primes q <= {top}, every p: parity of #(odd 3 <= N < q, first base-N digit of p/q odd) against (a/q), a = p or q-p odd"
    )
    print(f"tests {tests}, failures {bad}")
    print(
        f"sum of floor(pN/q) over odd N < 2q equals (2p-1)(q-1)/2 + p at all {tests}; its parity agrees with the symbol at {full_match}"
    )
    print(f"runtime {time.time() - t0:.2f} s")


# INTEGER COUNT


def in_k(k, n):
    h = (n - 1) // 2
    while k:
        if k % n > h:
            return False
        k //= n
    return True


def counts(top):
    t = isqrt(2 * top) + 1
    c = np.zeros(top + 1, dtype=np.int64)
    for n in range(3, t + 1, 2):
        h = (n - 1) // 2
        digits = np.arange(h + 1, dtype=np.int64)
        vals = digits[digits <= top]
        pw = n
        while pw <= top:
            step = digits * pw
            step = step[step <= top]
            vals = (vals[None, :] + step[:, None]).ravel()
            vals = vals[vals <= top]
            pw *= n
        vals = vals[vals >= (n + 1) // 2]
        c += np.bincount(vals, minlength=top + 1)
    diff = np.zeros(top + 2, dtype=np.int64)
    first = t + 1 if (t + 1) % 2 else t + 2
    j = 2
    while j * first <= 2 * top:
        ns = np.arange(first, 2 * top // j + 1, 2, dtype=np.int64)
        lo = j * ns // 2
        hi = np.minimum(((j + 1) * ns - 1) // 2, top)
        diff += np.bincount(lo, minlength=top + 2)
        diff -= np.bincount(hi + 1, minlength=top + 2)
        j += 2
    return c + np.cumsum(diff)[: top + 1]


def count(top):
    t0 = time.time()
    c = counts(top)
    t1 = time.time()
    small = 3000
    bad = sum(
        1
        for k in range(1, small + 1)
        if c[k] != sum(1 for n in range(3, 2 * k + 1, 2) if in_k(k, n))
    )
    print(
        f"c(k) = #(odd 3 <= N <= 2k : 2k in Z_N) for all k <= {top}: {t1 - t0:.2f} s; direct digit test at k <= {small}: {bad} failures"
    )
    k = np.arange(top + 1, dtype=np.float64)
    err = c - (1 - log(2)) * k
    ratio = np.zeros_like(err)
    ratio[1:] = err[1:] / np.sqrt(k[1:])
    i_max = int(np.argmax(np.abs(ratio[1:]))) + 1
    print(
        f"max over 1 <= k <= {top} of |c(k) - (1 - log 2) k| / sqrt(k): {math.ceil(abs(ratio[i_max]) * 10**6) / 10**6:.6f} at k = {i_max}"
    )
    for lo in [10**2, 10**3, 10**4, 10**5]:
        hi = min(10 * lo, top)
        seg = ratio[lo : hi + 1]
        print(
            f"  k in [{lo}, {hi}]: err/sqrt(k) min {math.floor(seg.min() * 10**6) / 10**6:.6f} max {math.ceil(seg.max() * 10**6) / 10**6:.6f} mean {seg.mean():.6f}"
        )
    for kk in [10**3, 10**4, 10**5, 10**6]:
        if kk <= top:
            print(
                f"  c({kk}) = {c[kk]}, (1 - log 2) k = {(1 - log(2)) * kk:.3f}, err/sqrt(k) = {ratio[kk]:.6f}"
            )
    bound = np.sqrt(2 * k[1:]) + 1
    print(
        f"largest |err| / (sqrt(2k) + 1) over 1 <= k <= {top}: {math.ceil(np.max(np.abs(err[1:]) / bound) * 10**6) / 10**6:.6f}"
    )
    print(f"runtime {time.time() - t0:.2f} s")


def fast(k):
    x = 2 * k
    s = isqrt(x)
    b = sum(1 for n in range(3, s + 1, 2) if in_k(k, n))
    a, j = 0, 2
    while x // j > s:
        lo, hi = max(x // (j + 1), s), x // j
        a += (hi + 1) // 2 - (lo + 1) // 2
        j += 2
    return a + b, b


def second(top):
    t0 = time.time()
    bad = sum(
        1
        for k in range(1, 3001)
        if fast(k)[0] != sum(1 for n in range(3, 2 * k + 1, 2) if in_k(k, n))
    )
    print(f"block count against the direct digit test at k <= 3000: {bad} failures")
    zeta_half = -1.4603545088095868
    kappa = -(2 - math.sqrt(2)) * zeta_half / 4
    beta = (math.sqrt(2) + (2 - math.sqrt(2)) * zeta_half) / 4
    print(
        f"kappa = -(2 - sqrt 2) zeta(1/2)/4 = {kappa:.6f}; beta = (sqrt 2 + (2 - sqrt 2) zeta(1/2))/4 = {beta:.6f}"
    )
    rng = np.random.default_rng(376)
    for e in range(6, top + 1):
        ks = [int(v) for v in rng.integers(10**e, 2 * 10**e, size=60)]
        r = []
        rb = []
        for k in ks:
            v, b = fast(k)
            r.append((v - (1 - log(2)) * k) / math.sqrt(k))
            rb.append(b / math.sqrt(k))
        print(
            f"  60 k in [10^{e}, 2 10^{e}): err/sqrt(k) mean {np.mean(r):.6f} min {min(r):.6f} max {max(r):.6f}; small-side part / sqrt(k) mean {np.mean(rb):.6f}"
        )
    print(f"runtime {time.time() - t0:.2f} s")


# FAMILY


def carries_ok(k, p, a):
    c = 0
    i = 0
    while k or c:
        s = 2 * (k % p) + c
        c = 1 if s >= p else 0
        if c and i % a == a - 1:
            return False
        k //= p
        i += 1
    return True


def family():
    t0 = time.time()
    meet = [
        m
        for m in range(3000)
        if all(m % 2 == 0 and in_k(m // 2, n) for n in range(3, m + 3, 2))
    ]
    print(f"integers below 3000 held by every odd side 3 <= N <= m+1: {meet}")
    miss = [
        (p, q)
        for q in range(2, 61)
        for p in range(1, q)
        if gcd(p, q) == 1 and held(p, q, 2 * q - 1)
    ]
    print(f"p/q in (0,1), q <= 60, held at side 2q-1: {len(miss)}")
    sub_int = all(
        in_k(k, n**e)
        for n in (3, 5, 7)
        for e in (2, 3)
        for k in range(10**5)
        if in_k(k, n)
    )
    sub_rat = all(
        held(p, q, pow(n, e, 2 * q))
        for q in range(2, 61)
        for p in range(q + 1)
        if gcd(p, q) == 1
        for n in range(3, 2 * q + 1, 2)
        for e in (2, 3, 4)
        if held(p, q, n)
    )
    print(
        f"E_N inside E_(N^e): integers k < 10^5 at N = 3, 5, 7, e = 2, 3: {sub_int}; rationals q <= 60, e = 2, 3, 4: {sub_rat}"
    )
    strict = all(
        not held(n + 1, n**e, n)
        and held(n + 1, n**e, n**e)
        and not in_k((n + 1) // 2, n)
        and in_k((n + 1) // 2, n**e)
        for n in range(3, 40, 2)
        for e in (2, 3, 4)
    )
    print(
        f"strict at e >= 2: (N+1)/N^e and N+1 lie at side N^e and not at side N, odd N < 40, e = 2, 3, 4: {strict}"
    )
    kum = all(
        (math.comb(2 * k, k) % p != 0) == in_k(k, p)
        for p in (3, 5, 7, 11, 13, 17, 19, 23)
        for k in range(1500)
    )
    print(
        f"Kummer: K_p = (k : p does not divide C(2k,k)) at odd primes p <= 23, k < 1500: {kum}"
    )
    pw = all(
        in_k(k, p**a) == carries_ok(k, p, a)
        for p, a in ((3, 2), (3, 3), (5, 2), (7, 2))
        for k in range(10**5)
    )
    print(
        f"K_(p^a) = no carry of k + k in base p out of a position = a-1 mod a, (p,a) in (3,2),(3,3),(5,2),(7,2), k < 10^5: {pw}"
    )
    v3 = next(k for k in range(10**5) if in_k(k, 9) and (math.comb(2 * k, k) % 9 == 0))
    print(
        f"least k in K_9 with 9 | C(2k,k): {v3}; least k with 9 not dividing C(2k,k) outside K_9: {next(k for k in range(10**5) if not in_k(k, 9) and math.comb(2 * k, k) % 9)}"
    )
    a = next(k for k in range(10**5) if in_k(k, 3) and in_k(k, 5) and not in_k(k, 15))
    b = next(k for k in range(10**5) if in_k(k, 15) and not in_k(k, 3))
    cc = next(k for k in range(10**5) if in_k(k, 15) and not in_k(k, 5))
    d = next(k for k in range(10**5) if in_k(k, 15) and not (math.comb(2 * k, k) % 15))
    print(
        f"least k in K_3 cap K_5 outside K_15: {a}; least in K_15 outside K_3: {b}; outside K_5: {cc}; least k in K_15 with 15 | C(2k,k): {d}"
    )
    print(f"runtime {time.time() - t0:.2f} s")


# FINITE INTERSECTIONS


def next_in(n, lo):
    h = (n - 1) // 2
    while True:
        x, i, bad = lo, 0, -1
        while x:
            if x % n > h:
                bad = i
            x //= n
            i += 1
        if bad < 0:
            return lo
        p = n ** (bad + 1)
        lo = (lo // p + 1) * p


def common(sides, top):
    b, rest = sides[0], sides[1:]
    h = (b - 1) // 2
    depth = 0
    while b**depth <= top:
        depth += 1
    pw = [b**i for i in range(depth + 1)]
    span = [h * (pw[i] - 1) // (b - 1) for i in range(depth + 1)]
    out, nodes, stack = [], 0, [(0, depth)]
    while stack:
        base, j = stack.pop()
        nodes += 1
        if base > top:
            continue
        hi = min(base + span[j], top)
        if any(next_in(n, base) > hi for n in rest):
            continue
        if j == 0:
            out.append(base)
            continue
        for d in range(h, -1, -1):
            stack.append((base + d * pw[j - 1], j - 1))
    return sorted(out), nodes


def a030979(path):
    terms = []
    with open(path) as f:
        for line in f:
            if line.startswith("A030979 "):
                return [int(t) for t in line.split(",")[1:] if t.strip()]
            if line[:2] in ("%S", "%T", "%U") and "A030979" in line:
                terms += [
                    int(t) for t in line.split(None, 2)[2].split(",") if t.strip()
                ]
    return terms


def inter(x, path=None):
    t0 = time.time()
    top = (x - 1) // 2
    small = 10**6
    for sides in ((3, 5), (3, 5, 7), (3, 5, 7, 11), (3, 5, 15)):
        brute = [k for k in range(small) if all(in_k(k, n) for n in sides)]
        got, _ = common(sides, small - 1)
        assert got == brute, sides
    print(
        f"control: pruned walk equals the direct digit test below k = {small} at four side sets"
    )
    oeis = a030979(path) if path else None
    if not oeis:
        print(
            f"A030979 comparison skipped: pass a copy of {A030979_URL} as the second argument"
        )
    for sides in (
        (3, 5),
        (3, 5, 15),
        (3, 5, 7),
        (3, 5, 7, 9),
        (3, 5, 7, 11),
        (3, 5, 7, 13),
        (3, 5, 7, 15),
        (3, 5, 7, 11, 13),
        (3, 5, 7, 11, 15),
    ):
        t = time.time()
        got, nodes = common(sides, top)
        ints = [2 * k for k in got]
        tail = f": {ints}" if len(ints) <= 20 else ""
        print(
            f"sides {sides}: {len(ints)} integers below {x:.0e}, {nodes} nodes, {time.time() - t:.2f} s{tail}"
        )
        if sides == (3, 5, 7) and oeis:
            print(
                f"  equals twice the A030979 terms below {x:.0e}: {got == [k for k in oeis if k <= top]}"
            )
    print(f"runtime {time.time() - t0:.2f} s")


def deep(e, sides):
    t0 = time.time()
    got, nodes = common(sides, (10**e - 1) // 2)
    print(
        f"sides {sides}: integers below 10^{e}: {[2 * k for k in got]}, {nodes} nodes, {time.time() - t0:.2f} s"
    )


# IRRATIONALS


def weyl(top, levels):
    t0 = time.time()
    scale = 10**90
    import mpmath

    mpmath.mp.dps = 110
    xs = {
        "sqrt(2) - 1": isqrt(2 * scale * scale) - scale,
        "(sqrt(5) - 1)/2": (isqrt(5 * scale * scale) - scale) // 2,
        "2^(1/3) - 1": int(mpmath.floor(mpmath.cbrt(2) * scale)) - scale,
        "pi - 3": int(mpmath.floor(mpmath.pi * scale)) - 3 * scale,
        "e - 2": int(mpmath.floor(mpmath.e * scale)) - 2 * scale,
    }
    sides = range(3, 2 * top + 2, 2)
    print(
        f"share of odd sides 3 <= N <= {2 * top + 1} holding x to level L (first L digits even), against 2^-L"
    )
    for name, v in xs.items():
        hits = [0] * (levels + 1)
        for n in sides:
            a = v
            for lev in range(1, levels + 1):
                a = a * n % (2 * scale)
                if a >= scale:
                    break
                hits[lev] += 1
        row = " ".join(
            f"{hits[lev] / len(sides) * 2**lev:.4f}" for lev in range(1, levels + 1)
        )
        print(f"  {name}: share times 2^L at L = 1..{levels}: {row}")
    print(f"runtime {time.time() - t0:.2f} s")


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else ""
    words = sys.argv[2:]
    args = [
        int(float(a)) for a in words if a.replace(".", "").replace("e", "").isdigit()
    ]
    if verb == "period":
        period(*(args or [200]))
    elif verb == "eisenstein":
        eisenstein(*(args or [200]))
    elif verb == "count":
        count(*(args or [10**6]))
    elif verb == "family":
        family()
    elif verb == "inter":
        files = [a for a in words if not a.replace(".", "").replace("e", "").isdigit()]
        inter(args[0] if args else 10**12, files[0] if files else None)
    elif verb == "second":
        second(*(args or [11]))
    elif verb == "deep":
        deep(args[0] if args else 30, (3, 5, 7, 11))
    elif verb == "weyl":
        weyl(*(args or [10**5, 6]))
    else:
        raise SystemExit(
            "verbs: period Q, eisenstein Q, count K, second E, family, inter X [FILE], deep D, weyl M L"
        )
