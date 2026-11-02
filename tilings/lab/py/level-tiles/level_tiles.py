import itertools
import math
import subprocess
import sys
import time
from collections import defaultdict
from functools import cache

# ARITHMETIC


@cache
def factor(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return tuple(sorted(out.items()))


def phi(n):
    r = n
    for p, _ in factor(n):
        r = r // p * (p - 1)
    return r


@cache
def divisors(n):
    ds = [1]
    for p, e in factor(n):
        ds = [d * p**i for d in ds for i in range(e + 1)]
    return tuple(sorted(ds))


def prime_power(n):
    f = factor(n)
    return f[0] if len(f) == 1 else None


def polydiv(num, den):
    num, q = list(num), [0] * max(len(num) - len(den) + 1, 1)
    for i in range(len(num) - len(den), -1, -1):
        c = num[i + len(den) - 1] // den[-1]
        q[i] = c
        if c:
            for j, d in enumerate(den):
                num[i + j] -= c * d
    rem = num[: len(den) - 1]
    while rem and rem[-1] == 0:
        rem.pop()
    return q, rem


@cache
def cyclo(m):
    poly = [-1] + [0] * (m - 1) + [1]
    for d in divisors(m)[:-1]:
        poly, _ = polydiv(poly, cyclo(d))
    return tuple(poly)


def mask(A):
    out = [0] * (max(A) + 1)
    for a in A:
        out[a] += 1
    return out


def valuation(poly, m):
    v, c = 0, cyclo(m)
    while True:
        q, r = polydiv(poly, c)
        if r:
            return v
        v, poly = v + 1, q


# ZERO SETS


def zero_set(F):
    F = [f - min(F) for f in F]
    deg, poly, out = max(F), mask(F), {}
    for m in range(2, 2 * deg * deg + 3):
        if phi(m) <= deg:
            red = [0] * m
            for f in F:
                red[f % m] += 1
            if not polydiv(red + [0], cyclo(m))[1]:
                out[m] = valuation(poly, m)
    return out


def level_zero_set(Z, b, n):
    out = defaultdict(int)
    for i in range(n):
        B = b**i
        for u, mu in Z.items():
            for g in divisors(B):
                if math.gcd(u, B // g) == 1:
                    out[u * g] += mu
    return dict(out)


def in_level(M, Z, b, n):
    return any(M // math.gcd(M, b**i) in Z for i in range(n))


def level_set(F, b, n):
    A = [0]
    for i in range(n):
        A = [a + f * b**i for a in A for f in F]
    return A


def exps(Z):
    out = defaultdict(set)
    for m in Z:
        pp = prime_power(m)
        if pp:
            out[pp[0]].add(pp[1])
    return out


def level_exps(Z, b, n):
    out = defaultdict(set)
    fb = dict(factor(b))
    for p, cs in exps(Z).items():
        e = fb.get(p, 0)
        for c in cs:
            for i in range(n if e else 1):
                out[p].add(c + e * i)
    return out


# COVEN-MEYEROWITZ


def t1(k, E):
    return k == math.prod(p ** len(cs) for p, cs in E.items())


def t2_fail(E, member):
    ps = sorted(E)
    for r in range(2, len(ps) + 1):
        for sub in itertools.combinations(ps, r):
            for choice in itertools.product(*[sorted(E[p]) for p in sub]):
                M = math.prod(p**a for p, a in zip(sub, choice))
                if not member(M):
                    return M
    return None


def omega(n):
    return len(factor(n))


def verdict(F, b, n, Z=None, size=None):
    Z = zero_set(F) if Z is None else Z
    k = (size or len(F)) ** n
    E = level_exps(Z, b, n)
    if not t1(k, E):
        return "no", "T1", E
    bad = t2_fail(E, lambda M: in_level(M, Z, b, n))
    if bad is None:
        return "tile", "", E
    return ("no" if omega(k) <= 2 else "open"), "T2 at %d" % bad, E


def depth(F, b, nmax, Z=None, size=None):
    Z = zero_set(F) if Z is None else Z
    for n in range(1, nmax + 1):
        v, why, _ = verdict(F, b, n, Z, size)
        if v != "tile":
            return n - 1, v, why
    return nmax, "tile", ""


def condition_p(F, b, Z):
    fb = dict(factor(b))
    E = exps(Z)
    if not t1(len(F), E):
        return False
    for p, cs in E.items():
        e = fb.get(p, 0)
        if e == 0 or len({c % e for c in cs}) < len(cs):
            return False
    return True


def horizon(Z, b):
    fb = dict(factor(b))
    E = exps(Z)
    C = defaultdict(int)
    for u in Z:
        for p, a in factor(u):
            C[p] = max(C[p], a)
    G = max((-(-C[p] // fb[p]) for p in E), default=0)
    return 1 + G * (2 * len(E) - 1)


def every_level(Z, size, b):
    E = exps(Z)
    fb = dict(factor(b))
    if not t1(size, E):
        return False, 0
    for p, cs in E.items():
        e = fb.get(p, 0)
        if e == 0 or len({c % e for c in cs}) < len(cs):
            return False, 0
    n0 = horizon(Z, b)
    for n in range(1, n0 + 1):
        if t2_fail(level_exps(Z, b, n), lambda M: in_level(M, Z, b, n)) is not None:
            return False, n0
    return True, n0


def t1_depth(F, b, Z):
    fb = dict(factor(b))
    E = exps(Z)
    if not t1(len(F), E):
        return 0
    best = math.inf
    for p, cs in E.items():
        e = fb.get(p, 0)
        if e == 0:
            return 1
        for c, d in itertools.combinations(sorted(cs), 2):
            if (d - c) % e == 0:
                best = min(best, (d - c) // e)
    return best


def tiles_residues(F, b):
    if b % len(F):
        return False
    return complement(F, b) is not None


def normalize(F):
    m = min(F)
    g = math.gcd(*[f - m for f in F])
    return tuple(sorted((f - m) // g for f in F))


# CERTIFICATES


def cm_complement(A, E):
    L = math.prod(p ** max(cs) for p, cs in E.items())
    B = [1]
    for p, cs in E.items():
        t = L // p ** max(cs)
        for a in range(1, max(cs) + 1):
            if a not in cs:
                c = cyclo(p**a)
                f = [0] * ((len(c) - 1) * t + 1)
                for j, x in enumerate(c):
                    f[j * t] = x
                B = polymul(B, f)
    return L, B


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def certify(A, E):
    L, B = cm_complement(A, E)
    if any(c not in (0, 1) for c in B):
        return L, None
    Bs = [i for i, c in enumerate(B) if c]
    seen = {(a + x) % L for a in A for x in Bs}
    return L, (Bs if len(seen) == L == len(A) * len(Bs) else None)


def complement(A, N, cap=2_000_000):
    R = sorted({a % N for a in A})
    if len(R) < len(A) or N % len(A):
        return None
    full = (1 << N) - 1
    base = sum(1 << r for r in R)
    rot = lambda s: ((base << s) | (base >> (N - s))) & full
    steps = [0]

    def rec(cov, B):
        if cov == full:
            return B
        steps[0] += 1
        if steps[0] > cap:
            raise TimeoutError
        r = ((cov + 1) & ~cov).bit_length() - 1
        for a in R:
            s = (r - a) % N
            m = rot(s)
            if not cov & m:
                got = rec(cov | m, B + [s])
                if got is not None:
                    return got
        return None

    return rec(0, [])


# SPECTRA


def spectrum(A, Zl, cap=3_000_000):
    L = math.lcm(*Zl)
    D = [d for d in range(1, L) if L // math.gcd(d, L) in Zl]
    idx = {d: i for i, d in enumerate(D)}
    nb = []
    for d in D:
        m = 0
        for e in D:
            if e != d and ((e - d) % L) in idx:
                m |= 1 << idx[e]
        nb.append(m)
    target = len(A) - 1
    steps = [0]

    def bk(size, cand, chosen):
        if size == target:
            return chosen
        steps[0] += 1
        if steps[0] > cap:
            raise TimeoutError
        while cand:
            if size + cand.bit_count() < target:
                return None
            v = cand.bit_length() - 1
            got = bk(size + 1, cand & nb[v], chosen + [D[v]])
            if got is not None:
                return got
            cand &= ~(1 << v)
        return None

    if target == 0:
        return L, [0]
    for m in sorted(Zl):
        v = idx[L // m]
        got = bk(1, nb[v], [L // m])
        if got is not None:
            return L, [0] + got
    return L, None


def check_spectrum(A, L, S):
    for x, y in itertools.combinations(S, 2):
        z = sum(
            complex(
                math.cos(2 * math.pi * a * (x - y) / L),
                math.sin(2 * math.pi * a * (x - y) / L),
            )
            for a in A
        )
        if abs(z) > 1e-6:
            return False
    return True


# PARI


def gp(lines):
    out = subprocess.run(
        ["gp", "-q", "-D", "parisizemax=2G"],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    return out.strip().splitlines()


def pari_cyclo(polys):
    cmd = []
    for P in polys:
        s = "+".join("x^%d" % a for a in P)
        cmd.append(
            "{f=factor(%s);print(vector(#f~,j,[poliscyclo(f[j,1]),f[j,2]]))}" % s
        )
    out = []
    for line in gp(cmd):
        pairs = eval(line.replace(";", ","))
        out.append({m: e for m, e in pairs if m})
    return out


def pari_unimodular(polys):
    cmd = ["default(realprecision,60);"]
    for P in polys:
        s = "+".join("x^%d" % a for a in P)
        cmd.append(
            "{f=factor(%s);g=1;for(j=1,#f~,if(!poliscyclo(f[j,1]),g*=f[j,1]^f[j,2]));h=gcd(g,polrecip(g));c=0;if(poldegree(h)>0,r=polroots(h);for(j=1,#r,if(abs(abs(r[j])-1)<1e-30,c++)));print(c)}"
            % s
        )
    return [int(x) for x in gp(cmd)]


# VERBS


def digit_sets(b, sizes=None):
    for r in range(1, b):
        if sizes and r + 1 not in sizes:
            continue
        for rest in itertools.combinations(range(1, b), r):
            yield (0,) + rest


def mirror_rep(F):
    G = tuple(sorted(max(F) - f for f in F))
    return min(F, G)


def ppfactor(F, p, a):
    q, cnt = p ** (a - 1), defaultdict(int)
    for f in F:
        cnt[f % (p * q)] += 1
    return all(cnt[r] == cnt[r + j * q] for r in range(q) for j in range(1, p))


def hand():
    F, b = (0, 2), 3
    A = level_set(F, b, 2)
    Z = zero_set(F)
    v, why, E = verdict(F, b, 2, Z)
    print("level 2 of {0,2} at base 3:", sorted(A))
    print("PARI cyclotomic factors of the mask:", pari_cyclo([A])[0])
    print("index lemma:", level_zero_set(Z, b, 2))
    print(
        "T1: |A_2| = %d against prod over S of Phi_s(1) = %d, S = %s"
        % (
            len(A),
            math.prod(p ** len(c) for p, c in E.items()),
            {p: sorted(c) for p, c in E.items()},
        )
    )
    print("verdict:", v, why)
    found = [N for N in range(4, 129, 4) if complement(A, N) is not None]
    print("complements in Z/N for 4 | N <= 128:", found or "none")


def index():
    polys, preds, t0 = [], [], time.time()
    for b in range(2, 9):
        for F in digit_sets(b):
            Z = zero_set(F)
            for n in range(1, (4 if b <= 4 else 3) + 1):
                A = level_set(F, b, n)
                polys.append(A)
                preds.append(level_zero_set(Z, b, n))
    got = pari_cyclo(polys)
    bad = sum(1 for g, p in zip(got, preds) if g != p)
    print(
        "masks factored by PARI: %d (bases 2..8, every digit set with 0, levels 1..3, level 4 at bases <= 4)"
        % len(polys)
    )
    print(
        "cyclotomic factorisation equal to the index lemma, with multiplicity: %d of %d"
        % (len(polys) - bad, len(polys))
    )
    print("time %.1f s" % (time.time() - t0))


def census(bmax=16, nmax=8):
    t0 = time.time()
    print(
        "base   sets  tileZb  every  T1every  raw  norm  hor  depth histogram (inf = every level)"
    )
    cxs, tot = [], defaultdict(int)
    for b in range(2, bmax + 1):
        row, hist = defaultdict(int), defaultdict(int)
        for F in digit_sets(b):
            Z = zero_set(F)
            d, _, _ = depth(F, b, nmax, Z)
            every, n0 = every_level(Z, len(F), b)
            P = condition_p(F, b, Z)
            T1d = t1_depth(F, b, Z)
            tz = tiles_residues(F, b)
            row["sets"] += 1
            row["tz"] += tz
            row["every"] += every
            row["P"] += P
            hist["inf" if every else (d if d < nmax else ">=%d" % nmax)] += 1
            assert P == (T1d == math.inf)
            T1read = next(
                (
                    n - 1
                    for n in range(1, 17)
                    if not t1(len(F) ** n, level_exps(Z, b, n))
                ),
                16,
            )
            assert T1read == min(T1d, 16)
            T2read = next(
                (
                    n - 1
                    for n in range(1, 17)
                    if t2_fail(level_exps(Z, b, n), lambda M: in_level(M, Z, b, n))
                    is not None
                ),
                16,
            )
            assert d <= T1d
            if omega(len(F)) <= 2:
                assert d == min(T1read, T2read, nmax)
                row["t1law"] += 1
                row["t2cut"] += d < min(T1d, nmax)
            assert (not every) or (P and d == nmax)
            assert (not tz) or every
            if P:
                far = all(
                    t2_fail(level_exps(Z, b, n), lambda M: in_level(M, Z, b, n)) is None
                    for n in range(1, n0 + 7)
                )
                assert far == every
                row["horizon"] = max(row["horizon"], n0)
            if prime_power(len(F)):
                assert d == min(T1d, nmax) and every == P
            if prime_power(b):
                assert every == tz
            if every:
                assert b % len(F) == 0
                if not prime_power(len(F)) and b % math.prod(
                    p ** max(c) for p, c in exps(Z).items()
                ):
                    row["wide"] += 1
            if every and not tz:
                row["raw"] += 1
                if not tiles_residues(normalize(F), b):
                    row["norm"] += 1
                    cxs.append((b, F))
        for key in row:
            tot[key] += row[key]
        h = " ".join(
            "%s:%d" % (k, hist[k])
            for k in sorted(hist, key=lambda x: (isinstance(x, str), str(x).rjust(4)))
        )
        print(
            "%4d %6d %7d %6d %8d %4d %5d %4d  %s"
            % (
                b,
                row["sets"],
                row["tz"],
                row["every"],
                row["P"],
                row["raw"],
                row["norm"],
                row["horizon"],
                h,
            )
        )
    print(
        "bases 2..%d, every digit set with 0 and at least two digits: %d sets; tile Z/b %d; T1 and T2 at every level %d; T1 at every level %d"
        % (bmax, tot["sets"], tot["tz"], tot["every"], tot["P"])
    )
    print(
        "T1 depth formula equal to the T1 depth read level by level to level 16 on all %d sets; tiling depth equal to the least of T1 depth, T2 depth and %d on the %d sets of size with at most two primes, T2 cutting below the T1 depth on %d"
        % (tot["sets"], nmax, tot["t1law"], tot["t2cut"])
    )
    print(
        "every level tiles but F does not tile Z/b: %d; and F/gcd does not either: %d; every-level sets with a composite non-prime-power size and lcm(S_F) not dividing b: %d"
        % (tot["raw"], tot["norm"], tot["wide"])
    )
    reps = sorted({(b, mirror_rep(F)) for b, F in cxs})
    print(
        "normalized counterexamples up to mirror (base, digits, zero set of hat F at roots of unity):"
    )
    for b, F in reps:
        print("  %d %s %s" % (b, F, sorted(zero_set(F))))
    print("time %.1f s" % (time.time() - t0))


def search(bmax=8, nmax=2, Nmax=512, cbmax=12, cnmax=3):
    t0 = time.time()
    agree, total, timeouts = 0, 0, 0
    for b in range(2, bmax + 1):
        for F in digit_sets(b):
            Z = zero_set(F)
            for n in range(1, nmax + 1):
                A = level_set(F, b, n)
                v, why, E = verdict(F, b, n, Z)
                found = None
                try:
                    for N in range(len(A), Nmax + 1, len(A)):
                        if complement(A, N) is not None:
                            found = N
                            break
                except TimeoutError:
                    timeouts += 1
                    continue
                total += 1
                ok = (found is not None) == (v == "tile")
                agree += ok
                if not ok:
                    print("  disagreement", b, F, n, v, why, found)
    print(
        "levels at bases <= %d, levels <= %d, checked by exhaustive complement search in Z/N for |A| | N <= %d: %d; agreeing with the T1 T2 verdict: %d; search capped: %d"
        % (bmax, nmax, Nmax, total, agree, timeouts)
    )
    certs, tiles = 0, 0
    for b in range(2, cbmax + 1):
        for F in digit_sets(b):
            Z = zero_set(F)
            for n in range(1, cnmax + 1):
                v, why, E = verdict(F, b, n, Z)
                if v == "tile":
                    tiles += 1
                    _, B = certify(level_set(F, b, n), E)
                    certs += B is not None
    print(
        "tile verdicts at bases <= %d, levels <= %d: %d; certified by the Coven-Meyerowitz complement with A + B checked to be all residues mod lcm(S): %d"
        % (cbmax, cnmax, tiles, certs)
    )
    print("time %.1f s" % (time.time() - t0))


def witness():
    for F, b, n in (
        ((0, 2), 6, 3),
        ((0, 1, 8, 9), 12, 3),
        ((0, 1, 2, 6, 7, 8), 18, 2),
        ((0, 1, 4, 5), 6, 3),
        ((0, 1, 8, 9), 10, 3),
    ):
        Z = zero_set(F)
        print(
            "digits %s base %d: zero set of hat F at roots of unity %s, normalized %s tiles Z/%d: %s"
            % (F, b, sorted(Z), normalize(F), b, tiles_residues(normalize(F), b))
        )
        Ls = [level_set(F, b, m) for m in range(1, n + 1)]
        facs = pari_cyclo(Ls)
        for m, (A, fac) in enumerate(zip(Ls, facs), 1):
            v, why, E = verdict(F, b, m, Z)
            line = (
                "  level %d: |A| = %d, PARI cyclotomic part %s, S = %s, verdict %s %s"
                % (
                    m,
                    len(A),
                    dict(sorted(fac.items())),
                    {p: sorted(c) for p, c in sorted(E.items())},
                    v,
                    why,
                )
            )
            if v == "tile":
                L, B = certify(A, E)
                line += ", complement %s mod %d" % (
                    B if len(B) <= 8 else "of size %d" % len(B),
                    L,
                )
            print(line)
    for F in ((0, 1, 8, 9), (0, 3, 8, 11)):
        hits = [
            v for v in range(1, 13) if complement([v * f for f in F], 12) is not None
        ]
        z2 = sorted(
            {
                a
                for v in range(1, 13)
                for a in exps(zero_set(tuple(v * f for f in F)))[2]
            }
        )
        print(
            "multiples v %s with v = 1..12 (all classes mod 12) tiling Z/12: %s; exponents of 2 in S_vF over these v: %s"
            % (F, hits or "none", z2)
        )
    for F, b in (((0, 1, 4, 5), 6), ((0, 1, 8, 9), 10), ((0, 1, 256, 257), 18)):
        Z = zero_set(F)
        print(
            "T1 depth of %s at base %d: %s; tiling depth read to level 10: %d"
            % (F, b, t1_depth(F, b, Z), depth(F, b, 10, Z)[0])
        )
    print(
        "digit sets that satisfy T1 at every level, by base and size with two primes:"
    )
    t2search()


def t2search(bmax=30, cap=200_000):
    t0 = time.time()
    for b in range(2, bmax + 1):
        for k in divisors(b):
            if omega(k) < 2 or k == b:
                continue
            if math.comb(b - 1, k - 1) > cap:
                print(
                    "base %d, |F| = %d: %d digit sets, skipped"
                    % (b, k, math.comb(b - 1, k - 1))
                )
                continue
            fb = dict(factor(b))
            cnt, passp, own, later = 0, 0, 0, []
            for F in digit_sets(b, {k}):
                cnt += 1
                E = defaultdict(set)
                for p in fb:
                    a = 1
                    while phi(p**a) <= b - 1:
                        if ppfactor(F, p, a):
                            E[p].add(a)
                        a += 1
                if not t1(k, E) or any(
                    len({c % fb[p] for c in cs}) < len(cs) for p, cs in E.items()
                ):
                    continue
                Z = zero_set(F)
                if not condition_p(F, b, Z):
                    continue
                passp += 1
                d, _, why = depth(F, b, 6, Z)
                every, _ = every_level(Z, k, b)
                assert every == (d == 6)
                if d == 0:
                    own += 1
                elif d < 6:
                    later.append((d, mirror_rep(F), why))
            reps = sorted(set(later))
            print(
                "base %d, |F| = %d: %d digit sets, %d with T1 at every level, %d fail T2 already at level 1, %d tile and fail T2 at a level <= 6 (%d up to mirror)%s"
                % (
                    b,
                    k,
                    cnt,
                    passp,
                    own,
                    len(later),
                    len(reps),
                    (
                        ": "
                        + ", ".join(
                            "%s at level %d (%s)" % (F, d + 1, w)
                            for d, F, w in reps[:6]
                        )
                    )
                    if reps
                    else "",
                )
            )
    print("time %.1f s" % (time.time() - t0))


def spectral(bmax=12, nmax=3, amax=36):
    t0 = time.time()
    sets = [(b, F) for b in range(2, bmax + 1) for F in digit_sets(b)]
    uni = pari_unimodular([list(F) for b, F in sets])
    print(
        "digit sets with a root on the unit circle that is not a root of unity: %d of %d"
        % (sum(1 for u in uni if u), len(sets))
    )
    tally, caps, odd = defaultdict(int), 0, []
    for (b, F), u in zip(sets, uni):
        Z = zero_set(F)
        for n in range(1, nmax + 1):
            if len(F) ** n > amax:
                continue
            A = level_set(F, b, n)
            v, _, _ = verdict(F, b, n, Z)
            Zl = level_zero_set(Z, b, n)
            S = None
            if Zl:
                try:
                    L, S = spectrum(A, set(Zl))
                except TimeoutError:
                    caps += 1
                    continue
            if S is not None:
                assert check_spectrum(A, L, S)
                w = "spectral"
            else:
                w = "undecided" if u else "not spectral"
            tally[(v, w)] += 1
            if (v == "tile") != (w == "spectral") and w != "undecided":
                odd.append((b, F, n, v, w))
    for key in sorted(tally):
        print("  tiling verdict %-5s  %-13s %d" % (key[0], key[1], tally[key]))
    for x in odd[:10]:
        print("  split:", x)
    print(
        "levels with |A| <= %d at bases <= %d, levels <= %d: %d; tile and spectral disagree on %d; search capped: %d; time %.1f s"
        % (amax, bmax, nmax, sum(tally.values()), len(odd), caps, time.time() - t0)
    )


def lwz(N, m, L, p):
    if p % N or p % L:
        return False
    d = max(i for i in range(200) if math.gcd(m * L // math.gcd(m * L, p**i), L) != 1)
    return (m // math.gcd(m, p**d)) % N == 0


def product_zero_set(N, m, L):
    Z = defaultdict(int)
    for u in divisors(N)[1:]:
        Z[u] += 1
    for v in divisors(L)[1:]:
        for g in divisors(m):
            if math.gcd(v, m // g) == 1:
                Z[v * g] += 1
    return dict(Z)


def measure(pmax=24, cmax=12):
    t0 = time.time()
    agree, total, opens, bad = 0, 0, defaultdict(int), []
    for p in range(2, pmax + 1):
        for N in range(2, cmax + 1):
            for L in range(2, cmax + 1):
                for m in range(N, p * p + 1):
                    Z = product_zero_set(N, m, L)
                    ours, _ = every_level(Z, N * L, p)
                    if (
                        not ours
                        and omega(N * L) > 2
                        and condition_p(range(N * L), p, Z)
                    ):
                        opens[lwz(N, m, L, p)] += 1
                        continue
                    theirs = lwz(N, m, L, p)
                    total += 1
                    agree += ours == theirs
                    if ours != theirs and len(bad) < 10:
                        bad.append((N, m, L, p, ours, theirs))
    print(
        "product-form digit sets D_N + m D_L at base p <= %d, 2 <= N, L <= %d, N <= m <= p^2: %d"
        % (pmax, cmax, total)
    )
    print(
        "T1 and T2 at every level  ==  the spectral condition of Liu-Wang-Zheng Theorem 1.3: %d of %d"
        % (agree, total)
    )
    print(
        "skipped, T1 at every level and T2 failing with three primes in N L: %d, of which the spectral condition holds on %d"
        % (sum(opens.values()), opens[True])
    )
    for x in bad:
        print("  disagreement N=%d m=%d L=%d p=%d ours=%s lwz=%s" % x)
    cons = [(N, p) for p in range(2, 65) for N in range(2, p + 1)]
    ok = sum(
        1
        for N, p in cons
        if every_level(zero_set(tuple(range(N))), N, p)[0] == (p % N == 0)
    )
    print(
        "consecutive digits {0..N-1} at base p <= 64: T1 and T2 at every level iff N | p on %d of %d"
        % (ok, len(cons))
    )
    print("time %.1f s" % (time.time() - t0))


VERBS = {
    "hand": hand,
    "index": index,
    "census": census,
    "search": search,
    "witness": witness,
    "spectral": spectral,
    "measure": measure,
}

if __name__ == "__main__":
    for v in sys.argv[1:] or list(VERBS):
        print("== %s" % v)
        VERBS[v]()
