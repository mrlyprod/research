import argparse
import time
from math import ceil, floor, gcd, log

import numpy as np


def up(x, d):
    return ceil(x * 10**d) / 10**d


def down(x, d):
    return floor(x * 10**d) / 10**d


# INCREMENT AUTOMATON


def level(z1, z2):
    start = z1 // 3
    if start == 0:
        return 0, 0
    dist = {start: 1}
    frontier = [start]
    d = 1
    while frontier:
        d += 1
        nxt = []
        for j in frontier:
            if j % 3 == 0:
                moves = (j // 3, (j + z1) // 3)
            elif (j - z2) % 3 == 0:
                moves = ((j - z2) // 3,)
            else:
                continue
            for k in moves:
                if k == 0:
                    return d, len(dist)
                if k not in dist:
                    dist[k] = d
                    nxt.append(k)
        frontier = nxt
    return 0, len(dist)


def reach_size(z1, z2):
    start = z1 // 3
    if start == 0:
        return 0
    seen = {0, start}
    stack = [start]
    while stack:
        j = stack.pop()
        if j % 3 == 0:
            moves = (j // 3, (j + z1) // 3)
        elif (j - z2) % 3 == 0:
            moves = ((j - z2) // 3,)
        else:
            continue
        for k in moves:
            if k not in seen:
                seen.add(k)
                stack.append(k)
    return len(seen)


def band_cap(z1, z2):
    return (z1 - 1) // 2 + (z2 - 1) // 2 + 1


def succ(j, z1, z2):
    if j % 3 == 0:
        return (j // 3, (j + z1) // 3)
    if (j - z2) % 3 == 0:
        return ((j - z2) // 3,)
    return ()


def reach(z1, z2):
    start = z1 // 3
    seen = {start}
    stack = [start]
    hit = 0
    while stack:
        for k in succ(stack.pop(), z1, z2):
            if k == 0:
                hit = 1
            if k not in seen:
                seen.add(k)
                stack.append(k)
    return hit, len(seen)


def core(z1, z2):
    lo = -((z2 - 1) // 2)
    hi = (z1 - 1) // 2
    deg = {}
    for j in range(lo, hi + 1):
        deg[j] = sum(1 for k in succ(j, z1, z2) if lo <= k <= hi)
    dead = [j for j in deg if j and not deg[j]]
    while dead:
        j = dead.pop()
        del deg[j]
        for p in (3 * j, 3 * j - z1, 3 * j + z2):
            if lo <= p <= hi and p in deg and p != 0 and j in succ(p, z1, z2):
                deg[p] -= 1
                if not deg[p]:
                    dead.append(p)
    live = set(deg)
    start = z1 // 3
    seen = {start}
    stack = [start]
    ok = start in live
    while stack:
        for k in succ(stack.pop(), z1, z2):
            if k in live:
                ok = True
            if k not in seen:
                seen.add(k)
                stack.append(k)
    return ok, len(live)


def base3(n):
    s = ""
    while n:
        s = str(n % 3) + s
        n //= 3
    return s or "0"


BETA = log(2) / log(3)


def binaries(k):
    out = [0]
    for i in range(k):
        out = out + [a + 3**i for a in out]
    return out


def pairs(cap):
    for z1 in range(3, cap + 1, 3):
        for z2 in range(1, cap + 1):
            if z2 % 3 and gcd(z1, z2) == 1:
                yield z1, z2


# PAIR CARRY AUTOMATON

OUT = ((0, 0), (0, 1), (1, 0))


def carry_occupied(z1, z2):
    seen = set()
    stack = []
    for d in (1, 2):
        s1 = d * z1
        s2 = d * z2
        if (s1 % 3, s2 % 3) in OUT:
            key = (s1 // 3) * z2 + (s2 // 3)
            if key not in seen:
                seen.add(key)
                stack.append((s1 // 3, s2 // 3))
    while stack:
        c1, c2 = stack.pop()
        e1 = c1 % 3
        if e1 == 2:
            continue
        for e2 in (0, 1) if e1 == 0 else (0,):
            d = (e2 - c2) * pow(z2 % 3, -1, 3) % 3
            n1 = (d * z1 + c1) // 3
            n2 = (d * z2 + c2) // 3
            if n1 == 0 and n2 == 0:
                return True, len(seen)
            key = n1 * z2 + n2
            if key not in seen:
                seen.add(key)
                stack.append((n1, n2))
    return False, len(seen)


# BRUTE GASKET


def brute(n, cap):
    out = {}
    for lab in range(3**n):
        u = v = 0
        w = lab
        for i in range(n):
            d = w % 3
            w //= 3
            if d == 1:
                u += 3**i
            elif d == 2:
                v += 3**i
        if u and v:
            g = gcd(u, v)
            z = (u // g, v // g)
            if max(z) <= cap:
                k = 1
                while 3**k <= max(u, v):
                    k += 1
                if z not in out or k < out[z]:
                    out[z] = k
    return out


# SWEEP


def cmd_sweep(args):
    cap = args.cap
    t = time.time()
    hs = []
    ws = []
    ls = []
    for z1, z2 in pairs(cap):
        d, _ = level(z1, z2)
        if d:
            hs.append(max(z1, z2))
            ws.append(z1 + z2)
            ls.append(d)
    secs = time.time() - t
    h = np.array(hs, dtype=np.int64)
    w = np.array(ws, dtype=np.int64)
    lv = np.array(ls, dtype=np.int64)
    print(
        "convention height max(z1,z2), ordered directions, occupied at some level, no level cap"
    )
    print("x A(x) theta local")
    ladder = []
    x = 32
    while x <= cap:
        ladder.append(x)
        x *= 2
    loc = []
    prev = None
    for x in ladder:
        a = 2 * int((h <= x).sum())
        th = log(a) / log(x)
        if prev:
            loc.append(log(a / prev[1]) / log(x / prev[0]))
        print(x, a, round(th, 4), round(loc[-1], 4) if prev else "-")
        prev = (x, a)
    print(
        "theta band",
        down(min(log(2 * (h <= x).sum()) / log(x) for x in ladder[-4:]), 4),
        up(max(log(2 * (h <= x).sum()) / log(x) for x in ladder[-4:]), 4),
        "local band",
        down(min(loc[-3:]), 4),
        up(max(loc[-3:]), 4),
        "over x",
        ladder[-4],
        ladder[-1],
    )
    print("W Zsum Zmean Zmax argmax logZmax/logW")
    cnt = np.bincount(w, minlength=cap + 1)[: cap + 1]
    zexp = []
    for x in ladder:
        seg = cnt[: x + 1]
        top = 2 * int(seg.max())
        zexp.append(log(top) / log(x))
        print(
            x,
            int(2 * seg.sum()),
            round(2 * seg.sum() / x, 4),
            top,
            int(seg.argmax()),
            round(zexp[-1], 4),
        )
    print(
        "Zmax exponent band",
        down(min(zexp), 4),
        up(max(zexp), 4),
        "over W",
        ladder[0],
        ladder[-1],
        "needed below",
        0.8073,
    )
    print("x meanlev maxlev meanc maxc share_lev_below_1.8073_log3_x")
    for x in ladder:
        sel = h <= x
        c = lv[sel] / (np.log(h[sel]) / log(3))
        print(
            x,
            round(float(lv[sel].mean()), 3),
            int(lv[sel].max()),
            round(float(c.mean()), 3),
            round(float(c.max()), 3),
            round(float((lv[sel] <= 1.8073 * log(x) / log(3)).mean()), 4),
        )
    print("secs", round(secs, 1))


# LEVELS


def cmd_levels(args):
    cap = args.cap
    t = time.time()
    tab = np.zeros(200, dtype=np.int64)
    for z1, z2 in pairs(cap):
        d, _ = level(z1, z2)
        if d:
            tab[d] += 2
    cum = np.cumsum(tab)
    print("cap", cap, "A(inf)", int(cum[-1]))
    print("n A(n,cap) new")
    for n in range(args.lo, args.hi + 1):
        print(n, int(cum[n]), int(tab[n]))
    print("secs", round(time.time() - t, 1))


# MISSING DIGIT MULTIPLES


def dcount(n, q):
    v = np.zeros(q, dtype=np.int64)
    v[0] = 1
    idx = np.arange(q)
    a = (3 * idx) % q
    b = (3 * idx + 1) % q
    for _ in range(n):
        nv = np.zeros(q, dtype=np.int64)
        nv[a] = v
        nv[b] += v
        v = nv
    return int(v[0])


def cmd_multiples(args):
    n = args.n
    print("uniform test D_n(q) q / 2^n over q <= Q coprime to 3, n =", n)
    worst = (0.0, 0)
    for q in range(2, args.qmax + 1):
        if q % 3 == 0:
            continue
        r = dcount(n, q) * q / 2**n
        if r > worst[0]:
            worst = (r, q)
    print("Q", args.qmax, "max ratio", round(worst[0], 4), "at q", worst[1])
    print("h q=1+3^h D_2h(q) 2^h 2^(2h)/q ratio")
    for h in range(1, args.hmax + 1):
        q = 1 + 3**h
        d = dcount(2 * h, q)
        print(h, q, d, 2**h, round(4**h / q, 2), round(d * q / 4**h, 3))


# BAND


def cmd_band(args):
    t = time.time()
    big = (0, None)
    frac = (0, 1, None)
    bad = 0
    for z1, z2 in pairs(args.cap):
        r = reach_size(z1, z2)
        cap = band_cap(z1, z2)
        if r > cap:
            bad += 1
        if r > big[0]:
            big = (r, (z1, z2))
        if r * frac[1] > frac[0] * cap:
            frac = (r, cap, (z1, z2))
    print("cap", args.cap, "violations of the halved cap", bad)
    print(
        "largest reachable set",
        big[0],
        "at",
        big[1],
        "halved cap",
        band_cap(*big[1]),
        "fraction",
        up(big[0] / band_cap(*big[1]), 4),
    )
    print(
        "fullest reachable set",
        frac[0],
        "of",
        frac[1],
        "at",
        frac[2],
        "fraction",
        up(frac[0] / frac[1], 4),
    )
    print("secs", round(time.time() - t, 1))


# WEIGHT LAYERS


def cmd_layers(args):
    cap = args.cap
    t = time.time()
    Z = [0] * (cap + 1)
    bad = 0
    for z1 in range(3, cap + 1, 3):
        for z2 in range(1, cap + 1):
            if z2 % 3 == 0 or gcd(z1, z2) != 1:
                continue
            w = z1 + z2
            if w > cap:
                break
            if level(z1, z2)[0]:
                Z[w] += 2
                if w <= 3 * z1 <= 2 * w:
                    bad += 1
    secs = time.time() - t
    print("cap", cap, "occupied directions with z1/w inside [1/3,2/3]", bad)
    print("octave argmax_Z Zmax base3 argmax_exp exp base3 argmax_const const base3")
    x = 32
    binary = True
    while 2 * x - 1 <= cap:
        seg = range(x, 2 * x)
        mz = max(seg, key=lambda w: Z[w])
        me = max(seg, key=lambda w: log(Z[w]) / log(w) if Z[w] else 0.0)
        mc = max(seg, key=lambda w: Z[w] / w**BETA)
        print(
            x,
            mz,
            Z[mz],
            base3(mz),
            me,
            up(log(Z[me]) / log(me), 4),
            base3(me),
            mc,
            up(Z[mc] / mc**BETA, 4),
            base3(mc),
        )
        binary = binary and all(set(base3(a)) <= {"0", "1"} for a in (mz, me, mc))
        x *= 2
    print("every argmax of all three columns binary base 3", binary)
    e = max(range(4, cap + 1), key=lambda w: log(Z[w]) / log(w) if Z[w] else 0.0)
    c = max(range(4, cap + 1), key=lambda w: Z[w] / w**BETA)
    print(
        "max exp",
        up(log(Z[e]) / log(e), 4),
        "at",
        e,
        "max const",
        up(Z[c] / c**BETA, 4),
        "at",
        c,
        "needed exp below",
        0.8073,
        "secs",
        round(secs, 1),
    )


def cmd_weights(args):
    print("w Z(w) phi3 logZ/logw Z/w^(log2/log3) meanreach/sqrt(w) secs")
    for w in args.w:
        t = time.time()
        z = n = r = 0
        for z1 in range(3, w, 3):
            if gcd(z1, w) != 1:
                continue
            n += 1
            hit, s = reach(z1, w - z1)
            z += 2 * hit
            r += s
        print(
            w,
            z,
            2 * n,
            up(log(z) / log(w), 4),
            up(z / w**BETA, 4),
            up(r / n / w**0.5, 4),
            round(time.time() - t, 1),
        )


def cmd_core(args):
    print("w Z(w) Zinf(w) Zinf/Z coremean band secs")
    for w in args.w:
        t = time.time()
        z = h = n = c = 0
        for z1 in range(3, w, 3):
            if gcd(z1, w) != 1:
                continue
            n += 1
            z += 2 * reach(z1, w - z1)[0]
            ok, sz = core(z1, w - z1)
            h += 2 * ok
            c += sz
        print(
            w,
            z,
            h,
            up(h / z, 4),
            round(c / n, 1),
            (w - 1) // 2,
            round(time.time() - t, 1),
        )


# SLOPE COVER


def tri(k, base):
    a = np.zeros(1, dtype=np.int64)
    b = np.zeros(1, dtype=np.int64)
    for i in range(k):
        p = 3 ** (base + i)
        a = np.concatenate([a, a + p, a])
        b = np.concatenate([b, b, b + p])
    return a, b


def cover(n, e, low=13):
    m = n + e
    if 3 * 3 ** (n + m) // 2 + 3**n >= 2**63:
        raise SystemExit("int64 overflow: reduce --extra or n")
    low = min(low, m)
    al, bl = tri(low, 0)
    ah, bh = tri(m - low, low)
    N = 3**n
    M = 3**m
    hit = np.zeros(N + 1, dtype=bool)
    u = np.empty(len(al), dtype=np.int64)
    v = np.empty(len(al), dtype=np.int64)
    c = np.empty(len(al), dtype=np.int64)
    d = np.empty(len(al), dtype=np.int64)
    for i in range(len(ah)):
        np.add(al, ah[i], out=u)
        np.add(bl, bh[i], out=v)
        np.add(u, M, out=c)
        np.add(c, v, out=d)
        c *= N
        c //= d
        hit[c] = True
        np.add(u, v, out=d)
        d += M
        np.multiply(u, N, out=c)
        c //= d
        hit[c] = True
    return int(hit[:N].sum())


def cmd_box(args):
    print(
        "cover counts both swap halves and is a lower estimate of the saturated cover"
    )
    print("saturation at n =", args.lo, "over extra digits 0 ..", args.extra)
    for e in range(args.extra + 1):
        print(args.lo, e, cover(args.lo, e))
    print("n cover(3^-n) 3^n cover*n/3^n log_3 cover / n local secs")
    prev = None
    for n in range(args.lo, args.hi + 1):
        t = time.time()
        c = cover(n, args.extra)
        print(
            n,
            c,
            3**n,
            down(c * n / 3**n, 4),
            down(log(c) / (n * log(3)), 4),
            "-" if prev is None else down(log(c / prev) / log(3), 4),
            round(time.time() - t, 1),
        )
        prev = c


# REPUNIT


def witness(z1, z2):
    start = z1 // 3
    par = {start: (None, z1)}
    frontier = [start]
    while frontier:
        nxt = []
        for j in frontier:
            for k in succ(j, z1, z2):
                a = k * 3 - j
                if k == 0:
                    path = [a]
                    cur = j
                    while cur is not None:
                        path.append(par[cur][1])
                        cur = par[cur][0]
                    e = sum(3**i for i, b in enumerate(reversed(path)) if b)
                    return e // (z1 + z2)
                if k not in par:
                    par[k] = (j, a)
                    nxt.append(k)
        frontier = nxt
    return 0


def factor(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def order3(q):
    d = 1
    x = 3 % q
    while x != 1:
        x = x * 3 % q
        d += 1
    return d


def count_multiples(k, q):
    v = np.zeros(q, dtype=np.int64)
    v[0] = 1
    idx = np.arange(q)
    for i in range(k):
        nv = v.copy()
        np.add.at(nv, (idx + pow(3, i, q)) % q, v)
        v = nv
    return int(v[0])


def fourier_multiples(k, q):
    d = order3(q)
    t = np.arange(q)
    P = np.ones(q, dtype=complex)
    for r in range(d):
        P *= 1 + np.exp(2j * np.pi * t * pow(3, r, q) / q)
    return (P ** (k // d)).sum() / q


def floor_formula(k):
    ps = factor((3**k - 1) // 2)
    exact = 0
    four = 0
    for mask in range(2 ** len(ps)):
        q = 1
        for i, p in enumerate(ps):
            if mask >> i & 1:
                q *= p
        sign = (-1) ** mask.bit_count()
        exact += sign * count_multiples(k, q)
        four += sign * (2**k if q == 1 else fourier_multiples(k, q))
    return ps, exact, four


def lift_union(k):
    w = (3**k - 1) // 2
    pw = np.array([3**i for i in range(k)], dtype=np.int64)
    S = np.arange(2**k)
    bits = ((S[:, None] >> np.arange(k)) & 1).astype(np.int64)
    occ = set()
    total = 0
    heur = 0.0
    for T in range(0, 2**k, 2):
        tb = np.array([(T >> i) & 1 for i in range(k)], dtype=np.int64)
        m = 1 + 2 * int(tb.dot(pw))
        heur += 2**k / m
        A = bits.dot(pw * (1 - tb) + pw * tb * 3**k)
        for a in A[A % m == 0]:
            z = int(a) // m
            if 0 < z < w and gcd(z, w) == 1:
                occ.add(z)
                total += 1
    return occ, total, heur


def cmd_repunit(args):
    print(
        "floor Phi_k = #{submask a of R_k : 0 < a < R_k, gcd(a, R_k) = 1}; Occ_T = directions with a witness K_T, "
        "T inside [1, k-1]; lift = |union of Occ_T|, liftsum = Sum |Occ_T|, heur = Sum 2^k / m_T; "
        "deep = Z - lift; delta = Phi / 2^k; X = (Z - Phi) / Phi; maxdigits = longest minimal witness m R_k, maxcol = most uses of one column mod k by a minimal witness"
    )
    print(
        "k R_k primes Phi_k formula fourier Z(R_k) excess lift liftsum heur deep Z/R_k^beta delta X maxdigits maxcol secs"
    )
    lists = {}
    for k in range(2, args.kmax + 1):
        t = time.time()
        w = (3**k - 1) // 2
        subs = {a for a in binaries(k) if 0 < a < w and gcd(a, w) == 1}
        ps, exact, four = floor_formula(k)
        occ = []
        for z1 in range(3, w, 3):
            if gcd(z1, w) != 1:
                continue
            m = witness(z1, w - z1)
            if m:
                occ.append((z1, m))
        Z = 2 * len(occ)
        if k <= args.liftmax:
            lu, liftsum, heur = lift_union(k)
            lift, deep, heur = len(lu), Z - len(lu), round(heur, 1)
        else:
            lift = deep = liftsum = heur = "-"
        digits = maxcol = 0
        for z1, m in occ:
            if z1 in subs:
                continue
            cols = [0] * k
            for i, ch in enumerate(reversed(base3(m * w))):
                cols[i % k] += ch == "1"
            digits = max(digits, len(base3(m * w)))
            maxcol = max(maxcol, max(cols))
        print(
            k,
            w,
            ps,
            len(subs),
            exact,
            round(four.real, 3) if abs(four.imag) < 1e-6 else four,
            Z,
            Z - len(subs),
            lift,
            liftsum,
            heur,
            deep,
            up(Z / w**BETA, 4),
            down(len(subs) / 2**k, 4),
            down((Z - len(subs)) / len(subs), 4),
            digits,
            maxcol,
            round(time.time() - t, 1),
        )
        if k in (7, 8, 9):
            lists[k] = [(z1, w - z1, m) for z1, m in occ if z1 not in subs]
    for k, rows in lists.items():
        print(
            "non-submask occupied directions of weight R_" + str(k),
            "=",
            (3**k - 1) // 2,
            "one of each swap pair, 3 | z1; columns z1 z2 m base3(z1) base3(z2) base3(m) base3(m z1) base3(m z2) base3(m w)",
        )
        for z1, z2, m in rows:
            print(
                z1,
                z2,
                m,
                base3(z1),
                base3(z2),
                base3(m),
                base3(m * z1),
                base3(m * z2),
                base3(m * (z1 + z2)),
            )


# LIFT UNION


def subset_sums(ws):
    out = np.zeros(1, dtype=np.int64)
    for x in ws:
        out = np.concatenate((out, out + x))
    return out


def submask_floor(k):
    w = (3**k - 1) // 2
    b = subset_sums([3**i for i in range(k)])
    return int(((b > 0) & (b < w) & (np.gcd(b, w) == 1)).sum())


def lift_scan(k):
    w = (3**k - 1) // 2
    h = max(k // 2, 1)
    v1 = {}
    for tl in range(0, 1 << h, 2):
        v1[tl] = subset_sums([3 ** (i + k) if tl >> i & 1 else 3**i for i in range(h)])
    v2 = {}
    for th in range(1 << (k - h)):
        v2[th] = subset_sums(
            [3 ** (i + k) if th >> (i - h) & 1 else 3**i for i in range(h, k)]
        )
    occ = set()
    total = 0
    model = 0.0
    best = (0.0, 0, 0)
    for th in range(1 << (k - h)):
        ah = sum(3**i for i in range(h, k) if th >> (i - h) & 1)
        V2 = v2[th]
        for tl in range(0, 1 << h, 2):
            al = sum(3**i for i in range(h) if tl >> i & 1)
            T = tl | (th << h)
            m = 1 + 2 * (al + ah)
            model += 2**k / m
            z = mod_match(v1[tl], V2, m, w)
            total += int(z.size)
            occ.update(z.tolist())
            ratio = z.size * m / 2**k
            if ratio > best[0]:
                best = (ratio, m, T)
    return occ, total, model, best


def occ_of(k, T):
    w = (3**k - 1) // 2
    m = 1 + 2 * sum(3**i for i in range(k) if T >> i & 1)
    A = subset_sums([3 ** (i + k) if T >> i & 1 else 3**i for i in range(k)])
    z = A[A % m == 0] // m
    z = z[(z > 0) & (z < w) & (np.gcd(z, w) == 1)]
    return m, set(z.tolist())


def cyclotomic_occ(t):
    k = 2 * t + 1
    w = (3**k - 1) // 2
    m = 3 ** (2 * t) - 3**t + 1
    out = set()
    for S in range(1, 1 << (t - 1)):
        c = sum(3 ** (1 + i) for i in range(t - 1) if S >> i & 1)
        z = c * (3**t + 1)
        if gcd(z, w) == 1:
            out.add(z)
            out.add(w - z)
    return m, out


def cmd_lifts(args):
    print(
        "Occ_T = {A / m_T : A submask of K_T = m_T R_k, m_T | A, 0 < A < K_T, gcd(A / m_T, R_k) = 1}, T inside [1, k-1], m_T = 1 + 2 a_T; "
        "U = |union of Occ_T|, sum = Sum_T |Occ_T|, model = Sum_T 2^k / m_T, peak = max_T |Occ_T| m_T / 2^k at multiplier mstar; "
        "Phi = submask floor, agg = sum 2^k / (Phi model) is the aggregate against the model after the coprime cut; "
        "cyc = |Occ_T| at T = [t, 2t-1] when k = 2t+1, pred its cyclotomic lower bound; deep = Z(R_k) - U where the band automaton reaches"
    )
    print("k R_k Phi U sum L U/2^k U/Phi U/sum agg peak mstar cyc pred Z deep secs")
    band = []
    for k in range(2, args.kmax + 1):
        t = time.time()
        w = (3**k - 1) // 2
        phi = submask_floor(k)
        occ, total, model, best = lift_scan(k)
        if k % 2 and k >= 5:
            j = k // 2
            _mstar, low = cyclotomic_occ(j)
            cyc = len(occ_of(k, ((1 << j) - 1) << j)[1])
            pred = len(low)
        else:
            cyc = pred = "-"
        if k <= args.zmax:
            z = 2 * sum(
                1 for z1 in range(3, w, 3) if gcd(z1, w) == 1 and witness(z1, w - z1)
            )
            deep = z - len(occ)
        else:
            z = deep = "-"
        agg = total * 2**k / (phi * model)
        band.append((len(occ) / phi, agg, model / 2**k, best[0], len(occ) / total))
        print(
            k,
            w,
            phi,
            len(occ),
            total,
            down(model / 2**k, 5),
            down(len(occ) / 2**k, 4),
            down(len(occ) / phi, 4),
            down(len(occ) / total, 4),
            down(agg, 4),
            down(best[0], 3),
            best[1],
            cyc,
            pred,
            z,
            deep,
            round(time.time() - t, 1),
        )
    tail = band[args.lo - 2 :]
    names = ("U/Phi", "agg", "L", "peak", "U/sum")
    print(
        "bands over k =",
        args.lo,
        "..",
        args.kmax,
        "".join(
            " "
            + n
            + " ["
            + str(down(min(r[i] for r in tail), 5))
            + ", "
            + str(up(max(r[i] for r in tail), 5))
            + "]"
            for i, n in enumerate(names)
        ),
    )


def raw_scan(k):
    h = max(k // 2, 1)
    v1 = {}
    for tl in range(0, 1 << h, 2):
        v1[tl] = subset_sums([3 ** (i + k) if tl >> i & 1 else 3**i for i in range(h)])
    v2 = {}
    for th in range(1 << (k - h)):
        v2[th] = subset_sums(
            [3 ** (i + k) if th >> (i - h) & 1 else 3**i for i in range(h, k)]
        )
    ms = []
    ns = []
    tops = []
    for th in range(1 << (k - h)):
        ah = sum(3**i for i in range(h, k) if th >> (i - h) & 1)
        V2 = v2[th]
        for tl in range(0, 1 << h, 2):
            al = sum(3**i for i in range(h) if tl >> i & 1)
            T = tl | (th << h)
            m = 1 + 2 * (al + ah)
            if m == 1:
                n = 1 << k
            else:
                r2s = np.sort(V2 % m)
                key = (m - v1[tl] % m) % m
                n = int(
                    (
                        np.searchsorted(r2s, key, side="right")
                        - np.searchsorted(r2s, key, side="left")
                    ).sum()
                )
            ms.append(m)
            ns.append(n)
            tops.append(T.bit_length() - 1)
    return np.array(ms, dtype=np.int64), np.array(ns, dtype=np.int64), np.array(tops)


def residue_dist(k, T):
    m = 1 + 2 * sum(3**i for i in range(k) if T >> i & 1)
    v = np.zeros(m, dtype=np.int64)
    v[0] = 1
    for i in range(k):
        p = i + k if T >> i & 1 else i
        v = v + np.roll(v, pow(3, p, m))
    return m, v


def fourier_direct(k, T, m):
    u = np.arange(m)
    F = np.ones(m, dtype=complex)
    for i in range(k):
        p = i + k if T >> i & 1 else i
        F *= 1 + np.exp(2j * np.pi * (u * pow(3, p, m) % m) / m)
    return F


def phi2p_T(p, t):
    k = (p - 1) * t + 1
    T = 0
    for i in range(1, p - 1, 2):
        for j in range(t):
            T |= 1 << (t * i + j)
    return k, T


def cyclotomic_poly(n, x):
    num = 1
    den = 1
    for d in range(1, n + 1):
        if n % d == 0:
            mu = mobius(n // d)
            if mu == 1:
                num *= x**d - 1
            elif mu == -1:
                den *= x**d - 1
    return num // den


def mobius(n):
    ps = factor(n)
    m = n
    for p in ps:
        if m % (p * p) == 0:
            return 0
    return (-1) ** len(ps)


def raw_count(k, T):
    m = 1 + 2 * sum(3**i for i in range(k) if T >> i & 1)
    A = subset_sums([3 ** (i + k) if T >> i & 1 else 3**i for i in range(k)])
    return m, int((A % m == 0).sum())


def no_one_digits(u, m, lo, hi):
    ok = np.ones(u.size, dtype=bool)
    for j in range(lo, hi + 1):
        ok &= (u * 3**j // m) % 3 != 1
    return ok


def cmd_agg(args):
    print(
        "N_T = #{A submask of K_T : m_T | A} with A = 0 and K_T kept, T over every subset of [1, k-1] with {} included; "
        "F_T(u) = Prod_{p in supp K_T} (1 + e(u 3^p / m_T)); identity N_T = (1/m_T) Sum_u F_T(u), the u = 0 term 2^k / m_T; "
        "M = Sum_T N_T, L = Sum_T 1/m_T, share = (M - 2^k L) / 2^k is the aggregate u != 0 part, "
        "agg' = Sum_T (N_T - 2) / (2^k L) the cut-free aggregate against the model, "
        "Abs = Sum_T (1/m_T) Sum_{u != 0} |F_T(u)| / 2^k, top = max_T max_{u != 0} F_T(u) / 2^k at T*, "
        "iderr = max over T and u of |F_T from the product - F_T from the FFT|, imag = max |Im F_T|"
    )
    print("k M M/2^k L share agg' Abs Abs/prev top T* iderr imag secs")
    prev = None
    for k in range(args.ulo, args.umax + 1):
        t0 = time.time()
        M = 0
        L = 0.0
        absum = 0.0
        top = (0.0, 0)
        iderr = 0.0
        imag = 0.0
        for T in range(0, 1 << k, 2):
            m, v = residue_dist(k, T)
            F = np.fft.fft(v)
            imag = max(imag, float(np.abs(F.imag).max()))
            Fr = F.real
            n = int(v[0])
            if abs(Fr.mean() - n) > 1e-6 * max(n, 1):
                raise SystemExit("identity fails at k %d T %d" % (k, T))
            if k <= args.directmax:
                D = fourier_direct(k, T, m)
                iderr = max(iderr, float(np.abs(D - np.conj(F)).max()))
            M += n
            L += 1 / m
            if m > 1:
                absum += (np.abs(Fr).sum() - 2**k) / m
                mx = float(Fr[1:].max()) / 2**k
                if mx > top[0]:
                    top = (mx, T)
        share = (M - 2**k * L) / 2**k
        aggp = (M - 2 * (1 << (k - 1))) / (2**k * L)
        absk = absum / 2**k
        Tset = [i for i in range(k) if top[1] >> i & 1]
        print(
            k,
            M,
            down(M / 2**k, 5),
            down(L, 5),
            down(share, 5),
            down(aggp, 5),
            down(absk, 4),
            "-" if prev is None else down(absk / prev, 4),
            down(top[0], 5),
            Tset,
            "%.1e" % iderr if k <= args.directmax else "-",
            "%.1e" % imag,
            round(time.time() - t0, 1),
        )
        prev = absk
    print(
        "mass at the cyclotomic T = [t, 2t-1], k = 2t+1, m = Phi_6(3^t): N_T against 2^t, F_T(1)/2^k against its bound 1 - 13 9^(-t) for t >= 2, "
        "share_no1 = the share of Sum_{u != 0} F_T(u) carried by the u whose digits d_2..d_t of u/m avoid 1, frac_no1 = how many such u per m"
    )
    print("t k m N_T 2^t F(1)/2^k bound share_no1 frac_no1")
    for t in range(2, args.cycmax + 1):
        k = 2 * t + 1
        T = ((1 << t) - 1) << t
        m, v = residue_dist(k, T)
        F = np.fft.fft(v).real
        u = np.arange(m)
        ok = no_one_digits(u, m, 2, t)
        ok[0] = False
        rest = F[1:].sum()
        print(
            t,
            k,
            m,
            int(v[0]),
            2**t,
            down(F[1] / 2**k, 6),
            down(1 - 13 * 9.0 ** (-t), 6),
            down(F[ok].sum() / rest, 4),
            down(ok.sum() / m, 5),
        )
    print(
        "the antipodal family: p odd, t >= 1, 0 <= s <= t, k = (p-1) t + s, T = Union_{i odd <= p-2} [ti, ti+t-1], "
        "m_T = (3^(pt) + 1) / (3^t + 1), N_T = 2^((p-1)(t-s)/2 + s) exactly, 2^t at s = t; the verb raises on any mismatch; "
        "printed per (p, t): the k range, m_T, whether m_T is that quotient, and N_T against pred at every s"
    )
    print("p t k_lo..k_hi m_T quotient N_T(s=0..t) pred(s=0..t) equal")
    for p in range(3, 20, 2):
        for t in range(1, 20):
            if (p - 1) * t > args.kmax:
                break
            T = phi2p_T(p, t)[1]
            m = 1 + 2 * sum(3**i for i in range(64) if T >> i & 1)
            ns = []
            preds = []
            ok = True
            for s in range(t + 1):
                k = (p - 1) * t + s
                if k > args.kmax:
                    break
                n = raw_count(k, T)[1]
                pred = 2 ** ((p - 1) * (t - s) // 2 + s)
                ns.append(n)
                preds.append(pred)
                ok &= n == pred
            print(
                p,
                t,
                "%d..%d" % ((p - 1) * t, (p - 1) * t + len(ns) - 1),
                m,
                m == (3 ** (p * t) + 1) // (3**t + 1),
                ns,
                preds,
                ok,
            )
            if not ok or m != (3 ** (p * t) + 1) // (3**t + 1):
                raise SystemExit("antipodal family fails at p %d t %d" % (p, t))
    print(
        "the cut-free aggregate to k = kmax by meet in the middle: M, M/2^k, L, share, agg', "
        "and the excess by top digit t = max T, X_t = Sum_{max T = t} (N_T - 2 - 2^k / m_T) / 2^k at its argmax"
    )
    print("k M M/2^k L share agg' t* X_t* secs")
    band = []
    for k in range(2, args.kmax + 1):
        t0 = time.time()
        ms, ns, tops = raw_scan(k)
        M = int(ns.sum())
        L = float((1.0 / ms).sum())
        share = (M - 2**k * L) / 2**k
        aggp = (M - 2 * (1 << (k - 1))) / (2**k * L)
        ex = (ns - 2 - 2.0**k / ms) / 2**k
        ex[ms == 1] = 0
        xt = np.array([ex[tops == t].sum() for t in range(k)])
        ts = int(xt.argmax())
        band.append((M / 2**k, share, aggp))
        print(
            k,
            M,
            down(M / 2**k, 5),
            down(L, 5),
            down(share, 5),
            down(aggp, 5),
            ts,
            down(xt[ts], 5),
            round(time.time() - t0, 1),
        )
    tail = band[args.lo - 2 :]
    names = ("M/2^k", "share", "agg'")
    print(
        "bands over k =",
        args.lo,
        "..",
        args.kmax,
        "".join(
            " "
            + n
            + " ["
            + str(down(min(r[i] for r in tail), 5))
            + ", "
            + str(up(max(r[i] for r in tail), 5))
            + "]"
            for i, n in enumerate(names)
        ),
    )


# DEEP TAIL


def col_sums(opts):
    out = np.zeros(1, dtype=np.int64)
    for o in opts:
        out = np.concatenate([out + v for v in o])
    return out


def block_half(k, cols, slots, key):
    opts = []
    g = 0
    for i, r in enumerate(cols):
        e = key // 3**i % 3
        if slots == 1:
            opts.append([0, 3 ** (r + k * e)])
        else:
            ps = [3 ** (r + k * j) for j in range(3) if j != e]
            opts.append([0, ps[0], ps[1], ps[0] + ps[1]])
        g += 3**r * (0 if e == 0 else 1 if e == 1 else 3**k + 1)
    return col_sums(opts), g


def mod_match(V1, V2, m, w):
    r2 = V2 % m
    order = np.argsort(r2, kind="stable")
    r2s = r2[order]
    key = (m - V1 % m) % m
    lo = np.searchsorted(r2s, key, side="left")
    hi = np.searchsorted(r2s, key, side="right")
    cnt = hi - lo
    nz = np.nonzero(cnt)[0]
    if nz.size == 0:
        return np.zeros(0, dtype=np.int64)
    reps = cnt[nz]
    left = np.repeat(nz, reps)
    base = np.repeat(np.cumsum(reps) - reps, reps)
    pos = np.repeat(lo[nz], reps) + np.arange(reps.sum()) - base
    z = (V1[left] + V2[order[pos]]) // m
    return z[(z > 0) & (z < w) & (np.gcd(z, w) == 1)]


def block3_scan(k, slots):
    if 3 * k * log(3) / log(2) > 62:
        raise ValueError("the depth-3 census is int64 bound to k <= 13")
    w = (3**k - 1) // 2
    X = 3**k
    h = max(k // 2, 1)
    low = [block_half(k, list(range(h)), slots, key) for key in range(3**h)]
    occ = set()
    for kh in range(3 ** (k - h)):
        V2, g2 = block_half(k, list(range(h, k)), slots, kh)
        for kl, (V1, g1) in enumerate(low):
            if (kl % 3 == 0) != (slots == 1):
                continue
            g = g1 + g2
            m = 1 + 2 * g if slots == 1 else X * (X + 1) - 2 * g
            occ.update(mod_match(V1, V2, m, w).tolist())
    return occ


def column_vectors(k, b):
    n = b * k
    S = np.arange(1, 1 << n)
    bits = ((S[:, None] >> np.arange(n)) & 1).astype(np.int64)
    val = bits.dot(np.array([3**i for i in range(n)], dtype=np.int64))
    keep = val % ((3**k - 1) // 2) == 0
    return bits[keep].reshape(-1, b, k).sum(axis=1)


def lift_lengths(k):
    w = (3**k - 1) // 2
    dep = {}
    for z1 in range(3, w, 3):
        if gcd(z1, w) != 1:
            continue
        m = witness(z1, w - z1)
        if m:
            d = len(base3(m * w))
            dep[z1] = d
            dep[w - z1] = d
    return dep


def cmd_tail(args):
    print(
        "d(z) is the base-3 length of m(z) R_k for the minimal witness m(z) of an occupied direction of weight R_k, "
        "and the depth b(z) = ceil(d(z) / k) is the number of k-blocks that lift fills; "
        "U = |Union_T Occ_T| is depth at most 2, V = #{b(z) <= 3} the depth-3 census, Z = Z(R_k), tail = Z - U the deep tail; "
        "cap = V - U is the tail captured at depth 3 and share = cap / tail; "
        "fam = |U union the depth-3 lift families|, new1 and new2 what the all-1 and the all-2 family add beyond U; "
        "the depth histogram of the tail is printed under each row"
    )
    print("k R_k U V Z tail cap share fam new1 new2 maxdepth secs")
    rows = []
    for k in range(2, args.kmax + 1):
        t = time.time()
        w = (3**k - 1) // 2
        occ, _, _, _ = lift_scan(k)
        dep = lift_lengths(k)
        deep = sorted(set(dep) - occ)
        hist = {}
        for z in deep:
            b = -(-dep[z] // k)
            hist[b] = hist.get(b, 0) + 1
        V = len(occ) + sum(1 for z in deep if dep[z] <= 3 * k)
        fam = new1 = new2 = "-"
        o1 = set()
        if k <= min(max(args.onemax, args.famax), 13):
            o1 = block3_scan(k, 1)
            new1 = len(o1 - occ)
        if k <= min(args.famax, 13):
            o2 = block3_scan(k, 2)
            new2 = len(o2 - occ)
            fam = len(occ | o1 | o2)
            if fam != V:
                raise ValueError(
                    "depth-3 family census "
                    + str(fam)
                    + " against the automaton "
                    + str(V)
                    + " at k = "
                    + str(k)
                )
        share = down((V - len(occ)) / len(deep), 4) if deep else "-"
        rows.append((k, len(deep), V - len(occ)))
        print(
            k,
            w,
            len(occ),
            V,
            len(dep),
            len(deep),
            V - len(occ),
            share,
            fam,
            new1,
            new2,
            max(hist) if hist else 0,
            round(time.time() - t, 1),
        )
        if hist:
            print("  depth histogram of the tail at k =", k, sorted(hist.items()))
    live = [r for r in rows if r[1]]
    print(
        "tail",
        [r[1] for r in live],
        "captured",
        [r[2] for r in live],
        "at k =",
        [r[0] for r in live],
    )
    print(
        "tail growth over two steps",
        [down(live[i + 2][1] / live[i][1], 4) for i in range(len(live) - 2)],
        "against 4 for 2^k, captured growth",
        [
            down(live[i + 2][2] / live[i][2], 4)
            for i in range(len(live) - 2)
            if live[i][2]
        ],
    )


# CHECKS


def cmd_check(args):
    cap = 60
    for n in (6, 9, 12):
        b = brute(n, cap)
        miss = [z for z in b if level(*(z if z[0] % 3 == 0 else z[::-1]))[0] == 0]
        levbad = [z for z in b if level(*(z if z[0] % 3 == 0 else z[::-1]))[0] != b[z]]
        print(
            "brute n",
            n,
            "rays",
            len(b),
            "not occupied by automaton",
            len(miss),
            "level mismatch",
            len(levbad),
        )
    b = brute(12, cap)
    extra = [
        (a, c)
        for a in range(1, cap + 1)
        for c in range(1, cap + 1)
        if gcd(a, c) == 1
        and (a % 3 == 0) != (c % 3 == 0)
        and level(*((a, c) if a % 3 == 0 else (c, a)))[0]
        and (a, c) not in b
    ]
    print("automaton occupied but absent from brute n=12", len(extra))
    bad = sym = gap = core_bad = 0
    for z1, z2 in pairs(120):
        d, _size = level(z1, z2)
        o, _ = carry_occupied(z1, z2)
        if (d > 0) != o:
            bad += 1
        if reach_size(z1, z2) > band_cap(z1, z2):
            sym += 1
        if d and z1 + z2 <= 3 * z1 <= 2 * (z1 + z2):
            gap += 1
        if d and not core(z1, z2)[0]:
            core_bad += 1
    print(
        "j vs pair-carry disagreements to 120",
        bad,
        "band violations",
        sym,
        "occupied with z1/w inside [1/3,2/3]",
        gap,
        "occupied not reaching the core",
        core_bad,
    )
    print(
        "gasket rays with u/(u+v) inside [1/3,2/3] at n=12",
        sum(1 for a, c in brute(12, 3**12) if a + c <= 3 * a <= 2 * (a + c)),
    )
    b = brute(12, 60)
    print(
        "gasket rays symmetric under swap",
        all((c, a) in b for a, c in b),
        "every ray has exactly one coordinate divisible by 3",
        all((a % 3 == 0) != (c % 3 == 0) for a, c in b),
        "so no direction with 3 dividing neither coordinate is occupied",
    )
    print("k w submask_floor Z(w) w^(log2/log3) max_lev floor_beats_w^0.6309")
    for k in range(2, 9):
        wt = (3**k - 1) // 2
        sub = [a for a in binaries(k) if 0 < a < wt and gcd(a, wt) == 1]
        got = [
            a
            for a in range(1, wt)
            if gcd(a, wt) == 1
            and (a % 3 == 0) != ((wt - a) % 3 == 0)
            and level(*((a, wt - a) if a % 3 == 0 else (wt - a, a)))[0]
        ]
        lv = max(level(*((a, wt - a) if a % 3 == 0 else (wt - a, a)))[0] for a in got)
        print(
            k,
            wt,
            len(sub),
            len(got),
            round(wt ** (log(2) / log(3)), 1),
            lv,
            len(sub) > wt ** (log(2) / log(3)),
        )
    R = {}
    for k in range(2, 7):
        s = set()
        for lab in range(3**k):
            u = v = 0
            m = lab
            for i in range(k):
                d = m % 3
                m //= 3
                if d == 1:
                    u += 3**i
                elif d == 2:
                    v += 3**i
            if u and v % 3:
                s.add(u * pow(v, -1, 3**k) % 3**k)
        R[k] = s
    Z = [0] * 729
    recov = miss = 0
    for z1 in range(3, 729, 3):
        for z2 in range(1, 729):
            if z2 % 3 == 0 or gcd(z1, z2) != 1 or z1 + z2 > 728:
                continue
            if not level(z1, z2)[0]:
                continue
            w = z1 + z2
            Z[w] += 2
            k = 1
            while 3**k <= w:
                k += 1
            r = z1 * pow(z2, -1, 3**k) % 3**k
            if r * w * pow(1 + r, -1, 3**k) % 3**k != z1:
                recov += 1
            if 2 <= k <= 6 and r not in R[k]:
                miss += 1
    over = 0
    for w in range(4, 729):
        k = 1
        while 3**k <= w:
            k += 1
        if Z[w] > 2 * len(R[k]):
            over += 1
    print(
        "|R_k| for k = 2..6",
        [len(R[k]) for k in range(2, 7)],
        "recovery failures",
        recov,
        "residues outside R_k",
        miss,
        "weights breaking Z(w) <= 2|R_k|",
        over,
    )
    for k in (7, 8):
        occ, total, _model, _best = lift_scan(k)
        direct = set()
        dtot = 0
        for T in range(0, 1 << k, 2):
            m, o = occ_of(k, T)
            direct |= o
            dtot += len(o)
        print(
            "lift meet-in-the-middle against direct submask enumeration at k",
            k,
            "union",
            len(occ),
            len(direct),
            "sum",
            total,
            dtot,
            "agree",
            occ == direct and total == dtot,
        )
    for t in range(2, 7):
        k = 2 * t + 1
        m, low = cyclotomic_occ(t)
        mt, full = occ_of(k, ((1 << t) - 1) << t)
        print(
            "cyclotomic lift t",
            t,
            "k",
            k,
            "m_T",
            m,
            "= 3^(2t)-3^t+1",
            m == mt,
            "divides 3^(3t)+1",
            (3 ** (3 * t) + 1) % m == 0,
            "cofactor 3^t+1 coprime to R_k",
            gcd(3**t + 1, (3**k - 1) // 2) == 1,
            "predicted",
            len(low),
            "actual",
            len(full),
            "equal",
            low == full,
        )
    for k in (7, 8):
        ms, ns, _tops = raw_scan(k)
        bad = 0
        for j, T in enumerate(range(0, 1 << k, 2)):
            m, v = residue_dist(k, T)
            F = fourier_direct(k, T, m)
            if m != ms[j] or int(v[0]) != ns[j] or abs(F.sum() / m - ns[j]) > 1e-6:
                bad += 1
        fam = all(
            raw_count((p - 1) * t + s, phi2p_T(p, t)[1])[1]
            == 2 ** ((p - 1) * (t - s) // 2 + s)
            for p, t in ((3, 2), (3, 3), (5, 2), (7, 1), (9, 1))
            for s in range(t)
            if (p - 1) * t + s <= 8
        )
        print(
            "Fourier identity N_T = (1/m_T) Sum_u Prod_p (1 + e(u 3^p / m_T)) against the residue DP and the meet-in-the-middle count at k",
            k,
            "over",
            len(ms),
            "sets T, failures",
            bad,
            "antipodal family N_T = 2^((p-1)(t-s)/2 + s) at k <= 8",
            fam,
        )
    rigid = sum(
        1 for k in range(2, 6) for c in column_vectors(k, 3) if len(set(c.tolist())) > 1
    )
    branch = [c.tolist() for c in column_vectors(3, 4) if len(set(c.tolist())) > 1]
    cen = []
    for k in (7, 8, 9):
        fam = len(lift_scan(k)[0] | block3_scan(k, 1) | block3_scan(k, 2))
        Z = len(lift_lengths(k))
        assert fam == Z, (k, fam, Z)
        cen.append((k, fam))
    print(
        "binary multiples of R_k below 3^(3k) carry a constant column vector at k = 2..5, exceptions",
        rigid,
        "and at four blocks the rigidity breaks, non-constant vectors at k = 3",
        len(branch),
        "first",
        branch[0],
        "; the depth-3 family census equals Z(R_k) at",
        cen,
    )
    worst = max(
        sum(
            1 / (1 + 2 * sum(3**i for i in range(k) if T >> i & 1))
            for T in range(0, 1 << k, 2)
        )
        for k in range(2, 13)
    )
    print(
        "Sum_T 1/m_T over T inside [1, k-1] at k = 2..12 stays below 3/2, largest",
        down(worst, 5),
    )
    print("regression A(inf, 3^5) >= 474 and A(9, 3^7) = 2818 are checked by levels")


def main():
    p = argparse.ArgumentParser()
    s = p.add_subparsers(dest="cmd", required=True)
    a = s.add_parser("sweep")
    a.add_argument("cap", type=int)
    a.set_defaults(fn=cmd_sweep)
    b = s.add_parser("levels")
    b.add_argument("cap", type=int)
    b.add_argument("--lo", type=int, default=2)
    b.add_argument("--hi", type=int, default=24)
    b.set_defaults(fn=cmd_levels)
    c = s.add_parser("multiples")
    c.add_argument("n", type=int)
    c.add_argument("--qmax", type=int, default=500)
    c.add_argument("--hmax", type=int, default=8)
    c.set_defaults(fn=cmd_multiples)
    e = s.add_parser("band")
    e.add_argument("cap", type=int)
    e.set_defaults(fn=cmd_band)
    f = s.add_parser("layers")
    f.add_argument("cap", type=int)
    f.set_defaults(fn=cmd_layers)
    g = s.add_parser("weights")
    g.add_argument("w", type=int, nargs="+")
    g.set_defaults(fn=cmd_weights)
    h = s.add_parser("core")
    h.add_argument("w", type=int, nargs="+")
    h.set_defaults(fn=cmd_core)
    i = s.add_parser("box")
    i.add_argument("lo", type=int)
    i.add_argument("hi", type=int)
    i.add_argument("--extra", type=int, default=3)
    i.set_defaults(fn=cmd_box)
    j = s.add_parser("repunit")
    j.add_argument("--kmax", type=int, default=13)
    j.add_argument("--liftmax", type=int, default=15)
    j.set_defaults(fn=cmd_repunit)
    l = s.add_parser("lifts")
    l.add_argument("--kmax", type=int, default=15)
    l.add_argument("--zmax", type=int, default=13)
    l.add_argument("--lo", type=int, default=11)
    l.set_defaults(fn=cmd_lifts)
    q = s.add_parser("tail")
    q.add_argument("--kmax", type=int, default=13)
    q.add_argument("--famax", type=int, default=12)
    q.add_argument("--onemax", type=int, default=13)
    q.set_defaults(fn=cmd_tail)
    n = s.add_parser("agg")
    n.add_argument("--kmax", type=int, default=19)
    n.add_argument("--umax", type=int, default=11)
    n.add_argument("--ulo", type=int, default=5)
    n.add_argument("--directmax", type=int, default=9)
    n.add_argument("--cycmax", type=int, default=5)
    n.add_argument("--lo", type=int, default=11)
    n.set_defaults(fn=cmd_agg)
    d = s.add_parser("check")
    d.set_defaults(fn=cmd_check)
    args = p.parse_args()
    args.fn(args)


main()
