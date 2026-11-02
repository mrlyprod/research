import argparse
import time
from math import ceil, floor, gcd, log

import numpy as np

LOW = 0.4475978
HIGH = 0.6402122


def up(x, d):
    return ceil(x * 10**d) / 10**d


def down(x, d):
    return floor(x * 10**d) / 10**d


# GASKET CHUNKS


def base_pair(m, fix0):
    u = np.zeros(1, dtype=np.int64)
    v = np.zeros(1, dtype=np.int64)
    start = 0
    if fix0:
        v = np.ones(1, dtype=np.int64)
        start = 1
    for i in range(start, m):
        p = 3**i
        u = np.concatenate([u, u + p, u])
        v = np.concatenate([v, v, v + p])
    return u, v


def top_offsets(lo, hi):
    out = [(0, 0)]
    for i in range(lo, hi):
        p = 3**i
        out = (
            [(a, b) for a, b in out]
            + [(a + p, b) for a, b in out]
            + [(a, b + p) for a, b in out]
        )
    return out


def chunks(n, fix0, m):
    m = min(m, n)
    bu, bv = base_pair(m, fix0)
    for du, dv in top_offsets(m, n):
        yield bu + du, bv + dv


# RATIO SET


def ratio_values(k, m):
    mod = 3**k
    for u, v in chunks(k, True, m):
        keep = u > 0
        u = u[keep]
        v = v[keep]
        vs = np.unique(v)
        iv = np.array([pow(int(x), -1, mod) for x in vs], dtype=np.int64)
        yield (u * iv[np.searchsorted(vs, v)]) % mod


def ratio_stats(k, m, full):
    mod = 3**k
    if full:
        hist = np.zeros(mod, dtype=np.int32)
        for r in ratio_values(k, m):
            hist += np.bincount(r, minlength=mod).astype(np.int32)
        size = int(np.count_nonzero(hist))
        h = hist.astype(np.int64)
        return size, int((h * h).sum()), int(h.max())
    parts = []
    for r in ratio_values(k, m):
        parts.append(np.unique(r))
        if len(parts) > 24:
            parts = [np.unique(np.concatenate(parts))]
    return int(np.unique(np.concatenate(parts)).size), 0, 0


def cmd_ratios(args):
    print("k R_k sigma_k growth c_k k*c_k M2/4^k mean CSlower top secs")
    prev = 0
    band = []
    for k in range(2, args.kmax + 1):
        t = time.time()
        size, second, top = ratio_stats(k, args.mem, k <= args.full)
        pairs = 3 ** (k - 1) - 2 ** (k - 1)
        growth = size / prev if prev else 0.0
        ck = 1 - log(growth) / log(3) if growth else 0.0
        if growth:
            band.append((k, ck))
        prev = size
        row = [
            k,
            size,
            round(size / 3**k, 6),
            round(growth, 4),
            round(ck, 7),
            round(k * ck, 7),
        ]
        if second:
            row += [
                round(second / 4**k, 4),
                round(pairs / size, 4),
                round(pairs * pairs / second / 3**k, 6),
                top,
            ]
        print(*row, round(time.time() - t, 1))
    if len(band) > 1:
        prod = [k * c for k, c in band[-6:]]
        print(
            "decay band k*c_k",
            down(min(prod), 4),
            up(max(prod), 4),
            "over k",
            band[-6:][0][0],
            band[-1][0],
        )


# CENSUS

POW3 = np.array([3**i for i in range(39)], dtype=np.int64)


def census(n, m, cut_pow, alphas, caps=()):
    thr = np.array([3.0 ** (a * n) for a in alphas] + [3.0**j for j in caps])
    keys = []
    pts = 0
    fibre = 0
    weighted = np.zeros(n + 2, dtype=np.int64)
    heavy = np.zeros(thr.size, dtype=np.int64)
    cut = max(3.0**cut_pow, float(thr.max()) if thr.size else 0.0)
    for u, v in chunks(n, False, m):
        ok = (u > 0) & (v > 0)
        u = u[ok]
        v = v[ok]
        fibre += int(ok.size - u.size)
        g = np.gcd(u, v)
        z1 = u // g
        z2 = v // g
        h = np.maximum(z1, z2)
        pts += h.size
        oct_ = np.searchsorted(POW3, h, side="right") - 1
        weighted += np.bincount(oct_, minlength=n + 2)[: n + 2]
        for i, t in enumerate(thr):
            heavy[i] += int(np.count_nonzero(h <= t))
        sel = h <= cut
        if sel.any():
            keys.append(np.unique((z1[sel] << np.int64(32)) + z2[sel]))
        if len(keys) > 48:
            keys = [np.unique(np.concatenate(keys))]
    keys = np.unique(np.concatenate(keys)) if keys else np.zeros(0, dtype=np.int64)
    kh = np.maximum(keys >> np.int64(32), keys & np.int64((1 << 32) - 1))
    occ = [int(np.count_nonzero(kh <= t)) for t in thr]
    return pts, fibre, weighted, heavy, occ, keys, kh


def cmd_census(args):
    alphas = args.alphas
    caps = [j for j in args.caps]
    print(
        "convention height max(z1,z2), window height <= 3^(alpha n), octave floor(log_3 h),"
        " ladder octave = octave + 1 above h = 1, ray totals exclude the two fibre rays"
    )
    print("n alpha A F meanM expA expF theta box A/box")
    seen = {}
    for n in args.levels:
        t = time.time()
        cutp = min(n, int(args.cut * n) + 1)
        pts, fibre, weighted, heavy, occ, keys, _kh = census(
            n, args.mem, cutp, alphas, caps
        )
        for i, a in enumerate([str(x) for x in alphas] + ["3^%d" % j for j in caps]):
            x = (
                3.0 ** (alphas[i] * n)
                if i < len(alphas)
                else 3.0 ** caps[i - len(alphas)]
            )
            e_a = log(occ[i]) / (n * log(3)) if occ[i] else 0.0
            th = log(occ[i]) / log(x) if occ[i] else 0.0
            box = int(x) ** 2
            seen.setdefault(a, []).append((n, e_a, th))
            print(
                n,
                a,
                occ[i],
                int(heavy[i]),
                round(heavy[i] / occ[i], 4) if occ[i] else 0.0,
                round(e_a, 4),
                round(log(int(heavy[i])) / (n * log(3)), 4) if heavy[i] else 0.0,
                round(th, 4),
                box,
                round(occ[i] / box, 6),
            )
        print(
            "level",
            n,
            "nonfibre",
            pts,
            "fibre",
            fibre,
            "cut",
            cutp,
            "keys",
            keys.size,
            "total_rays",
            keys.size if cutp >= n else 0,
            keys.size + 2 if cutp >= n else 0,
            round(time.time() - t, 1),
        )
        print("octaves", n, *[int(x) for x in weighted[: n + 1]])
    for a, rows in seen.items():
        print(
            "band",
            a,
            "expA",
            down(min(r[1] for r in rows), 4),
            up(max(r[1] for r in rows), 4),
            "theta",
            down(min(r[2] for r in rows), 4),
            up(max(r[2] for r in rows), 4),
            "over n",
            rows[0][0],
            rows[-1][0],
        )


# LARGE PRIME SUM


def spf_sieve(limit):
    spf = np.zeros(limit + 1, dtype=np.int32)
    spf[2::2] = 2
    for p in range(3, int(limit**0.5) + 1, 2):
        if spf[p] == 0:
            spf[p * p :: 2 * p] = np.where(
                spf[p * p :: 2 * p] == 0, p, spf[p * p :: 2 * p]
            )
    odd = np.arange(3, limit + 1, 2)
    spf[3::2] = np.where(spf[3::2] == 0, odd, spf[3::2])
    spf[1] = 1
    return spf


def gcd_hist(n, m):
    hist = np.zeros(3**n, dtype=np.int32)
    for u, v in chunks(n, False, m):
        g = np.gcd(u, v)
        hist += np.bincount(g, minlength=3**n).astype(np.int32)
    hist[0] = 0
    return hist


def cmd_sieve(args):
    print("n beta primesum F bound ratio")
    for n in args.levels:
        t = time.time()
        hist = gcd_hist(n, args.mem)
        spf = spf_sieve(3**n)
        vals = np.nonzero(hist)[0]
        mult = hist[vals].astype(np.int64)
        for beta in args.betas:
            thr = 3.0 ** (beta * n)
            work = vals.copy()
            big = np.zeros(work.size, dtype=np.int64)
            last = np.zeros(work.size, dtype=np.int64)
            while True:
                live = work > 1
                if not live.any():
                    break
                p = spf[work].astype(np.int64)
                new = live & (p > thr) & (p != last)
                big += new
                last = np.where(live, p, last)
                work = np.where(live, work // np.maximum(p, 1), work)
            total = int((big * mult).sum())
            _, _, _, heavy, _occ, _, _ = census(n, args.mem, 1, [1.0 - beta])
            bound = (int(heavy[0]) + 2 ** (n + 1)) / beta
            print(
                n, beta, total, int(heavy[0]), round(bound, 1), round(total / bound, 4)
            )
        print("level", n, "secs", round(time.time() - t, 1))


# CHECKS


def brute_ratios(k):
    mod = 3**k
    seen = set()
    for lab in range(3**k):
        u = v = 0
        w = lab
        for i in range(k):
            d = w % 3
            w //= 3
            if d == 1:
                u += 3**i
            elif d == 2:
                v += 3**i
        if u and v % 3:
            seen.add(u * pow(v, -1, mod) % mod)
    return len(seen)


def brute_rays(n, xcap):
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
            if max(z) <= xcap:
                out[z] = out.get(z, 0) + 1
    return out


def cmd_constants(args):
    c = log(4 / 3) / log(3)
    print("window edges", LOW, HIGH)
    print("first moment moves the window above alpha", up(1 - HIGH, 7))
    print("first moment closes the window at alpha", up(1 - LOW, 7))
    print("congruence decay cap c <=", up(c, 7))
    print("congruence alpha cap <=", up(1 / (2 - c), 6))
    print(
        "trivial box C threshold reading", 1, "octave reading", 9, "delta 1 - 2 alpha"
    )
    print("O needs theta <", down(1 / 0.5533, 4), "at alpha 0.5533, box theta 2")


def cmd_check(args):
    print("k brute numpy CSlower ok")
    for k in range(2, 10):
        a = brute_ratios(k)
        b, second, _ = ratio_stats(k, 13, True)
        floor = (3 ** (k - 1) - 2 ** (k - 1)) ** 2 / second
        print(k, a, b, round(floor, 2), a == b and b >= floor)
    reg = census(9, 13, min(9, int(0.62 * 9) + 1), [], [7])
    print(
        "regression A(9, 3^7)",
        reg[4][0],
        len(brute_rays(9, 3**7)),
        reg[4][0] == 2818 == len(brute_rays(9, 3**7)),
    )
    print("n bruteA A bruteF F onecoord3 weight3 nonfibre")
    for n in range(4, 11):
        cutp = min(n, int(0.62 * n) + 1)
        a = 0.5
        rays = brute_rays(n, int(3.0 ** (a * n)))
        pts, _fibre, _weighted, heavy, occ, _keys, _kh = census(n, 13, cutp, [a])
        div3 = all((z[0] % 3 == 0) != (z[1] % 3 == 0) for z in rays)
        w3 = all((z[0] + z[1]) % 3 for z in rays)
        print(
            n,
            len(rays),
            occ[0],
            sum(rays.values()),
            int(heavy[0]),
            div3,
            w3,
            pts == 3**n - 2 ** (n + 1) + 1,
        )


def main():
    p = argparse.ArgumentParser()
    s = p.add_subparsers(dest="cmd", required=True)
    r = s.add_parser("ratios")
    r.add_argument("kmax", type=int)
    r.add_argument("--full", type=int, default=15)
    r.add_argument("--mem", type=int, default=13)
    r.set_defaults(fn=cmd_ratios)
    c = s.add_parser("census")
    c.add_argument("levels", type=int, nargs="+")
    c.add_argument("--cut", type=float, default=0.62)
    c.add_argument("--mem", type=int, default=13)
    c.add_argument("--alphas", type=float, nargs="+", default=[0.45, 0.5, 0.5533, 0.6])
    c.add_argument("--caps", type=int, nargs="*", default=[5, 6, 7])
    c.set_defaults(fn=cmd_census)
    q = s.add_parser("sieve")
    q.add_argument("levels", type=int, nargs="+")
    q.add_argument("--betas", type=float, nargs="+", default=[0.45, 0.5])
    q.add_argument("--mem", type=int, default=13)
    q.set_defaults(fn=cmd_sieve)
    n = s.add_parser("constants")
    n.set_defaults(fn=cmd_constants)
    k = s.add_parser("check")
    k.set_defaults(fn=cmd_check)
    a = p.parse_args()
    a.fn(a)


main()
