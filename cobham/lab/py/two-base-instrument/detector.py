import math
import sys
import time

NGRID = 6000
WINDOW_FLOOR = 1000
OMEGA_LO, OMEGA_HI, OMEGA_N = 0.5, 14.0, 1351
TOL = 0.01
RATIO = 10.0
COMB = 8
OM3 = 2 * math.pi / math.log(3)
OM5 = 2 * math.pi / math.log(5)

INDEPENDENT = (3, (0, 1), 5, (0, 1, 2))
DEPENDENT = (3, (0, 1), 9, (0, 1, 3))

# ENUMERATION


def walk(L, p, pd, q, qd, sink):
    pw_p = [p**i for i in range(L + 1)]
    pw_q = [1]
    while pw_q[-1] <= pw_p[L]:
        pw_q.append(pw_q[-1] * q)
    allowed = frozenset(qd)
    stack = [(0, 0, len(pw_q) - 1)]
    nodes = 0
    while stack:
        j, a, s = stack.pop()
        nodes += 1
        hi = a + pw_p[L - j] - 1
        t = s
        while t > 0 and a // pw_q[t - 1] == hi // pw_q[t - 1]:
            t -= 1
        ok = True
        for u in range(t, s):
            if (a // pw_q[u]) % q not in allowed:
                ok = False
                break
        if not ok:
            continue
        if j == L:
            sink(a)
            continue
        step = pw_p[L - j - 1]
        for d in reversed(pd):
            stack.append((j + 1, a + d * step, t))
    return nodes


def in_base(n, b, ds):
    if n == 0:
        return 0 in ds
    while n:
        if n % b not in ds:
            return False
        n //= b
    return True


def brute(L, p, pd, q, qd):
    return [n for n in range(p**L) if in_base(n, p, pd) and in_base(n, q, qd)]


def digit_count(N, q, ds):
    if N <= 0:
        return 0
    k = len(ds)
    dig = []
    m = N
    while m:
        dig.append(m % q)
        m //= q
    total = 0
    for i in range(len(dig) - 1, -1, -1):
        total += sum(1 for d in ds if d < dig[i]) * k**i
        if dig[i] not in ds:
            break
    return total - (1 if 0 in ds else 0)


def histogram(L, p, pd, q, qd):
    h = L * math.log(p) / (NGRID - 1)
    cnt = [0] * (NGRID + 2)
    log = math.log

    def sink(n, cnt=cnt, h=h, log=log):
        if n:
            i = int(log(n) / h) + 1
            if i <= NGRID:
                cnt[i] += 1

    nodes = walk(L, p, pd, q, qd, sink)
    run = 0
    out = []
    for i in range(NGRID):
        run += cnt[i]
        out.append(run)
    return h, out, nodes


def one_base_ladder(L, p, b, ds):
    h = L * math.log(p) / (NGRID - 1)
    return h, [digit_count(int(math.exp(i * h)), b, ds) for i in range(NGRID)]


def crop(h, counts):
    start = 0
    while start < len(counts) and counts[start] < WINDOW_FLOOR:
        start += 1
    return [i * h for i in range(start, len(counts))], list(counts[start:])


# DETECTOR


def blackman(m):
    return [
        0.42
        - 0.5 * math.cos(2 * math.pi * j / (m - 1))
        + 0.08 * math.cos(4 * math.pi * j / (m - 1))
        for j in range(m)
    ]


def power_at(y, h, w):
    z = complex(math.cos(w * h), -math.sin(w * h))
    acc = 0j
    for v in reversed(y):
        acc = acc * z + v
    return abs(acc) ** 2


def golden(f, lo, hi):
    phi = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(60):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = f(d)
    return (a + b) / 2


def detrend(us, cs, deg):
    ys = [math.log(c) for c in cs]
    mid = (us[0] + us[-1]) / 2
    half = (us[-1] - us[0]) / 2
    ts = [(x - mid) / half for x in us]
    basis = []
    for k in range(deg + 1):
        v = [t**k for t in ts]
        for b in basis:
            d = sum(a * c for a, c in zip(v, b))
            v = [a - d * c for a, c in zip(v, b)]
        nrm = math.sqrt(sum(a * a for a in v))
        basis.append([a / nrm for a in v])
    g = list(ys)
    for b in basis:
        d = sum(a * c for a, c in zip(g, b))
        g = [a - d * c for a, c in zip(g, b)]
    n = len(ts)
    mt = sum(ts) / n
    my = sum(ys) / n
    slope = (
        sum((t - mt) * (y - my) for t, y in zip(ts, ys))
        / sum((t - mt) ** 2 for t in ts)
        / half
    )
    return slope, g


def spectrum(g, h):
    w = blackman(len(g))
    y = [a * b for a, b in zip(g, w)]
    m = sum(y) / len(y)
    y = [v - m for v in y]
    step = (OMEGA_HI - OMEGA_LO) / (OMEGA_N - 1)
    ws = [OMEGA_LO + i * step for i in range(OMEGA_N)]
    return y, ws, [power_at(y, h, x) for x in ws]


def maxima(y, h, ws, ps):
    out = []
    for i in range(1, len(ps) - 1):
        if ps[i] > ps[i - 1] and ps[i] >= ps[i + 1]:
            x = golden(lambda w: -power_at(y, h, w), ws[i - 1], ws[i + 1])
            out.append((x, power_at(y, h, x)))
    return out


def median(xs):
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


CANDIDATES = [
    ("2pi/ln3", 2 * math.pi / math.log(3)),
    ("2pi/ln5", 2 * math.pi / math.log(5)),
    ("2pi/ln9", 2 * math.pi / math.log(9)),
    ("2pi/ln15", 2 * math.pi / math.log(15)),
    ("beat", 2 * math.pi / math.log(3) - 2 * math.pi / math.log(5)),
    ("sum", 2 * math.pi / math.log(3) + 2 * math.pi / math.log(5)),
]

# READING


def nearest_comb(w):
    best = None
    for m in range(-COMB, COMB + 1):
        for n in range(-COMB, COMB + 1):
            if m == 0 and n == 0:
                continue
            v = m * OM3 + n * OM5
            if v <= 0:
                continue
            e = abs(v - w) / v * 100
            if best is None or e < best[0]:
                best = (e, m, n, v)
    return best


def analyse(us, cs, asked, deg):
    h = us[1] - us[0]
    slope, g = detrend(us, cs, deg)
    y, ws, ps = spectrum(g, h)
    peaks = maxima(y, h, ws, ps)
    med = median(ps)
    loud = [t for t in peaks if t[1] >= RATIO * med]
    cover = sum(2 * TOL * w for w, _ in loud) / (OMEGA_HI - OMEGA_LO)
    top = sorted(peaks, key=lambda t: -t[1])[:4]
    verdict = {}
    for label, target in CANDIDATES:
        if label not in asked:
            continue
        best = min(peaks, key=lambda t: abs(t[0] - target))
        err = abs(best[0] - target) / target * 100
        ratio = best[1] / med
        verdict[label] = (
            err < TOL * 100 and ratio >= RATIO,
            err,
            ratio,
            best[0],
            target,
        )
    return (
        slope,
        max(g) - min(g),
        len(loud),
        cover,
        [(w, p / med) for w, p in top],
        verdict,
    )


def reading(us, cs, asked, deg):
    slope, swing, nloud, cover, top, verdict = analyse(us, cs, asked, deg)
    print(
        "    detrend degree %d   fitted exponent %.6f   swing %.6f   maxima above %.0fx median %d   chance of a hit by position alone %.3f"
        % (deg, slope, swing, RATIO, nloud, cover)
    )
    print("      strongest maxima  " + "   ".join("%.4f (%.3g x)" % t for t in top))
    e, m, n, v = nearest_comb(top[0][0])
    print(
        "      strongest maximum %.6f: nearest m 2pi/ln3 + n 2pi/ln5 with |m|, |n| <= %d is (%d, %d) = %.6f, error %.3f%%"
        % (top[0][0], COMB, m, n, v, e)
    )
    for label, _ in CANDIDATES:
        if label not in verdict:
            continue
        ok, err, ratio, got, target = verdict[label]
        print(
            "      %-8s predicted %.6f  nearest max %.6f  error %.3f%%  power %.3g x median  present %s"
            % (label, target, got, err, ratio, ok)
        )
    return slope, verdict


def report(name, us, cs, asked):
    print("  %s" % name)
    print(
        "    window u in [%.4f, %.4f], span %.4f, %d points, count %s to %s, resolution 2pi/span %.4f"
        % (
            us[0],
            us[-1],
            us[-1] - us[0],
            len(us),
            fmt(cs[0]),
            fmt(cs[-1]),
            2 * math.pi / (us[-1] - us[0]),
        )
    )
    reading(us, cs, asked, 1)
    reading(us, cs, asked, 3)


def fmt(c):
    return "%d" % c if float(c).is_integer() else "%.6g" % c


# VERBS


def verb_control(M):
    p, pd, q, qd = INDEPENDENT
    print(
        "CONTROL: the pruned walk against direct digit filtering, base %d digits %s and base %d digits %s"
        % (p, list(pd), q, list(qd))
    )
    print("L  height        direct  walk  nodes  agree")
    for L in range(M + 1):
        ref = brute(L, p, pd, q, qd)
        got = []
        nodes = walk(L, p, pd, q, qd, got.append)
        print(
            "%-2d %-13d %-7d %-5d %-6d %s"
            % (L, p**L, len(ref), len(got), nodes, ref == got)
        )
        sys.stdout.flush()
    terms = []
    walk(10, p, pd, q, qd, terms.append)
    print(
        "first members of the cell below %d^10: %s"
        % (p, ", ".join(str(n) for n in terms[1:13]))
    )


def pair_counts(L, cell):
    p, pd, q, qd = cell
    t0 = time.time()
    h, counts, nodes = histogram(L, p, pd, q, qd)
    budget = math.log(len(pd), p) + math.log(len(qd), q) - 1
    print(
        "  base %d digits %s against base %d digits %s: budget %.6f rounded up, %d nodes, %d hits, %.1f s"
        % (
            p,
            list(pd),
            q,
            list(qd),
            math.ceil(budget * 1e6) / 1e6,
            nodes,
            counts[-1],
            time.time() - t0,
        )
    )
    sys.stdout.flush()
    return crop(h, counts)


def verb_cell(L):
    p, pd, q, qd = INDEPENDENT
    print("CELL: height %d^%d = %d" % (p, L, p**L))
    print(
        "  criterion: a local maximum of the Blackman periodogram of ln C(e^u) detrended in u, lying within %.0f%% of the prediction and carrying at least %.0fx the median power of the band [%.1f, %.1f]"
        % (TOL * 100, RATIO, OMEGA_LO, OMEGA_HI)
    )
    print("INDEPENDENT PAIR")
    us, cs = pair_counts(L, INDEPENDENT)
    report(
        "C(N) of the intersection",
        us,
        cs,
        ("2pi/ln3", "2pi/ln5", "2pi/ln15", "beat", "sum"),
    )
    half = len(us) // 2
    report(
        "C(N) on the upper half of the window",
        us[half:],
        cs[half:],
        ("2pi/ln3", "2pi/ln5"),
    )
    print("ONE-BASE CONTROLS")
    hp, cp = one_base_ladder(L, p, p, pd)
    up, vp = crop(hp, cp)
    report("base %d digits %s alone" % (p, list(pd)), up, vp, ("2pi/ln3", "2pi/ln5"))
    hq, cq = one_base_ladder(L, p, q, qd)
    uq, vq = crop(hq, cq)
    report("base %d digits %s alone" % (q, list(qd)), uq, vq, ("2pi/ln3", "2pi/ln5"))


def verb_ladder(A, B):
    p, pd, q, qd = INDEPENDENT
    print(
        "LADDER: the criterion at every height from %d^%d to %d^%d, full window, both detrends"
        % (p, A, p, B)
    )
    print(
        "  a height whose new hits are 0 is a decade of %d containing no member of the cell, on which C is exactly constant"
        % p
    )
    print(
        "L  hits      new       span   pts  loud cover  deg  2pi/ln3 err   power  present  2pi/ln5 err   power  present"
    )
    prev = None
    for L in range(A, B + 1):
        h, counts, _nodes = histogram(L, p, pd, q, qd)
        us, cs = crop(h, counts)
        new = "" if prev is None else "%d" % (counts[-1] - prev)
        prev = counts[-1]
        for deg in (1, 3):
            _, _, nloud, cover, _, v = analyse(us, cs, ("2pi/ln3", "2pi/ln5"), deg)
            a, b = v["2pi/ln3"], v["2pi/ln5"]
            print(
                "%-2d %-9d %-9s %6.2f %5d %4d %6.3f %3d %8.3f%% %7.3g %-8s %8.3f%% %7.3g %-8s"
                % (
                    L,
                    counts[-1],
                    new,
                    us[-1] - us[0],
                    len(us),
                    nloud,
                    cover,
                    deg,
                    a[1],
                    a[2],
                    a[0],
                    b[1],
                    b[2],
                    b[0],
                )
            )
            new = ""
            sys.stdout.flush()


def top_band(b, k, ds):
    return b**k, max(ds) * (b ** (k + 1) - 1) // (b - 1)


def verb_blocks(L):
    p, pd, q, qd = INDEPENDENT
    print("BLOCKS: the applicability hypothesis of the detector, height %d^%d" % (p, L))
    print(
        "  a member of the base %d design with k+1 digits lies in [%d^k, (%d^(k+1)-1)/2], so the design occupies the fraction ln(%d/2)/ln %d = %.6f of every decade of %d and the rest of the decade is empty"
        % (p, p, p, p, p, math.log(p / 2) / math.log(p), p)
    )
    print(
        "  the base %d design occupies ln(%d/2)/ln %d = %.6f of every decade of %d"
        % (q, q, q, math.log(q / 2) / math.log(q), q)
    )
    empty = []
    for j in range(L):
        lo3, hi3 = top_band(p, j, pd)
        hit = False
        i = 0
        while q**i <= hi3:
            lo5, hi5 = top_band(q, i, qd)
            if max(lo3, lo5) <= min(hi3, hi5):
                hit = True
                break
            i += 1
        if not hit:
            empty.append(j)
    print(
        "  decades [%d^j, %d^(j+1)) whose two design bands do not meet, so the cell has no member there and C is constant across the whole decade, j = %s"
        % (p, p, ", ".join(str(j) for j in empty))
    )
    print(
        "  the fraction of the exponents j < %d carrying no member is %d / %d = %.4f"
        % (L, len(empty), L, len(empty) / L)
    )
    print(
        "BLOCK MODEL: C3(N) C5(N) / N, the two band structures multiplied with no joint arithmetic"
    )
    h = L * math.log(p) / (NGRID - 1)
    model = []
    for i in range(NGRID):
        N = int(math.exp(i * h))
        model.append(digit_count(N, p, pd) * digit_count(N, q, qd) / N if N > 0 else 0)
    um, cm = crop(h, model)
    report("C3(N) C5(N) / N", um, cm, ("2pi/ln3", "2pi/ln5", "2pi/ln15", "beat", "sum"))


def verb_collapse(L):
    p, _pd, _q, _qd = DEPENDENT
    print("COLLAPSE: height %d^%d = %d" % (p, L, p**L))
    print(
        "  the dependent pair collapses to the base 9 design {0,1,3} of exponent log_9 3 = 0.5, one lattice of frequency 2pi/ln 9 whose first harmonic is 2pi/ln 3"
    )
    ud, cd = pair_counts(L, DEPENDENT)
    report(
        "C(N) of the dependent intersection", ud, cd, ("2pi/ln3", "2pi/ln5", "2pi/ln9")
    )


def main():
    args = sys.argv[1:]
    if len(args) == 3 and args[0] == "ladder":
        verb_ladder(int(args[1]), int(args[2]))
    elif len(args) == 2 and args[0] == "blocks":
        verb_blocks(int(args[1]))
    elif len(args) == 2 and args[0] == "collapse":
        verb_collapse(int(args[1]))
    elif len(args) == 2 and args[0] == "control":
        verb_control(int(args[1]))
    elif len(args) == 2 and args[0] == "cell":
        verb_cell(int(args[1]))
    else:
        print("verbs: control M | cell L | collapse M | ladder A B | blocks L")


main()
