import sys

import numpy as np
from mpmath import findroot, mp, mpc, mpf, power, zeta, zetazero
from scipy.ndimage import median_filter
from scipy.stats import spearmanr

SAMPLES = 32768
BAND = (4.0, 60.0)
SCORE = 8.0
FLOOR_WIDTH = 101
HEIGHT = 62.0
ZEROS = 20
LINES = 14
SHIFTS = 4000
SEED = 20260907
LADDER = (14, 16, 18, 20, 22)
LADDER_CENSUS = {
    14: (11, 105),
    16: (149, 173),
    18: (-30, 312),
    20: (496, 539),
    22: (1009, 1089),
}
DRAWS = 24
SIGNS = 8
DETAIL = ("3-01", "10-x9")
FIELD = ((100000000, 10, (0, 1, 2, 3, 4, 5, 6, 7, 8)), (129140163, 3, (0, 1)))
ESCAPE = ((10000000, 16, (0, 1)), (10000000, 10, (0, 1)))
ESCAPE_DECADES = (3, 4, 5, 6, 7)

DESIGNS = {
    "3-01": (3, (0, 1), 20, 3),
    "3-02": (3, (0, 2), 20, 3),
    "3-12": (3, (1, 2), 18, 3),
    "4-01": (4, (0, 1), 18, 4),
    "4-012": (4, (0, 1, 2), 12, 4),
    "5-01": (5, (0, 1), 16, 5),
    "5-012": (5, (0, 1, 2), 12, 5),
    "9-012": (9, (0, 1, 2), 11, 3),
    "9-0123": (9, (0, 1, 2, 3), 10, 3),
    "9-0134": (9, (0, 1, 3, 4), 10, 3),
    "10-x9": (10, (9,), 100000000, 5),
}

CENSUS = {
    "3-01": (496, 539),
    "3-02": (-382, 485),
    "3-12": (-1461, 1582),
    "10-x9": (2181, 5234),
}

QUADRATIC = {3: {1: 1, 2: -1}, 4: {1: 1, 3: -1}, 5: {1: 1, 2: -1, 3: -1, 4: 1}}

CONTROL = 100000000
CONTROL_ANCHOR = 1928


def primes_upto(limit):
    if limit < 2:
        return np.zeros(0, dtype=np.int64)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    return np.flatnonzero(sieve).astype(np.int64)


def is_prime(n):
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def mobius_values(values):
    mu = np.ones(values.size, dtype=np.int8)
    rem = values.astype(np.int64).copy()
    alive = np.flatnonzero(rem > 1)
    limit = round(float(values.max()) ** (1.0 / 3.0)) + 2
    for p in primes_upto(limit):
        if alive.size == 0:
            break
        hit = alive[rem[alive] % p == 0]
        if hit.size:
            rem[hit] //= p
            mu[hit] = -mu[hit]
            again = hit[rem[hit] % p == 0]
            if again.size:
                mu[again] = 0
                while again.size:
                    rem[again] //= p
                    again = again[rem[again] % p == 0]
        block = rem[alive]
        prime = alive[(block > 1) & (block < p * p)]
        mu[prime] = -mu[prime]
        alive = alive[rem[alive] >= p * p]
    for i in alive:
        c = int(rem[i])
        root = int(c**0.5)
        while root * root > c:
            root -= 1
        while (root + 1) * (root + 1) <= c:
            root += 1
        if root * root == c:
            mu[i] = 0
        elif is_prime(c):
            mu[i] = -mu[i]
    return mu


def elements(q, digits, length):
    lead = np.array([d for d in digits if d > 0], dtype=np.int64)
    tail = np.array(digits, dtype=np.int64)
    current = lead
    out = [lead]
    for _ in range(length - 1):
        current = (current[:, None] * q + tail[None, :]).ravel()
        out.append(current)
    values = np.concatenate(out)
    values.sort()
    return values


def design(q, digits, length):
    values = elements(q, digits, length)
    return values, mobius_values(values)


def mobius_sieve(n):
    rem = np.arange(n + 1, dtype=np.int32)
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes_upto(int(n**0.5) + 1):
        rem[p::p] //= p
        mu[p::p] = -mu[p::p]
        square = p * p
        if square <= n:
            mu[square::square] = 0
            power_of = square
            while power_of <= n:
                rem[power_of::power_of] //= p
                power_of *= p
    tail = rem > 1
    mu[tail] = -mu[tail]
    del rem
    return mu


def control_series(n, samples):
    mu = mobius_sieve(n)
    grid = np.exp(np.linspace(np.log(2.0), np.log(float(n)), samples))
    index = np.floor(grid).astype(np.int64)
    uniq, inverse = np.unique(index, return_inverse=True)
    starts = np.concatenate((np.ones(1, dtype=np.int64), uniq[:-1] + 1))
    running = np.cumsum(np.add.reduceat(mu, starts).astype(np.int64))
    del mu
    log_x = np.log(grid)
    return log_x, running[inverse] / np.exp(0.5 * log_x), int(running[-1])


def design_series(values, running, exponent, samples, window=None):
    lo = np.log(float(values[0])) if window is None else window[0]
    hi = np.log(float(values[-1])) if window is None else window[1]
    log_x = np.linspace(lo, hi, samples)
    slot = np.searchsorted(values, np.exp(log_x), side="right") - 1
    inside = slot >= 0
    out = np.zeros(samples, dtype=np.float64)
    out[inside] = running[slot[inside]]
    return log_x, out / np.exp(exponent * log_x)


def power_spectrum(log_x, series):
    step = (log_x[-1] - log_x[0]) / (log_x.size - 1)
    centred = series - series.mean()
    amplitude = np.fft.rfft(centred * np.hanning(log_x.size))
    gamma = np.fft.rfftfreq(log_x.size, d=step) * 2.0 * np.pi
    return gamma, np.abs(amplitude) ** 2


def scored(gamma, power):
    floor = median_filter(power, size=FLOOR_WIDTH, mode="nearest")
    return power / np.maximum(floor, np.finfo(float).tiny)


def window_max(score):
    return np.maximum.reduce([np.roll(score, 1), score, np.roll(score, -1)])


def band_slice(gamma):
    return np.flatnonzero((gamma > BAND[0]) & (gamma < BAND[1]))


def shift_test(gamma, score, targets, rng):
    inside = band_slice(gamma)
    lo, hi = float(gamma[inside[0]]), float(gamma[inside[-1]])
    span = hi - lo
    picked = np.array([t for t in targets if lo < t < hi], dtype=np.float64)
    if picked.size == 0:
        return 0, 0.0, 0.0, 1.0
    best = window_max(score)

    def statistic(offset):
        moved = lo + np.mod(picked - lo + offset, span)
        slot = np.clip(np.searchsorted(gamma, moved), 1, gamma.size - 2)
        return float(np.mean(np.log10(best[slot])))

    value = statistic(0.0)
    null = np.array([statistic(o) for o in rng.uniform(0.0, span, SHIFTS)])
    p = float((null >= value).mean())
    return picked.size, value, float(null.mean()), max(p, 1.0 / SHIFTS)


def peaks(gamma, score):
    inside = (gamma > BAND[0]) & (gamma < BAND[1])
    rise = np.zeros(score.size, dtype=bool)
    rise[1:-1] = (score[1:-1] > score[:-2]) & (score[1:-1] > score[2:])
    return np.flatnonzero(inside & rise & (score > SCORE))


def base_rate(gamma, score):
    inside = band_slice(gamma)
    return float((window_max(score)[inside] > SCORE).mean())


def zeta_zeros():
    mp.dps = 15
    out = []
    n = 1
    while True:
        g = float(zetazero(n).imag)
        if g > HEIGHT:
            return out
        out.append(g)
        n += 1


def dirichlet_l(s, q):
    table = QUADRATIC[q]
    return power(q, -s) * sum(c * zeta(s, mpf(a) / q) for a, c in table.items())


def l_zeros(q):
    mp.dps = 20
    grid = np.arange(0.5, HEIGHT, 0.05)
    mods = [abs(dirichlet_l(mpc(0.5, t), q)) for t in grid]
    out = []
    for i in range(1, grid.size - 1):
        if mods[i] < mods[i - 1] and mods[i] < mods[i + 1] and mods[i] < 0.3:
            root = findroot(lambda s: dirichlet_l(s, q), mpc(0.5, grid[i]))
            out.append(float(root.imag))
    return out


def lattice(q):
    return [2.0 * np.pi * j / np.log(q) for j in range(1, LINES + 1)]


def families(q, alpha, conductor, zeros, cache):
    if conductor not in cache:
        cache[conductor] = l_zeros(conductor)
    lz = cache[conductor]
    return (
        ("zeta", zeros),
        ("zeta+2", [g + 2.0 for g in zeros]),
        ("zetaflip", [64.0 - g for g in zeros]),
        ("alphazeta", [alpha * g for g in zeros]),
        (f"L{conductor}", lz),
        ("lattice", lattice(q)),
    )


def split_correlation(values, running, exponent):
    lo = np.log(float(values[0]))
    hi = np.log(float(values[-1]))
    mid = 0.5 * (lo + hi)
    halves = []
    for window in ((lo, mid), (mid, hi)):
        log_x, series = design_series(values, running, exponent, SAMPLES // 2, window)
        gamma, power = power_spectrum(log_x, series)
        halves.append((gamma, scored(gamma, power)))
    inside = band_slice(halves[0][0])
    rho, _ = spearmanr(halves[0][1][inside], halves[1][1][inside])
    return float(rho), float(halves[0][0][1])


def split_null(values, support, exponent, rng, draws):
    out = []
    for _ in range(draws):
        signs = np.zeros(values.size, dtype=np.int32)
        signs[support] = rng.choice([-1, 1], int(support.sum()))
        out.append(
            split_correlation(values, np.cumsum(signs, dtype=np.int64), exponent)[0]
        )
    return float(np.mean(out)), float(np.std(out))


def digit_mask(n, base, missing):
    mask = np.ones(n + 1, dtype=bool)
    place = 1
    while place <= n:
        block = base * place
        size = ((n + 1) // block) * block
        view = mask[:size].reshape(-1, base, place)
        view[:, missing, :] = False
        tail = mask[size:]
        if tail.size:
            offset = missing * place
            if offset < tail.size:
                tail[offset : offset + place] = False
        place = block
    mask[0] = False
    return mask


def digit_set_mask(n, base, digits):
    mask = np.ones(n + 1, dtype=bool)
    for d in range(base):
        if d not in digits:
            mask &= digit_mask(n, base, d)
    mask[0] = False
    return mask


def grid_running(weights, grid):
    index = np.floor(grid).astype(np.int64)
    uniq, inverse = np.unique(index, return_inverse=True)
    starts = np.concatenate((np.ones(1, dtype=np.int64), uniq[:-1] + 1))
    out = []
    for w in weights:
        out.append(np.cumsum(np.add.reduceat(w, starts).astype(np.int64))[inverse])
    return out


def mean_field(n, base, digits, samples):
    alpha = np.log(len(digits)) / np.log(base)
    mu = mobius_sieve(n)
    mask = digit_set_mask(n, base, digits)
    grid = np.exp(np.linspace(np.log(2.0), np.log(float(n)), samples))
    log_x = np.log(grid)
    uniq, inverse = np.unique(np.floor(grid).astype(np.int64), return_inverse=True)
    bounds = np.concatenate((np.ones(1, dtype=np.int64), uniq + 1))
    echo = np.empty(uniq.size)
    total = 0.0
    seen = 0
    for j in range(uniq.size):
        lo, hi = bounds[j], bounds[j + 1]
        run = seen + np.cumsum(mask[lo:hi], dtype=np.int64)
        line = np.arange(lo, hi, dtype=np.float64)
        total += float(np.dot(mu[lo:hi].astype(np.float64), run / line))
        seen = int(run[-1]) if run.size else seen
        echo[j] = total
    meter = np.cumsum(np.add.reduceat(mu * mask, bounds[:-1]))
    del mu, mask
    scale = np.exp(alpha * 0.5 * log_x)
    return (
        log_x,
        alpha,
        meter[inverse] / scale,
        echo[inverse] / scale,
        (meter[inverse] - echo[inverse]) / scale,
        int(meter[-1]),
        seen,
    )


def upper_rms(series):
    half = series[series.size // 2 :]
    return float(np.sqrt(np.mean(half * half)))


def kempner(n, base, missing):
    mu = mobius_sieve(n)
    mask = digit_mask(n, base, missing)
    values = np.flatnonzero(mask).astype(np.int64)
    weights = mu[values].copy()
    del mu, mask
    return values, weights


def line(tag, name, count, value, null, p):
    print(
        f"  {tag:8s} {name:8s} n {count:3d} logscore {value:6.3f} null {null:6.3f} p {p:.4f}"
    )


def verb_sieve():
    for name in DESIGNS:
        q, alpha, _, values, mu = sources(name)
        running = np.cumsum(mu, dtype=np.int64)
        top = int(running[-1])
        peak = int(np.abs(running).max())
        span = np.log(float(values[-1])) - np.log(float(values[0]))
        print(
            f"design {name:8s} q {q:2d} alpha {alpha:.6f} A {values.size:8d}"
            f" M {top:6d} Mmax {peak:6d} thetamax {np.log(peak) / np.log(values.size):.4f}"
            f" support {int((mu != 0).sum()):8d} logspan {span:.4f} bin {2 * np.pi / span:.4f}"
        )
        want = CENSUS.get(name)
        if want:
            print(
                f"anchor {name:8s} census {want} meter {(top, peak)} match {want == (top, peak)}"
            )
        del values, mu, running
    lhs, lmu = design(9, (0, 1, 3, 4), 10)
    rhs, rmu = design(3, (0, 1), 20)
    same = bool(np.array_equal(lhs, rhs)) and bool(np.array_equal(lmu, rmu))
    print(f"identity base 9 F 0134 equals base 3 F 01 at 3^20 {same}")
    _, _, anchor = control_series(CONTROL, 1024)
    print(
        f"control base 10 full set M(10^8) {anchor} A084237 {CONTROL_ANCHOR} match {anchor == CONTROL_ANCHOR}"
    )


def sources(name):
    q, digits, deep, conductor = DESIGNS[name]
    if q == 10:
        values, mu = kempner(deep, 10, digits[0])
        alpha = np.log(9.0) / np.log(10.0)
    else:
        values, mu = design(q, digits, deep)
        alpha = np.log(len(digits)) / np.log(q)
    return q, alpha, conductor, values, mu


def census_row(name, rng, zeros, cache, detail=False):
    q, alpha, conductor, values, mu = sources(name)
    support = mu != 0
    running = np.cumsum(mu, dtype=np.int64)
    top = int(running[-1])
    peak = int(np.abs(running).max())
    log_x, series = design_series(values, running, alpha * 0.5, SAMPLES)
    gamma, power = power_spectrum(log_x, series)
    score = scored(gamma, power)
    rho, _ = split_correlation(values, running, alpha * 0.5)
    mean, sigma = split_null(
        values, support, alpha * 0.5, rng, DRAWS if values.size < 5000000 else 6
    )
    z = (rho - mean) / sigma if sigma > 0 else 0.0
    cells = []
    for label, targets in families(q, alpha, conductor, zeros, cache):
        _, value, null, p = shift_test(gamma, score, targets, rng)
        cells.append(f"{label} {value:6.3f}/{null:5.3f} p {p:.4f}")
    print(
        f"census {name:8s} q {q:2d} alpha {alpha:.6f} A {values.size:8d}"
        f" M {top:6d} Mmax {peak:6d} bin {gamma[1]:.4f} rate {base_rate(gamma, score):.4f}"
        f" split {rho:7.4f} z {z:5.2f}"
    )
    print("  meter    " + "  ".join(cells))
    if detail:
        for tag, weight in (
            ("count", np.ones(values.size, dtype=np.int64)),
            ("random", None),
        ):
            if weight is None:
                weight = np.zeros(values.size, dtype=np.int32)
                weight[support] = rng.choice([-1, 1], int(support.sum()))
            run = np.cumsum(weight)
            exponent = alpha if tag == "count" else alpha * 0.5
            lx, sr = design_series(values, run, exponent, SAMPLES)
            g, pw = power_spectrum(lx, sr)
            sc = scored(g, pw)
            cells = []
            for label, targets in families(q, alpha, conductor, zeros, cache):
                _, value, null, p = shift_test(g, sc, targets, rng)
                cells.append(f"{label} {value:6.3f}/{null:5.3f} p {p:.4f}")
            print(f"  {tag:8s} " + "  ".join(cells))
            print(
                f"  {tag:8s} split {split_correlation(values, run, exponent)[0]:7.4f}"
            )
        found = peaks(gamma, score)
        for i in sorted(found[np.argsort(-score[found])][:6]):
            row = f"  peak {gamma[i]:8.4f} score {score[i]:8.2f}"
            for label, targets in families(q, alpha, conductor, zeros, cache):
                near = min(targets, key=lambda z: abs(z - gamma[i]))
                row += f" {label} {near:7.3f} off {abs(near - gamma[i]):5.3f}"
            print(row)
    return top, peak


def escape_line(n, base, digits):
    log_x, alpha, _, echo, _, top, count = mean_field(n, base, digits, SAMPLES)
    raw = echo * np.exp(alpha * 0.5 * log_x)
    print(
        f"escape base {base} F {''.join(str(d) for d in digits)} N {n} alpha {alpha:.6f}"
        f" A_F {count} M_F {top} echo rms {upper_rms(echo):.6f}"
    )
    for decade in ESCAPE_DECADES:
        j = min(int(np.searchsorted(log_x, decade * np.log(10.0))), log_x.size - 1)
        old_bound = float(np.exp((alpha - 0.5) * log_x[j]))
        print(
            f"  echo     x 1e{decade} sum mu A_F over n {float(raw[j]):8.4f}"
            f" old bound x^(alpha-1/2) {old_bound:7.4f} ratio {abs(float(raw[j])) / old_bound:6.2f}"
        )
    print()


def verb_spectrum():
    rng = np.random.default_rng(SEED)
    zeros = zeta_zeros()
    cache = {}
    log_x, series, anchor = control_series(CONTROL, SAMPLES)
    gamma, power = power_spectrum(log_x, series)
    score = scored(gamma, power)
    print(
        f"control base 10 full set N {CONTROL} M {anchor} A084237 {CONTROL_ANCHOR}"
        f" match {anchor == CONTROL_ANCHOR} logspan {log_x[-1] - log_x[0]:.4f}"
        f" bin {gamma[1]:.4f} rate {base_rate(gamma, score):.4f}"
    )
    line("control", "zeta", *shift_test(gamma, score, zeros, rng))
    found = peaks(gamma, score)
    order = found[np.argsort(-score[found])][:10]
    matched = 0
    for i in sorted(order):
        near = min(zeros, key=lambda z: abs(z - gamma[i]))
        matched += abs(near - gamma[i]) < gamma[1]
        print(
            f"  peak {gamma[i]:8.4f} score {score[i]:11.1f} zeta {near:8.4f} off {abs(near - gamma[i]):5.3f}"
        )
    print(f"  control top 10 peaks within one bin of a zeta zero {matched} of 10")
    print()
    for name in DESIGNS:
        got = census_row(name, rng, zeros, cache, detail=name in DETAIL)
        want = CENSUS.get(name)
        if want:
            print(f"anchor {name:8s} census {want} meter {got} match {want == got}")
        print()
    for n, base, digits in FIELD:
        log_x, alpha, meter, echo, residual, top, count = mean_field(
            n, base, digits, SAMPLES
        )
        print(
            f"field  base {base} F {''.join(str(d) for d in digits)} N {n} alpha {alpha:.6f}"
            f" A_F {count} M_F {top} bin {2 * np.pi / (log_x[-1] - log_x[0]):.4f}"
        )
        for tag, series in (("meter", meter), ("echo", echo), ("residual", residual)):
            g, pw = power_spectrum(log_x, series)
            sc = scored(g, pw)
            _, value, null, p = shift_test(g, sc, zeros, rng)
            top6 = peaks(g, sc)
            top6 = top6[np.argsort(-sc[top6])][:6]
            near = sum(1 for i in top6 if min(abs(z - g[i]) for z in zeros) < g[1])
            print(
                f"  {tag:9s} zeta {value:6.3f}/{null:5.3f} p {p:.4f}"
                f" rms {upper_rms(series):8.4f} top6 within one bin {near} of {top6.size}"
            )
        cut = log_x.size
        base_ratio = None
        for fraction in (0.6, 0.8, 1.0):
            take = int(cut * fraction)
            ratio = upper_rms(echo[:take]) / upper_rms(meter[:take])
            span = log_x[take - 1]
            if base_ratio is None:
                base_ratio, base_span = ratio, span
            rate = min(alpha, 1.0 - alpha)
            predicted = base_ratio * np.exp(-0.5 * rate * (span - base_span))
            print(
                f"  decay    logx {span:7.4f} echo over meter {ratio:.6f}"
                f" predicted {predicted:.6f} ratio {ratio / predicted:.4f}"
            )
        print()
    for n, base, digits in ESCAPE:
        escape_line(n, base, digits)
    for deep in LADDER:
        values, mu = design(3, (0, 1), deep)
        alpha = np.log(2.0) / np.log(3.0)
        running = np.cumsum(mu, dtype=np.int64)
        log_x, series = design_series(values, running, alpha * 0.5, SAMPLES)
        gamma, power = power_spectrum(log_x, series)
        score = scored(gamma, power)
        _, value, null, p = shift_test(gamma, score, zeros, rng)
        rho, _ = split_correlation(values, running, alpha * 0.5)
        mean, sigma = split_null(values, mu != 0, alpha * 0.5, rng, DRAWS)
        got = (int(running[-1]), int(np.abs(running).max()))
        want = LADDER_CENSUS.get(deep)
        print(
            f"ladder 3-01 L {deep:2d} A {values.size:8d} bin {gamma[1]:.4f}"
            f" zeta {value:6.3f}/{null:5.3f} p {p:.4f}"
            f" split {rho:7.4f} z {(rho - mean) / sigma if sigma > 0 else 0.0:5.2f}"
            f" census {want} meter {got} match {want == got}"
        )


def verb_family():
    rng = np.random.default_rng(SEED)
    zeros = zeta_zeros()
    cache = {}
    for left, right in (
        ("4-01", "9-012"),
        ("3-01", "9-0123"),
        ("3-01", "9-0134"),
        ("3-01", "3-02"),
        ("10-x9", "3-01"),
    ):
        rows = []
        for name in (left, right):
            q, alpha, conductor, values, mu = sources(name)
            log_x, series = design_series(
                values, np.cumsum(mu, dtype=np.int64), alpha * 0.5, SAMPLES
            )
            gamma, power = power_spectrum(log_x, series)
            score = scored(gamma, power)
            rows.append((name, q, alpha, gamma, score, peaks(gamma, score), values, mu))
        (an, aq, aa, ag, asc, ap, av, amu), (bn, bq, ba, bg, bsc, bp, bv, bmu) = rows
        width = max(float(ag[1]), float(bg[1]))
        shared = sum(
            1 for i in ap if bp.size and np.min(np.abs(bg[bp] - ag[i])) < width
        )
        cover = min(1.0, bp.size * 2.0 * width / (BAND[1] - BAND[0]))
        print(
            f"family {an:8s} alpha {aa:.6f} peaks {ap.size} against {bn:8s} alpha {ba:.6f}"
            f" peaks {bp.size} shared {shared} expect {cover * ap.size:.2f} bin {width:.4f}"
        )
        grid = np.linspace(BAND[0], BAND[1], 4096)
        rho, pv = spearmanr(np.interp(grid, ag, asc), np.interp(grid, bg, bsc))
        print(f"  spectra rho {float(rho):7.4f} p {float(pv):.3e}")
        for name, alpha, values, mu in ((an, aa, av, amu), (bn, ba, bv, bmu)):
            if values.size >= 5000000:
                continue
            support = mu != 0
            counts = []
            for _ in range(SIGNS):
                weight = np.zeros(values.size, dtype=np.int32)
                weight[support] = rng.choice([-1, 1], int(support.sum()))
                lx, sr = design_series(
                    values, np.cumsum(weight, dtype=np.int64), alpha * 0.5, SAMPLES
                )
                g, pw = power_spectrum(lx, sr)
                counts.append(int(peaks(g, scored(g, pw)).size))
            print(
                f"  random  {name:8s} peaks {min(counts)} to {max(counts)} over {SIGNS} sign draws"
            )
        for name, q, alpha, gamma, score in (
            (an, aq, aa, ag, asc),
            (bn, bq, ba, bg, bsc),
        ):
            conductor = DESIGNS[name][3]
            cells = []
            for label, targets in families(q, alpha, conductor, zeros, cache):
                _, value, null, p = shift_test(gamma, score, targets, rng)
                cells.append(f"{label} {value:6.3f}/{null:5.3f} p {p:.4f}")
            print(f"  {name:8s} " + "  ".join(cells))
        print()


def main():
    verb = sys.argv[1] if len(sys.argv) > 1 else "sieve"
    {"sieve": verb_sieve, "spectrum": verb_spectrum, "family": verb_family}[verb]()


if __name__ == "__main__":
    main()
