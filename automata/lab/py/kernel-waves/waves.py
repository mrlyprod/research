import itertools
import sys
from collections import Counter
from fractions import Fraction
from math import ceil, floor

import numpy as np

N = 256
SEED = 20260
SEEDS = (20260, 4093)
STEPS = 64
DENSITY = 0.5
KMAX = N // 2
DEAD = 0.02
FULL = 0.98
GAIN = 3.0
TOL = 1e-9

BUGS_B = (Fraction(34, 120), Fraction(45, 120))
BUGS_S = (Fraction(34, 120), Fraction(58, 120))
RULES = [
    ("bugs", BUGS_B, BUGS_S),
    ("wide-survive", BUGS_B, (Fraction(28, 100), Fraction(60, 100))),
    ("narrow-birth", (Fraction(30, 100), Fraction(34, 100)), BUGS_S),
]
GRID_B_LO = [Fraction(20, 100), Fraction(30, 100), Fraction(35, 100)]
GRID_B_HI = [
    Fraction(35, 100),
    Fraction(40, 100),
    Fraction(45, 100),
    Fraction(50, 100),
    Fraction(55, 100),
]
GRID_S_LO = [Fraction(0), Fraction(30, 100)]
GRID_S_HI = [
    Fraction(45, 100),
    Fraction(50, 100),
    Fraction(55, 100),
    Fraction(60, 100),
    Fraction(70, 100),
    Fraction(80, 100),
]


def residue_corners(d, base=2):
    return [
        [(i // base ** (d - 1 - j)) % base for j in range(d)] for i in range(base**d)
    ]


def tile(code, side=3):
    corners = residue_corners(2)
    filled = [(code >> i) & 1 for i in range(len(corners))]
    out = np.zeros((side, side), dtype=np.uint8)
    for x in range(side):
        for y in range(side):
            out[x, y] = filled[(x % 2) * 2 + (y % 2)]
    return out


def kron_power(t, level):
    out = t
    for _ in range(1, level):
        out = np.kron(out, t)
    return out


def design_mask(code, side, level):
    mask = kron_power(tile(code, side), level).copy()
    c = (side**level - 1) // 2
    mask[c, c] = 0
    return mask


def box_mask(r):
    mask = np.ones((2 * r + 1, 2 * r + 1), dtype=np.uint8)
    mask[r, r] = 0
    return mask


def torus_kernel(mask):
    side = mask.shape[0]
    c = (side - 1) // 2
    kernel = np.zeros((N, N))
    xs, ys = np.nonzero(mask)
    kernel[(xs - c) % N, (ys - c) % N] = 1
    return kernel


FX = np.fft.fftfreq(N) * N
KX, KY = np.meshgrid(FX, FX, indexing="ij")
RING = np.rint(np.sqrt(KX**2 + KY**2)).astype(int)
COUNT = np.bincount(RING.ravel())


def ring_sum(values):
    return np.bincount(RING.ravel(), weights=values.ravel(), minlength=COUNT.size)


def ring_mean(values):
    return ring_sum(values) / COUNT


def first_min(a, start):
    return next((k for k in range(start, KMAX) if a[k + 1] >= a[k]), None)


def first_max(a, start):
    return next((k for k in range(start, KMAX) if a[k + 1] < a[k]), None)


def mask_profile(full):
    a = ring_mean(np.abs(full))[: KMAX + 1]
    k_min = first_min(a, 1)
    k_2 = first_max(a, k_min + 1)
    k_min2 = first_min(a, k_2 + 1)
    half = a[k_2] / 2
    lobe = [k for k in range(k_min, k_min2 + 1) if a[k] >= half]
    g = ring_mean(np.real(full))[: KMAX + 1]
    k_neg = 1 + int(np.argmin(g[1:]))
    lo = k_neg
    while lo > 1 and g[lo - 1] < 0:
        lo -= 1
    hi = k_neg
    while hi < KMAX and g[hi + 1] < 0:
        hi += 1
    return a, k_min, k_2, k_min2, (lobe[0], lobe[-1]), g, k_neg, (lo, hi)


def riesz_product(code, side, level):
    t = tile(code, side)
    c = (side - 1) // 2
    xs, ys = np.nonzero(t)
    ex, ey = xs - c, ys - c
    out = np.ones((N, N), dtype=complex)
    for j in range(level):
        tx = side**j * KX / N
        ty = side**j * KY / N
        p = np.zeros((N, N), dtype=complex)
        for a, b in zip(ex, ey):
            p += np.exp(-2j * np.pi * (tx * a + ty * b))
        out *= p
    return out - int(t[c, c])


def thresholds(m, lo, hi):
    return ceil(lo * m), floor(hi * m)


def soup(seed):
    rng = np.random.default_rng(seed)
    return (rng.random((N, N)) < DENSITY).astype(np.uint8)


def evolve(half, m, birth, survive, seed=SEED):
    blo, bhi = thresholds(m, *birth)
    slo, shi = thresholds(m, *survive)
    grid = soup(seed)
    churn = None
    for _ in range(STEPS):
        count = np.rint(np.fft.irfft2(np.fft.rfft2(grid) * half, s=(N, N))).astype(int)
        born = (grid == 0) & (count >= blo) & (count <= bhi)
        kept = (grid == 1) & (count >= slo) & (count <= shi)
        new = (born | kept).astype(np.uint8)
        churn = int(np.count_nonzero(new != grid))
        grid = new
    return grid, churn


def spectrum(grid):
    f = np.fft.rfft2(grid - grid.mean())
    power = np.zeros((N, N))
    power[:, : N // 2 + 1] = np.abs(f) ** 2
    power[:, N // 2 + 1 :] = (
        np.abs(f[np.r_[0, N - 1 : 0 : -1], N // 2 - 1 : 0 : -1]) ** 2
    )
    s = ring_sum(power)[: KMAX + 1]
    d = s / COUNT[: KMAX + 1]
    if s[1:].sum() == 0:
        return 1, 0.0, 0.0
    k = 1 + int(np.argmax(d[1:]))
    share = s[k] / s[1:].sum()
    gain = d[k] / d[1:].mean()
    return k, share, gain


def classify(density, churn, k, gain, side):
    if density < DEAD:
        return "dead"
    if density > FULL:
        return "full"
    if churn > 0:
        return "active"
    if gain < GAIN:
        return "flat"
    if N / k > 2 * side:
        return "coarse"
    return "ring"


def fmt(fr):
    return f"{float(fr):.2f}"


def rule_label(birth, survive):
    return f"B[{fmt(birth[0])},{fmt(birth[1])}] S[{fmt(survive[0])},{fmt(survive[1])}]"


def main():
    t7 = tile(7)
    assert t7.sum() == 8 and t7[1, 1] == 0 and t7.shape == (3, 3)
    m7 = design_mask(7, 3, 2)
    assert m7.shape == (9, 9) and m7.sum() == 64
    assert design_mask(6, 3, 1).sum() == 4 and design_mask(6, 3, 1)[0, 1] == 1
    assert design_mask(9, 3, 1).sum() == 4 and design_mask(9, 3, 1)[0, 0] == 1
    print(
        f"domain: torus {N}x{N}, soup density {DENSITY}, seed {SEED}, {STEPS} steps, rings k = 1..{KMAX}"
    )
    print(
        f"classes: dead is density < {DEAD}, full is density > {FULL}, active is churn > 0 at the last step, flat is gain < {GAIN}, coarse is N/k* > 2 side, ring is the rest"
    )
    for seed in SEEDS:
        k0, share0, gain0 = spectrum(soup(seed))
        print(
            f"soup seed {seed} at step 0: k* {k0}, share {share0:.4f}, gain {gain0:.2f}"
        )
    print()
    print(
        "table 1: the Riesz product, max |DFT(mask) - (prod_j P(side^j xi / N) - tile centre)| over all frequencies"
    )
    print(
        f"{'code':>4} {'level':>5} {'side':>4} {'ones':>5} {'centre':>6} {'max err':>10} {'tag':>8}"
    )
    riesz_ok = True
    for code in (7, 6, 9):
        for level in (2, 3, 4):
            mask = design_mask(code, 3, level)
            err = np.max(
                np.abs(np.fft.fft2(torus_kernel(mask)) - riesz_product(code, 3, level))
            )
            ok = err < TOL
            riesz_ok &= ok
            print(
                f"{code:>4} {level:>5} {3**level:>4} {mask.sum():>5} {int(tile(code)[1, 1]):>6} {err:>10.2e} {'Verified' if ok else 'Refuted':>8}"
            )
    print(
        f"Riesz product at codes 7, 6, 9, levels 2..4, tolerance {TOL:.0e}: {'Verified' if riesz_ok else 'Refuted'}"
    )
    print()
    masks = [
        ("box r=5", box_mask(5), 11),
        ("box r=13", box_mask(13), 27),
        ("code 7 L2", design_mask(7, 3, 2), 9),
        ("code 7 L3", design_mask(7, 3, 3), 27),
        ("code 7 L4", design_mask(7, 3, 4), 81),
        ("code 6 L3", design_mask(6, 3, 3), 27),
        ("code 9 L3", design_mask(9, 3, 3), 27),
    ]
    print(
        "mask table: ones m, first minimum k_min, first secondary maximum k_2, second minimum, half-height lobe of the ring-mean |DFT|, then the most negative ring k_neg and the negative band of the ring-mean signed DFT"
    )
    print(
        f"{'mask':>10} {'side':>4} {'m':>5} {'k_min':>5} {'k_2':>5} {'k_min2':>6} {'lobe':>9} {'N/k_min':>8} {'N/k_2':>7} {'N/k_2/side':>10} {'A(k_min)/m':>10} {'A(k_2)/m':>9} {'k_neg':>5} {'neg band':>9} {'g(k_neg)/m':>10}"
    )
    profiles = {}
    for name, mask, side in masks:
        kernel = torus_kernel(mask)
        m = int(mask.sum())
        a, k_min, k_2, k_min2, lobe, g, k_neg, neg = mask_profile(np.fft.fft2(kernel))
        profiles[name] = (
            np.fft.rfft2(kernel),
            m,
            k_min,
            k_2,
            lobe,
            side,
            g,
            k_neg,
            neg,
        )
        print(
            f"{name:>10} {side:>4} {m:>5} {k_min:>5} {k_2:>5} {k_min2:>6} {f'{lobe[0]}..{lobe[1]}':>9} {N / k_min:>8.2f} {N / k_2:>7.2f} {N / k_2 / side:>10.3f} {a[k_min] / m:>10.4f} {a[k_2] / m:>9.4f} {k_neg:>5} {f'{neg[0]}..{neg[1]}':>9} {g[k_neg] / m:>10.4f}"
        )
    print()
    results = {}
    print(
        "run table: the state at the last step for every mask under the three named rules"
    )
    print(
        f"{'mask':>10} {'rule':>13} {'B':>9} {'S':>9} {'density':>7} {'churn':>6} {'k*':>3} {'N/k*':>7} {'share':>6} {'gain':>6} {'class':>6}"
    )
    for name, mask, side in masks:
        half, m, k_min, k_2, lobe, side, g, k_neg, neg = profiles[name]
        for rule, birth, survive in RULES:
            grid, churn = evolve(half, m, birth, survive)
            density = grid.mean()
            k, share, gain = spectrum(grid)
            cls = classify(density, churn, k, gain, side)
            results[(name, rule)] = (k, share, gain, density, churn, cls)
            b = "{}-{}".format(*thresholds(m, *birth))
            s = "{}-{}".format(*thresholds(m, *survive))
            print(
                f"{name:>10} {rule:>13} {b:>9} {s:>9} {density:>7.4f} {churn:>6} {k:>3} {N / k:>7.2f} {share:>6.3f} {gain:>6.1f} {cls:>6}"
            )
    print()
    print("table 2: same mask across the three named rules, k* of every ring still")
    print(f"{'mask':>10} {'bugs':>5} {'wide':>5} {'narrow':>6} {'rings':>5}")
    rings_named = 0
    for name, mask, side in masks:
        ks = [
            results[(name, rule)][0] if results[(name, rule)][5] == "ring" else None
            for rule, _, _ in RULES
        ]
        cells = " ".join(f"{k:>5}" if k is not None else f"{'-':>5}" for k in ks)
        kept = [k for k in ks if k is not None]
        rings_named += len(kept)
        print(f"{name:>10} {cells} {len(kept):>5}")
    print(
        f"ring stills under the three named rules: {rings_named} of {len(masks) * len(RULES)}; the wavelength laws are vacuous on the named rules, tested on the rule grid below"
    )
    print()
    print(
        "table 3: same named rule across the mask family, class of every cell against the mask's k_min and k_2"
    )
    print(
        f"{'rule':>13} {'mask':>10} {'side':>4} {'class':>6} {'k*':>3} {'k_min':>5} {'k_2':>5}"
    )
    for rule, _, _ in RULES:
        for name, mask, side in masks:
            k, share, gain, density, churn, cls = results[(name, rule)]
            half, m, k_min, k_2, lobe, side, g, k_neg, neg = profiles[name]
            print(
                f"{rule:>13} {name:>10} {side:>4} {cls:>6} {k if cls in ('ring', 'coarse', 'flat') else '-':>3} {k_min:>5} {k_2:>5}"
            )
    print()
    grid_rules = [
        (b_lo, b_hi, s_lo, s_hi)
        for b_lo, b_hi, s_lo, s_hi in itertools.product(
            GRID_B_LO, GRID_B_HI, GRID_S_LO, GRID_S_HI
        )
    ]
    print(
        f"rule grid: birth lower {[fmt(x) for x in GRID_B_LO]}, birth upper {[fmt(x) for x in GRID_B_HI]}, survive lower {[fmt(x) for x in GRID_S_LO]}, survive upper {[fmt(x) for x in GRID_S_HI]}, {len(grid_rules)} rules, seeds {list(SEEDS)}, {len(grid_rules) * len(masks) * len(SEEDS)} runs"
    )
    print()
    print(
        "census: class counts per mask over the rule grid, then the ring stills against the mask"
    )
    print(
        f"{'mask':>10} {'neg rings':>9} {'dead':>4} {'full':>4} {'active':>6} {'flat':>4} {'coarse':>6} {'ring':>4} {'k* range':>9} {'k*/k_2 range':>13} {'in lobe':>7} {'nearer k_2':>10} {'|k*-k_2|<=1':>11} {'|k*-k_min|<=1':>13} {'g(k*)<0':>7} {'|k*-k_neg|<=1':>13}"
    )
    census = {}
    coarse = {}
    totals = Counter()
    for name, mask, side in masks:
        half, m, k_min, k_2, lobe, side, g, k_neg, neg = profiles[name]
        tally = Counter()
        rings = []
        coarse[name] = Counter()
        for seed in SEEDS:
            for birth_lo, birth_hi, survive_lo, survive_hi in grid_rules:
                birth, survive = (birth_lo, birth_hi), (survive_lo, survive_hi)
                grid, churn = evolve(half, m, birth, survive, seed)
                density = grid.mean()
                k, share, gain = spectrum(grid)
                cls = classify(density, churn, k, gain, side)
                tally[cls] += 1
                if cls == "ring":
                    rings.append(
                        (rule_label(birth, survive), k, share, gain, density, seed)
                    )
                if cls == "coarse":
                    coarse[name][k] += 1
        census[name] = rings
        totals.update(tally)
        neg_rings = int(np.sum(g[1:] < 0))
        ks = [r[1] for r in rings]
        in_lobe = sum(lobe[0] <= k <= lobe[1] for k in ks)
        nearer = sum(abs(k - k_2) < abs(k - k_min) for k in ks)
        at_2 = sum(abs(k - k_2) <= 1 for k in ks)
        at_min = sum(abs(k - k_min) <= 1 for k in ks)
        negative = sum(g[k] < 0 for k in ks)
        at_neg = sum(abs(k - k_neg) <= 1 for k in ks)
        k_range = f"{min(ks)}..{max(ks)}" if ks else "-"
        ratio = f"{min(ks) / k_2:.2f}..{max(ks) / k_2:.2f}" if ks else "-"
        print(
            f"{name:>10} {neg_rings:>9} {tally['dead']:>4} {tally['full']:>4} {tally['active']:>6} {tally['flat']:>4} {tally['coarse']:>6} {tally['ring']:>4} {k_range:>9} {ratio:>13} {in_lobe:>7} {nearer:>10} {at_2:>11} {at_min:>13} {negative:>7} {at_neg:>13}"
        )
    print(
        f"{'total':>10} {'-':>9} {totals['dead']:>4} {totals['full']:>4} {totals['active']:>6} {totals['flat']:>4} {totals['coarse']:>6} {totals['ring']:>4}"
    )
    print()
    print(
        "ring stills per mask: k* histogram over the grid rules and seeds, the coarse stills' k* histogram, and the strongest ring still"
    )
    for name, mask, side in masks:
        rings = census[name]
        hist = Counter(r[1] for r in rings)
        line = " ".join(f"{k}x{hist[k]}" for k in sorted(hist))
        coarse_line = " ".join(f"{k}x{coarse[name][k]}" for k in sorted(coarse[name]))
        best = max(rings, key=lambda r: r[3]) if rings else None
        best_line = (
            f"{best[0]} seed {best[5]} k* {best[1]} N/k* {N / best[1]:.2f} share {best[2]:.3f} gain {best[3]:.1f} density {best[4]:.3f}"
            if best
            else "none"
        )
        print(
            f"{name:>10}: ring {line if line else '-'}; coarse {coarse_line if coarse_line else '-'}; strongest {best_line}"
        )
    print()
    print("every ring still: mask, rule, seed, k*, N/k*, g(k*)/m, share, gain, density")
    for name, mask, side in masks:
        half, m, k_min, k_2, lobe, side, g, k_neg, neg = profiles[name]
        for label, k, share, gain, density, seed in census[name]:
            print(
                f"{name:>10} {label} {seed:>5} {k:>3} {N / k:>6.2f} {g[k] / m:>8.4f} {share:>6.3f} {gain:>6.1f} {density:>6.3f}"
            )
    print()
    laws = [
        (
            "k* of every ring still lies in the mask's half-height first side lobe",
            lambda k, p: p[4][0] <= k <= p[4][1],
        ),
        (
            "k* of every ring still is nearer the mask's k_2 than its k_min",
            lambda k, p: abs(k - p[3]) < abs(k - p[2]),
        ),
        (
            "k* of every ring still sits within one ring of the mask's k_2",
            lambda k, p: abs(k - p[3]) <= 1,
        ),
        (
            "k* of every ring still sits within one ring of the mask's k_min",
            lambda k, p: abs(k - p[2]) <= 1,
        ),
        (
            "the ring-mean signed DFT of the mask is negative at k* of every ring still",
            lambda k, p: p[6][k] < 0,
        ),
        (
            "k* of every ring still sits within one ring of the mask's k_neg",
            lambda k, p: abs(k - p[7]) <= 1,
        ),
        (
            "k* of every ring still lies in the mask's negative band",
            lambda k, p: p[8][0] <= k <= p[8][1],
        ),
    ]
    total = sum(len(v) for v in census.values())
    for text, test in laws:
        failures = [
            (name, r[0], r[1])
            for name in census
            for r in census[name]
            if not test(r[1], profiles[name])
        ]
        verdict = "Verified" if total and not failures else "Refuted"
        witness = (
            f"; first witness {failures[0][0]} {failures[0][1]} k* {failures[0][2]}, {len(failures)} of {total} fail"
            if failures
            else f"; {total} ring stills"
        )
        print(f"law: {text}: {verdict}{witness}")
    print()
    print(
        "per-mask verdicts of the three sharpest laws over the ring stills, with the failing k* values"
    )
    print(
        f"{'mask':>10} {'rings':>5} {'in negative band':>28} {'nearer k_2 than k_min':>28} {'within one ring of k_min':>28}"
    )
    per_mask = [
        ("negative band", lambda k, p: p[8][0] <= k <= p[8][1]),
        ("nearer k_2", lambda k, p: abs(k - p[3]) < abs(k - p[2])),
        ("at k_min", lambda k, p: abs(k - p[2]) <= 1),
    ]
    verified_masks = []
    for name, mask, side in masks:
        cells = []
        for text, test in per_mask:
            bad = sorted({r[1] for r in census[name] if not test(r[1], profiles[name])})
            verdict = "Verified" if census[name] and not bad else "Refuted"
            if text == "negative band" and verdict == "Verified":
                verified_masks.append(name)
            cells.append(
                f"{verdict + (' k* ' + ' '.join(map(str, bad)) if bad else ''):>28}"
            )
        print(f"{name:>10} {len(census[name]):>5} " + " ".join(cells))
    print(
        f"negative band Verified on {len(verified_masks)} masks ({', '.join(verified_masks)}) holding {sum(len(census[n]) for n in verified_masks)} ring stills"
    )
    spreads = []
    for name, rows in census.items():
        ks = sorted({r[1] for r in rows})
        if len(ks) >= 2:
            spreads.append((name, ks[0], ks[-1]))
    fixed = "Verified" if total and not spreads else "Refuted"
    witness = (
        f"; first witness {spreads[0][0]} k* {spreads[0][1]}..{spreads[0][2]}"
        if spreads
        else ""
    )
    print(
        f"law: every ring still on one mask shares one k* across the grid rules: {fixed}{witness}"
    )
    print()
    print(
        "scaling: carpet family per level, k_2 and the median N/k* over side of its ring stills"
    )
    for name in ("code 7 L2", "code 7 L3", "code 7 L4"):
        half, m, k_min, k_2, lobe, side, g, k_neg, neg = profiles[name]
        ks = [r[1] for r in census[name]]
        med = f"{np.median([N / k / side for k in ks]):.3f}" if ks else "-"
        print(
            f"{name:>10} side {side:>3} k_2 {k_2:>3} N/k_2/side {N / k_2 / side:.3f} ring stills {len(ks):>3} median N/k*/side {med}"
        )
    if not riesz_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
