import math
import random
import sys
from itertools import combinations, product

import numpy as np

TOL = 1e-6
SEED = 20260930

EXHAUSTIVE = [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (4, 1), (4, 2)]
NEAR = [(2, 5), (2, 6), (3, 3), (3, 4), (4, 3)]

# WINDOWS


def all_windows(a, k):
    return list(product(range(a), repeat=k))


def from_code(code, wins):
    return frozenset(w for i, w in enumerate(wins) if (code >> i) & 1)


def letters(rule):
    return {c for w in rule for c in w}


def is_product(rule, a, k):
    return len(rule) == len(letters(rule)) ** k


# SPLICING


def straddles(u, v, rule, k):
    word = u + v
    return all(word[i : i + k] in rule for i in range(1, k))


def pair_count(rule, k):
    return sum(1 for u in rule for v in rule if straddles(u, v, rule, k))


def broken_pair(rule, k):
    for u in sorted(rule):
        for v in sorted(rule):
            if not straddles(u, v, rule, k):
                word = u + v
                bad = next(
                    word[i : i + k] for i in range(1, k) if word[i : i + k] not in rule
                )
                return u, v, bad
    return None


def closed(rule, k):
    return pair_count(rule, k) == len(rule) ** 2


# TRANSFER


def transfer(rule, a, k):
    if k == 1:
        return np.array([[float(len(rule))]])
    n = a ** (k - 1)
    m = np.zeros((n, n))
    for w in rule:
        s = 0
        for c in w[:-1]:
            s = s * a + c
        t = 0
        for c in w[1:]:
            t = t * a + c
        m[s, t] = 1.0
    return m


def rho_of(rule, a, k):
    return float(max(abs(np.linalg.eigvals(transfer(rule, a, k)))))


def counts_by_matrix(rule, a, k, n):
    m = transfer(rule, a, k)
    v = np.ones(m.shape[0])
    for _ in range(n - k + 1):
        v = m @ v
    return round(v.sum())


def counts_by_words(rule, a, k, n):
    total = 0
    for word in product(range(a), repeat=n):
        if all(word[i : i + k] in rule for i in range(n - k + 1)):
            total += 1
    return total


def kappa_of(rule, a, k):
    return math.log2(len(rule)) / k - math.log2(rho_of(rule, a, k))


def bound_of(rule, k):
    w = len(rule)
    return math.log2(w * w / pair_count(rule, k)) / (2 * k)


def gap_of(w, k):
    return math.log2(w * w / (w * w - 1)) / (2 * k)


# SWEEP


class Tally:
    def __init__(self):
        self.rules = 0
        self.live = 0
        self.products = 0
        self.mismatch = 0
        self.below_bound = 0
        self.below_gap = 0
        self.kappa_product = 0.0
        self.min_gap = None
        self.min_kappa = None
        self.min_bound = None

    def take(self, rule, a, k):
        if not rule:
            return
        self.rules += 1
        prod = is_product(rule, a, k)
        shut = closed(rule, k)
        self.products += prod
        if prod != shut:
            self.mismatch += 1
        rho = rho_of(rule, a, k)
        if rho < 0.5:
            return
        self.live += 1
        kappa = kappa_of(rule, a, k)
        if prod:
            self.kappa_product = max(self.kappa_product, abs(kappa))
            return
        bound = bound_of(rule, k)
        gap = gap_of(len(rule), k)
        if kappa < bound - TOL:
            self.below_bound += 1
        if kappa < gap - TOL:
            self.below_gap += 1
        if self.min_kappa is None or kappa < self.min_kappa:
            self.min_kappa, self.min_bound, self.min_gap = kappa, bound, gap

    def line(self, a, k, tag):
        low = "none" if self.min_kappa is None else f"{self.min_kappa:.6f}"
        return (
            f"{tag} a={a} k={k} rules={self.rules} products={self.products} "
            f"live={self.live} closed_vs_product_mismatch={self.mismatch} "
            f"kappa_below_pair_bound={self.below_bound} kappa_below_gap={self.below_gap} "
            f"max_abs_kappa_on_products={self.kappa_product:.2e} "
            f"min_kappa_on_live_nonproducts={low}"
        )


def sweep_all(a, k):
    wins = all_windows(a, k)
    tally = Tally()
    for code in range(1, 1 << len(wins)):
        tally.take(from_code(code, wins), a, k)
    return tally


def near_products(a, k):
    wins = all_windows(a, k)
    seen = set()
    for size in range(1, a + 1):
        for f in combinations(range(a), size):
            base = frozenset(w for w in wins if all(c in f for c in w))
            for w in wins:
                rule = base ^ {w}
                if rule and rule not in seen:
                    seen.add(rule)
                    yield rule


def sweep_near(a, k):
    tally = Tally()
    for rule in near_products(a, k):
        tally.take(rule, a, k)
    return tally


def sweep_random(a, k, density, count, rng):
    wins = all_windows(a, k)
    tally = Tally()
    for _ in range(count):
        rule = frozenset(w for w in wins if rng.random() < density)
        tally.take(rule, a, k)
    return tally


# VERBS


def verb_golden():
    a, k = 2, 2
    wins = all_windows(a, k)
    for name, rule in (
        ("golden, forbids 11", from_code(7, wins)),
        ("two constants, allows 00 and 11", frozenset({(0, 0), (1, 1)})),
    ):
        w = len(rule)
        n2 = pair_count(rule, k)
        bad = broken_pair(rule, k)
        by_words = counts_by_words(rule, a, k, 2 * k)
        by_matrix = counts_by_matrix(rule, a, k, 2 * k)
        rho = rho_of(rule, a, k)
        print(f"{name}: windows={sorted(rule)} w={w} w^2={w * w}")
        print(f"  N(2k) pairs={n2} words={by_words} matrix={by_matrix}")
        print(f"  first broken pair u={bad[0]} v={bad[1]} window outside W={bad[2]}")
        print(
            f"  rho={rho:.12f} rho^(2k)={rho ** (2 * k):.6f} <= N(2k)={n2} < w^2={w * w}"
        )
        print(
            f"  kappa={kappa_of(rule, a, k):.6f} "
            f"pair bound log2(w^2/N(2k))/(2k)={bound_of(rule, k):.6f} "
            f"gap bound log2(w^2/(w^2-1))/(2k)={gap_of(w, k):.6f}"
        )
        fekete = [(n, counts_by_words(rule, a, k, n), rho**n) for n in range(1, 13)]
        worst = min(c / r for _, c, r in fekete)
        print(f"  N(n) / rho^n over n=1..12, least {worst:.6f} (Fekete wants >= 1)")


def verb_generator():
    bad = 0
    total = 0
    for a, k in ((2, 2), (2, 3), (3, 2)):
        wins = all_windows(a, k)
        for code in range(1, 1 << len(wins)):
            rule = from_code(code, wins)
            for n in (k, k + 1, 2 * k):
                total += 1
                if counts_by_words(rule, a, k, n) != counts_by_matrix(rule, a, k, n):
                    bad += 1
            total += 1
            if pair_count(rule, k) != counts_by_words(rule, a, k, 2 * k):
                bad += 1
    print(
        f"counts by words, by the transfer matrix and by the pair test: {total} comparisons, {bad} disagreements"
    )


def verb_fekete():
    for a, k in ((2, 2), (2, 3), (3, 2)):
        wins = all_windows(a, k)
        live = 0
        worst = None
        broken = 0
        for code in range(1, 1 << len(wins)):
            rule = from_code(code, wins)
            rho = rho_of(rule, a, k)
            if rho < 0.5:
                continue
            live += 1
            for n in range(max(k - 1, 1), 15):
                ratio = counts_by_matrix(rule, a, k, n) / rho**n
                if ratio < 1 - 1e-9:
                    broken += 1
                if worst is None or ratio < worst:
                    worst = ratio
        print(
            f"a={a} k={k} live={live} N(n)/rho^n over n=k-1..14: least {worst:.6f}, below 1: {broken}"
        )


PRESETS = (
    ("golden", 2, 2, 7, 8),
    ("supergolden", 2, 3, 23, 8),
    ("plastic", 2, 3, 54, 8),
    ("tribonacci", 2, 3, 127, 8),
    ("no repeat", 4, 2, 31710, 4),
)


def verb_named():
    for name, a, k, code, levels in PRESETS:
        rule = from_code(code, all_windows(a, k))
        counts = [counts_by_words(rule, a, k, n) for n in range(1, levels + 1)]
        via = [
            counts_by_matrix(rule, a, k, n) for n in range(max(k - 1, 1), levels + 1)
        ]
        assert via == counts[len(counts) - len(via) :]
        print(
            f"{name}: alphabet={a} k={k} code={code} windows={len(rule)} "
            f"counts 1..{levels}={counts} rho={rho_of(rule, a, k):.9f} "
            f"kappa={kappa_of(rule, a, k):.6f}"
        )


def verb_splice():
    for a, k in EXHAUSTIVE:
        print(sweep_all(a, k).line(a, k, "exhaustive"))
    for a, k in NEAR:
        print(sweep_near(a, k).line(a, k, "near-product"))
    rng = random.Random(SEED)
    for a, k, density, count in (
        (2, 5, 0.5, 300),
        (2, 5, 0.9, 300),
        (2, 6, 0.97, 200),
        (3, 3, 0.9, 300),
        (4, 3, 0.95, 200),
    ):
        print(sweep_random(a, k, density, count, rng).line(a, k, f"random p={density}"))


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else ""
    if verb == "golden":
        verb_golden()
    elif verb == "generator":
        verb_generator()
    elif verb == "named":
        verb_named()
    elif verb == "fekete":
        verb_fekete()
    elif verb == "splice":
        verb_splice()
    else:
        verb_golden()
        verb_named()
        verb_generator()
        verb_fekete()
        verb_splice()
