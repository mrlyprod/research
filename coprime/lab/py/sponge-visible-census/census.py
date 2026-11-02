import sys
import time
from fractions import Fraction
from itertools import product

import numpy as np
from mpmath import mp, zeta

Q = 3
D = 3
PARITY = {(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)}
DIGITS = [v for v in product(range(Q), repeat=D) if tuple(c % 2 for c in v) in PARITY]
K = len(DIGITS)
WALK = 12
BASE_LOCAL = Fraction(sum(1 for v in DIGITS if any(v)), K)
BRACKET = BASE_LOCAL * Fraction(Q**D, Q**D - 1)
LEGAL = np.zeros((Q,) * D, dtype=bool)
for _v in DIGITS:
    LEGAL[_v] = True


def mobius(limit):
    mu = np.ones(limit + 1, dtype=np.int64)
    mu[0] = 0
    composite = np.zeros(limit + 1, dtype=bool)
    primes = []
    for i in range(2, limit + 1):
        if not composite[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > limit:
                break
            composite[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = mu[i] * mu[p]
    return mu


def divisible_count(d, n):
    counts = np.zeros((d, d, d), dtype=np.int64)
    counts[0, 0, 0] = 1
    place = 1 % d
    for _ in range(n):
        nxt = np.zeros_like(counts)
        for v in DIGITS:
            shift = (place * v[0] % d, place * v[1] % d, place * v[2] % d)
            nxt += np.roll(counts, shift, axis=(0, 1, 2))
        counts = nxt
        place = place * Q % d
    return int(counts[0, 0, 0])


def in_design(pts, n):
    ok = np.ones(len(pts), dtype=bool)
    rest = pts
    for _ in range(n):
        dig = rest % Q
        ok &= LEGAL[dig[:, 0], dig[:, 1], dig[:, 2]]
        rest = rest // Q
    return int(np.count_nonzero(ok))


def primitive_box(side):
    axis = np.arange(side + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
    pts = grid.reshape(-1, D)
    return pts[np.gcd.reduce(pts, axis=1) == 1]


def visible_count(n, cutoff):
    span = Q**n
    mu = mobius(cutoff)
    head = 0
    for d in range(1, cutoff + 1):
        if mu[d]:
            head += int(mu[d]) * (divisible_count(d, n) - 1)
    partial = np.zeros(span, dtype=np.int64)
    for d in range(1, cutoff + 1):
        if mu[d]:
            partial[d::d] += mu[d]
    tail = 0
    for g in range(cutoff + 1, span):
        side = (span - 1) // g
        weight = int(partial[g])
        if weight == 0:
            continue
        tail += weight * in_design(g * primitive_box(side), n)
    return head - tail


def brute_count(n):
    digits = np.array(DIGITS, dtype=np.int64)
    total = 0
    for lead in digits:
        pts = lead.reshape(1, D)
        for _ in range(n - 1):
            pts = (pts[:, None, :] * Q + digits[None, :, :]).reshape(-1, D)
        total += int(np.count_nonzero(np.gcd.reduce(pts, axis=1) == 1))
    return total


def ladder(n):
    return max(8, min(round(Q ** (n / 2)), 120))


def main():
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    mp.dps = 30
    delta = mp.mpf(BRACKET.numerator) / BRACKET.denominator / zeta(D)
    print(f"design q={Q} D={D} k={K} bracket={BRACKET} delta={mp.nstr(delta, 10)}")
    print(f"domain n=1..{top} coordinates 0..{Q**top - 1}")
    for n in range(1, top + 1):
        cut = ladder(n)
        clock = time.time()
        a = visible_count(n, cut)
        line = f"A({n}) = {a}  G={cut}  {time.time() - clock:.1f}s"
        if n <= 6:
            line += f"  brute={brute_count(n)}"
        if n >= 7:
            gap = (delta * mp.mpf(K) ** n - a) / mp.mpf(WALK) ** n
            line += f"  gap/{WALK}^n = {mp.nstr(gap, 5)}"
        print(line, flush=True)
    if top >= 9:
        low = visible_count(9, 100)
        high = visible_count(9, 150)
        print(f"cutoff check A(9) G=100 {low} G=150 {high} agree={low == high}")


main()
