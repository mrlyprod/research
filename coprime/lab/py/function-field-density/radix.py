import itertools
from math import log

import numpy as np

N = 12
SMALL = (2, 3, 5, 7, 11, 13)
SCHEDULES = [
    ("alternating base 2 {0,1} / base 3 {0,1}", ((2, (0, 1)), (3, (0, 1)))),
    ("alternating base 2 {0,1} / base 3 {0,2}", ((2, (0, 1)), (3, (0, 2)))),
    ("pure base 2 {0,1}", ((2, (0, 1)),)),
    ("pure base 3 {0,1}", ((3, (0, 1)),)),
    ("pure base 3 {0,2}", ((3, (0, 2)),)),
]


def build(schedule, n):
    vals = [0]
    place = 1
    for j in range(n):
        base, digits = schedule[j % len(schedule)]
        vals = [v + d * place for v in vals for d in digits]
        place *= base
    return sorted(vals)


def moebius_table(n):
    mu = np.ones(n + 1, dtype=np.int64)
    prime = np.ones(n + 1, dtype=bool)
    prime[:2] = False
    for p in range(2, int(n**0.5) + 1):
        if prime[p]:
            prime[p * p :: p] = False
    for p in np.nonzero(prime)[0]:
        mu[p::p] = -mu[p::p]
        mu[p * p :: p * p] = 0
    mu[0] = 0
    return mu


def present(vals, top):
    a = np.zeros(top + 1, dtype=np.int64)
    a[np.array(vals, dtype=np.int64)] = 1
    a[0] = 0
    return a


def coprime_pairs(vals, mu, top):
    a = present(vals, top)
    ordered = 0
    for m in range(1, top + 1):
        if mu[m]:
            c = int(a[m::m].sum())
            if c:
                ordered += int(mu[m]) * c * c
    unit = 1 if 1 in vals else 0
    diagonal = unit
    with_zero = unit
    return (ordered - diagonal) // 2 + with_zero


def divisible_fraction(vals, m):
    return sum(1 for v in vals if v % m == 0) / len(vals)


def joint_small(vals):
    arr = np.array(vals, dtype=np.int64)
    n = len(vals)
    total = 0
    for k in range(len(SMALL) + 1):
        for sub in itertools.combinations(SMALL, k):
            d = 1
            for p in sub:
                d *= p
            c = int((arr % d == 0).sum())
            total += (-1 if k % 2 else 1) * c * (c - 1) // 2
    return total / (n * (n - 1) // 2)


top = max(max(build(s, N)) for _, s in SCHEDULES)
mu = moebius_table(top)

print(
    "DOMAIN %d digit positions, unordered distinct pairs including the zero point" % N
)
print(
    "points per schedule %d, pairs per schedule %d"
    % (1 << N, (1 << N) * ((1 << N) - 1) // 2)
)

rows = {}
for name, sched in SCHEDULES:
    vals = build(sched, N)
    pairs = len(vals) * (len(vals) - 1) // 2
    c = coprime_pairs(vals, mu, max(vals))
    rows[name] = (vals, c, c / pairs)
    print("%-40s coprime %7d  density %.6f" % (name, c, c / pairs))

a, b = SCHEDULES[0][0], SCHEDULES[1][0]
print("alternating gap = %.6f" % (rows[b][2] - rows[a][2]))
print(
    "pure base 2 against 1/zeta(2) = 0.607927, difference %+.6f"
    % (rows[SCHEDULES[2][0]][2] - 0.6079271018540267)
)

vals = rows[a][0]
res = [sum(1 for v in vals if v % 6 == r) for r in range(6)]
print(
    "residues mod 6 of the first alternating schedule: %s"
    % (", ".join(str(x) for x in res))
)
for name in (a, b):
    vs = rows[name][0]
    print(
        "%-40s even fraction %.6f  divisible by 3 fraction %.6f"
        % (name, divisible_fraction(vs, 2), divisible_fraction(vs, 3))
    )

for name in (a, b):
    vs = rows[name][0]
    prod = 1.0
    for p in SMALL:
        prod *= 1 - divisible_fraction(vs, p) ** 2
    exact = joint_small(vs)
    print(
        "%-40s Euler product through 13 = %.6f  exact joint = %.6f  product minus exact = %+.6f"
        % (name, prod, exact, prod - exact)
    )

for p in SMALL:
    fa = divisible_fraction(rows[a][0], p)
    fb = divisible_fraction(rows[b][0], p)
    ra = 1 - fa * fa
    rb = 1 - fb * fb
    print(
        "prime %2d  factor %.6f -> %.6f  log advantage %+.6f"
        % (p, ra, rb, log(rb / ra))
    )
