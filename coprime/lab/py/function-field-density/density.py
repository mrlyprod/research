import itertools
from fractions import Fraction
from math import log

import numpy as np

Q = 3
NMAX = 14
FACTOR_DEG = 6
EULER_N = 10
EULER_DEG = 5
LEVELS = (10, 12, 14)


def trim(c):
    while c and c[-1] == 0:
        del c[-1]
    return c


def divmod_poly(a, b):
    r = list(a)
    db = len(b) - 1
    inv = pow(b[-1], Q - 2, Q)
    q = [0] * max(0, len(r) - db)
    while r and len(r) - 1 >= db:
        d = len(r) - 1 - db
        f = (r[-1] * inv) % Q
        q[d] = f
        for i, c in enumerate(b):
            r[i + d] = (r[i + d] - f * c) % Q
        trim(r)
    return tuple(trim(q)), tuple(r)


def monics(d):
    for tail in itertools.product(range(Q), repeat=d):
        yield tail + (1,)


def irreducibles(maxdeg):
    out = []
    for d in range(1, maxdeg + 1):
        for p in monics(d):
            if all(divmod_poly(p, s)[1] != () for s in out if len(s) - 1 <= d // 2):
                out.append(p)
    return out


def label(p):
    parts = []
    for k in range(len(p) - 1, -1, -1):
        c = p[k]
        if c == 0:
            continue
        if k == 0:
            parts.append(str(c))
        else:
            base = "t" if k == 1 else "t^%d" % k
            parts.append(base if c == 1 else "%d%s" % (c, base))
    return " + ".join(parts) if parts else "0"


def poly_of(i, n):
    return tuple(trim([(i >> k) & 1 for k in range(n)]))


def masks(primes, n):
    bits = ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(np.int64)
    out = {}
    for p in primes:
        dp = len(p) - 1
        rows = np.zeros((n, dp), dtype=np.int64)
        for k in range(n):
            r = divmod_poly(tuple([0] * k + [1]), p)[1]
            rows[k, : len(r)] = r
        out[p] = ~((bits @ rows) % Q).any(axis=1)
    return out


def factor_sets(primes, msk, n):
    small = [[] for _ in range(1 << n)]
    for p in primes:
        for i in np.nonzero(msk[p])[0]:
            if i:
                small[int(i)].append(p)
    ids = {p: k for k, p in enumerate(primes)}
    degs = [len(p) - 1 for p in primes]
    seen = {}
    sets = [()] * (1 << n)
    for i in range(1, 1 << n):
        h = poly_of(i, n)
        keys = []
        for p in small[i]:
            keys.append(ids[p])
            while True:
                q, r = divmod_poly(h, p)
                if r != ():
                    break
                h = q
        if len(h) > 1:
            if h not in seen:
                seen[h] = len(degs)
                degs.append(len(h) - 1)
            keys.append(seen[h])
        sets[i] = tuple(sorted(keys))
    return sets, degs


def sieve_count(sets, degs, n, degcap=None):
    cnt = {}
    for i in range(1, 1 << n):
        ps = sets[i]
        if degcap is not None:
            ps = tuple(k for k in ps if degs[k] <= degcap)
        for k in range(len(ps) + 1):
            for sub in itertools.combinations(ps, k):
                cnt[sub] = cnt.get(sub, 0) + 1
    total = 0
    for key, m in cnt.items():
        total += (-1 if len(key) % 2 else 1) * (m * m + 2 * m)
    return total


def brute(n):
    polys = [poly_of(i, n) for i in range(1 << n)]
    c = 0
    for a in polys:
        for b in polys:
            x, y = a, b
            while y != ():
                x, y = y, divmod_poly(x, y)[1]
            c += 1 if len(x) == 1 else 0
    return c


primes = irreducibles(FACTOR_DEG)
msk = masks(primes, NMAX)
sets, degs = factor_sets(primes, msk, NMAX)

print(
    "DOMAIN q = 3, S = {0, 1}, degree below n, n up to %d, ordered pairs over all %d^2 including the zero polynomial"
    % (NMAX, 1 << NMAX)
)
print("irreducibles of degree at most %d: %d" % (FACTOR_DEG, len(primes)))

for n in range(1, 5):
    assert sieve_count(sets, degs, n) == brute(n), n
print("cross-check against Euclid on every ordered pair, n = 1 to 4: passed")

for n in LEVELS:
    z = sieve_count(sets, degs, n)
    print(
        "n = %2d  Z = %d / %d  density = %.6f  gap from 9/16 = %+.6f"
        % (n, z, 4**n, z / 4**n, z / 4**n - 9 / 16)
    )

print("prediction (2/3) * (3/4) / (8/9) = 9/16 = %.6f" % (9 / 16))

pi = {
    p: Fraction(int(msk[p][: 1 << EULER_N].sum()), 1 << EULER_N)
    for p in primes
    if len(p) - 1 <= EULER_DEG
}
for p in primes:
    if len(p) - 1 == 1:
        print("pi(%s) at n = %d = %s = %.6f" % (label(p), EULER_N, pi[p], float(pi[p])))
for d in (2, 3):
    ps = [p for p in primes if len(p) - 1 == d]
    devs = [abs(pi[p] - Fraction(1, Q**d)) for p in ps]
    print(
        "degree %d, %d primes, mean |pi - 3^-%d| = %.6f, max = %.6f"
        % (d, len(ps), d, float(sum(devs) / len(devs)), float(max(devs)))
    )

prod = Fraction(1)
for d in range(1, EULER_DEG + 1):
    for p in primes:
        if len(p) - 1 == d:
            prod *= 1 - pi[p] ** 2
    print("marginal Euler product through degree %d = %.6f" % (d, float(prod)))

shared = sieve_count(sets, degs, EULER_N, degcap=EULER_DEG)
exact = Fraction(shared, 4**EULER_N)
print(
    "exact no shared prime of degree at most %d at n = %d = %d/%d = %.6f"
    % (EULER_DEG, EULER_N, shared, 4**EULER_N, float(exact))
)
print("marginal product minus exact = %+.6f" % float(prod - exact))

g = log(4) / log(3)
print(
    "gamma = log_3(4) = %.6f  gamma/2 = %.6f  window (gamma/2, 1/2] empty = %s"
    % (g, g / 2, g / 2 > 0.5)
)
