import sys
from itertools import pairwise
from math import log


def tail(x, d, s, t):
    lo, hi = 1, 2
    while x + hi * (d - x) + hi * (hi - 1) * s <= t:
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        m = (lo + hi) // 2
        if x + m * (d - x) + m * (m - 1) * s <= t:
            lo = m
        else:
            hi = m
    return lo - 1


def count(t):
    r = (-1, 2, 2, 3)
    n = 0
    work = []
    for i in range(4):
        a, b, c = r[:i] + r[i + 1 :]
        k = 2 * (a + b + c) - r[i]
        if k <= t:
            n += 1
            work += [(b, c, a, k), (a, c, b, k), (a, b, c, k)]
    while work:
        y, z, x, d = work.pop()
        s = y + z
        while True:
            e = 2 * (s + d) - x
            if e > t:
                break
            n += 1
            p = 2 * (y + d + e) - z
            q = 2 * (z + d + e) - y
            if p > t and q > t:
                n += tail(d, e, s, t)
                break
            if p <= t:
                work.append((y, d, z, e))
            if q <= t:
                work.append((z, d, y, e))
            x, d = d, e
    return n


if __name__ == "__main__":
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    row = [count(10**j) for j in range(1, top + 1)]
    print("root excluded", ",".join(map(str, row)))
    every = [1] + [v + 4 for v in row]
    print("every circle", ",".join(map(str, every)))
    print(
        "exponents",
        ",".join(f"{log(b / a) / log(10):.4f}" for a, b in pairwise(row)),
    )
    print(
        "exponents every",
        ",".join(f"{log(b / a) / log(10):.4f}" for a, b in pairwise(every)),
    )
    delta = 1.3056867280498771846
    print(
        "ratios", ",".join(f"{v / 10 ** (j * delta):.5f}" for j, v in enumerate(every))
    )
    print("cut mod 24", ",".join(str(10**j % 24) for j in range(1, top + 1)))
