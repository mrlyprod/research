import sys
from math import gcd


def new_terms(n):
    top = 3 * n * n
    norms = set()
    side = int((4 * top / 3) ** 0.5) + 2
    for a in range(-side, side + 1):
        for b in range(-side, side + 1):
            k = a * a - a * b + b * b
            if 0 < k <= top:
                norms.add(k)
    seen = set()
    out = []
    for m in range(1, n + 1):
        here = set()
        for k in norms:
            if k <= 3 * m * m:
                g = gcd(k, m * m)
                here.add((k // g, m * m // g))
        out.append(len(here - seen))
        seen |= here
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    print(", ".join(map(str, new_terms(n))))
