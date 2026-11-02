import json
import os
import sys
from time import perf_counter

from mpmath import log, mp, mpc, sqrt

mp.dps = 40
P = 14
SHIFT = 60
CUT = 30
T = [[1, 1], [1, 0]]
GAM = [[0, 0], [1, 0]]
LN2 = log(2)

low, mid = [], [[], []]
for n in range(1, 1 << P):
    if n & (n << 1):
        continue
    if n.bit_length() < P:
        low.append(n)
    else:
        mid[n & 1].append(n)


def mv(m, v):
    return [sum(mpc(m[i][j]) * v[j] for j in range(2)) for i in range(2)]


def solve(x, rhs):
    a = [[1 - x * T[0][0], -x * T[0][1]], [-x * T[1][0], 1 - x * T[1][1]]]
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    return [
        (a[1][1] * rhs[0] - a[0][1] * rhs[1]) / det,
        (a[0][0] * rhs[1] - a[1][0] * rhs[0]) / det,
    ]


def epoly(w):
    return [sum(mpc(n) ** (-w) for n in mid[u]) for u in range(2)]


def ladder(s):
    levels = int(mp.ceil(SHIFT - s.real))
    rows = levels + CUT + 2
    g = [[mpc(0), mpc(0)] for _ in range(rows)]
    num = None
    for j in range(levels - 1, -1, -1):
        w = s + j
        acc = epoly(w)
        c = mpc(1)
        for l in range(1, CUT + 1):
            c = c * (-w - (l - 1)) / l
            weight = c * mpc(2) ** (-w - l)
            lift = mv(GAM, g[j + l])
            acc = [acc[i] + weight * lift[i] for i in range(2)]
        if j == 0:
            num = acc
        g[j] = solve(mpc(2) ** (-w), acc)
    return g[0], num


def zeta(s):
    g, _ = ladder(s)
    return sum(mpc(n) ** (-s) for n in low) + g[0] + g[1]


def cofactor(s):
    _, num = ladder(s)
    x = mpc(2) ** (-s)
    det = 1 - x - x * x
    adj = [[mpc(1), x], [x, 1 - x]]
    return det * sum(mpc(n) ** (-s) for n in low) + sum(mv(adj, num))


CONTROL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "control.json")


def control():
    start = perf_counter()
    phi = (1 + sqrt(5)) / 2
    alpha = log(phi) / LN2
    period = 2 * mp.pi / LN2
    half = mp.pi / LN2
    rows = []
    for kind, s in [
        ("zeta", mpc(3)),
        ("zeta", mpc("0.8")),
        ("Z", mpc(3)),
        ("Z", mpc(2)),
        ("Z", mpc("0.8")),
        ("Z", mpc("-0.95", 20)),
        ("Z", mpc(0, 30)),
        ("Z", mpc("-0.693919320380", "9.492981163426")),
        ("Z", mpc("-0.737611737911", "14.343440171295")),
        ("Z", mpc("-0.708500086380", "27.443172842651")),
        ("Z", mpc("0.665353220372", "18.164471712166")),
        ("Z", mpc("-0.315485550754", "23.209873760385")),
        ("residue", mpc(alpha, period)),
        ("residue", mpc(-alpha, half * 5)),
    ]:
        v = zeta(s) if kind == "zeta" else cofactor(s) if kind == "Z" else residue(s)
        rows.append(
            {
                "kind": kind,
                "re": mp.nstr(s.real, 25),
                "im": mp.nstr(s.imag, 25),
                "value": [mp.nstr(v.real, 25), mp.nstr(v.imag, 25)],
            }
        )
    with open(CONTROL, "w") as f:
        json.dump(
            {
                "source": "lab/py/memory-zeta peel.py control",
                "dps": mp.dps,
                "peel": P,
                "shift": SHIFT,
                "cut": CUT,
                "rows": rows,
            },
            f,
            indent=1,
        )
    print("control.json", len(rows), "rows,", round(perf_counter() - start, 2), "s")


def residue(w0):
    _, num = ladder(w0)
    x = mpc(2) ** (-w0)
    adj = [[mpc(1), x], [x, 1 - x]]
    top = sum(mv(adj, num))
    den = -x * LN2 * (-1 - 2 * x)
    return top / den


def main():
    start = perf_counter()
    phi = (1 + sqrt(5)) / 2
    alpha = log(phi) / LN2
    period = 2 * mp.pi / LN2
    half = mp.pi / LN2
    print("code 7, k = 2, dps", mp.dps, "peel", P, "shift", SHIFT, "cut", CUT)
    print("alpha", mp.nstr(alpha, 25))
    for s in [mpc(3), mpc(2), mpc("0.8"), mpc("1.2", 9)]:
        print("zeta", s, mp.nstr(zeta(s), 25))
    for j in range(3):
        print("comb one j", j, mp.nstr(residue(mpc(alpha, period * j)), 25))
    for j in range(3):
        print("comb two j", j, mp.nstr(residue(mpc(-alpha, half * (2 * j + 1))), 25))
    print("runtime", round(perf_counter() - start, 2), "s")


if len(sys.argv) > 1 and sys.argv[1] == "control":
    control()
else:
    main()
