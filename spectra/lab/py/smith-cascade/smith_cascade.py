import time
from fractions import Fraction
from math import comb

PREC = 256
TOP = 511


def jacobsthal(k):
    return (2**k - (-1) ** k) // 3


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def digit_poly(D):
    base = [0] * (2 * D - 1)
    for k in range(D):
        base[2 * k] = comb(D - 1, k)
    return polymul(base, [1, D, 1])


def m_even(D):
    P = digit_poly(D)
    R = (D - 1) // 2
    n = R + 1

    def coefficient(x):
        return P[x] if 0 <= x < len(P) else 0

    M = [[0] * n for _ in range(n)]
    for cp in range(n):
        for c in range(n):
            v = coefficient(c + D - 3 * cp)
            if c:
                v += coefficient(-c + D - 3 * cp)
            M[cp][c] = v
    return M


def fill(D):
    return 2 ** (D - 1) * (D + 2)


def bareiss_det(A):
    n = len(A)
    M = [row[:] for row in A]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            p = next((r for r in range(k + 1, n) if M[r][k]), None)
            if p is None:
                return 0
            M[k], M[p] = M[p], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def v_2(x):
    return (x & -x).bit_length() - 1


def smith_valuations(A, prec=PREC):
    n = len(A)
    cur = prec
    M = [[x % (1 << cur) for x in row] for row in A]
    out = []
    for k in range(n):
        bv = None
        best = None
        for i in range(k, n):
            for j in range(k, n):
                x = M[i][j]
                if x:
                    v = v_2(x)
                    if bv is None or v < bv:
                        bv, best = v, (i, j)
            if bv == 0:
                break
        if best is None:
            raise ValueError("singular modulo 2^%d at step %d" % (cur, k))
        if bv > cur - 64:
            raise ValueError("precision exhausted")
        i0, j0 = best
        M[k], M[i0] = M[i0], M[k]
        if j0 != k:
            for row in M:
                row[k], row[j0] = row[j0], row[k]
        out.append(bv)
        cur -= bv
        m2 = 1 << cur
        u = (M[k][k] >> bv) % m2
        inv = pow(u, -1, m2)
        prow = M[k][k:]
        for i in range(k + 1, n):
            Mi = M[i]
            a = Mi[k] >> bv
            if a % m2:
                f = (a * inv) % m2
                M[i] = Mi[:k] + [(x - f * y) % m2 for x, y in zip(Mi[k:], prow)]
            elif bv:
                M[i] = [x % m2 for x in Mi]
            M[i][k] = 0
    return out


TROUGHS = sorted({2 * jacobsthal(k) + s for k in range(2, 20) for s in (1, 3)})


def tent(D):
    return min(abs(D - t) // 2 + 1 for t in TROUGHS)


def layer(a, j):
    return sum(1 for x in a if x >= j)


def profile(D, prec=PREC):
    a = smith_valuations(m_even(D), prec)
    n = len(a)
    return {
        "D": D,
        "n": n,
        "a": a,
        "v_2": sum(a),
        "L1": layer(a, 1),
        "L2": layer(a, 2),
        "L3": layer(a, 3),
        "L4": layer(a, 4),
        "L5": layer(a, 5),
        "amax": max(a),
        "X": sum(a) - layer(a, 1),
        "tent": tent(D),
    }


def octave(D):
    return D.bit_length() - 1


def cone_checks(rows, key):
    seq = [r[key] for r in rows]
    lipschitz = sum(1 for i in range(len(seq) - 1) if abs(seq[i + 1] - seq[i]) > 1)
    minima = [
        rows[i]["D"]
        for i in range(1, len(seq) - 1)
        if seq[i] <= seq[i - 1]
        and seq[i] <= seq[i + 1]
        and (seq[i] < seq[i - 1] or seq[i] < seq[i + 1])
        and seq[i] != 1
    ]
    plateaus = [
        rows[i]["D"]
        for i in range(1, len(seq) - 1)
        if seq[i - 1] == seq[i] == seq[i + 1] and seq[i] != 1
    ]
    return lipschitz, minima, plateaus


def argmax(rows, key):
    m = max(r[key] for r in rows)
    return m, [r["D"] for r in rows if r[key] == m]


def main():
    t0 = time.time()
    rows = [profile(D) for D in range(3, TOP + 1, 2)]
    print("odd D = 3..%d, rows %d, seconds %.0f" % (TOP, len(rows), time.time() - t0))
    print(
        "tent law nullity = tent(D): %d/%d, misses %s"
        % (
            sum(r["L1"] == r["tent"] for r in rows),
            len(rows),
            [r["D"] for r in rows if r["L1"] != r["tent"]],
        )
    )
    rows = [r for r in rows if r["D"] >= 5]
    print("rows at odd D = 5..%d: %d" % (TOP, len(rows)))
    for D in (5, 7):
        r = next(r for r in rows if r["D"] == D)
        print("D = %d profile %s v_2 = %d n = %d" % (D, r["a"], r["v_2"], r["n"]))
    print()
    print(
        "octave  D range   maxL2 J(k-2)  maxL3 J(k-4)  maxL4  maxL5  max_amax k+4  maxX  max(v_2-ceil(n/3))  at D"
    )
    for k in range(2, 9):
        W = [r for r in rows if octave(r["D"]) == k]
        l2, _ = argmax(W, "L2")
        l3, _ = argmax(W, "L3")
        l4, _ = argmax(W, "L4")
        l5, _ = argmax(W, "L5")
        am, _ = argmax(W, "amax")
        x, _ = argmax(W, "X")
        for r in W:
            r["slack"] = r["v_2"] - (-(-r["n"] // 3))
        sl, at = argmax(W, "slack")
        print(
            "%-7d %3d..%-4d %5d %6d  %5d %6d  %5d  %5d  %8d %3d  %4d  %17d  %s"
            % (
                k,
                W[0]["D"],
                W[-1]["D"],
                l2,
                jacobsthal(k - 2),
                l3,
                jacobsthal(k - 4) if k >= 4 else 0,
                l4,
                l5,
                am,
                k + 4,
                x,
                sl,
                at,
            )
        )
    print()
    for j, key in ((1, "L1"), (2, "L2"), (3, "L3")):
        lip, mins, plat = cone_checks(rows, key)
        print(
            "layer %d cones: Lipschitz breaks %d, interior local minima off 1 %s, plateaus off 1 %s"
            % (j, lip, mins, plat)
        )
    print()
    print(
        "v_2 <= ceil(n/3) + 9 at every row: %s, max slack %d at D = %s"
        % (all(r["slack"] <= 9 for r in rows), *argmax(rows, "slack"))
    )
    print("v_2 > n at D = %s" % [r["D"] for r in rows if r["v_2"] > r["n"]])
    print("v_2 = n at D = %s" % [r["D"] for r in rows if r["v_2"] == r["n"]])
    cls = [r for r in rows if r["D"] % 6 == 1 and r["D"] >= 13]
    ratio = max(Fraction(r["v_2"], r["D"] - 1) for r in cls)
    print(
        "class D = 1 mod 6, rows %d (%d..%d): max(v_2 - n) = %d, all v_2 < D - 1: %s, max v_2/(D-1) = %s at D = %s"
        % (
            len(cls),
            cls[0]["D"],
            cls[-1]["D"],
            max(r["v_2"] - r["n"] for r in cls),
            all(r["v_2"] < r["D"] - 1 for r in cls),
            ratio,
            [r["D"] for r in cls if Fraction(r["v_2"], r["D"] - 1) == ratio],
        )
    )
    last = rows[-1]
    print(
        "D = %d: v_2 = %d, ceil(n/3) + 13 = %d, ceil(n/3) + 9 = %d"
        % (last["D"], last["v_2"], -(-last["n"] // 3) + 13, -(-last["n"] // 3) + 9)
    )
    print(
        "max v_2/n over odd D >= 17: %.3f"
        % max(r["v_2"] / r["n"] for r in rows if r["D"] >= 17)
    )
    print()
    print(
        "first amax > 9 at D = %d, max amax %d at D = %s"
        % (next(r["D"] for r in rows if r["amax"] > 9), *argmax(rows, "amax"))
    )
    print(
        "first L3 != 1 at D = %d, max L3 %d at D = %s"
        % (next(r["D"] for r in rows if r["L3"] != 1), *argmax(rows, "L3"))
    )
    print(
        "first L2 > 5 at D = %d, max L2 %d at D = %s"
        % (next(r["D"] for r in rows if r["L2"] > 5), *argmax(rows, "L2"))
    )
    print("max X %d at D = %s" % argmax(rows, "X"))
    print(
        "X = v_2 - nullity at D = 255, 257: %s"
        % [r["X"] for r in rows if r["D"] in (255, 257)]
    )
    print(
        "second largest divisor valuation, max over rows %d; rows with L5 = 0: %s"
        % (
            max(sorted(r["a"])[-2] for r in rows),
            [r["D"] for r in rows if r["L5"] == 0],
        )
    )
    print(
        "v_2/(D-1) at D = 13, 19: %s"
        % [str(Fraction(r["v_2"], r["D"] - 1)) for r in rows if r["D"] in (13, 19)]
    )
    print()
    bad = [
        D
        for D in range(5, 62, 2)
        if v_2(bareiss_det(m_even(D))) != sum(smith_valuations(m_even(D)))
    ]
    print("exact determinant against Smith sum, odd D = 5..61: mismatches %s" % bad)
    t1 = time.time()
    stable = all(
        smith_valuations(m_even(D), 1024) == next(r["a"] for r in rows if r["D"] == D)
        for D in (255, 257, 511)
    )
    print(
        "profiles at D = 255, 257, 511 unchanged at precision 1024: %s, seconds %.0f"
        % (stable, time.time() - t1)
    )
    E = m_even(7)
    pencil = [[fill(7) * (i == j) - 3 * E[i][j] for j in range(4)] for i in range(4)]
    d7 = bareiss_det(pencil)
    print(
        "D = 7: det(fill I - 3 M_even) = %d, v_2 = %d, fill = %d"
        % (d7, v_2(d7), fill(7))
    )
    print("total seconds %.0f" % (time.time() - t0))


if __name__ == "__main__":
    main()
