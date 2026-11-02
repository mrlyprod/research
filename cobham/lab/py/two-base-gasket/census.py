import subprocess
import sys
import time
from fractions import Fraction
from math import ceil, floor, log, log2


def in_gasket(x, y, b):
    while x or y:
        u, v = x % b, y % b
        if u > 1 or v > 1 or (u and v):
            return False
        x //= b
        y //= b
    return True


def direct(m):
    n = 3**m
    c = 0
    for x in range(n):
        for y in range(n):
            if in_gasket(x, y, 2) and in_gasket(x, y, 3):
                c += 1
    return c


def count(m):
    pw = [3**i for i in range(m + 1)]
    bl = [p.bit_length() for p in pw]
    total = 0

    def walk(k, x, y):
        nonlocal total
        if k == 0:
            if x & y == 0:
                total += 1
            return
        b = bl[k]
        hx = x >> b
        hy = y >> b
        if (hx & hy) and (hx & (hy + 1)) and ((hx + 1) & hy) and ((hx + 1) & (hy + 1)):
            return
        i = k - 1
        p = pw[i]
        walk(i, x, y)
        walk(i, x + p, y)
        walk(i, x, y + p)

    walk(m, 0, 0)
    return total


def axis(m):
    return 2 ** (m + 1) - 1


def terms(lo, hi):
    print("m  C(3^m)        2^(m+1)-1   C/prev    log_3(C)/m  sec")
    prev = None
    for m in range(lo, hi + 1):
        t0 = time.time()
        c = count(m)
        assert c >= axis(m)
        r = c / prev if prev else float("nan")
        e = log(c, 3) / m if m else float("nan")
        print(
            "%-2d %-13d %-11d %-9.5f %-11.6f %.1f"
            % (m, c, axis(m), r, e, time.time() - t0)
        )
        sys.stdout.flush()
        prev = c


def control(hi):
    print("m  direct  pruned  axis  agree")
    for m in range(hi + 1):
        d = direct(m)
        c = count(m)
        print("%-2d %-7d %-7d %-5d %s" % (m, d, c, axis(m), d == c and c >= axis(m)))
        sys.stdout.flush()


def down(x, k=6):
    return floor(x * 10**k) / 10**k


def up(x, k=6):
    return ceil(x * 10**k) / 10**k


def budget():
    da = log2(3)
    db = 1.0
    ax = log(2, 3)
    print("dim A = log_2 3 = %.6f" % da)
    print("dim B = log_3 3 = %.6f" % db)
    print("budget dim A + dim B - 2 <= %.6f" % up(da + db - 2))
    print("axis exponent log_3 2 >= %.6f" % down(ax))
    print("excess >= %.6f" % down(ax - up(da + db - 2)))


def series(hi):
    return [count(m) for m in range(hi + 1)]


def bareiss(rows):
    n = len(rows)
    a = [row[:] for row in rows]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            p = next((i for i in range(k + 1, n) if a[i][k]), None)
            if p is None:
                return 0
            a[k], a[p] = a[p], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
            a[i][k] = 0
        prev = a[k][k]
    return sign * a[n - 1][n - 1]


def solve(rows, rhs, r):
    aug = [
        [Fraction(v) for v in rows[i]] + [Fraction(rhs[i])] for i in range(len(rows))
    ]
    piv = []
    row = 0
    for col in range(r):
        p = next((i for i in range(row, len(aug)) if aug[i][col]), None)
        if p is None:
            continue
        aug[row], aug[p] = aug[p], aug[row]
        d = aug[row][col]
        aug[row] = [v / d for v in aug[row]]
        for i in range(len(aug)):
            if i != row and aug[i][col]:
                f = aug[i][col]
                aug[i] = [aug[i][j] - f * aug[row][j] for j in range(r + 1)]
        piv.append(col)
        row += 1
    for i in range(row, len(aug)):
        if aug[i][r]:
            return "none", None
    c = [Fraction(0)] * r
    for i, col in enumerate(piv):
        c[col] = aug[i][r]
    return ("unique" if len(piv) == r else "free"), c


def charpoly_string(c):
    parts = ["x^%d" % len(c)]
    for i, v in enumerate(c):
        parts.append("- (%s)*x^%d" % (v, len(c) - 1 - i))
    return " ".join(parts)


def dominant_root(c):
    p = charpoly_string(c)
    out = subprocess.run(
        ["gp", "-q"],
        input="v=polroots(%s); vecmax(abs(v))\n" % p,
        capture_output=True,
        text=True,
        check=False,
    )
    return out.stdout.strip()


def hankel(hi, rmax):
    a = series(hi)
    print("terms C(3^m), m = 0..%d" % hi)
    print(",".join(str(v) for v in a))
    print()
    print("k  det H_k (entries a[i+j], i,j < k)")
    for k in range(1, rmax + 2):
        if 2 * k - 1 > len(a):
            break
        print(
            "%-2d %d" % (k, bareiss([[a[i + j] for j in range(k)] for i in range(k)]))
        )
    print()
    print("r  eqs spare fit    coefficients c_1..c_r of a(n) = sum c_i a(n-i), n >= r")
    found = []
    for r in range(1, rmax + 1):
        rows = [[a[n - 1 - j] for j in range(r)] for n in range(r, len(a))]
        rhs = [a[n] for n in range(r, len(a))]
        st, c = solve(rows, rhs, r)
        print(
            "%-2d %-3d %-5d %-6s %s"
            % (
                r,
                len(rows),
                len(rows) - r,
                st,
                "" if c is None else " ".join(str(v) for v in c),
            )
        )
        if st != "none":
            found.append((r, c))
    print()
    if not found:
        print(
            "no linear recurrence with constant coefficients of order r <= %d fits all %d terms"
            % (rmax, len(a))
        )
        print("last reading log_3 C(3^%d) / %d = %.6f" % (hi, hi, log(a[hi], 3) / hi))
        return
    for r, c in found:
        print("order %d characteristic polynomial %s" % (r, charpoly_string(c)))
        rho = dominant_root(c)
        print("order %d dominant root %s" % (r, rho))
        try:
            print("order %d log(root)/log(3) = %.6f" % (r, log(float(rho), 3)))
        except ValueError:
            print("order %d dominant root is not positive real" % r)
        print(
            "order %d last reading log_3 C(3^%d) / %d = %.6f"
            % (r, hi, hi, log(a[hi], 3) / hi)
        )


def first_order(a, rmax):
    for r in range(1, rmax + 1):
        rows = [[a[n - 1 - j] for j in range(r)] for n in range(r, len(a))]
        rhs = [a[n] for n in range(r, len(a))]
        st, c = solve(rows, rhs, r)
        if st != "none":
            return r, c
    return 0, None


def selftest():
    ok = True
    print("%-26s %-26s %-26s %s" % ("check", "expect", "got", "pass"))

    def report(name, want, got):
        nonlocal ok
        hit = str(want) == str(got)
        ok = ok and hit
        print("%-26s %-26s %-26s %s" % (name, want, got, hit))

    report("bareiss 1x1", 3, bareiss([[3]]))
    report("bareiss 2x2", -2, bareiss([[1, 2], [3, 4]]))
    report("bareiss 2x2 pivot swap", -6, bareiss([[0, 2], [3, 4]]))
    report("bareiss 2x2 singular", 0, bareiss([[1, 2], [2, 4]]))
    report("bareiss 3x3", 4, bareiss([[2, 1, 0], [1, 2, 1], [0, 1, 2]]))
    fib = [0, 1]
    while len(fib) < 20:
        fib.append(fib[-1] + fib[-2])
    r, c = first_order(fib, 6)
    report("fibonacci first order", 2, r)
    report("fibonacci charpoly", "x^2 - (1)*x^1 - (1)*x^0", charpoly_string(c))
    mix = [2**n + 3**n + 1 for n in range(20)]
    r, c = first_order(mix, 6)
    report("2^n + 3^n + 1 first order", 3, r)
    report("2^n + 3^n + 1 dominant root", "3.0000000", dominant_root(c)[:9])
    cube = [n**3 + 2**n for n in range(20)]
    r, c = first_order(cube, 8)
    report("n^3 + 2^n first order", 5, r)
    report(
        "n^3 + 2^n charpoly",
        "x^5 - (6)*x^4 - (-14)*x^3 - (16)*x^2 - (-9)*x^1 - (2)*x^0",
        charpoly_string(c),
    )
    print()
    print("all pass %s" % ok)


if __name__ == "__main__":
    v = sys.argv[1] if len(sys.argv) > 1 else "budget"
    if v == "terms":
        terms(int(sys.argv[2]), int(sys.argv[3]))
    elif v == "control":
        control(int(sys.argv[2]))
    elif v == "hankel":
        hankel(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 12)
    elif v == "selftest":
        selftest()
    else:
        budget()
