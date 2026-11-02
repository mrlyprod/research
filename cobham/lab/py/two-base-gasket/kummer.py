import sys
import time
from math import comb, gcd, log, log2

import numpy as np

PAIRS = ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1))


def coprime6(i, j):
    return gcd(comb(i + j, i), 6) == 1


def direct(m):
    n = 3**m
    return sum(1 for i in range(n) for j in range(n) if coprime6(i, j))


def count(m):
    pw = [3**i for i in range(m + 1)]
    bl = [p.bit_length() for p in pw]
    xs = np.zeros(1, dtype=np.int64)
    ys = np.zeros(1, dtype=np.int64)
    for k in range(m, 0, -1):
        p = pw[k - 1]
        b = bl[k - 1]
        nx = np.concatenate([xs + a * p for a, _ in PAIRS])
        ny = np.concatenate([ys + c * p for _, c in PAIRS])
        hx = nx >> b
        hy = ny >> b
        dead = (hx & hy) != 0
        dead &= (hx & (hy + 1)) != 0
        dead &= ((hx + 1) & hy) != 0
        dead &= ((hx + 1) & (hy + 1)) != 0
        keep = ~dead
        xs = nx[keep]
        ys = ny[keep]
    return int(np.count_nonzero((xs & ys) == 0))


def digits3(n, width):
    out = np.zeros((n, width), dtype=np.int8)
    v = np.arange(n, dtype=np.int64)
    for d in range(width):
        out[:, d] = v % 3
        v //= 3
    return out


def row(n):
    c = 0
    k = n
    while True:
        a, b = k, n
        ok = True
        while b:
            if a % 3 > b % 3:
                ok = False
                break
            a //= 3
            b //= 3
        if ok:
            c += 1
        if k == 0:
            break
        k = (k - 1) & n
    return c


def row_direct(n):
    return sum(1 for k in range(n + 1) if gcd(comb(n, k), 6) == 1)


def box(hi):
    width = 1
    while 3**width < hi + 1:
        width += 1
    d = digits3(hi + 1, width)
    idx = np.arange(hi + 1, dtype=np.int64)
    out = [0]
    total = 0
    for t in range(hi + 1):
        if t == 0:
            total = 1
        else:
            j = idx[:t]
            ok = ((j & t) == 0) & np.all(d[:t] + d[t] <= 2, axis=1)
            total += 2 * int(np.count_nonzero(ok))
        out.append(total)
    return out


def box_direct(n):
    return sum(1 for i in range(n) for j in range(n) if coprime6(i, j))


def terms(lo, hi):
    print("m  K(3^m)        K/prev    log_3(K)/m  sec")
    prev = None
    for m in range(lo, hi + 1):
        t0 = time.time()
        c = count(m)
        r = c / prev if prev else float("nan")
        e = log(c, 3) / m if m else float("nan")
        print("%-2d %-13d %-9.5f %-11.6f %.1f" % (m, c, r, e, time.time() - t0))
        sys.stdout.flush()
        prev = c


def control(hi):
    print("m  direct  walk  agree")
    for m in range(hi + 1):
        d = direct(m)
        c = count(m)
        print("%-2d %-9d %-9d %s" % (m, d, c, d == c))
        sys.stdout.flush()


def rows(hi):
    print("n  row  direct  agree")
    for n in range(hi + 1):
        a = row(n)
        b = row_direct(n)
        print("%-4d %-5d %-5d %s" % (n, a, b, a == b))


def rowlist(hi):
    for n in range(hi + 1):
        print(n, row(n))


def boxcheck(hi):
    b = box(hi)
    bad = 0
    for n in range(min(hi, 60) + 1):
        if box_direct(n) != b[n]:
            bad += 1
            print("mismatch", n, b[n], box_direct(n))
    print(
        "box N = 0..%d checked against binomial gcd, mismatches %d" % (min(hi, 60), bad)
    )
    m = 0
    while 3**m <= hi:
        print("K(3^%d) = %d box %d" % (m, count(m), b[3**m]))
        m += 1


def boxlist(hi):
    b = box(hi)
    for n in range(1, hi + 1):
        print(n, b[n])


def ident(hi):
    print("m  K(3^m)  sum of row over n < 3^m  agree")
    total = 0
    n = 0
    for m in range(hi + 1):
        while n < 3**m:
            total += row(n)
            n += 1
        c = count(m)
        print("%-2d %-9d %-9d %s" % (m, c, total, c == total))
        sys.stdout.flush()


def bfile(hi):
    for n in range(hi + 1):
        print(n, row(n))


def budget():
    da = log2(3)
    db = log(6, 3)
    print(
        "log_2 3 = %.6f  log_3 6 = %.6f  budget log_2 3 + log_3 6 - 2 = %.6f"
        % (da, db, da + db - 2)
    )


if __name__ == "__main__":
    v = sys.argv[1]
    a = [int(x) for x in sys.argv[2:]]
    if v == "budget":
        budget()
        sys.exit()
    {
        "terms": terms,
        "control": control,
        "rows": rows,
        "rowlist": rowlist,
        "boxcheck": boxcheck,
        "boxlist": boxlist,
        "ident": ident,
        "bfile": bfile,
    }[v](*a)
