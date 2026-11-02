import argparse
import json
import os
import resource
import sys
import time
from math import comb

from mpmath import cos, exp, floor, log, mp, mpf, nint, pi

# DIGIT POLYNOMIAL


def digit_poly(D):
    base = [0] * (2 * D - 1)
    for k in range(D):
        base[2 * k] = comb(D - 1, k)
    out = [0] * (2 * D + 1)
    for j, v in enumerate(base):
        out[j] += v
        out[j + 1] += D * v
        out[j + 2] += v
    return out


def fill_value(D):
    return 2 ** (D - 1) * (D + 2)


def half_width(D):
    return (D - 1) // 2


def coefficient(P, x):
    return P[x] if 0 <= x < len(P) else 0


def sign_vector(D):
    return [
        (3 if c % 3 == 0 else 0) - 1 for c in range(-half_width(D), half_width(D) + 1)
    ]


# STRUCTURAL ASSERTIONS


def assert_nonnegative(D):
    assert all(v >= 0 for v in digit_poly(D)), (
        "digit poly has a negative coefficient at D=%d" % D
    )


def assert_window_closed(D):
    P = digit_poly(D)
    H = half_width(D)
    for c in range(-H, H + 1):
        for s, val in enumerate(P):
            if not val or (c + D - s) % 3:
                continue
            assert abs((c + D - s) // 3) <= H, "carry window escapes at D=%d" % D


def assert_palindrome(D):
    P = digit_poly(D)
    assert all(P[j] == P[2 * D - j] for j in range(2 * D + 1)), (
        "digit poly not palindromic at D=%d" % D
    )


def assert_colsum(D):
    P = digit_poly(D)
    H = half_width(D)
    f = fill_value(D)
    for c in range(-H, H + 1):
        cs = sum(coefficient(P, c + D - 3 * cp) for cp in range(-H, H + 1))
        v = (3 if c % 3 == 0 else 0) - 1
        assert f - 3 * cs == (D - 1) * v, "colsum identity fails at D=%d c=%d" % (D, c)


# FOLDED CORE, COLUMNS OF THE TRANSPOSE


def folded_columns(D):
    P = digit_poly(D)
    H = half_width(D)
    cols = []
    for c in range(H + 1):
        col = []
        for cp in range(H + 1):
            v = coefficient(P, c + D - 3 * cp)
            if c:
                v += coefficient(P, -c + D - 3 * cp)
            if v:
                col.append((cp, v))
        cols.append(col)
    return cols


def functional(D):
    out = []
    for c in range(half_width(D) + 1):
        v = (3 if c % 3 == 0 else 0) - 1
        out.append(v if c == 0 else 2 * v)
    return out


# UNFOLDED CENSUS, THE INDEPENDENT ROUTE


def census_V(D, top):
    P = digit_poly(D)
    H = half_width(D)
    n = 2 * H + 1
    step = []
    for cp in range(-H, H + 1):
        step.append(
            [
                (c + H, coefficient(P, c + D - 3 * cp))
                for c in range(-H, H + 1)
                if coefficient(P, c + D - 3 * cp)
            ]
        )
    u = [0] * n
    u[H] = 1
    v = sign_vector(D)
    out = [sum(v[i] * u[i] for i in range(n))]
    for _ in range(top):
        w = [sum(x * u[j] for j, x in row) for row in step]
        u = w
        out.append(sum(v[i] * u[i] for i in range(n)))
    return out


# UNFOLDED ROW SWEEP, THE INDEPENDENT ROUTE


def unfolded_rows(D, top):
    P = digit_poly(D)
    H = half_width(D)
    cols = []
    for c in range(-H, H + 1):
        cols.append(
            [
                (cp + H, coefficient(P, c + D - 3 * cp))
                for cp in range(-H, H + 1)
                if coefficient(P, c + D - 3 * cp)
            ]
        )
    r = sign_vector(D)
    out = [r]
    for _ in range(top):
        r = [sum(x * r[j] for j, x in col) for col in cols]
        out.append(r)
    return out


# ROW SWEEP


def row_sweep(D, cap):
    assert_nonnegative(D)
    assert_window_closed(D)
    assert_palindrome(D)
    cols = folded_columns(D)
    r = functional(D)
    n = len(r)
    last_neg = -1
    for k in range(1, cap + 1):
        nxt = [0] * n
        for c in range(n):
            s = 0
            for cp, w in cols[c]:
                x = r[cp]
                if x:
                    s += w * x
            nxt[c] = s
        r = nxt
        if r[0] < 0:
            last_neg = k
        lo = min(r)
        if lo >= 0:
            return k, last_neg, lo > 0
    raise AssertionError("no nonnegative row within cap at D=%d" % D)


# CROSSING LEVEL


def k_star(D, levels=400):
    log_r = mpf(0)
    s = mpf(0)
    for i in range(2, levels):
        a = cos(pi / mpf(3) ** i)
        b = cos(2 * pi / mpf(3) ** i)
        log_r += log(a / b)
        s += log((D - 2 * a) / (D - 2)) - log((D + 2 * b) / (D + 2))
    return ((D - 1) * log_r + s) / log(mpf(D + 2) / mpf(D - 2))


def l_zero(D, guard):
    ks = k_star(D)
    near = abs(ks - nint(ks))
    assert near > guard, "K* within %s of an integer at D=%d" % (guard, D)
    L0 = int(floor(ks))
    return (L0 if L0 % 2 else L0 - 1), ks, near


# WINDOW


def envelope(D, K):
    e = 4 * (K - 1) * (
        mpf("0.7528157") ** (D - 1) + mpf("0.7052518") ** (D - 1)
    ) + 2 * mpf("0.2266816") ** (D - 1)
    return e * exp(2 * (K + 1) * mpf("0.8900159") ** (D - 1))


def window(lo, hi):
    inside = []
    print(
        "  D  K*              L0  E           A          w           dist(K*, odd)  inside"
    )
    for D in range(lo, hi + 1, 2):
        ks = k_star(D)
        L0 = int(floor(ks))
        L0 = L0 if L0 % 2 else L0 - 1
        A = log(mpf(D + 2) / (D - 2))
        E = envelope(D, int(floor(ks)) + 2)
        tail = (mpf("1.9") * (D - 1) + mpf(7) / (D - 2)) * mpf(9) ** (-L0)
        w = (-log(1 - E) + tail) / A
        dist = min(abs(ks - L0), abs(ks - L0 - 2))
        hit = dist <= w
        if hit:
            inside.append(D)
        print(
            "%3d  %-14s %4d  %-10s  %-9s  %-10s  %-13s  %s"
            % (
                D,
                mp.nstr(ks, 10),
                L0,
                mp.nstr(E, 4),
                mp.nstr(A, 5),
                mp.nstr(w, 4),
                mp.nstr(dist, 4),
                "YES" if hit else "no",
            )
        )
    print(
        "even D in %d..%d with an odd integer inside the window: %s"
        % (lo, hi, inside or "none")
    )


def values(D, lo, hi):
    assert_nonnegative(D)
    assert_window_closed(D)
    assert_palindrome(D)
    cols = folded_columns(D)
    r = functional(D)
    n = len(r)
    for k in range(1, hi + 1):
        r = [sum(w * r[cp] for cp, w in cols[c]) for c in range(n)]
        if k >= lo:
            v = r[0]
            print(
                "D=%d  V(%d) = %s%s... (%d digits)"
                % (D, k, "-" if v < 0 else "+", str(abs(v))[:12], len(str(abs(v))))
            )


# SELF TEST


def selftest():
    for D in range(4, 62, 2):
        assert_colsum(D)
    print(
        "colsum identity fill - 3 colsum(c) = (D-1) v_c bites at even D = 4..60; the theorem is prop:mass"
    )
    for D in range(4, 37, 2):
        H = half_width(D)
        cols = folded_columns(D)
        r = functional(D)
        n = len(r)
        top = 4 * D + 8
        mine = [r[0]]
        folded = [r]
        for _ in range(top):
            r = [sum(w * r[cp] for cp, w in cols[c]) for c in range(n)]
            folded.append(r)
            mine.append(r[0])
        assert mine == census_V(D, top), (
            "folded sweep disagrees with the census at D=%d" % D
        )
        full = unfolded_rows(D, top)
        for k in range(top + 1):
            a, b = folded[k], full[k]
            assert all(b[H - c] == b[H + c] for c in range(H + 1)), (
                "unfolded row not symmetric at D=%d k=%d" % (D, k)
            )
            assert all(a[c] == b[H + c] * (1 if c == 0 else 2) for c in range(H + 1)), (
                "folded row disagrees with the unfolded row at D=%d k=%d" % (D, k)
            )
    print(
        "folded sweep V(L) = 3 m0(L) - b(L) against the unfolded census: even D = 4..36, L <= 4D+8, exact"
    )
    print(
        "folded row against the unfolded row on all of S, entrywise: even D = 4..36, L <= 4D+8, exact"
    )


# RUN


def run(lo, hi, cap_slope, out_path, budget):
    rows = []
    if out_path and os.path.exists(out_path):
        with open(out_path) as fh:
            rows = json.load(fh)
    done = {row["D"] for row in rows}
    guard = mpf(10) ** -20
    start = time.time()
    for D in range(lo, hi + 1, 2):
        if D in done:
            continue
        L0, ks, near = l_zero(D, guard)
        t0 = time.time()
        t, last_neg, strict = row_sweep(D, int(cap_slope * D * D) + 200)
        dt = time.time() - t0
        rows.append(
            {
                "D": D,
                "Lstar": last_neg,
                "t": t,
                "L0": L0,
                "strict": strict,
                "kstar": mp.nstr(ks, 14),
                "gap": mp.nstr(near, 6),
                "agree": last_neg == L0,
                "seconds": round(dt, 3),
            }
        )
        print(
            "D=%3d  L*=%5d  L0=%5d  t=%5d  K*=%-16s  %s  %s  %.2fs"
            % (
                D,
                last_neg,
                L0,
                t,
                mp.nstr(ks, 12),
                "ok" if last_neg == L0 else "MISS",
                "strict" if strict else "SLACK",
                dt,
            ),
            flush=True,
        )
        if out_path:
            with open(out_path, "w") as fh:
                json.dump(rows, fh, indent=1)
        if budget and time.time() - start > budget:
            print("budget wall after D=%d" % D, flush=True)
            break
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lo", type=int, default=6)
    ap.add_argument("--hi", type=int, default=120)
    ap.add_argument("--cap-slope", type=float, default=0.08)
    ap.add_argument("--out", default="")
    ap.add_argument("--budget", type=float, default=0.0)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--window", nargs=2, type=int)
    ap.add_argument("--values", nargs=3, type=int)
    a = ap.parse_args()
    mp.dps = 60
    if a.window:
        window(*a.window)
        return 0
    if a.values:
        values(*a.values)
        return 0
    if a.selftest:
        selftest()
    rows = run(a.lo, a.hi, a.cap_slope, a.out, a.budget)
    bad = [r["D"] for r in rows if not r["agree"]]
    slack = [r["D"] for r in rows if not r["strict"]]
    gap = min(rows, key=lambda r: float(r["gap"]))
    steps = sorted({r["t"] - r["Lstar"] for r in rows})
    print(
        "rows %d  agree %d  misses %s"
        % (len(rows), sum(r["agree"] for r in rows), bad or "none")
    )
    print("t - L* values %s  non-strict rows %s" % (steps, slack or "none"))
    print("least K* distance to an integer: %s at D = %d" % (gap["gap"], gap["D"]))
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak = peak / 1048576.0 if sys.platform == "darwin" else peak / 1024.0
    print(
        "total sweep seconds %.1f, deepest row D = %d at %.1f s, peak resident %.1f MB"
        % (sum(r["seconds"] for r in rows), rows[-1]["D"], rows[-1]["seconds"], peak)
    )
    return 0 if not bad and not slack else 1


if __name__ == "__main__":
    sys.exit(main())
