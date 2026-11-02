import argparse
import math
import sys
import time
from fractions import Fraction

import numpy as np

A002487 = [
    0,
    1,
    1,
    2,
    1,
    3,
    2,
    3,
    1,
    4,
    3,
    5,
    2,
    5,
    3,
    4,
    1,
    5,
    4,
    7,
    3,
    8,
    5,
    7,
    2,
    7,
    5,
    8,
    3,
    7,
    4,
    5,
]
PHI = (1 + 5**0.5) / 2


def product(q, b, L):
    r = [1]
    for j in range(L):
        s = b**j
        out = [0] * (len(r) + s * (len(q) - 1))
        for k, w in enumerate(q):
            if w:
                for n, v in enumerate(r):
                    out[n + k * s] += w * v
        r = out
    return r


def digits(q, b, L, n, memo):
    if L == 0:
        return 1 if n == 0 else 0
    key = (L, n)
    if key not in memo:
        memo[key] = sum(
            w * digits(q, b, L - 1, (n - c) // b, memo)
            for c, w in enumerate(q)
            if w and c <= n and (n - c) % b == 0
        )
    return memo[key]


def stern(m):
    s = [0, 1]
    for n in range(2, m + 1):
        s.append(s[n // 2] if n % 2 == 0 else s[n // 2] + s[n // 2 + 1])
    return s


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def count(a):
    q = [1, 1, 1]
    s = stern(2 ** (a.top + 1) + 2)
    assert s[: len(A002487)] == A002487
    print(f"count: first {len(A002487)} terms of A002487 agree with the recursion")
    bad = 0
    for L in range(1, a.exact + 1):
        r = product(q, 2, L)
        memo = {}
        bad += sum(r[n] != digits(q, 2, L, n, memo) for n in range(len(r)))
    print(
        f"count: product formula against the digit recursion, Q = 1+z+z^2, b = 2, L = 1..{a.exact}: {bad} mismatches"
    )
    for Q, b in [([1, 1], 2), ([2, 1], 2), ([1, 1, 1, 1], 3), ([2, 2, 3, 2, 1], 2)]:
        memo = {}
        r = product(Q, b, a.exact - 4)
        print(
            f"count: Q = {Q}, b = {b}, L = {a.exact - 4}: {sum(r[n] != digits(Q, b, a.exact - 4, n, memo) for n in range(len(r)))} mismatches over {len(r)} lags"
        )
    print(
        "L lags total mean max F_(L+1) argmax single single_law stern_mismatch mirror_mismatch fold_mismatch peak/mean ratio"
    )
    prev = None
    for L in range(1, a.top + 1):
        r = product(q, 2, L)
        N = len(r)
        assert N == 2 ** (L + 1) - 1 and sum(r) == 3**L
        sm = sum(r[n] != s[n + 1] for n in range(2**L))
        mm = sum(r[n] != r[N - 1 - n] for n in range(N))
        top = max(r)
        arg = [n for n in range(N) if r[n] == top]
        m = (2 ** (L + 1) + (-1) ** L) // 3
        law = sorted({m - 1, 3 * 2 ** (L - 1) - m - 1, N - m, N - 3 * 2 ** (L - 1) + m})
        one = [n for n in range(N) if r[n] == 1]
        onelaw = sorted(
            {2**k - 1 for k in range(L + 1)}
            | {2 ** (L + 1) - 1 - 2**k for k in range(L + 1)}
        )
        fm = sum(
            s[2**L + j] != (r[j - 1] if j else 0) + r[j - 1 + 2**L] for j in range(2**L)
        )
        assert (
            top == fib(L + 1)
            and arg == law
            and one == onelaw
            and len(one) == 2 * L + 1
            and sm == 0
            and mm == 0
            and fm == 0
        )
        pm = Fraction(top * N, 3**L)
        ratio = float(pm / prev) if prev else float("nan")
        prev = pm
        shown = arg if len(arg) <= 4 else arg[:4]
        print(
            f"{L} {N} {3**L} {3**L / N:.4f} {top} {fib(L + 1)} {shown} {len(one)} {2 * L + 1} {sm} {mm} {fm} {float(pm):.4f} {ratio:.5f}"
        )
    print(f"count: 2 phi/3 = {2 * PHI / 3:.5f}")


CYCLO = {}


def cyclo(n):
    if n not in CYCLO:
        p = [Fraction(-1)] + [Fraction(0)] * (n - 1) + [Fraction(1)]
        for d in range(1, n):
            if n % d == 0:
                p = divide(p, cyclo(d))[0]
        CYCLO[n] = p
    return CYCLO[n]


def divide(p, d):
    p = list(p)
    out = [Fraction(0)] * max(1, len(p) - len(d) + 1)
    while len(p) >= len(d) and any(p):
        c = p[-1] / d[-1]
        k = len(p) - len(d)
        out[k] = c
        for i, v in enumerate(d):
            p[k + i] -= c * v
        p.pop()
    return out, p


def vanishes(q, n):
    return not any(divide(q, cyclo(n))[1])


def verdict(q, b, depth):
    q = [Fraction(x) for x in q]
    while q and q[-1] == 0:
        q.pop()
    for k in range(1, b**depth):
        if k % b == 0:
            continue
        if not any(vanishes(q, b**i // math.gcd(k, b**i)) for i in range(1, depth + 1)):
            return False, k
    return True, None


def fourier(q, b, t, terms=80):
    q = np.array([float(x) for x in q])
    z = 1 + 0j
    for i in range(1, terms):
        w = np.exp(2j * np.pi * t / b**i) ** np.arange(len(q))
        z *= (q @ w) / q.sum()
    return abs(z)


def family():
    rows = []
    for b in [2, 3, 4]:
        for K in range(2, 9):
            rows.append((f"uniform K={K}", [1] * K, b))
    for c in ["1/9", "1/4", "1/2", "1", "2", "4"]:
        cf = Fraction(c)
        rows.append((f"1+c(1+z) c={c}", [1 + cf, cf], 2))
    for c in ["1/81", "1/9", "1/4", "1/2", "1", "2", "4"]:
        cf = Fraction(c)
        rows.append((f"1+c(1+z+z^2)^2 c={c}", [1 + cf, 2 * cf, 3 * cf, 2 * cf, cf], 2))
    rows.append(("(1+z+z^2)^2", [1, 2, 3, 2, 1], 2))
    return rows


def limit(a):
    print("Q b verdict witness_k max|hat mu(k)|,k<b^3 |hat mu(1)|")
    agree = 0
    for name, q, b in family():
        ok, k = verdict(q, b, a.depth)
        top = max(fourier(q, b, t) for t in range(1, b**3) if t % b)
        print(
            f"{name} {b} {'absolutely continuous' if ok else 'singular'} {k if k else '-'} {top:.2e} {fourier(q, b, 1):.4f}"
        )
        assert ok == (top < 1e-9)
        agree += 1
    print(f"limit: {agree} rows, verdict agrees with the Fourier product on every row")
    s2 = 2**0.5
    q = np.polymul(np.polymul([1, -s2, 1], [1, 0, s2, 0, 1]), [1, 2, 3, 2, 1])[::-1]
    at = lambda t: abs(np.polyval(q[::-1], np.exp(2j * np.pi * t)))
    cover = [
        min(i for i in range(1, a.depth + 1) if at(k / 2**i) < 1e-12)
        for k in range(1, 2**a.depth, 2)
    ]
    off = [max(at(k / 2**i) for k in range(1, 2**i, 2)) for i in range(1, 5)]
    top = max(fourier(q, 2, t) for t in range(1, 256, 2))
    print(
        f"limit: real Q = (z^2 - sqrt2 z + 1)(z^4 + sqrt2 z^2 + 1)(1+z+z^2)^2, b = 2: min coefficient {q.min():.4f}, cover level of odd k < 2^{a.depth} {sorted(set(cover))}, max |Q| at primitive 2^i-th roots i = 1..4 {[round(float(x), 4) for x in off]}, max|hat mu(t)| over odd t < 256 {top:.1e}"
    )
    for name, q, b in [
        ("uniform K=3", [1, 1, 1], 2),
        ("uniform K=4", [1, 1, 1, 1], 2),
        ("1+(1+z+z^2)^2", [2, 2, 3, 2, 1], 2),
        ("1+(1+z)/4", [5, 1], 2),
    ]:
        cells = []
        for L in [8, 12, 16, 20]:
            r = np.array(product(q, b, L), dtype=float)
            p = np.sort(r)[::-1] / r.sum()
            cells.append(f"L={L}: {np.searchsorted(np.cumsum(p), 0.5) / len(p):.4f}")
        print(
            f"limit: {name}, share of lags carrying half the mass, {', '.join(cells)}"
        )


def stack(rng, C, L, b, K, sigma, resid):
    A = rng.normal(0, sigma, (L, K, C, C))
    if resid:
        A[:, 0] += np.eye(C)
    return A


def grad_filter(A, u, b, N):
    L, K, C, _ = A.shape
    g = np.zeros((N, C))
    g[0] = u
    for j in range(L - 1, -1, -1):
        s = b**j
        new = np.zeros_like(g)
        for k in range(K):
            new[k * s :] += g[: N - k * s] @ A[j, k]
        g = new
    return g


def gradient(a):
    rng = np.random.default_rng(1)
    C, L, b = a.channels, a.levels, 2
    for name, K, resid, cs in [
        ("pure K=3", 3, False, 1.0),
        ("residual 1+c(1+z)", 2, True, 0.5),
    ]:
        sigma = (cs / C) ** 0.5
        q = [1 + cs, cs] if resid else [cs] * K
        exact = np.array(
            product([Fraction(x).limit_denominator(1000) for x in q], b, L), dtype=float
        )
        N = len(exact)
        acc = np.zeros(N)
        sq = np.zeros(N)
        u = np.zeros(C)
        u[0] = 1.0
        for _ in range(a.draws):
            g = (grad_filter(stack(rng, C, L, b, K, sigma, resid), u, b, N) ** 2).sum(1)
            acc += g
            sq += g * g
        mean = acc / a.draws
        se = np.sqrt(np.maximum(sq / a.draws - mean**2, 0) / a.draws)
        z = np.abs(mean - exact) / np.maximum(se, 1e-300)
        rel = np.abs(mean - exact) / exact
        print(
            f"gradient: {name}, C={C}, L={L}, {a.draws} draws, {N} lags: max |z| {z.max():.2f}, share |z|>3 {np.mean(z > 3):.4f}, median rel err {np.median(rel):.4f}, corr {np.corrcoef(mean, exact)[0, 1]:.4f}"
        )


def train(A, u, v, n, lr, cap, hit):
    L, K, C, _ = A.shape
    b = 2
    N = 2 ** (L + 1) - 1
    for step in range(1, cap + 1):
        F = [np.zeros((N, C))]
        F[0][0] = v
        for j in range(L):
            s = b**j
            new = np.zeros((N, C))
            for k in range(K):
                new[k * s :] += F[j][: N - k * s] @ A[j, k].T
            F.append(new)
        f = F[L] @ u
        if f[n] >= hit:
            return step
        e = 2 * f
        e[n] -= 2
        G = np.outer(e, u)
        gu = e @ F[L]
        gA = np.zeros_like(A)
        for j in range(L - 1, -1, -1):
            s = b**j
            back = np.zeros((N, C))
            for k in range(K):
                gA[j, k] = G[k * s :].T @ F[j][: N - k * s]
                back[: N - k * s] += G[k * s :] @ A[j, k]
            G = back
        gv = G[0]
        A -= lr * gA
        u -= lr * gu
        v -= lr * gv
    return None


def copy(a):
    L, C = a.levels, a.channels
    r = product([1, 1, 1], 2, L)
    lags = []
    for want in [1, 2, 3, 5, 8, 13, 21, 34, 55]:
        pick = [n for n in range(2 ** (L - 1), 2**L) if r[n] == want]
        if pick:
            lags.append(pick[len(pick) // 2])
    print(
        f"copy: K=3, b=2, L={L}, C={C}, lags {lags}, r {[r[n] for n in lags]}, seeds {a.seeds}, lr {a.lr}, hit f(n) >= {a.hit}"
    )
    rows = []
    capped = 0
    for n in lags:
        steps = []
        for seed in range(a.seeds):
            rng = np.random.default_rng(100 + seed)
            A = rng.normal(0, a.scale * (1 / (3 * C)) ** 0.5, (L, 3, C, C))
            u = rng.normal(0, 1 / C**0.5, C)
            v = rng.normal(0, 1, C)
            steps.append(train(A, u, v, n, a.lr, a.cap, a.hit))
        med = float(np.median([x if x else math.inf for x in steps]))
        capped += steps.count(None)
        rows.append((n, r[n], med))
        print(f"copy: lag {n} r {r[n]} steps {steps} median {med}")
    x = np.log([1 / t[1] for t in rows])
    y = np.log([t[2] for t in rows])
    good = np.isfinite(y)
    beta, _alpha = np.polyfit(x[good], y[good], 1)
    print(
        f"copy: {capped} runs capped at {a.cap} steps, read as above the cap in each median"
    )
    print(
        f"copy: log steps against log 1/r over {good.sum()} lags: slope {beta:.3f}, corr {np.corrcoef(x[good], y[good])[0, 1]:.3f}, steps ratio r=1 over r=max {rows[0][2] / rows[-1][2]:.2f} against r ratio {rows[-1][1] / rows[0][1]}"
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("verb", choices=["count", "limit", "gradient", "copy", "all"])
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--exact", type=int, default=12)
    p.add_argument("--depth", type=int, default=6)
    p.add_argument("--channels", type=int, default=8)
    p.add_argument("--levels", type=int, default=8)
    p.add_argument("--draws", type=int, default=4000)
    p.add_argument("--seeds", type=int, default=5)
    p.add_argument("--lr", type=float, default=0.005)
    p.add_argument("--scale", type=float, default=1.0)
    p.add_argument("--cap", type=int, default=20000)
    p.add_argument("--hit", type=float, default=0.5)
    a = p.parse_args()
    verbs = ["count", "limit", "gradient", "copy"] if a.verb == "all" else [a.verb]
    for verb in verbs:
        t = time.time()
        globals()[verb](a)
        print(f"{verb}: {time.time() - t:.1f} seconds")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
