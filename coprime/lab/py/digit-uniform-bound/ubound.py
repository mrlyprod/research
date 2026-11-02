import itertools
import math
import sys
import time

import numpy as np
from mpmath import iv

FAILS = []

# CHECKS


def chk(value, want, tol):
    ok = abs(value - want) <= tol
    if not ok:
        FAILS.append((value, want))
    return "OK" if ok else f"MISMATCH want {want}"


def chk_le(value, cap, name):
    ok = value <= cap
    if not ok:
        FAILS.append((name, value, cap))
    return "OK" if ok else f"OVER {cap}"


def close(t0):
    if FAILS:
        raise SystemExit(f"{len(FAILS)} rows off: {FAILS[:4]}")
    print(f"all rows agree  ({time.time() - t0:.1f}s)")


# FLOAT TRANSFORM


def hat_missing(q, a0, t):
    th = np.asarray(t) % 1.0
    s = np.sin(np.pi * th)
    d = np.where(
        np.abs(s) < 1e-14,
        float(q),
        np.sin(q * np.pi * th) / np.where(np.abs(s) < 1e-14, 1.0, s),
    )
    ph = 2 * np.pi * (a0 - (q - 1) / 2.0) * th
    return np.sqrt(np.maximum(d * d - 2 * d * np.cos(ph) + 1.0, 0.0)) / (q - 1)


def dirichlet(q, t):
    th = np.asarray(t) % 1.0
    s = np.sin(np.pi * th)
    return np.where(
        np.abs(s) < 1e-14,
        float(q),
        np.abs(np.sin(q * np.pi * th) / np.where(np.abs(s) < 1e-14, 1.0, s)),
    )


def hat_excl(q, E, t):
    th = np.asarray(t) % 1.0
    acc = np.zeros(th.shape, dtype=complex)
    for d in range(q):
        if d not in E:
            acc = acc + np.exp(2j * np.pi * d * th)
    return np.abs(acc) / (q - len(E))


def domination(bases, m=1, pts=4001):
    t0 = time.time()
    t = (np.arange(pts) + 0.5) / pts
    for q in bases:
        u = np.minimum(1.0, (dirichlet(q, t) + float(m)) / (q - m))
        worst = 0.0
        gap = 1.0
        n = 0
        for E in itertools.combinations(range(q), m):
            h = hat_excl(q, E, t)
            worst = max(worst, float((h - u).max()))
            gap = min(gap, float((u - h).max()))
            n += 1
        print(
            f"base {q} at {m} excluded digits: all {n} sets, max(|hat F| - u_q) = {worst:.3e} {chk_le(worst, 1e-12, f'dom {q}')}, every set slack at least {gap:.6f}"
        )
    close(t0)


# THE SUBSET IDENTITY


def runs(E, N):
    out = []
    s = None
    for j in range(N + 1):
        if j < N and j in E:
            if s is None:
                s = j
        elif s is not None:
            out.append((s, j - s))
            s = None
    return out


def sigma_u(q, N, x, m=1):
    i = np.arange(q**N)
    t = x + i / q**N
    p = np.ones(len(i))
    for j in range(N):
        p *= (dirichlet(q, (q**j) * t) + float(m)) / (q - m)
    return float(p.sum())


def subset_sum(q, N, x, m=1):
    i = np.arange(q**N)
    t = x + i / q**N
    tot = 0.0
    for msk in range(1 << N):
        E = {j for j in range(N) if msk >> j & 1}
        p = np.ones(len(i))
        for s, l in runs(E, N):
            p *= dirichlet(q**l, (q**s) * t)
        tot += float(m) ** (N - len(E)) * float(p.sum())
    return tot / (q - m) ** N


def identity(rows, m=1):
    t0 = time.time()
    for q, N in rows:
        for x in (0.0, 0.1234567, 0.5, 1.0 / (3 * q**N), 0.4999):
            a, b = sigma_u(q, N, x, m), subset_sum(q, N, x, m)
            if abs(a - b) > 1e-9 * max(1.0, a):
                FAILS.append((q, N, x, a, b))
        print(
            f"base {q} levels {N} at {m} excluded digits: Sigma_N^u against the subset sum over {2**N} runs decompositions weighted {m}^(N - |E|) {chk(0.0, 0.0, 0.0)}, Sigma_N^u(0) = {sigma_u(q, N, 0.0, m):.6f}"
        )
    close(t0)


# THE LEBESGUE CONSTANT


def lam_at(M, th):
    j = np.arange(M)
    d = np.minimum((j + th) / M, 1.0 - (j + th) / M)
    d = np.where(d <= 0.0, 1.0 / (2 * M), d)
    return float((abs(math.sin(math.pi * th)) / np.sin(np.pi * d)).sum())


def lam_scan(M, grid=4000):
    ths = np.linspace(0.0, 0.5, grid + 1)[1:]
    vals = [lam_at(M, float(t)) for t in ths]
    k = int(np.argmax(vals))
    return vals[k], float(ths[k])


def lebesgue(Ms):
    t0 = time.time()
    c1 = 2.0 / math.pi
    for M in Ms:
        v, arg = lam_scan(M)
        half = lam_at(M, 0.5)
        cap = c1 * math.log(M) + G0 + c1 / M
        print(
            f"M = {M}: L_M = {v:.6f} at offset {arg:.6f}, half-offset value {half:.6f}, L_M/M - (2/pi) log M = {v / M - c1 * math.log(M):.6f} under gamma' + (2/pi)/M = {G0 + c1 / M:.7f} {chk_le(v / M - cap, 0.0, f'lam {M}')} {chk_le(v - half, 1e-9 * M, f'half {M}')}"
        )
    close(t0)


def least_uniform(lo, hi, nd, m=1):
    t0 = time.time()
    rows = [(q, uniform_alpha(q, nd, m)[0]) for q in range(lo, hi)]
    first = next((q for q, e in rows if e < 0.25), None)
    if first is None:
        best = min(rows, key=lambda r: r[1])
        print(
            f"the digit-uniform window bound at {nd} digits and {m} excluded digits clears 1/4 nowhere in [{lo}, {hi}): the best row is q = {best[0]} at alpha_1 < {best[1]:.6f}, so the {m}-digit uniform window ceiling sits above {hi - 1}"
        )
        close(t0)
        return
    bad = [q for q, e in rows if q > first and e >= 0.25]
    print(
        f"the digit-uniform window bound at {nd} digits and {m} excluded digits clears 1/4 first at q = {first}, and at every base above it in [{lo}, {hi}) {chk(float(len(bad)), 0.0, 0.0)}"
    )
    print(
        "  "
        + " ".join(f"q={q}:{e:.6f}" for q, e in rows if first - 4 <= q <= first + 3)
    )
    close(t0)


# THE PEELING BOUND


def lam_max(M, grid=20000):
    return lam_scan(M, grid)[0]


def a_seq(q, N, lams, m=1):
    a = {-1: 1.0, 0: 1.0}
    for n in range(1, N + 1):
        a[n] = (
            m * a[n - 1]
            + m * sum(lams[l] * a[n - 1 - l] for l in range(1, n))
            + lams[n]
        )
    return a


def peel(rows, m=1, grid=400):
    t0 = time.time()
    for q, N in rows:
        lams = {l: lam_max(q**l) / q**l for l in range(1, N + 1)}
        xs = (np.arange(grid) + 0.5) / (grid * q**N)
        i = np.arange(q**N)
        for msk in range(1, 1 << N):
            E = {j for j in range(N) if msk >> j & 1}
            rs = runs(E, N)
            best = 0.0
            for x in xs:
                t = x + i / q**N
                pr = np.ones(len(i))
                for sr, l in rs:
                    pr *= dirichlet(q**l, (q**sr) * t)
                best = max(best, float(pr.sum()))
            cap = q**N * math.prod(lams[l] for _, l in rs)
            if best > cap * (1 + 1e-9):
                FAILS.append((q, N, sorted(E), best, cap))
        a = a_seq(q, N, lams, m)
        sig = max(sigma_u(q, N, float(x), m) for x in xs)
        cap = a[N] * (q / (q - m)) ** N
        if sig > cap * (1 + 1e-9):
            FAILS.append(("chain", q, N, sig, cap))
        print(
            f"base {q} levels {N} at {m} excluded digits: every one of {2**N - 1} run terms under q^N prod lambda {chk(0.0, 0.0, 0.0)}, and max Sigma_N^u = {sig:.6f} under a_N (q/(q-m))^N = {cap:.6f}"
        )
    close(t0)


# THE THRESHOLD

C1 = 2.0 / math.pi
G0 = 0.9625229


def root_z(q, m=1, g0=G0):
    L = math.log(q)
    f = lambda z: (
        (z - m) * (z - 1) ** 2
        - m * (C1 * L * z + g0 * (z - 1) + C1 * (z - 1) ** 2 / (q * z - 1))
    )
    lo, hi = float(m), float(m) + 10.0
    while f(hi) < 0.0:
        hi *= 2.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return hi


def clears(q, m=1, thr=0.25, g0=G0):
    return root_z(q, m, g0) * q / (q - m) < q**thr


def threshold(m=1, thr=0.25, name="1/4", g0=G0, hi=100000):
    t0 = time.time()
    qu = next(
        q
        for q in range(m + 2, hi)
        if all(clears(r, m, thr, g0) for r in range(q, min(q + 2000, hi)))
    )
    bad = [q for q in range(qu, hi) if not clears(q, m, thr, g0)]
    z = root_z(qu, m, g0)
    print(
        f"{m} excluded digits against {name}, gamma' = {g0}: threshold q_u = {qu}, growth root z = {z:.6f} against q^({name})(1 - {m}/q) = {qu**thr * (1 - float(m) / qu):.6f}, alpha_1 < {math.ceil(math.log(z * qu / (qu - m)) / math.log(qu) * 1e6) / 1e6:.6f}"
    )
    print(
        f"  the chain fails at q = {qu - 1} and at no base in [{qu}, {hi}) {chk(float(len(bad)), 0.0, 0.0)}"
    )
    for q in (qu, 2 * qu, 10 * qu, 10**6):
        z = root_z(q, m, g0)
        print(
            f"  q = {q}: alpha_1 < {math.ceil(math.log(z * q / (q - m)) / math.log(q) * 1e6) / 1e6:.6f}"
        )
    close(t0)


# THE THRESHOLD, CERTIFIED


def certify_threshold(qu, hi, m=1, thr=(1, 4), name="1/4", g0=G0):
    t0 = time.time()
    iv.prec = 120
    c1 = 2 / iv.pi
    e = iv.mpf(thr[0]) / thr[1]
    worst = None
    gam = c1 * (iv.euler + iv.log(8 / iv.pi))
    if not float(gam.b) <= g0:
        raise SystemExit("gamma' constant is below the true Lebesgue constant")

    def margin(q):
        w = iv.mpf(q) ** e * (1 - iv.mpf(m) / q)
        return (w - m) * (w - 1) ** 2 - iv.mpf(m) * (
            c1 * iv.log(q) * w + iv.mpf(g0) * (w - 1) + c1 * (w - 1) ** 2 / (q * w - 1)
        )

    for q in range(qu, hi):
        g = float(margin(q).a)
        if g <= 0.0:
            FAILS.append((q, g))
        if worst is None or g < worst[1]:
            worst = (q, g)
    print(
        f"certified at 120 bits at {m} excluded digits, gamma' = {g0} above the true {float(gam.b):.10f}: (w - {m})(w - 1)^2 > {m}((2/pi) log q w + gamma'(w - 1) + (2/pi)(w-1)^2/(q w - 1)) at w = q^({name})(1 - {m}/q) for every q in [{qu}, {hi}), tightest margin {worst[1]:.6e} at q = {worst[0]} {chk(float(len(FAILS)), 0.0, 0.0)}"
    )
    print(
        f"  and the chain fails at q = {qu - 1}, margin {float(margin(qu - 1).b):.6e}, so {qu} is the least threshold this chain gives at {m} excluded digits"
    )
    close(t0)


# THE UNIFORM WINDOW BOUND


def dn(x):
    return np.nextafter(np.asarray(x, float), -np.inf)


def up(x):
    return np.nextafter(np.asarray(x, float), np.inf)


def sup_sin(A, B):
    hit = np.ceil(A - 0.5) <= np.floor(B - 0.5)
    ends = np.maximum(np.abs(np.sin(np.pi * A)), np.abs(np.sin(np.pi * B)))
    return np.where(hit, 1.0, up(up(ends) * (1.0 + 2.0**-45)))


def uniform_G(q, nd, m=1):
    W = q**nd
    w = np.arange(W, dtype=np.float64)
    delta = np.minimum(w, W - 1 - w) / W
    sn = dn(dn(np.sin(np.pi * delta)) * (1.0 - 2.0**-45))
    with np.errstate(divide="ignore"):
        S2 = np.where(delta <= 0.0, np.inf, up(1.0 / np.maximum(sn, 0.0)))
    S1 = sup_sin(w / q ** (nd - 1), (w + 1) / q ** (nd - 1))
    D = np.minimum(float(q), up(S1 * S2))
    return np.minimum(1.0, up(up(D + float(m)) / (q - m)))


def perron_up(G, q, nd, iters=20000):
    S = q ** (nd - 1)
    tgt = np.tile(np.arange(S), q)
    y = np.ones(S)
    lam = 0.0
    streak = 0
    for it in range(iters):
        z = (G * y[tgt]).reshape(S, q).sum(axis=1)
        nl = float(z.max())
        y = z / nl
        streak = streak + 1 if abs(nl - lam) <= 1e-13 * nl else 0
        lam = nl
        if streak >= 50 and it >= 300:
            break
    y = np.maximum(y, 1e-9 * float(y.max()))
    r = up(G * y[tgt]).reshape(S, q)
    acc = np.zeros(S)
    for c in range(q):
        acc = up(acc + r[:, c])
    return float(up(acc / y).max())


def uniform_alpha(q, nd, m=1):
    mu = perron_up(uniform_G(q, nd, m), q, nd)
    e = math.ceil(math.log(mu) / math.log(q) * 1e6) / 1e6
    if not mu < q**e:
        raise SystemExit("rounding unsafe")
    return e, mu


def window(rows, m=1):
    t0 = time.time()
    for q, nd in rows:
        e, mu = uniform_alpha(q, nd, m)
        v = "CLEARS 1/4" if e < 0.25 else "does not clear 1/4"
        print(
            f"base {q} uniform windows {nd} digits at {m} excluded digits: lambda < {mu:.6f} -> alpha_1 < {e:.6f} at every excluded set  [{v}]"
        )
    close(t0)


def compare(rows):
    t0 = time.time()
    for q, nd, a0, per in rows:
        e, _mu = uniform_alpha(q, nd)
        ok = e >= per
        if not ok:
            FAILS.append((q, e, per))
        print(
            f"base {q} missing {a0}: uniform {e:.6f} against the per-digit ladder {per:.6f}, uniform weaker by {e - per:.6f} {chk(1.0 if ok else 0.0, 1.0, 0.0)}"
        )
    close(t0)


# VERBS


def main():
    verb = sys.argv[1] if len(sys.argv) > 1 else "check"
    if verb == "check":
        domination([3, 4, 10, 20, 21, 34, 76, 99, 200])
        identity([(3, 4), (5, 3), (10, 3), (21, 2), (34, 2)])
        lebesgue([4, 21, 34, 76, 99, 200, 1000, 100000])
        peel([(3, 4), (5, 3), (10, 3)])
        threshold()
        certify_threshold(125, 3000)
        compare([(21, 4, 0, 0.250088), (34, 4, 16, 0.249371)])
        domination([4, 10, 21, 34], 2)
        identity([(4, 3), (10, 2), (21, 2)], 2)
        peel([(4, 3), (5, 3)], 2)
        threshold(2)
        certify_threshold(649, 3000, 2)
    if verb == "domination":
        domination([3, 4, 10, 20, 21, 34, 76, 99, 200])
    if verb == "peel":
        peel([(3, 4), (4, 4), (5, 3), (7, 3), (10, 3)])
    if verb == "least":
        least_uniform(35, 90, 3)
    if verb == "window":
        window(
            [
                (21, 2),
                (21, 3),
                (21, 4),
                (34, 3),
                (34, 4),
                (76, 3),
                (99, 3),
                (126, 3),
                (200, 3),
            ]
        )
    if verb == "compare":
        compare(
            [
                (9, 5, 0, 0.323432),
                (10, 5, 5, 0.350684),
                (20, 4, 6, 0.283414),
                (21, 4, 0, 0.250088),
                (33, 4, 15, 0.253000),
                (34, 4, 16, 0.249371),
            ]
        )
    if verb == "threshold":
        threshold()
        certify_threshold(125, 3000)
    if verb == "pairs":
        domination([4, 6, 10, 21, 34], 2)
        identity([(4, 3), (5, 3), (10, 2), (21, 2), (34, 2)], 2)
        peel([(4, 3), (5, 3), (7, 3)], 2)
        threshold(2)
        certify_threshold(649, 3000, 2)
        threshold(2, 1 / 3, "1/3")
        certify_threshold(105, 3000, 2, (1, 3), "1/3")
        threshold(1, 1 / 3, "1/3")
        certify_threshold(32, 3000, 1, (1, 3), "1/3")
        least_uniform(100, 300, 2, 2)
    if verb == "lebesgue":
        lebesgue([4, 9, 21, 34, 76, 99, 200, 1000, 5776, 100000])
    if verb == "identity":
        identity([(3, 4), (5, 3), (10, 3), (21, 2), (34, 2)])


if __name__ == "__main__":
    main()
