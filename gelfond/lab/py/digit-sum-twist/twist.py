import sys
import time
from math import gcd, log, log10, pi, sin, sqrt

import numpy as np

GAMMA = 0.5772156649015329

# DESIGNS


def digits_of(b, e):
    return np.array([d for d in range(b) if d not in e], dtype=np.int64)


def gabs(F, v):
    re = np.zeros_like(v)
    im = np.zeros_like(v)
    for d in F:
        re += np.cos(2 * pi * d * v)
        im += np.sin(2 * pi * d * v)
    return np.sqrt(re * re + im * im)


def lattice(b, alpha):
    x = (b - 1) * alpha
    return abs(x - round(x)) < 1e-12


def beta_of(b, alpha):
    x = (b - 1) * alpha
    return abs(x - round(x))


# LEMMA


def g_general(fill, x0):
    return fill - 4 * sin(pi * x0 / 2) ** 2


def kernel(b, x):
    return max(sin(pi * b * x) / sin(pi * x), 1 / sin(pi / b))


def variance(F):
    F = [float(d) for d in F]
    mean = sum(F) / len(F)
    return sum((d - mean) ** 2 for d in F) / len(F)


def g_parts(b, F, x0):
    fill = len(F)
    m = b - fill
    xf = 1 / (2 * (b - 1))
    gi = g_general(fill, x0)
    gii = m + kernel(b, x0)
    far = min(m + kernel(b, xf), g_general(fill, xf))
    giii = max(fill * sqrt(1 - 16 * x0 * x0 * variance(F)), far)
    return gi, gii, giii


def c_proved(b, F, beta):
    fill = len(F)
    gi, gii, giii = g_parts(b, F, beta / (b + 1))
    return 1 - min(gi, gii, giii) / fill, 1 - gi / fill, 1 - gii / fill, 1 - giii / fill


def harm(n):
    return log(n) + GAMMA + 1 / (2 * n)


def floor_sig(x, d):
    import mpmath

    if x <= 0:
        return mpmath.mpf(0)
    e = int(mpmath.floor(mpmath.log10(x)))
    s = mpmath.mpf(10) ** (e - d + 1)
    return mpmath.floor(x / s) * s


def margin_hp(b, chord):
    import mpmath

    mpmath.mp.dps = 40
    q = mpmath.mpf(b)
    p_ = mpmath.pi
    H = lambda n: mpmath.log(n) + mpmath.euler + 1 / (2 * n)
    if chord:
        P = mpmath.floor(q / 2)
        if b % 2 == 0:
            psi = (q / p_) * (2 * H(P - 1) - 1 + 1 / P) + (1 - 2 / p_) * q / 2
        else:
            psi = (q / p_) * (2 * H(P - 1) - 1 + 2 / P) + (1 - 2 / p_) * (
                q / 2 + 1 / (2 * q)
            )
        pb = ((4 / p_) * q + psi + q / 2 - mpmath.mpf(1) / 2) / q
    else:
        n = mpmath.ceil((q - 2) / 2)
        phi = (
            (4 / p_) * q
            + (2 * q / p_) * H(n)
            + (1 - 2 / p_) * (q - 2)
            + mpmath.mpf("0.727")
        )
        pb = 1 + phi / q
    gap = (q - 1) * q ** mpmath.mpf(-0.8) - pb
    return mpmath.log(1 + gap / pb) / mpmath.log(q), pb


def margin_iv(b):
    import mpmath
    from mpmath import iv

    iv.prec = 200
    q = iv.mpf(b)
    n = iv.mpf((b - 1) // 2)
    phi = (
        (4 / iv.pi) * q
        + (2 * q / iv.pi) * (iv.log(n) + iv.euler + 1 / (2 * n))
        + (1 - 2 / iv.pi) * (q - 2)
        + iv.mpf("0.727")
    )
    pb = 1 + phi / q
    m = iv.log(1 + ((q - 1) * q ** (iv.mpf(-4) / 5) - pb) / pb) / iv.log(q)
    return mpmath.mpf(m.a), mpmath.mpf(m.b)


def directed(x, d, up):
    import mpmath

    s = mpmath.mpf(10) ** (int(mpmath.floor(mpmath.log10(abs(x)))) - d + 1)
    return (mpmath.ceil(x / s) if up else mpmath.floor(x / s)) * s


def c_hp(b, m, beta):
    import mpmath

    mpmath.mp.dps = 40
    q = mpmath.mpf(b)
    fill = q - m
    x0 = mpmath.mpf(beta) / (q + 1)
    xf = 1 / (2 * (q - 1))
    sn = mpmath.sin
    P = mpmath.pi
    gen = lambda x: fill - 4 * sn(P * x / 2) ** 2
    ker = lambda x: max(sn(P * q * x) / sn(P * x), 1 / sn(P / q))
    e0 = b - 1
    mean = (q * (q - 1) / 2 - e0) / fill
    var = ((q - 1) * q * (2 * q - 1) / 6 - e0**2) / fill - mean**2
    g3 = max(fill * mpmath.sqrt(1 - 16 * x0**2 * var), min(m + ker(xf), gen(xf)))
    return 1 - min(gen(x0), m + ker(x0), g3) / fill


def delta_hp(b, m, beta, chord):
    import mpmath

    mpmath.mp.dps = 40
    marg, _pb = margin_hp(b, chord)
    a1 = mpmath.mpf(1) / 5 - marg
    c = c_hp(b, m, beta)
    gam = -mpmath.log(1 - c) / (2 * mpmath.log(b))
    sig = gam * (1 - 4 * a1) / (5 - 4 * a1)
    ab = mpmath.log(b - m) / mpmath.log(b)
    eta = min(marg, sig)
    return marg, c, gam, sig, ab * eta / 4


def verb_lemma():
    b = 10
    F = digits_of(b, {7})
    print(
        "two-position constant c = 1 - G(beta/(b+1))/fill at base 10 missing 7, G the least of the three forms"
    )
    print("alpha  (b-1)alpha  beta  c_pair  c_kernel  c_variance  c  gamma")
    for name, a in ALPHAS:
        if lattice(b, a):
            print(f"{name:>10}  {(b - 1) * a:.6f}  0  lattice")
            continue
        be = beta_of(b, a)
        c, c1, c2, c3 = c_proved(b, F, be)
        gam = -log(1 - c) / (2 * log(b))
        print(
            f"{name:>10}  {(b - 1) * a:.6f}  {be:.6f}  {c1:.6e}  {c2:+.6f}  {c3:+.6f}  {c:.6e}  {gam:.6e}"
        )
    print()
    print(
        "delta_0 = (alpha_b/4) min(1/5 - alpha_1, gamma (1 - 4 alpha_1)/(5 - 4 alpha_1)), each printed floored at four digits"
    )
    print("base  form  beta  1/5-alpha_1  c  gamma  sigma  delta_0")
    for b2, chord in ((92317, False), (39363, True), (10**6, False)):
        for be in (0.5, 0.2, 1e-3):
            marg, c, gam, sig, dl = delta_hp(b2, 1, be, chord)
            if marg <= 0:
                raise SystemExit(f"margin not positive at {b2}")
            form = "chord" if chord else "step3"
            print(
                f"{b2}  {form}  {be}  {floor_sig(marg, 4)}  {floor_sig(c, 4)}  {floor_sig(gam, 4)}  {floor_sig(sig, 4)}  {floor_sig(dl, 4)}"
            )
    marg, _ = margin_hp(92316, False)
    print(
        "the variance form at large base is read at the extreme excluded digit, the least variance"
    )
    for be in (0.5, 0.2):
        lo, hi = 92317, 10**6
        bind = lambda q: delta_hp(q, 1, be, False)[3] < delta_hp(q, 1, be, False)[0]
        if bind(lo) or not bind(hi):
            raise SystemExit("binding switch not bracketed")
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if bind(mid) else (mid, hi)
        print(
            f"step3, beta {be}: region A's 1/5 - alpha_1 binds at 92317 and at {lo}, sigma at {hi} and at 10^6 (bisection)"
        )
    F2 = np.array([0, 2])
    vals = [exact_at(5, F2, 1 / 8, k, 3 / 8) for k in (4, 8, 12)]
    if any(abs(v - 1) > 1e-9 for v in vals):
        raise SystemExit("sharpness example not full")
    print(
        f"sharpness: base 5, F = {{0,2}}, alpha = 1/8 (beta = {beta_of(5, 1 / 8)}): abs(hat w_k(3/8))/fill^k at k = 4, 8, 12 = {[round(v, 12) for v in vals]}"
    )
    print(f"step3 margin at 92316: {float(marg):.4e} (negative: the wall is 92317)")
    lo, _ = margin_iv(92317)
    _, hi = margin_iv(92316)
    if not (lo > 0 and hi < 0):
        raise SystemExit("step3 wall not certified at 92317")
    print(
        f"step3 wall in interval arithmetic at 200 bits: 1/5 - alpha_1 >= {directed(lo, 9, False)} at 92317, <= {directed(hi, 9, True)} at 92316"
    )


# SUP


def sup_dp(b, F, alpha, kmax, W):
    fill = len(F)
    n = b**W
    h = 1.0 / n
    v = (np.arange(n) + 0.5) * h + alpha
    lip = 2 * pi * float(F.sum())
    up = np.minimum(1.0, (gabs(F, v) + lip * h / 2) / fill)
    lu = np.log(up)
    out = []
    best = lu.reshape(b, n // b).max(0)
    args = [lu.reshape(b, n // b).argmax(0)]
    for k in range(1, kmax + 1):
        if k > 1:
            val = np.repeat(best, b) + lu
            r = val.reshape(b, n // b)
            args.append(r.argmax(0))
            best = r.max(0)
        s = int(best.argmax())
        words = []
        for j in range(k - 1, -1, -1):
            w = int(args[j][s]) * (n // b) + s
            words.append(w)
            s = w // b
        words.reverse()
        T = words[0]
        for w in words[1:]:
            T = T * b + (w % b)
        L = k + W - 1
        tl = local_sup(b, F, alpha, k, T, L)
        out.append((float(np.exp(best.max())), tl))
    return out


def local_sup(b, F, alpha, k, T, L):
    fill = len(F)
    offs = np.linspace(-1.0, 2.0, 3001)
    vals = np.ones_like(offs)
    for j in range(k):
        mod = b ** (L - j)
        v = ((T % mod) + offs) / mod + alpha
        vals *= gabs(F, v) / fill
    return float(vals.max())


def exact_at(b, F, alpha, k, t):
    fill = len(F)
    p = 1.0
    for j in range(k):
        p *= float(gabs(F, np.array([(b**j) * t + alpha]))[0]) / fill
    return p


def verb_sup():
    for b, e, kmax, W in ((10, {7}, 8, 6), (37, {18}, 4, 4)):
        F = digits_of(b, e)
        print(
            f"sup_t abs(hat w_k(t))/fill^k at base {b} missing {sorted(e)}: DP upper and attained lower, window {W}"
        )
        for name, a in ALPHAS:
            lat = lattice(b, a)
            res = sup_dp(b, F, a, kmax, W)
            if lat:
                at = exact_at(b, F, a, kmax, (-a) % 1.0)
                if abs(at - 1) > 1e-9:
                    raise SystemExit(f"lattice alpha {name} not full at t = -alpha")
                row = " ".join(f"{lo:.6f}" for _, lo in res)
                print(
                    f"  {name:>10} lattice   attained {row}   value at t=-alpha {at:.12f}"
                )
                continue
            be = beta_of(b, a)
            c = c_proved(b, F, be)[0]
            row = []
            for k, (hi, lo) in enumerate(res, 1):
                pk = (1 - c) ** (k // 2)
                if lo > pk * (1 + 1e-12):
                    raise SystemExit(
                        f"REFUTED: alpha {name} k {k} attained {lo} above proved {pk}"
                    )
                if k >= 2 and hi > pk:
                    raise SystemExit(f"upper bracket above proved: alpha {name} k {k}")
                row.append(f"{lo:.6f}/{hi:.6f}/{pk:.6f}")
            rate = res[-1][0] ** (1 / kmax)
            print(
                f"  {name:>10} beta {be:.4f}  attained/upper/proved: {' '.join(row)}  per-digit upper {rate:.6f} proved {(1 - c) ** 0.5:.6f}"
            )
        print()


# PEEL


def grid_sum(b, F, alpha, i, s):
    a = np.arange(b**i, dtype=np.int64)
    prod = np.ones(b**i)
    for j in range(i):
        frac = ((b**j) * s) % 1.0
        v = frac + ((a * (b**j)) % (b**i)) / (b**i) + alpha
        prod *= gabs(F, v)
    return prod


def gsum(b, F, T):
    return sum(gabs(F, (T + r) / b) for r in range(b))


def verb_peel():
    b, e = 10, {7}
    F = digits_of(b, e)
    N = 200000
    Tg = np.arange(N) / N
    G = gsum(b, F, Tg)
    lip = 2 * pi * float(F.sum())
    bup = float(G.max()) + lip / (2 * N)
    phi = (
        (4 / pi) * b
        + (2 * b / pi) * harm(-(-(b - 2) // 2))
        + (1 - 2 / pi) * (b - 2)
        + 0.727
    )
    print(
        f"B_10(F) at missing 7: grid max {G.max():.6f}, certified upper {np.ceil(bup * 1e6) / 1e6:.6f}, Parseval floor 2(b-1) = {2 * (b - 1)}, step 3 bound {b * (1 + phi / b):.4f}"
    )
    rng = np.random.default_rng(20250)
    worst_ratio = 0.0
    worst_step = 0.0
    worst_id = 0.0
    nchk = 0
    for name, a in ALPHAS:
        for _ in range(12):
            s = float(rng.random())
            for i in range(1, 6):
                cur = grid_sum(b, F, a, i, s)
                Si = float(cur.sum())
                prev = grid_sum(b, F, a, i - 1, (b * s) % 1.0) if i > 1 else np.ones(1)
                ap = np.arange(b ** (i - 1))
                T = (b * (s + a)) % 1.0 + ap / (b ** (i - 1))
                step = float((prev * gsum(b, F, T)).sum())
                worst_id = max(worst_id, abs(step - Si) / Si)
                worst_ratio = max(worst_ratio, Si / bup**i)
                worst_step = max(worst_step, Si / (bup * float(prev.sum())))
                nchk += 1
    if worst_ratio > 1 or worst_step > 1 or worst_id > 1e-9:
        raise SystemExit("REFUTED: twisted peel")
    print(f"{nchk} checks over {len(ALPHAS)} alpha classes, 12 shifts, depths 1..5")
    print(f"  one-step identity, largest relative error {worst_id:.2e}")
    print(
        f"  Sigma_i(s)/B^i at most {worst_ratio:.6f};  Sigma_i(s)/(B Sigma_(i-1)(b s)) at most {worst_step:.6f}"
    )
    for name, a in (("0", 0.0), ("1/2", 0.5), ("1/7", 1 / 7)):
        m0 = float(grid_sum(b, F, a, 5, 0.0).sum())
        print(
            f"  alpha {name}: Sigma_5(0) = {m0:.1f}, Sigma_5(0)^(1/5) = {m0**0.2:.4f}"
        )


# METER


def sieve(N):
    isp = np.ones(N + 1, dtype=bool)
    isp[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if isp[p]:
            isp[p * p :: p] = False
    primes = np.nonzero(isp)[0]
    mu = np.ones(N + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes:
        mu[p::p] *= -1
        if p * p <= N:
            mu[p * p :: p * p] = 0
    lam = np.zeros(N + 1)
    for p in primes:
        q = int(p)
        lp = log(q)
        while q <= N:
            lam[q] = lp
            q *= int(p)
    return isp, mu, lam


def design_mask(N, b, e):
    n = np.arange(N + 1, dtype=np.int64)
    ok = n > 0
    s = np.zeros(N + 1, dtype=np.int64)
    x = n.copy()
    while True:
        live = x > 0
        if not live.any():
            break
        d = x % b
        for ex in e:
            ok &= ~(live & (d == ex))
        s += d
        x //= b
    return ok, s


CACHE = {}


def arrays(N):
    if N not in CACHE:
        CACHE[N] = sieve(N)
    return CACHE[N]


def verb_meter():
    N = 10**7
    _isp, mu, lam = arrays(N)
    for e in ({7}, {0}):
        ok, s = design_mask(N, 10, e)
        idx = np.nonzero(ok)[0]
        A = np.cumsum(ok)
        fill = 10 - len(e)
        kappa = 2.5 * len([f for f in (1, 3, 7, 9) if f not in e]) / fill
        print(f"base 10 missing {sorted(e)}: A_F(10^7) = {A[N]}, kappa_F = {kappa:.6f}")
        print(
            "  alpha  class  abs(M)/A at 10^4..10^7  max abs(M)/A at 10^7  psi^alpha/psi^0 at 10^5..10^7  limit  model at 10^7"
        )
        psi0 = np.cumsum(lam[idx])
        cop = np.gcd(idx, 30) == 1
        top = [0.0, 0.0, 0.0, 0.0]
        for name, a in ALPHAS:
            ph = np.exp(2j * pi * a * s[idx])
            M = np.cumsum(mu[idx] * ph)
            P = np.cumsum(lam[idx] * ph)
            cols = []
            for j in range(4, 8):
                x = 10**j
                pos = np.searchsorted(idx, x, side="right") - 1
                cols.append(f"{abs(M[pos]) / A[x]:.5f}")
            mx = float(np.abs(M).max()) / A[N]
            pc = []
            for j in range(5, 8):
                pos = np.searchsorted(idx, 10**j, side="right") - 1
                pc.append(
                    f"{(P[pos] / psi0[pos]).real:+.4f}{(P[pos] / psi0[pos]).imag:+.4f}i"
                )
            if lattice(10, a):
                r = round(9 * a) % 9
                q = 9 // gcd(r, 9)
                pred = {1: 1.0, 3: -0.5, 9: 0.0}[q]
                cls = "lattice"
            else:
                pred = 0.0
                cls = "off"
            mod = ph[cop].mean() if cls == "off" else pred
            ratio = P[-1] / psi0[-1]
            top[0] = max(top[0], abs(M[-1]) / A[N])
            top[1] = max(top[1], mx)
            k = 2 if cls == "off" else 3
            top[k] = max(top[k], abs(ratio - mod))
            print(
                f"  {name:>10} {cls:>7}  {' '.join(cols)}  {mx:.5f}  {' '.join(pc)}  {pred:+.2f}  {mod.real:+.4f}{mod.imag:+.4f}i"
            )
        print(
            f"  over all classes at 10^7: abs(M)/A at most {top[0]:.5f}, running max at most {top[1]:.5f}"
        )
        print(
            f"  psi^alpha/psi^0 at 10^7 against its model: off the lattice within {top[2]:.4f}, on it within {top[3]:.4f} of the Ramanujan value"
        )
        print()


def verb_gelfond():
    N = 10**7
    isp, _mu, _lam = arrays(N)
    for e in ({7}, {0}):
        ok, s = design_mask(N, 10, e)
        sel = ok & isp
        sp = s[sel]
        tot = int(sel.sum())
        n = np.arange(N + 1)
        sc = s[ok & (np.gcd(n, 30) == 1)]
        print(
            f"primes of base 10 missing {sorted(e)} below 10^7: {tot}; model: the elements of the design coprime to b (b - 1) = 90"
        )
        for m in (2, 4, 5, 7, 8, 10, 3, 9):
            cnt = np.bincount(sp % m, minlength=m)
            mc = np.bincount(sc % m, minlength=m) / len(sc)
            dev = float(np.abs(cnt / (tot / m) - 1).max())
            mdev = float(np.abs(mc * m - 1).max())
            gap = float(np.abs(cnt / tot - mc).max() * m)
            tag = "coprime to 9" if gcd(m, 9) == 1 else "shares a factor with 9"
            print(
                f"  m {m:>2} ({tag}): counts {cnt.tolist()}  deviation from 1/m: primes {dev:.5f}, model {mdev:.5f}; primes against model {gap:.5f}"
            )
        print()


# LATTICE


def odd_primes(N):
    h = (N + 1) // 2
    odd = np.ones(h, dtype=bool)
    odd[0] = False
    for i in range(1, (int(N**0.5) - 1) // 2 + 1):
        if odd[i]:
            p = 2 * i + 1
            odd[p * p // 2 :: p] = False
    return np.concatenate(([2], 2 * np.nonzero(odd)[0] + 1)).astype(np.int64)


def digits_ok(v, b, e):
    ok = np.ones(len(v), dtype=bool)
    s = np.zeros(len(v), dtype=np.int64)
    x = v.copy()
    while (x > 0).any():
        live = x > 0
        d = x % b
        for ex in e:
            ok &= ~(live & (d == ex))
        s += d
        x //= b
    return ok, s


def sieved_model(b, F, k, P):
    M = P * 9 // gcd(P, 9)
    cnt = np.zeros((M, 2), dtype=np.int64)
    cnt[0, 0] = 1
    acc = np.zeros_like(cnt)
    r = np.arange(M)
    for _ in range(k):
        new = np.zeros((M, 2), dtype=np.int64)
        for d in F:
            to = (b * r + d) % M
            np.add.at(new[:, (d % 2)], to, cnt[:, 0])
            np.add.at(new[:, 1 - (d % 2)], to, cnt[:, 1])
        cnt = new
        if 0 not in F:
            acc += cnt
    if 0 not in F:
        cnt = acc
    keep = np.gcd(r, P) == 1
    return {
        (int(n), par): int(cnt[n, par])
        for n in np.nonzero(keep)[0]
        for par in (0, 1)
        if cnt[n, par]
    }


def tests(u9):
    return [
        (
            "p mod 9, law 1/6 on the units",
            9,
            lambda n, s: n % 9,
            lambda r: 1 / 6 if r in u9 else 0.0,
        ),
        (
            "s mod 3, law 1/2 on 1, 2",
            3,
            lambda n, s: s % 3,
            lambda r: 0.5 if r % 3 else 0.0,
        ),
        ("s mod 2, law 1/2", 2, lambda n, s: s % 2, lambda r: 0.5),
        (
            "s mod 6, law 1/4 on r prime to 3",
            6,
            lambda n, s: s % 6,
            lambda r: 0.25 if r % 3 else 0.0,
        ),
        (
            "s mod 18, law 1/12 on r prime to 9",
            18,
            lambda n, s: s % 18,
            lambda r: 1 / 12 if gcd(r, 9) == 1 else 0.0,
        ),
        (
            "(p mod 9, s mod 6), law 1/12 where s = p mod 3",
            54,
            lambda n, s: (n % 9) * 6 + s % 6,
            lambda r: 1 / 12 if r // 6 in u9 and (r % 6) % 3 == (r // 6) % 3 else 0.0,
        ),
        (
            "(p mod 7, s mod 3), law 1/12 on the units mod 7 and s in 1, 2",
            21,
            lambda n, s: (n % 7) * 3 + s % 3,
            lambda r: 1 / 12 if r // 3 and r % 3 else 0.0,
        ),
    ]


def verb_lattice():
    N = 10**8
    pr = odd_primes(N)
    u9 = [r for r in range(9) if gcd(r, 9) == 1]
    P = 2 * 3 * 5 * 7 * 11 * 13
    for b, e in ((10, {7}), (10, {0})):
        F = [d for d in range(b) if d not in e]
        ok, s = digits_ok(pr, b, e)
        p0, s0 = pr[ok], s[ok]
        model = sieved_model(b, F, 8, P)
        mn = np.array([k[0] for k in model])
        mpar = np.array([k[1] for k in model])
        mpar = mn % 9 + 9 * ((mpar - mn % 9) % 2)
        mw = np.array(list(model.values()), dtype=float)
        print(
            f"base {b} missing {sorted(e)}: each class as (count - law)/sqrt(law) against the limit law, then against the design's 8-digit strings prime to {P}, read exactly by a digit recursion"
        )
        worst_law = worst_model = 0.0
        for X in (10**6, 10**7, 10**8):
            sel = (p0 > 3 * b) & (p0 <= X)
            p, sp = p0[sel], s0[sel]
            tot = len(p)
            print(f"  below 10^{round(log10(X))}: {tot} primes of S_F above {3 * b}")
            for name, size, key, law in tests(u9):
                kp = key(p, sp)
                cells = list(range(size))
                cp = {c: int((kp == c).sum()) for c in cells}
                zs, zm = [], []
                if X == N:
                    km = key(mn, mpar)
                    mt = mw.sum()
                for c in cells:
                    if law(c) == 0:
                        if cp[c]:
                            raise SystemExit("a prime in a class of law 0")
                        continue
                    pred = law(c) * tot
                    zs.append((cp[c] - pred) / sqrt(pred))
                    if X == N:
                        mp = mw[km == c].sum() / mt * tot
                        zm.append((cp[c] - mp) / sqrt(mp))
                wl = max(abs(z) for z in zs)
                line = f"    {name}: worst against the law {wl:.2f}"
                if X == N:
                    worst_law = max(worst_law, wl)
                    wm = max(abs(z) for z in zm)
                    worst_model = max(worst_model, wm)
                    line += f", against the sieved model {wm:.2f}"
                print(line)
            even = int((sp % 2 == 0).sum())
            print(f"    share of even digit sums: {even / tot:.5f}")
        print(
            f"  at 10^8: worst against the law {worst_law:.2f}, worst against the sieved model {worst_model:.2f}"
        )
        for Q in (30, 210, 330, 390, 2310, P):
            mq = sieved_model(b, F, 8, Q)
            ev = sum(v for (n, par), v in mq.items() if par == 0) / sum(mq.values())
            print(
                f"  share of even digit sums on the 8-digit strings prime to {Q}: {ev:.5f}"
            )
        ratio = abs(sum(np.exp(2j * pi * d / (2 * (b + 1))) for d in F)) / len(F)
        print(f"  abs(g_F(1/(2(b + 1))))/fill = {ratio:.4f}")
        bias = {}
        for k in (6, 7, 8, 10, 12):
            sh = []
            for Q in (30, 330):
                mq = sieved_model(b, F, k, Q)
                sh.append(
                    sum(v for (n, par), v in mq.items() if par == 0) / sum(mq.values())
                )
            bias[k] = sh[1] - sh[0]
        print(
            "  bias from 11, even share prime to 330 less prime to 30, at k digits: "
            + ", ".join(f"k {k}: {v:+.5f}" for k, v in bias.items())
        )
        print(
            "  its ratio over two digits at even k: "
            + ", ".join(
                f"{k} to {k + 2}: {bias[k + 2] / bias[k]:.2f}" for k in (6, 8, 10)
            )
            + f", against abs(g_F(1/22))^2/fill^2 = {ratio**2:.3f}"
        )
        if worst_model > 4:
            raise SystemExit("a class drifts from the sieved model")
        print()


ALPHAS = [
    ("0", 0.0),
    ("1/9", 1 / 9),
    ("1/3", 1 / 3),
    ("2/3", 2 / 3),
    ("4/9", 4 / 9),
    ("1/2", 1 / 2),
    ("1/4", 1 / 4),
    ("1/5", 1 / 5),
    ("1/7", 1 / 7),
    ("1/10", 1 / 10),
    ("1/11", 1 / 11),
    ("1/27", 1 / 27),
    ("1/90", 1 / 90),
    ("1/900", 1 / 900),
    ("sqrt2-1", sqrt(2) - 1),
    ("golden", (sqrt(5) - 1) / 2),
]

VERBS = {
    "lemma": verb_lemma,
    "sup": verb_sup,
    "peel": verb_peel,
    "meter": verb_meter,
    "gelfond": verb_gelfond,
    "lattice": verb_lattice,
}

if __name__ == "__main__":
    for v in sys.argv[1:] or list(VERBS):
        t0 = time.time()
        VERBS[v]()
        print(f"[{v}: {time.time() - t0:.1f} s]")
        print()
