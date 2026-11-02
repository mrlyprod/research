import math
import sys
import time

import numpy as np
from mpmath import iv

G0 = 0.9625229
C1F = 2 / math.pi
GATES = ((1, 5), (1, 6), (1, 8))

# FLOAT TRANSFORM


def kern(b, t):
    t = np.asarray(t, dtype=float)
    n = np.round(t)
    th = t - n
    s = np.sin(np.pi * th)
    small = np.abs(s) < 1e-13
    sign = np.where(((b - 1) * n.astype(np.int64)) % 2 == 0, 1.0, -1.0)
    return sign * np.where(
        small, float(b), np.sin(b * np.pi * th) / np.where(small, 1.0, s)
    )


def classes(b):
    h = [(2 * a - (b - 1)) / 2 for a in range(b)]
    mags = sorted({abs(x) for x in h})
    return [(m1, m2) for i, m1 in enumerate(mags) for m2 in mags[: i + 1]]


def rep(b, m1, m2):
    return (round(m1 + (b - 1) / 2), round(m2 + (b - 1) / 2))


def fgrid(b, m1, m2, t1, t2):
    K1, K2 = kern(b, t1), kern(b, t2)
    e1, e2 = np.exp(2j * np.pi * m1 * t1), np.exp(2j * np.pi * m2 * t2)
    return np.abs(np.outer(K1, K2) - np.outer(e1, e2)) / (b * b - 1)


def fdirect(b, a, t1, t2):
    acc = np.zeros((len(t1), len(t2)), dtype=complex)
    for x in range(b):
        for y in range(b):
            if (x, y) != a:
                acc += np.outer(
                    np.exp(2j * np.pi * x * t1), np.exp(2j * np.pi * y * t2)
                )
    return np.abs(acc) / (b * b - 1)


# CENSUS


def sigma1(b, m1, m2, G):
    s = np.arange(G) / (G * b)
    i = np.arange(b) / b
    t = (s[:, None] + i[None, :]).ravel()
    F = fgrid(b, m1, m2, t, t).reshape(G, b, G, b)
    S = F.sum(axis=(1, 3))
    return float(S.min()), float(S.max()), float(S[0, 0])


def perron(M, iters=4000, tol=1e-12):
    y = np.ones(M.shape[0])
    lam = 0.0
    for _ in range(iters):
        z = M @ y
        nl = float(z.max())
        if nl <= 0.0:
            return 0.0, y
        z /= nl
        if abs(nl - lam) <= tol * nl and np.abs(z - y).max() < 1e-10:
            return nl, z
        y, lam = z, nl
    return lam, y


def window2(b, m1, m2, r):
    n = b * b * r
    t = np.arange(n + 1) / n
    F = fgrid(b, m1, m2, t, t)
    hi = np.full((b * b, b * b), 0.0)
    lo = np.full((b * b, b * b), np.inf)
    for u in range(r + 1):
        for v in range(r + 1):
            blk = F[u : u + b * b * r : r, v : v + b * b * r : r]
            hi = np.maximum(hi, blk)
            lo = np.minimum(lo, blk)

    def mat(C):
        return C.reshape(b, b, b, b).transpose(0, 2, 1, 3).reshape(b * b, b * b)

    return perron(mat(hi))[0], perron(mat(lo))[0]


def sig_at(b, m1, m2, S):
    i = np.arange(b) / b
    t1 = S[:, 0, None] + i[None, :]
    t2 = S[:, 1, None] + i[None, :]
    P = kern(b, t1)[:, :, None] * kern(b, t2)[:, None, :]
    E = (
        np.exp(2j * np.pi * m1 * t1)[:, :, None]
        * np.exp(2j * np.pi * m2 * t2)[:, None, :]
    )
    return np.abs(P - E).sum(axis=(1, 2)) / (b * b - 1)


def local_min(b, m1, m2, n=13, width=2.5):
    R = width / b**2
    u = np.linspace(-R, R, n) / b
    S = np.stack(np.meshgrid(u, u, indexing="ij"), axis=-1).reshape(-1, 2)
    v = sig_at(b, m1, m2, S)
    k = int(v.argmin())
    best, x = float(v[k]), S[k].copy()
    step = 2 * R / (n - 1) / b
    while step > 1e-7 / b**3:
        moved = False
        for d in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)):
            y = x + step * np.array(d)
            w = float(sig_at(b, m1, m2, y[None, :])[0])
            if w < best:
                best, x, moved = w, y, True
        if not moved:
            step /= 2
    return best, x


def floor1():
    t0 = time.time()
    lo, hi = [int(x) for x in sys.argv[2:4]] if len(sys.argv) > 3 else [3, 31]
    print(
        f"THE ONE-DIGIT FLOOR: min over shifts of Sigma_1 near the grid, a 13 x 13 scan of |s_i| <= 2.5/b^3 refined by a shrinking pattern search, every orbit of every base {lo}..{hi}, float"
    )
    w = float(sig_at(7, 0.0, 0.0, np.array([[0.019396 / 7, 0.980604 / 7]]))[0])
    assert abs(w - 1.997526) < 1e-6
    print(f"  base 7 missing (3, 3) at s = (0.019396/7, 0.980604/7): Sigma_1 = {w:.6f}")
    below, flat = [], []
    for b in range(lo, hi + 1):
        for m1, m2 in classes(b):
            v, x = local_min(b, m1, m2)
            (below if v < 2.0 - 1e-9 else flat).append(
                (b, rep(b, m1, m2), (m1, m2), v, x)
            )
    for b, a, m, v, x in below:
        print(
            f"  base {b} missing {a}, m = ({m[0]:g}, {m[1]:g}): min Sigma_1 {v:.6f} at b s = ({b * x[0]:.6f}, {b * x[1]:.6f})"
        )
    zero = all(0 in r[2] for r in below)
    odd = all(r[0] % 2 for r in below)
    every = all(not (r[0] % 2 and 0 in r[2]) for r in flat)
    least = min(below, key=lambda r: r[3]) if below else None
    print(
        f"  {len(below)} orbits dip below 2 and {len(flat)} stay at 2; every dipping orbit has a zero coordinate of m: {zero}; all at odd bases: {odd}; every odd-base orbit with a zero coordinate dips: {every}"
    )
    if least:
        print(f"  least {least[3]:.6f} at base {least[0]} missing {least[1]}")
    print(f"runtime {time.time() - t0:.1f} s")


def census():
    t0 = time.time()
    rng = [int(x) for x in sys.argv[2:4]] if len(sys.argv) > 3 else [3, 44]
    print(
        "THE CENSUS: base b on both coordinates, F = {0..b-1}^2 minus one vector a, classes m = a - ((b-1)/2, (b-1)/2) up to sign and swap"
    )
    print(
        "  float readings, not certificates: Sigma_1 over a shift grid, and the two-digit window matrix with cell sup and inf read on a sub-grid"
    )
    t = (np.arange(7) + 0.37) / 7.3
    for b in (3, 4, 5):
        for a in [(x, y) for x in range(b) for y in range(b)]:
            m = (a[0] - (b - 1) / 2, a[1] - (b - 1) / 2)
            dev = np.abs(fgrid(b, m[0], m[1], t, t) - fdirect(b, a, t, t)).max()
            assert dev < 1e-12, (b, a, dev)
    print(
        "  the closed form |K(t1) K(t2) - e(m . t)|/(b^2 - 1) matches the digit sum at every vector of bases 3, 4, 5"
    )
    summary = {g: [None, None] for g in GATES}
    for b in range(rng[0], rng[1] + 1):
        r = 2 if b <= 20 else 1
        G = 24 if b <= 20 else 12
        rows = []
        for m1, m2 in classes(b):
            smin, smax, s0 = sigma1(b, m1, m2, G)
            assert abs(s0 - 2.0) < 1e-9
            wh, wl = window2(b, m1, m2, r)
            L2 = math.log(b * b)
            rows.append(
                (
                    m1,
                    m2,
                    math.log(smin) / L2,
                    math.log(smax) / L2,
                    math.log(wl) / L2 if wl > 0 else float("nan"),
                    math.log(wh) / L2,
                )
            )
        hi_all = max(x[5] for x in rows)
        lo_min = min(x[5] for x in rows)
        for g in GATES:
            gv = g[0] / g[1]
            if summary[g][0] is None and lo_min < gv:
                summary[g][0] = b
            if hi_all < gv:
                if summary[g][1] is None:
                    summary[g][1] = b
            else:
                summary[g][1] = None
        print(
            f"base {b}: {len(rows)} classes, window upper reading from {lo_min:.4f} to {hi_all:.4f}"
        )
        for m1, m2, fl, fh, wl, wh in rows:
            a = rep(b, m1, m2)
            print(
                f"  a={a} m=({m1:g},{m2:g}) shift floor N=1 {fl:.4f}  max Sigma_1 {fh:.4f}  window [{wl:.4f}, {wh:.4f}]"
            )
    for g in GATES:
        print(
            f"gate {g[0]}/{g[1]}: first base with some class read below {summary[g][0]}, every class read below from {summary[g][1]} in the range"
        )
    print(f"runtime {time.time() - t0:.1f} s")


# THE SQUARE CHAIN


def ivc():
    c1 = 2 / iv.pi
    g = iv.mpf(G0)
    return c1, g


def phi_iv(w, Lhi, blo):
    c1, g = ivc()
    A = (
        c1 * c1 * Lhi * Lhi * w * (w + 1) / (w - 1) ** 3
        + 2 * c1 * Lhi * g * w / (w - 1) ** 2
        + g * g / (w - 1)
    )
    T = (
        2 * c1 * c1 * iv.log(blo) * blo * w / (blo * w - 1) ** 2
        + 2 * c1 * g / (blo * w - 1)
        + c1 * c1 / (blo * blo * w - 1)
    )
    return A + T - (w - 1)


def target(b, p, q):
    b = iv.mpf(b)
    return b ** (iv.mpf(2 * p) / q) * (1 - 1 / (b * b))


def phi_f(z, b):
    L = math.log(b)
    A = (
        C1F * C1F * L * L * z * (z + 1) / (z - 1) ** 3
        + 2 * C1F * L * G0 * z / (z - 1) ** 2
        + G0 * G0 / (z - 1)
    )
    T = (
        2 * C1F * C1F * L * b * z / (b * z - 1) ** 2
        + 2 * C1F * G0 / (b * z - 1)
        + C1F * C1F / (b * b * z - 1)
    )
    return A + T - (z - 1)


def phi1_f(z, B):
    L = math.log(B)
    return C1F * L * z / (z - 1) ** 2 + G0 / (z - 1) + C1F / (B * z - 1) - (z - 1)


def root(phi, b):
    lo, hi = 1.0 + 1e-12, 1e6
    for _ in range(200):
        mid = (lo + hi) / 2
        if phi(mid, b) > 0:
            lo = mid
        else:
            hi = mid
    return hi


def clears_f(b, p, q):
    return phi_f(b ** (2 * p / q) * (1 - 1 / b**2), b) < 0


def cap_iv(b):
    c1, g = ivc()
    L = iv.log(b)
    tau = (4 * c1 * c1 * L + 2 * c1 * g) / b + c1 * c1 / (b * b)
    return 1 + iv.sqrt(3 * c1 * c1 * L * L + 3 * c1 * L * g + g * g + tau)


def slope_iv(b, p, q):
    c1, g = ivc()
    d = iv.mpf(2 * p) / q
    return d * target(b, p, q) - iv.sqrt(3) * c1 - iv.sqrt(3) * g / (2 * iv.log(b))


def block_ok(lo, hi, p, q):
    return float(phi_iv(target(lo, p, q), iv.log(hi), iv.mpf(lo)).b) < 0.0


def certify(lo, top, p, q):
    blocks = 0
    b = lo
    step = 1
    while b <= top:
        hi = min(b + step - 1, top)
        if block_ok(b, hi, p, q):
            blocks += 1
            b = hi + 1
            step *= 2
        elif step > 1:
            step //= 2
        else:
            raise SystemExit(f"block at base {b} does not clear for gate {p}/{q}")
    return blocks


def chain():
    t0 = time.time()
    iv.prec = 120
    gam = (2 / iv.pi) * (iv.euler + iv.log(8 / iv.pi))
    assert float(gam.b) <= G0
    print(
        "THE SQUARE CHAIN: |hat F| <= |D_b(t1)| |D_b(t2)| + 1 at every missing vector; runs telescope in each coordinate, so lambda_l enters squared"
    )
    print(
        "  z_c - 1 = sum_l (lambda_l^+)^2 z_c^(-l), lambda_l^+ = (2/pi) l log b + gamma' + (2/pi) b^(-l); alpha_2 <= log_(b^2)(z_c b^2/(b^2 - 1)) in units of b^(2k)"
    )
    for p, q in GATES:
        w0 = next(
            b
            for b in range(3, 10**6)
            if all(clears_f(r, p, q) for r in range(b, b + 2000))
        )
        assert (
            float(phi_iv(target(w0 - 1, p, q), iv.log(w0 - 1), iv.mpf(w0 - 1)).a) > 0.0
        )
        bc = w0
        while not (
            float((target(bc, p, q) - cap_iv(bc)).a) > 0.0
            and float((target(bc, p, q) - 3).a) > 0.0
            and float(slope_iv(bc, p, q).a) > 0.0
        ):
            bc = int(bc * 1.05) + 1
        nb = certify(w0, bc, p, q)
        z = root(phi_f, w0)
        al = (
            math.ceil(math.log(z * w0 * w0 / (w0 * w0 - 1)) / math.log(w0 * w0) * 1e6)
            / 1e6
        )
        zf = root(phi_f, w0 - 1)
        af = (
            math.floor(
                math.log(zf * (w0 - 1) ** 2 / ((w0 - 1) ** 2 - 1))
                / math.log((w0 - 1) ** 2)
                * 1e6
            )
            / 1e6
        )
        print(
            f"  gate {p}/{q}: wall {w0}, certified at 120 bits on [{w0}, {bc}] in {nb} interval blocks; cap max(3, 1 + sqrt(3 c_pi^2 L_b^2 + 3 c_pi L_b gamma' + gamma'^2 + tau_b)) below the target from {bc} with the slope test positive there"
        )
        print(
            f"    z_c = {z:.6f} and alpha_2 < {al:.6f} at {w0}; at {w0 - 1} the root equation fails at the target, alpha_2 bound {af:.6f} from the float root"
        )
    print(
        "THE BOX AGAINST THE MORTON LINE: the same set read at base b^2 in one dimension, chain (z_M-1)^3 = (2/pi) log(b^2) z_M + gamma'(z_M-1) + (2/pi)(z_M-1)^2/(b^2 z_M - 1)"
    )
    for e in range(1, 7):
        b = 10**e
        z2 = root(phi_f, b)
        z1 = root(phi1_f, b * b)
        a2 = math.log(z2 * b * b / (b * b - 1)) / math.log(b * b)
        a1 = math.log(z1 * b * b / (b * b - 1)) / math.log(b * b)
        print(
            f"  b = 10^{e}: box z_c = {z2:.4f} alpha {a2:.5f}, Morton z_M = {z1:.4f} alpha {a1:.5f}, ratio {a2 / a1:.4f}, z_c / ((2/pi) log b) = {z2 / (C1F * math.log(b)):.4f}"
        )
    print(f"runtime {time.time() - t0:.1f} s")


# THE SQUARE UNIFORM WINDOW

MARGIN = 1e-12


def ksup(b, r):
    n = b * b * r
    br = b * r
    j = np.arange(n + 1)
    num = np.abs(np.sin(np.pi * ((b * j) % (2 * n)) / n)) + MARGIN
    den = np.sin(np.pi * j / n) - MARGIN
    jl = j[:-1]
    o = -(-(2 * jl) // br)
    o = o + (o % 2 == 0)
    hi = np.where(o * br <= 2 * jl + 2, 1.0 + MARGIN, np.maximum(num[:-1], num[1:]))
    lo = np.minimum(den[:-1], den[1:])
    s = np.full(n, float(b))
    ok = lo > 0
    s[ok] = hi[ok] / lo[ok]
    left = (jl + 1 <= br) & (jl >= 1)
    s[left] = num[:-1][left] / den[:-1][left]
    right = (jl >= n - br) & (jl + 1 <= n - 1)
    s[right] = num[1:][right] / den[1:][right]
    s[(jl == 0) | (jl == n - 1)] = float(b)
    return np.minimum(s, float(b)).reshape(b * b, r).max(axis=1)


def uwin_rho(b, r):
    A = ksup(b, r).reshape(b, b)
    Y = np.ones((b, b))
    lam = 0.0
    for _ in range(3000):
        Z = A @ Y @ A.T + Y.sum()
        nl = float(Z.max())
        Z /= nl
        if abs(nl - lam) <= 1e-13 * nl and np.abs(Z - Y).max() < 1e-11:
            Y = Z
            break
        Y, lam = Z, nl
    Y = np.maximum(Y, 1e-300)
    ratio = (A @ Y @ A.T + Y.sum()) / Y
    return (
        float(ratio.max()) * (1 + 1e-9),
        float(ratio.min()) * (1 - 1e-9),
        float(Y.max() / Y.min()),
    )


def uwin_alpha(b, r):
    rho, _rlo, spread = uwin_rho(b, r)
    return math.log(rho / (b * b - 1)) / math.log(b * b), spread


def uwin():
    t0 = time.time()
    from mpmath import mp, mpf
    from mpmath import pi as mpi
    from mpmath import sin as msin

    mp.prec = 80
    rs = np.random.default_rng(7)
    worst = 0.0
    for _ in range(400):
        b = int(rs.integers(3, 3000))
        n = b * b * 2
        jj = int(rs.integers(1, n))
        for x, y in (
            (np.sin(np.pi * ((b * jj) % (2 * n)) / n), msin(mpi * b * mpf(jj) / n)),
            (np.sin(np.pi * jj / n), msin(mpi * mpf(jj) / n)),
        ):
            worst = max(worst, abs(float(mpf(float(x)) - y)))
    assert worst < MARGIN / 100
    print(
        "THE SQUARE UNIFORM WINDOW: Psi(t) = (|D_b(t1)| |D_b(t2)| + 1)/(b^2 - 1) dominates |hat F|/fill at every missing vector, its cell sup is (Delta[k1] Delta[k2] + 1)/(b^2 - 1)"
    )
    print(
        f"  Delta[k] a 1D sup of |D_b| on the cell [k/b^2, (k+1)/b^2], read on r sub-cells from endpoint sines with margin {MARGIN:g} (float sines within {worst:.1e} of 80-bit values at 800 sampled points)"
    )
    print(
        "  the window matrix at two digit-vectors is (A (x) A + 1 1^T)/(b^2 - 1), A[c, c'] = Delta[c b + c'], and its Collatz-Wielandt bound gives alpha_2 <= log_(b^2) of that bound at every vector at once"
    )
    for b in (6, 8, 10):
        for a in [(0, 0), (b // 2, b // 3), (b - 1, 1)]:
            t = (np.arange(b * b * 2) + 0.5) / (b * b * 2)
            assert np.all(
                fdirect(b, a, t[::7], t[::5])
                <= (
                    (np.outer(np.abs(kern(b, t[::7])), np.abs(kern(b, t[::5]))) + 1)
                    / (b * b - 1)
                )
                + 1e-12
            )
    print(
        "  the domination |hat F|/fill <= u is checked on a grid at three vectors of bases 6, 8, 10"
    )
    walls = {}
    for p, q in GATES[:2]:
        z = next(
            b
            for b in range(3, 2000)
            if clears_f(b, p, q) and all(clears_f(r, p, q) for r in range(b, b + 500))
        )
        b = z - 1
        rows = {}
        while True:
            a, sp = uwin_alpha(b, 2)
            rows[b] = (a, sp)
            if a >= p / q:
                break
            b -= 1
        lo = b + 1
        walls[(p, q)] = lo
        top = max(range(lo, z), key=lambda x: rows[x][0])
        print(
            f"  gate {p}/{q}: alpha < {p}/{q} at every base {lo}..{z - 1}, largest {math.ceil(rows[top][0] * 1e6) / 1e6:.6f} at {top}; base {b} reads {math.floor(rows[b][0] * 1e6) / 1e6:.6f} and fails; the square chain holds from {z}"
        )
        print(
            f"    alpha < {math.ceil(rows[lo][0] * 1e6) / 1e6:.6f} at {lo}, Perron vector spread {rows[lo][1]:.2f}"
        )
    for b in (300, 600, 1000, 2000, 3000):
        a, sp = uwin_alpha(b, 1)
        print(
            f"  base {b}: square uniform window alpha < {math.ceil(a * 1e6) / 1e6:.6f}, square chain {math.log(root(phi_f, b) * b * b / (b * b - 1)) / math.log(b * b):.6f}"
        )
    print(f"runtime {time.time() - t0:.1f} s")


def window1(B, a0, r):
    n = B * B * r
    t = np.arange(n + 1) / n
    g = np.abs(kern(B, t) - np.exp(2j * np.pi * (a0 - (B - 1) / 2) * t)) / (B - 1)
    hi = np.zeros(B * B)
    for u in range(r + 1):
        hi = np.maximum(hi, g[u : u + B * B * r : r])
    return perron(hi.reshape(B, B))[0]


def morton():
    t0 = time.time()
    print(
        "THE SAME SET READ TWICE: the box reading at base b on Z[i] and the Morton reading at base b^2 on Z, both two-digit window matrices with sampled cell sups, both in units of b^(2k)"
    )
    for b in (8, 12, 16, 24, 32):
        r = 2 if b <= 16 else 1
        for a in [(0, 0), (b // 2, b // 2), (b - 1, 0), (b // 3, b - 2)]:
            m1, m2 = a[0] - (b - 1) / 2, a[1] - (b - 1) / 2
            w2 = math.log(window2(b, m1, m2, r)[0]) / math.log(b * b)
            a0 = a[0] + b * a[1]
            w1 = math.log(window1(b * b, a0, r)) / math.log(b * b)
            print(
                f"  b = {b}, a = {a}: box reading {w2:.4f}, Morton digit {a0} at base {b * b} reading {w1:.4f}, difference {w2 - w1:+.4f}"
            )
    print(f"runtime {time.time() - t0:.1f} s")


# THE COUNT


def design(b, a, k):
    F = [(x, y) for x in range(b) for y in range(b) if (x, y) != a]
    fx = np.array([v[0] for v in F], dtype=np.int64)
    fy = np.array([v[1] for v in F], dtype=np.int64)
    X, Y = np.zeros(1, dtype=np.int64), np.zeros(1, dtype=np.int64)
    for _ in range(k):
        X = (X[:, None] * b + fx[None, :]).ravel()
        Y = (Y[:, None] * b + fy[None, :]).ravel()
    return X, Y, F


def gprimes(b):
    out = []
    for p in range(2, b + 1):
        if b % p == 0 and all(p % q for q in range(2, p)):
            if p == 2:
                out.append((1, 1))
            elif p % 4 == 3:
                out.append((p, 0))
            else:
                u = next(u for u in range(1, p) if (u * u + 1) % p == 0)
                out += [(p, u), (p, p - u)]
    return out


def coprime(v, P):
    x, y = v
    for p, u in P:
        if (p, u) == (1, 1):
            if (x + y) % 2 == 0:
                return False
        elif u == 0:
            if x % p == 0 and y % p == 0:
                return False
        elif (x + u * y) % p == 0:
            return False
    return True


def count():
    t0 = time.time()
    print(
        "THE COUNT: sum of the Gaussian von Mangoldt function over S_k against (4/pi) kappa |S_k|, kappa = prod_(pi | b) N(pi)/(N(pi) - 1) times the share of F coprime to b"
    )
    print(
        "  a reading of the main term's shape at small bases, never a bound; 4/pi is the density of Lambda per lattice point, 1 over the residue pi/4 of the Dedekind zeta of Q(i)"
    )
    rows = [
        (3, (0, 0), 7),
        (3, (1, 1), 7),
        (3, (2, 0), 7),
        (4, (0, 0), 6),
        (4, (1, 0), 6),
        (4, (3, 3), 6),
        (5, (2, 2), 5),
        (5, (1, 2), 5),
        (7, (3, 3), 4),
        (7, (0, 1), 4),
    ]
    top = max(2 * (b**k) ** 2 for b, a, k in rows)
    sv = np.ones(top + 1, dtype=bool)
    sv[:2] = False
    for p in range(2, math.isqrt(top) + 1):
        if sv[p]:
            sv[p * p :: p] = False
    for b, a, k in rows:
        X, Y, F = design(b, a, k)
        P = gprimes(b)
        cop = sum(coprime(v, P) for v in F)
        kap = cop / len(F)
        for p, u in P:
            nrm = 2 if (p, u) == (1, 1) else (p * p if u == 0 else p)
            kap *= nrm / (nrm - 1)
        N = X * X + Y * Y
        inner = (X > 0) & (Y > 0)
        lam = np.where(inner & sv[N], np.log(np.maximum(N, 1)), 0.0).sum()
        ax = np.where(X == 0, Y, np.where(Y == 0, X, 0))
        lam += np.where(
            (ax > 0) & sv[ax] & (ax % 4 == 3), 2 * np.log(np.maximum(ax, 1)), 0.0
        ).sum()
        pred = 4 / math.pi * kap * len(X)
        print(
            f"  base {b} missing {a} level {k}: |S_k| = {len(X)}, kappa = {kap:.6f}, ratio {lam / pred:.4f}"
        )
    print(f"runtime {time.time() - t0:.1f} s")


# THE MINOR-ARC LEMMA, READ


def minor():
    t0 = time.time()
    lo, hi = 1, 2048
    X = 2 * (hi - 1) ** 2
    L = math.log(X)
    pp = np.zeros(X + 1, dtype=np.int64)
    sv = np.ones(X + 1, dtype=bool)
    sv[:2] = False
    for p in range(2, X + 1):
        if p * p > X:
            break
        if sv[p]:
            sv[p * p :: p] = False
    pp[sv] = np.nonzero(sv)[0]
    for p in np.nonzero(sv[: math.isqrt(X) + 1])[0]:
        q = int(p) * int(p)
        while q <= X:
            pp[q] = p
            q *= int(p)
    xs = np.arange(lo, hi, dtype=np.int64)
    Z1, Z2 = np.meshgrid(xs, xs, indexing="ij")
    N = Z1 * Z1 + Z2 * Z2
    P = pp[N]
    g = np.gcd(Z1, Z2)
    keep = (P > 0) & ((P == 2) | (g % np.maximum(P, 1) != 0))
    z1, z2 = Z1[keep].astype(float), Z2[keep].astype(float)
    lam = np.log(P[keep].astype(float))
    npr = int(((P > 0) & (P == N) & keep).sum())
    mass = float(lam.sum())
    least = min(
        X**0.875 + X / math.sqrt(n) + X ** (11 / 16) * math.sqrt(n) + math.sqrt(X) * n
        for n in range(1, X)
    )
    print(
        f"THE MINOR-ARC LEMMA, READ: S(xi) = sum Lambda(z) e(Re(z xi)) over the box [{lo}, {hi})^2, X = {X}, prime powers included: {len(lam)} prime-power elements, {npr} of them prime"
    )
    print(
        "  bound(xi) = X^(7/8) + X N(q)^(-1/2) + X^(11/16) N(q)^(1/2) + X^(1/2) N(q), constant 1 and no (log X)^4, at a sampled approximation lambda/q, half the samples within 4/N(q) of lambda/q at random and half within 1/X"
    )
    print(
        f"  the box's whole mass sum Lambda = {mass:.0f} = {mass / least:.4f} of the least value of the bound over N(q), so the reading can fail"
    )
    assert mass > least
    rs = np.random.default_rng(11)
    worst = 0.0
    for ell in (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45):
        rows = []
        for _ in range(100):
            while True:
                r = math.sqrt(X**ell) * math.sqrt(rs.uniform(0.7, 1.4))
                th = rs.uniform(0, 2 * math.pi)
                q = complex(round(r * math.cos(th)), round(r * math.sin(th)))
                nq = int(q.real) ** 2 + int(q.imag) ** 2
                if nq >= 1:
                    break
            while True:
                lm = complex(
                    int(rs.integers(-3 * nq - 3, 3 * nq + 3)),
                    int(rs.integers(-3 * nq - 3, 3 * nq + 3)),
                )
                w = lm * q.conjugate()
                if nq == 1 or all(
                    (
                        int(w.real) % p
                        or int(w.imag) % p
                        or (int(lm.real) ** 2 + int(lm.imag) ** 2) % p
                    )
                    for p in range(2, nq + 1)
                    if nq % p == 0 and sv[p]
                ):
                    break
            rad = (
                4.0 / nq * math.sqrt(rs.uniform(0, 1))
                if _ % 2 == 0
                else rs.uniform(0, 1) / X
            )
            ph = rs.uniform(0, 2 * math.pi)
            xi = lm / q + rad * complex(math.cos(ph), math.sin(ph))
            val = abs(np.sum(lam * np.exp(2j * np.pi * (z1 * xi.real - z2 * xi.imag))))
            bnd = (
                X**0.875
                + X / math.sqrt(nq)
                + X ** (11 / 16) * math.sqrt(nq)
                + math.sqrt(X) * nq
            )
            rows.append((val, bnd, nq))
        mv = max(r[0] for r in rows)
        rb = max(r[0] / r[1] for r in rows)
        worst = max(worst, rb)
        print(
            f"  N(q) near X^{ell:.2f}: max |S| = X^{math.log(mv) / L:.4f}, lemma exponent {max(0.875, 1 - ell / 2, 11 / 16 + ell / 2, 0.5 + ell):.4f}, max |S|/bound {rb:.4f}"
        )
    assert worst < 1.0
    print(
        f"  largest |S|/bound over all 1000 samples {worst:.4f} < 1 with no logarithm; the lemma allows (log X)^4 = {L**4:.0f} times more"
    )
    print(f"runtime {time.time() - t0:.1f} s")


# THE GATE AT ONE FIFTH, READ


def ggcd_int(a, b):
    a, b = (int(a.real), int(a.imag)), (int(b.real), int(b.imag))
    while b != (0, 0):
        n = b[0] * b[0] + b[1] * b[1]
        p0, p1 = a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1]
        q0, q1 = (2 * p0 + n) // (2 * n), (2 * p1 + n) // (2 * n)
        a, b = b, (a[0] - (q0 * b[0] - q1 * b[1]), a[1] - (q0 * b[1] + q1 * b[0]))
    return a[0] * a[0] + a[1] * a[1]


def lam_box(y):
    X = 2 * (y - 1) ** 2
    sv = np.ones(X + 1, dtype=bool)
    sv[:2] = False
    for p in range(2, math.isqrt(X) + 1):
        if sv[p]:
            sv[p * p :: p] = False
    small = np.nonzero(sv[: math.isqrt(X) + 1])[0]
    pw, pb = [], []
    for p in small:
        q = int(p) * int(p)
        while q <= X:
            pw.append(q)
            pb.append(int(p))
            q *= int(p)
    o = np.argsort(pw)
    pw, pb = np.array(pw, dtype=np.int64)[o], np.array(pb, dtype=np.int64)[o]
    xs = np.arange(1, y, dtype=np.int64)
    Z1, Z2, LM = [], [], []
    for s in range(0, len(xs), 256):
        a = xs[s : s + 256]
        A, B = np.meshgrid(a, xs, indexing="ij")
        N = A * A + B * B
        P = np.where(sv[N], N, 0)
        j = np.clip(np.searchsorted(pw, N), 0, len(pw) - 1)
        hit = (pw[j] == N) & (P == 0)
        P = np.where(hit, pb[j], P)
        g = np.gcd(A, B)
        keep = (P > 0) & ((P == 2) | (g % np.maximum(P, 1) != 0))
        Z1.append(A[keep])
        Z2.append(B[keep])
        LM.append(np.log(P[keep].astype(float)))
    return (
        X,
        np.concatenate(Z1).astype(float),
        np.concatenate(Z2).astype(float),
        np.concatenate(LM),
    )


def sample_q(rs, nq_target):
    while True:
        r = math.sqrt(nq_target) * math.sqrt(rs.uniform(0.7, 1.4))
        th = rs.uniform(0, 2 * math.pi)
        q = complex(round(r * math.cos(th)), round(r * math.sin(th)))
        nq = int(q.real) ** 2 + int(q.imag) ** 2
        if nq >= 2:
            break
    while True:
        lm = complex(int(rs.integers(0, nq)), int(rs.integers(0, nq)))
        if ggcd_int(lm, q) == 1:
            return q, lm, nq


def gate():
    t0 = time.time()
    print(
        "THE GATE AT ONE FIFTH, READ: S(xi) = sum Lambda(z) e(Re(z xi)) over the box [1, y)^2, y = 3^6, 3^7, 3^8, prime powers included, against y^(8/5)"
    )
    print(
        "  xi = lambda/q + beta with N(q) near y^theta, theta in 0.8..1.2, the band of the gate 1/5; half the samples with abs(beta) <= 4/N(q) at random, half with abs(beta) <= 1/X"
    )
    rs = np.random.default_rng(5)
    rows = []
    for e in (6, 7, 8):
        y = 3**e
        X, z1, z2, lam = lam_box(y)
        mass = float(lam.sum())
        g = y**1.6
        assert mass > 10 * g
        out = []
        for th in (0.8, 0.9, 1.0, 1.1, 1.2):
            best = 0.0
            for i in range(40):
                q, lm, nq = sample_q(rs, y**th)
                rad = (
                    4.0 / nq * math.sqrt(rs.uniform(0, 1))
                    if i % 2 == 0
                    else rs.uniform(0, 1) / X
                )
                ph = rs.uniform(0, 2 * math.pi)
                xi = lm / q + rad * complex(math.cos(ph), math.sin(ph))
                v = abs(
                    np.sum(lam * np.exp(2j * np.pi * (z1 * xi.real - z2 * xi.imag)))
                )
                best = max(best, v)
            out.append(best / g)
        rows.append((y, out))
        print(
            f"  y = {y}: {len(lam)} prime-power elements, mass {mass / g:.2f} y^(8/5); max abs(S)/y^(8/5) by band "
            + ", ".join(
                f"y^{th:.1f}: {r:.4f}" for th, r in zip((0.8, 0.9, 1.0, 1.1, 1.2), out)
            )
        )
    for j in range(5):
        assert rows[-1][1][j] < rows[0][1][j]
    print(
        "  in every band the ratio falls from y = 729 to y = 6561: the reading does not grow with y"
    )
    print(f"runtime {time.time() - t0:.1f} s")


# THE COUNT LEMMA AND THE BILINEAR BOUND, CHECKED


def tdist(w1, w2):
    return np.hypot(w1 - np.round(w1), w2 - np.round(w2))


def disc(H):
    r = math.isqrt(int(H))
    a = np.arange(-r, r + 1, dtype=np.int64)
    A, B = np.meshgrid(a, a, indexing="ij")
    k = A * A + B * B <= H
    return A[k].astype(float), B[k].astype(float)


def hcount():
    t0 = time.time()
    H = 10**6
    h1, h2 = disc(H)
    print(
        f"THE COUNT LEMMA, CHECKED: #{{h : N(h) <= H, ||h xi|| <= eta}} against (H/N(q) + 1)(N(q) eta^2 + 1) at H = 10^6, {len(h1)} points, abs(xi - lambda/q) <= 4/N(q), ||w|| the Euclidean distance from w to Z[i]"
    )
    rs = np.random.default_rng(3)
    worst, wrow, top = 0.0, None, 0
    for e in np.linspace(0, 5.5, 23):
        for i in range(8):
            q, lm, nq = sample_q(rs, 10**e)
            rad = 0.0 if i == 0 else 4.0 / nq * math.sqrt(rs.uniform(0, 1))
            ph = rs.uniform(0, 2 * math.pi)
            xi = lm / q + rad * complex(math.cos(ph), math.sin(ph))
            d = tdist(h1 * xi.real - h2 * xi.imag, h1 * xi.imag + h2 * xi.real)
            ds = np.sort(d)
            top = max(top, nq)
            for eta in (
                0.1 / math.sqrt(nq),
                0.5 / math.sqrt(nq),
                2 / math.sqrt(nq),
                1 / math.sqrt(H),
                10 / math.sqrt(H),
                0.05,
                0.3,
                0.7,
            ):
                c = int(np.searchsorted(ds, eta, side="right"))
                r = c / ((H / nq + 1) * (nq * eta * eta + 1))
                if r > worst:
                    worst, wrow = r, (nq, eta, c, rad)
    print(
        f"  184 frequencies, N(q) from 2 to {top}, eight eta each: largest ratio {worst:.4f} at N(q) = {wrow[0]}, eta = {wrow[1]:.3g}, count {wrow[2]}, abs(beta) = {wrow[3]:.3g}; the proof allows 2^18"
    )
    q, lm, Hc, e2 = complex(9, 4), complex(53, 77), 1940, 5 / 97
    assert ggcd_int(lm, q) == 1
    c1, c2 = disc(Hc)
    xi = lm / q
    c = int(
        (
            tdist(c1 * xi.real - c2 * xi.imag, c1 * xi.imag + c2 * xi.real) ** 2
            <= e2 + 1e-12
        ).sum()
    )
    r = c / ((Hc / 97 + 1) * (97 * e2 + 1))
    print(
        f"  the cell q = 9 + 4i, lambda = 53 + 77i, xi = lambda/q, H = 1940, eta^2 = 5/97: count {c}, ratio {r:.4f}; readings, not the supremum"
    )
    assert worst < 2**18 and r < 2**18
    print(f"runtime {time.time() - t0:.1f} s")


def annulus(U0):
    h1, h2 = disc(2 * U0 - 1)
    k = h1 * h1 + h2 * h2 >= U0
    return h1[k], h2[k]


def bilinear():
    t0 = time.time()
    print(
        "THE BILINEAR BOUND, CHECKED: B = sum_v sum_u alpha_v beta_u e(Re(v u xi)) over N(v) in [V0, 2 V0), N(u) in [U0, 2 U0), against norm2(alpha) norm2(beta) (U0 + V0 + N(q) + U0 V0/N(q))^(1/2)"
    )
    print(
        "  alpha_v chosen as the conjugate of sum_u beta_u e(Re(v u xi)), the extremal choice for Cauchy-Schwarz, beta unimodular at random or the von Mangoldt function"
    )
    rs = np.random.default_rng(7)
    worst, wrow = 0.0, None
    for U0, V0 in ((400, 400), (1500, 1500), (300, 6000), (6000, 300), (2500, 2500)):
        u1, u2 = annulus(U0)
        v1, v2 = annulus(V0)
        for e in np.linspace(0.2, math.log10(4 * U0 * V0), 10):
            for kind in (0, 1):
                q, lm, nq = sample_q(rs, 10**e)
                rad = 4.0 / nq * math.sqrt(rs.uniform(0, 1))
                ph = rs.uniform(0, 2 * math.pi)
                xi = lm / q + rad * complex(math.cos(ph), math.sin(ph))
                if kind == 0:
                    beta = np.exp(2j * np.pi * rs.uniform(0, 1, len(u1)))
                else:
                    beta = np.array(
                        [
                            math.log(max(2, int(a * a + b * b)))
                            if all(
                                (int(a * a + b * b)) % p
                                for p in range(2, math.isqrt(int(a * a + b * b)) + 1)
                            )
                            else 0.0
                            for a, b in zip(u1, u2)
                        ]
                    )
                    if not beta.any():
                        continue
                ur, ui = u1 * xi.real - u2 * xi.imag, u1 * xi.imag + u2 * xi.real
                tot = 0.0
                for c in range(0, len(v1), 800):
                    F = (
                        np.exp(
                            2j
                            * np.pi
                            * (
                                np.outer(v1[c : c + 800], ur)
                                - np.outer(v2[c : c + 800], ui)
                            )
                        )
                        @ beta
                    )
                    tot += float(np.sum(np.abs(F) ** 2))
                val = math.sqrt(tot) / np.linalg.norm(beta)
                r = val / math.sqrt(U0 + V0 + nq + U0 * V0 / nq)
                if r > worst:
                    worst, wrow = r, (U0, V0, nq, kind)
    print(
        f"  on the annuli, 5 shapes, 10 bands of N(q) from 2 to 4 U0 V0, two kinds of beta: largest ratio {worst:.4f} at U0 = {wrow[0]}, V0 = {wrow[1]}, N(q) = {wrow[2]}, beta {'random' if wrow[3] == 0 else 'Lambda'}; the proof allows 2 * 10^4"
    )
    top = 0.0
    for U0, q, lm in (
        (200, complex(1, 0), complex(0, 0)),
        (100, complex(1, 1), complex(1, 0)),
    ):
        d1, d2 = disc(2 * U0 - 1)
        xi = lm / q
        nq = int(q.real) ** 2 + int(q.imag) ** 2
        Mx = np.exp(
            2j
            * np.pi
            * (
                np.outer(d1, d1 * xi.real - d2 * xi.imag)
                - np.outer(d2, d1 * xi.imag + d2 * xi.real)
            )
        )
        s1 = float(np.linalg.norm(Mx, 2))
        r = s1 / math.sqrt(2 * U0 + nq + U0 * U0 / nq)
        top = max(top, r)
        print(
            f"  operator norm on the discs N(v), N(u) < {2 * U0}, {len(d1)} points each, q = {int(q.real)} + {int(q.imag)}i, xi = lambda/q: sigma_1 = {s1:.1f}, {r:.4f} times the bracket"
        )
    assert worst < 2e4 and top < 2e4
    print(f"runtime {time.time() - t0:.1f} s")


# THE SEPARATION COST OF A SQUARE, READ


def sstep(t):
    t = np.clip(t, 0.0, 1.0)
    with np.errstate(divide="ignore", over="ignore"):
        a = np.where(t > 0, np.exp(-1.0 / np.maximum(t, 1e-300)), 0.0)
        b = np.where(t < 1, np.exp(-1.0 / np.maximum(1 - t, 1e-300)), 0.0)
    return a / (a + b)


def anorm(f):
    A = np.abs(np.fft.rfft2(f))
    return float((A[:, 0].sum() + A[:, -1].sum() + 2 * A[:, 1:-1].sum()) / f.size)


def sqnorm(h, n):
    rho = 0.25 + np.arange(n) / n
    th = 0.35 + 0.9 * np.arange(n) / n
    er = np.exp(rho)[:, None]
    z1 = er * np.cos(th)[None, :]
    sq = sstep((z1 - 1) / h) * sstep((2 - z1) / h)
    del z1
    z2 = er * np.sin(th)[None, :]
    sq *= sstep((z2 - 1) / h) * sstep((2 - z2) / h)
    del z2
    assert (
        sq[0, :].max() == 0
        and sq[:, 0].max() == 0
        and sq[-1, :].max() == 0
        and sq[:, -1].max() == 0
    )
    return anorm(sq)


def sep():
    t0 = time.time()
    print(
        "THE SEPARATION COST OF A SQUARE, READ: the Fourier l^1 norm in (log abs(z), arg z) of the square [1, 2]^2 smoothed at scale h_s, which is the L^1 cost of writing W(v u) through the characters abs(z)^(i t_M) (z/abs z)^l, against a polar box smoothed at the same scale"
    )
    ra, pr, ta, pt = 0.25, 1.0, 0.35, 0.9
    rows = []
    for k in (4, 5, 6, 7, 8, 9):
        h = 2.0**-k
        n = 2 ** (k + 5) if k < 8 else 2 ** (k + 4)
        rho = ra + pr * np.arange(n) / n
        th = ta + pt * np.arange(n) / n
        a = sqnorm(h, n)
        dc = h / 2
        po = (
            sstep((rho - 0.45) / dc)[:, None]
            * sstep((1.05 - rho) / dc)[:, None]
            * (sstep((th - 0.5) / dc) * sstep((1.1 - th) / dc))[None, :]
        )
        rows.append((h, a, anorm(po)))
        del po
    for i, (h, a, b) in enumerate(rows):
        g = (
            ""
            if i == 0
            else f", growth {a / rows[i - 1][1]:.4f} against {b / rows[i - 1][2]:.4f}"
        )
        print(
            f"  h_s = 1/{round(1 / h)}: square {a:.3f}, polar box {b:.3f}, square (h_s)^(1/2) {a * math.sqrt(h):.4f}{g}"
        )
    for k, n in ((7, 2**11), (8, 2**13)):
        a = sqnorm(2.0**-k, n)
        b = rows[k - 4][1]
        print(
            f"  grid check at h_s = 1/{2**k}: {a:.3f} on {n}^2 points against {b:.3f}, relative change {abs(a - b) / b:.4f}"
        )
        assert abs(a - b) / b < 0.005
    assert rows[-1][1] / rows[-2][1] > 1.3 and rows[-1][2] / rows[-2][2] < 1.2
    print(f"runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    {
        "chain": chain,
        "census": census,
        "floor": floor1,
        "uwin": uwin,
        "morton": morton,
        "count": count,
        "minor": minor,
        "gate": gate,
        "hcount": hcount,
        "bilinear": bilinear,
        "sep": sep,
    }[sys.argv[1]]()
