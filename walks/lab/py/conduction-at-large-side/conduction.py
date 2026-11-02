import math
import sys
import time
from fractions import Fraction

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

R3 = math.sqrt(3.0)
SIGMA = 1.0 / R3

# RENDER


def odd_digits(n, level, k=1):
    i = np.arange(n**level)
    out = []
    for _ in range(level):
        out.append((i % n) % 2 == 1)
        i = i // n
    out = np.array(out)
    if k > 1:
        out = np.repeat(out, k, axis=1)
    return out


def render(n, level, k=1):
    o = odd_digits(n, level, k)
    void = np.zeros((o.shape[1], o.shape[1]), bool)
    for row in o:
        void |= np.outer(row, row)
    return ~void


def obnosov(z):
    return math.sqrt((1 + 3 * z) / (3 + z))


def fill(n):
    return n * n - ((n - 1) // 2) ** 2


# NETWORK


def pair(a, b):
    s = a + b
    return np.where(s > 0, 2 * a * b / np.where(s > 0, s, 1), 0.0)


def solve(n_nodes, rows, cols, w, diag, rhs):
    A = sp.coo_matrix((-w, (rows, cols)), shape=(n_nodes, n_nodes))
    A = (A + A.T + sp.diags(diag)).tocsc()
    return spl.spsolve(A, rhs, permc_spec="COLAMD")


def network(a, wx, wy, left, right, right_value):
    _nx, _ny = a.shape
    live = a > 0
    idx = -np.ones(a.shape, np.int64)
    idx[live] = np.arange(live.sum())
    n_nodes = int(live.sum())
    rows, cols, ws = [], [], []
    gx = pair(a[:-1, :], a[1:, :]) * wx[:-1, :]
    gy = pair(a[:, :-1], a[:, 1:]) * wy[:, :-1]
    for g, i0, i1 in ((gx, idx[:-1, :], idx[1:, :]), (gy, idx[:, :-1], idx[:, 1:])):
        m = g > 0
        rows.append(i0[m])
        cols.append(i1[m])
        ws.append(g[m])
    rows, cols, ws = np.concatenate(rows), np.concatenate(cols), np.concatenate(ws)
    diag = np.bincount(rows, ws, n_nodes) + np.bincount(cols, ws, n_nodes)
    rhs = np.zeros(n_nodes)
    lm = live[0, :]
    diag[idx[0, lm]] += left[lm]
    rm = live[-1, :] & (right > 0)
    diag[idx[-1, rm]] += right[rm]
    rhs[idx[-1, rm]] += right[rm] * right_value
    u = solve(n_nodes, rows, cols, ws, diag, rhs)
    return float(np.sum(left[lm] * u[idx[0, lm]]))


def conductance_full(a):
    a = np.asarray(a, float)
    one = np.ones(a.shape)
    return network(a, one, one, 2 * a[0, :], 2 * a[-1, :], 1.0)


def conductance(a):
    a = np.asarray(a, float)
    s = a.shape[0]
    h = (s + 1) // 2
    q = a[:h, :h].copy()
    wx = np.ones(q.shape)
    wy = np.ones(q.shape)
    rw = np.ones(h)
    if s % 2:
        wx[:, -1] = 0.5
        rw[-1] = 0.5
        center = q[-1, :] > 0
        left = 2 * q[0, :] * rw
        g = pair(q[-2, :], q[-1, :]) * rw
        q2 = q[:-1, :]
        wx2, wy2 = wx[:-1, :], wy[:-1, :]
        right = np.where(center, g, 0.0)
        return 2 * network(q2, wx2, wy2, left, right, 0.5)
    left = 2 * q[0, :]
    right = 2 * q[-1, :]
    return 2 * network(q, wx, wy, left, right, 0.5)


def rich(g1, g2, r=2.0, p=4.0 / 3.0):
    return g2 + (g2 - g1) / (r**p - 1)


def aitken(g1, g2, g3, r=2.0):
    e1, e2 = g2 - g1, g3 - g2
    p = math.log(e1 / e2) / math.log(r)
    return rich(g2, g3, r, p), p


def level(n, lv, k):
    return conductance(render(n, lv, k).astype(float))


# CELL


def cell(m, z):
    a = np.ones((2 * m, 2 * m))
    a[m:, m:] = z
    return conductance_full(a)


def cell3(m, rtol=1e-11):
    t = np.arange(2 * m) >= m
    up = t[:, None, None].astype(int) + t[None, :, None] + t[None, None, :]
    live = up <= 1
    s = 2 * m
    idx = -np.ones(live.shape, np.int64)
    n_nodes = int(live.sum())
    idx[live] = np.arange(n_nodes)
    rows, cols = [], []
    for ax in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[ax] = slice(0, -1)
        b[ax] = slice(1, None)
        both = live[tuple(a)] & live[tuple(b)]
        rows.append(idx[tuple(a)][both])
        cols.append(idx[tuple(b)][both])
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    w = np.ones(len(rows))
    diag = np.bincount(rows, w, n_nodes) + np.bincount(cols, w, n_nodes)
    left = idx[0][live[0]]
    right = idx[-1][live[-1]]
    diag[left] += 2.0
    diag[right] += 2.0
    rhs = np.zeros(n_nodes)
    rhs[right] = 2.0
    A = sp.coo_matrix((-w, (rows, cols)), shape=(n_nodes, n_nodes))
    A = (A + A.T + sp.diags(diag)).tocsr()
    x0 = (np.indices(live.shape)[0][live] + 0.5) / s
    u, info = spl.cg(A, rhs, x0=x0, M=sp.diags(1.0 / diag), rtol=rtol, maxiter=20000)
    assert info == 0
    return 2.0 * float(np.sum(u[left])) / s


def verb_cell():
    ms = (64, 128, 256, 512)
    g = [cell(m, 0.0) for m in ms]
    for i, m in enumerate(ms):
        line = f"insulating M={m} sigma={g[i]:.10f} err={g[i] - SIGMA:+.3e}"
        if i:
            line += f" order={math.log((g[i - 1] - SIGMA) / (g[i] - SIGMA)) / math.log(2):.4f} rich={rich(g[i - 1], g[i]):.10f}"
        print(line)
    e1, e2 = rich(g[1], g[2]), rich(g[2], g[3])
    bar = abs(e2 - e1)
    print(
        f"insulating sigma={e2:.10f} bar={bar:.1e} target={SIGMA:.10f} off={e2 - SIGMA:+.1e}"
    )
    assert abs(e2 - SIGMA) < 2 * bar + 1e-9
    zs = (0.01, 1 / 9, 1 / 3, 3.0, 9.0, 100.0)
    val = {}
    for z in zs:
        h = [cell(m, z) for m in (64, 128, 256)]
        val[z], p = aitken(*h)
        print(
            f"contrast z={z:.6g} sigma={val[z]:.9f} order={p:.4f} obnosov={obnosov(z):.9f} off={val[z] - obnosov(z):+.1e}"
        )
        assert abs(val[z] - obnosov(z)) < 6e-8
    for z in (0.01, 1 / 9, 1 / 3):
        prod = val[z] * val[1 / z if z != 1 / 3 else 3.0]
        print(f"keller z={z:.6g} sigma(z)*sigma(1/z)={prod:.9f}")
        assert abs(prod - 1) < 7e-8
    big = 1e8
    h = [cell(m, big) for m in (64, 128, 256)]
    s_inf, p = aitken(*h)
    print(
        f"conductor z=1e8 sigma={s_inf:.9f} order={p:.4f} sqrt3={R3:.9f} keller={e2 * s_inf:.9f}"
    )


def verb_side():
    for n in (11, 21, 41, 81):
        row = []
        for k in (4, 8, 16):
            row.append(n * (level(n, 1, k) - cell(k // 2, 0.0)))
        print(
            f"level1 N={n} N*(sigma_k(N,1)-sigma_k(inf,1)) k=4,8,16: "
            + " ".join(f"{x:.5f}" for x in row)
            + f" sigma(N,1)~{SIGMA + row[-1] / n:.6f}"
        )
    for n in (5, 7, 9, 11, 15, 21, 31, 41):
        row = []
        for k in (1, 2, 4):
            if k * n * n > 1800 and k > 1:
                break
            row.append(level(n, 1, k) / level(n, 2, k))
        print(
            f"ratio N={n} sigma(N,1)/sigma(N,2) k=1,2,4: "
            + " ".join(f"{x:.5f}" for x in row)
            + f" N*(sqrt3-ratio)={n * (R3 - row[-1]):.4f}"
        )


def up(x):
    e = math.floor(math.log10(x))
    return f"{math.ceil(x / 10**e)}e{e}" if math.ceil(x / 10**e) < 10 else f"1e{e + 1}"


def bounds(n):
    beta = (n + 1) / (2 * n) + (n - 1) / (n + 1)
    alpha = 2 * n / (n + 1)
    return beta, alpha


def ratios(n, k, cap):
    g = [1.0]
    lv = 1
    while k * n ** (lv + 1) <= cap:
        if lv == 1:
            g.append(level(n, 1, k))
        g.append(level(n, lv + 1, k))
        lv += 1
    return [g[i] / g[i + 1] for i in range(1, len(g) - 1)]


def verb_scale():
    rows = []
    for n in (3, 5, 7, 9, 11):
        beta, alpha = bounds(n)
        assert Fraction(3 * n * n + 1, 2 * n * (n + 1)) == Fraction(
            n + 1, 2 * n
        ) + Fraction(n - 1, n + 1)
        seq = {k: ratios(n, k, 3200 if k == 1 else 2700) for k in (1, 2)}
        for k in (1, 2):
            print(
                f"scale N={n} k={k} rho(N,L) L=1..{len(seq[k])}: "
                + " ".join(f"{x:.5f}" for x in seq[k])
            )
            assert all(beta <= x <= alpha for x in seq[k])
        j = len(seq[2]) - 1
        mesh = abs(seq[1][j] - seq[2][j])
        kk = min((1, 2), key=lambda k: abs(seq[k][-1] - seq[k][-2]))
        drift = abs(seq[kk][-1] - seq[kk][-2])
        rho = seq[kk][-1]
        bar = max(mesh, drift)
        m = fill(n)
        dw = math.log(m * rho) / math.log(n)
        dbar = math.log((rho + bar) / rho) / math.log(n)
        law = 2 + math.log(3 * R3 / 4) / math.log(n)
        archie = math.log(rho) / math.log(n * n / m)
        print(
            f"table N={n} k={kk} rho={rho:.5f} mesh={mesh:.1e} drift={drift:.1e} bar={bar:.1e} dw_bar={dbar:.1e} bounds=[{beta:.5f}, {alpha:.5f}] dw={dw:.4f} law={law:.4f} (dw-2)logN={(dw - 2) * math.log(n):.4f} archie={archie:.4f} drift_k1={seq[1][0] - seq[1][-1]:.2e} drift_k2={seq[2][0] - seq[2][-1]:.2e} lower={(math.log(m * beta / n**2)):.4f} upper={(math.log(m * alpha / n**2)):.4f} archie_gap={math.log(3) / math.log(16 / 9) - archie:.4f} excess={(n - 1) ** 3 / (8 * n**3):.5f} dw-law={dw - law:.4f}"
        )
        rows.append(
            f"| {n} | {kk} | "
            + ", ".join(f"{x:.5f}" for x in seq[kk])
            + f" | `{rho:.5f} +- {up(bar)}` | `[{Fraction(3 * n * n + 1, 2 * n * (n + 1))}, {Fraction(2 * n, n + 1)}]` | `{dw:.4f} +- {up(dbar)}` | {law:.4f} |"
        )
    print("\n".join(rows))
    print(
        f"limits sqrt3={R3:.5f} bounds=[1.5, 2] (dw-2)logN->{math.log(3 * R3 / 4):.4f} in [{math.log(9 / 8):.4f}, {math.log(1.5):.4f}] archie->{math.log(3) / math.log(16 / 9):.4f}"
    )


def verb_dim3():
    ms = (8, 16, 32, 64)
    g = [cell3(m) for m in ms]
    for m, x in zip(ms, g):
        print(f"dim3 M={m} sigma={x:.7f}")
    a1, p1 = aitken(*g[:3])
    a2, p2 = aitken(*g[1:])
    r2 = rich(g[2], g[3])
    bar = max(abs(a2 - a1), abs(a2 - r2))
    expo = math.log(1 / a2) / math.log(2)
    print(
        f"dim3 aitken(8,16,32)={a1:.6f} order={p1:.4f} aitken(16,32,64)={a2:.6f} order={p2:.4f} rich(32,64)={r2:.6f}"
    )
    print(
        f"dim3 sigma={a2:.5f} bar={bar:.1e} exponent=log(1/sigma)/log2={expo:.4f} band=[{math.log(1 / (a2 - bar)) / math.log(2):.4f}, {math.log(1 / (a2 + bar)) / math.log(2):.4f}]"
    )


def verb_rate():
    ns = (9, 11, 15, 21, 31, 41)
    y = np.array([level(n, 1, 2) / level(n, 2, 2) for n in ns])
    x = np.array(ns, float)
    A = np.stack([np.ones_like(x), -1 / x, 1 / x ** (4 / 3)], axis=1)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ c
    print(
        "rate N="
        + ",".join(map(str, ns))
        + " rho(N,1) k=2: "
        + " ".join(f"{v:.5f}" for v in y)
    )
    print(
        "rate N*(sqrt3-rho(N,1)): "
        + " ".join(f"{n * (R3 - v):.4f}" for n, v in zip(ns, y))
    )
    print(
        f"rate fit rho(N,1) = r - a/N + b/N^(4/3): r={c[0]:.4f} a={c[1]:.3f} b={c[2]:.3f} max residual={np.abs(res).max():.1e} sqrt3={R3:.4f}"
    )
    for q in (5 / 4, 4 / 3, 3 / 2):
        A2 = np.stack([-1 / x, 1 / x**q], axis=1)
        c2, *_ = np.linalg.lstsq(A2, y - R3, rcond=None)
        print(
            f"rate fit rho(N,1) = sqrt3 - a/N + b/N^{q:.4g}: a={c2[0]:.3f} b={c2[1]:.3f} max residual={np.abs(y - R3 - A2 @ c2).max():.1e}"
        )


def verb_dual():
    for n, lv in ((5, 1), (3, 2)):
        ext, bars = {}, {}
        for z in (0.0, 1e8):
            g = []
            for k in (8, 16, 32, 64):
                m = render(n, lv, k)
                g.append(conductance(np.where(m, 1.0, z)))
            ext[z], p = aitken(*g[1:])
            bars[z] = abs(ext[z] - aitken(*g[:3])[0])
            print(
                f"dual N={n} L={lv} z={z:.0e} k=8,16,32,64: "
                + " ".join(f"{x:.7f}" for x in g)
                + f" aitken={ext[z]:.7f} order={p:.3f} bar={bars[z]:.1e}"
            )
        prod = ext[0.0] * ext[1e8]
        bar = bars[0.0] * ext[1e8] + bars[1e8] * ext[0.0]
        print(f"dual N={n} L={lv} sigma(0)*sigma(inf)={prod:.7f} bar={bar:.1e}")
        assert abs(prod - 1) < bar


VERBS = {
    "cell": verb_cell,
    "dual": verb_dual,
    "side": verb_side,
    "rate": verb_rate,
    "scale": verb_scale,
    "dim3": verb_dim3,
}

if __name__ == "__main__":
    names = sys.argv[1:] or list(VERBS)
    for name in names:
        t = time.time()
        VERBS[name]()
        print(f"{name} {time.time() - t:.1f} s", flush=True)
