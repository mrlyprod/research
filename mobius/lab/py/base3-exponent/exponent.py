import importlib.util
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np
from mpmath import cos as mpcos
from mpmath import log as mplog
from mpmath import mp, mpf, quad
from mpmath import pi as mppi

HERE = os.path.dirname(os.path.abspath(__file__))
NORMS = os.path.join(
    HERE,
    "..",
    "..",
    "..",
    "..",
    "coprime",
    "lab",
    "py",
    "digit-transform-norms",
    "norms.py",
)


def load():
    spec = importlib.util.spec_from_file_location("norms", NORMS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return vars(mod)


LIB = load()
dn, up = LIB["dn"], LIB["up"]
B, A0 = 3, 2
DESIGN = LIB["base_missing"](B, A0)
LN3 = math.log(B)
FAILS = []


def lo7(x):
    return math.floor(x * 1e7) / 1e7


def hi7(x):
    return math.ceil(x * 1e7) / 1e7


def expect(ok, name):
    if not ok:
        FAILS.append(name)
    return "ok" if ok else "MISMATCH"


def one_digit():
    total = Fraction(1) + 2 * Fraction(1, 2)
    top = hi7(math.log(2) / LN3)
    print(
        f"one digit: sum over t in Z/3 of |hat F(t/3)| = 1 + 1/2 + 1/2 = {total} exactly, so max_x Sigma_1 >= 2 at x = 0 and alpha_1 <= log_3 2 < {top:.7f} once the max is certified {expect(total == 2, 'digit sum')}"
    )


def entropy_floor():
    mp.dps = 30
    jensen = quad(lambda t: mplog(abs(mpcos(mppi * t))), [0, mpf(1) / 2, 1])
    exact = -mplog(2)
    floor = 1 + exact / mplog(3)
    fl = lo7(float(floor))
    print(
        f"entropy floor: int_0^1 log|cos(pi t)| dt = {float(jensen):.12f} against -log 2 = {float(exact):.12f} (Jensen's formula, int log|1 + e(u)| = 0), so I_level >= 2^(-level) and alpha_1 >= 1 - log_3 2 = log(3/2)/log 3 > {fl:.7f} {expect(abs(jensen - exact) < mpf(10) ** -9, 'jensen')}"
    )
    return fl


def sandwich(N, m, chunk=2048):
    t0 = time.time()
    F = DESIGN[2]
    Q = B**N * m
    LIB["costable"](Q)
    dif, k = LIB["folded"](F), len(F)
    n = B**N
    grid = np.arange(n, dtype=np.int64) * m
    slack = up(
        up(up(LIB["lip_up"](F) * ((B**N - 1) // (B - 1))) * B**N) / (2 * B**N * m)
    )
    best_lo, best_hi = np.inf, 0.0
    for s0 in range(0, m, chunk):
        sh = np.arange(s0, min(s0 + chunk, m), dtype=np.int64)
        pl = np.ones((len(sh), n))
        ph = np.ones((len(sh), n))
        for j in range(N):
            kc = [((B**j) * (grid[None, :] + sh[:, None])) % Q]
            pl = dn(pl * LIB["hat_iv"](kc, Q, dif, k, False))
            ph = up(ph * LIB["hat_iv"](kc, Q, dif, k, True))
        al = np.zeros(len(sh))
        ah = np.zeros(len(sh))
        for c in range(n):
            al = dn(al + pl[:, c])
            ah = up(ah + ph[:, c])
        best_lo = min(best_lo, float(al.min()))
        best_hi = max(best_hi, float(ah.max()))
    cl, ch = float(dn(best_lo - slack)), float(up(best_hi + slack))
    el = lo7(math.log(cl) / (N * LN3))
    eh = hi7(math.log(ch) / (N * LN3))
    assert cl > B ** (N * el) and ch < B ** (N * eh)
    print(
        f"CVY sandwich N={N} m={m}: min_x Sigma_N > {cl:.6f}, max_x Sigma_N < {ch:.6f}, slack {slack:.2e} -> alpha_1 in ({el:.7f}, {eh:.7f})  ({time.time() - t0:.1f}s)"
    )
    return el, eh, cl, ch


def masked_perron(G, nd, keep, iters=6000, tol=1e-14):
    S, U = B ** (nd - 1), B ** (nd - 2)
    G3 = G.reshape(B, U, B)
    y = np.where(keep, 1.0, 0.0)
    z = np.empty((B, U))
    prev = 0.0
    for it in range(iters):
        Y = y.reshape(U, B)
        for c in range(B):
            z[c] = (G3[c] * Y).sum(axis=1)
        top = float(z.max())
        if top <= 0.0:
            return y
        y = np.where(keep, (z / top).reshape(S), 0.0)
        if it > 300 and abs(top - prev) <= tol * top:
            break
        prev = top
    return y


def window(nd, m):
    t0 = time.time()
    Ghi, Glo, slack = LIB["base_windows"](B, A0, nd, m, 1)
    _l, yh = LIB["perron_red"](Ghi, B, nd)
    mu_hi = LIB["cw_red"](Ghi, yh, B, nd, True)
    eh = hi7(math.log(mu_hi) / LN3)
    assert mu_hi < B**eh
    _l, yl = LIB["perron_red"](Glo, B, nd)
    keep = yl > 1e-9 * yl.max()
    mu_tree = LIB["cw_red"](Glo, yl, B, nd, False, keep)
    tree = lo7(math.log(mu_tree) / LN3) if mu_tree > 0 else None
    rounds = 0
    for rounds in range(1, 41):
        y = masked_perron(Glo, nd, keep)
        nxt = y > 1e-6 * y.max()
        if (nxt == keep).all():
            break
        keep = nxt
    mu_lo = LIB["cw_red"](Glo, np.where(keep, y, 1.0), B, nd, False, keep)
    el = lo7(math.log(mu_lo) / LN3) if mu_lo > 0 else None
    if el is not None:
        assert mu_lo > B**el
    lo_txt = "none" if el is None else f"{el:.7f}"
    tr_txt = "none" if tree is None else f"{tree:.7f}"
    print(
        f"window {nd} digits sub-scan {m}: slack {slack:.2e} mu_hi {mu_hi:.7f} alpha_1 < {eh:.7f}; masked lower {lo_txt} on {int(keep.sum())} of {B ** (nd - 1)} states after {rounds} rounds, unmasked restriction {tr_txt}  ({time.time() - t0:.1f}s)"
    )
    return el, eh, tree


def grid_reading(top):
    prev = None
    for L in range(1, top + 1):
        a = np.arange(B**L, dtype=np.float64)
        p = np.ones(B**L)
        for j in range(L):
            p *= np.abs(np.cos(np.pi * a * B**j / B**L))
        s = float(p.sum())
        step = "" if prev is None else f"  log_3 step {math.log(s / prev) / LN3:.6f}"
        print(f"  reading Sigma_{L}(0) = {s:.6f}{step}")
        prev = s


def check():
    t0 = time.time()
    LIB["hat_missing_check"](B, A0)
    one_digit()
    fl = entropy_floor()
    _a, hi1, _cl, ch1 = sandwich(1, 4096)
    _a, hi2, _cl2, _ch2 = sandwich(2, 20000)
    print(
        f"  {expect(hi1 >= hi2 and ch1 >= 2.0 - 1e-9, 'one step bound')} N=1 bound {hi1:.7f} against N=2 bound {hi2:.7f}"
    )
    el, eh, _tree = window(6, 16)
    print(
        f"  {expect(fl < el < eh < hi2, 'order')} floor {fl:.7f} < lower {el:.7f} < upper {eh:.7f} < N=2 bound {hi2:.7f}"
    )
    print(f"runtime {time.time() - t0:.1f}s")


def enclosure():
    t0 = time.time()
    LIB["hat_missing_check"](B, A0)
    one_digit()
    fl = entropy_floor()
    sandwich(1, 100000)
    _a, _hi2, cl2, ch2 = sandwich(2, 500000)
    norm = float(up(4.0 * ch2))
    step = math.ceil(math.log(ch2) / (2 * LN3) * 1e7) / 1e7
    print(
        f"two-step norm: ||T^2 1|| = fill^2 max_x Sigma_2 < {norm:.6f} (max_x Sigma_2 in ({cl2:.6f}, {ch2:.6f}) so 4 max_x Sigma_2 = 15.5175 to six figures), alpha_1 <= (1/2) log_3 max_x Sigma_2 < {step:.7f}"
    )
    sandwich(4, 4000)
    sandwich(6, 1500)
    lows, highs = [fl], [step]
    for nd, m in ((6, 16), (8, 16), (9, 16), (10, 16), (11, 8), (12, 16)):
        el, eh, _tree = window(nd, m)
        if el is not None:
            lows.append(el)
        highs.append(eh)
    lo, hi = max(lows), min(highs)
    one = hi7(math.log(2) / LN3)
    ok = expect(fl < lo <= hi < step < one, "enclosure")
    print(
        f"certified enclosure alpha_1(3, {{0,1}}) in [{lo:.7f}, {hi:.7f}], width {hi - lo:.7f}; floor {fl:.7f}, two-step bound {step:.7f}, one-step bound {one:.7f} {ok}"
    )
    print(
        f"hat kappa_1 = 1 - alpha_1 in [{1 - hi:.7f}, {1 - lo:.7f}] against 1/2 (Chow, Varju and Yu Remark 6.1: below 1/2 at base 3), floor against the bars 1/4 and 1/5: {expect(fl > 0.25, 'bars')}"
    )
    grid_reading(12)
    print(f"runtime {time.time() - t0:.1f}s")


def main():
    verb = sys.argv[1] if len(sys.argv) > 1 else "check"
    if verb == "check":
        check()
    elif verb == "enclosure":
        enclosure()
    else:
        raise SystemExit("verbs: check enclosure")
    if FAILS:
        raise SystemExit(f"{len(FAILS)} rows off target: {FAILS}")
    print("every asserted row matches")


if __name__ == "__main__":
    main()
