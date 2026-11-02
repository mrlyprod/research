import sys
import time
from fractions import Fraction
from itertools import combinations, product

import numpy as np

CHECK_TOP = 8
READ_TOP = 13
MIDS = [(5, 100), (7, 1000), (8, 2000), (8, 5000)]
SMALL = [(5, 2, 1), (5, 2, 2), (5, 2, 3), (5, 2, 4), (4, 2, 3)]

# ARITHMETIC


def mobius(n):
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, n + 1):
        if not sieve[p]:
            continue
        sieve[p * p :: p] = False
        mu[p::p] *= -1
        mu[p * p :: p * p] = 0
    return mu


def digits(values, base, level):
    out = np.zeros((len(values), level), dtype=np.int64)
    v = np.array(values, dtype=np.int64)
    for j in range(level):
        out[:, j] = v % base
        v //= base
    return out


def level_set(base, digs, level):
    vals = np.zeros(1, dtype=np.int64)
    for j in range(level):
        vals = (
            vals[:, None] + np.array(digs, dtype=np.int64)[None, :] * base**j
        ).ravel()
    return np.sort(vals)


def rep(values, n):
    out = values % n
    out[out == 0] = n
    return out


def residue_vector(mu, x, n):
    g = np.zeros(n, dtype=np.int64)
    idx = np.arange(1, x + 1)
    np.add.at(g, idx % n, mu[1 : x + 1])
    return g


# FAMILIES


def translates(base, digs, level, mu, x):
    n = base**level
    d = level_set(base, digs, level)
    total = 0
    for s0 in range(0, n, 512):
        s = np.arange(s0, min(n, s0 + 512), dtype=np.int64)
        members = rep(d[None, :] + s[:, None], n)
        vals = np.where(members <= x, mu[members], 0).sum(axis=1)
        total += int((vals * vals).sum())
    return Fraction(total, n)


def shifts(base, digs, level, mu, x):
    n = base**level
    fd = digits(level_set(base, digs, level), base, level)
    total = 0
    for c0 in range(0, n, 512):
        cs = digits(np.arange(c0, min(n, c0 + 512)), base, level)
        members = np.zeros((len(cs), len(fd)), dtype=np.int64)
        for j in range(level):
            members += ((fd[None, :, j] + cs[:, None, j]) % base) * base**j
        members = rep(members, n)
        vals = np.where(members <= x, mu[members], 0).sum(axis=1)
        total += int((vals * vals).sum())
    return Fraction(total, n)


def subsets(base, fill, level, mu, x):
    n = base**level
    subs = list(combinations(range(base), fill))
    total = 0
    count = 0
    for choice in product(subs, repeat=level):
        members = np.zeros(1, dtype=np.int64)
        for j, fj in enumerate(choice):
            members = (
                members[:, None] + np.array(fj, dtype=np.int64)[None, :] * base**j
            ).ravel()
        members = rep(members, n)
        v = int(np.where(members <= x, mu[members], 0).sum())
        total += v * v
        count += 1
    return Fraction(total, count)


# PAIR FORMS


def squarefree(mu, x):
    idx = np.arange(1, x + 1)
    keep = mu[1 : x + 1] != 0
    return idx[keep], mu[1 : x + 1][keep]


def pair_translates(base, digs, level, mu, x):
    n = base**level
    d = level_set(base, digs, level)
    r = np.bincount(((d[:, None] - d[None, :]) % n).ravel(), minlength=n)
    ns, ms = squarefree(mu, x)
    total = 0
    for i0 in range(0, len(ns), 256):
        a, sa = ns[i0 : i0 + 256], ms[i0 : i0 + 256]
        total += int(
            (sa[:, None] * ms[None, :] * r[(a[:, None] - ns[None, :]) % n]).sum()
        )
    return Fraction(total, n)


def pair_shifts(base, digs, level, mu, x):
    n = base**level
    rf = np.zeros(base, dtype=np.int64)
    for f in digs:
        for g in digs:
            rf[(f - g) % base] += 1
    ns, ms = squarefree(mu, x)
    dg = digits(ns % n, base, level)
    total = 0
    for i0 in range(0, len(ns), 128):
        delta = (dg[i0 : i0 + 128, None, :] - dg[None, :, :]) % base
        p = rf[delta].prod(axis=2)
        total += int((ms[i0 : i0 + 128, None] * ms[None, :] * p).sum())
    return Fraction(total, n)


def pair_subsets(base, fill, level, mu, x):
    n = base**level
    ns, ms = squarefree(mu, x)
    dg = digits(ns % n, base, level)
    total = 0
    for i0 in range(0, len(ns), 128):
        dist = (dg[i0 : i0 + 128, None, :] != dg[None, :, :]).sum(axis=2)
        p = fill**level * (fill - 1) ** dist * (base - 1) ** (level - dist)
        total += int((ms[i0 : i0 + 128, None] * ms[None, :] * p).sum())
    return Fraction(total, n * (base - 1) ** level)


# FOURIER FORMS


def omega_mul(a, b, k):
    if k % 3 == 0:
        return a, b
    if k % 3 == 1:
        return -b, a - b
    return b - a, -a


def vilenkin3(g, level):
    a = g.reshape((3,) * level).astype(np.int64)
    b = np.zeros_like(a)
    for ax in range(level):
        a0, a1, a2 = (np.take(a, i, axis=ax) for i in range(3))
        b0, b1, b2 = (np.take(b, i, axis=ax) for i in range(3))
        outs = []
        for k in range(3):
            ra, rb = a0.copy(), b0.copy()
            for v, (pa, pb) in ((1, (a1, b1)), (2, (a2, b2))):
                qa, qb = omega_mul(pa, pb, k * v)
                ra, rb = ra + qa, rb + qb
            outs.append((ra, rb))
        a = np.stack([o[0] for o in outs], axis=ax)
        b = np.stack([o[1] for o in outs], axis=ax)
    return a, b


def weight_count(level):
    w = np.zeros((3,) * level, dtype=np.int64)
    for ax in range(level):
        shape = [1] * level
        shape[ax] = 3
        w = w + (np.arange(3) != 0).reshape(shape)
    return w


def top_position(level):
    t = np.full((3,) * level, -1, dtype=np.int64)
    for ax in range(level):
        shape = [1] * level
        shape[ax] = 3
        nz = np.broadcast_to((np.arange(3) != 0).reshape(shape), t.shape)
        pos = level - 1 - ax
        t = np.where(nz & (t < pos), pos, t)
    return t


def fourier_digits3(mu, x, level):
    n = 3**level
    a, b = vilenkin3(residue_vector(mu, x, n), level)
    norm = a * a - a * b + b * b
    w = weight_count(level)
    num = sum(
        int(norm[w == k].sum(dtype=np.uint64)) * 4 ** (level - k)
        for k in range(level + 1)
    )
    return Fraction(num, 9**level), norm, w


def fourier_translates(base, digs, level, mu, x):
    n = base**level
    ind = np.zeros(n)
    ind[level_set(base, digs, level)] = 1
    hd = np.abs(np.fft.fft(ind)) ** 2
    hs = np.abs(np.fft.fft(residue_vector(mu, x, n).astype(float))) ** 2
    return float((hd * hs).sum()) / n**2, hd, hs


def twists3(mu, x):
    idx = np.arange(1, x + 1)
    cls = [int(mu[1 : x + 1][idx % 3 == r].sum()) for r in range(3)]
    s = []
    for c in range(3):
        ca = cb = 0
        for r in range(3):
            pa, pb = omega_mul(cls[r], 0, c * r)
            ca, cb = ca + pa, cb + pb
        s.append(ca * ca - ca * cb + cb * cb)
    return cls, s


# VERBS


def check():
    t0 = time.time()
    mu = mobius(3**CHECK_TOP + 1)
    q = Fraction(1, 4)
    fails = 0
    print("flat mean, base 3, fill 2")
    for level in range(1, CHECK_TOP + 1):
        n = 3**level
        for digs in ((0, 1), (0, 2), (1, 2)):
            d = level_set(3, digs, level)
            cnt = np.zeros(n + 1, dtype=np.int64)
            for s in range(n):
                np.add.at(cnt, rep(d + s, n), 1)
            fails += int(np.any(cnt[1:] != 2**level))
        fd = digits(level_set(3, (0, 1), level), 3, level)
        cnt = np.zeros(n + 1, dtype=np.int64)
        for c in digits(np.arange(n), 3, level):
            np.add.at(cnt, rep(((fd + c) % 3) @ (3 ** np.arange(level)), n), 1)
        fails += int(np.any(cnt[1:] != 2**level))
    print(
        f"  translates of {{0,1}}, {{0,2}}, {{1,2}} and digit shifts, every n <= 3^level, level <= {CHECK_TOP}: {fails} counts off 2^level"
    )
    d2 = level_set(3, (0, 1), 2)
    print(
        f"  integer translates without wrap, level 2: n = 1 in {int((d2 <= 1).sum())} of 9, n = 8 in {int((d2 <= 8).sum())} of 9"
    )

    def in_design(m, digs):
        while m:
            if m % 3 not in digs:
                return False
            m //= 3
        return True

    subs = list(combinations(range(3), 2))
    print(
        f"  one digit set at every position: n = 1 in {sum(in_design(1, f) for f in subs)} of 3, n = 15 in {sum(in_design(15, f) for f in subs)} of 3"
    )
    print("level x | E translates {0,1} | E digit shifts | M(x) | ok")
    jobs = [(lv, 3**lv) for lv in range(1, CHECK_TOP + 1)] + MIDS
    for level, x in jobs:
        n = 3**level
        et = translates(3, (0, 1), level, mu, x)
        pt = pair_translates(3, (0, 1), level, mu, x)
        ft, _hd, _hs = fourier_translates(3, (0, 1), level, mu, x)
        et2 = translates(3, (0, 2), level, mu, x)
        pt2 = pair_translates(3, (0, 2), level, mu, x)
        ev = shifts(3, (0, 1), level, mu, x)
        pv = pair_shifts(3, (0, 1), level, mu, x)
        er = subsets(3, 2, level, mu, x)
        pr = pair_subsets(3, 2, level, mu, x)
        fv, _, _ = fourier_digits3(mu, x, level)
        m = int(mu[1 : x + 1].sum())
        cls, s = twists3(mu, x)
        mean2 = Fraction(4, 9) ** level
        low = mean2 * (s[0] + q * (s[1] + s[2]))
        idx = np.arange(1, x + 1)
        sq9 = sum(int(mu[1 : x + 1][idx % 9 == r].sum()) ** 2 for r in range(9))
        hf = lambda t: abs(1 + np.exp(2j * np.pi * t)) ** 2
        kappa2 = min(hf(c / 9) * hf(3 * c / 9) for c in range(9)) / 16
        nine = level < 2 or (
            q * q * mean2 * 9 * sq9 <= ev
            and kappa2 * float(mean2) * 9 * sq9 <= float(et) * (1 + 1e-12)
        )
        ok = (
            et == pt
            and et2 == pt2
            and ev == pv == pr == fv == er
            and abs(ft - float(et)) <= 1e-9 * float(et)
            and s[0] == m * m
            and s[0] + s[1] + s[2] == 3 * sum(c * c for c in cls)
            and mean2 * m * m <= min(et, et2, ev)
            and low <= min(et, et2, ev)
            and q * mean2 * 3 * sum(c * c for c in cls) <= min(et, ev)
            and nine
        )
        fails += int(not ok)
        print(f"{level} {x} | {et} | {ev} | {m} | {ok}")
    print("other bases, pair form against direct average")
    for base, fill, level in SMALL:
        n = base**level
        x = n
        mub = mobius(n + 1)
        digs = tuple(range(fill)) if base != 4 else (0, 2)
        et, pt = (
            translates(base, digs, level, mub, x),
            pair_translates(base, digs, level, mub, x),
        )
        ev, pv = (
            shifts(base, digs, level, mub, x),
            pair_shifts(base, digs, level, mub, x),
        )
        er, pr = (
            subsets(base, fill, level, mub, x),
            pair_subsets(base, fill, level, mub, x),
        )
        ft, _, _ = fourier_translates(base, digs, level, mub, x)
        ok = (
            et == pt
            and ev == pv
            and er == pr
            and abs(ft - float(et)) <= 1e-9 * float(et)
        )
        fails += int(not ok)
        print(f"  base {base} digits {digs} level {level}: T {et} V {ev} R {er} {ok}")
    print(f"{fails} failures, {time.time() - t0:.1f} s")


def read():
    t0 = time.time()
    mu = mobius(3**READ_TOP + 1)
    alpha = np.log(2) / np.log(3)
    print(f"half plane edge 1 - alpha/2 = {1 - alpha / 2:.9f}")
    print("level | E_R/mass | E_T/mass | a=0 share R | low share R | low share T")
    for level in range(4, READ_TOP + 1):
        x = 3**level
        e, norm, w = fourier_digits3(mu, x, level)
        mass = (2 / 3) ** level * x
        m = int(mu[1 : x + 1].sum())
        top = top_position(level)
        k = int(np.floor((1 - alpha) * level))
        wts = (4.0 ** (level - w)) * norm / 9.0**level
        low_r = float(wts[top < k].sum()) / float(e)
        ft, hd, hs = fourier_translates(3, (0, 1), level, mu, x)
        a = np.arange(3**level)
        low_t = float(((hd * hs)[a % 3 ** (level - k) == 0]).sum()) / 9.0**level / ft
        share = float(Fraction(4, 9) ** level * m * m / e)
        print(
            f"{level} | {float(e) / mass:.6f} | {ft / mass:.6f} | {share:.6f} | {low_r:.6f} (k={k}) | {low_t:.6f}"
        )
    print(f"{time.time() - t0:.1f} s")


if __name__ == "__main__":
    verbs = {"check": check, "read": read}
    for name in sys.argv[1:] or list(verbs):
        verbs[name]()
