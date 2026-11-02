import itertools
import subprocess
import sys
import time

import mpmath
import numpy as np
import sympy

NB = 1.45
TWO_PI = 2 * np.pi
CONTRASTS = (1.6, 2.3, 3.5)
DEEP = 0.02


def lay(n, ph):
    c, s = np.cos(ph), np.sin(ph)
    return (c + 0j, 1j * s / n, 1j * n * s, c + 0j)


def mul(p, q):
    a, b, c, d = p
    e, f, g, h = q
    return (a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)


def num_den(m):
    a, b, c, d = m
    x, y = NB * a + NB * NB * b, c + NB * d
    return x - y, x + y


def peak(m):
    return np.maximum.reduce([np.abs(z) for z in m])


def levels(base, digits, top, delta, na):
    m = lay(na, delta)
    dead = np.zeros(delta.shape, bool)
    out = [(m, dead.copy())]
    for k in range(top):
        spacer = lay(NB, np.mod(float(base) ** k * delta, TWO_PI))
        c = np.where(dead, 1.0, np.maximum(peak(m), 1e10) / 1e10)
        unit = tuple(z / c for z in m)
        n = None
        for d in range(base):
            x = unit if d in digits else spacer
            n = x if n is None else mul(n, x)
        big = peak(n)
        grow = len(digits) * np.log(c)
        dead |= np.log(big) + grow >= np.log(1e100)
        scale = np.where(dead, 1 / big, np.exp(np.minimum(grow, 700)))
        m = tuple(z * scale for z in n)
        out.append((m, dead.copy()))
    return out


def log_rt(m, dead):
    num, _ = num_den(m)
    v = np.log(np.maximum(np.abs(num) / (2 * NB), 1e-14))
    return np.where(dead, np.inf, v)


def phase_t(m):
    _, den = num_den(m)
    return -np.angle(den)


def word(base, digits, level):
    w = [1]
    for _ in range(level):
        w = [
            x if d in digits else 0
            for d in range(base)
            for x in (w if d in digits else [0] * len(w))
        ]
    return w


def brute(base, digits, level, delta, na):
    mpmath.mp.dps = 40
    out = []
    for dl in delta:
        dl = mpmath.mpf(float(dl))
        c, s = mpmath.cos(dl), mpmath.sin(dl)
        mats = {}
        for n in (mpmath.mpf(na), mpmath.mpf(NB)):
            mats[n] = (c, 1j * s / n, 1j * n * s, c)
        m = (1, 0, 0, 1)
        for x in word(base, digits, level):
            m = mul(m, mats[mpmath.mpf(na) if x else mpmath.mpf(NB)])
        out.append([complex(z) for z in m])
    return tuple(np.array(z) for z in zip(*out))


def mahler_p(digits):
    c = np.zeros(max(digits) + 1)
    c[list(digits)] = 1
    r = np.roots(c[::-1])
    return float(np.sum(np.log(np.maximum(1, np.abs(r)))))


def mahler_q(digits, n=1 << 15):
    y = np.exp(2j * np.pi * (np.arange(n) + 0.5) / n)
    deg = len(digits) - 1
    if deg == 0:
        return 0.0
    coef = np.stack([y ** (d - i) for i, d in enumerate(digits)], axis=1)
    lead = coef[:, -1]
    comp = np.zeros((n, deg, deg), complex)
    comp[:, 0, :] = -coef[:, -2::-1] / lead[:, None]
    if deg > 1:
        comp[:, np.arange(1, deg), np.arange(deg - 1)] = 1
    r = np.linalg.eigvals(comp)
    return float(np.mean(np.sum(np.log(np.maximum(1, np.abs(r))), axis=1)))


def q_value(digits, x, y):
    return sum(x**i * y ** (d - i) for i, d in enumerate(digits))


X, Y, Z = sympy.symbols("x y z")


def cyclotomic_one(coeffs):
    poly = sympy.Poly(sum(c * Z**e for e, c in coeffs), Z)
    return poly.is_cyclotomic or (-poly).is_cyclotomic


def boyd_class(digits):
    q = sum(X**i * Y ** (d - i) for i, d in enumerate(digits))
    _, fs = sympy.factor_list(q, X, Y)
    for f, _ in fs:
        terms = sympy.Poly(f, X, Y).terms()
        if len(terms) == 1:
            continue
        e0 = np.array(terms[0][0])
        diffs = [np.array(t[0]) - e0 for t in terms]
        v = next(d for d in diffs if d.any())
        u = v // np.gcd(abs(int(v[0])), abs(int(v[1])))
        steps = []
        for d in diffs:
            if u[0] * d[1] - u[1] * d[0] != 0:
                return False
            steps.append(int(d @ u) // int(u @ u))
        lo = min(steps)
        if not cyclotomic_one([(s - lo, int(c)) for s, (_, c) in zip(steps, terms)]):
            return False
    return True


def p_cyclotomic(digits):
    _, fs = sympy.factor_list(sum(Z**d for d in digits), Z)
    return all(sympy.Poly(f, Z).is_cyclotomic for f, _ in fs)


def fmt(xs, p=3):
    return " ".join(f"{x:.{p}f}" for x in xs)


def check():
    rng = np.random.default_rng(1)
    delta = rng.uniform(0, np.pi, 8)
    worst = 0.0
    for base, digits, level in [
        (3, {0, 2}, 6),
        (4, {0, 1, 3}, 5),
        (5, {0, 1, 2, 4}, 4),
        (7, {0, 2, 3, 4, 6}, 4),
        (9, {0, 1, 4, 7, 8}, 3),
    ]:
        err = 0.0
        for na in CONTRASTS:
            a = brute(base, digits, level, delta, na)
            b = levels(base, digits, level, delta, na)[level][0]
            scale = np.maximum.reduce([np.abs(z) for z in a])
            err = max(
                err, max(float(np.max(np.abs(x - y) / scale)) for x, y in zip(a, b))
            )
        worst = max(worst, err)
        print(
            f"check b={base} D={sorted(digits)} level={level} cells={base**level} max rel err {err:.1e}"
        )
    assert worst < 1e-12
    h = 1e-8
    for base, digits in [(4, [0, 1, 3]), (7, [0, 2, 3, 4, 6]), (3, [0, 2])]:
        delta = rng.uniform(0, np.pi, 64)
        lv = levels(base, set(digits), 4, delta, NB * (1 + h))
        r0 = num_den(lv[0][0])[0]
        errs = []
        for k in range(1, 5):
            born = np.prod(
                [
                    q_value(
                        digits,
                        np.exp(-2j * (base**j) * delta),
                        np.exp(-2j * (base**j) * delta),
                    )
                    for j in range(k)
                ],
                axis=0,
            )
            got = num_den(lv[k][0])[0] / r0
            errs.append(
                float(
                    np.max(
                        np.abs(np.abs(got) - np.abs(born)) / np.maximum(1, np.abs(born))
                    )
                )
            )
        print(
            f"born derivative b={base} D={digits} levels 1-4 max rel err at h=1e-8: {fmt(errs, 9)}"
        )
        assert max(errs) < 1e-5


def born():
    rng = np.random.default_rng(2)
    delta = rng.uniform(0, np.pi, 200000)
    na = 1.455
    gaps, ses = [], []
    for base, digits, top in [
        (4, [0, 1, 3], 7),
        (5, [0, 1, 2, 4], 6),
        (7, [0, 1, 5], 6),
        (3, [0, 2], 8),
        (6, [0, 1, 3, 4], 6),
        (7, [0, 2, 3, 4, 6], 6),
    ]:
        lv = levels(base, set(digits), top, delta, na)
        lr = [log_rt(*x) for x in lv]
        steps = [lr[k + 1] - lr[k] for k in range(top)]
        per = [float(np.mean(s)) for s in steps]
        se = float(np.std(lr[top] - lr[0]) / np.sqrt(delta.size) / top)
        slope = float(np.mean(lr[top] - lr[0]) / top)
        mp = mahler_p(digits)
        gaps.append(abs(slope - mp))
        ses.append(se)
        print(
            f"born nA={na} b={base} D={digits} m(P)={mp:.5f} m(Q)={mahler_q(digits):.5f} per-level {fmt(per, 4)} mean {slope:.4f} se {se:.4f} |mean - m(P)| {abs(slope - mp):.5f}"
        )
    print(
        f"born: |mean - m(P)| at most {max(gaps[:3]):.5f} on the three with m(P) > 0 and {max(gaps):.5f} on all six designs, se at most {max(ses):.4f}"
    )


def deep_run(base, digits, top, delta, na):
    lv = levels(base, set(digits), top, delta, na)
    lr = [log_rt(*x) for x in lv]
    r = {
        key: []
        for key in ("drift", "se", "lift", "lemma", "coup", "torus", "tse", "dead")
    }
    stop = False
    for k in range(1, top):
        a_all = phase_t(lv[k][0])
        b_all = phase_t(lay(NB, np.mod(float(base) ** k * delta, TWO_PI)))
        lt = np.log(
            np.maximum(
                np.abs(q_value(digits, np.exp(2j * a_all), np.exp(2j * b_all))), 1e-300
            )
        )
        r["torus"].append(float(np.mean(lt)))
        r["tse"].append(float(np.std(lt) / np.sqrt(lt.size)))
        r["dead"].append(float(np.mean(lv[k][1])))
        sel = (lr[k] < np.log(DEEP)) & np.isfinite(lr[k + 1])
        stop = stop or sel.sum() < 300
        if stop:
            continue
        al = phase_t(lv[k][0])[sel]
        be = phase_t(lay(NB, np.mod(float(base) ** k * delta[sel], TWO_PI)))
        lq = np.log(np.abs(q_value(digits, np.exp(2j * al), np.exp(2j * be))))
        dd = (lr[k + 1] - lr[k])[sel]
        r["drift"].append(float(np.mean(dd)))
        r["se"].append(float(np.std(dd) / np.sqrt(dd.size)))
        r["lift"].append(float(np.mean(lq)))
        tiny = lr[k][sel] < np.log(1e-3)
        if tiny.sum() >= 50:
            r["lemma"].append(float(np.median(np.abs(dd - lq)[tiny])))
        r["coup"].append(float(np.abs(np.mean(np.exp(2j * (al - be))))))
    return r


DEEP_DESIGNS = [
    (4, [0, 1, 3], 12),
    (5, [0, 1, 2, 4], 10),
    (7, [0, 1, 5], 10),
    (7, [0, 2, 3, 4, 6], 10),
    (9, [0, 1, 4, 7, 8], 9),
    (8, [0, 2, 3, 4, 5, 7], 9),
    (3, [0, 2], 14),
    (7, [0, 1, 3, 4], 10),
    (7, [0, 2, 4, 6], 10),
]

TABLE = DEEP_DESIGNS[2:]


def deep():
    rng = np.random.default_rng(3)
    delta = rng.uniform(0, np.pi, 200000)
    cells = lemmas = 0
    worst = track = 0.0
    zero, pos = [], []
    for base, digits, top in DEEP_DESIGNS:
        mq = mahler_q(digits)
        print(f"deep b={base} D={digits} m(P)={mahler_p(digits):.4f} m(Q)={mq:.4f}")
        for na in CONTRASTS:
            r = deep_run(base, digits, top, delta, na)
            tail, tse = r["drift"][-4:], r["se"][-4:]
            mean, se = float(np.mean(tail)), float(np.mean(tse))
            gap = max(abs(a - b) for a, b in zip(r["drift"][-4:], r["lift"][-4:]))
            cells += 1
            track = max(track, gap)
            lemmas += bool(r["lemma"])
            worst = max([worst] + r["lemma"])
            law = (
                f"{max(r['lemma']):.1e}"
                if r["lemma"]
                else "n/a, fewer than 50 frequencies below 1e-3"
            )
            print(
                f"  nA={na} deep drift {fmt(r['drift'])} | tail mean {mean:.3f} se {se:.3f} | (tail - m(Q))/se {(mean - mq) / se:.1f}"
            )
            print(
                f"    per-level se {fmt(r['se'])} | deep lift {fmt(r['lift'])} | tail max |drift - lift| {gap:.4f}"
            )
            print(f"    step law median err {law} | coupling {fmt(r['coup'], 2)}")
            print(
                f"    torus mean of ln|Q_D| over all frequencies {fmt(r['torus'])} | last {r['torus'][-1]:.3f} se {r['tse'][-1]:.3f} | dark fraction {fmt(r['dead'], 2)}"
            )
            if (base, digits, top) in TABLE:
                (zero if mq < 1e-6 else pos).append((mean, se))
    print(
        f"deep step law: {lemmas} of {cells} design-contrast cells have 50 or more frequencies below 1e-3, worst median err {worst:.1e}; deep set, last four levels: |drift - lift| at most {track:.4f}"
    )
    print(
        f"deep zero set: m(Q) = 0 tail |mean| at most {max(abs(m) for m, _ in zero):.3f} (se at most {max(e for _, e in zero):.3f}); m(Q) > 0 tail mean at least {min(m for m, _ in pos):.3f} (se at most {max(e for _, e in pos):.3f})"
    )


def smyth():
    s = subprocess.run(
        ["gp", "-q"],
        input="print(3*sqrt(3)/(4*Pi)*lfun(-3,2))",
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    mq = mahler_q([0, 1, 3], 1 << 17)
    print(
        f"smyth constant from PARI {s[:14]}  m(Q_(0,1,3)) numeric {mq:.10f}  m(P_(0,1,3)) {mahler_p([0, 1, 3]):.10f}"
    )
    assert abs(mq - float(s)) < 1e-8


def passband(base, digits, top, delta, na):
    lv = levels(base, set(digits), top, delta, na)
    return np.array([float(np.mean(log_rt(*x) < 0)) for x in lv])


def refute():
    rng = np.random.default_rng(4)
    delta = rng.uniform(0, np.pi, 200000)
    for base, digits in [
        (7, [0, 2, 3, 4, 6]),
        (9, [0, 1, 4, 7, 8]),
        (8, [0, 2, 3, 4, 5, 7]),
        (7, [0, 2, 4, 6]),
        (7, [0, 1, 3, 4]),
    ]:
        f = passband(base, digits, 12, delta, 2.3)
        ks = np.arange(1, 13)
        pc, bc = p_cyclotomic(digits), boyd_class(digits)
        print(
            f"refute b={base} D={digits} P cyclotomic {pc} m(P)={mahler_p(digits):.1e} Q generalised cyclotomic {bc} m(Q)={mahler_q(digits):.4f}"
        )
        print(f"  f_k sqrt(k) k=1..12 at nA=2.3: {fmt(f[1:] * np.sqrt(ks))}")
        if digits in ([0, 2, 3, 4, 6], [0, 1, 4, 7, 8], [0, 2, 3, 4, 5, 7]):
            assert pc and not bc


def designs(lo, hi):
    seen = set()
    for base in range(lo, hi + 1):
        for size in range(2, base):
            for rest in itertools.combinations(range(1, base), size - 1):
                d = (0,) + rest
                mirror = tuple(sorted(max(d) - x for x in d))
                key = (base, min(d, mirror))
                if key not in seen:
                    seen.add(key)
                    yield key


def census():
    rng = np.random.default_rng(5)
    delta = rng.uniform(0, np.pi, 100000)
    top, half = 12, 6
    agree = total = 0
    flats = {True: [], False: []}
    for base, digits in designs(3, 6):
        digits = list(digits)
        bc = boyd_class(digits)
        mq = mahler_q(digits)
        f = passband(base, digits, top, delta, 2.3)
        flat = f[top] * np.sqrt(top) / (f[half] * np.sqrt(half)) if f[half] > 0 else 0.0
        critical = flat > 0.9
        total += 1
        agree += critical == bc
        flats[bc].append(flat)
        tag = "critical" if critical else "decays"
        print(
            f"census b={base} D={digits} Q gen-cyclotomic {bc!s:5} m(Q)={mq:.4f} m(P)={mahler_p(digits):.4f} f6={f[half]:.4f} f12={f[top]:.5f} flatness {flat:.3f} {tag}"
        )
    print(
        f"census designs {total} (bases 3 to 6, up to shift and mirror), flatness verdict matches the Boyd class on {agree}"
    )
    for bc, name in ((True, "Boyd-cyclotomic"), (False, "other")):
        print(
            f"census {name}: {len(flats[bc])} designs, flatness {min(flats[bc]):.3f} to {max(flats[bc]):.3f}"
        )


VERBS = {
    "check": check,
    "born": born,
    "smyth": smyth,
    "deep": deep,
    "refute": refute,
    "census": census,
}

if __name__ == "__main__":
    import warnings

    warnings.filterwarnings("ignore")
    for v in sys.argv[1:] or list(VERBS):
        t0 = time.time()
        VERBS[v]()
        print(f"[{v} {time.time() - t0:.1f}s]")
