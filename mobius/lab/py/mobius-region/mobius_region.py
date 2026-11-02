import itertools
import math
import sys
from fractions import Fraction

import numpy as np
from mpmath import mp

# DESIGNS


def missing(q, a0):
    return q, tuple(d for d in range(q) if d != a0)


def proper_sets(q):
    out = []
    for m in range(2, q):
        for F in itertools.combinations(range(q), m):
            out.append((q, F))
    return out


def census():
    rows = []
    for q in (3, 4, 5):
        rows.extend(proper_sets(q))
    for a0 in range(10):
        rows.append(missing(10, a0))
    rows.append(missing(21, 0))
    return rows


DEPTH = {3: (9, 8), 4: (7, 8), 5: (6, 8), 10: (5, 8), 21: (5, 8)}

TGRID = [
    Fraction(1),
    Fraction(11, 10),
    Fraction(6, 5),
    Fraction(13, 10),
    Fraction(7, 5),
    Fraction(235, 154),
    Fraction(3, 2),
    Fraction(8, 5),
    Fraction(17, 10),
    Fraction(9, 5),
    Fraction(19, 10),
    Fraction(39, 20),
]

# SAFE ROUNDING


def fdown(x, n=7):
    return math.floor(x * 10**n) / 10**n


def fup(x, n=7):
    return math.ceil(x * 10**n) / 10**n


FAILS = []


def want(cond, name):
    if not cond:
        FAILS.append(name)
    return cond


# THE MASS EXPONENT


def alpha_bracket(q, k, den=10**5):
    mp.dps = 40
    a = int(mp.floor(den * mp.log(k) / mp.log(q)))
    kb = k**den
    while q ** (a + 1) <= kb:
        a += 1
    while q**a > kb:
        a -= 1
    return Fraction(a, den), Fraction(a + 1, den)


# THE WINDOW FACTORS


def differences(F):
    k = len(F)
    mult = {}
    for f1 in F:
        for f2 in F:
            d = abs(f1 - f2)
            mult[d] = mult.get(d, 0) + 1
    return sorted(mult.items()), k


def window_factors(q, F, nd, m, chunk=1 << 17):
    dif, k = differences(F)
    lip = 2.0 * math.pi * sum(F) / k
    slack = lip / (2.0 * m * q**nd)
    tol = 1e-12
    total = q**nd
    gup = np.empty(total)
    glo = np.empty(total)
    off = (np.arange(m) + 0.5) / (m * q**nd)
    for lo in range(0, total, chunk):
        hi = min(lo + chunk, total)
        t = np.arange(lo, hi, dtype=np.float64)[:, None] / q**nd + off[None, :]
        acc = np.zeros_like(t)
        for d, c in dif:
            if d == 0:
                acc += c
            else:
                acc += c * np.cos(2.0 * math.pi * d * t)
        acc /= k * k
        gup[lo:hi] = np.sqrt(np.clip(acc + tol, 0.0, None)).max(axis=1) + slack
        glo[lo:hi] = np.clip(
            np.sqrt(np.clip(acc - tol, 0.0, None)).min(axis=1) - slack, 0.0, None
        )
    return gup, glo


# THE TRANSFER ROOT


def roots(G, q, nd, power=1.0, iters=400):
    S = q ** (nd - 1)
    W = (G.reshape(S, q) ** power) if power != 1.0 else G.reshape(S, q)
    tgt = ((np.arange(S, dtype=np.int64)[:, None] * q + np.arange(q)) % S).astype(
        np.int32
    )
    y = np.ones(S)
    for _ in range(iters):
        z = (W * y[tgt]).sum(axis=1)
        top = z.max()
        if top <= 0.0:
            return 0.0, 0.0
        y = z / top + 1e-30
    z = (W * y[tgt]).sum(axis=1)
    r = z / y
    return float(r.min()) * (1 - 1e-12), float(r.max()) * (1 + 1e-12)


def exponent(mu, q):
    if mu <= 0.0:
        return None
    return math.log(mu) / math.log(q)


# THE THRESHOLD FROM BELOW


def beta_lower(glo, q, nd, ahi, target=0.25, tail=1.99, cap_cells=600):
    tailbound = (1.0 - ahi) / (2.0 - tail)
    work = [(1.0, tail)]
    done = []
    cells = 0
    while work:
        t0, t1 = work.pop()
        mu, _ = roots(glo, q, nd, power=t1)
        e = exponent(mu, q)
        cells += 1
        b = -1.0 if e is None else e / (2.0 - t0)
        if b > target or cells > cap_cells:
            done.append((t0, t1, b))
            if b <= target:
                return None, cells, tailbound, (t0, t1, b)
        else:
            mid = 0.5 * (t0 + t1)
            if mid - t0 < 1e-4:
                return None, cells, tailbound, (t0, t1, b)
            work.append((t0, mid))
            work.append((mid, t1))
    return min(min(r[2] for r in done), tailbound), cells, tailbound, None


# THE THREE PARAMETERS


def params(q, F, nd=None, m=None, tgrid=None):
    k = len(F)
    if nd is None:
        nd, m = DEPTH.get(q, (5, 8))
    if tgrid is None:
        tgrid = TGRID if q <= 10 else TGRID[:1]
    alo, ahi = alpha_bracket(q, k)
    gup, glo = window_factors(q, F, nd, m)
    _, muup = roots(gup, q, nd)
    mulo, _ = roots(glo, q, nd)
    a1hi = exponent(muup, q)
    a1lo = exponent(mulo, q)
    floor = 1.0 - float(ahi)
    if a1lo is None or a1lo < floor:
        a1lo = floor
    mts = []
    for t in tgrid:
        _, mu = roots(gup, q, nd, power=float(t))
        mt = exponent(mu, q)
        mts.append((t, mt, mt / float(2 - t)))
    tbest, mtbest, bhi = min(mts, key=lambda r: r[2])
    return {
        "q": q,
        "F": F,
        "k": k,
        "nd": nd,
        "m": m,
        "alo": float(alo),
        "ahi": float(ahi),
        "a1lo": a1lo,
        "a1hi": a1hi,
        "bhi": bhi,
        "tbest": tbest,
        "mtbest": mtbest,
        "bfloor": 1.0 - float(ahi),
        "mts": mts,
        "gup": gup,
        "glo": glo,
    }


def name_of(q, F):
    if len(F) == q - 1:
        return "missing " + str(next(d for d in range(q) if d not in F))
    return "{" + ",".join(str(d) for d in F) + "}"


# THE CRITERION

CONDS = ["L1", "C1", "C2", "L2", "L3", "L4", "L5"]


def cap(name, alpha, a1):
    if name == "C1":
        return Fraction(1, 4)
    if name == "C2":
        return Fraction(2, 5) * (1 - a1)
    if name == "L2":
        return (1 - a1) / (2 * (2 - alpha))
    if name == "L4":
        return (1 + alpha / 2) / 5
    if name == "L5":
        return (1 - a1) * (1 - a1 + alpha / 2) / 2
    if name == "L3":
        if a1 < Fraction(1, 3):
            return None
        u = min(Fraction(1), 2 * a1 / alpha)
        return u * alpha / (4 * (3 * a1 - 1 + u * (1 - a1)))
    raise ValueError(name)


def beta_cap(alpha, a1):
    best = None
    who = None
    for name in CONDS[1:]:
        c = cap(name, alpha, a1)
        if c is None:
            continue
        if best is None or c < best:
            best, who = c, name
    return best, who


def verdict(alpha, a1, beta):
    rows = []
    rows.append(("L1", 2 * a1, alpha, a1 < alpha / 2))
    for name in CONDS[1:]:
        c = cap(name, alpha, a1)
        if c is None:
            rows.append((name, beta, None, True))
        else:
            ok = beta <= c if name in ("C1", "C2") else beta < c
            rows.append((name, beta, c, ok))
    ok = all(r[3] for r in rows) and a1 < Fraction(1, 2)
    binder = None
    slack = None
    for name, left, right, good in rows:
        if right is None:
            continue
        s = float(right) - float(left)
        if slack is None or s < slack:
            slack, binder = s, name
    return rows, ok, binder, slack


# VERB PARAMS


def verb_params(argv):
    q = int(argv[0])
    F = tuple(int(c, 36) for c in argv[1])
    nd = int(argv[2]) if len(argv) > 2 else None
    m = int(argv[3]) if len(argv) > 3 else None
    p = params(q, F, nd, m)
    print(
        f"q = {p['q']}  F = {sorted(p['F'])}  k = {p['k']}  window {p['nd']} digits, sub-scan {p['m']}"
    )
    print(
        f"alpha in [{fdown(p['alo'])}, {fup(p['ahi'])}]  (integer comparison, k^b against q^a, b = 10^5)"
    )
    print(
        f"alpha_1 in [{fdown(p['a1lo'])}, {fup(p['a1hi'])}]  (lower from the inf window and the l^1 floor, upper from the sup window)"
    )
    print(
        f"beta <= {fup(p['bhi'])} at t = {p['tbest']}, m_t <= {fup(p['mtbest'])}; beta >= {fdown(p['bfloor'])} by Parseval"
    )


# VERB CRITERION


def frac_up(x, den=10**7):
    return Fraction(math.ceil(x * den), den)


def frac_down(x, den=10**7):
    return Fraction(math.floor(x * den), den)


def decide(p, betalo):
    pes = verdict(frac_down(p["alo"]), frac_up(p["a1hi"]), frac_up(p["bhi"]))
    opt = verdict(frac_up(p["ahi"]), frac_down(p["a1lo"]), frac_down(betalo))
    holds = pes[1]
    dead = [r[0] for r in opt[0] if not r[3]]
    if not (frac_down(p["a1lo"]) < Fraction(1, 2)):
        dead.append("LAT")
    return holds, dead, pes


def verb_criterion(_argv):
    print(
        "q  design         k  alpha        alpha_1 in                beta in                   verdict   binds                     GRH theta <="
    )
    out = []
    for q, F in census():
        p = params(q, F)
        betalo = 1.0 - float(p["ahi"])
        cert = ""
        if betalo <= 0.25 and 2.0 * p["a1lo"] < float(p["ahi"]):
            lo, cells, _, _ = beta_lower(p["glo"], q, p["nd"], float(p["ahi"]))
            if lo is not None:
                betalo, cert = lo, f" ({cells} cells)"
        holds, dead, pes = decide(p, betalo)
        tag = "HOLDS" if holds else ("REFUTED" if dead else "open")
        grh = f"{fup(max(1.0 - (0.25 - p['a1hi']) / p['alo'], 1.0 - (0.25 - p['a1hi']) / p['ahi']))}"
        binds = ",".join(dead) if dead else (pes[2] or "-")
        print(
            f"{q:<3}{name_of(q, F):<15}{p['k']:<3}{fdown(p['alo'], 6):<13}"
            f"[{fdown(p['a1lo'], 6)}, {fup(p['a1hi'], 6)}]{'':<6}"
            f"[{fdown(betalo, 6)}, {fup(p['bhi'], 6)}]{'':<5}{tag:<10}{binds:<26}{grh}{cert}"
        )
        out.append((q, F, p, tag, dead))
    print()
    print(
        "Every alpha and every lower bound truncates down, every upper bound rounds up."
    )
    print(
        "HOLDS uses the pessimistic corner (alpha low, alpha_1 high, beta high); REFUTED uses the optimistic one; open is neither."
    )
    print(
        "A beta lower bound is 1 - alpha by Parseval, or the certified cell chain where that does not decide."
    )
    for tag in ("HOLDS", "REFUTED", "open"):
        rows = [r for r in out if r[3] == tag]
        print(
            f"{tag}: {len(rows)} of {len(out)} designs"
            + (
                ""
                if tag != "HOLDS"
                else ": " + ", ".join(f"q={r[0]} {name_of(r[0], r[1])}" for r in rows)
            )
        )
    counts = {}
    for _, _, _, tag, dead in out:
        if tag == "REFUTED":
            counts[",".join(dead)] = counts.get(",".join(dead), 0) + 1
    for key in sorted(counts):
        print(f"  refuted on {key}: {counts[key]} designs")
    best = min((r for r in out if r[0] == 10), key=lambda r: r[2]["bhi"])
    worst = max((r for r in out if r[0] == 10), key=lambda r: r[2]["bhi"])
    print(
        f"base 10 cheapest column {name_of(10, best[1])}: beta <= {fup(best[2]['bhi'])}, over the bar by {fdown(best[2]['bhi'] - 0.25)}"
    )
    print(
        f"base 10 dearest column {name_of(10, worst[1])}: beta <= {fup(worst[2]['bhi'])}, over the bar by {fdown(worst[2]['bhi'] - 0.25)}"
    )


# VERB THRESHOLD


def verb_threshold(argv):
    target = 0.25
    if argv and len(argv) == 1:
        target = float(Fraction(argv[0]))
        argv = []
    elif len(argv) >= 3:
        target = float(Fraction(argv[2]))
    fams = (
        [(10, tuple(d for d in range(10) if d != a0)) for a0 in range(10)]
        if not argv
        else [(int(argv[0]), tuple(int(c, 36) for c in argv[1]))]
    )
    print(
        "The exceptional-set threshold from below: m_t is non-increasing in t and m_2 = 1 - alpha exactly,"
    )
    print(
        "so on a cell [t0, t1] every t has m_t/(2 - t) >= m_(t1)/(2 - t0), and above the cut the Parseval value alone decides."
    )
    print(f"target {fdown(target)}")
    print(
        "design            alpha_1 <=   beta <=      beta >=      cells  tail cut  verdict"
    )
    for q, F in fams:
        p = params(q, F)
        lo, cells, tb, bad = beta_lower(
            p["glo"], q, p["nd"], float(p["ahi"]), target=target
        )
        if lo is None:
            note = f"undecided at the target, t in [{bad[0]:.4f}, {bad[1]:.4f}]"
            shown = "-"
        else:
            note = (
                f"REFUTED, no admissible beta clears {fdown(target)}"
                if lo > target
                else "open"
            )
            shown = str(fdown(lo))
        print(
            f"{q} {name_of(q, F):<15} {fup(p['a1hi']):<12} {fup(p['bhi']):<12} {shown:<12} {cells:<6} {tb:<9.4f} {note}"
        )
    print()
    print(
        "A lower bound truncates down and an upper bound rounds up; the lower bound comes from the infimum window,"
    )
    print("which is a lower bound on the grid moment and so on the sup-over-shift one.")


# VERB REGION


def cross(alpha, lo, hi, name, target=Fraction(1, 4), steps=60):
    for _ in range(steps):
        mid = (lo + hi) / 2
        c = cap(name, alpha, mid)
        if c is not None and c < target:
            hi = mid
        else:
            lo = mid
    return hi


def verb_region(argv):
    alphas = [
        Fraction(3, 4),
        Fraction(4, 5),
        Fraction(9, 10),
        Fraction(19, 20),
        Fraction(98397, 10**5),
    ]
    if argv:
        alphas = [Fraction(argv[0])]
    print(
        "The region at fixed alpha, in the (alpha_1, beta) plane: the wall L1 puts a ceiling on alpha_1 and the six caps put one on beta."
    )
    print(
        "alpha        L1 wall      floor 1-alpha  C1 cut at    C2 cut at    L2 cut at    L3 cut at    L5 cut at    beta at the wall"
    )
    for alpha in alphas:
        top = min(alpha / 2, Fraction(1, 2))
        cuts = {}
        for name in ("C2", "L2", "L3", "L5"):
            ct = cap(name, alpha, top)
            cuts[name] = (
                "never"
                if (ct is None or ct >= Fraction(1, 4))
                else fdown(float(cross(alpha, Fraction(0), top, name)), 6)
            )
        wall, who = beta_cap(alpha, top)
        print(
            f"{fdown(float(alpha), 6):<12} {fdown(float(top), 6):<12} {fdown(1 - float(alpha), 6):<14} "
            f"{'0':<12} {cuts['C2']!s:<12} {cuts['L2']!s:<12} {cuts['L3']!s:<12} {cuts['L5']!s:<12} {fdown(float(wall), 6)} by {who}"
        )
    print()
    print("The boundary itself, printed at nine stations from the origin to the wall:")
    for alpha in alphas:
        top = min(alpha / 2, Fraction(1, 2))
        print(f"alpha = {fdown(float(alpha), 6)}")
        for i in range(9):
            a1 = top * i / 8
            c, who = beta_cap(alpha, a1)
            reach = "design" if a1 >= 1 - alpha and c >= 1 - alpha else "-"
            print(
                f"  alpha_1 {fdown(float(a1), 6):<11} beta <= {fdown(float(c), 6):<11} {who:<4} {reach}"
            )
    print()
    print(
        "A cap is an upper limit on beta and prints truncated down; a design also obeys alpha_1 >= 1 - alpha and beta >= 1 - alpha."
    )


# VERB BOUNDARY


def verb_boundary(_argv):
    print(
        "Which of the seven inequalities can bind. Each cap is monotone in alpha_1, so its minimum over the region sits at the wall alpha_1 = alpha/2."
    )
    print(
        "alpha        C1       C2 at wall   L2 at wall   L3 at wall   L4          L5 at wall"
    )
    for na in (68, 70, 75, 80, 85, 90, 95, 98397 / 1000, 99, 999 / 10):
        alpha = Fraction(round(na * 1000), 100000)
        top = min(alpha / 2, Fraction(1, 2))
        line = f"{fdown(float(alpha), 6):<12} {'0.25':<8}"
        for name in ("C2", "L2", "L3", "L4", "L5"):
            c = cap(name, alpha, top)
            line += f" {fdown(float(c), 6):<12}"
        print(line)
    print()
    print(
        "L2 and L3 read exactly 1/4 at the wall at every alpha, L5 reads (2 - alpha)/4 and L4 reads (1 + alpha/2)/5,"
    )
    print(
        "so inside the wall only C1 and C2 ever cut below 1/4, and C2 cuts exactly from alpha_1 = 3/8."
    )
    bad = []
    exact = []
    for na in range(670, 1000):
        alpha = Fraction(na, 1000)
        top = min(alpha / 2, Fraction(1, 2))
        for i in range(1, 201):
            a1 = top * i / 200
            for name in ("L2", "L3", "L4", "L5"):
                c = cap(name, alpha, a1)
                if c is not None and c < Fraction(1, 4):
                    bad.append((alpha, a1, name, c))
        for name in ("L2", "L3"):
            c = cap(name, alpha, top)
            exact.append(c == Fraction(1, 4))
    print(
        f"sweep of alpha in [67/100, 999/1000] by 1/1000 and alpha_1 in (0, alpha/2] by alpha/400: "
        f"{len(bad)} cells where L2, L3, L4 or L5 falls below 1/4, out of {330 * 200 * 4}"
    )
    print(
        f"L2 and L3 equal 1/4 at the wall in {sum(exact)} of {len(exact)} exact rational tests"
    )
    want(len(bad) == 0, "no inactive cap dips")
    want(all(exact), "L2 and L3 meet the corner")
    print()
    print(
        "So the region is exactly alpha_1 < alpha/2 and beta <= min(1/4, (2/5)(1 - alpha_1)), with alpha_1 >= 1 - alpha and beta >= 1 - alpha at a design:"
    )
    print("alpha        design alpha_1 window        beta window            non-empty")
    for na in (60, 66, 67, 70, 75, 80, 90, 95, 98397 / 1000):
        alpha = Fraction(round(na * 1000), 100000)
        top = min(alpha / 2, Fraction(1, 2))
        fl = 1 - alpha
        c, _ = beta_cap(alpha, fl)
        ok = fl < top and fl <= c
        print(
            f"{fdown(float(alpha), 6):<12} [{fdown(float(fl), 6)}, {fdown(float(top), 6)})       "
            f"[{fdown(float(fl), 6)}, {fdown(float(c), 6)}]      {'yes' if ok else 'no'}"
        )
    print()
    print(
        "The single-window branch, read with the step beta <= m_1 <= alpha_1 that t = 1 in the infimum and one shift of the supremum give:"
    )
    bmax = Fraction(1) / (1 + Fraction(13, 4))
    amin = 1 - bmax
    print(
        f"  alpha_1 <= 1 - (13/4) beta against beta <= alpha_1 asks beta <= {bmax} = {fup(float(bmax))}, under the {Fraction(1, 4)} it replaces"
    )
    print(
        f"  with the Parseval floor beta >= 1 - alpha that asks alpha >= {amin} = {fdown(float(amin), 6)}, over the {Fraction(3, 4)} the window cap asks: the branch narrows the route"
    )
    want(bmax < Fraction(1, 4), "branch cap under a quarter")
    want(amin > Fraction(3, 4), "branch gate over three quarters")
    alpha, a1, beta = Fraction(9, 10), Fraction(77, 500), Fraction(13, 50)
    rows = [
        ("floor alpha_1", a1 >= 1 - alpha),
        ("floor beta", beta >= 1 - alpha),
        ("L1", 2 * a1 < alpha),
        ("branch", a1 <= 1 - Fraction(13, 4) * beta),
    ]
    for name in ("C2", "L2", "L3", "L4", "L5"):
        c = cap(name, alpha, a1)
        rows.append(
            (name, True if c is None else (beta <= c if name == "C2" else beta < c))
        )
    print(
        f"  the step is load-bearing: without it alpha = {fdown(float(alpha), 6)}, alpha_1 = {fdown(float(a1), 6)}, beta = {fdown(float(beta), 6)} > 1/4 passes "
        + ", ".join(n for n, ok in rows if ok)
    )
    want(
        all(ok for _, ok in rows) and beta > Fraction(1, 4), "the load-bearing witness"
    )
    print()
    print("FAILS: " + (", ".join(FAILS) if FAILS else "none"))
    if FAILS:
        raise SystemExit(1)


# VERB CHECK


def verb_check(_argv):
    print(
        "Parseval anchor, its expectation from unique base-q expansion and not from the sweep:"
    )
    for q, F in (
        (3, (0, 1)),
        (5, (0, 1, 2, 3)),
        (10, tuple(d for d in range(10) if d != 5)),
    ):
        k = len(F)
        for L in (2, 3):
            x = q**L
            a = np.arange(x)
            t = a / x
            v = np.ones(x)
            for j in range(L):
                s = np.zeros(x, dtype=complex)
                for f in F:
                    s += np.exp(2j * math.pi * f * (q**j) * t)
                v *= np.abs(s) / k
            got = float((v**2).sum())
            wantv = (q / k) ** L
            print(
                f"  q = {q}, k = {k}, L = {L}: sum of F^2 = {got:.9f} against (q/k)^L = {wantv:.9f}"
            )
            want(abs(got - wantv) < 1e-6 * wantv, f"parseval {q} {L}")
    print()
    print(
        "The two floors, each an upper bound on nothing and a lower bound on both parameters:"
    )
    for q, F in (
        (3, (0, 1)),
        (5, (0, 1, 2, 3)),
        (10, tuple(d for d in range(10) if d != 5)),
        (21, tuple(range(1, 21))),
    ):
        p = params(q, F, tgrid=TGRID[:6] if q <= 10 else TGRID[:1])
        fl = 1.0 - float(p["ahi"])
        print(
            f"  q = {q}, k = {p['k']}: 1 - alpha <= {fup(fl)}, alpha_1 <= {fup(p['a1hi'])}, beta <= {fup(p['bhi'])}"
        )
        want(p["a1hi"] >= fl - 1e-9, f"l1 floor {q}")
        want(p["bhi"] >= fl - 1e-9, f"beta floor {q}")
    print()
    print(
        "Both window bounds against an exact grid sum, two inequalities the construction forces term by term:"
    )
    for q, F, L in (
        (3, (0, 1), 6),
        (5, (0, 1, 2, 3), 5),
        (10, tuple(d for d in range(10) if d != 5), 4),
    ):
        k = len(F)
        x = q**L
        t = np.arange(x) / x
        v = np.ones(x)
        for j in range(L):
            s2 = np.zeros(x, dtype=complex)
            for f in F:
                s2 += np.exp(2j * math.pi * f * (q**j) * t)
            v *= np.abs(s2) / k
        exact = float(v.sum())
        nd, m = DEPTH[q]
        gup, glo = window_factors(q, F, nd, m)
        a = np.arange(x, dtype=np.int64)
        up = np.ones(x)
        dn2 = np.ones(x)
        for j in range(L):
            idx = (a * q**j % x) * q**nd // x
            up *= gup[idx]
            dn2 *= glo[idx]
        path = float(up.sum())
        pathlo = float(dn2.sum())
        print(
            f"  q = {q}, k = {k}, L = {L}: infimum path {pathlo:.6f} at most exact grid {exact:.6f} at most supremum path {path:.6f}, ratios {pathlo / exact:.6f} and {exact / path:.6f}"
        )
        want(exact <= path, f"grid under sup window {q}")
        want(pathlo <= exact, f"inf window under grid {q}")
    print()
    print("The q = 21 recompute from scratch, against the quarter bar:")
    p = params(21, tuple(range(1, 21)), tgrid=TGRID[:1])
    print(
        f"  q = 21 missing 0, window {p['nd']} digits, sub-scan {p['m']}: alpha_1 in [{fdown(p['a1lo'])}, {fup(p['a1hi'])}]"
    )
    print(
        f"  clears 1/4 by {fdown(0.25 - p['a1hi'])}; beta <= alpha_1 at t = 1, so the exceptional-set threshold clears too"
    )
    want(p["a1hi"] < 0.25, "q21 clears")
    print()
    print("FAILS: " + (", ".join(FAILS) if FAILS else "none"))
    if FAILS:
        raise SystemExit(1)


VERBS = {
    "params": verb_params,
    "criterion": verb_criterion,
    "region": verb_region,
    "boundary": verb_boundary,
    "threshold": verb_threshold,
    "check": verb_check,
}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in VERBS:
        print("verbs: " + " ".join(sorted(VERBS)))
        raise SystemExit(2)
    VERBS[sys.argv[1]](sys.argv[2:])
