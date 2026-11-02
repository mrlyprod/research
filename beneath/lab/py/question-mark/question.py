import argparse
import math
import time
import urllib.error
import urllib.request
from fractions import Fraction

import numpy as np
from scipy.special import zeta as hurwitz

OEIS = "https://raw.githubusercontent.com/oeis/oeisdata/main/seq/A002/A002487.seq"
E2 = "0.5312805062772051416"
LN2 = math.log(2.0)
PHI = (1.0 + 5.0**0.5) / 2.0
HOLDER = LN2 / (2.0 * math.log(PHI))
M0 = ((1, 1), (0, 1))
M1 = ((1, 0), (1, 1))

# STERN


def stern(top):
    s = [0, 1]
    for n in range(2, top + 1):
        s.append(s[n // 2] if n % 2 == 0 else s[n // 2] + s[n // 2 + 1])
    return s


def mul(v, m):
    return (v[0] * m[0][0] + v[1] * m[1][0], v[0] * m[0][1] + v[1] * m[1][1])


def mm(a, b):
    return tuple(
        tuple(sum(a[i][t] * b[t][j] for t in range(2)) for j in range(2))
        for i in range(2)
    )


def oeis():
    request = urllib.request.Request(OEIS, headers={"User-Agent": "curl/8"})
    try:
        with urllib.request.urlopen(request, timeout=10) as handle:
            text = handle.read().decode()
    except (urllib.error.URLError, OSError, TimeoutError):
        return None
    terms = []
    for line in text.splitlines():
        if line[:2] in ("%S", "%T", "%U"):
            terms += [
                int(x) for x in line.split(" ", 2)[2].strip().strip(",").split(",")
            ]
    return terms


def brocot_row(depth, unit):
    frontier = [((0, 1), (1, 1) if unit else (1, 0))]
    for _ in range(depth + 1):
        nxt = []
        row = []
        for lo, hi in frontier:
            mid = (lo[0] + hi[0], lo[1] + hi[1])
            row.append(mid)
            nxt += [(lo, mid), (mid, hi)]
        frontier = nxt
    return [Fraction(p, q) for p, q in row]


def carry(span, depth):
    t0 = time.time()
    s = stern((1 << (span + 1)) + 2)
    bad = 0
    for n in range(1 << span):
        v = (s[n], s[n + 1])
        if mul(v, M0) != (s[2 * n], s[2 * n + 1]) or mul(v, M1) != (
            s[2 * n + 1],
            s[2 * n + 2],
        ):
            bad += 1
    print(f"carry: v(2n+b) = v(n) M_b mismatches below 2^{span}: {bad}")
    bad = 0
    for n in range(1, 1 << span):
        v = (0, 1)
        for bit in bin(n)[2:]:
            v = mul(v, M1 if bit == "1" else M0)
        if v != (s[n], s[n + 1]):
            bad += 1
    print(f"carry: v(n) = (0,1) M_(d_1) ... M_(d_L) mismatches below 2^{span}: {bad}")
    ref = oeis()
    if ref is None:
        print("carry: A002487 fetch skipped, offline")
    else:
        print(f"carry: A002487 listed terms {len(ref)}, agree: {ref == s[: len(ref)]}")
    print(
        f"carry: M_0 M_1 = {mm(M0, M1)}, M_1 M_0 = {mm(M1, M0)}, equal: {mm(M0, M1) == mm(M1, M0)}"
    )
    print(
        f"carry: s(3), s(5), s(6) = {s[3]}, {s[5]}, {s[6]} on the words 011, 101, 110"
    )
    for lead, low in (("without", 1), ("with", 0)):
        least = None
        for L in range(2, 12):
            lo = (1 << (L - 1)) if low else 0
            for n in range(lo, 1 << L):
                for m in range(n + 1, 1 << L):
                    if n.bit_count() == m.bit_count() and s[n] != s[m]:
                        least = (
                            L,
                            format(n, f"0{L}b"),
                            format(m, f"0{L}b"),
                            s[n],
                            s[m],
                        )
                        break
                if least:
                    break
            if least:
                break
        print(
            f"carry: least pair of one length and one weight with different s, {lead} leading zeros: {least}"
        )
    print(f"carry: level 2 cells read {s[:4]}")
    worst = 0
    for d in range(depth + 1):
        row = brocot_row(d, False)
        cw = [Fraction(s[n], s[n + 1]) for n in range((1 << d), (1 << (d + 1)))]
        rev = [int(format(i, f"0{d}b")[::-1], 2) if d else 0 for i in range(1 << d)]
        worst = max(worst, sum(1 for i in range(1 << d) if row[i] != cw[rev[i]]))
    print(
        f"carry: Stern-Brocot row against Calkin-Wilf under bit reversal, depths 0..{depth}, mismatches: {worst}"
    )
    print(f"carry: {time.time() - t0:.2f} s")


# QUESTION MARK


def quotients(x):
    out = []
    while x:
        q = 1 / x
        a = q.numerator // q.denominator
        out.append(a)
        x = q - a
    return out


def denjoy(x):
    y = Fraction(0)
    total = 0
    for i, a in enumerate(quotients(x)):
        total += a
        y += (-1) ** i * Fraction(1, 2 ** (total - 1))
    return y


def farey(x):
    return x / (1 - x) if x <= Fraction(1, 2) else (1 - x) / x


def tent(y):
    return 2 * y if y <= Fraction(1, 2) else 2 - 2 * y


def runs_word(qs, length):
    word = "0" * (qs[0] - 1)
    bit = "1"
    for a in qs[1:]:
        word += bit * a
        bit = "1" if bit == "0" else "0"
    return word[:length] if len(word) >= length else None


def question(depth, level, alphabets):
    t0 = time.time()
    q = {Fraction(0): Fraction(0), Fraction(1): Fraction(1)}
    frontier = [(Fraction(0), Fraction(1))]
    for _ in range(depth + 1):
        nxt = []
        for a, b in frontier:
            m = Fraction(a.numerator + b.numerator, a.denominator + b.denominator)
            q[m] = (q[a] + q[b]) / 2
            nxt += [(a, m), (m, b)]
        frontier = nxt
    bad_denjoy = sum(1 for x in q if denjoy(x) != q[x])
    print(
        f"question: mediant recursion against Denjoy on {len(q)} nodes to depth {depth}: mismatches {bad_denjoy}"
    )
    bad_tent = sum(1 for x in q if 0 < x < 1 and denjoy(farey(x)) != tent(q[x]))
    print(f"question: ?(F(x)) = T(?(x)) on the same nodes: mismatches {bad_tent}")
    bad_branch = 0
    for x, qx in q.items():
        for a in range(1, 6):
            if denjoy(1 / (a + x)) != Fraction(2 - qx, 2**a):
                bad_branch += 1
    print(
        f"question: ?(1/(a+x)) = 2^(-a)(2 - ?(x)) for a = 1..5 on the same nodes: mismatches {bad_branch}"
    )
    bad_addr = 0
    for d in range(depth + 1):
        for i, x in enumerate(brocot_row(d, True)):
            if q[x] != Fraction(2 * i + 1, 2 ** (d + 1)):
                bad_addr += 1
    print(
        f"question: ?(node i of row d) = (2i+1)/2^(d+1), depths 0..{depth}: mismatches {bad_addr}"
    )
    for k, alpha in alphabets:
        code = code_of(k, alpha)
        accepted = {w for w in words(level) if accepts(code, k, w)}
        prefixes = prefixes_of(alpha, k, level)
        extra = accepted - prefixes
        missing = prefixes - accepted
        heads = sorted({(w[0], len(w) - len(w.lstrip(w[0]))) for w in extra})
        print(
            f"question: width {k} code {code} alphabet {name(alpha)}: accepted {len(accepted)}, prefixes of ?(E_A) {len(prefixes)}, missing {len(missing)}, extra {len(extra)} with leading runs {heads}"
        )
    print(f"question: {time.time() - t0:.2f} s")


# RULES


def words(length):
    return [format(i, f"0{length}b") for i in range(1 << length)]


def accepts(code, k, w):
    return all((code >> int(w[i : i + k], 2)) & 1 for i in range(len(w) - k + 1))


def forbidden(alpha, k):
    pats = []
    if alpha[0] == "finite":
        m = max(alpha[1])
        pats += ["0" * (m + 1), "1" * (m + 1)]
        gaps = [j for j in range(1, m) if j not in alpha[1]]
    else:
        gaps = sorted(alpha[1])
    for j in gaps:
        pats += ["1" + "0" * j + "1", "0" + "1" * j + "0"]
    assert all(len(p) <= k for p in pats)
    return pats


def code_of(k, alpha):
    pats = forbidden(alpha, k)
    code = 0
    for w in range(1 << k):
        ws = format(w, f"0{k}b")
        if not any(p in ws for p in pats):
            code |= 1 << w
    return code


def name(alpha):
    if alpha[0] == "finite":
        return "{" + ",".join(str(a) for a in alpha[1]) + "}"
    return "N" if not alpha[1] else "N\\{" + ",".join(str(a) for a in alpha[1]) + "}"


def alphabets_at(k):
    out = []
    for mask in range(1, 1 << (k - 1)):
        out.append(("finite", tuple(j + 1 for j in range(k - 1) if (mask >> j) & 1)))
    for mask in range(1 << (k - 2)):
        out.append(("cofinite", tuple(j + 1 for j in range(k - 2) if (mask >> j) & 1)))
    return out


def prefixes_of(alpha, k, level):
    if alpha[0] == "finite":
        letters = list(alpha[1])
    else:
        letters = [a for a in range(1, level + 2) if a not in alpha[1]]
    out = set()
    stack = [((), 0)]
    while stack:
        qs, total = stack.pop()
        if total >= level + 1 and len(qs) >= 2:
            w = runs_word(qs, level)
            if w is not None:
                out.add(w)
            continue
        for a in letters:
            stack.append((qs + (a,), total + a))
    return out


def transfer(code, k):
    n = 1 << (k - 1)
    m = np.zeros((n, n))
    for s in range(n):
        for c in range(2):
            w = s * 2 + c
            if (code >> w) & 1:
                m[s, w % n] = 1.0
    return m


def perron(code, k):
    ev = np.linalg.eigvals(transfer(code, k))
    return float(max(ev.real))


def nacci(alpha):
    if alpha[0] == "finite":
        m = max(alpha[1])
        p = [0] * (m + 1)
        p[m] = 1
        for a in alpha[1]:
            p[m - a] -= 1
        return p
    f = max(alpha[1]) if alpha[1] else 0
    p = [0] * (f + 2)
    p[f + 1] += 1
    p[f] -= 2
    for j in alpha[1]:
        p[f - j + 1] += 1
        p[f - j] -= 1
    return p


def largest_root(p):
    r = np.roots(p[::-1])
    return float(max(z.real for z in r if abs(z.imag) < 1e-9))


def charpoly_int(m):
    n = m.shape[0]
    a = np.array(np.rint(m).astype(int), dtype=object)
    cs = [1]
    mj = a.copy()
    ident = np.array(
        [[1 if i == j else 0 for j in range(n)] for i in range(n)], dtype=object
    )
    for j in range(1, n + 1):
        tr = sum(mj[i, i] for i in range(n))
        c = -tr // j
        cs.append(c)
        mj = a.dot(mj + c * ident)
    return cs


def divides(p, q):
    q = list(q)
    p = list(p)
    while q and q[-1] == 0:
        q.pop()
    while len(q) >= len(p):
        f = q[-1]
        if f % p[-1]:
            return False
        f //= p[-1]
        d = len(q) - len(p)
        for i, c in enumerate(p):
            q[d + i] -= f * c
        while q and q[-1] == 0:
            q.pop()
    return not q


# TRANSFER OPERATOR


def nodes(n):
    x = (1.0 - np.cos(np.pi * np.arange(n) / (n - 1))) / 2.0
    w = (-1.0) ** np.arange(n)
    w[0] /= 2.0
    w[-1] /= 2.0
    return x, w


def bary_rows(t, x, w):
    d = t[:, None] - x[None, :]
    hit = np.abs(d) < 1e-15
    r = w[None, :] / np.where(hit, 1.0, d)
    r = np.where(hit.any(axis=1)[:, None], hit.astype(float), r)
    return r / r.sum(axis=1, keepdims=True)


def diff_matrix(x, w):
    n = len(x)
    d = np.zeros((n, n))
    for j in range(n):
        for i in range(n):
            if i != j:
                d[j, i] = (w[i] / w[j]) / (x[j] - x[i])
        d[j, j] = -d[j].sum()
    return d


def operator(s, alpha, x, w, cut, taylor):
    n = len(x)
    m = np.zeros((n, n))
    if alpha[0] == "finite":
        letters = list(alpha[1])
    else:
        letters = [a for a in range(1, cut + 1) if a not in alpha[1]]
    for a in letters:
        t = 1.0 / (a + x)
        m += ((a + x) ** (-2.0 * s))[:, None] * bary_rows(t, x, w)
    if alpha[0] == "cofinite":
        d = diff_matrix(x, w)
        dk = np.eye(n)
        fact = 1.0
        for k in range(taylor + 1):
            m += hurwitz(2.0 * s + k, cut + 1.0 + x)[:, None] * (dk[0][None, :] / fact)
            dk = dk.dot(d)
            fact *= k + 1
    return m


def leading(s, alpha, x, w, cut, taylor):
    ev = np.linalg.eigvals(operator(s, alpha, x, w, cut, taylor))
    i = int(np.argmax(np.abs(ev)))
    return ev[i].real, abs(ev[i].imag)


def pressure_zero(alpha, modes, cut, taylor):
    x, w = nodes(modes)
    if alpha[0] == "finite" and len(alpha[1]) == 1:
        lam, _ = leading(0.0, alpha, x, w, cut, taylor)
        return 0.0, lam
    lo, hi = (0.0, 1.0) if alpha[0] == "finite" else (0.51, 1.1)
    g = lambda s: math.log(leading(s, alpha, x, w, cut, taylor)[0])
    glo, ghi = g(lo), g(hi)
    assert glo > 0 > ghi, (alpha, glo, ghi)
    for _ in range(70):
        mid = (lo + hi) / 2.0
        gm = g(mid)
        if gm > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-16:
            break
    s = (lo + hi) / 2.0
    return s, leading(s, alpha, x, w, cut, taylor)[0]


def control(modes, cut, taylor):
    a12 = ("finite", (1, 2))
    s, lam = pressure_zero(a12, modes, cut, taylor)
    ref = float(E2)
    digits = -math.log10(abs(s - ref)) if s != ref else 17
    print(
        f"control: A = {{1,2}} pressure zero {s:.16f} against {E2}, gap {abs(s - ref):.1e}, {digits:.1f} digits, leading eigenvalue {lam:.16f}"
    )
    assert abs(s - ref) < 1e-10, "E_2 control failed"
    for m in (24, 32, 48, 56):
        s2, _ = pressure_zero(a12, m, cut, taylor)
        print(
            f"control: A = {{1,2}} at {m} modes: {s2:.16f}, gap to {modes} modes {abs(s2 - s):.1e}"
        )
    sn, lam = pressure_zero(("cofinite", ()), modes, cut, taylor)
    print(
        f"control: A = N pressure zero {sn:.16f} against 1, gap {abs(sn - 1.0):.1e}, leading eigenvalue {lam:.16f}"
    )
    assert abs(sn - 1.0) < 1e-10, "Gauss control failed"
    for c in (500, 1000, 4000):
        s3, _ = pressure_zero(("cofinite", (1,)), modes, c, taylor)
        print(f"control: A = N\\{{1}} at cut {c}: {s3:.16f}")
    s4, _ = pressure_zero(("cofinite", (1,)), modes, cut, taylor)
    print(f"control: A = N\\{{1}} at cut {cut}: {s4:.16f}")
    return s


def table(widths, modes, cut, taylor):
    t0 = time.time()
    control(modes, cut, taylor)
    named = {
        (2, 7): 1.618033988749,
        (3, 23): 1.465571231876,
        (3, 54): 1.324717957244,
        (3, 127): 1.839286755214,
    }
    for (k, code), rho in named.items():
        assert abs(perron(code, k) - rho) < 1e-9, (k, code)
    print(
        "table: transfer matrix convention pinned on the four named codes of the census: ok"
    )
    print("| k | code | A | P_A | rho | log_2 rho | dim_CF | alpha log_2 rho |")
    print("|---|---|---|---|---|---|---|---|")
    rows = []
    for k in widths:
        for alpha in alphabets_at(k):
            code = code_of(k, alpha)
            rho = perron(code, k)
            p = nacci(alpha)
            r = largest_root(p)
            cp = charpoly_int(transfer(code, k))
            div = divides(p, cp[::-1])
            assert abs(rho - r) < 1e-9 and div, (k, code, rho, r, div)
            s, _ = pressure_zero(alpha, modes, cut, taylor)
            dy = math.log(rho) / LN2
            bound = HOLDER * dy
            assert s + 1e-12 >= bound, (k, code, s, bound)
            rows.append((k, code, name(alpha), poly_text(p), rho, dy, s, bound))
            print(
                f"| {k} | {code} | `{name(alpha)}` | `{poly_text(p)}` | {trunc(rho)} | {trunc(dy)} | {trunc(s)} | {trunc(bound)} |"
            )
    seen = {}
    for k, code, a, p, rho, dy, s, b in rows:
        seen.setdefault(a, []).append((k, code, s))
    print(
        f"table: {len(rows)} codes on {len(seen)} alphabets; the same alphabet at two widths prints the same dim_CF: {all(max(v)[2] - min(v)[2] == 0 for v in seen.values())}"
    )
    print(f"table: {time.time() - t0:.2f} s")


def trunc(x):
    return f"{math.floor(x * 1e12 + 1e-3) / 1e12:.12f}"


def poly_text(p):
    terms = []
    for i in range(len(p) - 1, -1, -1):
        c = p[i]
        if c == 0:
            continue
        mono = "" if i == 0 else ("x" if i == 1 else f"x^{i}")
        if abs(c) == 1 and i > 0:
            body = mono
        else:
            body = f"{abs(c)}{mono}" if i > 0 else f"{abs(c)}"
        terms.append(("- " if c < 0 else "+ ") + body)
    text = " ".join(terms)
    return text[2:] if text.startswith("+ ") else "-" + text[2:]


def symmetric_orphans(k):
    codes = {code_of(k, a) for a in alphabets_at(k)}
    total = 1 << (1 << k)
    sym = []
    for code in range(total):
        flip = 0
        rev = 0
        for w in range(1 << k):
            if (code >> w) & 1:
                flip |= 1 << ((1 << k) - 1 - w)
                rev |= 1 << int(format(w, f"0{k}b")[::-1], 2)
        if flip == code and rev == code:
            sym.append(code)
    orphans = [c for c in sym if c not in codes]
    live = [c for c in orphans if perron(c, k) > 1.0 + 1e-9]
    return codes, sym, orphans, live


def obstruction(widths):
    for k in widths:
        codes, sym, orphans, live = symmetric_orphans(k)
        total = 1 << (1 << k)
        print(
            f"obstruction: width {k}: {len(codes)} run-length codes of {total}, {len(sym)} codes fixed by both the digit flip and reversal, {len(orphans)} of those with no alphabet"
        )
        print(
            f"obstruction: width {k}: {len(live)} of the {len(orphans)} carry rho > 1"
        )
        for c in orphans if k <= 3 else live[:1]:
            allowed = [format(w, f"0{k}b") for w in range(1 << k) if (c >> w) & 1]
            forbid = [format(w, f"0{k}b") for w in range(1 << k) if not (c >> w) & 1]
            print(
                f"obstruction: width {k} code {c} allows {allowed if k <= 3 else len(allowed)} forbids {forbid}, rho {trunc(perron(c, k))}"
            )


# GRAPH


def graph_of(code, k):
    n = 1 << (k - 1)
    edges = {}
    tails = {}
    for u in range(n):
        b = 1 - (u & 1)
        s = u
        for a in range(1, k + 1):
            w = 2 * s + b
            if not (code >> w) & 1:
                break
            s = w % n
            if a < k:
                edges.setdefault((u, s), []).append(a)
            else:
                tails[u] = s
    return edges, tails


def components(n, edges, tails):
    succ = {u: set() for u in range(n)}
    for u, v in edges:
        succ[u].add(v)
    for u, v in tails.items():
        succ[u].add(v)
    reach = [[False] * n for _ in range(n)]
    for u in range(n):
        stack = list(succ[u])
        while stack:
            v = stack.pop()
            if not reach[u][v]:
                reach[u][v] = True
                stack += list(succ[v])
    comp = {}
    for u in range(n):
        if reach[u][u]:
            comp[u] = frozenset(v for v in range(n) if reach[u][v] and reach[v][u])
    return comp


def branch_table(x, w, cut):
    return np.array([bary_rows(1.0 / (a + x), x, w) for a in range(1, cut + 1)])


def taylor_rows(x, w, taylor):
    d = diff_matrix(x, w)
    dk = np.eye(len(x))
    fact = 1.0
    rows = []
    for j in range(taylor + 1):
        rows.append(dk[0] / fact)
        dk = dk.dot(d)
        fact *= j + 1
    return np.array(rows)


def graph_operator(s, k, comp, edges, tails, x, B, rows, cut):
    n = len(x)
    states = sorted(comp)
    idx = {u: i for i, u in enumerate(states)}
    m = np.zeros((n * len(states), n * len(states)))
    W = (np.arange(1, cut + 1)[:, None] + x[None, :]) ** (-2.0 * s)
    small = {a: W[a - 1][:, None] * B[a - 1] for a in range(1, k)}
    tail = None
    for (u, v), labels in edges.items():
        if u in comp and v in comp[u]:
            for a in labels:
                m[idx[u] * n : (idx[u] + 1) * n, idx[v] * n : (idx[v] + 1) * n] += (
                    small[a]
                )
    for u, v in tails.items():
        if u in comp and v in comp[u]:
            if tail is None:
                tail = np.einsum("ai,aij->ij", W[k - 1 :], B[k - 1 :])
                for j in range(rows.shape[0]):
                    tail += (
                        hurwitz(2.0 * s + j, cut + 1.0 + x)[:, None] * rows[j][None, :]
                    )
            m[idx[u] * n : (idx[u] + 1) * n, idx[v] * n : (idx[v] + 1) * n] += tail
    return m


def graph_shape(code, k):
    n = 1 << (k - 1)
    edges, tails = graph_of(code, k)
    comp = components(n, edges, tails)
    live = any(u in comp and v in comp[u] for u, v in tails.items())
    out = {u: 0 for u in comp}
    for (u, v), labels in edges.items():
        if u in comp and v in comp[u]:
            out[u] += len(labels)
    for u, v in tails.items():
        if u in comp and v in comp[u]:
            out[u] += 2
    return edges, tails, comp, live, max(out.values(), default=0)


def graph_zero(code, k, x, B, rows, cut):
    edges, tails, comp, live, fan = graph_shape(code, k)
    g = lambda s: math.log(
        max(
            np.linalg.eigvals(
                graph_operator(s, k, comp, edges, tails, x, B, rows, cut)
            ).real
        )
    )
    if fan <= 1:
        return 0.0, math.exp(g(0.0)) if comp else 0.0
    lo, hi = (0.51, 1.1) if live else (0.0, 1.0)
    glo, ghi = g(lo), g(hi)
    assert glo > 0 > ghi, (code, glo, ghi)
    for _ in range(70):
        mid = (lo + hi) / 2.0
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-16:
            break
    s = (lo + hi) / 2.0
    return s, math.exp(g(s))


def markov_prefixes(letters, pairs, level):
    out = set()
    stack = [((), 0)]
    while stack:
        qs, total = stack.pop()
        if total >= level + 1 and len(qs) >= 2:
            w = runs_word(qs, level)
            if w is not None:
                out.add(w)
            continue
        for a in letters:
            if qs and (qs[-1], a) in pairs:
                continue
            stack.append((qs + (a,), total + a))
    return out


def graph(widths, modes, cut, taylor, level):
    t0 = time.time()
    x, w = nodes(modes)
    B = branch_table(x, w, cut)
    rows = taylor_rows(x, w, taylor)
    worst = 0.0
    for k in widths:
        for alpha in alphabets_at(k):
            code = code_of(k, alpha)
            s_graph, _ = graph_zero(code, k, x, B, rows, cut)
            s_alpha, _ = pressure_zero(alpha, modes, cut, taylor)
            worst = max(worst, abs(s_graph - s_alpha))
            assert abs(s_graph - s_alpha) < 1e-12, (k, code, s_graph, s_alpha)
    print(
        f"graph: the graph form recovers the pressure zero of all 18 run-length codes, largest gap {worst:.1e}"
    )
    k = 4
    code = 11892
    accepted = {v for v in words(level) if accepts(code, k, v)}
    prefixes = markov_prefixes((1, 2), {(2, 2)}, level)
    extra = accepted - prefixes
    heads = sorted({(v[0], len(v) - len(v.lstrip(v[0]))) for v in extra})
    print(
        f"graph: code {code} at level {level}: accepted {len(accepted)}, prefixes of ?(M) for M = {{1,2}} without the pair 22: {len(prefixes)}, missing {len(prefixes - accepted)}, extra {len(extra)} with leading runs {heads}"
    )
    edges, _, comp, live, _ = graph_shape(code, k)
    print(
        f"graph: code {code} graph: recurrent states {sorted(format(u, '03b') for u in comp)}, edges {sorted((format(u, '03b'), format(v, '03b'), tuple(l)) for (u, v), l in edges.items() if u in comp and v in comp[u])}, tail edges {live}"
    )
    _, _, _, live_codes = symmetric_orphans(k)
    order = [code] + [c for c in live_codes if c != code]
    print(
        "| code | forbids | states | tail | rho | log_2 rho | dim_CF | alpha log_2 rho |"
    )
    print("|---|---|---|---|---|---|---|---|")
    for c in order:
        rho = perron(c, k)
        dy = math.log(rho) / LN2
        s, _ = graph_zero(c, k, x, B, rows, cut)
        bound = HOLDER * dy
        assert s + 1e-12 >= bound, (c, s, bound)
        _, _, comp, live, _ = graph_shape(c, k)
        forbid = " ".join(
            format(v, f"0{k}b") for v in range(1 << k) if not (c >> v) & 1
        )
        print(
            f"| {c} | `{forbid}` | {len(comp)} | {'yes' if live else 'no'} | {trunc(rho)} | {trunc(dy)} | {trunc(s)} | {trunc(bound)} |"
        )
    print(
        "| code | rule | CF constraint | rho | log_2 rho | dim_CF | alpha log_2 rho |"
    )
    print("|---|---|---|---|---|---|---|")
    named = (
        (2, 7, "no 11", "even quotients 1"),
        (
            3,
            23,
            "at most one 1 per 3",
            "even quotients 1, odd quotients at least 2 past the first",
        ),
        (
            3,
            54,
            "no 11, no 000",
            "even quotients 1, odd quotients 1 or 2 past the first",
        ),
        (3, 127, "no 111", "even quotients 1 or 2"),
    )
    for k2, c, rule, cf in named:
        rho = perron(c, k2)
        dy = math.log(rho) / LN2
        s, _ = graph_zero(c, k2, x, B, rows, cut)
        bound = HOLDER * dy
        assert s + 1e-12 >= bound, (c, s, bound)
        print(
            f"| {c} | {rule} | {cf} | {trunc(rho)} | {trunc(dy)} | {trunc(s)} | {trunc(bound)} |"
        )
    print(f"graph: {time.time() - t0:.2f} s")


# SUBLEADING


def complex_parts(modes, letters):
    x, w = nodes(modes)
    return [np.log(a + x) for a in letters], [
        bary_rows(1.0 / (a + x), x, w) for a in letters
    ]


def determinant(s, logs, B):
    n = len(logs[0])
    m = np.zeros((n, n), dtype=complex)
    for lg, b in zip(logs, B):
        m += np.exp(-2.0 * s * lg)[:, None] * b
    return np.linalg.det(np.eye(n) - m), m


def newton(s, logs, B):
    h = 1e-6
    for _ in range(60):
        f, _ = determinant(s, logs, B)
        fp = (determinant(s + h, logs, B)[0] - determinant(s - h, logs, B)[0]) / (2 * h)
        step = f / fp
        s = s - step
        if abs(step) < 1e-15:
            break
    return s


def winding(logs, B, s0, s1, t0, t1, n=3000):
    pts = [complex(s0 + (s1 - s0) * k / n, t0) for k in range(n)]
    pts += [complex(s1, t0 + (t1 - t0) * k / n) for k in range(n)]
    pts += [complex(s1 - (s1 - s0) * k / n, t1) for k in range(n)]
    pts += [complex(s0, t1 - (t1 - t0) * k / n) for k in range(n)]
    ph = np.angle(np.array([determinant(p, logs, B)[0] for p in pts]))
    d = np.diff(np.concatenate([ph, ph[:1]]))
    d = (d + np.pi) % (2 * np.pi) - np.pi
    return round(d.sum() / (2 * np.pi))


def zeros_in(logs, B, s0, s1, ds, t0, t1, dt):
    sig = np.arange(s0, s1 + ds / 2, ds)
    tau = np.arange(t0, t1 + dt / 2, dt)
    F = np.array(
        [[abs(determinant(complex(a, b), logs, B)[0]) for b in tau] for a in sig]
    )
    mins = sorted(
        (F[i, j], sig[i], tau[j])
        for i in range(1, len(sig) - 1)
        for j in range(1, len(tau) - 1)
        if F[i, j] == F[i - 1 : i + 2, j - 1 : j + 2].min()
    )
    zeros = []
    for _, a, b in mins:
        z = newton(complex(a, b), logs, B)
        inside = s0 <= z.real <= s1 and t0 <= z.imag <= t1
        if (
            inside
            and abs(determinant(z, logs, B)[0]) < 1e-10
            and not any(abs(z - q) < 1e-8 for q in zeros)
        ):
            zeros.append(z)
    return sorted(zeros, key=lambda z: -z.real)


def census_of(m, jmax):
    import importlib.util
    import os

    path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "..",
        "..",
        "apollonian",
        "lab",
        "py",
        "ford-horocycle",
        "ford_horocycle.py",
    )
    spec = importlib.util.spec_from_file_location("ford_horocycle", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.census_walk(m, 2**jmax)


def fit_zero(cum, delta, z, jlo, jmax, per):
    js = np.arange(jlo * per, jmax * per + 1) / per
    Q = np.array([round(2.0**j) for j in js])
    R = np.array([int(cum[q]) for q in Q], dtype=float) / Q ** (2 * delta)
    lnQ = np.log(Q)
    env = Q ** (2 * (z.real - delta))
    X0 = np.ones((len(Q), 1))
    X1 = np.column_stack(
        [
            np.ones(len(Q)),
            env * np.cos(2 * z.imag * lnQ),
            env * np.sin(2 * z.imag * lnQ),
        ]
    )
    out = []
    for X in (X0, X1):
        coef, *_ = np.linalg.lstsq(X, R, rcond=None)
        res = R - X @ coef
        out.append((coef, math.sqrt(float(np.mean(res**2))), float(np.abs(res).max())))
    return Q, R, out


def subleading(modes, jmax, jlo, per):
    t0 = time.time()
    delta = float(E2)
    letters = (1, 2)
    logs, B = complex_parts(modes, letters)
    small = zeros_in(logs, B, 0.0, 0.53, 0.01, 0.2, 14.0, 0.05)
    assert len(small) == 1 and winding(logs, B, 0.0, 0.53, 0.2, 14.0) == 1, small
    z = small[0]
    _, m = determinant(z, logs, B)
    ev = np.linalg.eigvals(m)
    near = ev[int(np.argmin(np.abs(ev - 1.0)))]
    print(
        f"subleading: A = {{1,2}}, one zero of det(1 - L_s) in 0 <= sigma <= 0.53, 0.2 <= tau <= 14 at {modes} modes, the winding number agreeing: s_1 = {z.real:.12f} + {z.imag:.12f} i, eigenvalue nearest 1 there {near.real:.12f} {near.imag:+.1e} i"
    )
    for m2 in (60, 100, 140):
        logs2, B2 = complex_parts(m2, letters)
        z2 = newton(z, logs2, B2)
        print(
            f"subleading: s_1 at {m2} modes {z2.real:.15f} + {z2.imag:.15f} i, gap to {modes} modes {abs(z2 - z):.1e}"
        )
    logs2, B2 = complex_parts(100, letters)
    wide = zeros_in(logs2, B2, 0.0, 0.53, 0.01, 0.2, 80.0, 0.1)
    count = winding(logs2, B2, 0.0, 0.53, 0.2, 80.0)
    assert len(wide) == count and abs(wide[0] - z) < 1e-9, (len(wide), count)
    print(
        f"subleading: {len(wide)} zeros in 0 <= sigma <= 0.53, 0.2 <= tau <= 80 at 100 modes, the winding number agreeing; the three of largest real part: "
        + ", ".join(f"{q.real:.9f} + {q.imag:.9f} i" for q in wide[:3])
    )
    period = math.pi / (z.imag * LN2)
    ratio = 2.0 ** (2.0 * (z.real - delta))
    alias = 2.0 * z.imag * LN2 - 3.0 * math.pi
    print(
        f"subleading: delta_2 - sigma_1 = {delta - z.real:.6f}; a zero s = sigma + i tau puts Q^(2 sigma) cos(2 tau log Q) into N(Q), so the period in octaves is pi/(tau log 2) = {period:.4f} and the amplitude ratio per octave 2^(2(sigma - delta_2)) = {ratio:.4f}; per octave the phase advances 2 tau log 2 = 3 pi + {alias:.4f}, so at integer octaves the wave reads as a sign alternation under an envelope of period {2 * math.pi / alias:.1f} octaves"
    )
    t1 = time.time()
    cum = census_of(2, jmax)
    print(
        f"subleading: E_2 census walk to 2^{jmax} by lab/py/ford-horocycle census_walk, N_2(2^{jmax}) = {int(cum[2**jmax])} ({time.time() - t1:.1f}s)"
    )
    Q, R, (flat, one) = fit_zero(cum, delta, z, jlo, jmax, per)
    C, b1, b2 = one[0]
    print(
        f"subleading: N_2(Q)/Q^(2 delta_2) on {len(Q)} points, {per} per octave from 2^{jlo} to 2^{jmax}: constant fit C = {flat[0][0]:.6f} leaves rms {flat[1]:.2e}, max {flat[2]:.2e}; constant plus the wave of s_1 with only C, amplitude and phase fitted: C = {C:.6f}, amplitude {math.hypot(b1, b2):.6f}, phase {math.atan2(-b2, b1):.4f}, rms {one[1]:.2e}, max {one[2]:.2e}"
    )
    fit = lambda q: (
        C
        + q ** (2 * (z.real - delta))
        * (
            b1 * math.cos(2 * z.imag * math.log(q))
            + b2 * math.sin(2 * z.imag * math.log(q))
        )
    )
    rows = []
    for j in range(jlo, jmax):
        q = 2**j
        obs = math.log2(int(cum[2 * q]) / int(cum[q])) - 2 * delta
        rows.append((j, obs, math.log2(fit(2 * q) / fit(q))))
    print(
        "subleading: octave exponent minus 2 delta_2, observed against the fitted wave of s_1: "
        + ", ".join(f"j = {j}: {o:+.4f} / {p:+.4f}" for j, o, p in rows)
    )
    assert all(abs(o - p) < 4e-3 for _, o, p in rows), rows
    from scipy.optimize import least_squares

    lnQ = np.log(Q)
    resid = lambda p: (
        p[2]
        + Q ** (2 * (p[0] - delta))
        * (p[3] * np.cos(2 * p[1] * lnQ) + p[4] * np.sin(2 * p[1] * lnQ))
        - R
    )
    sol = least_squares(resid, [z.real, z.imag, C, b1, b2])
    print(
        f"subleading: sigma and tau fitted freely on the same points: {sol.x[0]:.4f}, {sol.x[1]:.4f} against the operator's {z.real:.4f}, {z.imag:.4f}, rms {math.sqrt(float(np.mean(sol.fun**2))):.2e}"
    )
    delta3 = 0.705660908028
    logs3, B3 = complex_parts(100, (1, 2, 3))
    three = zeros_in(logs3, B3, 0.0, 0.71, 0.01, 0.2, 80.0, 0.05)
    count3 = winding(logs3, B3, 0.0, 0.71, 0.2, 80.0)
    assert len(three) == count3, (len(three), count3)
    z3 = three[0]
    cum3 = census_of(3, 18)
    _, _, (flat3, one3) = fit_zero(cum3, delta3, z3, 10, 18, per)
    print(
        f"subleading: control A = {{1,2,3}}: {len(three)} zeros in 0 <= sigma <= 0.71, 0.2 <= tau <= 80 at 100 modes, the winding number agreeing, the largest real part at {z3.real:.12f} + {z3.imag:.12f} i, delta_3 - sigma = {delta3 - z3.real:.4f}, ratio per octave {2.0 ** (2.0 * (z3.real - delta3)):.4f}; N_3(Q)/Q^(2 delta_3) from 2^10 to 2^18: constant fit rms {flat3[1]:.2e}, with the wave {one3[1]:.2e}"
    )
    print(f"subleading: {time.time() - t0:.2f} s")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "verb",
        choices=[
            "carry",
            "question",
            "table",
            "obstruction",
            "graph",
            "subleading",
            "all",
        ],
    )
    parser.add_argument("--span", type=int, default=16)
    parser.add_argument("--depth", type=int, default=12)
    parser.add_argument("--level", type=int, default=14)
    parser.add_argument("--modes", type=int, default=40)
    parser.add_argument("--cut", type=int, default=2000)
    parser.add_argument("--taylor", type=int, default=4)
    parser.add_argument("--jmax", type=int, default=24)
    parser.add_argument("--jlo", type=int, default=12)
    parser.add_argument("--per", type=int, default=16)
    args = parser.parse_args()
    widths = (2, 3, 4)
    if args.verb in ("carry", "all"):
        carry(args.span, args.depth)
    if args.verb in ("question", "all"):
        alphas = [(k, a) for k in widths for a in alphabets_at(k)]
        question(args.depth, args.level, alphas)
    if args.verb in ("table", "all"):
        table(widths, args.modes, args.cut, args.taylor)
    if args.verb in ("obstruction", "all"):
        obstruction(widths)
    if args.verb in ("graph", "all"):
        graph(widths, args.modes, args.cut, args.taylor, args.level)
    if args.verb in ("subleading", "all"):
        subleading(args.modes, args.jmax, args.jlo, args.per)


if __name__ == "__main__":
    main()
