import math
import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "py", "mrly-pairing"
    ),
)

import pairing

# MONOID

B = 32

DEEP = 18

WIDE = 14

BASES = [(3, 16), (4, 12), (5, 11), (6, 9), (7, 8)]

RHO3 = 2.207512


def gens(L):
    out = []
    for m in range(2, 1 << L):
        v = 0
        for i in range(L):
            if (m >> i) & 1:
                v |= 1 << (i * B)
        out.append((m.bit_length() - 1, v))
    return out


def deg(p):
    return (p.bit_length() - 1) // B


def coeffs(p):
    c = []
    while p:
        c.append(p & ((1 << B) - 1))
        p >>= B
    return c


def value(p, q):
    v = 0
    for a in reversed(coeffs(p)):
        v = v * q + a
    return v


def show(p):
    t = []
    for i, a in enumerate(coeffs(p)):
        if a == 0:
            continue
        h = "" if a == 1 and i else str(a)
        b = "" if i == 0 else ("x" if i == 1 else "x^" + str(i))
        t.append(h + b)
    return " + ".join(t) if t else "0"


def monoid(L):
    g = gens(L)
    seen = {1}
    cur = [1]
    while cur:
        nxt = []
        for p in cur:
            d = deg(p)
            for dd, v in g:
                if dd + d >= L:
                    break
                r = p * v
                if r not in seen:
                    seen.add(r)
                    nxt.append(r)
        cur = nxt
    return sorted(seen), g


# INVERSE


def nustar(L):
    els, g = monoid(L)
    nu = dict.fromkeys(els, 0)
    nu[1] = 1
    for p in els:
        v = nu[p]
        if v == 0:
            continue
        d = deg(p)
        for dd, w in g:
            if dd + d >= L:
                break
            nu[p * w] -= v
    return els, nu


def power(l):
    return 1 if l == 0 else (-1 if l == 1 else 0)


def ladder(L):
    els, nu = nustar(L)
    for l in range(L):
        assert nu[1 << (l * B)] == power(l)
    run = 0
    best = 0
    mx, sums, cens = [], [], [0] * L
    top = 0
    d = 0
    for p in els:
        while deg(p) > d:
            sums.append(run + power(d + 1))
            mx.append(max(best, abs(run + power(d + 1))))
            d += 1
        run += nu[p]
        cens[d] += 1
        top = max(top, max(coeffs(p)))
        best = max(best, abs(run))
    sums.append(run + power(L))
    mx.append(max(best, abs(run + power(L))))
    return mx, sums, cens, top


# BOUND


def qset(L, els):
    top = [p for p in els if deg(p) < L]
    for q in range(2, 4096):
        hi = max(value(p, q) for p in top)
        if hi < q**L:
            return q, max(top, key=lambda p: value(p, q - 1))
    return None, None


def closed(L):
    return next(q for q in range(2, 4096) if (q + 1) ** (L - 1) < q**L)


def digits(n, q):
    d = []
    while n:
        d.append(n % q)
        n //= q
    return "".join(str(a) for a in reversed(d))


def qord(L, els):
    lst = [p for p in els if deg(p) <= L]
    for q in range(2, 4096):
        prev = -1
        bad = None
        for i, p in enumerate(lst):
            v = value(p, q)
            if v <= prev:
                bad = (lst[i - 1], p)
                break
            prev = v
        if bad is None:
            return q, wit(lst, q - 1)
    return None, None


def wit(lst, q):
    prev = -1
    for i, p in enumerate(lst):
        v = value(p, q)
        if v <= prev:
            return (lst[i - 1], p)
        prev = v
    return None


def window(q):
    return max(l for l in range(2, 64) if all(closed(j) <= q for j in range(2, l + 1)))


def ordwindow(tab, q):
    return max([l for l in sorted(tab) if tab[l] <= q] or [1])


# VERBS


def line(*a):
    print(*a, flush=True)


def lemma():
    els, _ = monoid(WIDE + 1)
    line("THE CARRY BOUND for F = {0,1}, the monoid M* of products of 0/1 polynomials")
    line("")
    line(
        "  max coefficient over degree < L, against the exact cap binomial(L-1, floor((L-1)/2)) and the crude cap 2^(L-1)"
    )
    cmax = [0] * (WIDE + 2)
    for p in els:
        for a in coeffs(p):
            cmax[deg(p)] = max(cmax[deg(p)], a)
    run = 0
    for L in range(2, WIDE + 2):
        run = max(run, cmax[L - 1])
        line(
            "   L %2d  max coef %6d  exact cap %6d  crude cap %8d"
            % (L, run, math.comb(L - 1, (L - 1) // 2), 2 ** (L - 1))
        )
    line("")
    line(
        "  q_set(L): least q with every P of degree < L below q^L, and the last P over the cut"
    )
    line("  q_ord(L): least q with evaluation increasing on all of M* to degree L")
    tset, tord = {}, {}
    for L in range(3, WIDE + 1):
        qs, ps = qset(L, els)
        qo, w = qord(L, els)
        tset[L], tord[L] = qs, qo
        line(
            "   L %2d  q_set %2d  q_ord %2d  closed %2d %s"
            % (L, qs, qo, closed(L), "" if closed(L) == qs else "SPLIT")
        )
        line(
            "          cut at q = %d: %s = %d over q^L = %d, %d digits base %d"
            % (
                qs - 1,
                show(ps),
                value(ps, qs - 1),
                (qs - 1) ** L,
                len(digits(value(ps, qs - 1), qs - 1)),
                qs - 1,
            )
        )
        line(
            "          order at q = %d: %s = %d before %s = %d"
            % (qo - 1, show(w[1]), value(w[1], qo - 1), show(w[0]), value(w[0], qo - 1))
        )
    line("")
    mx, sums, _, _ = ladder(WIDE + 2)
    line(
        "  the pushforward term by term, nu_F(n) against the sum of nu*(P) over P(q) = n"
    )
    for q, L in [(2, 15), (3, 10), (4, 8), (5, 7)]:
        bad, coll, dis = push(q, L)
        line(
            "   q %2d  to n <= q^L = %10d  mismatches %d  colliding values %d  distinct values %d"
            % (q, q**L, bad, coll, dis)
        )
    line("")
    line(
        "  the prediction against nu_F, sibling generator lab/py/mrly-pairing verb inverse"
    )
    line("  base-free ladder: M(q^L) %s" % sums)
    line("  base-free ladder: maxima %s" % mx)
    for q, L in BASES:
        nu = pairing.dirichlet_inverse(q, 0b11, L)
        sig, run = pairing.ladder_stats(nu, q, L)
        ds = next((l + 1 for l in range(len(sig)) if sig[l] != sums[l]), None)
        dm = next((l + 1 for l in range(len(run)) if run[l] != mx[l]), None)
        line(
            "   q %2d  to level %2d  carry-free window %2d  order window %2d  sum departs %s  maxima depart %s"
            % (q, L, window(q), ordwindow(tord, q), ds, dm)
        )
        line("          M(q^L) %s" % sig)
        line("          maxima %s" % run)
        del nu


def push(q, L):
    els, nu = nustar(L + 1)
    g = {}
    for p in els:
        if nu[p] == 0:
            continue
        v = value(p, q)
        if v <= q**L:
            g[v] = g.get(v, 0) + nu[p]
    nf = pairing.dirichlet_inverse(q, 0b11, L)
    bad = 0
    coll = 0
    seen = {}
    for p in els:
        v = value(p, q)
        if v <= q**L:
            seen[v] = seen.get(v, 0) + 1
            coll += seen[v] > 1
    for n in range(1, q**L + 1):
        if int(nf[n]) != g.get(n, 0):
            bad += 1
    return bad, coll, len(seen)


def grep(seq):
    path = os.environ.get("OEIS_STRIPPED")
    if not path or not os.path.exists(path):
        return "no local dump named by OEIS_STRIPPED"
    key = "," + ",".join(str(v) for v in seq) + ","
    hits = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for row in f:
            if key in row:
                hits.append(row.split()[0])
    return "absent" if not hits else ", ".join(hits[:8])


def sequence():
    mx, sums, cens, top = ladder(DEEP)
    line("THE BASE-FREE LADDER of nu*, levels 1..%d" % DEEP)
    line("")
    line("   level  running max  M(q^L)  monoid census by degree")
    for l in range(1, DEEP + 1):
        line("   %2d  %8d  %6d  %8d" % (l, mx[l - 1], sums[l - 1], cens[l - 1]))
    line("")
    line("  running maxima  %s" % mx)
    line("  M(q^L)          %s" % sums)
    line("  monoid census   %s" % cens)
    line("  partial census  %s" % [sum(cens[:l]) for l in range(1, DEEP + 1)])
    line("")
    line(
        "  max coefficient at degree %d is %d, under the packing width 2^%d and the crude cap 2^%d"
        % (DEEP - 1, top, B, DEEP - 1)
    )
    line("")
    line("  OEIS running maxima  %s" % grep(mx[:12]))
    line("  OEIS monoid census   %s" % grep(cens[:12]))
    line("  OEIS partial census  %s" % grep([sum(cens[:l]) for l in range(1, 13)]))


def exponent():
    mx, _, _, _ = ladder(DEEP)
    line("THE BASE-FREE MERTENS EXPONENT against the design mass 2^L")
    line("")
    line("   L   max      2^L      max/2^L   log2(max)/L  step   log2 step")
    for l in range(1, DEEP + 1):
        st = mx[l - 1] / mx[l - 2] if l > 1 and mx[l - 2] else 0.0
        line(
            "   %2d %8d %8d  %9.6f  %9.6f  %6.4f  %9.6f"
            % (
                l,
                mx[l - 1],
                2**l,
                mx[l - 1] / 2**l,
                math.log2(mx[l - 1]) / l if mx[l - 1] else 0.0,
                st,
                math.log2(st) if st else 0.0,
            )
        )
    line("")
    for r in (0, 1):
        cl = [l for l in range(1, DEEP + 1) if l % 2 == r]
        line(
            "  L = %d mod 2   max/2^L  %s"
            % (r, ["%.6f" % (mx[l - 1] / 2**l) for l in cl])
        )
    line("")
    for w in (4, 6, 8):
        g = (mx[DEEP - 1] / mx[DEEP - 1 - w]) ** (1.0 / w)
        line(
            "  geometric mean step over the last %2d levels  %.6f   log2 %.6f"
            % (w, g, math.log2(g))
        )
    line("")
    line(
        "  the design mass rate is 2 and q^(Re rho) at base 3 {0,1} is %.6f, lab/py/mrly-pairing verb box"
        % RHO3
    )


VERBS = {"lemma": lemma, "sequence": sequence, "exponent": exponent}

if __name__ == "__main__":
    VERBS[sys.argv[1]]()
