import collections
import itertools
import sys
import time
from fractions import Fraction

import numpy as np

SPAN = 13
SIZE = 6
ORDERS = (2, 3)
DIAL = ([0, 1, 2], [0, 1, 4, 6], [0, 1, 2, 3, 7], [0, 1, 4, 7, 9], [0, 1, 2, 3, 7, 11])
CAP = 17000
LAGS = 8_000_000


def coarray(S):
    return {x - y for x in S for y in S}


def weights(S):
    w = collections.Counter()
    for x in S:
        for y in S:
            w[x - y] += 1
    return w


def design(G, b, r):
    F = [0]
    for i in range(r):
        F = sorted({f + g * b**i for f in F for g in G})
    return F


def digits(s, b, r):
    return [(s // b**i) % b for i in range(r)]


def brute(S):
    D = coarray(S)
    return {s for s in S if coarray([x for x in S if x != s]) != D}


def fast(S):
    pairs = collections.defaultdict(list)
    for x in S:
        for y in S:
            if x != y and len(pairs[x - y]) < 3:
                pairs[x - y].append((x, y))
    E = set()
    for p in pairs.values():
        if len(p) == 1:
            E.update(p[0])
        elif len(p) == 2:
            E.update(set(p[0]) & set(p[1]))
    return E


def paired(G):
    w = weights(G)
    return {g for g in G if any(w[g - h] == 1 for h in G if h != g)}


def holefree(G):
    a = max(G) - min(G)
    return coarray(G) == set(range(-a, a + 1))


def generators():
    for a in range(1, SPAN + 1):
        for L in range(2, SIZE + 1):
            for mid in itertools.combinations(range(1, a), L - 2):
                G = [0, *mid, a]
                if tuple(a - x for x in reversed(G)) < tuple(G):
                    continue
                if holefree(G):
                    yield G


def exact():
    t0 = time.time()
    gens = list(generators())
    cases = bad = tight_bad = 0
    loose = []
    for G in gens:
        L, a = len(G), max(G)
        M = 2 * a + 1
        U = paired(G)
        E = brute(G)
        assert U <= E
        for r in ORDERS:
            F = design(G, M, r)
            assert len(F) == L**r and coarray(F) == set(
                range(-(M**r - 1) // 2, (M**r + 1) // 2)
            )
            got = brute(F)
            want = {s for s in F if all(d in U for d in digits(s, M, r))}
            cases += 1
            bad += got != want
            tight_bad += (len(got) == len(E) ** r) != (E == U)
            assert fast(F) == got
        if E != U:
            loose.append((G, L, len(E), len(U)))
    econ = [x for x in loose if x[2] == x[1]]
    print(
        f"generators {len(gens)} cases {cases} mismatches {bad} tightness-criterion failures {tight_bad}"
    )
    print(
        f"loose generators (F_G > u/L) {len(loose)}, of them maximally economic (F_G = 1) {len(econ)}"
    )
    byL = collections.Counter(x[1] for x in econ)
    print(
        f"  maximally economic loose by L {dict(sorted(byL.items()))}; card E(G) - u over the loose {sorted({e - u for _, _, e, u in loose})}"
    )
    for G, L, e, u in loose:
        print(f"  {G} L {L} essential {e} paired {u}{' economic' if e == L else ''}")
    for name, G in (
        ("ula3", [0, 1, 2]),
        ("mra4", [0, 1, 4, 6]),
        ("ula4", [0, 1, 2, 3]),
    ):
        M = 2 * max(G) + 1
        row = [(len(brute(design(G, M, r))), len(G) ** r) for r in (1, 2, 3)]
        print(f"  {name} {G} M {M} (essential, sensors) at r = 1, 2, 3: {row}")
    print(f"exact {time.time() - t0:.1f}s")


def law(row, u):
    if all(e == u ** (i + 1) for i, e in enumerate(row) if i):
        return "u^r"
    for k in (2, 3):
        tail = row[k - 1 :]
        if len(tail) < 3:
            continue
        d0, d1 = tail[1] - tail[0], tail[2] - tail[1]
        lam = Fraction(d1, d0) if d0 else Fraction(1)
        c = tail[1] - lam * tail[0]
        if lam >= 1 and all(
            tail[j + 1] == lam * tail[j] + c for j in range(len(tail) - 1)
        ):
            return f"e_(r+1) = {lam} e_r + {c} from r = {k}, fitted on 2 steps, checked on {len(tail) - 3}"
    return "no affine law"


def ess(F, w, A):
    F = np.asarray(F, dtype=np.int64)
    inside = np.zeros(A + 1, dtype=bool)
    inside[F] = True
    out = np.zeros(len(F), dtype=bool)
    ch = max(1, 4_000_000 // len(F))
    for i in range(0, len(F), ch):
        s = F[i : i + ch, None]
        t = s - F[None, :]
        W = w[t + A]
        k = s + t
        two = (W == 2) & (t != 0) & (k >= 0) & (k <= A) & inside[np.clip(k, 0, A)]
        out[i : i + ch] = ((W == 1) & (t != 0)).any(axis=1) | two.any(axis=1)
    return set(F[out].tolist())


def dial():
    t0 = time.time()
    for G in DIAL:
        a, L, u = max(G), len(G), len(paired(G))
        wG = weights(G)
        for b in range(a + 1, 2 * a + 2):
            row, lags, w, A = [], 0, np.ones(1, dtype=np.int64), 0
            for r in range(1, 12):
                An = a * (b**r - 1) // (b - 1)
                if L**r > CAP or 2 * An + 1 > LAGS:
                    break
                wn = np.zeros(2 * An + 1, dtype=np.int64)
                for d, m in wG.items():
                    lo = An + d * b ** (r - 1) - A
                    wn[lo : lo + 2 * A + 1] += m * w
                w, A = wn, An
                F = design(G, b, r)
                assert (
                    len(F) == L**r
                    and F[-1] == A
                    and (w > 0).all()
                    and w.sum() == len(F) ** 2
                )
                E = ess(F, w, A)
                assert len(F) > 250 or (
                    E == brute(F) == fast(F)
                    and weights(F) == {t - A: int(x) for t, x in enumerate(w) if x}
                )
                row.append(len(E))
                lags = 2 * A + 1
            print(
                f"  {G} b {b} e_r {row} {law(row, u)}; {L ** len(row)} sensors, {lags} lags at the last r"
            )
    print(f"dial {time.time() - t0:.1f}s")


if __name__ == "__main__":
    verbs = sys.argv[1:] or ["exact", "dial"]
    for v in verbs:
        {"exact": exact, "dial": dial}[v]()
