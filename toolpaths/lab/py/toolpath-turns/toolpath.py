import sys
import time
from itertools import pairwise

U = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
F, R, M, MR = 0, 1, 2, 3
REV = [R, F, MR, M]
MIR = [M, MR, F, R]
NAMES = ["F", "R", "M", "MR"]


def norm(c):
    a, b = c
    return a * a + a * b + b * b


def mul(p, q):
    x, y = p
    s, t = q
    return (x * s - y * t, x * t + y * s + y * t)


def mirror_index(c):
    c2 = mul(c, c)
    n = norm(c)
    for m in range(6):
        if (U[m][0] * n, U[m][1] * n) == c2:
            return m
    return None


def turn(d1, d2):
    t = (d2 - d1) % 6
    return t - 6 if t > 3 else t


def walks(n, c, gentle):
    out = []
    ta, tb = c

    def dfs(path, seen, a, b, d):
        k = len(path)
        if k == n:
            if (a, b) == (ta, tb):
                out.append(tuple(path))
            return
        rem = n - k - 1
        for nd in range(6):
            if k:
                t = abs(turn(d, nd))
                if t == 3 or (gentle and t > 1):
                    continue
            x, y = a + U[nd][0], b + U[nd][1]
            dx, dy = ta - x, tb - y
            if max(abs(dx), abs(dy), abs(dx + dy)) > rem:
                continue
            if (x, y) in seen:
                continue
            seen.add((x, y))
            path.append(nd)
            dfs(path, seen, x, y, nd)
            path.pop()
            seen.discard((x, y))

    dfs([], {(0, 0)}, 0, 0, None)
    return out


def copy(d, f, g, fl, m):
    n = len(g)
    if f == F:
        return [((d + g[j]) % 6, fl[j]) for j in range(n)]
    if f == R:
        return [((d + g[n - 1 - j]) % 6, REV[fl[n - 1 - j]]) for j in range(n)]
    if f == M:
        return [((d + m - g[j]) % 6, MIR[fl[j]]) for j in range(n)]
    return [((d + m - g[n - 1 - j]) % 6, REV[MIR[fl[n - 1 - j]]]) for j in range(n)]


def expand(seq, g, fl, m):
    out = []
    for d, f in seq:
        out.extend(copy(d, f, g, fl, m))
    return out


def points(seq):
    a = b = 0
    P = [(0, 0)]
    for d, _ in seq:
        a, b = a + U[d][0], b + U[d][1]
        P.append((a, b))
    return P


def avoids(seq):
    a = b = 0
    seen = {(0, 0)}
    for d, _ in seq:
        a, b = a + U[d][0], b + U[d][1]
        if (a, b) in seen:
            return False
        seen.add((a, b))
    return True


def counts(seq):
    c = [0, 0, 0, 0]
    for (d1, _), (d2, _) in pairwise(seq):
        c[abs(turn(d1, d2))] += 1
    return c


def ends(g, fl, m):
    out = {}
    for f in (F, R, M, MR):
        if m is None and f >= M:
            continue
        s = copy(0, f, g, fl, m)
        out[f] = (s[0][0], s[0][1], s[-1][0], s[-1][1])
    return out


def orbit(t, f, h, E):
    seen = set()
    x = (t % 6, f, h)
    worst = 0
    while x not in seen:
        seen.add(x)
        worst = max(worst, abs(turn(0, x[0])))
        t, f, h = x
        x = ((t + E[h][0] - E[f][2]) % 6, E[f][3], E[h][1])
    return worst, seen


def closure(g, fl, m):
    E = ends(g, fl, m)
    types = set()
    flags = set(fl)
    todo = [(g[j + 1] - g[j], fl[j], fl[j + 1]) for j in range(len(g) - 1)]
    fdone = set()
    while todo or flags - fdone:
        while todo:
            x = todo.pop()
            x = (x[0] % 6, x[1], x[2])
            if x in types:
                continue
            types.add(x)
            t, f, h = x
            todo.append(((t + E[h][0] - E[f][2]) % 6, E[f][3], E[h][1]))
        for f in list(flags - fdone):
            fdone.add(f)
            s = copy(0, f, g, fl, m)
            flags.update(x[1] for x in s)
            todo.extend((b[0] - a[0], a[1], b[1]) for a, b in pairwise(s))
    return max(abs(turn(0, t)) for t, _, _ in types)


def flag_search(g, flags, m, gentle=True):
    n = len(g)
    shapes = {}
    for d in range(6):
        for f in flags:
            steps = copy(d, f, g, [F] * n, m)
            shapes[(d, f)] = (points(steps)[1:], steps[0][0], steps[-1][0])
    tj = [(g[j + 1] - g[j]) % 6 for j in range(n - 1)]
    res = []
    for f0 in flags:
        for fz in flags:
            proto = [f0] + [F] * (n - 2) + [fz]
            E = ends(g, proto, m)
            good = {
                (t, f, h): not gentle or orbit(t, f, h, E)[0] <= 1
                for t in range(6)
                for f in flags
                for h in flags
            }

            def dfs(j, fl, used, o, last):
                if j == n:
                    res.append(tuple(fl))
                    return
                cand = [f0] if j == 0 else ([fz] if j == n - 1 else flags)
                for f in cand:
                    if j and not good[(tj[j - 1], fl[-1], f)]:
                        continue
                    K, first, end = shapes[(g[j], f)]
                    if last is not None and abs(turn(last, first)) > (
                        1 if gentle else 2
                    ):
                        continue
                    P = [(o[0] + p[0], o[1] + p[1]) for p in K]
                    if any(p in used for p in P):
                        continue
                    used.update(P)
                    fl.append(f)
                    dfs(j + 1, fl, used, P[-1], end)
                    fl.pop()
                    used.difference_update(P)

            dfs(0, [], {(0, 0)}, (0, 0), None)
    return res


def canon(g, fl, mm):
    orb = [(g, fl), (g[::-1], fl[::-1])]
    if mm is not None:
        h = tuple((mm - d) % 6 for d in g)
        orb += [(h, fl), (h[::-1], fl[::-1])]
    return min(orb)


A229214 = (
    1,
    2,
    -1,
    3,
    1,
    1,
    -3,
    1,
    2,
    2,
    -1,
    -2,
    3,
    2,
    3,
    -1,
    -1,
    -3,
    1,
    -2,
    -1,
    3,
    -1,
    -3,
    -2,
    3,
    3,
    2,
    1,
    2,
    -1,
    3,
    1,
    1,
    -3,
    1,
    2,
    -1,
    3,
    1,
    1,
    -3,
    -2,
    -3,
    -3,
    2,
    3,
    1,
    -3,
)

CHORDS = [
    (3, (1, 1)),
    (4, (2, 0)),
    (7, (2, 1)),
    (9, (3, 0)),
    (12, (2, 2)),
    (13, (3, 1)),
    (16, (4, 0)),
    (19, (3, 2)),
]


def census(orders, top, stop, mirrors=True):
    for n, c in CHORDS:
        if n not in orders:
            continue
        t0 = time.time()
        mm = mirror_index(c)
        m = mm if mirrors else None
        flags = [F, R] if m is None else [F, R, M, MR]
        G = walks(n, c, True)
        W = G[:stop] if stop else G
        lv2 = sharp = 0
        cross = {}
        alive = []
        kept = set()
        for g in W:
            for fl in flag_search(g, flags, m):
                lv2 += 1
                if closure(g, fl, m) > 1:
                    sharp += 1
                    continue
                kept.add(canon(g, fl, mm))
                seq = expand([(d, f) for d, f in zip(g, fl)], g, fl, m)
                L = 2
                while L < top:
                    seq = expand(seq, g, fl, m)
                    L += 1
                    if not avoids(seq):
                        cross[L] = cross.get(L, 0) + 1
                        break
                else:
                    alive.append((g, fl))
        mir = "F R" if m is None else "F R M MR"
        sym = "reversal" if mm is None else "reversal and mirror"
        if len(W) == len(G):
            tail = f"classes up to {sym} {len(kept)}, exhausted"
        else:
            tail = f"CUT after walk {len(W)} of {len(G)}"
        print(
            f"order {n} chord {c} flags {mir}: gentle walks {len(G)}, (generator, flags) pairs gentle and point-avoiding at level 2 by the typed search {lv2}, sharp at a later level by the junction closure {sharp}, gentle at every level {lv2 - sharp}, first crossing level {sorted(cross.items())}, alive at level {top} {len(alive)}, {tail}, {time.time() - t0:.1f}s"
        )
        for g, fl in alive[:4]:
            print("  alive", g, [NAMES[f] for f in fl])
        sys.stdout.flush()


def curve(g, fl, m, k):
    seq = [(0, F)]
    for _ in range(k):
        seq = expand(seq, g, fl, m)
    return seq


def harmonic(seq):
    c = [0, 0, 0]
    for d, _ in seq:
        c[d % 3] += 1
    return c, c[0] ** 2 + c[1] ** 2 + c[2] ** 2 - c[0] * c[1] - c[1] * c[2] - c[0] * c[
        2
    ]


def gosper():
    g, fl = (0, 5, 3, 4, 0, 0, 1), (F, R, R, F, F, F, R)
    for k in range(1, 7):
        seq = curve(g, fl, None, k)
        tc = counts(seq)
        want = [7 ** (k - 1), 4 * 7 ** (k - 1) - 1, 2 * 7 ** (k - 1), 0]
        assert tc == want, (k, tc, want)
        assert k > 5 or avoids(seq)
        ax, h2 = harmonic(seq)
        assert h2 == 7**k
        if k == 2:
            assert [(-d) % 6 for d, _ in seq] == [
                x - 1 if x > 0 else 2 - x for x in A229214
            ]
        print(
            f"gosper level {k}: steps {len(seq)}, turns 0/60/120 {tc[:3]}, axis counts {ax}, |h2|^2 {h2}"
        )


def lsys(rules, axiom, k):
    s = axiom
    for _ in range(k):
        s = "".join(rules.get(ch, ch) for ch in s)
    return s


def square(rules, axiom, k):
    d = 0
    h = v = 0
    tc = [0, 0, 0]
    last = None
    for ch in lsys(rules, axiom, k):
        if ch == "+":
            d = (d + 1) % 4
        elif ch == "-":
            d = (d - 1) % 4
        elif ch == "F":
            if last is not None:
                tc[min((d - last) % 4, (last - d) % 4)] += 1
            last = d
            if d % 2:
                v += 1
            else:
                h += 1
    return h, v, tc


def rate(orders):
    for n, c in CHORDS:
        if n not in orders:
            continue
        t0 = time.time()
        m = mirror_index(c)
        flags = [F, R] if m is None else [F, R, M, MR]
        pairs = 0
        rows = {}
        for g in walks(n, c, False):
            for fl in flag_search(g, flags, m, False):
                seq = curve(g, fl, m, 4)
                if not avoids(seq) or counts(seq)[3]:
                    continue
                pairs += 1
                s = [counts(curve(g, fl, m, k))[2] for k in range(1, 6)]
                rows[canon(g, fl, m)] = s
        sym = "reversal" if m is None else "reversal and mirror"
        print(
            f"order {n} chord {c}: (generator, flags) pairs avoiding points through level 4 {pairs}, classes up to {sym} {len(rows)}, {time.time() - t0:.1f}s"
        )
        for (g, fl), s in sorted(rows.items(), key=lambda x: (x[1], x[0])):
            print(
                f"  {g} {[NAMES[f] for f in fl]}: 120-degree turns at levels 1..5 {s}, per unit step at level 5 {s[-1]}/{n**5}"
            )
        sys.stdout.flush()


def bead():
    hil = {"A": "+BF-AFA-FB+", "B": "-AF+BFB+FA-"}
    pea = {"X": "XFYFX+F+YFXFY-F-XFYFX", "Y": "YFXFY-F-XFYFX+F+YFXFY"}
    for k in range(1, 7):
        h, v, tc = square(hil, "A", k)
        assert h + v == 4**k - 1
        print(
            f"hilbert level {k}: steps {h + v}, horizontal {h}, vertical {v}, h2 {h - v}, turns 0/90 {tc[:2]}"
        )
    for k in range(1, 5):
        h, v, tc = square(pea, "X", k)
        assert h + v == 9**k - 1
        print(
            f"peano level {k}: steps {h + v}, horizontal {h}, vertical {v}, h2 {h - v}, turns 0/90 {tc[:2]}"
        )


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "census"
    orders = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else None
    stop = int(sys.argv[3]) if len(sys.argv) > 3 else None
    if verb == "census":
        census(orders or [3, 4, 7, 9, 13], 4, stop)
    if verb == "plain":
        census(orders or [12], 4, stop, False)
    if verb == "gosper":
        gosper()
    if verb == "bead":
        bead()
    if verb == "rate":
        rate(orders or [7])
