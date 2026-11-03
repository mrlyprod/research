import json
import math
import os
import sys

import numpy as np

KBRUTE = 20
KMAX = 64
KFRONT = 20
TOP = 255
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", ".."))
DATA_DIR = os.path.join(REPO, "data", os.path.relpath(HERE, REPO))
CENSUS = os.path.join(REPO, "data", "research", "random", "lab", "rs", "mess")
COLUMNS = (
    "degree",
    "ones",
    "nonlinearity",
    "walsh_max",
    "autocorrelation_max",
    "differential",
    "immunity",
)
WEIGHTS = {
    "carpet": (0, 1, 1),
    "net": (1, 1, 2),
    "tree": (0, 1),
    "void": (0, 2),
    "point": (2,),
    "dust": (0,),
    "line": (1, 2),
    "star": (1, 1),
}


def primes_to(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = bytearray(len(sieve[p * p :: p]))
    return {n for n in range(limit + 1) if sieve[n]}


def fibonacci_to(limit):
    out, a, b = set(), 0, 1
    while a <= limit:
        out.add(a)
        a, b = b, a + b
    return out


def fill(corners, side):
    even, odd = (side + 1) // 2, side // 2
    return sum(even ** (2 - w) * odd**w for w in corners)


def design_set(corners, void, limit):
    out, side = set(), 1
    while side <= limit + 3:
        value = side * side - fill(corners, side) if void else fill(corners, side)
        if value > limit:
            break
        out.add(value)
        side += 2
    return out


def sources(limit):
    out = {
        "evens": set(range(0, limit + 1, 2)),
        "odds": set(range(1, limit + 1, 2)),
        "primes": primes_to(limit),
        "binary": {1 << i for i in range(limit.bit_length()) if (1 << i) <= limit},
        "fibonacci": fibonacci_to(limit),
        "grid_squares": {s * s for s in range(1, limit + 2, 2) if s * s <= limit},
    }
    for name, corners in WEIGHTS.items():
        out[f"{name}_fills"] = design_set(corners, False, limit)
        out[f"{name}_voids"] = design_set(corners, True, limit)
    return out


def indicator(terms, limit):
    return [1 if j in terms else 0 for j in range(limit + 1)]


def lucas_anf(g):
    return [
        sum(g[j] for j in range(w + 1) if math.comb(w, j) % 2) % 2
        for w in range(len(g))
    ]


def degree(c):
    return max((w for w, bit in enumerate(c) if bit), default=-1)


def fwht(values):
    a = values.copy()
    n = len(values)
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        a = np.stack((a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]), axis=1)
        h *= 2
    return a.reshape(-1)


def popcounts(k):
    pc = np.zeros(1, dtype=np.int64)
    for _ in range(k):
        pc = np.concatenate((pc, pc + 1))
    return pc


def brute(g, k, pc):
    table = np.array(g[: k + 1], dtype=np.int64)[pc]
    walsh = fwht(1 - 2 * table)
    top = int(np.abs(walsh).max())
    squares = fwht(walsh * walsh)
    autocorrelation = int(np.abs(squares[1:] // (1 << k)).max())
    anf = table.copy()
    for bit in range(k):
        index = np.arange(1 << k)
        mask = (index >> bit) & 1 == 1
        anf[mask] ^= anf[index[mask] ^ (1 << bit)]
    support = pc[anf == 1]
    nonzero = pc[(walsh != 0) & (pc > 0)]
    return {
        "degree": int(support.max()) if len(support) else -1,
        "ones": int(table.sum()),
        "nonlinearity": (1 << (k - 1)) - top // 2,
        "walsh_max": top,
        "autocorrelation_max": autocorrelation,
        "differential": ((1 << k) + autocorrelation) // 2,
        "immunity": int(nonzero.min()) - 1 if len(nonzero) else k,
    }


def krawtchouk(k):
    rows = []
    for v in range(k + 1):
        poly = [1] + [0] * k
        for _ in range(k - v):
            poly = [poly[i] + (poly[i - 1] if i else 0) for i in range(k + 1)]
        for _ in range(v):
            poly = [poly[i] - (poly[i - 1] if i else 0) for i in range(k + 1)]
        rows.append(poly)
    return rows


def exact(g, k, kraw):
    walsh = [
        sum(kraw[v][w] * (1 - 2 * g[w]) for w in range(k + 1)) for v in range(k + 1)
    ]
    parseval = sum(math.comb(k, v) * walsh[v] ** 2 for v in range(k + 1)) == 1 << (
        2 * k
    )
    corr = [sum(walsh[v] ** 2 * kraw[u][v] for v in range(k + 1)) for u in range(k + 1)]
    assert all(x % (1 << k) == 0 for x in corr)
    corr = [x // (1 << k) for x in corr]
    top = max(abs(x) for x in walsh)
    ac = max(abs(x) for x in corr[1:])
    ones = sum(math.comb(k, w) * g[w] for w in range(k + 1))
    cells = {
        "degree": degree(lucas_anf(g[: k + 1])),
        "ones": ones,
        "nonlinearity": (1 << (k - 1)) - top // 2,
        "walsh_max": top,
        "autocorrelation_max": ac,
        "differential": ((1 << k) + ac) // 2,
        "immunity": next((v for v in range(1, k + 1) if walsh[v] != 0), k + 1) - 1,
    }
    return cells, parseval


def ceiling(k):
    return (1 << (k - 1)) - (1 << (k // 2 - 1) if k % 2 == 0 else 1 << ((k - 1) // 2))


def aggregates(terms_by_name, kraw):
    classes = {}
    for name, terms in terms_by_name.items():
        classes.setdefault(tuple(indicator(terms, KMAX)), []).append(name)
    affine = tight = total = full = cells = 0
    best = (0.0, "", 0)
    for g, names in classes.items():
        g = list(g)
        if degree(lucas_anf(g)) <= 1:
            affine += 1
            continue
        for k in range(3, KMAX + 1):
            row, _ = exact(g, k, kraw[k])
            cells += 1
            full += row["degree"] == k
            if k >= 8:
                total += 1
                tight += row["nonlinearity"] == min(row["ones"], (1 << k) - row["ones"])
            if k >= 6 and row["nonlinearity"] / ceiling(k) > best[0]:
                best = (row["nonlinearity"] / ceiling(k), names[0], k)
    return len(classes), affine, tight, total, full, cells, best


def frontier(k, kraw):
    size = 1 << (k + 1)
    matrix = np.array(kraw, dtype=np.int64).T
    lucas = np.array(
        [[math.comb(w, j) % 2 for w in range(k + 1)] for j in range(k + 1)],
        dtype=np.int64,
    )
    rows = [[0, None, 0, 0, None, 0] for _ in range(k + 1)]
    for start in range(0, size, 1 << 16):
        index = np.arange(start, min(size, start + (1 << 16)), dtype=np.int64)
        g = (index[:, None] >> np.arange(k + 1)) & 1
        walsh = (1 - 2 * g) @ matrix
        top = np.abs(walsh).max(axis=1)
        c = (g @ lucas) % 2
        deg = np.where(c.any(axis=1), k - np.argmax(c[:, ::-1], axis=1), 0)
        even = walsh[:, 0] == 0
        for d in range(k + 1):
            pick = deg == d
            if pick.any():
                rows[d][0] += int(pick.sum())
                best = int(top[pick].min())
                if rows[d][1] is None or best < rows[d][1]:
                    rows[d][1], rows[d][2] = best, 0
                if best == rows[d][1]:
                    rows[d][2] += int((top[pick] == best).sum())
            both = pick & even
            if both.any():
                rows[d][3] += int(both.sum())
                best = int(top[both].min())
                if rows[d][4] is None or best < rows[d][4]:
                    rows[d][4], rows[d][5] = best, 0
                if best == rows[d][4]:
                    rows[d][5] += int((top[both] == best).sum())
    return rows


def main():
    with open(os.path.join(CENSUS, "census.json")) as f:
        census = json.load(f)
    with open(os.path.join(CENSUS, "frontier.json")) as f:
        fronts = json.load(f)
    mine = sources(KMAX)
    bad = []
    for entry in census:
        theirs = indicator(mine[entry["name"]], KMAX)
        if theirs != entry["g"]:
            bad.append(entry["name"])
    print(
        f"the sources: {len(census)} read from the library, {len(census) - len(bad)} rebuilt here from definitions and the fill law, mismatches {bad or 'none'}"
    )
    classes = {}
    for entry in census:
        classes.setdefault(tuple(entry["g"]), []).append(entry["name"])
    print(f"the sources fall into {len(classes)} classes")
    twins = [names for names in classes.values() if len(names) > 1]
    for names in twins:
        group = [e for e in census if e["name"] in names]
        assert all(
            e["rows"] == group[0]["rows"] and e["c"] == group[0]["c"] for e in group
        )
    print(f"the {len(twins)} classes with twins carry identical rows and ANF vectors")
    krawt = {k: krawtchouk(k) for k in range(1, KMAX + 1)}
    pcs = {k: popcounts(k) for k in range(1, KBRUTE + 1)}
    cells, wrong, classes_wrong = 0, [], set()
    for names in classes.values():
        entry = next(e for e in census if e["name"] == names[0])
        for k in range(1, KBRUTE + 1):
            got = brute(entry["g"], k, pcs[k])
            row = dict(zip(("k",) + COLUMNS, entry["rows"][k - 1]))
            for column in COLUMNS:
                cells += 1
                if got[column] != row[column]:
                    wrong.append((names[0], k, column, got[column], row[column]))
                    classes_wrong.add(names[0])
    print(
        f"brute-force truth table, fast Walsh-Hadamard transform and Moebius transform at k <= {KBRUTE}: {cells} cells over {len(classes)} classes, disagreements {len(wrong)}"
    )
    print(f"classes and k where they disagree: {wrong[:12] or 'none'}")
    cells, wrong, off = 0, [], 0
    for names in classes.values():
        entry = next(e for e in census if e["name"] == names[0])
        anf_here = lucas_anf(entry["g"])
        off += anf_here != entry["c"]
        for k in range(1, KMAX + 1):
            got, parseval = exact(entry["g"], k, krawt[k])
            row = dict(zip(("k",) + COLUMNS, entry["rows"][k - 1]))
            off += not parseval
            for column in COLUMNS:
                cells += 1
                if got[column] != row[column]:
                    wrong.append((names[0], k, column, got[column], row[column]))
    print(
        f"exact Krawtchouk, Wiener-Khinchin and binomial-parity ANF at k <= {KMAX}: {cells} cells, disagreements {len(wrong)}, ANF or Parseval failures {off}"
    )
    print(f"first disagreements: {wrong[:12] or 'none'}")
    for label, low in (("keep 0 and 1", 0), ("drop 0", 1), ("drop 0 and 1", 2)):
        kept = {name: {t for t in terms if t >= low} for name, terms in mine.items()}
        sets, affine, tight, total, full, cells, best = aggregates(kept, krawt)
        print(
            f"{label}, rebuilt here: {sets} sets, {affine} affine, NL = min(wt, 2^k - wt) on {tight} of {total} cells with k >= 8, full degree on {full} of {cells} cells, best nl/ceiling {best[0]:.4f} by {best[1]} at k = {best[2]}"
        )
    cells = 0
    for front in fronts:
        k = front["k"]
        if k > KFRONT:
            continue
        mine_rows = frontier(k, krawt[k] if k in krawt else krawtchouk(k))
        for d, (row, got) in enumerate(zip(front["rows"], mine_rows)):
            want = [
                row[0],
                row[1] if row[0] else None,
                row[2],
                row[3],
                row[4] if row[3] else None,
                row[5],
            ]
            cells += len(want)
            assert want == got, f"frontier k={k} degree {d}: library {want} here {got}"
    print(
        f"the exhaustive symmetric frontier at k <= {KFRONT}, all 2^(k+1) functions rebuilt from a matrix product: {cells} cells agree"
    )
    wide = sources(TOP)
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(os.path.join(DATA_DIR, "sequences.txt"), "w") as out:
        for names in classes.values():
            name = names[0]
            g = indicator(wide[name], TOP)
            c = lucas_anf(g)
            deg = [degree(c[: k + 1]) for k in range(1, KMAX + 1)]
            out.write(f"{name}.c ,{','.join(map(str, c))}\n")
            out.write(f"{name}.degree ,{','.join(map(str, deg))}\n")
            if name == "primes":
                ones = [
                    sum(math.comb(k, w) for w in range(k + 1) if g[w])
                    for k in range(1, 40)
                ]
                out.write(f"{name}.ones ,{','.join(map(str, ones))}\n")
    print(
        f"sequences written for the OEIS lookup: {len(classes)} classes, c(w) to w = {TOP} and degree(k) to k = {KMAX}"
    )


if __name__ == "__main__":
    sys.exit(main())
