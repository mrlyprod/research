import csv
import math
import subprocess
import sys
import tempfile
import time
from itertools import permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent

CASES = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2)]
A000616 = {1: 3, 2: 6, 3: 22, 4: 402}
NAMED = {
    (1, -1, -1): "golden",
    (1, -1, 0, -1): "supergolden",
    (1, 0, -1, -1): "plastic",
    (1, -1, -1, -1): "tribonacci",
}

# WINDOWS


def digits_of(w, d, k):
    a = 1 << d
    return [(w >> (d * (k - 1 - j))) & (a - 1) for j in range(k)]


def word_of(ds, d):
    w = 0
    for c in ds:
        w = (w << d) | c
    return w


# GROUP


def bd_apply(c, perm, flip, d):
    out = 0
    for i in range(d):
        b = ((c >> i) & 1) ^ ((flip >> i) & 1)
        out |= b << perm[i]
    return out


def window_tables(d, k, reversal):
    nw = 1 << (d * k)
    revs = (False, True) if reversal else (False,)
    tables = set()
    for perm in permutations(range(d)):
        for flip in range(1 << d):
            for rev in revs:
                t = []
                for w in range(nw):
                    ds = [bd_apply(c, perm, flip, d) for c in digits_of(w, d, k)]
                    if rev:
                        ds = ds[::-1]
                    t.append(word_of(ds, d))
                tables.add(tuple(t))
    return sorted(tables)


def code_array(table, nbits):
    n = 1 << nbits
    arr = [0] * n
    for code in range(1, n):
        low = code & -code
        arr[code] = arr[code ^ low] | (1 << table[low.bit_length() - 1])
    return arr


def orbit_walk(maps, ncodes):
    seen = bytearray(ncodes)
    reps, sizes = [], []
    for c in range(ncodes):
        if seen[c]:
            continue
        orb = {m[c] for m in maps}
        for x in orb:
            seen[x] = 1
        reps.append(c)
        sizes.append(len(orb))
    return reps, sizes


# TRANSFER MATRIX


def transfer(code, d, k):
    a = 1 << d
    if k == 1:
        return [[code.bit_count()]]
    s_n = a ** (k - 1)
    m = [[0] * s_n for _ in range(s_n)]
    for s in range(s_n):
        for c in range(a):
            w = s * a + c
            if (code >> w) & 1:
                m[s][w % s_n] = 1
    return m


def matmul(x, y):
    n = len(x)
    return [
        [sum(x[i][t] * y[t][j] for t in range(n)) for j in range(n)] for i in range(n)
    ]


def charpoly(m):
    n = len(m)
    cs = [1]
    mj = [row[:] for row in m]
    for j in range(1, n + 1):
        tr = sum(mj[i][i] for i in range(n))
        cj = -tr // j
        cs.append(cj)
        if j < n:
            t = [row[:] for row in mj]
            for i in range(n):
                t[i][i] += cj
            mj = matmul(m, t)
    return tuple(cs)


def poly_str(cs):
    n = len(cs) - 1
    parts = []
    for i, c in enumerate(cs):
        if c == 0:
            continue
        e = n - i
        mono = "1" if e == 0 else ("x" if e == 1 else f"x^{e}")
        if abs(c) != 1 or e == 0:
            mono = f"{abs(c)}*{mono}" if e else f"{abs(c)}"
        parts.append(("- " if c < 0 else "+ ") + mono)
    s = " ".join(parts)
    return s[2:] if s.startswith("+ ") else "-" + s[2:]


def divides(p, q):
    p = list(p)
    q = list(q)
    while len(q) >= len(p):
        if q[0] % p[0]:
            return False
        f = q[0] // p[0]
        for i, c in enumerate(p):
            q[i] -= f * c
        if q[0] != 0:
            return False
        q.pop(0)
    return all(c == 0 for c in q)


# PARI


def pari_roots(polys, scratch):
    body = (
        "for(i = 1, #V, p = V[i]; f = factor(p); best = -1; bg = 0; "
        "for(j = 1, #f~, g = f[j, 1]; if(poldegree(g) > 0, r = polrootsreal(g); "
        "if(#r > 0, m = vecmax(r); if(m > best, best = m; bg = g)))); "
        "s = 0; r = polroots(bg); for(t = 1, #r, if(abs(r[t] - best) > 1e-25, s = max(s, abs(r[t])))); "
        'print(i, "|", Vec(bg), "|", best, "|", s))'
    )
    lines = [f"V = [{', '.join(poly_str(p) for p in polys)}];", body, "quit()"]
    src = scratch / "perron.gp"
    src.write_text("\n".join(lines) + "\n")
    out = subprocess.run(
        ["gp", "-q", str(src)],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
        check=False,
    )
    if out.returncode != 0:
        sys.exit(out.stderr)
    got = {}
    for line in out.stdout.strip().splitlines():
        idx, vec, val, sec = line.split("|")
        cs = tuple(int(t) for t in vec.strip().strip("[]").split(","))
        got[int(idx) - 1] = (cs, val.strip(), float(sec))
    return [got[i] for i in range(len(polys))]


# CENSUS


def burnside(d, k, reversal):
    total = 0
    tables = window_tables(d, k, reversal)
    for t in tables:
        seen = [False] * len(t)
        cycles = 0
        for w in range(len(t)):
            if seen[w]:
                continue
            cycles += 1
            x = w
            while not seen[x]:
                seen[x] = True
                x = t[x]
        total += 1 << cycles
    return total // len(tables)


def product_codes(d, k):
    a = 1 << d
    out = set()
    for f in range(1 << a):
        full = 0
        for w in range(1 << (d * k)):
            if all((f >> c) & 1 for c in digits_of(w, d, k)):
                full |= 1 << w
        out.add(full)
    return out


def main():
    scratch = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.gettempdir())
    t_all = time.time()
    census, table = [], []
    for d, k in CASES:
        t0 = time.time()
        nbits = 1 << (d * k)
        ncodes = 1 << nbits
        maps_g = [code_array(t, nbits) for t in window_tables(d, k, True)]
        maps_b = [code_array(t, nbits) for t in window_tables(d, k, False)]
        reps, sizes = orbit_walk(maps_g, ncodes)
        reps_b, _ = orbit_walk(maps_b, ncodes)
        orbit_of = dict(zip(reps, sizes))
        products = product_codes(d, k)
        t_orbit = time.time() - t0

        t0 = time.time()
        polys = {}
        for c in reps:
            polys.setdefault(charpoly(transfer(c, d, k)), []).append(c)
        t_poly = time.time() - t0

        t0 = time.time()
        keys = sorted(polys)
        roots = pari_roots(keys, scratch)
        t_pari = time.time() - t0

        rows = []
        strict = {}
        for key, (minp, val, sec) in zip(keys, roots):
            strict[minp] = sec
            for c in polys[key]:
                w = c.bit_count()
                rho = float(val)
                zero = w > 0 and divides(minp, [1] + [0] * (k - 1) + [-w])
                rows.append((c, w, minp, val, rho, zero))
        n_zero = sum(1 for r in rows if r[5])
        bad = [r for r in rows if r[5] and r[0] not in products]
        live = [r for r in rows if r[4] > 0]
        kappas = [math.log2(r[1]) / k - math.log2(r[4]) for r in live]
        k_max = max(kappas)
        tied = [r for kap, r in zip(kappas, live) if kap > k_max - 1e-12]
        rho1 = [r for r in live if abs(r[4] - 1.0) < 1e-12]
        by_poly = {}
        for c, w, minp, val, rho, zero in rows:
            e = by_poly.setdefault(
                minp,
                {"rho": val, "classes": 0, "rules": 0, "code": c, "w": w, "ks": []},
            )
            e["classes"] += 1
            e["rules"] += orbit_of[c]
            if rho > 0 and w > 0:
                e["ks"].append(math.log2(w) / k - math.log2(rho))
            if c < e["code"]:
                e["code"], e["w"] = c, w
        weak = {
            m: e
            for m, e in by_poly.items()
            if float(e["rho"]) > 0 and strict[m] >= float(e["rho"]) - 1e-20
        }
        n_weak = len(weak)
        roots_seen = {
            m
            for m in by_poly
            if float(by_poly[m]["rho"]) > 0
            and strict[m] < float(by_poly[m]["rho"]) - 1e-20
        }
        n_radical = 0
        for m in weak:
            step = 0
            for i, c in enumerate(m):
                if c:
                    step = math.gcd(step, len(m) - 1 - i)
            if step > 1 and tuple(m[i] for i in range(0, len(m), step)) in roots_seen:
                n_radical += 1
        entry = {
            "D": d,
            "k": k,
            "rules": ncodes,
            "classes_G": len(reps),
            "classes_B": len(reps_b),
            "group_G": len(maps_g),
            "group_B": len(maps_b),
            "charpolys": len(keys),
            "perron_polys": len(by_poly),
            "classes_kappa0": n_zero,
            "kappa0_nonproduct": len(bad),
            "kappa_min": round(min(kappas), 6),
            "kappa_max": round(k_max, 6),
            "kappa_max_code": min(r[0] for r in tied),
            "kappa_max_ties": len(tied),
            "kappa_max_windows": max(r[1] for r in tied),
            "rho1_windows": max((r[1] for r in rho1), default=0),
            "weak_perron_polys": n_weak,
            "weak_are_radicals": n_radical,
            "burnside_G": burnside(d, k, True),
            "burnside_B": burnside(d, k, False),
            "dead_classes": len(rows) - len(live),
            "a000616": A000616.get(d * k, ""),
            "t_orbit": round(t_orbit, 2),
            "t_poly": round(t_poly, 2),
            "t_pari": round(t_pari, 2),
        }
        census.append(entry)
        for minp, e in sorted(by_poly.items(), key=lambda kv: -float(kv[1]["rho"])):
            table.append(
                {
                    "D": d,
                    "k": k,
                    "minpoly": poly_str(minp),
                    "degree": len(minp) - 1,
                    "rho": e["rho"][:14],
                    "classes": e["classes"],
                    "rules": e["rules"],
                    "code": e["code"],
                    "windows": e["w"],
                    "kappa_min": round(min(e["ks"]), 6) if e["ks"] else "",
                    "kappa_max": round(max(e["ks"]), 6) if e["ks"] else "",
                    "strict": int(strict[minp] < float(e["rho"]) - 1e-20)
                    if float(e["rho"]) > 0
                    else "",
                    "name": NAMED.get(minp, ""),
                }
            )
        print(
            f"D={d} k={k} rules={ncodes} classes(G)={len(reps)} classes(B)={len(reps_b)} "
            f"charpoly={len(keys)} perron={len(by_poly)} weak={n_weak} kappa0={n_zero} nonproduct={len(bad)} "
            f"kappa_max={round(k_max, 6)} ties={len(tied)} "
            f"orbit={t_orbit:.2f}s poly={t_poly:.2f}s pari={t_pari:.2f}s"
        )
        for r in bad:
            print(
                f"  NONPRODUCT kappa=0 code={r[0]} windows={r[1]} minpoly={poly_str(r[2])}"
            )

    with (HERE / "census.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(census[0]))
        w.writeheader()
        w.writerows(census)
    with (HERE / "classes.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    print(f"live classes {sum(r['classes_G'] - r['dead_classes'] for r in census)}")
    seen_polys = {}
    for r in table:
        seen_polys.setdefault((r["D"], r["k"]), set()).add(r["minpoly"])
    for a, b in (
        ((1, 3), (1, 4)),
        ((1, 2), (1, 3)),
        ((1, 1), (1, 2)),
        ((2, 1), (2, 2)),
    ):
        miss = seen_polys[a] - seen_polys[b]
        print(
            f"nest {a} -> {b}: {len(seen_polys[a])} of {len(seen_polys[b])}, missing {len(miss)}"
        )
    t0 = time.time()
    print(
        "burnside D=1 G  " + ", ".join(str(burnside(1, k, True)) for k in range(1, 9))
    )
    print(
        "burnside D=1 B  " + ", ".join(str(burnside(1, k, False)) for k in range(1, 9))
    )
    print(
        "burnside D=2 G  " + ", ".join(str(burnside(2, k, True)) for k in range(1, 5))
    )
    print(
        "burnside D=2 B  " + ", ".join(str(burnside(2, k, False)) for k in range(1, 5))
    )
    print(f"burnside extension {time.time() - t0:.2f}s")
    print(f"total {time.time() - t_all:.2f}s")


# CLOSED FORM


def classes_closed(k):
    if k % 2 == 0:
        m = k // 2
        return (
            (1 << (2 ** (2 * m) - 2))
            + (1 << (2 ** (2 * m - 1) - 2))
            + (1 << (2 ** (2 * m - 1) + 2 ** (m - 1) - 1))
        )
    m = (k - 1) // 2
    return (
        (1 << (2 ** (2 * m + 1) - 2))
        + (1 << (2 ** (2 * m) - 1))
        + (1 << (2 ** (2 * m) + 2**m - 2))
    )


def word_orbits(d, k):
    tables = window_tables(d, k, True)
    seen = [False] * (1 << (d * k))
    n = 0
    for w in range(len(seen)):
        if seen[w]:
            continue
        n += 1
        for t in tables:
            seen[t[w]] = True
    return n


def verb_burnside(top=11):
    t0 = time.time()
    cycle = [burnside(1, k, True) for k in range(1, top + 1)]
    closed = [classes_closed(k) for k in range(1, top + 1)]
    words = [word_orbits(1, k) for k in range(1, top + 1)]
    print("classes G_(1,k) k=1..8  " + ", ".join(str(v) for v in cycle[:8]))
    for k in range(9, top + 1):
        print(f"classes G_(1,{k})  {cycle[k - 1]}")
    agree = sum(1 for a, b in zip(cycle, closed) if a == b)
    print(f"closed form against the cycle index k=1..{top}: {agree} of {top} agree")
    for k in range(1, top + 1):
        if cycle[k - 1] != closed[k - 1]:
            print(f"  MISMATCH k={k}")
    print(
        "digits k=9,10,11  "
        + ", ".join(str(len(str(cycle[k - 1]))) for k in (9, 10, 11))
    )
    print("words under the same group k=1..8  " + ", ".join(str(v) for v in words[:8]))
    print(f"burnside {time.time() - t0:.2f}s")


# DETERMINANT


def det_bareiss(m):
    n = len(m)
    a = [row[:] for row in m]
    sign = 1
    prev = 1
    for i in range(n - 1):
        if a[i][i] == 0:
            p = -1
            for r in range(i + 1, n):
                if a[r][i]:
                    p = r
                    break
            if p < 0:
                return 0
            a[i], a[p] = a[p], a[i]
            sign = -sign
        for r in range(i + 1, n):
            for c in range(i + 1, n):
                a[r][c] = (a[r][c] * a[i][i] - a[r][i] * a[i][c]) // prev
            a[r][i] = 0
        prev = a[i][i]
    return sign * a[n - 1][n - 1]


def perm_sign(p):
    seen = [False] * len(p)
    s = 1
    for i in range(len(p)):
        if seen[i]:
            continue
        j, c = i, 0
        while not seen[j]:
            seen[j] = True
            j = p[j]
            c += 1
        if c % 2 == 0:
            s = -s
    return s


def block_frame(k):
    s_n = 1 << (k - 1)
    h = 1 << (k - 2)
    rows, cols = [], []
    for s in range(h):
        rows += [s, s + h]
        cols += [(2 * s) % s_n, (2 * s + 1) % s_n]
    return rows, cols, perm_sign(rows) * perm_sign(cols)


def block_det(m, rows, cols, sign):
    d = 1
    for i in range(0, len(rows), 2):
        d *= (
            m[rows[i]][cols[i]] * m[rows[i + 1]][cols[i + 1]]
            - m[rows[i]][cols[i + 1]] * m[rows[i + 1]][cols[i]]
        )
        if d == 0:
            return 0
    return sign * d


def block_support(k):
    s_n = 1 << (k - 1)
    h = 1 << (k - 2)
    m = transfer((1 << (1 << k)) - 1, 1, k)
    cover = []
    for s in range(h):
        pair = {(2 * s) % s_n, (2 * s + 1) % s_n}
        for r in (s, s + h):
            if {c for c in range(s_n) if m[r][c]} - pair:
                return False
        cover += sorted(pair)
    return sorted(cover) == list(range(s_n))


# PERRON


def pari_minpolys(polys, scratch):
    body = (
        "for(i = 1, #V, p = V[i]; rho = vecmax(abs(polroots(p))); "
        "f = factor(p); g = 0; n = 0; w = 0; "
        "for(j = 1, #f~, h = f[j, 1]; if(poldegree(h) > 0, s = polrootsreal(h); "
        "if(#s > 0 && abs(vecmax(s) - rho) < 1e-20, g = h; n = n + 1; "
        "w = vecmax(abs(polroots(h)))))); "
        'print(i, "|", Vec(g), "|", n, "|", rho, "|", w))'
    )
    lines = [f"V = [{', '.join(poly_str(p) for p in polys)}];", body, "quit()"]
    src = scratch / "lemmas.gp"
    src.write_text("\n".join(lines) + "\n")
    out = subprocess.run(
        ["gp", "-q", str(src)],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
        check=False,
    )
    if out.returncode != 0:
        sys.exit(out.stderr)
    got = {}
    for line in out.stdout.strip().splitlines():
        idx, vec, hits, rho, wide = line.split("|")
        cs = tuple(int(t) for t in vec.strip().strip("[]").split(","))
        got[int(idx) - 1] = (cs, int(hits), rho.strip(), float(wide) - float(rho))
    return [got[i] for i in range(len(polys))]


def verb_lemmas(scratch):
    t0 = time.time()
    keys, dets = {}, {}
    for k in range(1, 5):
        ncodes = 1 << (1 << k)
        seen = set()
        if k < 2:
            for code in range(ncodes):
                seen.add(charpoly(transfer(code, 1, k)))
        else:
            rows, cols, sign = block_frame(k)
            small, agree = 0, 0
            for code in range(ncodes):
                m = transfer(code, 1, k)
                seen.add(charpoly(m))
                d = det_bareiss(m)
                small += d in (-1, 0, 1)
                agree += block_det(m, rows, cols, sign) == d
            dets[k] = (ncodes, small, agree, 1 << (k - 2), block_support(k))
        keys[k] = sorted(seen)
    for k in sorted(dets):
        ncodes, small, agree, blocks, support = dets[k]
        print(
            f"det k={k} rules={ncodes} blocks={blocks} det in [-1,0,1] {small} of {ncodes} "
            f"block product agrees {agree} of {ncodes} support={int(support)}"
        )
    minpolys, radii = [], []
    for k in range(1, 5):
        got = pari_minpolys(keys[k], scratch)
        ones = sum(1 for _, n, _, _ in got if n == 1)
        tight = sum(1 for _, _, _, gap in got if abs(gap) < 1e-20)
        mp = {cs for cs, _, _, _ in got}
        rr = {rho for _, _, rho, _ in got}
        minpolys.append(len(mp))
        radii.append(len(rr))
        print(
            f"perron k={k} charpolys={len(keys[k])} minpolys={len(mp)} distinct_rho={len(rr)} "
            f"unique factor {ones} of {len(got)} no conjugate above rho {tight} of {len(got)}"
        )
    print(
        "minpolys "
        + ", ".join(str(v) for v in minpolys)
        + " and distinct rho "
        + ", ".join(str(v) for v in radii)
    )
    print(f"lemmas {time.time() - t0:.2f}s")


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else ""
    if verb == "burnside":
        verb_burnside()
    elif verb == "lemmas":
        verb_lemmas(
            Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.gettempdir())
        )
    else:
        main()
