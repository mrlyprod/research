import re
import sys
from collections import Counter, defaultdict

from sympy import Poly, Symbol, factor_list

x = Symbol("x")
SUPERGOLDEN = 1.4655712318768


def parse(path):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.startswith("ray "):
                continue
            m = re.match(
                r"ray (\d+) (\d+) states (\d+) mu ([\d.]+) width ([\d.e+-]+) match (.*?) exact \[(.*?)\] returns (\S+)( renewal-\w+)?( charpoly \[(.*?)\])?",
                line,
            )
            rows.append(
                {
                    "a": int(m.group(1)),
                    "b": int(m.group(2)),
                    "n": int(m.group(3)),
                    "rho": float(m.group(4)),
                    "match": m.group(6),
                    "exact": m.group(7).split(),
                    "returns": m.group(8),
                    "cp": [int(v) for v in m.group(11).split(",")]
                    if m.group(11)
                    else None,
                }
            )
    return rows


def v3(z):
    k = 0
    while z % 3 == 0:
        z //= 3
        k += 1
    return k


def digits3(z):
    out = []
    while z:
        out.append(z % 3)
        z //= 3
    return out


def binary3(z):
    return all(d < 2 for d in digits3(z))


def in_gasket(u, v):
    return binary3(u) and binary3(v) and binary3(u + v)


def least_multiplier(a, b, cap=3**13):
    m = 1
    while m < cap:
        if in_gasket(m * a, m * b):
            return m
        m += 1
    return None


def perron_factor(cp):
    poly = Poly(list(reversed(cp)), x)
    roots = []
    for fac, _ in factor_list(poly.as_expr())[1]:
        fp = Poly(fac, x)
        roots += [(float(r), r, fp) for r in fp.real_roots()]
    roots.sort(key=lambda t: -t[0])
    if (
        len(roots) > 1
        and roots[0][0] - roots[1][0] < 1e-6
        and not roots[0][1] > roots[1][1]
    ):
        roots[0], roots[1] = roots[1], roots[0]
    return roots[0][1], roots[0][2]


def trinomial(tag):
    k = int(tag[1:])
    return Poly(x**k - (x ** (k - 1) if tag[0] == "A" else x) - 1, x)


def profile(r):
    if r["returns"] == "cyclic":
        return None
    return tuple(tuple(map(int, s.split(":"))) for s in r["returns"].split(","))


def main(path):
    rows = parse(path)
    print(
        "occupied",
        len(rows),
        "live-state distribution",
        sorted(Counter(r["n"] for r in rows).items()),
    )
    print("no exact charpoly:", [(r["a"], r["b"]) for r in rows if r["cp"] is None])
    mismatch = 0
    for r in rows:
        top, fp = perron_factor(r["cp"])
        r["mp"] = fp
        r["top"] = float(top)
        if abs(r["top"] - r["rho"]) > 1e-9:
            mismatch += 1
            print("  LARGEST REAL ROOT DIFFERS", (r["a"], r["b"]), r["top"], r["rho"])
        if (
            r["match"] != "neither"
            and not trinomial(r["match"].split()[0]).rem(fp).is_zero
        ):
            mismatch += 1
            print(
                "  TRINOMIAL IS NOT THE MINIMAL POLYNOMIAL",
                (r["a"], r["b"]),
                r["match"],
                fp.as_expr(),
            )
    print(
        "every mu is the largest real root of its exact charpoly to 1e-9 and every matched trinomial is the minimal polynomial of that root:",
        mismatch == 0,
    )
    print(
        "distinct rates",
        len({round(r["rho"], 10) for r in rows}),
        "minimal-polynomial degrees",
        sorted(Counter(r["mp"].degree() for r in rows).items()),
    )
    neither = [r for r in rows if r["match"] == "neither"]
    print(
        "rays in neither family",
        len(neither),
        "least by height",
        min(((r["b"], r["a"]) for r in neither), default=None),
    )
    for r in sorted(neither, key=lambda r: (r["b"], r["a"])):
        print(
            "  (%d,%d) states %d rho %.10f v3(q) %d returns %s minpoly %s"
            % (
                r["a"],
                r["b"],
                r["n"],
                r["rho"],
                v3(r["a"] if r["a"] % 3 == 0 else r["b"]),
                r["returns"],
                r["mp"].as_expr(),
            )
        )
    print(
        "family B rays",
        [
            ((r["a"], r["b"]), r["match"].split()[0], r["returns"])
            for r in rows
            if r["match"][0] == "B"
        ],
    )
    single = agree = capped = 0
    for r in rows:
        m0 = least_multiplier(r["a"], r["b"])
        r["l0"] = len(digits3((r["a"] + r["b"]) * m0)) if m0 else None
        p = profile(r)
        if p and len(p) == 2:
            single += 1
            if r["l0"] is None:
                capped += 1
            elif p[1][0] == r["l0"]:
                agree += 1
    print(
        "single-return rays",
        single,
        "with k equal to the ternary length of (a + b) m0 on",
        agree,
        "and m0 past the search cap 3^13 on",
        capped,
    )
    acyclic = [r for r in rows if profile(r)]
    print(
        "acyclic rays",
        len(acyclic),
        "with the shortest return equal to that length on",
        sum(1 for r in acyclic if r["l0"] == profile(r)[1][0]),
    )
    prof = Counter(tuple(j for j, _ in profile(r)) for r in acyclic)
    print(
        "return-length profiles",
        sorted(prof.items(), key=lambda kv: (len(kv[0]), kv[0])),
    )
    grid = defaultdict(Counter)
    for r in rows:
        if r["match"] != "neither" and r["match"][0] == "A":
            p, q = (r["b"], r["a"]) if r["a"] % 3 == 0 else (r["a"], r["b"])
            q1 = q // 3 ** v3(q)
            grid[(v3(q), v3(q1 - p) if q1 != p else None)][
                int(r["match"].split()[0][1:])
            ] += 1
    print("family A: k against (v3(q), v3(q1 - p))")
    for key in sorted(grid, key=str):
        print("  ", key, dict(sorted(grid[key].items())))
    off = [r for r in rows if not (r["a"] == 1 and r["b"] == 3 ** v3(r["b"]))]
    print(
        "attainers of the supergolden",
        sorted(
            ((r["a"], r["b"]), r["a"] + r["b"])
            for r in rows
            if abs(r["rho"] - SUPERGOLDEN) < 1e-9
        ),
    )
    print("top off-shift rates")
    seen = set()
    for r in sorted(off, key=lambda r: -r["rho"]):
        key = round(r["rho"], 10)
        if key in seen:
            continue
        seen.add(key)
        if len(seen) > 6:
            break
        print(
            "  %.13f %s on %s"
            % (
                r["rho"],
                r["mp"].as_expr(),
                sorted((s["a"], s["b"]) for s in off if round(s["rho"], 10) == key),
            )
        )


if __name__ == "__main__":
    main(sys.argv[1])
