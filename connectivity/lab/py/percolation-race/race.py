import statistics
import time
from collections import deque
from itertools import product
from random import Random

import numpy as np

NUMPY_SEED = 20260725
PYTHON_SEED_BASE = 1000


def corners(base, dimension, drop):
    return [c for c in product(range(base), repeat=dimension) if not drop(c)]


DESIGNS = {
    "gasket": (2, 2, [(0, 0), (0, 1), (1, 0)]),
    "diagonal": (2, 2, [(0, 0), (1, 1)]),
    "antidiagonal": (2, 2, [(0, 1), (1, 0)]),
    "seven-of-eight": (2, 3, corners(2, 3, lambda c: c == (1, 1, 1))),
    "carpet": (3, 2, corners(3, 2, lambda c: c == (1, 1))),
    "sponge": (3, 3, corners(3, 3, lambda c: c.count(1) > 1)),
}

SAME = {"antidiagonal": "diagonal"}

PASS_A = [
    ("gasket", 5, 400),
    ("gasket", 6, 400),
    ("antidiagonal", 5, 400),
    ("antidiagonal", 6, 400),
    ("seven-of-eight", 3, 400),
    ("seven-of-eight", 4, 400),
    ("carpet", 3, 400),
    ("carpet", 4, 200),
    ("sponge", 3, 200),
    ("sponge", 4, 25),
]

PASS_B = [
    ("gasket", 5, 400),
    ("gasket", 6, 400),
    ("gasket", 7, 200),
    ("gasket", 8, 100),
    ("diagonal", 5, 400),
    ("diagonal", 6, 400),
    ("diagonal", 7, 200),
    ("seven-of-eight", 3, 400),
    ("seven-of-eight", 4, 200),
    ("carpet", 3, 400),
    ("carpet", 4, 200),
    ("sponge", 3, 100),
    ("sponge", 4, 20),
]


def kron_power(base, dimension, cells, level):
    tile = np.zeros((base,) * dimension, dtype=np.uint8)
    for cell in cells:
        tile[cell] = 1
    out = tile
    for _ in range(level - 1):
        out = np.kron(out, tile)
    return out.astype(bool)


def row_strides(shape):
    out = [1] * len(shape)
    for axis in range(len(shape) - 2, -1, -1):
        out[axis] = out[axis + 1] * shape[axis + 1]
    return out


def edge_ranks(occupied, grid, cells):
    shape = occupied.shape
    stride = row_strides(shape)
    lows = []
    for axis in range(len(shape)):
        lo = [slice(None)] * len(shape)
        hi = [slice(None)] * len(shape)
        lo[axis] = slice(0, shape[axis] - 1)
        hi[axis] = slice(1, shape[axis])
        both = occupied[tuple(lo)] & occupied[tuple(hi)]
        left = grid[tuple(lo)][both]
        lows.append((left, left + stride[axis]))
    a = np.concatenate([pair[0] for pair in lows])
    b = np.concatenate([pair[1] for pair in lows])
    return np.searchsorted(cells, a), np.searchsorted(cells, b)


def union_stats(n, ea, eb):
    parent = list(range(n))
    weight = [1] * n
    count = n
    largest = 1
    for i, j in zip(ea, eb):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        while parent[j] != j:
            parent[j] = parent[parent[j]]
            j = parent[j]
        if i == j:
            continue
        if weight[i] < weight[j]:
            i, j = j, i
        parent[j] = i
        weight[i] += weight[j]
        count -= 1
        largest = max(largest, weight[i])
    return count, largest


def kron_stats(occupied, grid, dimension):
    cells = np.flatnonzero(occupied.ravel())
    n = int(cells.size)
    ranks_a, ranks_b = edge_ranks(occupied, grid, cells)
    edges = int(ranks_a.size)
    count, largest = union_stats(n, ranks_a.tolist(), ranks_b.tolist())
    return n, count, largest / n, (2 * dimension * n - 2 * edges) / n


def substitution(base, dimension, cells, level):
    out = [(0,) * dimension]
    for _ in range(level):
        out = [
            tuple(c[a] * base + f[a] for a in range(dimension))
            for c in out
            for f in cells
        ]
    return out


def digit_rule(base, dimension, cells, level):
    keep = set(cells)
    side = base**level
    out = []
    for index in range(side**dimension):
        coord = []
        left = index
        for _ in range(dimension):
            coord.append(left % side)
            left //= side
        ok = True
        for position in range(level):
            if tuple((c // base**position) % base for c in coord) not in keep:
                ok = False
                break
        if ok:
            out.append(index)
    return out


def padded_steps(side, dimension):
    return [(side + 2) ** a for a in range(dimension)]


def place(coords, side, dimension, step):
    occupied = bytearray((side + 2) ** dimension)
    flat = [sum((c[a] + 1) * step[a] for a in range(dimension)) for c in coords]
    for index in flat:
        occupied[index] = 1
    return occupied, flat


def spread_out(picks, side, dimension, step):
    out = []
    for index in picks:
        left = index
        total = 0
        for a in range(dimension):
            total += (left % side + 1) * step[a]
            left //= side
        out.append(total)
    return out


def bfs_stats(occupied, flat, step, dimension):
    n = len(flat)
    edges = 0
    for index in flat:
        for s in step:
            if occupied[index + s]:
                edges += 1
    moves = step + [-s for s in step]
    seen = bytearray(len(occupied))
    components = 0
    largest = 0
    for start in flat:
        if seen[start]:
            continue
        components += 1
        seen[start] = 1
        queue = deque([start])
        size = 0
        while queue:
            index = queue.popleft()
            size += 1
            for s in moves:
                neighbour = index + s
                if occupied[neighbour] and not seen[neighbour]:
                    seen[neighbour] = 1
                    queue.append(neighbour)
        largest = max(largest, size)
    return n, components, largest / n, (2 * dimension * n - 2 * edges) / n


def run_a(name, level, seeds):
    base, dimension, cells = DESIGNS[name]
    design = kron_power(base, dimension, cells, level)
    grid = np.arange(design.size, dtype=np.int64).reshape(design.shape)
    n, dc, df, db = kron_stats(design, grid, dimension)
    rng = np.random.default_rng(NUMPY_SEED)
    comps, frac, bound = [], [], []
    for _ in range(seeds):
        flat = np.zeros(design.size, dtype=bool)
        flat[rng.choice(design.size, n, replace=False)] = True
        _, c, f, b = kron_stats(flat.reshape(design.shape), grid, dimension)
        comps.append(float(c))
        frac.append(f)
        bound.append(b)
    return pack(
        name,
        level,
        base**level,
        dimension,
        n,
        design.size,
        seeds,
        dc,
        df,
        db,
        comps,
        frac,
        bound,
        lambda v: (float(np.mean(v)), float(np.std(v, ddof=1))),
    )


def draw_b(name, level, seeds):
    base, dimension, cells = DESIGNS[name]
    side = base**level
    total = side**dimension
    step = padded_steps(side, dimension)
    design = substitution(base, dimension, cells, level)
    occupied, flat = place(design, side, dimension, step)
    n, dc, df, db = bfs_stats(occupied, flat, step, dimension)
    comps, frac, bound = [], [], []
    for seed in range(seeds):
        picks = Random(PYTHON_SEED_BASE + seed).sample(range(total), n)
        holes = spread_out(picks, side, dimension, step)
        blank = bytearray((side + 2) ** dimension)
        for index in holes:
            blank[index] = 1
        _, c, f, b = bfs_stats(blank, holes, step, dimension)
        comps.append(float(c))
        frac.append(f)
        bound.append(b)
    return pack(
        name,
        level,
        side,
        dimension,
        n,
        total,
        seeds,
        dc,
        df,
        db,
        comps,
        frac,
        bound,
        lambda v: (statistics.mean(v), statistics.stdev(v)),
    )


def pack(
    name,
    level,
    side,
    dimension,
    n,
    total,
    seeds,
    dc,
    df,
    db,
    comps,
    frac,
    bound,
    moments,
):
    return {
        "name": name,
        "level": level,
        "side": side,
        "dimension": dimension,
        "cells": n,
        "density": n / total,
        "seeds": seeds,
        "design": (dc, df, db),
        "comps": moments(comps),
        "frac": moments(frac),
        "bound": moments(bound),
        "tie_comps": sum(1 for v in comps if v <= dc),
        "tie_bound": sum(1 for v in bound if v >= db),
    }


def show(result):
    grid = "x".join([str(result["side"])] * result["dimension"])
    print(
        "  {:<15} L={} grid {:<12} cells {:>7} density {:.4f} seeds {}".format(
            result["name"],
            result["level"],
            grid,
            result["cells"],
            result["density"],
            result["seeds"],
        )
    )
    dc, df, db = result["design"]
    print(
        "    components  design {:>10}   random {:12.4f} +/- {:.4f}   ties {}/{}".format(
            dc,
            result["comps"][0],
            result["comps"][1],
            result["tie_comps"],
            result["seeds"],
        )
    )
    print(
        "    largest     design {:>10.4f}   random {:12.4f} +/- {:.4f}".format(
            df, result["frac"][0], result["frac"][1]
        )
    )
    print(
        "    boundary    design {:>10.4f}   random {:12.4f} +/- {:.4f}   ties {}/{}".format(
            db,
            result["bound"][0],
            result["bound"][1],
            result["tie_bound"],
            result["seeds"],
        )
    )


def main():
    start = time.time()
    print("DOMAIN")
    print("  pass A: 2D grids to 81x81, 3D to 81^3, 400 seeds thinning to 25")
    print("  pass B: 2D grids to 256x256, 3D to 81^3, 400 seeds thinning to 20")
    print(
        "  ties: random draws at or below the design count, at or above the design boundary"
    )
    print()

    print("SUBSTITUTION AGAINST DIGIT RULE")
    agree = True
    for name in ("gasket", "diagonal", "seven-of-eight", "carpet", "sponge"):
        base, dimension, cells = DESIGNS[name]
        for level in (1, 2, 3):
            side = base**level
            a = sorted(digit_rule(base, dimension, cells, level))
            raw = substitution(base, dimension, cells, level)
            b = sorted({sum(c[x] * side**x for x in range(dimension)) for c in raw})
            same = a == b and len(a) == len(cells) ** level
            agree = agree and same
            print(f"  {name:<15} L={level} cells {len(a):>7}  routes agree {same}")
    print(f"  every route agrees and the fill count is multiplicative: {agree}")
    print()

    print("PASS A - KRONECKER POWERS, UNION-FIND, NUMPY PCG64")
    a_results = {}
    for name, level, seeds in PASS_A:
        result = run_a(name, level, seeds)
        show(result)
        a_results[(SAME.get(name, name), level)] = result
    print()

    print("PASS B - SUBSTITUTION, BREADTH-FIRST SEARCH, PYTHON MERSENNE TWISTER")
    b_results = {}
    for name, level, seeds in PASS_B:
        result = draw_b(name, level, seeds)
        show(result)
        b_results[(SAME.get(name, name), level)] = result
    print()

    print("DISPERSING EXTREME - COMPONENTS ARE NOT THE METRIC")
    for level in (5, 6, 7):
        result = b_results[("diagonal", level)]
        print(
            "  diagonal L={} side {:>4} cells {:>5} seeds {:>3}  design comps {:>5}   random comps {:10.4f} +/- {:.4f}".format(
                level,
                result["side"],
                result["cells"],
                result["seeds"],
                result["design"][0],
                result["comps"][0],
                result["comps"][1],
            )
        )
    print()

    print("MAXIMUM BOUNDARY")
    for key in (
        ("diagonal", 6),
        ("gasket", 6),
        ("carpet", 4),
        ("sponge", 3),
        ("seven-of-eight", 4),
    ):
        result = b_results[key]
        print(
            "  {:<15} L={} boundary per cell {:.4f} of the maximum {}".format(
                key[0], key[1], result["design"][2], 2 * result["dimension"]
            )
        )
    print()

    print("THE TWO PASSES AGAINST EACH OTHER")
    worst = 0.0
    for key in sorted(set(a_results) & set(b_results)):
        left, right = a_results[key], b_results[key]
        for field in ("comps", "frac", "bound"):
            gap = abs(left[field][0] - right[field][0])
            sd = max(left[field][1], right[field][1])
            worst = max(worst, gap / sd if sd else 0.0)
            print(
                f"  {key[0]:<15} L={key[1]} {field:<6} pass A {left[field][0]:12.4f}   pass B {right[field][0]:12.4f}   gap {gap / sd if sd else 0.0:.4f} sd"
            )
    print(f"  every comparison agrees to within one standard deviation: {worst <= 1.0}")
    print(f"  the widest gap is {worst:.4f} standard deviations")
    print()

    print("WITNESSES")
    gasket8 = b_results[("gasket", 8)]
    gasket5 = b_results[("gasket", 5)]
    sponge4 = b_results[("sponge", 4)]
    print(
        "  gasket 256x256   design comps {}   random {:.2f} +/- {:.2f}".format(
            gasket8["design"][0], gasket8["comps"][0], gasket8["comps"][1]
        )
    )
    print(
        "  gasket 256x256   design boundary {:.4f}   random {:.4f} +/- {:.4f}".format(
            gasket8["design"][2], gasket8["bound"][0], gasket8["bound"][1]
        )
    )
    print(
        "  sponge 81^3      design comps {}   random {:.2f} +/- {:.2f}".format(
            sponge4["design"][0], sponge4["comps"][0], sponge4["comps"][1]
        )
    )
    print(
        "  gasket 32x32     random draws reaching one component {}/{}".format(
            gasket5["tie_comps"], gasket5["seeds"]
        )
    )
    print()
    print(f"  every check passed: {agree and worst <= 1.0}")
    print(f"  wall {time.time() - start:.1f} s")
    if not agree or worst > 1.0:
        raise SystemExit(1)


main()
