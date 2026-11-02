import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ZERO = 1e-9
MULT = 1e-8
SHOW = 12
BINS = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])


def carpet(level):
    tile = np.ones((3, 3), dtype=np.uint8)
    tile[1, 1] = 0
    grid = np.ones((1, 1), dtype=np.uint8)
    for _ in range(level):
        grid = np.kron(grid, tile)
    assert int(grid.sum()) == 8**level
    return grid


def edges(grid):
    index = -np.ones(grid.shape, dtype=np.int64)
    filled = np.argwhere(grid == 1)
    index[filled[:, 0], filled[:, 1]] = np.arange(len(filled))
    right = (grid[:, :-1] == 1) & (grid[:, 1:] == 1)
    down = (grid[:-1, :] == 1) & (grid[1:, :] == 1)
    a = np.concatenate([index[:, :-1][right], index[:-1, :][down]])
    b = np.concatenate([index[:, 1:][right], index[1:, :][down]])
    return len(filled), a, b


def laplacian(n, a, b):
    lap = np.zeros((n, n))
    lap[a, b] = -1.0
    lap[b, a] = -1.0
    lap[np.arange(n), np.arange(n)] = -lap.sum(axis=1)
    return lap


def nodal(n, a, b, v):
    sign = np.where(v > ZERO, 1, np.where(v < -ZERO, -1, 0))
    keep = (sign[a] == sign[b]) & (sign[a] != 0)
    graph = coo_matrix((np.ones(int(keep.sum())), (a[keep], b[keep])), shape=(n, n))
    parts, _ = connected_components(graph, directed=False)
    zeros = int((sign == 0).sum())
    return parts - zeros, zeros


def clusters(values):
    n = len(values)
    last = np.zeros(n, dtype=np.int64)
    mult = np.zeros(n, dtype=np.int64)
    start = 0
    for i in range(1, n + 1):
        if i == n or values[i] - values[i - 1] > MULT:
            last[start:i] = i
            mult[start:i] = i - start
            start = i
    return last, mult


def report(name, n, a, b, values, vectors):
    k = np.arange(1, n + 1)
    nu = np.zeros(n, dtype=np.int64)
    zeros = np.zeros(n, dtype=np.int64)
    for i in range(n):
        nu[i], zeros[i] = nodal(n, a, b, vectors[:, i])
    last, mult = clusters(values)
    ratio = nu / k
    hist = np.histogram(np.minimum(ratio, 1.0), bins=BINS)[0]
    violations = nu > k
    degenerate_k = int((mult > 1).sum())
    degenerate_classes = len(set(last[mult > 1].tolist()))
    classes = len(set(last.tolist()))
    tail = ratio[1:]
    print(
        f"{name}: {n} nodes, {len(a)} edges, lambda_2 = {values[1]:.10f}, lambda_max = {values[-1]:.10f}"
    )
    print(f"  nu_k, k = 1..{SHOW}: {nu[:SHOW].tolist()}")
    print(f"  r_k, k = 1..{SHOW}: {mult[:SHOW].tolist()}")
    print(f"  zero cells, k = 1..{SHOW}: {zeros[:SHOW].tolist()}")
    print(
        f"  max nu_k = {nu.max()} at k = {int(k[nu.argmax()])}; over k >= 2 max nu_k/k = {tail.max():.6f} at k = {int(tail.argmax()) + 2}, mean nu_k/k = {tail.mean():.6f}"
    )
    print(
        f"  equality nu_k = k at k = {k[nu == k][:SHOW].tolist()}, {int((nu == k).sum())} times"
    )
    print(
        f"  Courant nu_k <= k: {int((~violations).sum())} of {n}, fraction {(~violations).mean():.6f}; violations {int(violations.sum())}, of which at degenerate k {int((mult[violations] > 1).sum())}"
    )
    print(
        f"  DGLS nu_k <= k + r_k - 1: {int((nu <= last).sum())} of {n}, fraction {(nu <= last).mean():.6f}"
    )
    print(
        f"  eigenvalue classes {classes}, degenerate classes {degenerate_classes}, degenerate indices {degenerate_k} of {n}, fraction {degenerate_k / n:.6f}, max multiplicity {int(mult.max())}"
    )
    print(f"  eigenvectors with a zero cell: {int((zeros > 0).sum())}")
    top = vectors[:, -1]
    sign = np.sign(top) * (np.abs(top) > ZERO)
    print(
        f"  top eigenvector: zero cells {int(zeros[-1])}, min |v| = {np.abs(top).min():.3e}, edges joining two nonzero cells of one sign: {int(((sign[a] == sign[b]) & (sign[a] != 0)).sum())}"
    )
    print(
        f"  nu_k/k in [0,0.2) [0.2,0.4) [0.4,0.6) [0.6,0.8) [0.8,1]: {hist.tolist()}, above 1: {int((ratio > 1).sum())}"
    )
    return nu, mult, last


def mixing(n, a, b, vectors, mult):
    found = []
    for i in np.flatnonzero(mult == 2)[::2]:
        u, v = vectors[:, i], vectors[:, i + 1]
        counts = [
            nodal(n, a, b, w)[0]
            for w in (u, v, (u + v) / np.sqrt(2), (u - v) / np.sqrt(2))
        ]
        if len(set(counts)) > 1:
            found.append((i + 1, counts))
    print(
        f"  double eigenvalues where the two returned vectors, their sum and their difference disagree in nu: {len(found)} of {int((mult == 2).sum()) // 2}; first at k = {found[0][0] if found else None}, counts {found[0][1] if found else None}"
    )


def separable(side):
    i = np.arange(side)
    lam = 2.0 - 2.0 * np.cos(np.pi * np.arange(side) / side)
    cos = np.cos(np.pi * np.outer(np.arange(side), i + 0.5) / side)
    items = sorted((lam[p] + lam[q], p, q) for p in range(side) for q in range(side))
    vectors = np.empty((side * side, side * side))
    expect = np.empty(side * side, dtype=np.int64)
    for col, (_, p, q) in enumerate(items):
        vectors[:, col] = np.outer(cos[p], cos[q]).ravel()
        expect[col] = (p + 1) * (q + 1)
    return np.array([t[0] for t in items]), vectors, expect


def run_carpet(level):
    t = time.time()
    n, a, b = edges(carpet(level))
    values, vectors = np.linalg.eigh(laplacian(n, a, b))
    assert (
        connected_components(
            coo_matrix((np.ones(len(a)), (a, b)), shape=(n, n)), directed=False
        )[0]
        == 1
    )
    _nu, mult, _last = report(f"carpet L = {level}", n, a, b, values, vectors)
    mixing(n, a, b, vectors, mult)
    print(f"  {time.time() - t:.1f} s")
    print()


def run_grid(side):
    t = time.time()
    n, a, b = edges(np.ones((side, side), dtype=np.uint8))
    values, vectors = np.linalg.eigh(laplacian(n, a, b))
    _nu, mult, _last = report(
        f"grid {side} x {side}, numpy basis", n, a, b, values, vectors
    )
    mixing(n, a, b, vectors, mult)
    svalues, svectors, expect = separable(side)
    assert np.abs(svalues - values).max() < 1e-9
    snu, _smult, _slast = report(
        f"grid {side} x {side}, separable basis", n, a, b, svalues, svectors
    )
    print(
        f"  separable nu_k equals (p + 1)(q + 1) for all k: {bool((snu == expect).all())}"
    )
    print(f"  {time.time() - t:.1f} s")
    print()


def main():
    t = time.time()
    print(
        f"domain: carpet L = 3, 4 on 4-neighbour adjacency, combinatorial Laplacian D - A, numpy.linalg.eigh; zero tolerance {ZERO:g}, multiplicity tolerance {MULT:g}; controls grid 22 x 22 and 64 x 64"
    )
    print(
        "nu_k counts strong nodal domains of the k-th returned eigenvector; on a degenerate eigenvalue the count depends on the basis, so every nu_k below is a fact about numpy's returned basis only"
    )
    print()
    for level in (3, 4):
        run_carpet(level)
    for side in (22, 64):
        run_grid(side)
    print(f"total {time.time() - t:.1f} s")


if __name__ == "__main__":
    main()
