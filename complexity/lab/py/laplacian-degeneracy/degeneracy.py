import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

CODE = 7
BASE = 2
TOL = 1e-9
LEVELS = list(range(1, 9))
SWEEP = [1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05]
SECOND = 1.0 - np.sqrt(30.0) / 6.0
THIRD = 1.0 - 0.988332421566
NEAR = 1e-6


def tile(code, q):
    cells = np.array([(code >> i) & 1 for i in range(q * q)], dtype=np.uint8)
    return cells.reshape(q, q)


def fractal(code, q, level):
    unit = tile(code, q)
    out = unit
    for _ in range(1, level):
        out = np.kron(out, unit)
    return out


def and_form(grid):
    side = grid.shape[0]
    rows = np.arange(side).reshape(-1, 1)
    cols = np.arange(side).reshape(1, -1)
    return np.array_equal(grid != 0, (rows & cols) == 0)


def graph(grid):
    flat = grid.reshape(-1)
    index = np.full(flat.size, -1, dtype=np.int64)
    filled = np.flatnonzero(flat)
    index[filled] = np.arange(filled.size)
    index = index.reshape(grid.shape)
    pairs = []
    for axis in range(grid.ndim):
        low = np.take(index, np.arange(index.shape[axis] - 1), axis=axis)
        high = np.take(index, np.arange(1, index.shape[axis]), axis=axis)
        keep = (low >= 0) & (high >= 0)
        pairs.append(np.stack([low[keep], high[keep]], axis=1))
    return filled.size, np.concatenate(pairs)


def one_component(nodes, edges):
    ones = np.ones(edges.shape[0])
    adj = coo_matrix((ones, (edges[:, 0], edges[:, 1])), shape=(nodes, nodes))
    count, _ = connected_components(adj, directed=False)
    return count == 1


def laplacian(nodes, edges):
    degree = np.zeros(nodes)
    np.add.at(degree, edges[:, 0], 1.0)
    np.add.at(degree, edges[:, 1], 1.0)
    root = 1.0 / np.sqrt(degree)
    out = np.zeros((nodes, nodes))
    np.fill_diagonal(out, 1.0)
    weight = root[edges[:, 0]] * root[edges[:, 1]]
    out[edges[:, 0], edges[:, 1]] = -weight
    out[edges[:, 1], edges[:, 0]] = -weight
    return out


def spectrum(level):
    grid = fractal(CODE, BASE, level)
    if not and_form(grid):
        raise SystemExit("kronecker fractal is not the AND set at level %d" % level)
    nodes, edges = graph(grid)
    if not one_component(nodes, edges):
        raise SystemExit("graph is disconnected at level %d" % level)
    return np.sort(np.linalg.eigvalsh(laplacian(nodes, edges)))


def classes(mu, tol):
    starts = np.concatenate(([0], np.flatnonzero(np.diff(mu) > tol) + 1))
    ends = np.concatenate((starts[1:], [mu.size]))
    counts = ends - starts
    values = np.add.reduceat(mu, starts) / counts
    spreads = mu[ends - 1] - mu[starts]
    return values, counts, spreads


def locate(values, target):
    at = int(np.argmin(np.abs(values - target)))
    if abs(values[at] - target) < NEAR:
        return at
    return -1


def multiplicity(values, counts, target):
    at = locate(values, target)
    return int(counts[at]) if at >= 0 else 0


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def main():
    print(
        "object: design code %d, base %d, dimension 2, normalised Laplacian"
        % (CODE, BASE)
    )
    print("clustering tolerance: %g" % TOL)
    print()
    print("level nodes distinct degenerate repeated mult_1 mult_second max_spread")
    table = {}
    top = None
    for level in LEVELS:
        mu = spectrum(level)
        values, counts, spreads = classes(mu, TOL)
        nodes = mu.size
        big = counts > 1
        distinct = int(values.size)
        degenerate = int(big.sum())
        repeated = float(counts[big].sum()) / nodes
        unit = multiplicity(values, counts, 1.0)
        second = multiplicity(values, counts, SECOND)
        spread = float(spreads[big].max()) if degenerate else 0.0
        print(
            "%d %d %d %d %.4f %d %d %.2e"
            % (level, nodes, distinct, degenerate, repeated, unit, second, spread)
        )
        third_at = locate(values, THIRD)
        table[level] = {
            "nodes": nodes,
            "distinct": distinct,
            "degenerate": degenerate,
            "unit": unit,
            "second": second,
            "second_value": values[locate(values, SECOND)] if second else None,
            "third": int(counts[third_at]) if third_at >= 0 else 0,
            "third_value": float(values[third_at]) if third_at >= 0 else None,
            "spread": spread,
        }
        if level == LEVELS[-1]:
            top = (mu, values, counts)
    print()

    print("multiplicity of eigenvalue 1, levels 1 to 8")
    print("  measured: %s" % ", ".join(str(table[l]["unit"]) for l in LEVELS))
    print("  3^(L-1):  %s" % ", ".join(str(3 ** (l - 1)) for l in LEVELS))
    print()

    print("multiplicity of the 1 -/+ sqrt(30)/6 family, levels 3 to 8")
    print("  measured:   %s" % ", ".join(str(table[l]["second"]) for l in LEVELS[2:]))
    print("  3^(L-3)+1:  %s" % ", ".join(str(3 ** (l - 3) + 1) for l in LEVELS[2:]))
    worst = 0.0
    for level in LEVELS[1:]:
        worst = max(worst, abs(table[level]["second_value"] - SECOND))
    print("  1 - sqrt(30)/6 = %.12f" % SECOND)
    print("  1 + sqrt(30)/6 = %.12f" % (2.0 - SECOND))
    print("  worst deviation over levels 2 to 8: %.2e" % worst)
    print()

    print("third family, levels 4 to 8")
    print("  eigenvalue sought: 1 - %.12f" % (1.0 - THIRD))
    print("  measured at level 8: 1 - %.12f" % (1.0 - table[8]["third_value"]))
    print("  measured:   %s" % ", ".join(str(table[l]["third"]) for l in LEVELS[3:]))
    print("  3^(L-4)+1:  %s" % ", ".join(str(3 ** (l - 4) + 1) for l in LEVELS[3:]))
    _, values8, counts8 = top
    print("  classes at level 8 with multiplicity 82: %d" % int((counts8 == 82).sum()))
    print()

    print("counting fits")
    print(
        "  distinct:        %s" % ", ".join(str(table[l]["distinct"]) for l in LEVELS)
    )
    print("  2*Fib(2L)+1:     %s" % ", ".join(str(2 * fib(2 * l) + 1) for l in LEVELS))
    print(
        "  degenerate:      %s"
        % ", ".join(str(table[l]["degenerate"]) for l in LEVELS[1:])
    )
    print(
        "  2*Fib(2L-3)-1:   %s"
        % ", ".join(str(2 * fib(2 * l - 3) - 1) for l in LEVELS[1:])
    )
    print()

    mu8, values8, counts8 = top
    print("reading the integers off a floating-point spectrum, level 8")
    print("  widest degenerate class spans %.2e" % table[8]["spread"])
    at = locate(values8, 1.0)
    print("  nearest distinct eigenvalue below 1: %.10f" % values8[at - 1])
    print("  nearest distinct eigenvalue above 1: %.10f" % values8[at + 1])
    print(
        "  isolation of the class at 1: %.3f below, %.3f above"
        % (1.0 - values8[at - 1], values8[at + 1] - 1.0)
    )
    print()

    stop = 3**9
    print("why level 9 stops under this method")
    print("  nodes at level 9: %d" % stop)
    print(
        "  dense matrix before eigensolver workspace: %.2f GB" % (stop * stop * 8 / 1e9)
    )
    print()

    print("tolerance sweep at level 8")
    print("  tol distinct mult_1 mult_second")
    for tol in SWEEP:
        values, counts, _ = classes(mu8, tol)
        print(
            "  %.0e %d %d %d"
            % (
                tol,
                values.size,
                multiplicity(values, counts, 1.0),
                multiplicity(values, counts, SECOND),
            )
        )


main()
