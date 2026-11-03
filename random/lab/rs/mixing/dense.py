import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

N = 3
BITS = N * N
M = 1 << BITS
HOODS = {"carpet": 495, "net": 186, "void": 341}
SOURCES = {
    "primes": [2, 3, 5, 7],
    "binary": [1, 2, 4, 8],
    "fibonacci": [0, 1, 2, 3, 5, 8],
    "odds": [1, 3, 5, 7],
}
SPREADS = {
    "id": lambda i, j: (i, j),
    "shift": lambda i, j: (i, (j + 1) % N),
    "drop": lambda i, j: ((i + 1) % N, j),
    "slide": lambda i, j: ((i + 1) % N, (j + 1) % N),
    "turn": lambda i, j: (j, N - 1 - i),
    "flip": lambda i, j: (j, i),
}
MATS = [
    (1, 0, 0, 1),
    (0, -1, 1, 0),
    (-1, 0, 0, -1),
    (0, 1, -1, 0),
    (0, 1, 1, 0),
    (0, -1, -1, 0),
    (-1, 0, 0, 1),
    (1, 0, 0, -1),
]


def cell_perm(f):
    return [N * f(p // N, p % N)[0] + f(p // N, p % N)[1] for p in range(BITS)]


def symmetries():
    out = []
    for a, b, c, d in MATS:
        for ti in range(N):
            for tj in range(N):
                out.append(
                    [
                        ((a * (p // N) + b * (p % N) + ti) % N) * N
                        + (c * (p // N) + d * (p % N) + tj) % N
                        for p in range(BITS)
                    ]
                )
    return out


def image(perm, mask):
    return sum(1 << perm[p] for p in range(BITS) if mask >> p & 1)


def reps(spread):
    sp = cell_perm(SPREADS[spread])
    central = [
        s for s in symmetries() if all(s[sp[p]] == sp[s[p]] for p in range(BITS))
    ]
    seen, out = set(), []
    for mask in range(1, M - 1):
        canon = min(image(s, mask) for s in central)
        if canon not in seen:
            seen.add(canon)
            out.append(mask)
    return out


def table(fill, hood, source, spread):
    offs = [(b // 3 - 1, b % 3 - 1) for b in range(9) if hood >> b & 1 and b != 4]
    void = (M - 1) & ~fill
    reads = []
    for a in range(BITS):
        if fill >> a & 1:
            i, j = a // N, a % N
            reads.append(
                (
                    a,
                    sum(1 << (((i + di) % N) * N + (j + dj) % N) for di, dj in offs)
                    & void,
                )
            )
    x = np.arange(M)
    bits = (x[:, None] >> np.arange(BITS)) & 1
    mess = np.zeros(M, dtype=np.int64)
    for a, mask in reads:
        cols = [q for q in range(BITS) if mask >> q & 1]
        count = bits[:, cols].sum(axis=1) if cols else np.zeros(M, dtype=np.int64)
        mess |= np.isin(count, source).astype(np.int64) << a
    y = x ^ mess
    sp = cell_perm(SPREADS[spread])
    out = np.zeros(M, dtype=np.int64)
    for p in range(BITS):
        out |= ((y >> p) & 1) << sp[p]
    return out


def kernel(h):
    x = np.arange(M)
    d = h[x[None, :]] ^ h[x[None, :] ^ x[:, None]]
    k = np.bincount((x[:, None] * M + d).ravel(), minlength=M * M).reshape(M, M)
    return k[1:, 1:] / M


def sigma2(q):
    return np.linalg.svd(q, compute_uv=False)[1]


def curve(p, eps):
    n = p.shape[0]
    d, power = [1 - 1 / n], np.eye(n)
    while d[-1] > eps and len(d) < 400:
        power = power @ p
        d.append(0.5 * np.abs(power - 1 / n).sum(axis=1).max())
    return d


def first(d, eps):
    return next(i + 1 for i, v in enumerate(d) if v < eps)


def main():
    counts, rows = {"nonlinear": 0, "affine": 0}, []
    total = 0
    for spread in SPREADS:
        masks = reps(spread)
        for hood, code in HOODS.items():
            for source, values in SOURCES.items():
                for fill in masks:
                    total += 1
                    p = kernel(table(fill, code, values, spread))
                    comps, _ = connected_components(
                        csr_matrix(p > 0), directed=True, connection="strong"
                    )
                    if comps != 1:
                        continue
                    kind = "affine" if source == "odds" else "nonlinear"
                    counts[kind] += 1
                    if kind == "affine":
                        continue
                    d = curve(p, 0.01)
                    if counts[kind] % 7 != 1:
                        rows.append(
                            (
                                None,
                                None,
                                None,
                                first(d, 0.25),
                                first(d, 0.01),
                                f"{fill:#06x} {hood} {source} {spread}",
                            )
                        )
                        continue
                    lam = np.sort(np.abs(np.linalg.eigvals(p)))[::-1]
                    powers, q = [], np.eye(p.shape[0])
                    while not powers or powers[-1] > 1 - 1e-9:
                        q = q @ p
                        powers.append(sigma2(q))
                    rows.append(
                        (
                            lam[1],
                            len(powers),
                            powers[-1],
                            first(d, 0.25),
                            first(d, 0.01),
                            f"{fill:#06x} {hood} {source} {spread}",
                        )
                    )
    print(
        f"  n 3: {total} cells, irreducible nonlinear {counts['nonlinear']}, affine {counts['affine']}"
    )
    print(
        f"  spectra on {sum(r[0] is not None for r in rows)} of {len(rows)} cells (every 7th), curves on all"
    )
    for name, col in (
        ("|lambda2|", 0),
        ("k*", 1),
        ("sigma2(P^k*)", 2),
        ("rounds to 1/4", 3),
        ("rounds to 0.01", 4),
    ):
        v = sorted(r[col] for r in rows if r[col] is not None)
        print(f"  {name}: min {v[0]:.6f} median {v[len(v) // 2]:.6f} max {v[-1]:.6f}")
    slow = max(rows, key=lambda r: r[4])
    print(f"  slowest: {slow[5]}: rounds to 1/4 {slow[3]}, to 0.01 {slow[4]}")
    for r in rows:
        if r[0] is not None or len(sys.argv) > 1:
            print(
                f"  {r[5]}: |lambda2| {r[0]}, k* {r[1]}, sigma2(P^k*) {r[2]}, rounds to 1/4 {r[3]}, to 0.01 {r[4]}"
            )


main()
