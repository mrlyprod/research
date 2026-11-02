import sys
from collections import deque

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, maximum_bipartite_matching

TOP = 10
DENSE = 8
EXACT = 4
TOL = 1e-9


def nodes(level):
    side = 1 << level
    grid = np.arange(side)
    keep = (grid[:, None] & grid[None, :]) == 0
    rows, cols = np.nonzero(keep)
    return rows.astype(np.int64), cols.astype(np.int64)


def edges(rows, cols, side):
    index = {int(r) * side + int(c): k for k, (r, c) in enumerate(zip(rows, cols))}
    out = []
    for k, (r, c) in enumerate(zip(rows, cols)):
        for dr, dc in ((1, 0), (0, 1)):
            nr, nc = int(r) + dr, int(c) + dc
            if nr < side and nc < side:
                other = index.get(nr * side + nc)
                if other is not None:
                    out.append((k, other))
    return np.array(out, dtype=np.int64)


def build(level):
    rows, cols = nodes(level)
    return rows, cols, edges(rows, cols, 1 << level)


def adjacency(count, pairs):
    ones = np.ones(len(pairs))
    half = coo_matrix((ones, (pairs[:, 0], pairs[:, 1])), shape=(count, count))
    return half.tocsr() + half.T.tocsr()


def is_tree(count, pairs):
    if len(pairs) != count - 1:
        return False
    parts, _ = connected_components(adjacency(count, pairs), directed=False)
    return parts == 1


def stripped(count, pairs, gone=()):
    nbr = [[] for _ in range(count)]
    for a, b in pairs:
        nbr[a].append(b)
        nbr[b].append(a)
    dead = [False] * count
    for g in gone:
        dead[g] = True
    deg = [sum(not dead[x] for x in nbr[v]) for v in range(count)]
    queue = deque(v for v in range(count) if not dead[v] and deg[v] == 1)
    size = 0
    while queue:
        u = queue.popleft()
        if dead[u]:
            continue
        w = next((x for x in nbr[u] if not dead[x]), None)
        if w is None:
            continue
        size += 1
        dead[u] = dead[w] = True
        for x in nbr[w]:
            if not dead[x]:
                deg[x] -= 1
                if deg[x] == 1:
                    queue.append(x)
    return size


def hopcroft_karp(count, pairs, parity, gone=()):
    skip = set(gone)
    even = [v for v in range(count) if parity[v] == 0 and v not in skip]
    odd = [v for v in range(count) if parity[v] == 1 and v not in skip]
    row = {v: k for k, v in enumerate(even)}
    col = {v: k for k, v in enumerate(odd)}
    hits = []
    for a, b in pairs:
        if a in skip or b in skip:
            continue
        if parity[a] == 1:
            a, b = b, a
        hits.append((row[a], col[b]))
    if not hits:
        return 0
    data = np.ones(len(hits))
    mat = coo_matrix(
        (data, ([h[0] for h in hits], [h[1] for h in hits])),
        shape=(len(even), len(odd)),
    ).tocsr()
    match = maximum_bipartite_matching(mat, perm_type="column")
    return int((match >= 0).sum())


def blocks(level, rows, cols, pairs):
    h = 1 << (level - 1)
    block = (rows >= h).astype(int) * 2 + (cols >= h).astype(int)
    cross = {}
    for a, b in pairs:
        if block[a] != block[b]:
            key = (int(min(block[a], block[b])), int(max(block[a], block[b])))
            cross.setdefault(key, []).append((int(a), int(b)))
    return h, block, cross


def lattice(rows, cols, pairs):
    return {
        tuple(sorted(((int(rows[a]), int(cols[a])), (int(rows[b]), int(cols[b])))))
        for a, b in pairs
    }


def copies(level, rows, cols, pairs, below):
    h = 1 << (level - 1)
    shifts = {0: (0, 0), 1: (0, h), 2: (h, 0)}
    inside = {0: set(), 1: set(), 2: set()}
    for p, q in lattice(rows, cols, pairs):
        pb = (p[0] >= h) * 2 + (p[1] >= h)
        qb = (q[0] >= h) * 2 + (q[1] >= h)
        if pb == qb and pb in inside:
            dr, dc = shifts[pb]
            inside[pb].add(
                tuple(sorted(((p[0] - dr, p[1] - dc), (q[0] - dr, q[1] - dc))))
            )
    return all(inside[k] == below for k in inside)


def recursion(level, rows, cols, pairs, below):
    h, _block, cross = blocks(level, rows, cols, pairs)
    if not copies(level, rows, cols, pairs, below):
        return False
    if set(cross) != {(0, 1), (0, 2)}:
        return False
    if any(len(v) != 1 for v in cross.values()):
        return False
    link_b = cross[(0, 1)][0]
    link_c = cross[(0, 2)][0]
    ends_b = {(int(rows[v]), int(cols[v])) for v in link_b}
    ends_c = {(int(rows[v]), int(cols[v])) for v in link_c}
    return ends_b == {(0, h - 1), (0, h)} and ends_c == {(h - 1, 0), (h, 0)}


def matching_polynomial(count, pairs):
    nbr = [[] for _ in range(count)]
    for a, b in pairs:
        nbr[a].append(b)
        nbr[b].append(a)
    parent = [-1] * count
    order = []
    seen = [False] * count
    seen[0] = True
    stack = [0]
    while stack:
        v = stack.pop()
        order.append(v)
        for w in nbr[v]:
            if not seen[w]:
                seen[w] = True
                parent[w] = v
                stack.append(w)

    def mul(p, q):
        out = [0] * (len(p) + len(q) - 1)
        for i, x in enumerate(p):
            for j, y in enumerate(q):
                out[i + j] += x * y
        return out

    def add(p, q):
        size = max(len(p), len(q))
        return [
            (p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
            for i in range(size)
        ]

    free = [None] * count
    bare = [None] * count
    for v in reversed(order):
        kids = [w for w in nbr[v] if parent[w] == v]
        g = [1]
        for w in kids:
            g = mul(g, free[w])
        total = g
        for w in kids:
            rest = [1]
            for x in kids:
                if x != w:
                    rest = mul(rest, free[x])
            total = add(total, [0] + mul(bare[w], rest))
        bare[v] = g
        free[v] = total
    return free[order[0]]


def characteristic(count, pairs):
    mat = np.zeros((count, count), dtype=object)
    mat[:] = 0
    for a, b in pairs:
        mat[a, b] = 1
        mat[b, a] = 1
    ident = np.zeros((count, count), dtype=object)
    ident[:] = 0
    for i in range(count):
        ident[i, i] = 1
    coef = [1]
    cur = ident.copy()
    for k in range(1, count + 1):
        prod = mat.dot(cur)
        c = -int(sum(prod[i, i] for i in range(count))) // k
        coef.append(c)
        cur = prod + c * ident
    return coef


def exact_check(count, pairs):
    poly = matching_polynomial(count, pairs)
    char = characteristic(count, pairs)
    want = [0] * (count + 1)
    for k, m in enumerate(poly):
        want[2 * k] = (-1) ** k * m
    return char == want, poly


def nullity_adjacency(count, pairs):
    mat = adjacency(count, pairs).toarray()
    vals = np.linalg.eigvalsh(mat)
    zero = np.abs(vals) < TOL
    gap = float(np.abs(vals[~zero]).min()) if (~zero).any() else float("nan")
    return int(zero.sum()), gap, mat


def multiplicity_one(mat):
    deg = mat.sum(axis=1)
    root = 1.0 / np.sqrt(deg)
    lap = np.eye(len(deg)) - root[:, None] * mat * root[None, :]
    vals = np.linalg.eigvalsh(lap)
    near = np.abs(vals - 1.0)
    gap = float(near[near >= TOL].min())
    return int((near < TOL).sum()), gap


def main():
    print("object: design code 7, base 2, dimension 2; cells with i AND j = 0")
    print("graph: axis-nearest neighbours, no diagonals; levels 1 to %d" % TOP)
    print()
    print(
        "level nodes edges tree recursion nu_strip nu_hk nu_minus_corner "
        "n-2nu 3^(L-1) null_A mult_1 gap_A gap_1"
    )
    bad = []
    below = None
    for level in range(1, TOP + 1):
        rows, cols, pairs = build(level)
        count = len(rows)
        tree = is_tree(count, pairs)
        rec = level == 1 or recursion(level, rows, cols, pairs, below)
        below = lattice(rows, cols, pairs)
        parity = ((rows + cols) % 2).tolist()
        pl = [(int(a), int(b)) for a, b in pairs]
        nu_s = stripped(count, pl)
        nu_h = hopcroft_karp(count, pl, parity)
        origin = int(np.flatnonzero((rows == 0) & (cols == 0))[0])
        nu_c = hopcroft_karp(count, pl, parity, gone=(origin,))
        nu_c2 = stripped(count, pl, gone=(origin,))
        law = 3 ** (level - 1)
        null = unit = "-"
        gap_a = gap_n = "-"
        if level <= DENSE:
            null, gap_a, mat = nullity_adjacency(count, pairs)
            unit, gap_n = multiplicity_one(mat)
            gap_a = "%.4f" % gap_a
            gap_n = "%.4f" % gap_n
        print(
            "%d %d %d %s %s %d %d %d %d %d %s %s %s %s"
            % (
                level,
                count,
                len(pairs),
                tree,
                rec,
                nu_s,
                nu_h,
                nu_c,
                count - 2 * nu_h,
                law,
                null,
                unit,
                gap_a,
                gap_n,
            )
        )
        sys.stdout.flush()
        if not (tree and rec):
            bad.append("level %d: structure" % level)
        if nu_s != nu_h or nu_c != nu_c2:
            bad.append("level %d: two matching routines disagree" % level)
        if nu_c != nu_h - 1:
            bad.append("level %d: the corner is not essential" % level)
        if count - 2 * nu_h != law:
            bad.append(
                "level %d: n - 2 nu = %d, not %d" % (level, count - 2 * nu_h, law)
            )
        if level <= DENSE and (null != law or unit != law):
            bad.append("level %d: dense nullity %s, mult_1 %s" % (level, null, unit))
    print()
    print("matching polynomial against the characteristic polynomial, exact integers")
    for level in range(1, EXACT + 1):
        rows, cols, pairs = build(level)
        pl = [(int(a), int(b)) for a, b in pairs]
        same, poly = exact_check(len(rows), pl)
        top = max(k for k, m in enumerate(poly) if m)
        print(
            "  level %d: det(xI - A) = sum (-1)^k m_k x^(n-2k): %s; max matching size %d; m_nu = %d"
            % (level, same, top, poly[top])
        )
        if not same:
            bad.append("level %d: Sachs expansion fails" % level)
    print()
    if bad:
        print("FAILED")
        for line in bad:
            print("  " + line)
        raise SystemExit(1)
    print(
        "all levels: tree, n - 2 nu = 3^(L-1), corner essential, nullity and mult_1 agree"
    )


main()
