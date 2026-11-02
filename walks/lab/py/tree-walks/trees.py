import math
import sys
from fractions import Fraction
from heapq import heapify, heappop, heappush

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, csc_matrix, identity
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.sparse.linalg import spsolve

SUBJECTS = {"79": (3, 79, 6), "7": (2, 7, 9)}
MASK = (1 << 64) - 1
OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
PORTS = ("E", "W", "N", "S")
OPPOSITE = {"E": "W", "W": "E", "N": "S", "S": "N"}

# RENDER


def tile_cells(code, base):
    return [
        (a, c) for a in range(base) for c in range(base) if code >> (base * a + c) & 1
    ]


def render(code, base, level):
    tile = tile_cells(code, base)
    a = np.array([t[0] for t in tile], dtype=np.int64)
    c = np.array([t[1] for t in tile], dtype=np.int64)
    x0 = np.zeros(1, dtype=np.int64)
    x1 = np.zeros(1, dtype=np.int64)
    for step in range(level):
        scale = base**step
        x0 = (a[:, None] * scale + x0[None, :]).ravel()
        x1 = (c[:, None] * scale + x1[None, :]).ravel()
    side = base**level
    return np.sort(x0 * side + x1), side


def grid_edges(flat, side):
    index = np.full(side * side, -1, dtype=np.int64)
    index[flat] = np.arange(len(flat))
    right = (flat % side != side - 1) & (
        index[np.minimum(flat + 1, side * side - 1)] >= 0
    )
    down = flat + side < side * side
    down[down] = index[flat[down] + side] >= 0
    u = np.concatenate([np.arange(len(flat))[right], np.arange(len(flat))[down]])
    v = np.concatenate([index[flat[right] + 1], index[flat[down] + side]])
    order = np.lexsort((flat[v], flat[u]))
    return u[order], v[order]


def adjacency(n, u, v):
    data = np.ones(len(u), dtype=np.int8)
    half = coo_matrix((data, (u, v)), shape=(n, n))
    return (half + half.T).tocsr()


def fold(h, word):
    return ((h ^ word) * PRIME) & MASK


def digest(flat, u, v):
    cells = OFFSET
    for f in flat.tolist():
        cells = fold(cells, f)
    edges = OFFSET
    for a, b in zip(flat[u].tolist(), flat[v].tolist()):
        edges = fold(fold(edges, a), b)
    return cells, edges


def is_tree(n, u, v):
    return (
        len(u) == n - 1
        and connected_components(adjacency(n, u, v), directed=False)[0] == 1
    )


def corners(side, flat):
    index = {f: i for i, f in enumerate(flat.tolist())}
    return index[0], index[side * (side - 1)], index[side - 1]


# EXACT RESISTANCE


def schur_resistance(n, u, v, s, t):
    nbr = [{} for _ in range(n)]
    one = Fraction(1)
    for a, b in zip(u.tolist(), v.tolist()):
        nbr[a][b] = nbr[a].get(b, 0) + one
        nbr[b][a] = nbr[b].get(a, 0) + one
    alive = [True] * n
    heap = [(len(nbr[i]), i) for i in range(n) if i != s and i != t]
    heapify(heap)
    while heap:
        degree, node = heappop(heap)
        if not alive[node]:
            continue
        if degree != len(nbr[node]):
            heappush(heap, (len(nbr[node]), node))
            continue
        total = sum(nbr[node].values())
        items = list(nbr[node].items())
        for a, _ in items:
            del nbr[a][node]
        for i, (a, ca) in enumerate(items):
            for b, cb in items[i + 1 :]:
                w = ca * cb / total
                nbr[a][b] = nbr[a].get(b, 0) + w
                nbr[b][a] = nbr[b].get(a, 0) + w
        alive[node] = False
        nbr[node] = {}
        for a, _ in items:
            if a != s and a != t:
                heappush(heap, (len(nbr[a]), a))
    return 1 / nbr[s][t]


def float_resistance(n, u, v, s, t):
    adj = adjacency(n, u, v).astype(np.float64)
    degree = np.asarray(adj.sum(axis=1)).ravel()
    lap = (
        coo_matrix((degree, (np.arange(n), np.arange(n))), shape=(n, n)) - adj
    ).tocsc()
    keep = np.array([i for i in range(n) if i != t])
    rhs = np.zeros(n)
    rhs[s] = 1.0
    sol = spsolve(lap[keep][:, keep].tocsc(), rhs[keep])
    return float(sol[np.searchsorted(keep, s)])


def solve_exact(matrix, rhs):
    n = len(rhs)
    rows = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next(r for r in range(col, n) if rows[r][col] != 0)
        rows[col], rows[pivot] = rows[pivot], rows[col]
        head = rows[col][col]
        rows[col] = [x / head for x in rows[col]]
        for r in range(n):
            if r != col and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [x - factor * y for x, y in zip(rows[r], rows[col])]
    return [rows[i][n] for i in range(n)]


def hitting_time(n, u, v, source, target):
    nbrs = [[] for _ in range(n)]
    for a, b in zip(u.tolist(), v.tolist()):
        nbrs[a].append(b)
        nbrs[b].append(a)
    keep = [i for i in range(n) if i != target]
    where = {node: k for k, node in enumerate(keep)}
    matrix = [[Fraction(0)] * len(keep) for _ in keep]
    rhs = []
    for node in keep:
        row = where[node]
        matrix[row][row] = Fraction(len(nbrs[node]))
        for other in nbrs[node]:
            if other != target:
                matrix[row][where[other]] -= 1
        rhs.append(Fraction(len(nbrs[node])))
    return solve_exact(matrix, rhs)[where[source]]


def show(x):
    return str(x.numerator) if x.denominator == 1 else str(x)


def distances(adj, sources):
    return dijkstra(adj, directed=False, unweighted=True, indices=sources).astype(
        np.int64
    )


# VERBS


def verb_render():
    for base, code, top in SUBJECTS.values():
        for level in range(1, top + 1):
            flat, side = render(code, base, level)
            u, v = grid_edges(flat, side)
            cells, edges = digest(flat, u, v)
            print(
                f"base {base} code {code} level {level} nodes {len(flat)} edges {len(u)} "
                f"cells {cells:016x} edge-list {edges:016x}"
            )
            assert is_tree(len(flat), u, v)


def verb_resist():
    for base, code, top in SUBJECTS.values():
        print(f"base {base} code {code}")
        for level in range(1, top + 1):
            flat, side = render(code, base, level)
            u, v = grid_edges(flat, side)
            n = len(flat)
            o, a, b = corners(side, flat)
            adj = adjacency(n, u, v)
            far = schur_resistance(n, u, v, a, b)
            near = [schur_resistance(n, u, v, o, a), schur_resistance(n, u, v, o, b)]
            dist = distances(adj, [a, o])
            want_far, want_near = 2 * base**level - 2, base**level - 1
            assert far == dist[0][b] == want_far and near == [want_near] * 2
            assert dist[1][a] == dist[1][b] == want_near
            drift = abs(float_resistance(n, u, v, a, b) - want_far)
            twist = ""
            if level >= 2:
                half = base ** (level - 1)
                where = {f: i for i, f in enumerate(flat.tolist())}
                p, q = where[(half - 1) * side + half], where[half * side + half - 1]
                gap = int(distances(adj, [p])[0][q])
                assert gap == 4 * half - 2
                twist = f" diagonal neighbours at graph distance {gap} = 4 b^(n-1) - 2"
            line = (
                f"level {level} nodes {n} R(A,B) {show(far)} d(A,B) {dist[0][b]} formula {want_far} "
                f"R(O,A) {show(near[0])} R(O,B) {show(near[1])} formula {want_near} float gap {drift:.1e}{twist}"
            )
            if n <= 130:
                hab = hitting_time(n, u, v, a, b)
                hba = hitting_time(n, u, v, b, a)
                assert hab + hba == 2 * (n - 1) * far
                line += f" commute {show(hab + hba)} = 2 (nodes - 1) R"
            print(line)


def verb_volume():
    plan = {"79": 5, "7": 8}
    for key, (base, code, top) in SUBJECTS.items():
        level = plan[key]
        flat, side = render(code, base, level + 1)
        u, v = grid_edges(flat, side)
        n = len(flat)
        fill = len(tile_cells(code, base))
        alpha = math.log(fill) / math.log(base)
        adj = adjacency(n, u, v)
        x0, x1 = flat // side, flat % side
        inner = base**level
        centres = np.where((x0 < inner) & (x1 < inner))[0]
        reach = 2 * (inner - 1)
        counts = np.zeros((len(centres), reach + 1), dtype=np.int32)
        blocks = [0] * (level + 1)
        ratio_origin = 0.0
        origin = int(np.where(flat == 0)[0][0])
        for start in range(0, len(centres), 256):
            chunk = centres[start : start + 256]
            dist = distances(adj, chunk.tolist())
            for row, centre in enumerate(chunk.tolist()):
                counts[start + row] = np.cumsum(
                    np.bincount(dist[row], minlength=reach + 1)[: reach + 1]
                )
                for k in range(1, level + 1):
                    near = dist[row] <= base**k - 2
                    ids = np.unique((x0[near] // base**k) * side + x1[near] // base**k)
                    blocks[k] = max(blocks[k], len(ids))
                if centre == origin:
                    far = np.maximum(x0, x1)
                    ratio_origin = float((dist[row][far > 0] / far[far > 0]).max())
        print(
            f"base {base} code {code}: graph of level {level + 1} ({n} nodes), "
            f"{len(centres)} centres in its corner copy of level {level}, alpha {alpha:.6f}"
        )
        for k in range(level + 1):
            low = counts[:, 2 * (base**k - 1)].min()
            cap = max(base**k - 2, 0)
            high = counts[:, cap].max()
            assert low >= fill**k and high <= 7 * fill**k and blocks[k] <= 7
            print(
                f"k {k}: min V(x, {2 * (base**k - 1)}) = {low} vs fill^k = {fill**k}; "
                f"max V(x, {cap}) = {high} vs 7 fill^k = {7 * fill**k}; "
                f"level-{k} blocks met by the ball of radius {cap}: at most {blocks[k]}"
            )
        radii = np.arange(1, inner - 1)
        for k in range(level):
            window = radii[(radii >= base**k) & (radii < base ** (k + 1))]
            ratios = counts[:, window] / window[None, :] ** alpha
            print(
                f"r in [{base**k}, {base ** (k + 1)}): V/r^alpha min {ratios.min():.4f} max {ratios.max():.4f}"
            )
        print(f"max over cells y of d(O, y)/max(y0, y1): {ratio_origin:.4f}")


def low_spectrum(code, base, level, modes=4):
    flat, side = render(code, base, level)
    u, v = grid_edges(flat, side)
    n = len(flat)
    adj = adjacency(n, u, v).toarray().astype(np.float64)
    lap = np.diag(adj.sum(axis=1)) - adj
    return eigh(lap, subset_by_index=[0, modes])[0]


def exit_time(adj, dist, radius):
    inside = np.where(dist < radius)[0]
    degree = np.asarray(adj.sum(axis=1)).ravel().astype(np.float64)
    sub = adj[inside][:, inside].astype(np.float64)
    system = (
        identity(len(inside)) - csc_matrix(sub.multiply(1.0 / degree[inside][:, None]))
    ).tocsc()
    return float(
        spsolve(system, np.ones(len(inside)))[
            np.searchsorted(inside, np.where(dist == 0)[0][0])
        ]
    )


def verb_scale():
    for key, (base, code, top) in SUBJECTS.items():
        fill = len(tile_cells(code, base))
        exact = math.log(fill * base) / math.log(base)
        print(f"base {base} code {code}: log({fill * base})/log({base}) = {exact:.6f}")
        pairs = [(3, 4), (4, 5)] if key == "79" else [(5, 6), (6, 7)]
        for lo, hi in pairs:
            a, b = low_spectrum(code, base, lo), low_spectrum(code, base, hi)
            est = [math.log(a[k] / b[k]) / math.log(base) for k in range(1, 5)]
            print(
                f"spectral levels {lo} to {hi}, modes 1..4: "
                + " ".join(f"{x:.4f}" for x in est)
            )
        level = 6
        flat, side = render(code, base, level if key == "79" else 9)
        u, v = grid_edges(flat, side)
        n = len(flat)
        adj = adjacency(n, u, v)
        dist = distances(adj, [int(np.where(flat == 0)[0][0])])[0]
        last = None
        for k in range(1, (5 if key == "79" else 8) + 1):
            r = base**k
            t = exit_time(adj, dist, r)
            ratio = (
                ""
                if last is None
                else f" ratio {t / last:.4f} exponent {math.log(t / last) / math.log(base):.4f}"
            )
            print(f"mean exit time from the corner at radius {r}: {t:.4f}{ratio}")
            last = t


def exit_time_exact(adj, dist, radius):
    indptr, indices = adj.indptr, adj.indices
    inside = np.where(dist < radius)[0]
    order = inside[np.argsort(-dist[inside], kind="stable")]
    slope, shift = {}, {}
    for node in order.tolist():
        nbrs = indices[indptr[node] : indptr[node + 1]].tolist()
        kids = [c for c in nbrs if dist[c] == dist[node] + 1 and dist[c] < radius]
        degree = len(nbrs)
        sum_shift = sum((shift[c] for c in kids), Fraction(0))
        sum_slope = sum((slope[c] for c in kids), Fraction(0))
        denom = degree - sum_slope
        if dist[node] == 0:
            return (degree + sum_shift) / denom
        shift[node] = (degree + sum_shift) / denom
        slope[node] = 1 / denom
    raise ValueError("no root")


def verb_exact():
    for key, (base, code, top) in SUBJECTS.items():
        fill = len(tile_cells(code, base))
        flat, side = render(code, base, top)
        u, v = grid_edges(flat, side)
        assert is_tree(len(flat), u, v)
        adj = adjacency(len(flat), u, v)
        dist = distances(adj, [int(np.where(flat == 0)[0][0])])[0]
        print(f"base {base} code {code}")
        for k in range(1, (5 if key == "79" else 8) + 1):
            t = exit_time_exact(adj, dist, base**k)
            want = Fraction(base**k * (fill**k + 1), 2)
            assert t == want
            print(
                f"k {k}: E_O tau at radius {base**k} = {show(t)} = {base}^k ({fill}^k + 1)/2"
            )


def symmetries(cells, base):
    last = base - 1
    maps = [
        lambda a, c: (a, c),
        lambda a, c: (last - a, c),
        lambda a, c: (a, last - c),
        lambda a, c: (last - a, last - c),
        lambda a, c: (c, a),
        lambda a, c: (last - c, a),
        lambda a, c: (c, last - a),
        lambda a, c: (last - c, last - a),
    ]
    return [frozenset(f(a, c) for a, c in cells) for f in maps]


def code_of(cells, base):
    return sum(1 << (base * a + c) for a, c in cells)


def tile_tree(cells, base):
    cells = sorted(cells)
    have = set(cells)
    steps = [
        (a, c, axis, 1)
        for a, c in cells
        for axis in (0, 1)
        if ((a + 1, c) if axis == 0 else (a, c + 1)) in have
    ]
    if len(cells) < 2 or len(steps) != len(cells) - 1:
        return None
    seen, stack = {cells[0]}, [cells[0]]
    while stack:
        a, c = stack.pop()
        for da, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nxt = (a + da, c + dc)
            if nxt in have and nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    if len(seen) != len(cells):
        return None
    last = base - 1
    ports = {}
    axes = {axis for _, _, axis, _ in steps}
    if 0 in axes:
        both = {c for a, c in cells if a == last} & {c for a, c in cells if a == 0}
        if len(both) != 1:
            return None
        row = both.pop()
        ports["P+"], ports["P-"] = (last, row), (0, row)
    if 1 in axes:
        both = {a for a, c in cells if c == last} & {a for a, c in cells if c == 0}
        if len(both) != 1:
            return None
        col = both.pop()
        ports["Q+"], ports["Q-"] = (col, last), (col, 0)
    return ports


def port_name(axis, sign):
    return ("P" if axis == 0 else "Q") + ("+" if sign > 0 else "-")


def opposite(name):
    return name[0] + ("-" if name[1] == "+" else "+")


def tile_path(cells, start, goal):
    have = set(cells)
    back, queue = {start: None}, [start]
    for node in queue:
        for axis, sign in ((0, 1), (0, -1), (1, 1), (1, -1)):
            nxt = (node[0] + sign * (axis == 0), node[1] + sign * (axis == 1))
            if nxt in have and nxt not in back:
                back[nxt] = (node, axis, sign)
                queue.append(nxt)
    path, node = [], goal
    while back[node] is not None:
        prev, axis, sign = back[node]
        path.append((axis, sign))
        node = prev
    return path[::-1]


def port_system(cells, ports, base):
    names = sorted(ports)
    pairs = [(x, y) for i, x in enumerate(names) for y in names[i + 1 :]]
    level_one = {}
    for x, y in pairs:
        level_one[(x, y)] = len(tile_path(cells, ports[x], ports[y]))
    active = [p for p in pairs if level_one[p] > 0]
    cell_of = ports
    rows, constants = {}, {}
    for x, y in pairs:
        steps = tile_path(cells, cell_of[x], cell_of[y])
        enter, row = x, {}
        for axis, sign in steps:
            leave = port_name(axis, sign)
            if enter != leave:
                key = tuple(sorted((enter, leave)))
                row[key] = row.get(key, 0) + 1
            enter = opposite(leave)
        if enter != y:
            key = tuple(sorted((enter, y)))
            row[key] = row.get(key, 0) + 1
        rows[(x, y)] = row
        constants[(x, y)] = len(steps)
    return pairs, active, level_one, rows, constants


def components(matrix):
    n = len(matrix)
    reach = (matrix > 0) | np.eye(n, dtype=bool)
    for k in range(n):
        reach = reach | (reach[:, [k]] & reach[[k], :])
    labels, classes = [-1] * n, []
    for i in range(n):
        if labels[i] < 0:
            members = [j for j in range(n) if reach[i, j] and reach[j, i]]
            for j in members:
                labels[j] = len(classes)
            classes.append(members)
    return reach, labels, classes


def analyse(cells, ports, base):
    pairs, active, level_one, rows, constants = port_system(cells, ports, base)
    index = {p: i for i, p in enumerate(active)}
    matrix = np.zeros((len(active), len(active)), dtype=np.int64)
    for p in active:
        for key, count in rows[p].items():
            if key in index:
                matrix[index[p], index[key]] += count
    values = {p: level_one[p] for p in pairs}
    history = [dict(values)]
    for _ in range(40):
        values = {
            p: constants[p] + sum(count * values[key] for key, count in rows[p].items())
            for p in pairs
        }
        history.append(dict(values))
    if not active:
        return False, 0.0, history
    reach, labels, classes = components(matrix)
    radius = []
    for members in classes:
        block = matrix[np.ix_(members, members)].astype(float)
        radius.append(float(max(abs(np.linalg.eigvals(block)))))
    lam = max(radius)
    chain = {}

    def deepest(c):
        if c not in chain:
            below = [
                labels[j]
                for i in classes[c]
                for j in range(len(active))
                if matrix[i, j] and labels[j] != c
            ]
            chain[c] = (abs(radius[c] - lam) < 1e-9) + max(
                (deepest(d) for d in set(below)), default=0
            )
        return chain[c]

    uniform = lam > 1 + 1e-9
    for i in range(len(active)):
        reachable = [labels[j] for j in range(len(active)) if reach[i, j]]
        uniform &= (
            abs(max(radius[c] for c in reachable) - lam) < 1e-9
            and deepest(labels[i]) <= 1
        )
    return bool(uniform), lam, history


def verify_history(code, base, ports, history, levels):
    for level in range(1, levels + 1):
        flat, side = render(code, base, level)
        u, v = grid_edges(flat, side)
        n = len(flat)
        assert is_tree(n, u, v)
        index = {f: i for i, f in enumerate(flat.tolist())}
        adj = adjacency(n, u, v)
        names = sorted(ports)
        for i, x in enumerate(names):
            sx = port_cell(ports, x, base, level)
            dist = distances(adj, [index[sx[0] * side + sx[1]]])[0]
            for y in names[i + 1 :]:
                sy = port_cell(ports, y, base, level)
                assert (
                    dist[index[sy[0] * side + sy[1]]] == history[level - 1][(x, y)]
                ), (code, base, level, x, y)


def port_cell(ports, name, base, level):
    last = base**level - 1
    rep = (base**level - 1) // (base - 1)
    if name[0] == "P":
        row = ports[name][1] * rep
        return (last if name[1] == "+" else 0, row)
    col = ports[name][0] * rep
    return (col, last if name[1] == "+" else 0)


def brute(base, levels):
    total = base**2
    agree, trees = 0, 0
    for code in range(1, 1 << total):
        cells = frozenset((i // base, i % base) for i in range(total) if code >> i & 1)
        ok = len(cells) >= 2
        for level in range(1, levels + 1):
            if not ok:
                break
            flat, side = render(code, base, level)
            u, v = grid_edges(flat, side)
            ok = is_tree(len(flat), u, v)
        passes = len(cells) >= 2 and tile_tree(cells, base) is not None
        assert ok == passes, (base, code)
        agree += 1
        trees += ok
    return agree, trees


def verb_siblings():
    for base, levels in ((2, 4), (3, 3), (4, 2)):
        agree, trees = brute(base, levels)
        print(
            f"base {base}: {agree} codes, {trees} trees at levels 1..{levels} (tile test agrees on every code)"
        )
    for base in (2, 3, 4):
        total = base**2
        seen, classes, rows = 0, {}, []
        for code in range(1, 1 << total):
            cells = frozenset(
                (i // base, i % base) for i in range(total) if code >> i & 1
            )
            ports = tile_tree(cells, base)
            if ports is None:
                continue
            seen += 1
            canon = min(code_of(f, base) for f in symmetries(cells, base))
            if canon not in classes:
                classes[canon] = (cells, ports)
            rows.append(code)
        print(
            f"base {base}: {seen} codes pass the tile test, {len(classes)} classes under the 8 symmetries of the square"
        )
        uniform, mixed = [], []
        for canon in sorted(classes):
            cells, ports = classes[canon]
            uniform_flag, lam, history = analyse(cells, ports, base)
            fill = len(cells)
            verify_history(canon, base, ports, history, 3)
            (uniform if uniform_flag else mixed).append((canon, fill, lam, history))
        groups = {}
        for canon, fill, lam, history in uniform:
            groups.setdefault((fill, round(lam, 6)), []).append(canon)
        print(
            f"  uniform port growth: {len(uniform)} classes; not uniform: {len(mixed)}"
        )
        for (fill, lam), members in sorted(groups.items()):
            dw_graph = 1 + math.log(fill) / math.log(lam)
            dw_side = math.log(fill * lam) / math.log(base)
            shown = " ".join(map(str, members[:8])) + (
                " ..." if len(members) > 8 else ""
            )
            print(
                f"  fill {fill} lambda {lam:.6f} d_w graph {dw_graph:.6f} d_w by side {dw_side:.6f}: {len(members)} classes, codes {shown}"
            )
        for canon, fill, lam, history in mixed:
            print(f"  not uniform: code {canon} fill {fill}")


VERBS = {
    "render": verb_render,
    "resist": verb_resist,
    "volume": verb_volume,
    "scale": verb_scale,
    "exact": verb_exact,
    "siblings": verb_siblings,
}

if __name__ == "__main__":
    names = sys.argv[1:] or list(VERBS)
    for name in names:
        print(f"== {name}")
        VERBS[name]()
