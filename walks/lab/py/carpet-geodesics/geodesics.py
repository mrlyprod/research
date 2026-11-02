import heapq
import math
import sys
import time

import numpy as np

R2 = math.sqrt(2.0)
R5 = math.sqrt(5.0)
LEVEL1 = 2 * R5 - 3 * R2

# DESIGN


def code(n):
    return sum(
        1 << (n * i + j) for i in range(n) for j in range(n) if not (i % 2 and j % 2)
    )


def voids(n, level):
    out = []

    def walk(x0, y0, size, depth):
        if depth == level:
            return
        c = size // n
        for i in range(n):
            for j in range(n):
                if i % 2 and j % 2:
                    out.append((x0 + i * c, y0 + j * c, c))
                else:
                    walk(x0 + i * c, y0 + j * c, c, depth + 1)

    walk(0, 0, n**level, 0)
    return np.array(out, dtype=float).reshape(-1, 3)


# VISIBILITY


def blocked(p, Q, H):
    if len(H) == 0:
        return np.zeros(len(Q), bool)
    out = np.zeros(len(Q), bool)
    dx = Q[:, 0] - p[0]
    dy = Q[:, 1] - p[1]
    bx0 = np.minimum(p[0], Q[:, 0])
    bx1 = np.maximum(p[0], Q[:, 0])
    by0 = np.minimum(p[1], Q[:, 1])
    by1 = np.maximum(p[1], Q[:, 1])
    step = max(1, 4_000_000 // max(1, len(Q)))
    for s in range(0, len(H), step):
        a0 = H[s : s + step, 0][None, :]
        b0 = H[s : s + step, 1][None, :]
        a1 = a0 + H[s : s + step, 2][None, :]
        b1 = b0 + H[s : s + step, 2][None, :]
        near = (
            (a1 > bx0[:, None])
            & (a0 < bx1[:, None])
            & (b1 > by0[:, None])
            & (b0 < by1[:, None])
        )
        if not near.any():
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            t1 = (a0 - p[0]) / dx[:, None]
            t2 = (a1 - p[0]) / dx[:, None]
            u1 = (b0 - p[1]) / dy[:, None]
            u2 = (b1 - p[1]) / dy[:, None]
        zx = dx[:, None] == 0
        zy = dy[:, None] == 0
        inx = (a0 < p[0]) & (p[0] < a1)
        iny = (b0 < p[1]) & (p[1] < b1)
        lox = np.where(zx, np.where(inx, -np.inf, np.inf), np.minimum(t1, t2))
        hix = np.where(zx, np.where(inx, np.inf, -np.inf), np.maximum(t1, t2))
        loy = np.where(zy, np.where(iny, -np.inf, np.inf), np.minimum(u1, u2))
        hiy = np.where(zy, np.where(iny, np.inf, -np.inf), np.maximum(u1, u2))
        lo = np.maximum(np.maximum(lox, loy), 0.0)
        hi = np.minimum(np.minimum(hix, hiy), 1.0)
        out |= (near & (hi - lo > 1e-12)).any(1)
    return out


def shortest(H, A, B, bound):
    A = np.array(A, float)
    B = np.array(B, float)
    c = H[:, :2] + H[:, 2:3] / 2
    r = H[:, 2] * R2 / 2
    keep = np.hypot(*(c - A).T) + np.hypot(*(c - B).T) - 2 * r <= bound
    H = H[keep]
    C = np.concatenate(
        [
            H[:, :2],
            H[:, :2] + H[:, 2:3] * [1, 0],
            H[:, :2] + H[:, 2:3] * [0, 1],
            H[:, :2] + H[:, 2:3],
        ]
    )
    C = np.unique(C, axis=0)
    C = C[np.hypot(*(C - A).T) + np.hypot(*(C - B).T) <= bound + 1e-9]
    P = np.concatenate([[A, B], C])
    n = len(P)
    tail = np.hypot(*(P - B).T)
    dist = np.full(n, np.inf)
    dist[0] = 0.0
    done = np.zeros(n, bool)
    heap = [(tail[0], 0)]
    while heap:
        _f, u = heapq.heappop(heap)
        if done[u]:
            continue
        done[u] = True
        if u == 1:
            return dist[1], len(H), n
        cand = np.where(~done)[0]
        step = np.hypot(*(P[cand] - P[u]).T)
        cand = cand[dist[u] + step + tail[cand] <= bound + 1e-9]
        if len(cand) == 0:
            continue
        ok = ~blocked(P[u], P[cand], H)
        v = cand[ok]
        nd = dist[u] + np.hypot(*(P[v] - P[u]).T)
        for w, d in zip(v[nd < dist[v] - 1e-13], nd[nd < dist[v] - 1e-13]):
            dist[w] = d
            heapq.heappush(heap, (d + tail[w], w))
    return np.inf, len(H), n


def corner(n, level):
    S = n**level
    H = voids(n, level)
    slack = 0.05
    while True:
        bound = (R2 + slack) * S
        d, _nh, nv = shortest(H, (0, 0), (S, S), bound)
        if d < np.inf:
            return d / S, len(H), nv
        slack *= 2


def envelope(n, level):
    return 2 - (2 - R2) * (1 - 1 / n) ** level


# VERBS


def verb_code():
    print("code of the side-N tile, bit N*i + j, void iff both digits odd")
    for n in (3, 5, 7, 9, 11):
        print(f"  N = {n:<2}  code {code(n)}  fill {n * n - ((n - 1) // 2) ** 2}")
    assert code(3) == 495


def verb_corner():
    print("corner distance D(N, L) = d((0,0), (1,1)), exact visibility graph")
    cases = [(3, 1), (3, 2), (3, 3), (3, 4), (5, 1), (5, 2), (7, 1), (7, 2), (9, 2)]
    if len(sys.argv) > 2:
        cases = [tuple(map(int, a.split(","))) for a in sys.argv[2:]]
    last = {}
    for n, level in cases:
        t = time.time()
        d, nh, nv = corner(n, level)
        env = envelope(n, level)
        assert R2 <= d <= env + 1e-12
        assert d >= last.get(n, 0.0) - 1e-12
        last[n] = d
        if n == 3:
            assert abs(d - 2 * R5 / 3) < 1e-12
        print(
            f"  N = {n:<2} L = {level}  D = {d:.12f}  N(D - sqrt2) = {n * (d - R2):.6f}  2 - D = {2 - d:.6f}  envelope {env:.6f}  voids {nh}  vertices {nv}  {time.time() - t:.1f} s"
        )


def verb_level1():
    print("level 1: D(N, 1) against sqrt2 + (2 sqrt5 - 3 sqrt2)/N")
    worst = 0.0
    for n in range(3, 42, 2):
        d, _, _ = corner(n, 1)
        f = R2 + LEVEL1 / n
        worst = max(worst, abs(d - f))
        print(f"  N = {n:<2}  D = {d:.12f}  formula {f:.12f}  diff {d - f:+.1e}")
    print(f"  largest difference {worst:.1e}; 2 sqrt5 - 3 sqrt2 = {LEVEL1:.9f}")


# MAP


def periodic_edges(K=3):
    corners = np.array([(1.0, 1.0), (2.0, 1.0), (1.0, 2.0), (2.0, 2.0)])
    H = np.array(
        [
            (1.0 + 2 * i, 1.0 + 2 * j, 1.0)
            for i in range(-K - 2, K + 3)
            for j in range(-K - 2, K + 3)
        ]
    )
    ea, eb, ev = [], [], []
    for a, pa in enumerate(corners):
        Q, idx = [], []
        for b, pb in enumerate(corners):
            for i in range(-K, K + 1):
                for j in range(-K, K + 1):
                    q = pb + 2 * np.array([i, j])
                    if np.any(q != pa):
                        Q.append(q)
                        idx.append(b)
        Q = np.array(Q)
        ok = ~blocked(pa, Q, H)
        for q, b in zip(Q[ok], np.array(idx)[ok]):
            ea.append(a)
            eb.append(b)
            ev.append(q - pa)
    return np.array(ea), np.array(eb), np.array(ev)


def gauge(h, P, V):
    return np.max((V @ P.T) / h[None, :], axis=1)


def homogenise(h, P, E):
    ea, eb, ev = E
    w = gauge(h, P, ev)
    g0 = P @ ev.T
    lo = np.zeros(len(P))
    hi = np.full(len(P), 2.0)
    flat = ea * 4 + eb
    for _ in range(48):
        lam = (lo + hi) / 2
        g = g0 - lam[:, None] * w[None, :]
        G = np.full((len(P), 16), -np.inf)
        for k in range(16):
            sel = flat == k
            if sel.any():
                G[:, k] = g[:, sel].max(1)
        G = G.reshape(len(P), 4, 4)
        for k in range(4):
            G = np.maximum(G, G[:, :, k : k + 1] + G[:, k : k + 1, :])
        pos = (np.diagonal(G, axis1=1, axis2=2) > 1e-12).any(1)
        lo = np.where(pos, lam, lo)
        hi = np.where(pos, hi, lam)
    return (lo + hi) / 2


def octagon(V):
    a = np.abs(V)
    big = np.maximum(a[:, 0], a[:, 1])
    small = np.minimum(a[:, 0], a[:, 1])
    return big - small + R2 * small


def node_cycles():
    from itertools import combinations, permutations

    out = [[(a, a)] for a in range(4)]
    for k in (2, 3, 4):
        for sub in combinations(range(4), k):
            for rest in permutations(sub[1:]):
                seq = (sub[0],) + rest
                if k == 2 and rest[0] < sub[0]:
                    continue
                out.append([(seq[i], seq[(i + 1) % k]) for i in range(k)])
    return out


def hull_gauge(W, V):
    from scipy.spatial import ConvexHull

    eq = ConvexHull(W).equations
    return np.max((V @ eq[:, :2].T) / (-eq[:, 2])[None, :], axis=1)


def cycle_points(cost, P, E, cycles):
    ea, eb, ev = E
    g0 = P @ ev.T
    pair = ea * 4 + eb
    groups = [np.where(pair == k)[0] for k in range(16)]
    lo = np.zeros(len(P))
    hi = np.full(len(P), 2.0)

    def best(lam):
        g = g0 - lam[:, None] * cost[None, :]
        G = np.full((len(P), 16), -np.inf)
        A = np.zeros((len(P), 16), int)
        for k, idx in enumerate(groups):
            if len(idx):
                j = g[:, idx].argmax(1)
                G[:, k] = g[np.arange(len(P)), idx[j]]
                A[:, k] = idx[j]
        S = np.stack([sum(G[:, a * 4 + b] for a, b in c) for c in cycles], 1)
        return S, A

    for _ in range(48):
        lam = (lo + hi) / 2
        S, _ = best(lam)
        pos = S.max(1) > 1e-12
        lo = np.where(pos, lam, lo)
        hi = np.where(pos, hi, lam)
    S, A = best(lo)
    pick = S.argmax(1)
    W = np.zeros((len(P), 2))
    for t in range(len(P)):
        idx = [A[t, a * 4 + b] for a, b in cycles[pick[t]]]
        W[t] = ev[idx].sum(0) / cost[idx].sum()
    return W


def verb_map():
    ndir = 720
    t = time.time()
    th = np.linspace(0, 2 * np.pi, ndir, endpoint=False)
    P = np.stack([np.cos(th), np.sin(th)], 1)
    U = np.stack([np.cos(th + np.pi / ndir), np.sin(th + np.pi / ndir)], 1)
    probe = np.array(
        [
            [1, 0],
            [math.cos(math.pi / 8), math.sin(math.pi / 8)],
            [2 / R5, 1 / R5],
            [1 / R2, 1 / R2],
        ]
    )
    names = ["0", "22.5", "atan(1/2)", "45"]
    cycles = node_cycles()
    print(
        f"homogenisation map: upper bounds on nu_L from cycle points, {ndir} support directions, {len(cycles)} node cycles"
    )
    for K in (3, 6):
        E = periodic_edges(K)
        cost = np.hypot(*E[2].T)
        print(f"  window {K} periods, {len(E[0])} edges")
        print(
            "  L    "
            + "  ".join(f"{s:>10}" for s in names)
            + "   gap to octagon at least"
        )
        for level in range(1, 21):
            W = cycle_points(cost, P, E, cycles)
            cost = hull_gauge(W, E[2])
            nu = hull_gauge(W, probe)
            gap = np.max(
                1
                - hull_gauge(W, np.concatenate([P, U]))
                / octagon(np.concatenate([P, U]))
            )
            if level <= 4 or level % 5 == 0:
                print(
                    f"  {level:<4} "
                    + "  ".join(f"{x:10.6f}" for x in nu)
                    + f"   {gap:.2e}"
                )
    print("  oct  " + "  ".join(f"{x:10.6f}" for x in octagon(probe)))
    two = np.array([[4.0, 2.0]])
    print(
        f"  explicit path (0,0)-(2,1)-(3,2)-(4,2): {R5 + R2 + 1:.6f} against octagon {octagon(two)[0]:.6f} for displacement (4,2)"
    )
    print(f"  {time.time() - t:.1f} s")


def verb_bridge():
    E = periodic_edges()
    ndir = 1440
    th = np.linspace(0, 2 * np.pi, ndir, endpoint=False)
    P = np.stack([np.cos(th), np.sin(th)], 1)
    h = homogenise(np.ones(ndir), P, E)
    print(
        "bridge at level 1: d_N((0,0), (1, (N-1)/(2N))) against the level-1 norm of the same vector"
    )
    for n in (11, 21, 31, 41):
        H = voids(n, 1)
        B = (n, (n - 1) // 2)
        nu = gauge(h, P, np.array([[1.0, (n - 1) / (2 * n)]]))[0]
        d, _, _ = shortest(H, (0, 0), B, 1.2 * n)
        print(
            f"  N = {n:<2}  d = {d / n:.6f}  level-1 norm {nu:.6f}  euclid {math.hypot(1, (n - 1) / (2 * n)):.6f}  N(d - norm) = {n * (d / n - nu):.4f}"
        )


def verb_hull():
    from itertools import combinations, product

    from scipy.spatial import ConvexHull

    V = []
    for i in range(3):
        for s in (1, -1):
            e = np.zeros(3)
            e[i] = s
            V.append(e)
    for i, j in combinations(range(3), 2):
        for s, t in product((1, -1), repeat=2):
            e = np.zeros(3)
            e[i], e[j] = s, t
            V.append(e / R2)
    hull = ConvexHull(np.array(V))
    faces = len(np.unique(np.round(hull.equations, 9), axis=0))
    edges = len({tuple(sorted(p)) for s in hull.simplices for p in combinations(s, 2)})
    print(
        f"dim 3 ball: hull of the 18 free unit vectors has {len(hull.vertices)} vertices, {edges} edges, {faces} faces, all triangles: {faces == len(hull.simplices)}"
    )
    assert (len(hull.vertices), edges, faces) == (18, 48, 32)


VERBS = {
    "code": verb_code,
    "corner": verb_corner,
    "level1": verb_level1,
    "map": verb_map,
    "bridge": verb_bridge,
    "hull": verb_hull,
}

if __name__ == "__main__":
    names = sys.argv[1:2] or list(VERBS)
    for name in names:
        t = time.time()
        VERBS[name]()
        print(f"[{name}] {time.time() - t:.1f} s", flush=True)
