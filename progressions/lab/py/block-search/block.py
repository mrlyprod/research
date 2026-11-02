import argparse
import json
import math
import os
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import coo_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()

# GRID


class Grid:
    def __init__(self, m, k):
        self.m, self.k, self.n = m, k, m**k
        idx = np.arange(self.n)
        self.coords = np.stack(
            [(idx // m ** (k - 1 - i)) % m for i in range(k)], axis=1
        )
        half = m // 2
        cls_key = (2 * self.coords) % m
        self.cls = np.unique(cls_key @ (m ** np.arange(k)[::-1]), return_inverse=True)[
            1
        ]
        self.ncls = self.cls.max() + 1
        x, z = np.triu_indices(self.n, 1)
        s = self.coords[x] + self.coords[z]
        if m % 2:
            inv2 = pow(2, -1, m)
            y = ((s * inv2) % m) @ (m ** np.arange(k)[::-1])
            X, Z, Y = x, z, y
        else:
            ok = np.all(s % 2 == 0, axis=1)
            x, z, s = x[ok], z[ok], s[ok] // 2
            X, Z, Y = [], [], []
            for e in range(2**k):
                shift = np.array([(e >> (k - 1 - i)) & 1 for i in range(k)]) * half
                y = ((s + shift) % m) @ (m ** np.arange(k)[::-1])
                X.append(x)
                Z.append(z)
                Y.append(y)
            X, Z, Y = np.concatenate(X), np.concatenate(Z), np.concatenate(Y)
        self.X, self.Z, self.Y = X, Z, Y
        d = self.coords[X] - self.coords[Z]
        self.D = np.sum(d * d, axis=1)

    def point(self, v):
        return int(np.dot(v, self.m ** np.arange(self.k)[::-1]))

    def rows(self, mask):
        sel = mask[self.X] & mask[self.Z] & mask[self.Y]
        return self.X[sel], self.Z[sel], self.Y[sel], self.D[sel]

    def core(self, mask, mode):
        rx, rz, ry, _ = self.rows(mask)
        alive = mask.copy()
        while alive.any():
            keep = alive[rx] & alive[rz]
            rx, rz, ry = rx[keep], rz[keep], ry[keep]
            cnt = np.bincount(ry, minlength=self.n)
            if mode == "point":
                rem = alive & (cnt == 0)
            else:
                bad = np.bincount(self.cls[alive & (cnt > 0)], minlength=self.ncls) > 0
                rem = alive & ~bad[self.cls]
            if not rem.any():
                break
            alive &= ~rem
        return alive

    def rounds(self, mask, mode):
        rx, rz, ry, _ = self.rows(mask)
        alive = mask.copy()
        level = np.full(self.n, -1)
        r = 0
        while alive.any():
            keep = alive[rx] & alive[rz]
            rx, rz, ry = rx[keep], rz[keep], ry[keep]
            cnt = np.bincount(ry, minlength=self.n)
            if mode == "point":
                rem = alive & (cnt == 0)
            else:
                bad = np.bincount(self.cls[alive & (cnt > 0)], minlength=self.ncls) > 0
                rem = alive & ~bad[self.cls]
            if not rem.any():
                return None
            level[rem] = r
            alive &= ~rem
            r += 1
        return level, r

    def peelable(self, mask, mode):
        if (
            mode == "point"
            and np.bincount(self.cls[mask], minlength=self.ncls).max() > 1
        ):
            return False
        return not self.core(mask, mode).any()

    def minimal_core(self, mask, mode):
        core = self.core(mask, mode).copy()
        for v in np.flatnonzero(core):
            core[v] = False
            if self.peelable(core, mode):
                core[v] = True
        return core

    def certificate(self, mask, mode):
        out = self.rounds(mask, mode)
        if out is None:
            return None
        level, r = out
        g = {int(v): 3 ** (r - 1 - int(level[v])) for v in np.flatnonzero(mask)}
        rx, rz, ry, _ = self.rows(mask)
        for x, z, y in zip(rx.tolist(), rz.tolist(), ry.tolist()):
            assert g[x] + g[z] - 2 * g[y] > 0
        if mode == "class":
            for c in range(self.ncls):
                vals = {g[int(v)] for v in np.flatnonzero(mask & (self.cls == c))}
                assert len(vals) <= 1
        return g

    def lp(self, mask, mode, margin=True):
        pts = np.flatnonzero(mask)
        loc = -np.ones(self.n, dtype=int)
        loc[pts] = np.arange(len(pts))
        rx, rz, ry, rd = self.rows(mask)
        nv = len(pts)
        rows, cols, vals, rhs = [], [], [], []
        for i, (x, z, y, d) in enumerate(zip(loc[rx], loc[rz], loc[ry], rd)):
            rows += [i, i, i]
            cols += [x, z, y]
            vals += [-1.0, -1.0, 2.0]
            rhs.append(-float(d) / self.m**2 if margin else -1.0)
        i = len(rhs)
        if mode == "class":
            for c in range(self.ncls):
                mem = pts[self.cls[pts] == c]
                for a in mem:
                    for b in mem:
                        if a != b:
                            rows += [i, i]
                            cols += [loc[a], loc[b]]
                            vals += [-1.0, 1.0]
                            rhs.append(0.0)
                            i += 1
        rows += [i + j for j in range(nv) for _ in range(1)]
        cols += list(range(nv))
        vals += [1.0] * nv
        rows_m = rows + list(range(i, i + nv))
        cols_m = cols + [nv] * nv
        vals_m = vals + [-1.0] * nv
        rhs += [0.0] * nv
        A = coo_matrix((vals_m, (rows_m, cols_m)), shape=(i + nv, nv + 1)).tocsr()
        c = np.zeros(nv + 1)
        c[nv] = 1.0
        res = linprog(
            c, A_ub=A, b_ub=np.array(rhs), bounds=[(0, None)] * (nv + 1), method="highs"
        )
        if res.status != 0:
            return None
        f = res.x[:nv]
        slack = (
            f[loc[rx]]
            + f[loc[rz]]
            - 2 * f[loc[ry]]
            - (rd / self.m**2 if margin else 1.0)
        )
        return (
            {int(v): float(fv) for v, fv in zip(pts, f)},
            float(res.x[nv]),
            float(slack.min()) if len(slack) else 0.0,
        )

    def farkas(self, mask, mode):
        pts = np.flatnonzero(mask)
        loc = -np.ones(self.n, dtype=int)
        loc[pts] = np.arange(len(pts))
        rx, rz, ry, rd = self.rows(mask)
        nv, nr = len(pts), len(rx)
        rows, cols, vals = [], [], []
        for i, (x, z, y, d) in enumerate(zip(loc[rx], loc[rz], loc[ry], rd)):
            rows += [i, i, i, i]
            cols += [x, z, y, nv]
            vals += [-1.0, -1.0, 2.0, float(d) / self.m**2]
        extra = []
        if mode == "class":
            i = nr
            for c in range(self.ncls):
                mem = pts[self.cls[pts] == c]
                for a in mem:
                    for b in mem:
                        if a != b:
                            rows += [i, i]
                            cols += [loc[a], loc[b]]
                            vals += [-1.0, 1.0]
                            extra.append((int(a), int(b)))
                            i += 1
        A = coo_matrix((vals, (rows, cols)), shape=(nr + len(extra), nv + 1)).tocsr()
        c = np.zeros(nv + 1)
        c[nv] = -1.0
        res = linprog(
            c,
            A_ub=A,
            b_ub=np.zeros(nr + len(extra)),
            bounds=[(0, 1)] * nv + [(None, None)],
            method="highs",
        )
        t = -res.fun
        lam = -res.ineqlin.marginals
        return t, lam[:nr], lam[nr:], (rx, rz, ry, rd), extra


# BLOCKS


def ehps(a, b, eps):
    one, half = Fraction(1), Fraction(1, 2)
    s = a + b
    if half <= a < one and Fraction(2, 3) < s <= Fraction(7, 6):
        return 1
    if half <= a < one and b < half and Fraction(7, 6) + eps <= s <= Fraction(17, 12):
        return 2
    if (
        a < half <= b < one
        and Fraction(7, 6) + eps <= s <= Fraction(17, 12)
        and 2 * a + b >= Fraction(3, 2) + eps
    ):
        return 3
    return 0


def ehps_f(a, b, eps):
    g = a * a if a < Fraction(1, 2) else (a - Fraction(1, 2)) ** 2
    return 24 / eps**2 * (a + b) ** 2 + 6 * g


def rasterise(grid, eps, shift):
    m = grid.m
    mask = np.zeros(grid.n, dtype=bool)
    which = np.zeros(grid.n, dtype=int)
    for v in range(grid.n):
        a = Fraction(int(grid.coords[v, 0]), m) + shift[0]
        b = Fraction(int(grid.coords[v, 1]), m) + shift[1]
        which[v] = ehps(a, b, eps)
    mask = which > 0
    return mask, which


def bitmap(grid, mask, which=None):
    m = grid.m
    lines = []
    for b in range(m - 1, -1, -1):
        row = ""
        for a in range(m):
            v = grid.point(np.array([a, b]))
            if which is not None and mask[v]:
                row += str(which[v])
            else:
                row += "#" if mask[v] else "."
        lines.append(row)
    return "\n".join(lines)


def C(alpha, k):
    return 2 * math.sqrt((2 / k) * math.log2(1 / alpha))


def stamp(msg):
    print(f"{msg}  [{time.time() - T0:.1f}s]")


# CONTROLS


def controls(args):
    print("CONTROL 1: the EHPS block rasterised, EHPS form (classes), k = 2")
    for m in (12, 24, 36):
        grid = Grid(m, 2)
        eps = Fraction(1, 4 * m)
        best = None
        for i in range(4):
            for j in range(4):
                shift = (Fraction(i, 4 * m), Fraction(j, 4 * m))
                mask, which = rasterise(grid, eps, shift)
                if best is None or mask.sum() > best[0]:
                    best = (int(mask.sum()), shift, mask, which)
        cnt, shift, mask, which = best
        peel = grid.peelable(mask, "class")
        lp = grid.lp(mask, "class")
        fx = {}
        for v in np.flatnonzero(mask):
            a = Fraction(int(grid.coords[v, 0]), m) + shift[0]
            b = Fraction(int(grid.coords[v, 1]), m) + shift[1]
            fx[int(v)] = ehps_f(a, b, eps)
        rx, rz, ry, rd = grid.rows(mask)
        worst = min(
            (
                fx[x] + fx[z] - 2 * fx[y] - Fraction(int(d), m * m)
                for x, z, y, d in zip(
                    rx.tolist(), rz.tolist(), ry.tolist(), rd.tolist()
                )
            ),
            default=Fraction(0),
        )
        print(
            f"  m = {m}: {cnt} points, density {cnt / m**2:.5f} against 7/24 = {7 / 24:.5f}, eps = {eps}, shift = ({shift[0]}, {shift[1]}), {len(rx)} triples; peelable {peel}; LP {'feasible, max f = %.1f, min slack = %.2e' % (lp[1], lp[2]) if lp else 'infeasible'}; EHPS f worst slack {float(worst):.3e} {'>= 0 ok' if worst >= 0 else 'NEGATIVE'}"
        )
        assert peel and lp is not None and worst >= 0
        if m == 12:
            print(bitmap(grid, mask, which))
    stamp("control 1 passed")
    print(
        "CONTROL 2: k = 1, exact optimum by the lazy-core solver against (m + 1)/2 at odd m, m/2 + 1 at even m in EHPS form, m/2 in HS form"
    )
    for m in range(2, 21):
        grid = Grid(m, 1)
        for mode in ("class", "point"):
            if mode == "point" and m % 2:
                continue
            size, mask, bound, done = exact(grid, mode, 20)
            want = (
                (m + 1) // 2 if m % 2 else (m // 2 + 1 if mode == "class" else m // 2)
            )
            lp = grid.lp(mask, mode)
            print(
                f"  m = {m} {mode}: optimum {size} of {m}, formula {want}, {'exact' if done else 'bound ' + str(bound)}, T = {sorted(np.flatnonzero(mask).tolist())}, LP {'ok' if lp else 'FAIL'}"
            )
            assert done and size == want and lp is not None
    stamp("control 2 passed")
    print("CONTROL 3: peel oracle against the LP oracle on random sets, k = 2")
    rng = np.random.default_rng(1)
    agree = 0
    for m in (5, 6, 7, 8, 9, 12):
        grid = Grid(m, 2)
        for _ in range(12):
            mask = rng.random(grid.n) < rng.uniform(0.2, 0.5)
            for mode in ("class", "point"):
                if mode == "point":
                    for c in range(grid.ncls):
                        mem = np.flatnonzero(mask & (grid.cls == c))
                        mask[mem[1:]] = False
                p = grid.peelable(mask, mode)
                lp = grid.lp(mask, mode)
                t, lam, lam2, _rows, extra = grid.farkas(mask, mode)
                ok = (lp is not None) == p and ((t > 1e-9) == p)
                agree += ok
                if not ok:
                    print(
                        f"  DISAGREE m = {m} {mode}: peel {p}, lp {lp is not None}, margin t = {t:.3e}, |T| = {mask.sum()}"
                    )
                if p:
                    g = grid.certificate(mask, mode)
                    assert g is not None
    print(f"  {agree} of 144 cases agree (peel, LP feasibility, max-margin sign)")
    assert agree == 144
    grid = Grid(6, 2)
    mask = np.zeros(grid.n, dtype=bool)
    rng = np.random.default_rng(3)
    while True:
        mask = rng.random(grid.n) < 0.45
        if not grid.peelable(mask, "class"):
            break
    core = grid.minimal_core(mask, "class")
    t, lam, lam2, (rx, rz, ry, rd), extra = grid.farkas(mask, "class")
    bal = np.zeros(grid.n)
    np.add.at(bal, rx, lam)
    np.add.at(bal, rz, lam)
    np.add.at(bal, ry, -2 * lam)
    for (a, b), l in zip(extra, lam2):
        bal[a] += l
        bal[b] -= l
    print(
        f"  witness at m = 6: |T| = {mask.sum()}, minimal core {sorted(np.flatnonzero(core).tolist())} ({core.sum()} points), max margin t = {t:.3e}, Farkas weights on {int((lam > 1e-9).sum())} triples, balance max |.| = {np.abs(bal).max():.2e}, weighted distance = {float(np.dot(lam, rd)) / 36:.4f}"
    )
    print(bitmap(grid, core))
    stamp("control 3 passed")


# EXACT


def symmetries(grid):
    m, k = grid.m, grid.k
    maps = []
    co = grid.coords
    for perm in ([0, 1], [1, 0]) if k == 2 else ([0],):
        for signs in range(2**k):
            img = co[:, perm].copy()
            for i in range(k):
                if (signs >> i) & 1:
                    img[:, i] = (-img[:, i]) % m
            for t in range(grid.n):
                sh = (img + co[t]) % m
                maps.append(sh @ (m ** np.arange(k)[::-1]))
    return np.unique(np.stack(maps), axis=0)


def exact(grid, mode, cap, seed_cores=None, verbose=False):
    n = grid.n
    syms = symmetries(grid)
    cuts = set()
    if mode == "point":
        for c in range(grid.ncls):
            mem = np.flatnonzero(grid.cls == c)
            for i in range(len(mem)):
                for j in range(i + 1, len(mem)):
                    cuts.add((int(mem[i]), int(mem[j])))
    start = time.time()
    best_mask, bound, it = None, n, 0
    while True:
        it += 1
        if cuts:
            rows, cols = [], []
            rhs = []
            for r, cut in enumerate(cuts):
                rows += [r] * len(cut)
                cols += list(cut)
                rhs.append(len(cut) - 1)
            A = coo_matrix(
                (np.ones(len(rows)), (rows, cols)), shape=(len(cuts), n)
            ).tocsr()
            cons = [LinearConstraint(A, -np.inf, np.array(rhs, dtype=float))]
        else:
            cons = []
        res = milp(
            -np.ones(n),
            constraints=cons,
            integrality=np.ones(n),
            bounds=Bounds(0, 1),
            options={"time_limit": max(1.0, cap - (time.time() - start))},
        )
        if res.x is None:
            return (
                (int(best_mask.sum()) if best_mask is not None else 0),
                best_mask,
                bound,
                False,
            )
        mask = res.x > 0.5
        bound = min(bound, round(-res.fun))
        if grid.peelable(mask, mode):
            return int(mask.sum()), mask, int(mask.sum()), True
        core = grid.minimal_core(mask, mode)
        cv = np.flatnonzero(core)
        for s in syms:
            cuts.add(tuple(sorted(int(x) for x in s[cv])))
        if verbose and it % 20 == 0:
            stamp(
                f"    iteration {it}: {len(cuts)} cuts, bound {bound}, last core {len(cv)} points"
            )
        if time.time() - start > cap:
            return 0, None, bound, False


def exact_cmd(args):
    k = args.k
    for m in range(args.lo, args.hi + 1):
        grid = Grid(m, k)
        for mode in ("class", "point"):
            if mode == "point" and m % 2:
                continue
            size, mask, bound, done = exact(grid, mode, args.cap, verbose=True)
            label = "EHPS form" if mode == "class" else "HS form"
            if done:
                alpha = size / m**k
                print(
                    f"m = {m} k = {k} {label}: exact optimum {size} of {m**k}, alpha = {alpha:.5f} = {size}/{m**k}, 7/24 = {7 / 24:.5f}, C = {C(alpha, k):.4f}"
                )
                if k == 2:
                    print(bitmap(grid, mask))
            else:
                print(
                    f"m = {m} k = {k} {label}: not closed in {args.cap}s, upper bound {bound} of {m**k}, alpha <= {bound / m**k:.5f}"
                )
        stamp(f"m = {m} done")


# SEARCH


def strip(grid, mask, mode):
    mask = mask.copy()
    if mode == "point":
        for c in range(grid.ncls):
            mem = np.flatnonzero(mask & (grid.cls == c))
            mask[mem[1:]] = False
    while True:
        core = grid.core(mask, mode)
        if not core.any():
            return mask
        _rx, _rz, ry, _ = grid.rows(core)
        cnt = np.bincount(ry, minlength=grid.n)
        mask[int(np.argmax(cnt))] = False


def extend(grid, mask, mode, rng, deadline, ok=None):
    ok = ok or (lambda mk: grid.peelable(mk, mode))
    mask = mask.copy()
    out = np.flatnonzero(~mask)
    rng.shuffle(out)
    for v in out:
        if time.time() > deadline:
            break
        if mode == "point" and mask[grid.cls == grid.cls[v]].any():
            continue
        mask[v] = True
        if not ok(mask):
            mask[v] = False
    return mask


def perturb(grid, mask, mode, rng, kick, core_fn=None):
    core_fn = core_fn or (lambda mk: grid.core(mk, mode))
    mask = mask.copy()
    out = np.flatnonzero(~mask)
    for v in rng.choice(out, size=min(kick, len(out)), replace=False):
        if mode == "point":
            mask[grid.cls == grid.cls[v]] = False
        mask[v] = True
        while True:
            core = core_fn(mask)
            if not core.any():
                break
            cand = np.flatnonzero(core)
            cand = cand[cand != v]
            if len(cand) == 0:
                mask[v] = False
                break
            mask[rng.choice(cand)] = False
    return mask


def greedy(grid, mode, order, rng):
    m, k, n = grid.m, grid.k, grid.n
    mask = np.zeros(n, dtype=bool)
    mids = np.zeros(n, dtype=int)
    co = grid.coords
    weights = m ** np.arange(k)[::-1]
    half = m // 2
    shifts = (
        [
            np.array([(e >> (k - 1 - i)) & 1 for i in range(k)]) * half
            for e in range(2**k)
        ]
        if m % 2 == 0
        else [np.zeros(k, dtype=int)]
    )
    inv2 = pow(2, -1, m) if m % 2 else None

    def midpoints(p, q):
        s = co[p] + co[q]
        if m % 2:
            return [int(((s * inv2) % m) @ weights)]
        if np.any(s % 2):
            return []
        return [int((((s // 2) + sh) % m) @ weights) for sh in shifts]

    def add(p):
        for q in np.flatnonzero(mask):
            for y in midpoints(p, q):
                mids[y] += 1
        mask[p] = True

    def addable(S):
        for y in S:
            if mids[y]:
                return False
            for z in S:
                if z != y:
                    x = int(((2 * co[y] - co[z]) % m) @ weights)
                    if mask[x] or x in S:
                        return False
        return True

    for c in order:
        mem = [int(v) for v in np.flatnonzero(grid.cls == c)]
        if mode == "point":
            rng.shuffle(mem)
            subsets = [[v] for v in mem]
        else:
            subsets = [
                [mem[i] for i in range(len(mem)) if (b >> i) & 1]
                for b in range(1, 2 ** len(mem))
            ]
            rng.shuffle(subsets)
            subsets.sort(key=len, reverse=True)
        for S in subsets:
            if addable(S):
                for v in S:
                    add(v)
                break
    return mask


def ehps_seed(grid):
    m = grid.m
    best = None
    for den in (4,):
        eps = Fraction(1, den * m)
        for i in range(den):
            for j in range(den):
                mask, _ = rasterise(
                    grid, eps, (Fraction(i, den * m), Fraction(j, den * m))
                )
                if best is None or mask.sum() > best.sum():
                    best = mask
    return best


def search(grid, mode, seconds, rng, seeds, log=None):
    deadline = time.time() + seconds
    best = None
    starts = []
    for name, mask in seeds:
        starts.append((name, strip(grid, mask, mode)))
    n_rand = 0
    while True:
        if starts:
            name, mask = starts.pop(0)
        else:
            if time.time() > deadline:
                break
            n_rand += 1
            order = np.arange(grid.ncls)
            rng.shuffle(order)
            name, mask = f"random greedy {n_rand}", greedy(grid, mode, order, rng)
        mask = extend(grid, mask, mode, rng, deadline)
        local_best = mask.copy()
        stale = 0
        while time.time() < deadline and stale < 30:
            cand = extend(
                grid,
                perturb(grid, mask, mode, rng, 1 + int(rng.integers(3))),
                mode,
                rng,
                deadline,
            )
            if cand.sum() >= mask.sum():
                stale = stale + 1 if cand.sum() == mask.sum() else 0
                mask = cand
                if mask.sum() > local_best.sum():
                    local_best = mask.copy()
            else:
                stale += 1
        if best is None or local_best.sum() > best[1].sum():
            best = (name, local_best)
            if log:
                log(f"    seed {name}: {int(local_best.sum())}, new best")
        elif log and n_rand % 25 == 0:
            log(f"    {n_rand} random seeds, best {int(best[1].sum())}")
    return best


def search_cmd(args):
    k = args.k
    rng = np.random.default_rng(args.seed)
    results = {}
    for m in args.m:
        grid = Grid(m, k)
        modes = ["class"] if m % 2 else ["class", "point"]
        if args.hs_only:
            modes = ["point"]
        if args.ehps_only:
            modes = ["class"]
        for mode in modes:
            label = "EHPS form" if mode == "class" else "HS form"
            seeds = []
            if k == 2 and mode == "class":
                seeds.append(("EHPS rasterised", ehps_seed(grid)))
            name, mask = search(
                grid,
                mode,
                args.seconds,
                rng,
                seeds,
                log=stamp if args.verbose else None,
            )
            size = int(mask.sum())
            alpha = size / m**k
            assert grid.peelable(mask, mode)
            g = grid.certificate(mask, mode)
            lp = grid.lp(mask, mode)
            assert g is not None and lp is not None
            levels = max(int(math.log(v, 3)) for v in g.values()) + 1
            print(
                f"m = {m} k = {k} {label}: best {size} of {m**k}, alpha = {alpha:.5f}, 7 m^2/24 = {7 * m**k / 24:.2f}, EHPS raster {int(ehps_seed(grid).sum()) if k == 2 and mode == 'class' else '-'}, C = {C(alpha, k):.4f}, c = alpha^(1/k) = {alpha ** (1 / k):.4f}, peel levels {levels}, LP max f = {lp[1]:.2f}, from {name}"
            )
            if k == 2:
                print(bitmap(grid, mask))
            results[f"{m}-{mode}"] = {
                "m": m,
                "k": k,
                "form": label,
                "size": size,
                "alpha": alpha,
                "T": sorted(int(v) for v in np.flatnonzero(mask)),
                "coords": [
                    [int(c) for c in grid.coords[v]] for v in np.flatnonzero(mask)
                ],
                "f_lp": {str(v): fv for v, fv in lp[0].items()},
                "lp_max_f": lp[1],
                "peel_levels": levels,
                "g_peel": {str(v): str(gv) for v, gv in g.items()},
                "seed": name,
            }
        stamp(f"m = {m} done")
    if args.save:
        path = os.path.join(HERE, args.save)
        with open(path, "w") as fh:
            json.dump(results, fh, indent=1)
        stamp(f"saved {path}")


# REFINE


def dilate(grid, mask, j):
    big = Grid(grid.m * j, grid.k)
    out = np.zeros(big.n, dtype=bool)
    for v in np.flatnonzero(mask):
        base = grid.coords[v] * j
        for e in range(j**grid.k):
            off = np.array([(e // j ** (grid.k - 1 - i)) % j for i in range(grid.k)])
            out[big.point(base + off)] = True
    return big, out


def refine_cmd(args):
    with open(os.path.join(HERE, args.load)) as fh:
        saved = json.load(fh)
    rng = np.random.default_rng(args.seed)
    for rec in saved.values():
        if rec["form"] != "EHPS form" or rec["k"] != 2:
            continue
        m = rec["m"]
        grid = Grid(m, 2)
        mask = np.zeros(grid.n, dtype=bool)
        mask[rec["T"]] = True
        assert grid.peelable(mask, "class")
        print(
            f"m = {m}: saved set of {int(mask.sum())} points, alpha = {mask.sum() / m**2:.5f}"
        )
        for j in args.factors:
            big, dil = dilate(grid, mask, j)
            ok = big.peelable(dil, "class")
            core = big.minimal_core(dil, "class") if not ok else None
            line = f"  dilated by {j} to m = {m * j}: {int(dil.sum())} points, density {dil.sum() / big.n:.5f}, peelable {ok}"
            if not ok:
                line += f", minimal core of {int(core.sum())} points at {[tuple(int(c) for c in big.coords[v]) for v in np.flatnonzero(core)]}"
            print(line)
            if not ok and args.repair > 0:
                deadline = time.time() + args.repair
                rep = extend(big, strip(big, dil, "class"), "class", rng, deadline)
                local = rep.copy()
                while time.time() < deadline:
                    cand = extend(
                        big, perturb(big, rep, "class", rng, 2), "class", rng, deadline
                    )
                    if cand.sum() >= rep.sum():
                        rep = cand
                        if rep.sum() > local.sum():
                            local = rep.copy()
                print(
                    f"    repaired at m = {m * j} in {args.repair}s: {int(local.sum())} points, density {local.sum() / big.n:.5f}"
                )
        stamp(f"m = {m} refined")


# STABLE


class Stable:
    def __init__(self, grid, factors):
        self.grid = grid
        self.levels = [(j, Grid(grid.m * j, grid.k)) for j in factors]

    def lift(self, mask, j, big):
        out = np.zeros(big.n, dtype=bool)
        pts = np.flatnonzero(mask)
        for e in range(j**self.grid.k):
            off = np.array(
                [(e // j ** (self.grid.k - 1 - i)) % j for i in range(self.grid.k)]
            )
            out[
                ((self.grid.coords[pts] * j + off) @ (big.m ** np.arange(big.k)[::-1]))
            ] = True
        return out

    def ok(self, mask):
        if not self.grid.peelable(mask, "class"):
            return False
        return all(
            big.peelable(self.lift(mask, j, big), "class") for j, big in self.levels
        )

    def core(self, mask):
        core = self.grid.core(mask, "class")
        if core.any():
            return core
        for j, big in self.levels:
            c = big.core(self.lift(mask, j, big), "class")
            if c.any():
                base = np.zeros(self.grid.n, dtype=bool)
                base[
                    (big.coords[np.flatnonzero(c)] // j)
                    @ (self.grid.m ** np.arange(self.grid.k)[::-1])
                ] = True
                return base
        return core

    def strip(self, mask):
        mask = mask.copy()
        while True:
            core = self.core(mask)
            if not core.any():
                return mask
            mask[
                np.random.default_rng(int(mask.sum())).choice(np.flatnonzero(core))
            ] = False


def stable_cmd(args):
    rng = np.random.default_rng(args.seed)
    results = {}
    for m in args.m:
        grid = Grid(m, 2)
        st = Stable(grid, args.factors)
        deadline = time.time() + args.seconds
        seeds = [("EHPS rasterised", ehps_seed(grid))]
        best = None
        n_rand = 0
        while True:
            if seeds:
                name, mask = seeds.pop(0)
            else:
                if time.time() > deadline:
                    break
                n_rand += 1
                order = np.arange(grid.ncls)
                rng.shuffle(order)
                name, mask = (
                    f"random greedy {n_rand}",
                    greedy(grid, "class", order, rng),
                )
            mask = extend(grid, st.strip(mask), "class", rng, deadline, ok=st.ok)
            local = mask.copy()
            stale = 0
            while time.time() < deadline and stale < 20:
                cand = extend(
                    grid,
                    perturb(
                        grid,
                        mask,
                        "class",
                        rng,
                        1 + int(rng.integers(3)),
                        core_fn=st.core,
                    ),
                    "class",
                    rng,
                    deadline,
                    ok=st.ok,
                )
                if cand.sum() >= mask.sum():
                    stale = stale + 1 if cand.sum() == mask.sum() else 0
                    mask = cand
                    if mask.sum() > local.sum():
                        local = mask.copy()
                else:
                    stale += 1
            if best is None or local.sum() > best[1].sum():
                best = (name, local)
                stamp(f"    seed {name}: {int(local.sum())}, new best")
        name, mask = best
        assert st.ok(mask)
        size = int(mask.sum())
        lp = grid.lp(mask, "class")
        print(
            f"m = {m} stable under dilation by {args.factors}: best {size} of {m**2}, alpha = {size / m**2:.5f}, 7 m^2/24 = {7 * m**2 / 24:.2f}, C = {C(size / m**2, 2):.4f}, LP max f = {lp[1]:.2f}, from {name}"
        )
        print(bitmap(grid, mask))
        for j in args.test:
            big = Grid(m * j, 2)
            dil = st.lift(mask, j, big)
            ok = big.peelable(dil, "class")
            line = f"  dilated by {j} to m = {m * j}: {int(dil.sum())} points, peelable {ok}"
            if not ok:
                core = big.minimal_core(dil, "class")
                line += f", minimal core of {int(core.sum())} points at {[tuple(int(c) for c in big.coords[v]) for v in np.flatnonzero(core)]}"
            print(line)
        results[str(m)] = {
            "m": m,
            "factors": args.factors,
            "size": size,
            "alpha": size / m**2,
            "T": sorted(int(v) for v in np.flatnonzero(mask)),
            "coords": [[int(c) for c in grid.coords[v]] for v in np.flatnonzero(mask)],
            "f_lp": {str(v): fv for v, fv in lp[0].items()},
            "lp_max_f": lp[1],
            "seed": name,
        }
        stamp(f"m = {m} done")
    if args.save:
        path = os.path.join(HERE, args.save)
        with open(path, "w") as fh:
            json.dump(results, fh, indent=1)
        stamp(f"saved {path}")


# LIFT


def lift_cmd(args):
    import itertools

    for m in args.m:
        grid = Grid(m, 2)
        saved = None
        if args.load and os.path.exists(os.path.join(HERE, args.load)):
            with open(os.path.join(HERE, args.load)) as fh:
                saved = json.load(fh).get(f"{m}-class")
        if m <= 6:
            _size, mask, _bound, done = exact(grid, "class", args.cap)
            assert done
            origin = "exact optimum"
        elif saved:
            mask = np.zeros(grid.n, dtype=bool)
            mask[saved["T"]] = True
            origin = f"saved search best from {args.load}"
        else:
            rng = np.random.default_rng(args.seed)
            name, mask = search(grid, "class", args.seconds, rng, [])
            origin = f"search best ({name})"
        assert m % 2 == 1
        g = grid.certificate(mask, "class")
        pts = sorted(g)
        G = {v: 2 * g[v] for v in pts}
        maxG = max(G.values())
        t = len(pts)
        print(
            f"m = {m}: block of {t} points, {origin}, certificate levels {max(g.values())} = 3^{int(math.log(max(g.values()), 3) + 0.5)}, bound r_3(Z_{m}^n) >= {t}^(n/2) / ({maxG // 2} n + 1) = {math.sqrt(t):.4f}^n / O(n) against the product base (m + 1)/2 = {(m + 1) / 2}"
        )
        for n in args.n:
            blocks = n // 2
            tot = {}
            for combo in itertools.product(pts, repeat=blocks):
                s = sum(G[v] for v in combo)
                tot.setdefault(s, []).append(combo)
            s_best = max(tot, key=lambda s: len(tot[s]))
            S = tot[s_best]
            vecs = np.array(
                [[int(c) for v in combo for c in grid.coords[v]] for combo in S]
            )
            key = {tuple(row): i for i, row in enumerate(vecs.tolist())}
            aps = 0
            for i in range(len(S)):
                for j in range(i + 1, len(S)):
                    y = tuple(((vecs[i] + vecs[j]) * pow(2, -1, m)) % m)
                    if y in key:
                        aps += 1
            floor = t**blocks / ((maxG // 2) * n + 1)
            print(
                f"  n = {n}: {t**blocks} lifted points in {len(tot)} slices, largest slice {len(S)} points at level {s_best}, progressions found by brute force {aps}, guaranteed slice size {floor:.1f}, ((m + 1)/2)^n = {((m + 1) // 2) ** n}"
            )
            assert aps == 0 and len(S) >= floor
        stamp(f"m = {m} lifted")


# BRUTE


def brute_cmd(args):
    for m in args.m:
        grid = Grid(m, 2)
        n = grid.n
        mids = {}
        for x, z, y in zip(grid.X.tolist(), grid.Z.tolist(), grid.Y.tolist()):
            mids.setdefault((x, z), []).append(y)
        cls = grid.cls.tolist()

        def peelable(T):
            alive = set(T)
            while alive:
                cnt = {v: 0 for v in alive}
                pts = sorted(alive)
                for i in range(len(pts)):
                    for j in range(i + 1, len(pts)):
                        for y in mids.get((pts[i], pts[j]), ()):
                            if y in alive:
                                cnt[y] += 1
                bad = {cls[v] for v in alive if cnt[v]}
                rem = [v for v in alive if cls[v] not in bad]
                if not rem:
                    return False
                alive.difference_update(rem)
            return True

        best, count = [0], [0]

        def dfs(T, start):
            count[0] += 1
            if len(T) > len(best):
                best[:] = T
            for v in range(start, n):
                T.append(v)
                if peelable(T):
                    dfs(T, v + 1)
                T.pop()

        dfs([0], 1)
        print(
            f"m = {m}: every peelable set containing the origin enumerated by depth-first search with hereditary pruning, {count[0]} sets, largest {len(best)} of {n}, alpha = {len(best) / n:.5f}, witness {sorted(best)}"
        )
        stamp(f"m = {m} brute")


# PROGRESSION-FREE


def r3_cmd(args):
    for m in args.m:
        grid = Grid(m, args.k)
        n = grid.n
        rows = np.repeat(np.arange(len(grid.X)), 3)
        cols = np.stack([grid.X, grid.Y, grid.Z], axis=1).ravel()
        A = coo_matrix(
            (np.ones(len(rows)), (rows, cols)), shape=(len(grid.X), n)
        ).tocsr()
        res = milp(
            -np.ones(n),
            constraints=[LinearConstraint(A, -np.inf, 2)],
            integrality=np.ones(n),
            bounds=Bounds(0, 1),
            options={"time_limit": args.cap},
        )
        mask = res.x > 0.5
        assert not np.any(mask[grid.X] & mask[grid.Y] & mask[grid.Z])
        size = int(mask.sum())
        closed = res.status == 0
        tail = "" if closed else f", <= {int(np.floor(-res.mip_dual_bound + 1e-6))}"
        print(
            f"m = {m} k = {args.k}: r_3(Z_{m}^{args.k}) {'=' if closed else '>='} {size}{tail}, {len(grid.X)} triples, the set checked progression-free, constant {size ** (1 / args.k):.4f}"
        )


def eps_def34(a, b):
    half = Fraction(1, 2)
    s = a + b
    if a < half and b >= half and Fraction(7, 12) <= s <= Fraction(4, 3):
        return True
    if a < half and b < half and s > Fraction(5, 6):
        return True
    return (
        a >= half
        and b < half
        and Fraction(7, 12) <= s < Fraction(5, 6)
        and 2 * a + b < Fraction(3, 2)
    )


def region_mask(grid, region, c, n, i, j):
    m = grid.m
    scale = 12 * n * m * c
    a = 12 * c * (n * grid.coords[:, 0] + i)
    b = 12 * c * (n * grid.coords[:, 1] + j)
    s = a + b
    half, one = scale // 2, scale
    if region == "EPS":
        a_lo = a < half
        b_lo = b < half
        t1 = a_lo & ~b_lo & (12 * s >= 7 * scale) & (3 * s <= 4 * scale)
        t1b = a_lo & b_lo & (6 * s > 5 * scale)
        t2 = (
            ~a_lo
            & b_lo
            & (12 * s >= 7 * scale)
            & (6 * s < 5 * scale)
            & (2 * (2 * a + b) < 3 * scale)
        )
        return t1 | t1b | t2
    e = 12 * n
    r1 = (a >= half) & (3 * s > 2 * one) & (6 * s <= 7 * one)
    r2 = (a >= half) & (b < half) & (6 * s >= 7 * one + 6 * e) & (12 * s <= 17 * one)
    r3 = (
        (a < half)
        & (b >= half)
        & (6 * s >= 7 * one + 6 * e)
        & (12 * s <= 17 * one)
        & (2 * (2 * a + b) >= 3 * one + 2 * e)
    )
    return r1 | r2 | r3


def raster_cmd(args):
    n = args.grid
    print(
        f"RASTER: the published blocks counted on the grid, k = 2, EHPS form: EHPS Definition 4.1 at eps = 1/(4m), 1/(16m), 1/(64m) and EPS Definition 3.4, at the {n} x {n} shifts (i, j)/({n} m), i, j < {n}, first maximiser in that order"
    )
    for m in args.m:
        grid = Grid(m, 2)
        found = []
        for region, c in (("EHPS", 4), ("EHPS", 16), ("EHPS", 64), ("EPS", 1)):
            top = None
            for i in range(n):
                for j in range(n):
                    cnt = int(region_mask(grid, region, c, n, i, j).sum())
                    if top is None or cnt > top[0]:
                        top = (cnt, i, j)
            found.append((top[0], region, c, top[1], top[2]))
        cnt, region, c, i, j = max(found, key=lambda f: f[0])
        mask = region_mask(grid, region, c, n, i, j)
        shift = (Fraction(i, n * m), Fraction(j, n * m))
        if region == "EPS":
            oracle = np.array(
                [
                    eps_def34(
                        Fraction(int(x), m) + shift[0], Fraction(int(y), m) + shift[1]
                    )
                    for x, y in grid.coords
                ]
            )
        else:
            oracle = rasterise(grid, Fraction(1, c * m), shift)[0]
        assert (oracle == mask).all()
        assert grid.peelable(mask, "class")
        g = grid.certificate(mask, "class")
        assert g is not None
        rounds = grid.rounds(mask, "class")[1]
        bound = -(-7 * m * m // 24)
        assert cnt >= bound
        root = math.floor(math.sqrt(cnt) * 1e4) / 1e4
        ehps = max(f[0] for f in found if f[1] == "EHPS")
        eps = max(f[0] for f in found if f[1] == "EPS")
        which = (
            "EPS Definition 3.4"
            if region == "EPS"
            else f"EHPS Definition 4.1, eps = 1/({c} m)"
        )
        print(
            f"  m = {m}: {cnt} points, {which}, shift ({i}, {j})/({n} m), ceil(7 m^2/24) = {bound}, excess {cnt - bound}; best EHPS {ehps}, best EPS {eps}; peelable, {rounds} peel rounds, certificate g <= {max(g.values())} checked on every triple, sqrt(card T) = {root:.4f}.."
        )


# MAIN


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "verb",
        choices=[
            "controls",
            "exact",
            "search",
            "refine",
            "stable",
            "lift",
            "brute",
            "r3",
            "raster",
        ],
    )
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--lo", type=int, default=2)
    ap.add_argument("--hi", type=int, default=8)
    ap.add_argument("--cap", type=float, default=60.0)
    ap.add_argument("--m", type=int, nargs="+", default=[12])
    ap.add_argument("--seconds", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--save", default=None)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--hs-only", action="store_true")
    ap.add_argument("--ehps-only", action="store_true")
    ap.add_argument("--load", default="best.json")
    ap.add_argument("--factors", type=int, nargs="+", default=[2, 3])
    ap.add_argument("--repair", type=float, default=0.0)
    ap.add_argument("--test", type=int, nargs="+", default=[4, 6])
    ap.add_argument("--n", type=int, nargs="+", default=[4, 6])
    ap.add_argument("--grid", type=int, default=12)
    args = ap.parse_args()
    {
        "controls": controls,
        "exact": exact_cmd,
        "search": search_cmd,
        "refine": refine_cmd,
        "stable": stable_cmd,
        "lift": lift_cmd,
        "brute": brute_cmd,
        "r3": r3_cmd,
        "raster": raster_cmd,
    }[args.verb](args)
    stamp("done")


if __name__ == "__main__":
    main()
