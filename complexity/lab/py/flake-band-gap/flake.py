import sys
import time
from fractions import Fraction
from itertools import pairwise
from math import lcm, sqrt

import numpy as np

CORNERS = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=np.int64)
TINY = 1e-13


def cells(level):
    pts = np.zeros((1, 3), dtype=np.int64)
    for _ in range(level):
        pts = (2 * pts[:, None, :] + CORNERS[None, :, :]).reshape(-1, 3)
    return pts


def face_pairs(pts, level):
    bits = level + 1
    key = (pts[:, 0] << (2 * bits)) | (pts[:, 1] << bits) | pts[:, 2]
    rank = np.argsort(key)
    ordered = key[rank]
    src, dst = [], []
    for step in (1 << (2 * bits), 1 << bits, 1):
        want = key + step
        slot = np.clip(np.searchsorted(ordered, want), 0, key.size - 1)
        hit = ordered[slot] == want
        src.append(np.nonzero(hit)[0])
        dst.append(rank[slot[hit]])
    return np.concatenate(src), np.concatenate(dst)


def layers(count, indptr, nbr):
    seen = np.zeros(count, dtype=bool)
    seen[0] = True
    front = np.array([0], dtype=np.int64)
    groups, parents = [front], [np.array([-1], dtype=np.int64)]
    while True:
        span = indptr[front + 1] - indptr[front]
        total = int(span.sum())
        if total == 0:
            break
        base = np.repeat(indptr[front], span)
        step = np.arange(total) - np.repeat(np.cumsum(span) - span, span)
        cand = nbr[base + step]
        came = np.repeat(front, span)
        fresh = ~seen[cand]
        cand, came = cand[fresh], came[fresh]
        if cand.size == 0:
            break
        seen[cand] = True
        groups.append(cand)
        parents.append(came)
        front = cand
    return groups, parents, int(seen.sum())


def rooted(level):
    pts = cells(level)
    count = pts.shape[0]
    src, dst = face_pairs(pts, level)
    tail = np.concatenate([src, dst])
    deg = np.bincount(tail, minlength=count)
    indptr = np.zeros(count + 1, dtype=np.int64)
    np.cumsum(deg, out=indptr[1:])
    nbr = np.concatenate([dst, src])[np.argsort(tail, kind="stable")]
    groups, parents, reached = layers(count, indptr, nbr)
    label = np.concatenate(groups)
    place = np.empty(count, dtype=np.int64)
    place[label] = np.arange(count)
    offset = np.zeros(len(groups) + 1, dtype=np.int64)
    np.cumsum([g.size for g in groups], out=offset[1:])
    local = [np.array([], dtype=np.int64)]
    for d in range(1, len(groups)):
        local.append(place[parents[d]] - offset[d - 1])
    return {
        "count": count,
        "deg": deg[label],
        "offset": offset,
        "local": local,
        "tree": reached == count and int(src.size) == count - 1,
        "edges": (place[src], place[dst]),
    }


def below(tree, shift):
    offset, degree = tree["offset"], tree["deg"].astype(np.float64)
    acc = np.zeros(tree["count"])
    neg = flat = 0
    for d in range(offset.size - 2, -1, -1):
        lo, hi = offset[d], offset[d + 1]
        piv = degree[lo:hi] - shift - acc[lo:hi]
        small = np.abs(piv) < TINY
        hits = int(small.sum())
        if hits:
            flat += hits
            piv = np.where(small, np.where(piv >= 0.0, TINY, -TINY), piv)
        neg += int((piv < 0.0).sum())
        if d:
            acc[offset[d - 1] : offset[d]] += np.bincount(
                tree["local"][d], weights=1.0 / piv, minlength=int(lo - offset[d - 1])
            )
    return neg, flat


def edge(tree, target, lo, hi):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if below(tree, mid)[0] >= target:
            hi = mid
        else:
            lo = mid
    return hi


def exact_pivots(tree, shift):
    offset, degree = tree["offset"], tree["deg"]
    acc = [Fraction(0)] * tree["count"]
    piv = [Fraction(0)] * tree["count"]
    neg = 0
    for d in range(offset.size - 2, -1, -1):
        lo, hi = int(offset[d]), int(offset[d + 1])
        up, home = int(offset[d - 1]), tree["local"][d]
        for v in range(lo, hi):
            here = Fraction(int(degree[v])) - shift - acc[v]
            piv[v] = here
            if here < 0:
                neg += 1
            if d:
                if here == 0:
                    return neg, piv, False
                acc[up + int(home[v - lo])] += Fraction(1) / here
    return neg, piv, True


def null_vector(tree, piv):
    offset = tree["offset"]
    val = [Fraction(0)] * tree["count"]
    val[0] = Fraction(1)
    for d in range(1, offset.size - 1):
        lo, hi = int(offset[d]), int(offset[d + 1])
        up, home = int(offset[d - 1]), tree["local"][d]
        for v in range(lo, hi):
            val[v] = val[up + int(home[v - lo])] / piv[v]
    scale = 1
    for x in val:
        scale = lcm(scale, x.denominator)
    return [int(x * scale) for x in val]


def residual(tree, vec):
    res = [(int(tree["deg"][v]) - 4) * x for v, x in enumerate(vec)]
    for a, b in zip(*tree["edges"]):
        res[a] -= vec[b]
        res[b] -= vec[a]
    return max(abs(x) for x in res)


def main():
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    exact_top = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    vector_top = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    print(f"flake band gap, code 23 base 2, levels 1..{top}")
    print("\n    L           N  tree  #(eig<2)  #(eig<4)    target  flat  lo(L)")
    los, his, above = {}, {}, {}
    for level in range(1, top + 1):
        clock = time.time()
        tree = rooted(level)
        target = 3 * 4 ** (level - 1)
        b2, t2 = below(tree, 2.0)
        b4, t4 = below(tree, 4.0)
        los[level] = edge(tree, target, 1.0, 3.0)
        his[level] = edge(tree, tree["count"], 4.0, 8.0)
        above[level] = tree["count"] - b4
        good = tree["tree"] and b2 == b4 == target and t2 == 0 and t4 == 1
        print(
            f"{level:5}  {tree['count']:10}  {'yes' if tree['tree'] else 'NO':>4}  {b2:8}  "
            f"{b4:8}  {target:8}  {t2 + t4:4}  {los[level]:.15f}  "
            f"[{time.time() - clock:.1f}s] {'PASS' if good else 'FAIL'}"
        )
    print("\n    L              hi(L)  #(eig>=4)    4^(L-1)")
    for level, upper in his.items():
        split = above[level] == 4 ** (level - 1)
        print(
            f"{level:5}  {upper:.10f}  {above[level]:9}  {4 ** (level - 1):9}  "
            f"{'PASS' if split else 'FAIL'}"
        )
    if 2 in his:
        print(f"hi(2) against 3+sqrt(5): {his[2] - (3.0 + sqrt(5.0)):.1e}")

    print("\n    L         2-lo(L)   ratio   (2-lo)*8^L")
    prev, ratios = None, []
    for level, lo in los.items():
        defect = 2.0 - lo
        ratio = None if prev is None else prev / defect
        ratios += [] if ratio is None else [ratio]
        shown = "-" if ratio is None else f"{ratio:.4f}"
        print(f"{level:5}  {defect:.12f}  {shown:>6}  {defect * 8.0**level:11.6f}")
        prev = defect
    climb = all(a < b < 8.0 for a, b in pairwise(ratios))
    print(f"ratios rise monotonically and stay under 8: {'PASS' if climb else 'FAIL'}")
    print(f"(2-lo)*8^L at L={top}: {(2.0 - los[top]) * 8.0**top:.6f}")
    print("\nexact rational elimination")
    print(
        "    L      N  root pivot at 4  earlier zeros  #(eig<2)  #(eig<4-d)  #(eig<4+d)  mult"
    )
    delta = Fraction(1, 10**9)
    for level in range(1, exact_top + 1):
        tree = rooted(level)
        piv, clean = exact_pivots(tree, Fraction(4))[1:]
        n2 = exact_pivots(tree, Fraction(2))[0]
        low = exact_pivots(tree, Fraction(4) - delta)[0]
        high = exact_pivots(tree, Fraction(4) + delta)[0]
        good = (
            clean
            and piv[0] == 0
            and n2 == low == 3 * 4 ** (level - 1)
            and high - low == 1
        )
        print(
            f"{level:5}  {tree['count']:5}  {piv[0]!s:>15}  {'none' if clean else 'SOME':>13}"
            f"  {n2:8}  {low:10}  {high:10}  {high - low:4}  {'PASS' if good else 'FAIL'}"
        )
    print("\ninteger null vector at 4")
    print("    L      N  digits  max|(Lap-4I)v|")
    for level in range(1, vector_top + 1):
        tree = rooted(level)
        vec = null_vector(tree, exact_pivots(tree, Fraction(4))[1])
        worst = residual(tree, vec)
        digits = len(str(max(abs(x) for x in vec)))
        print(
            f"{level:5}  {tree['count']:5}  {digits:6}  {worst:14}  "
            f"{'PASS' if worst == 0 else 'FAIL'}"
        )


main()
