import resource
import sys
import time
from fractions import Fraction
from itertools import product
from math import log

# CELLS

CELLS = [
    ("G", ((4, (0, 1)), (8, (0, 1, 2, 3)))),
    ("H", ((4, (0, 1)), (16, (0, 1, 4, 5)))),
    ("N", ((9, (0, 1, 2)), (27, tuple(range(9))))),
    ("T", ((4, (0, 1)), (8, (0, 1, 2, 3)), (16, tuple(range(8))))),
]

HEIGHT = 10**13

SHARP = ((2, (1,)), (4, (0, 1, 2, 3)))

IRRATIONAL = ((4, (0, 1, 2)), (16, tuple(range(16))))

# FRAME


def root_of(bases):
    for cand in range(2, min(bases) + 1):
        ok = True
        for b in bases:
            x = b
            while x % cand == 0:
                x //= cand
            if x != 1:
                ok = False
                break
        if ok:
            return cand
    return None


def exponent_of(b, r):
    e = 0
    while b % r == 0:
        b //= r
        e += 1
    return e


def lcm(values):
    out = 1
    for v in values:
        g, w = out, v
        while w:
            g, w = w, g % w
        out = out * v // g
    return out


def frame(cell):
    bases = [b for b, _ in cell]
    r = root_of(bases)
    exps = [exponent_of(b, r) for b in bases]
    return r, exps, lcm(exps)


# BLOCKS


def block_digits(cell):
    r, exps, m = frame(cell)
    keeps = [set(a) for _, a in cell]
    out = []
    for word in product(range(r), repeat=m):
        ok = True
        for keep, e in zip(keeps, exps):
            for g in range(m // e):
                v = 0
                for t in range(e - 1, -1, -1):
                    v = v * r + word[g * e + t]
                if v not in keep:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            value = 0
            for t in range(m - 1, -1, -1):
                value = value * r + word[t]
            out.append(value)
    return r, exps, m, tuple(sorted(out))


def in_design(n, b, keep):
    while n:
        if n % b not in keep:
            return False
        n //= b
    return True


def bottom_block(cell):
    r, _, m = frame(cell)
    keeps = [(b, set(a)) for b, a in cell]
    return tuple(
        c for c in range(r**m) if all(in_design(c, b, keep) for b, keep in keeps)
    )


# DIMENSION


def log_exact(value, r):
    t, x = 0, value
    while x % r == 0:
        x //= r
        t += 1
    return Fraction(t) if x == 1 else None


def design_dim(b, allowed, r):
    t = log_exact(len(allowed), r)
    return None if t is None else t / exponent_of(b, r)


def dim_text(count, r, m):
    t = log_exact(count, r)
    if t is None:
        return "log_%d(%d) / %d" % (r, count, m), log(count) / (m * log(r))
    return str(t / m), float(t / m)


# ENUMERATION


def walk(base, digits, limit):
    ds = sorted(digits)
    level = [d for d in ds if d and d < limit]
    while level:
        for v in level:
            yield v
        nxt = []
        for v in level:
            head = v * base
            if head >= limit:
                break
            for d in ds:
                w = head + d
                if w < limit:
                    nxt.append(w)
        level = nxt


def thin_side(cell, r):
    best = None
    for b, allowed in cell:
        d = design_dim(b, allowed, r)
        if best is None or d < best[0]:
            best = (d, b, allowed)
    return best[1], best[2]


# VERBS


def verb_blocks():
    for name, cell in CELLS:
        r, exps, m, blocks = block_digits(cell)
        assert blocks == bottom_block(cell), name
        shown = list(blocks) if len(blocks) <= 16 else list(blocks[:16]) + ["..."]
        print(
            "%s  bases %s  root %d  exponents %s  M %d  card A %d  A %s"
            % (name, [b for b, _ in cell], r, exps, m, len(blocks), shown)
        )
    r, exps, m = frame(SHARP)
    blocks = bottom_block(SHARP)
    keep = set(SHARP[0][1])
    bad = [n for n in walk(r**m, blocks, 10**6) if not in_design(n, SHARP[0][0], keep)]
    print(
        "sharpness  bases %s  exponents %s  M %d  A %s  in the collapse and not in the joint set %s"
        % ([b for b, _ in SHARP], exps, m, list(blocks), bad[:4])
    )


def verb_dim():
    for name, cell in CELLS + [("I", IRRATIONAL)]:
        r, _exps, m, blocks = block_digits(cell)
        text, value = dim_text(len(blocks), r, m)
        dims = [design_dim(b, allowed, r) for b, allowed in cell]
        if None in dims:
            print(
                "%s  bases %s  root %d  M %d  card A %d  exact dim %s = %.6f  card A is not a power of the root, no rational part and no budget"
                % (name, [b for b, _ in cell], r, m, len(blocks), text, value)
            )
            continue
        raw = sum(dims) - (len(cell) - 1)
        budget = max(Fraction(0), raw)
        print(
            "%s  exact dim %s = %.6f  parts %s  raw sum %s  naive budget max(0, raw) %s = %.6f  gap %s"
            % (
                name,
                text,
                value,
                " + ".join(str(d) for d in dims),
                raw,
                budget,
                float(budget),
                Fraction(text) - budget,
            )
        )


def verb_check():
    for name, cell in CELLS:
        r, _exps, m, blocks = block_digits(cell)
        big = r**m
        start = time.time()
        keeps = [(b, set(a)) for b, a in cell]
        inside = 0
        for n in walk(big, blocks, HEIGHT):
            assert all(in_design(n, b, keep) for b, keep in keeps), (name, n)
            inside += 1
        base, allowed = thin_side(cell, r)
        others = [(b, set(a)) for b, a in cell if b != base]
        marks = []
        power = big
        while power <= HEIGHT:
            marks.append(power)
            power *= big
        tally = [0] * len(marks)
        joint = 0
        for n in walk(base, allowed, HEIGHT):
            if all(in_design(n, b, keep) for b, keep in others):
                joint += 1
                for i, cut in enumerate(marks):
                    if n < cut:
                        tally[i] += 1
        law = [len(blocks) ** (i + 1) for i in range(len(marks))]
        seen = [t + 1 for t in tally]
        assert inside == joint, (name, inside, joint)
        assert seen == law, (name, seen, law)
        print(
            "%s  height %d  collapse inside the joint set %d  joint count %d  counts %s at %s  %.1f s"
            % (name, HEIGHT, inside, joint, law, marks, time.time() - start)
        )
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    print("peak resident %.0f MB" % peak)


def main():
    verbs = sys.argv[1:] or ["blocks", "dim", "check"]
    table = {"blocks": verb_blocks, "dim": verb_dim, "check": verb_check}
    for verb in verbs:
        print("-- %s" % verb)
        table[verb]()


main()
