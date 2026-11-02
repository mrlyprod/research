import hashlib
import time
from math import gcd

FAREY_N = 55
Q = 9973
M = 2 * Q
PROBES = 48
PROBE_NS = [55, 5555, 10**18]
SWEEP_NS = [1, 2, 55, 5555, 19945, 19946, 19947, 40001]
R = 512
D = 2 * R
IMAGE_N = 55
ODDS = list(range(1, IMAGE_N + 1, 2))
K = len(ODDS)


def farey_nodes(n_max):
    out = []
    for b in range(1, n_max + 1):
        for a in range(1, b + 1):
            if gcd(a, b) == 1:
                out.append((a, b))
    return out


def farey_by_form(n_max):
    return ["%d/%d:%d" % (a, b, n_max // b) for (a, b) in farey_nodes(n_max)]


def farey_by_membership(n_max):
    scales = []
    for n in range(1, n_max + 1):
        scales.append({(k // gcd(k, n), n // gcd(k, n)) for k in range(1, n + 1)})
    lines = []
    for a, b in farey_nodes(n_max):
        lines.append("%d/%d:%d" % (a, b, sum(1 for s in scales if (a, b) in s)))
    return lines


def digest(lines):
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def brightness_sum(lines):
    return sum(int(t.split(":")[1]) for t in lines)


def probe_points():
    return [((137 * i + 61) % Q, (211 * i + 97) % Q) for i in range(PROBES)]


def bad_residues(a, b):
    return [r for r in range(1, M, 2) if (r * a) % M >= Q and (r * b) % M >= Q]


def brightness_closed(bad, n_max):
    dark = 0
    for r in bad:
        if r <= n_max:
            dark += (n_max - r) // M + 1
    return (n_max + 1) // 2 - dark


def good_residues(a, b):
    return [r for r in range(1, M, 2) if (r * a) % M < Q or (r * b) % M < Q]


def brightness_direct(good, n_max):
    lit = 0
    for r in good:
        if r <= n_max:
            lit += (n_max - r) // M + 1
    return lit


def brightness_literal(a, b, n_max):
    c = 0
    for n in range(1, n_max + 1, 2):
        if not (((n * a) // Q) % 2 == 1 and ((n * b) // Q) % 2 == 1):
            c += 1
    return c


def column_masks():
    masks = []
    for j in range(R):
        u = 2 * j + 1
        m = 0
        for k, n in enumerate(ODDS):
            if ((n * u) // D) % 2 == 1:
                m |= 1 << k
        masks.append(m)
    return masks


def render_masks():
    masks = column_masks()
    buf = bytearray(R * R)
    for i in range(R):
        my = masks[i]
        row = i * R
        for j in range(R):
            buf[row + j] = K - (my & masks[j]).bit_count()
    return buf


def render_layers():
    xs = [2 * j + 1 for j in range(R)]
    img = [[0] * R for _ in range(R)]
    for n in ODDS:
        odd_col = [((n * v) // D) % 2 == 1 for v in xs]
        for i in range(R):
            row = img[i]
            if ((n * xs[i]) // D) % 2 == 1:
                for j in range(R):
                    if not odd_col[j]:
                        row[j] += 1
            else:
                for j in range(R):
                    row[j] += 1
    buf = bytearray()
    for i in range(R):
        buf.extend(bytes(img[i]))
    return buf


def render_pixels():
    buf = bytearray()
    for i in range(R):
        y = 2 * i + 1
        for j in range(R):
            x = 2 * j + 1
            c = 0
            for n in ODDS:
                if not (((n * x) // D) % 2 == 1 and ((n * y) // D) % 2 == 1):
                    c += 1
            buf.append(c)
    return buf


def main():
    form = farey_by_form(FAREY_N)
    member = farey_by_membership(FAREY_N)
    print("farey N", FAREY_N)
    print("farey nodes", len(form))
    print("farey brightness sum", brightness_sum(form))
    print("farey N(N+1)/2", FAREY_N * (FAREY_N + 1) // 2)
    print("farey digest closed form", digest(form))
    print("farey digest membership", digest(member))
    print("farey digests equal", digest(form) == digest(member))

    pts = probe_points()
    t0 = time.perf_counter()
    tables = [bad_residues(a, b) for (a, b) in pts]
    t_tables = time.perf_counter() - t0
    print("probe count", len(pts))
    print("probe modulus 2q", M)
    print("probe residue tables seconds", round(t_tables, 4))

    values = {}
    for n_max in PROBE_NS:
        t0 = time.perf_counter()
        values[n_max] = [brightness_closed(t, n_max) for t in tables]
        dt = time.perf_counter() - t0
        print("closed form N", n_max, "seconds", round(dt, 4))
        print("closed form N", n_max, "first four", *values[n_max][:4])

    t0 = time.perf_counter()
    twins = [good_residues(a, b) for (a, b) in pts]
    print("direct count residue tables seconds", round(time.perf_counter() - t0, 4))
    for n_max in PROBE_NS:
        t0 = time.perf_counter()
        other = [brightness_direct(t, n_max) for t in twins]
        dt = time.perf_counter() - t0
        print("direct count N", n_max, "seconds", round(dt, 4))
        print("direct count N", n_max, "equals closed form", other == values[n_max])

    for n_max in (55, 5555):
        agree = sum(
            1
            for i, (a, b) in enumerate(pts)
            if values[n_max][i] == brightness_literal(a, b, n_max)
        )
        print("probes closed equal literal N", n_max, "%d/%d" % (agree, len(pts)))

    bad = 0
    for i, (a, b) in enumerate(pts):
        for n_max in SWEEP_NS:
            if brightness_closed(tables[i], n_max) != brightness_literal(a, b, n_max):
                bad += 1
    print("sweep N values", *SWEEP_NS)
    print("sweep comparisons", len(pts) * len(SWEEP_NS))
    print("sweep mismatches", bad)

    a = render_masks()
    b = render_layers()
    c = render_pixels()
    print("image side", R)
    print("image layers", K)
    print("image sha256 masks", hashlib.sha256(bytes(a)).hexdigest())
    print("image sha256 layers", hashlib.sha256(bytes(b)).hexdigest())
    print("image sha256 pixels", hashlib.sha256(bytes(c)).hexdigest())
    print("image three routes equal", a == b == c)


main()
