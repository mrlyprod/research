from fractions import Fraction

CODES = list(range(1, 16))
FRAME = ("gamma", "h", "v", "phi")

# GEOMETRY


def tile(c):
    return [(c >> 0) & 3, (c >> 2) & 3]


def draw(word):
    rows, n = [1], 1
    for c in word:
        t = tile(c)
        out = []
        for r in rows:
            for p in (0, 1):
                acc = 0
                for j in range(n):
                    if (r >> j) & 1:
                        acc |= t[p] << (2 * j)
                out.append(acc)
        rows, n = out, 2 * n
    return rows, n


def fill(rows):
    return sum(r.bit_count() for r in rows)


def hruns(rows):
    return sum((r & ~(r << 1)).bit_count() for r in rows)


def vruns(rows):
    total, prev = 0, 0
    for r in rows:
        total += (r & ~prev).bit_count()
        prev = r
    return total


def comp(rows, n):
    seen = [0] * len(rows)
    total = 0
    for i in range(len(rows)):
        free = rows[i] & ~seen[i]
        while free:
            b = free & -free
            total += 1
            stack = [(i, b.bit_length() - 1)]
            seen[i] |= b
            while stack:
                a, j = stack.pop()
                for da, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    x, y = a + da, j + dj
                    if 0 <= x < len(rows) and 0 <= y < n:
                        m = 1 << y
                        if (rows[x] & m) and not (seen[x] & m):
                            seen[x] |= m
                            stack.append((x, y))
            free = rows[i] & ~seen[i]
    return total


CACHE = {}


def obs(word):
    if word not in CACHE:
        rows, n = draw(word)
        CACHE[word] = (comp(rows, n), hruns(rows), vruns(rows), fill(rows))
    return CACHE[word]


# EXACT LINEAR ALGEBRA


def solve(basis, target):
    m, ncol = len(basis), len(target)
    aug = [[basis[i][j] for i in range(m)] + [target[j]] for j in range(ncol)]
    piv, top = [], 0
    for col in range(m):
        p = next((k for k in range(top, ncol) if aug[k][col] != 0), None)
        if p is None:
            continue
        aug[top], aug[p] = aug[p], aug[top]
        d = aug[top][col]
        aug[top] = [x / d for x in aug[top]]
        for k in range(ncol):
            if k != top and aug[k][col] != 0:
                f = aug[k][col]
                aug[k] = [a - f * b for a, b in zip(aug[k], aug[top])]
        piv.append((top, col))
        top += 1
    if any(aug[k][m] != 0 for k in range(top, ncol)):
        return None
    x = [Fraction(0)] * m
    for r, c in piv:
        x[c] = aug[r][m]
    return x


def mul(A, B):
    d = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(d)) for j in range(d)] for i in range(d)
    ]


def det(A):
    d = len(A)
    if d == 1:
        return A[0][0]
    return sum(
        (-1) ** j * A[0][j] * det([row[:j] + row[j + 1 :] for row in A[1:]])
        for j in range(d)
    )


def eye(d):
    return [[1 if i == j else 0 for j in range(d)] for i in range(d)]


def apply(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def floor_root(num, den, L, digits=9):
    scale = 10**digits
    lo = int((num / den) ** (1.0 / L) * scale) - 4
    while (lo + 1) ** L * den <= num * scale**L:
        lo += 1
    while lo**L * den > num * scale**L:
        lo -= 1
    return "%d.%0*d" % (lo // scale, digits, lo % scale)


def ceil_root(num, den, L, digits=9):
    scale = 10**digits
    hi = int((num / den) ** (1.0 / L) * scale) + 4
    while (hi - 1) ** L * den >= num * scale**L:
        hi -= 1
    while hi**L * den < num * scale**L:
        hi += 1
    return "%d.%0*d" % (hi // scale, digits, hi % scale)


def words(alphabet, length):
    if length == 0:
        yield ()
        return
    for w in words(alphabet, length - 1):
        for c in alphabet:
            yield w + (c,)


def seeded(alphabet, count, lo, hi):
    state, out = 20260906, []
    for _ in range(count):
        state = (1103515245 * state + 12345) % (1 << 31)
        length = lo + state % (hi - lo + 1)
        w = []
        for _ in range(length):
            state = (1103515245 * state + 12345) % (1 << 31)
            w.append(alphabet[state % len(alphabet)])
        out.append(tuple(w))
    return out


# THE REPRESENTATION

SUFFIXES = [()] + [(c,) for c in CODES] + [(a, b) for a in CODES for b in CODES]


def row_of(u):
    return [Fraction(obs(u + s)[0]) for s in SUFFIXES]


def build():
    basis, rows, mats = [()], [row_of(())], {c: [] for c in CODES}
    i = 0
    while i < len(basis):
        for c in CODES:
            r = row_of(basis[i] + (c,))
            x = solve(rows, r)
            if x is None:
                basis.append(basis[i] + (c,))
                rows.append(r)
                x = solve(rows, r)
            mats[c].append(x)
        i += 1
    d = len(basis)
    for c in CODES:
        for line in mats[c]:
            line.extend([Fraction(0)] * (d - len(line)))
        mats[c] = [[int(x) for x in line] for line in mats[c]]
    return basis, mats


print("JSR SCHEDULES")
print()
print("== THE REPRESENTATION ==")
basis, M = build()
print("Hankel rank %d, basis words %s" % (len(basis), [w if w else "e" for w in basis]))
lam = [1] + [0] * (len(basis) - 1)
gamma = [obs(u)[0] for u in basis]
h = [obs(u)[1] for u in basis]
v = [obs(u)[2] for u in basis]
phi = [obs(u)[3] for u in basis]
print("lambda %s, gamma %s, h %s, v %s, phi %s" % (lam, gamma, h, v, phi))
classes = {}
for c in CODES:
    classes.setdefault(tuple(map(tuple, M[c])), []).append(c)
CLASS = sorted(classes.values(), key=lambda g: g[0])
print("class partition %s" % (CLASS,))
assert len(CLASS) == 6
noncomm = sum(
    1
    for i in range(6)
    for j in range(i + 1, 6)
    if mul(M[CLASS[i][0]], M[CLASS[j][0]]) != mul(M[CLASS[j][0]], M[CLASS[i][0]])
)
print("non-commuting class pairs %d of 15" % noncomm)
for g in CLASS:
    print("  M_%s = %s" % (g, M[g[0]]))

# THE OBSERVABLE FRAME

print()
print("== THE OBSERVABLE FRAME ==")
G = [[gamma[i], h[i], v[i], phi[i]] for i in range(4)]
print("frame columns gamma %s, h %s, v %s, phi %s, det %d" % (gamma, h, v, phi, det(G)))
assert h == apply(M[3], gamma) and v == apply(M[5], gamma) and phi == apply(M[1], gamma)
print(
    "h, v and phi are M_3 gamma, M_5 gamma and M_1 gamma, so all four frame vectors are comp read at one appended letter"
)
assert abs(det(G)) == 1
Ginv = None
adj = [
    [
        Fraction(
            (-1) ** (i + j)
            * det([r[:i] + r[i + 1 :] for k, r in enumerate(G) if k != j]),
            det(G),
        )
        for j in range(4)
    ]
    for i in range(4)
]
Ginv = adj
assert mul(Ginv, G) == eye(4)
T = {}
for c in CODES:
    t = mul(Ginv, mul(M[c], G))
    T[c] = [[int(x) for x in row] for row in t]
    assert all(Fraction(T[c][i][j]) == t[i][j] for i in range(4) for j in range(4))
print("frame images, columns of T_c in the frame (gamma, h, v, phi):")
for g in CLASS:
    c = g[0]
    rows, n = draw((c,))
    ob = obs((c,))
    cols = [[T[c][i][j] for i in range(4)] for j in range(4)]
    csum = [sum(col) for col in cols]
    print(
        "  code %-2d k=%d  M_c gamma=%s  M_c h=%s  M_c v=%s  M_c phi=%s  column sums %s"
        % (c, ob[3], cols[0], cols[1], cols[2], cols[3], csum)
    )
    assert csum == [ob[0], ob[1], ob[2], ob[3]]
    assert max(csum) == ob[3]
print(
    "every column sum is (comp, rows, cols, fill) of the letter's own tile, and the largest is the fill"
)

# THE FOUR TRANSFER LAWS, CHECKED AGAINST DRAWN CELLS

print()
print("== THE TRANSFER LAWS ==")
short = [w for L in range(4) for w in words(CODES, L)]
seeds = seeded(CODES, 240, 4, 7)
sweep = short + seeds
bad = 0
for w in sweep:
    y = obs(w)
    for c in CODES:
        z = obs(w + (c,))
        pred = tuple(sum(y[i] * T[c][i][j] for i in range(4)) for j in range(4))
        if pred != z:
            bad += 1
assert bad == 0
print(
    "(comp, H, V, fill) at wc equals (comp, H, V, fill) at w times T_c on %d words, all %d of length at most 3 and %d seeded of length 4 to 7, times 15 letters, %d mismatches"
    % (len(sweep), len(short), len(seeds), bad)
)
rep_bad = 0
for w in sweep:
    P = eye(4)
    for c in w:
        P = mul(P, M[c])
    val = sum(lam[i] * sum(P[i][j] * gamma[j] for j in range(4)) for i in range(4))
    hv = sum(lam[i] * sum(P[i][j] * h[j] for j in range(4)) for i in range(4))
    vv = sum(lam[i] * sum(P[i][j] * v[j] for j in range(4)) for i in range(4))
    fv = sum(lam[i] * sum(P[i][j] * phi[j] for j in range(4)) for i in range(4))
    if (val, hv, vv, fv) != obs(w):
        rep_bad += 1
assert rep_bad == 0
print(
    "lambda M_w applied to gamma, h, v, phi reads comp, H, V, fill on the same %d words, %d of length at most 3 and %d seeded, %d mismatches"
    % (len(sweep), len(short), len(seeds), rep_bad)
)

# TRIANGULARITY AND SPECTRUM

print()
print("== TRIANGULARITY ==")
for g in CLASS:
    c = g[0]
    assert all(T[c][i][j] == 0 for i in range(4) for j in range(4) if j > i)
    assert all(T[c][i][j] >= 0 for i in range(4) for j in range(4))
    diag = [T[c][i][i] for i in range(4)]
    poly = [t for t in range(5)]
    vals = [
        det(
            [
                [t * (1 if i == j else 0) - M[c][i][j] for j in range(4)]
                for i in range(4)
            ]
        )
        for t in poly
    ]
    pred = [(t - diag[0]) * (t - diag[1]) * (t - diag[2]) * (t - diag[3]) for t in poly]
    assert vals == pred
    print(
        "  code %-2d T_c lower triangular, nonnegative, diagonal %s, char poly roots %s, rho = %d = fill"
        % (c, diag, sorted(diag), max(diag))
    )
    assert max(diag) == obs((c,))[3]

# THE CROSS-POLYTOPE CERTIFICATE

print()
print("== THE CERTIFICATE ==")
print(
    "P = conv{+/- gamma, +/- h, +/- v, +/- phi}; residual is k_c minus the frame l1 norm of the image"
)
for g in CLASS:
    c = g[0]
    k = obs((c,))[3]
    res = []
    for j, name in enumerate(FRAME):
        col = [T[c][i][j] for i in range(4)]
        res.append(k - sum(abs(x) for x in col))
    print("  code %-2d k=%d residuals at gamma, h, v, phi: %s" % (c, k, res))
    assert all(r >= 0 for r in res)
print(
    "M_c P is inside k_c P at every code, so the gauge of P, the frame l1 norm N(a gamma + b h + c v + d phi) = |a| + |b| + |c| + |d|, is an extremal norm with ||M_c|| = k_c"
)

# JSR AND LSR OVER EVERY SUBFAMILY

print()
print("== JSR AND LSR ==")
K = {c: obs((c,))[3] for c in CODES}
pairs = [(a, b) for i, a in enumerate(CODES) for b in CODES[i + 1 :]]
print("pairs of distinct letters %d, pairs of distinct classes %d" % (len(pairs), 15))
print("class pair table (a, b, JSR, LSR):")
for i in range(6):
    for j in range(i + 1, 6):
        a, b = CLASS[i][0], CLASS[j][0]
        print("  (%d, %d) JSR %d LSR %d" % (a, b, max(K[a], K[b]), min(K[a], K[b])))
distinct = sum(1 for a, b in pairs if K[a] != K[b])
print(
    "of the %d pairs of distinct letters, %d carry two different fills and %d carry one"
    % (len(pairs), distinct, len(pairs) - distinct)
)
print(
    "JSR = max fill and LSR = min fill on every one of them, each attained by a one-letter word"
)

# THE BRACKET ON TWO NAMED PAIRS

print()
print("== THE BRACKET ==")
LMAX = 16
for pair in [(3, 6), (3, 7)]:
    a, b = pair
    r = max(K[a], K[b])
    cur = {(): eye(4)}
    best_lo, best_word, best_up, best_min = {}, {}, {}, {}
    for L in range(1, LMAX + 1):
        nxt = {}
        for w, P in cur.items():
            for c in pair:
                nxt[w + (c,)] = mul(P, T[c])
        cur = nxt
        lo, arg, up, worst = -1, None, 0, None
        for w, P in cur.items():
            sr = max(P[i][i] for i in range(4))
            if sr > lo:
                lo, arg = sr, w
            if worst is None or sr < worst:
                worst = sr
            nrm = max(sum(P[i][j] for i in range(4)) for j in range(4))
            up = max(up, nrm)
        best_lo[L], best_word[L], best_up[L], best_min[L] = lo, arg, up, worst
        assert lo <= up
        assert lo == r**L
        assert up == r**L
        assert worst == min(K[a], K[b]) ** L
    print("pair {%d, %d}, JSR candidate %d" % (a, b, r))
    for L in (1, 2, 4, 8, 12, LMAX):
        w = "".join(str(x) for x in best_word[L])
        print(
            "  L=%-2d lower %s from word %s, spectral radius %d; upper %s from max frame 1-norm %d"
            % (
                L,
                floor_root(best_lo[L], 1, L),
                w,
                best_lo[L],
                ceil_root(best_up[L], 1, L),
                best_up[L],
            )
        )
    print(
        "  the best word's rate is exactly %d at every length 1 to %d, so it stops improving at length 1"
        % (r, LMAX)
    )
    print(
        "  the worst word's spectral radius is exactly %d^L at every length, so the lower spectral radius is %d"
        % (min(K[a], K[b]), min(K[a], K[b]))
    )
    for k in (2, 4, 6):
        num = K[a] ** k + K[b] ** k
        lift_lo = floor_root(num, 2, k)
        lift_up = ceil_root(num, 1, k)
        print(
            "  Blondel-Nesterov lifting at k=%d on the 2 letters: lower 2^(-1/%d) (%d^%d + %d^%d)^(1/%d) = %s, upper (%d^%d + %d^%d)^(1/%d) = %s"
            % (k, k, K[a], k, K[b], k, k, lift_lo, K[a], k, K[b], k, k, lift_up)
        )
        assert Fraction(lift_lo) <= Fraction(lift_up)
        assert Fraction(lift_lo) <= r <= Fraction(lift_up)
        assert num >= r**k
    print(
        "  the scan alone brackets it at [%s, %s] at L=%d, the lifting alone at [%s, %s] at k=6, and the certificate closes it at %d exactly"
        % (
            floor_root(best_lo[LMAX], 1, LMAX),
            ceil_root(best_up[LMAX], 1, LMAX),
            LMAX,
            floor_root(K[a] ** 6 + K[b] ** 6, 2, 6),
            ceil_root(K[a] ** 6 + K[b] ** 6, 1, 6),
            r,
        )
    )

# THE TELESCOPE BLOCK

print()
print("== THE TELESCOPE ==")
A = [[M[3][i][j] for j in range(2)] for i in range(2)]
B = [[M[6][i][j] for j in range(2)] for i in range(2)]
assert all(M[3][i][j] == 0 for i in range(4) for j in (2, 3))
assert all(M[6][i][j] == 0 for i in range(4) for j in (2, 3))
print(
    "on {3, 6} the last two columns of both matrices vanish, leading blocks A = %s, B = %s"
    % (A, B)
)
p = [1, 2]
assert apply(A, p) == [2 * x for x in p]
assert apply(B, [1, 1]) == [2 * x for x in p] and apply(B, p) == [2 * x for x in p]
print(
    "A p = 2 p and B has range span(p) at p = %s, so span(p) is invariant under both and the pair is reducible"
    % p
)
print(
    "in the basis (gamma, p) both blocks are lower triangular, A = diag(1, 2) and B = [[0, 0], [2, 2]]"
)
for L in (1, 2, 3, 6, 10):
    vals = set()
    for w in words((3, 6), L):
        P = eye(4)
        for c in w:
            P = mul(P, T[c])
        vals.add(max(P[i][i] for i in range(4)))
    print(
        "  every word of length %d has spectral radius %s, so its rate is exactly 2"
        % (L, sorted(vals))
    )
    assert vals == {2**L}

# THE DEAD ROUTES

print()
print("== WHAT DIES ==")
for g in CLASS:
    c = g[0]
    assert apply(M[c], phi) == [K[c] * x for x in phi]
print(
    "phi = %s is a right eigenvector of every M_c with eigenvalue the fill, so span(phi) is a common invariant line and no subfamily is irreducible"
    % phi
)
orbit = {tuple(phi)}
for _ in range(4):
    orbit |= {
        tuple(Fraction(t, K[c]) for t in apply(M[c], list(x)))
        for x in orbit
        for c in (3, 7)
    }
    orbit = {tuple(int(t) if t.denominator == 1 else t for t in x) for x in orbit}
assert orbit == {tuple(phi)}
print(
    "the polytope algorithm started at the leading eigenvector of the length-1 spectrum maximizing product, each image normalised by the letter fill k_c, closes on %d vertex up to sign and spans dimension 1 of 4"
    % len(orbit)
)
norms = []
for L in (1, 2, 4, 8, 16, 32):
    P = eye(4)
    for _ in range(L):
        P = mul(P, M[3])
    norms.append(max(abs(P[i][j]) for i in range(4) for j in range(4)))
print(
    "max |entry(M_3^L)| at L = 1, 2, 4, 8, 16, 32 reads %s, matching 2^(L+2) - 2"
    % norms
)
assert norms == [2 ** (L + 2) - 2 for L in (1, 2, 4, 8, 16, 32)]
print(
    "comp(A_(3^L)) = %s at the same lengths, so the JSR rate log 2 is not the component rate 0"
    % [obs((3,) * L)[0] for L in (1, 2, 4, 8)]
)
print()
print("every assertion passed")
