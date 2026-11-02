import sys

# BITMASK POLYNOMIALS OVER F2


def polymul2(a, b):
    r = 0
    while b:
        lb = b & -b
        r ^= a * lb
        b ^= lb
    return r


def pow2mask(m):
    r = 1
    s = 1
    while m:
        if m & 1:
            r ^= r << s
        s <<= 1
        m >>= 1
    return r


def divmod2(a, b):
    db = b.bit_length() - 1
    q = 0
    while a and a.bit_length() - 1 >= db:
        sh = a.bit_length() - 1 - db
        q ^= 1 << sh
        a ^= b << sh
    return q, a


def bits(v):
    out = []
    while v:
        lb = v & -v
        out.append(lb.bit_length() - 1)
        v ^= lb
    return out


def echelon(vecs):
    B = []
    for v in vecs:
        for u in B:
            v = min(v, v ^ u)
        if v:
            B.append(v)
            B.sort(reverse=True)
    return B


def reduce_by(v, B):
    for u in B:
        v = min(v, v ^ u)
    return v


def pivots(vecs):
    piv = {}
    for v in vecs:
        while v:
            t = v.bit_length() - 1
            if t in piv:
                v ^= piv[t]
            else:
                piv[t] = v
                break
    return piv


def residue(v, piv, keys):
    for t in keys:
        if (v >> t) & 1:
            v ^= piv[t]
    return v


def spans(v, piv):
    while v:
        t = v.bit_length() - 1
        if t not in piv:
            return False
        v ^= piv[t]
    return True


def nullspace_right(rows, n):
    piv = {}
    for r in rows:
        for p in sorted(piv, reverse=True):
            if (r >> p) & 1:
                r ^= piv[p]
        if r:
            q = r.bit_length() - 1
            for pp in list(piv):
                if (piv[pp] >> q) & 1:
                    piv[pp] ^= r
            piv[q] = r
    out = []
    for f in [c for c in range(n) if c not in piv]:
        v = 1 << f
        for p, r in piv.items():
            if (r >> f) & 1:
                v |= 1 << p
        out.append(v)
    return out


# JACOBSTHAL AND THE SLOT


def jac(k):
    return 0 if k < 0 else (2**k - (-1) ** k) // 3


def slot(R):
    b = 0
    while (1 << b) < 3 * R - 1:
        b += 1
    g = abs(2 * R - (1 << (b - 1)) - 1)
    e = 1
    while jac(e) < (g + 1) // 2:
        e += 1
    return b, g, e, b - 1 - e


def window(R, b):
    A = 1 << b
    i0 = max(2, 4 * R - A)
    hi = (6 * R + 2 - A) // 3
    hi -= hi % 2
    return i0, (hi - i0) // 2


def fibpoly(t):
    if t == 0:
        return 1
    a, b = 0b11, 1
    for _ in range(t - 1):
        a, b = polymul2(a, 2) ^ b, a
    return a


def law_e(D):
    R = (D - 1) // 2
    b, _, e, k = slot(R)
    _, K = window(R, b)
    tt = (jac(k) - 1) // 2
    C = K - (2 * jac(e - 1) if k % 2 == 0 else 0)
    w = C - tt * (1 << e)
    ce = 2 * jac(e - 2) - 1 if e >= 3 else 1
    m = max(0, 2 * w - ce)
    g = 0
    for a in bits(fibpoly(tt)):
        g ^= 1 << (a * (1 << e))
    return K, C, g << m, min(w + 1, ce + 1 - w)


def tent(R):
    b, gslot, e, _ = slot(R)
    N = jac(e) - jac(e - 1)
    u = (gslot + 1) // 2 - jac(e - 1) - 1
    p = u if R > (1 << (b - 2)) else N - 1 - u
    return p, min(p, N - 1 - p)


def reach_law(R):
    b, _, e, k = slot(R)
    p, w = tent(R)
    if e == 1 and k % 2:
        return 1 if R == 1 << (b - 2) else (3 if R == 3 else 5)
    return 3 * w + (2 if e % 2 == 0 else 0) + ((1 + p % 2) if k % 2 else 0)


# THE MOD-4 SYMBOL AND THE KERNEL FAMILY

FW = 20
FM = (1 << FW) - 1


def symbol4(D):
    out = [0] * (2 * D + 1)
    c = 1
    for k in range(D):
        v = c % 4
        if v:
            out[2 * k] = (out[2 * k] + v) % 4
            out[2 * k + 1] = (out[2 * k + 1] + v * D) % 4
            out[2 * k + 2] = (out[2 * k + 2] + v) % 4
        c = c * (D - 1 - k) // (k + 1)
    return out


def symbol4pack(D):
    p = 0
    for i, v in enumerate(symbol4(D)):
        if v:
            p |= v << (FW * i)
    return p


def family(D):
    R = (D - 1) // 2
    b, _, _, _ = slot(R)
    i0, K = window(R, b)
    A = 1 << b
    G = polymul2(pow2mask(4 * R), 0b111)
    out = []
    for j in range(K + 1):
        i = i0 + 2 * j
        q = 1
        for a in bits(i):
            q = polymul2(q, 1 ^ (1 << (3 << a)))
        f = polymul2(1 << ((6 * R + 2 - 3 * i - A) // 2), q)
        f ^= f << A
        h, r = divmod2(f, G)
        assert r == 0, "family element off the ideal at D=%d" % D
        out.append(h)
    return out, K


def colblock(R, G, wide):
    dec = [0, 0, 0]
    for a in bits(G):
        dec[a % 3] |= 1 << (a // 3)
    m = (1 << (R + 1)) - 1

    def T(a):
        c = 1 - R + a
        r = c % 3
        q = (c - r) // 3
        v = dec[r]
        return ((v >> q) if q >= 0 else (v << -q)) & m

    out = [T(0)]
    for j in range(1, wide):
        out.append(T(-j) ^ T(j))
    return out


def columns(D):
    R = (D - 1) // 2
    return colblock(R, polymul2(pow2mask(4 * R), 0b111), R + 1)


def obstruction(D, X, Ppack):
    R = (D - 1) // 2
    out = []
    for h in X:
        hp = 0
        for a in bits(h):
            hp |= 1 << (FW * a)
        pr = hp * Ppack
        o = 0
        for nu in range(R + 1):
            s = (pr >> (FW * (3 * nu + 1))) & FM
            assert s % 2 == 0, "family element not in the mod-2 kernel at D=%d" % D
            if s % 4 == 2:
                o |= 1 << nu
        out.append(o)
    return out


# THE LAYER-2 WINDOW


def layer2(D):
    X, K = family(D)
    cols = columns(D)
    piv = pivots(cols)
    keys = sorted(piv, reverse=True)
    raw = obstruction(D, X, symbol4pack(D))
    red = [residue(o, piv, keys) for o in raw]
    nb = max((o.bit_length() for o in red), default=0)
    rows = []
    for nu in range(nb):
        r = 0
        for j in range(K + 1):
            if (red[j] >> nu) & 1:
                r |= 1 << j
        if r:
            rows.append(r)
    V2 = echelon(nullspace_right(rows, K + 1))
    return K, V2, raw, cols, X, piv


def shift_law(D, X, K, raw, piv):
    R = (D - 1) // 2
    G = polymul2(pow2mask(4 * R), 0b111)
    box = (1 << (2 * R + 1)) - 1
    for j in range(K):
        H = X[j]
        hi, lo = H << 3, H >> 3
        both = hi & lo
        assert X[j + 1] == (hi ^ lo), "family shift breaks at D=%d j=%d" % (D, j)
        Z = H ^ both
        assert (Z | box) == box, "Z leaves the box at D=%d j=%d" % (D, j)
        assert Z == rev_box(Z, R), "Z not palindromic at D=%d j=%d" % (D, j)
        ZG = polymul2(Z, G)
        az = 0
        for nu in range(R + 1):
            if (ZG >> (3 * nu + 1)) & 1:
                az |= 1 << nu
        assert spans(az, piv), "A(Z) off the image at D=%d j=%d" % (D, j)
        assert raw[j + 1] == fold_shift(raw[j], R) ^ az, (
            "Lambda law breaks at D=%d j=%d" % (D, j)
        )
    return True


def rev_box(v, R):
    out = 0
    for a in bits(v):
        out |= 1 << (2 * R - a)
    return out


def fold_shift(f, R):
    out = 0
    for nu in range(R + 1):
        a = f >> (nu - 1) & 1 if nu >= 1 else 0
        b = (
            (f >> (nu + 1)) & 1
            if nu + 1 <= R
            else ((f >> (R - 1)) & 1 if nu == R else 0)
        )
        if a ^ b:
            out |= 1 << nu
    return out


def corrector_index(cols, target):
    piv = {}
    for J, c in enumerate(cols):
        v = c
        while v:
            t = v.bit_length() - 1
            if t in piv:
                v ^= piv[t]
            else:
                piv[t] = v
                break
        if spans(target, piv):
            return J
    return None


# THE SWEEP


def sweep(lo, hi):
    n = bw = bg = bc = bm = br = bt = bf = ne = fr = 0
    cap, flo, tie = [], [], []
    r4 = [0] * 4
    r8 = [0] * 8
    for D in range(lo, hi + 1, 2):
        R = (D - 1) // 2
        b, _, e, k = slot(R)
        K, V2, raw, cols, X, piv = layer2(D)
        assert V2, "empty layer-2 window at D=%d" % D
        assert X[K].bit_length() - 1 <= 2 * R, "family top leaves the box at D=%d" % D
        gen = min(V2)
        dg = gen.bit_length() - 1
        C = max(v.bit_length() - 1 for v in V2)
        n += 1
        r4[R % 4] += 1
        r8[R % 8] += 1
        if echelon([polymul2(gen, 1 << s) for s in range(C - dg + 1)]) != V2:
            bw += 1
            print("window fails at D=%d" % D)
        Kp, Cp, gp, L2p = law_e(D)
        if (Kp, Cp, gp, L2p) != (K, C, gen, len(V2)):
            print(
                "law E fails at D=%d: %r against %r"
                % (D, (Kp, Cp, gp, L2p), (K, C, gen, len(V2)))
            )
            if gp != gen:
                bg += 1
            if Cp != C or L2p != len(V2):
                bc += 1
        shift_law(D, X, K, raw, piv)
        om = 0
        for j in bits(gen):
            om ^= raw[j]
        J = corrector_index(cols, om)
        if k % 2 or e == 1:
            assert C == K, "upper half not free at D=%d" % D
            fr += 1
        else:
            assert K - C == 2 * jac(e - 1), "ceiling deficit is not 2J(e-1) at D=%d" % D
        if tent(R)[1] != C - dg:
            bt += 1
            print("tent identity fails at D=%d" % D)
        rl = reach_law(R)
        if rl != R - J:
            br += 1
            print("reach law fails at D=%d: %d against %d" % (D, rl, R - J))
        a, bq = K - dg, rl // 3
        if e == 1 and k % 2 and R > (1 << (b - 2)):
            ne += 1
            if (bq, C - dg, a) != (1, 0, 0):
                bf += 1
                print("escaping row misread at D=%d" % D)
        elif bq != C - dg + (1 if (k % 2 and e % 2 == 0) else 0):
            bf += 1
            print("floor identity fails at D=%d" % D)
        if C - dg != min(a, bq):
            bm += 1
            print("corrector law fails at D=%d" % D)
        elif a < bq:
            cap.append(D)
        elif bq < a:
            flo.append(D)
        else:
            tie.append(D)
    print("odd D = %d..%d: %d rows" % (lo, hi, n))
    print("R mod 4 classes %r, R mod 8 classes %r" % (r4, r8))
    print("window structure V2 = g F_2[z]_(<= C - deg g): %d/%d" % (n - bw, n))
    print("Law E generator z^m c_t(z^(2^e)): %d/%d" % (n - bg, n))
    print("ceiling C = K - 2J(e-1)[k even] and L_2: %d/%d" % (n - bc, n))
    print(
        "the family shift law H^(j+1) = psi H^(j) - 2 Z_j and ob(X_(j+1)) = Lambda ob(X_j) + A(Z_j): %d/%d"
        % (n, n)
    )
    print("the slot tent w = min(p, N - 1 - p) = C - deg g: %d/%d" % (n - bt, n))
    print(
        "the reach law reach = 3w + 2[e even] + [k odd](1 + p mod 2), D = 4^m + 3 apart: %d/%d"
        % (n - br, n)
    )
    print(
        "floor(reach/3) = C - deg g + [k odd and e even] off the %d rows D = 4^m + 3, where it reads 1 against C - deg g = K - deg g = 0: %d/%d"
        % (ne, n - bf, n)
    )
    print(
        "corrector law C - deg g = min(K - deg g, floor(reach/3)) from the closed form: %d/%d"
        % (n - bm, n)
    )
    print(
        "upper half free where C = K, that is k odd or e = 1, else deficit K - C = 2J(e-1): free %d, open %d"
        % (fr, n - fr)
    )
    print(
        "branches: floor strict %d, K cap strict %d, tie %d"
        % (len(flo), len(cap), len(tie))
    )
    print(
        "floor-strict rows are exactly the C < K rows: %s"
        % (
            sorted(flo)
            == sorted(d for d in range(lo, hi + 1, 2) if law_e(d)[1] < law_e(d)[0])
        )
    )


def unboxed(lo, hi):
    n = bad = co = 0
    for D in range(lo, hi + 1, 2):
        R = (D - 1) // 2
        G = polymul2(pow2mask(4 * R), 0b111)
        piv = pivots(colblock(R, G, 4 * R + 12))
        X, _ = family(D)
        n += 1
        if len(piv) != R:
            co += 1
        if any(not spans(o, piv) for o in obstruction(D, X, symbol4pack(D))):
            bad += 1
    print("unboxed corrector image has corank exactly 1: %d/%d" % (n - co, n))
    print(
        "every family obstruction meets that image, so the unboxed layer-2 window is all of the mod-2 kernel: %d/%d"
        % (n - bad, n)
    )


if __name__ == "__main__":
    a = sys.argv[1:]
    lo = int(a[0]) if a else 5
    hi = int(a[1]) if len(a) > 1 else 1601
    sweep(lo, hi)
    unboxed(5, min(hi, 601))
