import math
import sys
import time
from fractions import Fraction

import numpy as np
from mpmath import iv

iv.prec = 128
BOX = iv.mpc(iv.mpf([-1, 1]), iv.mpf([-1, 1]))

# EXACT PRINTING


def to_frac(t):
    sign, man, exp, bc = t
    if man == 0 and bc < 0:
        return None
    v = Fraction(man) * (Fraction(2) ** exp)
    return -v if sign else v


def lo(x):
    return to_frac(x._mpi_[0])


def hi(x):
    return to_frac(x._mpi_[1])


def dec(n, d):
    neg = n < 0
    n = abs(n)
    ip, fp = divmod(n, 10**d)
    return ("-" if neg else "") + f"{ip}.{fp:0{d}d}"


def floor_str(fr, d):
    return dec((fr.numerator * 10**d) // fr.denominator, d)


def ceil_str(fr, d):
    return dec(-((-fr.numerator * 10**d) // fr.denominator), d)


def fmt_r(x, d):
    return f"[{floor_str(lo(x), d)}, {ceil_str(hi(x), d)}]"


def fmt_c(z, d):
    return fmt_r(z.real, d) + " + i " + fmt_r(z.imag, d)


def fmt_up(x, d=3):
    v = hi(x)
    if v <= 0:
        return f"{dec(0, d)}e+00"
    e = 0
    while v >= 10:
        v /= 10
        e += 1
    while v < 1:
        v *= 10
        e -= 1
    m = -((-v.numerator * 10**d) // v.denominator)
    if m >= 10 ** (d + 1):
        m //= 10
        e += 1
    return f"{dec(m, d)}e{e:+03d}"


def dist_zero(z, d):
    dx = max(Fraction(0), lo(z.real), -hi(z.real))
    dy = max(Fraction(0), lo(z.imag), -hi(z.imag))
    r2 = dx * dx + dy * dy
    n = math.isqrt((r2.numerator * 10 ** (2 * d)) // r2.denominator)
    return dec(n, d), r2 > 0


def meets(z1, z2):
    ok_re = max(lo(z1.real), lo(z2.real)) <= min(hi(z1.real), hi(z2.real))
    ok_im = max(lo(z1.imag), lo(z2.imag)) <= min(hi(z1.imag), hi(z2.imag))
    return ok_re and ok_im


def contains(z, cre, cim):
    return lo(z.real) <= cre <= hi(z.real) and lo(z.imag) <= cim <= hi(z.imag)


def width(z):
    return max(float(hi(z.real) - lo(z.real)), float(hi(z.imag) - lo(z.imag)))


def rat(p, q):
    return iv.mpf(p) / iv.mpf(q)


# DESIGN


class Design:
    def __init__(self, q, F):
        assert 0 in F and len(F) >= 2
        self.q = q
        self.F = sorted(F)
        self.N = len(F)
        self.F1 = [a for a in self.F if a > 0]
        self.N1 = len(self.F1)
        self.amin = min(self.F1)
        self.amax = max(self.F1)
        self.logq = iv.log(q)
        self.s0 = iv.log(self.N) / self.logq
        self.name = f"base {q} digits {{{','.join(map(str, self.F))}}}"

    def gamma(self, l):
        return sum(a**l for a in self.F)

    def s(self, k):
        return iv.mpc(self.s0, 2 * iv.pi * k / self.logq)

    def lev0(self, w):
        acc = iv.mpc(0)
        for a in self.F1:
            acc = acc + (iv.mpc(1) if a == 1 else iv.exp(-w * iv.log(a)))
        return acc

    def amin_pow(self, i):
        return (iv.mpf(self.amin) ** (-(self.s0 + i))).b

    def j_start(self):
        j0 = 0
        while 3 * self.amax > self.amin * self.q ** (j0 + 1):
            j0 += 1
        return j0

    def level(self, j):
        lev = list(self.F1)
        for _ in range(j):
            lev = [self.q * n + a for n in lev for a in self.F]
        return lev


# TAILS


def poch_row(a, T):
    pr = [iv.mpf(1)]
    for l in range(T + 1):
        pr.append(pr[-1] * (a + l) / (l + 1))
    return pr


def hyp_tail(a, rho, T, pr=None):
    if pr is None:
        pr = poch_row(a, T)
    part = iv.mpf(0)
    rp = iv.mpf(1)
    for l in range(T + 1):
        part = part + pr[l] * rp
        rp = rp * rho
    return ((1 - rho) ** (-a) - part).b


def choose_I(a, rho, eps):
    T = 8
    while True:
        if hi(hyp_tail(a, rho, T)) < eps:
            return T + 2
        T += 2


# ENGINE ONE: THE LEVEL RECURSION


def level_limit(D, k, L, I):
    s = D.s(k)
    q, N = D.q, D.N
    coeff, poch, absw = [], [], []
    for i in range(I + 1):
        w = s + i
        aw = abs(w).b
        T = I - i
        row, P = [], iv.mpc(1)
        for l in range(T + 1):
            row.append(P * rat(((-1) ** l) * D.gamma(l), q**l))
            P = P * (w + l) / (l + 1)
        coeff.append(row)
        poch.append(poch_row(aw, T))
        absw.append(aw)
    j0 = D.j_start()
    lev = D.level(j0)
    z = [iv.exp(-s * iv.log(n)) if n > 1 else iv.mpc(1) for n in lev]
    pw = [iv.mpf(1) for n in lev]
    V = []
    for i in range(I + 1):
        acc = iv.mpc(0)
        for t, n in enumerate(lev):
            acc = acc + z[t] * pw[t]
            pw[t] = pw[t] / n
        V.append(acc)
    mult = [rat(1, N * q**i) for i in range(I + 1)]
    apow = [D.amin_pow(i) for i in range(I + 1)]
    for j in range(j0, L):
        rho = rat(D.amax, D.amin * q ** (j + 1))
        newV = []
        for i in range(I + 1):
            T = I - i
            row = coeff[i]
            acc = iv.mpc(0)
            for l in range(T + 1):
                acc = acc + row[l] * V[i + l]
            tail = hyp_tail(absw[i], rho, T, poch[i])
            tau = (N * D.N1 * apow[i] * rat(1, q ** (j * i)) * tail).b
            newV.append(mult[i] * (acc + tau * BOX))
        V = newV
    rhoL = rat(D.amax, D.amin * q ** (L + 1))
    a = absw[0]
    E = (D.N1 * apow[0] * a * (1 - rhoL) ** (-(a + 1)) * rhoL * rat(q, q - 1)).b
    R = V[0] + E * BOX
    return R, R / D.logq, E


# ENGINE TWO: THE FUNCTIONAL EQUATION WITH DIRECT SUMS


def direct_enclosure(D, k, Lp, M):
    s = D.s(k)
    q, N = D.q, D.N
    nums = [n for j in range(Lp) for n in D.level(j)]
    z = [iv.exp(-s * iv.log(n)) if n > 1 else iv.mpc(1) for n in nums]
    inv = [rat(1, n) for n in nums]
    pw = [iv.mpf(1) for n in nums]
    a = abs(s).b
    R = D.lev0(s)
    P = iv.mpc(1)
    for m in range(1, M + 1):
        P = P * (s + m - 1) / m
        Km = iv.mpc(0)
        for t in range(len(nums)):
            pw[t] = pw[t] * inv[t]
            Km = Km + z[t] * pw[t]
        tailK = (D.N1 * D.amin_pow(m) * rat(1, q ** (Lp * m)) / (1 - rat(1, q**m))).b
        Km = Km + tailK * BOX
        R = R + P * rat(((-1) ** m) * D.gamma(m), N * q**m) * Km
    rho0 = rat(D.amax, D.amin * q)
    tailM = (D.N1 * D.amin_pow(0) / (1 - rat(1, q ** (M + 1))) * hyp_tail(a, rho0, M)).b
    R = R + tailM * BOX
    return R, R / D.logq


# THE COLUMN


def column(D, k, lam0, mmax):
    s = D.s(k)
    q, N = D.q, D.N
    lam = [lam0]
    for m in range(1, mmax + 1):
        acc = iv.mpc(0)
        fall = iv.mpc(1)
        for j in range(1, m + 1):
            fall = fall * (m - j + 1 - s)
            acc = acc + fall / math.factorial(j) * rat(D.gamma(j), q**j) * lam[m - j]
        lam.append(-acc / N / (1 - rat(1, q**m)))
    return lam


def column_closed(D, k, lam0):
    s = D.s(k)
    mean = Fraction(sum(D.F), D.N)
    var = Fraction(sum(a * a for a in D.F), D.N) - mean * mean
    m1 = mean / (D.q - 1)
    m2 = var / (D.q * D.q - 1) + m1 * m1
    b1 = -m1
    b2 = m1 * m1 - m2 / 2
    c1 = -(s - 1) * lam0 * rat(b1.numerator, b1.denominator)
    c2 = (s - 1) * (s - 2) * lam0 * rat(b2.numerator, b2.denominator)
    return [c1, c2], [b1, b2]


# CONTROLS IN FLOATS


def dist_box(z, cre, cim):
    dx = max(Fraction(0), lo(z.real) - cre, cre - hi(z.real))
    dy = max(Fraction(0), lo(z.imag) - cim, cim - hi(z.imag))
    return math.hypot(float(dx), float(dy))


def level_bound(D, k, J):
    s = math.log(D.N) / math.log(D.q) + 2j * math.pi * k / math.log(D.q)
    a = abs(s)
    rho = D.amax / (D.amin * D.q ** (J + 1))
    return (
        D.N1
        * D.amin ** (-math.log(D.N) / math.log(D.q))
        * a
        * (1 - rho) ** (-a - 1)
        * rho
        * D.q
        / (D.q - 1)
        / math.log(D.q)
    )


def control_levels(D, k, J):
    s = math.log(D.N) / math.log(D.q) + 2j * math.pi * k / math.log(D.q)
    arr = np.array(D.F1, dtype=np.int64)
    out = []
    J = min(J, int(20 * math.log(2) / math.log(D.N)))
    for j in range(J + 1):
        if j > 0:
            arr = np.concatenate([D.q * arr + a for a in D.F])
        out.append(complex(np.sum(np.exp(-s * np.log(arr.astype(np.float64))))))
    return out


def count_le(xs, D):
    q, F, N = D.q, D.F, D.N
    x = xs.copy()
    nd = 1 + int(math.log(float(xs.max())) / math.log(q)) + 1
    digits = np.zeros((len(xs), nd), dtype=np.int64)
    for p in range(nd):
        digits[:, p] = x % q
        x = x // q
    cnt = np.zeros(len(xs), dtype=np.int64)
    tight = np.ones(len(xs), dtype=bool)
    for p in reversed(range(nd)):
        dp = digits[:, p]
        less = sum((dp > a).astype(np.int64) for a in F)
        cnt += np.where(tight, less * (N**p), 0)
        tight &= np.isin(dp, F)
    cnt += tight
    return cnt - 1


def control_fourier(D, k, J, M):
    u = (np.arange(M) + 0.5) / M
    x = np.floor(float(D.q) ** (u + J)).astype(np.int64)
    A = count_le(x, D)
    Phi = A / (float(D.N) ** (u + J))
    ck = complex(np.mean(Phi * np.exp(-2j * np.pi * k * u)))
    s0 = math.log(D.N) / math.log(D.q)
    s = s0 + 2j * math.pi * k / math.log(D.q)
    top = u >= 1 - s0
    dev = (
        float(np.abs(Phi[top] - D.N ** (1 - u[top])).max())
        if D.F == [0, 1] and D.q == 3
        else float("nan")
    )
    return ck, s * ck * math.log(D.q), float(Phi.min()), float(Phi.max()), dev


# REPORT


def report(D, k, L, second=False, controls=False, col=0, d=15):
    s = D.s(k)
    a = abs(s).b
    j0 = D.j_start()
    I = choose_I(a, rat(D.amax, D.amin * D.q ** (j0 + 1)), Fraction(1, 10**24))
    t0 = time.time()
    R, lam, E = level_limit(D, k, L, I)
    t1 = time.time()
    print(
        f"k = {k}: s = {fmt_c(s, 12)}, recursion from level {j0} ({len(D.level(j0))} terms summed directly) to L = {L}, I = {I}, {t1 - t0:.1f} s"
    )
    print(f"  R_k = lim_j Lev_j(s) = {fmt_c(R, d)}")
    print(f"  lambda_(0,k) = R_k/log q = {fmt_c(lam, d)}")
    dz, excl = dist_zero(lam, 9)
    print(
        f"  |lambda_(0,k)| >= {dz}, excludes zero: {excl}, width {width(lam):.1e}, level tail E_L = {fmt_up(E)}"
    )
    out = {"R": R, "lam": lam, "I": I, "excl": excl, "dz": dz}
    if second:
        t0 = time.time()
        _, lam2 = direct_enclosure(D, k, 13, 36)
        print(
            f"  engine two (direct sums to level 13, M = 36, {time.time() - t0:.1f} s): lambda = {fmt_c(lam2, 8)}, width {width(lam2):.1e}, meets engine one: {meets(lam, lam2)}"
        )
    if controls:
        lv = control_levels(D, k, 20)
        J = len(lv) - 1
        lamf = [v / math.log(D.q) for v in lv]
        dist = dist_box(lam, Fraction(lamf[J].real), Fraction(lamf[J].imag))
        print(
            f"  control, level sums by enumeration (Burnol 5.1), lambda from Lev_{J - 2}, Lev_{J - 1}, Lev_{J}: "
            + ", ".join(f"{v.real:.9f}{v.imag:+.9f}i" for v in lamf[J - 2 : J + 1])
        )
        print(
            f"  successive differences |Lev_(j+1) - Lev_j|, j = {J - 4}..{J - 1}: "
            + ", ".join(f"{abs(lv[j + 1] - lv[j]):.2e}" for j in range(J - 4, J))
            + f"; distance of Lev_{J}/log q from engine one {dist:.1e}, its own bound E_{J}/log q = {level_bound(D, k, J):.1e}, within: {dist <= level_bound(D, k, J)}"
        )
        for M in (2**15, 2**17):
            ck, Rf, pmin, pmax, dev = control_fourier(D, k, 30, M)
            lf = Rf / math.log(D.q)
            print(
                f"  control, Fourier coefficient of Phi on {M} midpoints, J = 30: c_k = {ck.real:.9f}{ck.imag:+.9f}i, "
                f"lambda = s c_k = {lf.real:.9f}{lf.imag:+.9f}i, distance from engine one {dist_box(lam, Fraction(lf.real), Fraction(lf.imag)):.1e}; "
                f"Phi in [{pmin:.6f}, {pmax:.6f}], max deviation of Phi from N^(1-u) on [1 - s_0, 1]: {dev:.1e}"
            )
    if col:
        lams = column(D, k, lam, col)
        closed, bs = column_closed(D, k, lam)
        for m in range(1, col + 1):
            line = f"  lambda_({m},k) = {fmt_c(lams[m], 10)}"
            if m <= 2:
                line += f", closed form via Theorem 7.4 with [t^{m}](1/E) = {bs[m - 1]}: {fmt_c(closed[m - 1], 10)}, meets: {meets(lams[m], closed[m - 1])}"
            print(line)
    return out


def main():
    band = "--band" in sys.argv
    t_all = time.time()
    print(
        "BURNOL RESIDUE: certified enclosures of Res K at s_(0,k) = log_q N + 2 pi i k/log q, mpmath.iv at 128 bits"
    )
    print(
        "rests on: Burnol 2026 Proposition 5.1 (lambda_(0,k) log q = lim_j Lev_j(s_(0,k))), the level recursion, the two tail lemmas, outward-rounded interval arithmetic"
    )
    D = Design(3, [0, 1])
    print(
        f"\n{D.name}: N = {D.N}, s_0 = {fmt_r(D.s0, 15)}, log 3 = {fmt_r(D.logq, 15)}, 2^(s_0) = {fmt_r(iv.mpf(2) ** D.s0, 9)}, 2^(1 - s_0) = {fmt_r(iv.mpf(2) ** (1 - D.s0), 9)}"
    )
    L = 40
    res = {}
    res[0] = report(D, 0, L)
    res[1] = report(D, 1, L, second=True, controls=True, col=3)
    ks = list(range(2, 11)) if band else []
    for k in ks:
        res[k] = report(D, k, L, controls=True)
    print(
        "\nBAND (printed by the generator; endpoints floored and ceiled at 12 decimals; |lambda| truncated at 9)"
    )
    print(
        "| `k` | `Re lambda_(0,k)` | `Im lambda_(0,k)` | `abs(lambda_(0,k)) >=` | zero excluded |"
    )
    print("|---|---|---|---|---|")
    for k in sorted(res):
        z = res[k]["lam"]
        print(
            f"| {k} | `{fmt_r(z.real, 12)}` | `{fmt_r(z.imag, 12)}` | `{res[k]['dz']}` | {'yes' if res[k]['excl'] else 'NO'} |"
        )

    D2 = Design(3, [0, 2])
    print(
        f"\n{D2.name}: N = {D2.N}; every element is twice one of {{0,1}}, so lambda scales by 2^(-s_(0,k))"
    )
    r2 = report(D2, 1, L, controls=True)
    scaled = res[1]["lam"] * iv.exp(-D.s(1) * iv.log(2))
    print(
        f"  2^(-s) times the {{0,1}} enclosure = {fmt_c(scaled, 15)}, meets the {{0,2}} enclosure: {meets(scaled, r2['lam'])}"
    )

    D3 = Design(3, [0, 1, 2])
    print(
        f"\n{D3.name}: K is the Riemann zeta function, s_0 = 1; the residue at 1 is 1 and every off-real lattice point is regular"
    )
    r30 = report(D3, 0, L)
    print(f"  contains 1: {contains(r30['lam'], Fraction(1), Fraction(0))}")
    r31 = report(D3, 1, L, controls=True)
    print(f"  contains 0: {contains(r31['lam'], Fraction(0), Fraction(0))}")
    print(f"\ntotal {time.time() - t_all:.1f} s")


if __name__ == "__main__":
    main()
