import math
import subprocess
import time

import numpy as np

N = 30_000_000
JMIN = 8.0
JMAX = 23.5
PER_OCTAVE = 16
ZERO_HEIGHT = 300
PI2 = math.pi**2
T0 = time.time()


def clock(label):
    print(f"[{time.time() - T0:6.1f}s] {label}")


def totients(n):
    phi = np.arange(n + 1, dtype=np.int32)
    rem = np.arange(n + 1, dtype=np.int32)
    r = math.isqrt(n)
    small = np.ones(r + 1, dtype=bool)
    small[:2] = False
    for p in range(2, math.isqrt(r) + 1):
        if small[p]:
            small[p * p :: p] = False
    for p in np.flatnonzero(small):
        p = int(p)
        phi[p::p] -= phi[p::p] // p
        pk = p
        while pk <= n:
            rem[pk::pk] //= p
            pk *= p
    big = rem > 1
    phi[big] -= phi[big] // rem[big]
    return phi


def totients_brute(n):
    return np.array(
        [sum(1 for a in range(1, k + 1) if math.gcd(a, k) == 1) for k in range(n + 1)]
    )


def totient_sum_mobius(x):
    mu = np.ones(x + 1, dtype=np.int64)
    mu[0] = 0
    prime = np.ones(x + 1, dtype=bool)
    prime[:2] = False
    for p in range(2, x + 1):
        if prime[p]:
            prime[p * p :: p] = False
            mu[p::p] *= -1
            mu[p * p :: p * p] = 0
    total = 0
    for d in range(1, x + 1):
        if mu[d]:
            m = x // d
            total += int(mu[d]) * m * (m + 1) // 2
    return total


def bump_c2(u):
    v = (u - 1.0) * (2.0 - u)
    return np.where((u > 1.0) & (u < 2.0), 64.0 * v**3, 0.0)


def bump_cinf(u):
    v = (u - 1.0) * (2.0 - u)
    inside = (u > 1.0) & (u < 2.0)
    safe = np.where(inside, v, 1.0)
    return np.where(inside, np.exp(4.0 - 1.0 / safe), 0.0)


C2_POLY = 64.0 * np.polynomial.polynomial.polypow([-2.0, 3.0, -1.0], 3)


def mellin_c2(s):
    s = np.asarray(s, dtype=complex)
    total = np.zeros_like(s)
    for k, a in enumerate(C2_POLY):
        total += a * (2.0 ** (s + k) - 1.0) / (s + k)
    return total


GL_X, GL_W = np.polynomial.legendre.leggauss(2000)
GL_U = 1.5 + 0.5 * GL_X
GL_F = bump_cinf(GL_U) * 0.5 * GL_W


def mellin_cinf(s):
    s = np.asarray(s, dtype=complex)
    return np.array(
        [np.sum(GL_F * np.exp((z - 1.0) * np.log(GL_U))) for z in s.ravel()]
    ).reshape(s.shape)


def zeros_from_pari(height):
    script = f"""default(realprecision, 30);
z = lfunzeros(1, {height});
for(k = 1, #z, r = 1/2 + I*z[k]; a = zeta(r - 1); b = lfun(1, r, 1); print(z[k], " ", real(a), " ", imag(a), " ", real(b), " ", imag(b)))
"""
    out = subprocess.run(
        ["gp", "-q"], input=script, capture_output=True, text=True, check=True
    ).stdout
    rows = np.array(
        [[float(t) for t in line.split()] for line in out.strip().splitlines()]
    )
    gamma = rows[:, 0]
    zeta_left = rows[:, 1] + 1j * rows[:, 2]
    zeta_prime = rows[:, 3] + 1j * rows[:, 4]
    return gamma, zeta_left, zeta_prime


def smoothed_sum(phi, y, bump):
    lo = math.ceil(1.0 / y)
    hi = math.floor(2.0 / y)
    n = np.arange(lo, hi + 1, dtype=np.float64)
    return float(np.sum(phi[lo : hi + 1].astype(np.float64) * bump(n * y)))


def indicator_sum(prefix, y):
    lo = math.ceil(1.0 / y)
    hi = math.floor(2.0 / y)
    return int(prefix[hi] - prefix[lo - 1])


def fit(logq, logerr):
    a, b = np.polyfit(logq, logerr, 1)
    resid = logerr - (a * logq + b)
    return a, float(np.sqrt(np.mean(resid**2)))


def octave_rms(j, err):
    octaves = np.floor(j).astype(int)
    rows = []
    for k in sorted(set(octaves)):
        sel = octaves == k
        if sel.sum() >= PER_OCTAVE // 2:
            rows.append((k, float(np.sqrt(np.mean(err[sel] ** 2)))))
    return np.array(rows)


def slopes(name, j, err):
    logq = -2.0 * j * math.log(2.0)
    rows = octave_rms(j, err)
    ok = np.abs(err) > 0
    a_all, r_all = fit(logq[ok], np.log(np.abs(err[ok])))
    lq = -2.0 * (rows[:, 0] + 0.5) * math.log(2.0)
    le = np.log(rows[:, 1])
    half = len(rows) // 2
    a_lo, r_lo = fit(lq[:half], le[:half])
    a_hi, r_hi = fit(lq[half:], le[half:])
    a_oct, r_oct = fit(lq, le)
    print(
        f"{name}: octave rms slope in q {a_oct:.4f} (residual {r_oct:.3f}) over j in [{rows[0, 0]:.0f}, {rows[-1, 0] + 1:.0f}]"
    )
    print(
        f"{name}: lower window j in [{rows[0, 0]:.0f}, {rows[half - 1, 0] + 1:.0f}] slope {a_lo:.4f} (residual {r_lo:.3f}); upper window j in [{rows[half, 0]:.0f}, {rows[-1, 0] + 1:.0f}] slope {a_hi:.4f} (residual {r_hi:.3f})"
    )
    print(
        f"{name}: all-sample slope of log|E| in q {a_all:.4f} (residual {r_all:.3f}) on {int(ok.sum())} samples"
    )
    return a_oct


def main():
    clock(f"sieve to {N}")
    phi = totients(N)
    prefix = np.cumsum(phi, dtype=np.int64)
    clock("sieve done")
    brute = totients_brute(2000)
    print(
        f"brute-force totients agree to 2000: {bool(np.array_equal(brute, phi[:2001]))}"
    )
    x = 10**6
    print(
        f"totient sum at {x}: sieve {int(prefix[x])}, mobius route {totient_sum_mobius(x)}, equal {int(prefix[x]) == totient_sum_mobius(x)}"
    )
    j = np.arange(JMIN, JMAX + 1e-9, 1.0 / PER_OCTAVE)
    y = 2.0 ** (-j)
    c_ind = 6.0 / PI2 * 1.5
    c_c2 = 6.0 / PI2 * (24.0 / 35.0)
    f2_cinf = float(np.sum(GL_F * GL_U))
    c_cinf = 6.0 / PI2 * f2_cinf
    print(
        f"mellin transforms at 2: indicator 3/2, c2 bump 24/35 = {24 / 35:.12f} (closed form {mellin_c2(2.0).real:.12f}), cinf bump {f2_cinf:.12f}"
    )
    e_ind = np.array([yy * yy * indicator_sum(prefix, yy) - c_ind for yy in y])
    clock("indicator sums done")
    e_c2 = np.array([yy * yy * smoothed_sum(phi, yy, bump_c2) - c_c2 for yy in y])
    clock("c2 bump sums done")
    e_cinf = np.array([yy * yy * smoothed_sum(phi, yy, bump_cinf) - c_cinf for yy in y])
    clock("cinf bump sums done")
    yy = y[-1]
    lo = math.ceil(1.0 / yy)
    hi = math.floor(2.0 / yy)
    n = np.arange(lo, hi + 1, dtype=np.float64)
    prod = phi[lo : hi + 1].astype(np.float64) * bump_cinf(n * yy)
    noise = abs(float(np.sum(prod)) - math.fsum(prod.tolist())) * yy * yy
    print(
        f"summation noise at j = {JMAX}: |pairwise - fsum| scaled by y^2 is {noise:.3e} against |E_cinf| = {abs(e_cinf[-1]):.3e}"
    )
    print()
    print("SLOPES")
    a_ind = slopes("indicator", j, e_ind)
    a_c2 = slopes("c2 bump", j, e_c2)
    a_cinf = slopes("cinf bump", j, e_cinf)
    print(
        f"verdict: indicator {a_ind:.3f} against 1/2, c2 bump {a_c2:.3f} against 3/4, cinf bump {a_cinf:.3f} against 3/4"
    )
    print()
    print("EXPLICIT FORMULA")
    gamma, zeta_left, zeta_prime = zeros_from_pari(ZERO_HEIGHT)
    rho = 0.5 + 1j * gamma
    print(
        f"{len(gamma)} zeros to height {ZERO_HEIGHT} from PARI, first {gamma[0]:.6f}, last {gamma[-1]:.6f}"
    )
    for name, mellin, err in (
        ("c2 bump", mellin_c2, e_c2),
        ("cinf bump", mellin_cinf, e_cinf),
    ):
        coef = mellin(rho) * zeta_left / zeta_prime
        print(
            f"{name}: |c_rho| at zeros 1, 2, 10, 50, {len(gamma)}: "
            + ", ".join(f"{abs(coef[k]):.3e}" for k in (0, 1, 9, 49, len(gamma) - 1))
        )
        power = np.polyfit(np.log(gamma), np.log(np.abs(coef)), 1)[0]
        print(
            f"{name}: least-squares power of |c_rho| against gamma on {len(gamma)} zeros: {power:.2f}"
        )
        for K in (1, 10, 30, len(gamma)):
            model = np.array(
                [2.0 * np.sum((coef[:K] * yy ** (2.0 - rho[:K])).real) for yy in y]
            )
            resid = np.max(np.abs(err - model)) / np.max(np.abs(err))
            print(f"{name}: first {K:3d} zeros, max |E - sum| / max |E| = {resid:.3e}")
        amp = np.abs(err) / y**1.5
        print(
            f"{name}: |E|/y^(3/2) over the grid: min {amp.min():.4e}, max {amp.max():.4e}, bound 2 sum |c_rho| = {2 * np.sum(np.abs(coef)):.4e}"
        )
    print()
    print("SHARP CUTOFF")
    primes = [p for p in (10**6 + 3, 10**7 + 19) if phi[p] == p - 1]
    for p in primes:
        before = prefix[p - 1] - 3.0 / PI2 * (p - 1) ** 2
        after = prefix[p] - 3.0 / PI2 * p**2
        print(
            f"prime {p}: totient error jumps from {before:.1f} to {after:.1f}, jump {after - before:.1f} against phi(p) - 3(2p - 1)/pi^2 = {p - 1 - 3.0 / PI2 * (2 * p - 1):.1f}"
        )
    clock("done")


if __name__ == "__main__":
    main()
