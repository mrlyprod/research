import mpmath as mp
import numpy as np
import sympy as sp

TOP = 7
MC_DIMS = (3, 4, 5)
MC_SAMPLES = 10**8
BATCH = 10**6


def beta_sym(x, y):
    return sp.gamma(x) * sp.gamma(y) / sp.gamma(x + y)


def exact_miss(d):
    th, ps = sp.symbols("th ps", positive=True)
    alpha = sp.Rational(d * d - 1, 2)
    inner = sp.integrate(sp.sin(th) ** (d * d), (th, 0, sp.pi / 2 - ps))
    body = sp.expand_trig(sp.expand(inner * sp.sin(ps) ** (d - 2)))
    num = 2 * sp.integrate(body, (ps, 0, sp.pi / 2))
    denom = beta_sym(sp.Rational(1, 2), sp.Rational(d - 1, 2)) * beta_sym(
        sp.Rational(1, 2), alpha + 1
    )
    return sp.radsimp(sp.expand(sp.nsimplify(sp.simplify(2**d * num / denom))))


def quad_miss(d, dps):
    saved = mp.mp.dps
    mp.mp.dps = dps
    alpha = mp.mpf(d * d - 1) / 2
    half = mp.mpf(1) / 2

    def shell(h):
        return mp.betainc(mp.mpf(d - 1) / 2, half, 0, h * h, regularized=True)

    num = mp.quad(lambda h: (1 - h * h) ** alpha * shell(h), [0, half, 1])
    out = mp.mpf(2) ** d * num / mp.beta(half, alpha + 1)
    mp.mp.dps = saved
    return out


def half_ball(rng, count, dim):
    pts = rng.standard_normal((count, dim, dim))
    pts /= np.linalg.norm(pts, axis=2, keepdims=True)
    pts *= rng.random((count, dim, 1)) ** (1.0 / dim)
    pts[:, :, dim - 1] = np.abs(pts[:, :, dim - 1])
    return pts


def flat_face_rate(rng, dim, samples):
    hits = 0
    done = 0
    while done < samples:
        take = min(BATCH, samples - done)
        pts = half_ball(rng, take, dim)
        mat = np.concatenate([pts[:, :, : dim - 1], np.ones((take, dim, 1))], axis=2)
        rhs = pts[:, :, dim - 1][:, :, None]
        sol = np.linalg.solve(mat, rhs)[:, :, 0]
        slope = np.linalg.norm(sol[:, : dim - 1], axis=1)
        hits += int(np.count_nonzero(np.abs(sol[:, dim - 1]) <= slope))
        done += take
    return hits / samples


def main():
    mp.mp.dps = 30
    print(
        "VERSION H: d uniform points, the hyperplane through them misses the flat face"
    )
    print("  reduction: unoriented normal n with n_d = t >= 0 and signed offset h")
    print("  the section is a (d-1)-ball of radius sqrt(1-h^2) and stays above the")
    print("  flat face exactly when h >= sqrt(1-t^2)")
    print("  the within-hyperplane simplex integral scales as radius^(d^2-1)")
    print("  the half-ball Cartesian d-tuple measure is 2^-d of the full-ball one")

    forms = {}
    print("  exact values by symbolic reduction")
    for d in range(2, TOP + 1):
        p = exact_miss(d)
        forms[d] = p
        print(f"  d = {d}  {p}")

    print("  decimals to ten places")
    for d in range(2, TOP + 1):
        print(f"  d = {d}  {float(sp.N(forms[d], 30)):.10f}")

    print("  quadrature against the exact values, 60 and 80 working digits")
    for d in range(2, TOP + 1):
        lo = quad_miss(d, 60)
        hi = quad_miss(d, 80)
        mp.mp.dps = 80
        target = mp.mpf(str(sp.N(forms[d], 90)))
        print(
            f"  d = {d}  60 against 80 digits {mp.nstr(abs(lo - hi), 3)}"
            f"  quadrature against exact {mp.nstr(abs(hi - target), 3)}"
        )
        mp.mp.dps = 30

    print("  parity read off the exact values")
    for d in range(2, TOP + 1):
        powers = sorted(
            int(term.as_coeff_exponent(sp.pi)[1]) for term in sp.Add.make_args(forms[d])
        )
        shape = " + ".join("Q" if e == 0 else f"Q Pi^{e}" for e in powers)
        print(f"  d = {d}  pi powers {powers}  {shape}")

    rng = np.random.default_rng(20250828)
    print(f"  random-point check, {MC_SAMPLES} samples per dimension")
    for d in MC_DIMS:
        flat_exact = 1 - mp.mpf(str(sp.N(forms[d], 30)))
        est = mp.mpf(flat_face_rate(rng, d, MC_SAMPLES))
        se = mp.sqrt(est * (1 - est) / MC_SAMPLES)
        print(
            f"  d = {d}  flat-face estimate {mp.nstr(est, 8)}"
            f"  exact {mp.nstr(flat_exact, 8)}"
            f"  deviation {mp.nstr(est - flat_exact, 3)}  one sigma {mp.nstr(se, 3)}"
        )


main()
