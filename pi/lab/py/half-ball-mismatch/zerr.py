import mpmath as mp
import sympy as sp

DPS = 50
SAMPLES = 50


def symbolic_residuals():
    a, u = sp.symbols("a u", positive=True)
    r = 1 - a**2 + a**2 * u**2
    cube = (-a * u + sp.sqrt(r)) ** 3
    t1 = -(a**3) * u**3
    t2 = 3 * a**2 * u**2 * sp.sqrt(r)
    t3 = -3 * a * u * r
    t4 = r ** sp.Rational(3, 2)
    even = sp.sqrt(r) * (1 - a**2 + 4 * a**2 * u**2)
    prim = u * r ** sp.Rational(3, 2)
    res = {
        "cube minus its four terms": sp.simplify(sp.expand(cube - (t1 + t2 + t3 + t4))),
        "even part minus sqrt(R)(1-a^2+4a^2u^2)": sp.simplify(even - (t2 + t4)),
        "d/du[u R^(3/2)] minus even part": sp.simplify(sp.diff(prim, u) - even),
    }
    boundary = sp.simplify(prim.subs(u, 1) - prim.subs(u, -1))
    closed = sp.simplify(sp.integrate(cube, (u, -1, 1)))
    return res, boundary, closed


def inner_integral(a):
    return mp.quad(
        lambda u: (-a * u + mp.sqrt(1 - a**2 + a**2 * u**2)) ** 3, [-1, 0, 1]
    )


def chord(a, phi):
    return -a * mp.cos(phi) + mp.sqrt(1 - a**2 * mp.sin(phi) ** 2)


def i_diam_direct():
    half = mp.pi / 2
    return mp.quad(
        lambda a, phi: chord(a, phi) ** 3 * mp.sin(phi),
        [-1, 0, 1],
        [0, half, mp.pi],
    )


def main():
    mp.mp.dps = DPS
    print("ZERR DECOMPOSITION")
    res, boundary, closed = symbolic_residuals()
    for name, value in res.items():
        print(f"  residual {name} = {value}")
    print(f"  boundary [u R^(3/2)] from -1 to 1 = {boundary}")
    print(f"  closed form of the inner integral = {closed}")

    worst = mp.mpf(0)
    for i in range(SAMPLES):
        a = mp.mpf(i + 1) / (SAMPLES + 1)
        worst = max(worst, abs(inner_integral(a) - 2))
    print(f"  inner integral at {SAMPLES} values of a, working digits {DPS}")
    print(f"  max deviation from 2 = {mp.nstr(worst, 5)}")

    i_diam = mp.quad(inner_integral, [-1, 0, 1])
    print(f"  I_diam by the constant inner integral = {mp.nstr(i_diam, 30)}")

    mp.mp.dps = 25
    direct = i_diam_direct()
    print(f"  I_diam by direct chord-cube quadrature = {mp.nstr(direct, 20)}")

    mp.mp.dps = DPS
    area = mp.pi / 2
    prob = mp.mpf(4) / (3 * area**2)
    exact = 16 / (3 * mp.pi**2)
    print(f"  Area(H) = pi/2, Area(H)^2 = pi^2/4 = {mp.nstr(area**2, 20)}")
    print(
        f"  P(chord crosses the diameter) = I_diam/(3 Area(H)^2) = {mp.nstr(prob, 20)}"
    )
    print(f"  16/(3 Pi^2) = {mp.nstr(exact, 20)}")
    print(f"  16/(3 Pi^2) to 16 places = {mp.nstr(exact, 16)}")
    print(f"  1 - 16/(3 Pi^2) to 16 places = {mp.nstr(1 - exact, 16)}")
    print(f"  |P - 16/(3 Pi^2)| = {mp.nstr(abs(prob - exact), 5)}")


main()
