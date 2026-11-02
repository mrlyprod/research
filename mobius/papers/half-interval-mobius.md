---
title: An Unconditional Mertens Bound on the Half Interval
lead: At every odd base `base >= 94939`, on the integers whose digits all lie in the lower half `{0, ..., (base-1)/2}`, the Mobius function cancels, `|M_F(x)| <= C A_F(x) exp(-c sqrt(log x))`, and the primes are counted, the von Mangoldt sum being `(base/(base+1)) A_F(x)` within the same error, both proved here with no hypothesis and computable constants; at a prime base `p` that set is `{k : p does not divide C(2k,k)}`, and a transfer operator whose cost per digit is `(2/pi) log base + O(1)` sets the base.
date: 2026-10-02
figure: paper-half-interval-mobius
---

Write the integers in an odd base and keep those whose every digit is at most `(base-1)/2`, the lower half of the digits kept as an interval. They form a set `S_F` with `A_F(x) = x^(alpha_base + o(1))` members up to `x`, where `fill = (base+1)/2` and `alpha_base = log fill/log base`: thinner than the integers by a power of `x`, with no multiplicative structure, and at a prime base `p` exactly the set of `k >= 1` for which `p` does not divide the central binomial coefficient `C(2k,k)`. This paper proves two things on it with no hypothesis, at every odd base `base >= 94939` and every `x >= 2`. The Mobius function cancels, `|M_F(x)| <= C A_F(x) exp(-c sqrt(log x))`, where `M_F(x)` sums `mu(n)` over the set up to `x`. And the primes are counted: the von Mangoldt sum over the set is `kappa_F A_F(x)` within the same error, `kappa_F = base/(base+1)`, which is `p/(p+1)` at a prime base `p`. `C` and `c` are computable from the base alone. The proof is the dissection of the earlier paper [An Unconditional Mertens Bound at Large Base](../../coprime/papers/unconditional-mertens-at-large-base.md), which asks the digit set for two consecutive digits and one number: a bound on the `l^1` mass of the digit transform over every shifted grid, growing like `fill^i base^(i alpha_1)` with `alpha_1 < 1/5`. The earlier paper runs it for `mu`, and Section 7 runs it for `Lambda`. On the lower half the transform is a Dirichlet kernel of length `fill`, and a transfer operator with the explicit weight `1 + |sin(pi t)|/2` proves that bound with `base^(alpha_1) = (2/pi) log base + O(1)`, the constant `2/pi` sharp from both sides, where the proved one-step bound of the earlier paper pays `sqrt(m)` per digit at `m` avoided digits, a power of the base at the `fill - 1` avoided here. The certificate clears `1/5` at every odd base from `94939`, by `1.3678 * 10^(-8)` at the wall, certified at `120` bits up to a monotone tail, and clears `1/4` from `3789`, where under the generalized Riemann hypothesis it gives a power saving. The prime asymptotic on such sets at sufficiently large base, with a log-power error, is a theorem of Maynard (2022), which remarks at one missing digit that its error could be made effective; what is new here is the written proof of that effective shape on the half interval, each step that differs from the Mobius case written out, the Mobius side, and the explicit base.

## Introduction

![The horizontal axis is `log10 base`, from `2` to `7`, a hairline at each decade; the vertical axis is the margin `bar - alpha_1` on a linear scale from `-0.181` to `0.102`, the zero rule at `0`. Two rising curves over the odd bases from 101 to 10^7 on that log axis, faint hairlines at the decades 10^3 to 10^6 and a dim zero rule across the middle: the blue curve is the margin 1/5 - alpha_1 of the certificate of Theorem 4.1, rising from -0.170 at base 101 to 0.042 at 10^7 and crossing the zero rule at 94939, just left of the hairline at 10^5; the orange curve above it is the margin 1/4 - alpha_1, the same curve lifted by 1/20, crossing at 3789; a disc marks each crossing.](paper-half-interval-mobius)

An unconditional Mertens bound on a digit set spends its base on one number. Expand the indicator of the strings of `n` digits over `F` in additive characters modulo `base^n`: every frequency `a/base^n` carries the Mobius exponential sum `sum mu(v) e(va/base^n)`, weighted by the digit transform `|hat F_n(a/base^n)|`, and far from every fraction of small denominator that sum is bounded by `x^(4/5 + eps)`, the minor-arc bound of Basak, Robles and Zaharescu (2023). Paid against the whole `l^1` mass of the transform, that bound beats the count of the set exactly when the mass grows like `fill^n base^(n alpha_1)` with `alpha_1 < 1/5`. The figure plots, for the lower half `F = {0, ..., (base-1)/2}`, the margin `1/5 - alpha_1` of the certificate this paper proves, `base^(alpha_1) = lambda/fill` with `lambda` the closed form of Theorem 4.1. At small bases the mass is too large; the margin climbs like `1/5 - log((2/pi) log base)/log base` and crosses zero at `94939` (Proposition 6.2). The orange curve is the same certificate against `1/4`, the bar the generalized Riemann hypothesis asks, and crosses at `3789`. Both curves are closed forms; the figure's census writer recomputes them in floating point and asserts both walls and both margins against the generator, `lab/py/interval-digits`.

**Theorem 1.1.** Let `base >= 94939` be odd, `fill = (base+1)/2` and `F = {0, 1, ..., fill-1}`. Let `S_F` be the set of positive integers whose every base-`base` digit lies in `F`, `A_F(x) = #{n in S_F : n <= x}`, `M_F(x) = sum_(n in S_F, n <= x) mu(n)` and `kappa_F = (base/phi(base)) #{f in F : gcd(f, base) = 1}/fill`. Then there are `C > 0` and `c > 0`, depending on `base` alone and effectively computable, such that for every `x >= 2`

```
|M_F(x)| <= C A_F(x) exp(-c sqrt(log x)) ,
|sum_(n in S_F, n <= x) Lambda(n) - kappa_F A_F(x)| <= C A_F(x) exp(-c sqrt(log x)) .
```

**Corollary 1.2 (the central binomial coefficients).** Let `p >= 94939` be prime and `A(x)` the number of integers `1 <= k <= x` such that `p` does not divide `C(2k,k)`. Then there are `C > 0` and `c > 0`, depending on `p` alone and effectively computable, such that for every `x >= 2`

```
|sum_(1 <= k <= x, p does not divide C(2k,k)) mu(k)| <= C A(x) exp(-c sqrt(log x)) ,
|sum_(1 <= k <= x, p does not divide C(2k,k)) Lambda(k) - (p/(p+1)) A(x)|
    <= C A(x) exp(-c sqrt(log x)) .
```

**Theorem 1.3 (under the generalized Riemann hypothesis).** Assume that `L(s, chi)` has no zero in `Re s > 1/2` for every Dirichlet character `chi`. Let `base >= 3789` be odd and `F`, `S_F`, `A_F`, `M_F` as in Theorem 1.1, put `alpha_base = log fill/log base`, and let `alpha_1 = log_base(lambda/fill)` with `lambda` the explicit constant of Theorem 4.1. Then `delta = (1/4 - alpha_1)/alpha_base > 0`, and for every `eps > 0` and every `x >= 2`, `|M_F(x)| <<_(base, eps) A_F(x)^(1 - delta + eps)`. As the base grows, `1 - delta` tends to `3/4`.

The saving in Theorem 1.1 is of zero-free-region shape, not a power, and `c` is tiny; the constants are computable from the base but not printed, because the implied constant of the minor-arc input and the zero-free constants of the character bounds are unstated at their sources (Section 7). Below `94939` the theorem says nothing. The two bounds share the base `94939`, the wall of the one certificate, and the dissection; they differ in the minor-arc input and the arithmetic input, so each has its own `C` and `c`, and Theorem 1.1 takes the larger `C` and the smaller `c`. The prime count adds the main-term constant `kappa_F = base/(base + 1)` (Lemma 7.3).

Why this set. The dissection of [An Unconditional Mertens Bound at Large Base](../../coprime/papers/unconditional-mertens-at-large-base.md), called the earlier paper below, its Theorem 1.1, holds at any number of avoided digits and reads the digit set only through a shifted-grid `l^1` certificate below `1/5` and two consecutive digits. The lower half keeps the pair `{0, 1}`, so everything rests on the certificate. There the proved one-step bound of the earlier paper, its Lemma 2.5, fails at every base: it is `sqrt(m) + Phi/base` per digit at `m` avoided digits, the `sqrt(m)` coming from Cauchy-Schwarz on the avoided digits, and at `m = fill - 1` that term alone puts the exponent above `1/2` at every odd base from `5`. The exact one-step constant `max G` is at most `(3/2) lambda` by Theorem 4.1, so it is logarithmic too, but a reading puts it a factor `sqrt 2` above `lambda` in the leading term (Conjecture 9.2), and that factor is what the weight below buys. The lower half is also the set where the cost should be smallest, because its transform is a single Dirichlet kernel, `|hat F(t)| = |sin(pi fill t)/sin(pi t)|`, whose `l^1` mass over one period is logarithmic in its length. The step from that heuristic to a certificate is a transfer operator `T`, Section 3, whose iterates are exactly the shifted grid sums, and a weight `phi(t) = 1 + |sin(pi t)|/2` with `T phi <= lambda phi` everywhere, Theorem 4.1. The weight is chosen so that its extra term is summed in closed form: `|hat F(u)| sin(pi u) = |sin(pi fill u)|`, and `fill` is the inverse of `2` modulo the base, so the weighted sum over the grid is one geometric series (Lemma 3.3). The growth constant obeys `lambda/fill <= (2/pi) log base + 2.60043004 + 3.4/base` from base `101`, and no shifted-grid certificate can do better than `(2/pi) log base - 1.31` (Section 5); the bar asks `lambda/fill < base^(1/5)`, and the closed form meets that power at `94939` (Section 6).

The prime asymptotic on the half interval at sufficiently large base, with a log-power error, is contained in Maynard (2022), Theorem 1.3, which takes consecutive excluded digits with `q - s >= q^(4/5 + eps)`, and that paper remarks after its Theorem 1.1, at one missing digit, that Siegel zeros play no role at its moduli, so its errors could be made effective of the shape proved here, a remark its Section 9 does not restate for Theorem 1.3; Maynard (2019) gives the order of magnitude of the primes on such sets. Neither carries a Mobius sum. What this paper adds is the written proof of the effective shape on the half interval (Section 7), the Mobius bound, and the explicit base `94939` for both sums (Section 6), placed against the literature in Section 8.

Section 2 fixes the set and proves the binomial reading. Section 3 sets up the transfer operator, Section 4 proves the certificate, Section 5 the constant `2/pi` from both sides, and Section 6 the walls with their margins. Section 7 restates what is reused from the earlier paper, proves the prime count by the same dissection with `Lambda` in place of `mu`, each changed step written out, and deduces Theorems 1.1 and 1.3 and Corollary 1.2. Section 8 places the result against the literature, Section 9 is the falsification, the limits and a conjecture on the true growth rate, and Section 10 says what is left open.

## The half interval

**Definition 2.1 (the half interval).** Fix an odd base `base >= 3` and put `fill = (base+1)/2`, so `base = 2 fill - 1`; the two words keep these meanings throughout, except in a statement made for any digit set, where `fill` is the number of its digits. The half interval is `F = {0, 1, ..., fill-1}`: `m = fill - 1` avoided digits, all consecutive, the upper half. `S_F`, `A_F(x)` and `M_F(x)` are as in Theorem 1.1, and `alpha_base = log fill/log base` is the dimension of `S_F`. For `n >= 0`, `D_n` is the set of the `fill^n` integers `0 <= u < base^n` whose `n` padded digits lie in `F`, so `S_F` below `base^n` is `D_n` less `{0}`. With `e(t) = exp(2 pi i t)`, the digit transform is `hat F(t) = sum_(a in F) e(at)` and its level form is `hat F_n(t) = prod_(i < n) hat F(base^i t) = sum_(u in D_n) e(ut)`, unnormalised as in the earlier paper, so `|hat F_n| <= fill^n` and

```
|hat F(t)| = |sin(pi fill t)/sin(pi t)| ,   read as fill at the integers .
```

`|hat F|` has period `1` and is even. Write `||z||` for the distance from `z` to the nearest integer. The set is the subject of the section The half interval of the note [The half interval and the even digits](../notes/half-interval-even-digits.md), which records every statement proved here.

**Lemma 2.2 (the central binomial coefficients).** At a prime base `p`, `v_p(C(2k,k))` is the number of carries in the base-`p` addition `k + k`, so `S_F = {k >= 1 : p does not divide C(2k,k)}`.

**Proof.** Legendre's formula `v_p(n!) = (n - s_p(n))/(p - 1)`, with `s_p` the base-`p` digit sum, gives `v_p(C(2k,k)) = (2 s_p(k) - s_p(2k))/(p - 1)`. In the addition `k + k` each carry removes `p` from one position and adds `1` at the next, so `s_p(2k) = 2 s_p(k) - (p - 1) c` with `c` the number of carries, and `v_p(C(2k,k)) = c`. No carry enters the lowest position, so a first carry occurs at the lowest position whose digit `d` has `2d >= p`; hence `c = 0` exactly when every digit of `k` is at most `(p-1)/2 = fill - 1`. □

This is the case `k + k` of the theorem of Kummer (1852) that `v_p(C(a + b, a))` counts the carries in `a + b`. So Corollary 1.2 is Theorem 1.1 read at a prime base, with `A(x) = A_F(x)`.

**Remark 2.3 (doubling).** Every digit `d` of a member of `S_F` has `2d <= base - 1`, so doubling carries nowhere and `2 S_F` is the set of positive integers whose digits lie in the even digits `{0, 2, ..., base-1}`: a digit set with no two consecutive digits, all of whose members are even. Theorem 1.1 does not reach it; Section 10 says what is missing.

## The transfer operator

**Definition 3.1 (the operator and the grid sums).** For any digit set `F` at the base `base` and any bounded `phi` of period `1`, put

```
(T phi)(t) = sum_(r < base) |hat F((t+r)/base)| phi((t+r)/base) ,   G = T 1 ,
Sigma_i(s) = sum_(a < base^i) |hat F_i(s + a/base^i)| .
```

`T` is positive and keeps period `1`: replacing `t` by `t + 1` sends the term `r` to the term `r + 1`, and the term `base - 1` to the term `0` moved by `1`.

**Lemma 3.2 (the transfer identity).** For every digit set, every `i >= 0`, every `phi` of period `1` and every real `t`,

```
(T^i phi)(t) = sum_(a < base^i) |hat F_i((t+a)/base^i)| phi((t+a)/base^i) ,
so   Sigma_i(s) = (T^i 1)(base^i s) .
```

Hence if `1 <= phi <= Phi` and `T phi <= lambda phi` everywhere, then `Sigma_i(s) <= Phi lambda^i` at every level `i` and every real shift `s`: a shifted-grid certificate in the sense of Theorem 1.1 of the earlier paper, with `C_F = Phi` and `base^(alpha_1) = lambda/fill`.

**Proof.** Induct on `i`, the case `i = 0` being trivial. By the definition of `T` and the case `i - 1` at the point `(t + r)/base`,

```
(T^i phi)(t) = sum_(r < base) |hat F((t+r)/base)| sum_(a' < base^(i-1))
                 |hat F_(i-1)((t + r + base a')/base^i)| phi((t + r + base a')/base^i) .
```

Put `a = r + base a'`, which runs once over `0 <= a < base^i`. Since `hat F` has period `1`, `hat F((t+r)/base) = hat F((t+a)/base) = hat F(base^(i-1) (t+a)/base^i)`, the factor at position `i - 1` of `hat F_i((t+a)/base^i)`, and `hat F_(i-1)((t+a)/base^i)` holds the positions below it; together they are `hat F_i((t+a)/base^i)`. At `phi = 1` and `t = base^i s` the right side is `Sigma_i(s)`. `T` is positive, so `T phi <= lambda phi` gives `T^i phi <= lambda^i phi` by induction, and `Sigma_i(s) = (T^i 1)(base^i s) <= (T^i phi)(base^i s) <= Phi lambda^i`. The certificate of the earlier paper reads `Sigma_i(s) <= C_F fill^i base^(i alpha_1)`, which is this with `C_F = Phi` and `fill base^(alpha_1) = lambda`. □

At `phi = 1` the lemma is the peel of Lemma 2.3(i) of the earlier paper, with `lambda = max G`, the one-step constant `B`. A weight that is not constant lets a point where `G` is large borrow from points where `phi` is small; the certificate below spends that freedom on one sine.

**Lemma 3.3 (the sine sum).** On the half interval, with `sigma` the fractional part of `fill t`,

```
S(t) = sum_(r < base) |hat F((t+r)/base)| |sin(pi (t+r)/base)|
     = cos(pi (sigma - 1/2)/base)/sin(pi/(2 base)) <= csc(pi/(2 base)) .
```

**Proof.** `|hat F(u)| |sin(pi u)| = |sin(pi fill u)|` at every `u`, the integers included, so `S(t) = sum_(r < base) |sin(pi (fill t + fill r)/base)|`. Since `2 fill = base + 1`, `fill` is the inverse of `2` modulo `base`, and `fill r` runs once over the residues modulo `base` as `r` does; `|sin(pi z/base)|` is unchanged when `z` moves by a multiple of `base`, and moving `fill t` to `sigma` by an integer permutes the residues again, so `S(t) = sum_(j < base) sin(pi (sigma + j)/base)`, every argument in `[0, pi)`. That is the imaginary part of `e(sigma/(2 base)) sum_(j < base) e(j/(2 base)) = e(sigma/(2 base)) 2/(1 - e(1/(2 base)))`, and `2/(1 - e(1/(2 base))) = i e(-1/(4 base))/sin(pi/(2 base))`. □

## The certificate

**Theorem 4.1 (the certificate).** At every odd base `base >= 3` put `c_P = sqrt 2 - 4/pi`, `gamma` Euler's constant,

```
X_0 = (base/pi)(log(base + 3) + gamma + log tan(3 pi/8 + pi/(4 base))) + c_P (base + 1)^2/(8 base) ,
lambda = fill + X_0 + csc(pi/(2 base))/2 ,                                               (4.1)
```

and `phi(t) = 1 + |sin(pi t)|/2`. Then `T phi <= lambda phi` everywhere on the half interval, so `Sigma_i(s) <= (3/2) lambda^i` at every level `i` and every real shift `s`: the half interval carries the shifted-grid certificate `C_F = 3/2`, `base^(alpha_1) = lambda/fill`.

**Proof.** `T phi` and `lambda phi` have period `1` and are even, because `|hat F|` and `phi` are and `r -> base - 1 - r` carries the points `(1 - t + r)/base` to `1 - (t + r)/base`; so take `t in [0, 1/2]`, and put `a = sin(pi t/2)` and `b = cos(pi t/2)`. By Lemma 3.3, `T phi = G + S/2 <= G + csc(pi/(2 base))/2`, so it suffices to prove

```
G(t) <= fill + X_0 + sin(pi t) (X_0/2 + K_2 base/(4 pi)) ,                               (4.2)
K_2 = 2 (4/3 + (36/35)(7 zeta(3)/8 - 1)) < 2.7733 ,
```

since then `T phi <= lambda + sin(pi t)(X_0/2 + K_2 base/(4 pi)) <= lambda + (lambda/2) sin(pi t) = lambda phi(t)`, the last step because `K_2 base/(2 pi) < 0.45 base < fill`.

The points. The `base` points `u_r = (t + r)/base` sit at distance `delta = ||u_r||` from the integers: the low points `r < fill` at `delta = (t + r)/base`, and the high points `r = base - 1 - r'`, `0 <= r' <= fill - 2`, at `delta = (1 - t + r')/base`, every `delta` at most `1/2`. There `|hat F(u_r)| = |sin(pi fill delta)|/sin(pi delta)`, and `fill delta = (t + r + delta)/2` at a low point and `(1 - t + r' + delta)/2` at a high one. Expanding the sine of `pi fill delta` in the angles `pi t/2` and `pi delta/2` and dividing by `sin(pi delta) = 2 sin(pi delta/2) cos(pi delta/2)`, every point weighs a combination of

```
P(delta) = 1/(2 sin(pi delta/2)) ,   Q(delta) = 1/(2 cos(pi delta/2)) ,
```

as in the table, with `J = floor((fill-2)/2)` and `K = floor((fill-1)/2)`.

| points | `delta` | range | weight |
| --- | --- | --- | --- |
| low, `r = 0` | `t/base` | one point | at most `fill` |
| low, `r = 2k` | `(2k + t)/base` | `1 <= k <= K` | `aP + bQ` |
| low, `r = 2j + 1` | `(2j + 1 + t)/base` | `0 <= j <= J` | `bP - aQ` |
| high, `r' = 2j` | `(2j + 1 - t)/base` | `0 <= j <= J` | `bP + aQ` |
| high, `r' = 2k - 1` | `(2k - t)/base` | `1 <= k <= K` | `abs(bQ - aP)` |

For instance at a low odd point `pi fill delta = pi j + pi/2 + pi t/2 + pi delta/2`, whose sine is `cos(pi t/2 + pi delta/2) = b cos(pi delta/2) - a sin(pi delta/2)` up to sign, and the angle `pi t/2 + pi delta/2` is at most `pi/2`, so the weight is `bP - aQ >= 0`. At a high odd point `pi fill delta = pi k + pi delta/2 - pi t/2`, so the weight is `abs(bQ - aP)`, which is `bQ - aP` where `delta > t` and at most `aP` where `delta <= t`. Every `|hat F|` is at most `fill`, which pays the point `r = 0`. So

```
G(t) <= fill + b Sigma_P + b Sigma_Q + a U + a Sigma_Q' ,
```

where `Sigma_P` sums `P` over the low odd and high even points, `Sigma_Q` sums `Q` over the low even and all the high odd points, `U = sum_(low even) P + sum_(high odd, delta <= t) P - sum_(high odd, delta > t) P`, and `Sigma_Q' = sum_(j <= J) (Q((2j + 1 - t)/base) - Q((2j + 1 + t)/base))` collects the `aQ` terms; the high odd points with `delta <= t` are bounded by `aP`, their `-bQ` dropped.

The `aQ` terms. `Q` increases on `[0, 1/2]`, and the high even point `(2j + 1 - t)/base` sits below the low odd point `(2j + 1 + t)/base`, so every term of `Sigma_Q'` is at most `0`.

Two inequalities. On `(0, 1/2]`, `1/(pi delta) <= P(delta) <= 1/(pi delta) + c_P delta`: with `z = pi delta/2` in `(0, pi/4]`, `2 P(delta) - 2/(pi delta) = csc z - 1/z`, which is positive, and convex because its Taylor series at `0` has positive coefficients, so it lies under its chord, `csc z - 1/z <= (sqrt 2 - 4/pi)(4/pi) z = 2 c_P delta`. And `Q` is convex on `[0, 1)`, since `sec` is convex on `[0, pi/2)`.

The `bP` terms. The low odd and high even points pair as `(m + t)/base` and `(m - t)/base` with `m = 2j + 1`, so

```
Sigma_P <= (base/pi) sum_(j <= J) (1/(m + t) + 1/(m - t)) + (c_P/base) sum_(j <= J) 2m ,
1/(m + t) + 1/(m - t) = 2/m + 2 t^2/(m (m^2 - t^2)) .
```

First, `sum_(j <= J) 2/(2j + 1) = psi(J + 3/2) + gamma + 2 log 2 <= log(4J + 6) + gamma <= log(base + 3) + gamma`, by `psi(x) <= log x - 1/(2x)` and `4J + 6 <= 2 fill + 2`. Second, with `t <= 1/2`, `1/(1 - t^2) <= 4/3` at `m = 1` and `m^2/(m^2 - t^2) <= 36/35` at `m >= 3`, so the `t^2` terms sum to at most `2 t^2 (4/3 + (36/35) sum_(m odd, m >= 3) m^(-3)) = K_2 t^2`, the odd cubes summing to `7 zeta(3)/8 - 1`. Third, `sum_(j <= J) 2m = 2 (J + 1)^2 <= (base + 1)^2/8`, since `J + 1 <= fill/2`. So `Sigma_P <= (base/pi)(log(base + 3) + gamma + K_2 t^2) + c_P (base + 1)^2/(8 base)`.

The `bQ` terms. The low even points `(2k + t)/base` and the high odd points `(2k - t)/base`, `1 <= k <= K`, are two progressions of spacing `2/base`. By the Hermite-Hadamard inequality on the convex `Q`, `Q(x) <= (base/2) int_(x - 1/base)^(x + 1/base) Q`, and the intervals of one progression are disjoint and lie in `[0, (2K + 1 + t)/base]`, inside `[0, 1/2 + 1/base]` because `2K + 1 <= fill`. With `int_0^y Q = (1/pi) log tan(pi/4 + pi y/4)`, `Sigma_Q <= 2 (base/2)(1/pi) log tan(3 pi/8 + pi/(4 base)) = (base/pi) log tan(3 pi/8 + pi/(4 base))`. So `b Sigma_P + b Sigma_Q <= b X_0 + (base/pi) K_2 t^2`.

The `aP` terms. This is the one step where the sign of the high odd weights is used. The low even points, `delta = (2k + t)/base` with `1 <= k <= K <= (base - 1)/4`, have `delta <= 1/2`, so by the chord, `1/(2k + t) <= 1/(2k)` and `H_K <= log K + 1`,

```
sum_(low even) P <= (base/pi) sum_(k <= K) 1/(2k) + c_P K/2
                 <= (base/(2 pi))(log((base - 1)/4) + 1) + c_P (base - 1)/8 ,
```

the right side positive also at `K = 0`, where the sum is empty. The high odd points, `delta = (2k - t)/base`, have `delta <= t` exactly when `k <= X = t (base + 1)/2`; put `k* = min(floor(X), K)`. If `k* = 0` every high odd point is subtracted and their signed sum is at most `0`. Otherwise the points `k <= k*` carry at most `(base/pi)(2/3 + (1/2) log((4k* - 1)/3)) + c_P X t`: the term `k = 1` has `1/(2 - t) <= 2/3`, every `k >= 2` has `1/(2k - t) <= 1/(2k - 1/2) <= (1/2) int_(k-1)^k dx/(x - 1/4)`, and the `k*` values of `delta` are each at most `t`, with `k* <= X`. The points `k > k*` all have `delta > t` and are subtracted, and they carry at least `(base/(2 pi)) log((K + 1)/(k* + 1))`, from `P(delta) >= 1/(pi delta)` and `1/(2k - t) >= 1/(2k) >= (1/2) int_k^(k+1) dx/x`. Now `(4k* - 1)(k* + 1) <= 4 X^2 + 3 X <= (base + 1)(base + 4)/4`, since `X <= (base + 1)/4`, and `K + 1 >= fill/2 = (base + 1)/4`, so with `X t <= (base + 1)/8` the signed high odd sum is at most

```
(base/pi)(2/3 + (1/2) log((base + 4)/3)) + c_P (base + 1)/8 .
```

Adding the low even bound,

```
U <= (base/pi)((1/2) log((base - 1)(base + 4)/12) + 7/6) + c_P base/4
  <= (base/pi)(log(base + 4) - 0.0757) + c_P base/4 ,
```

since `(base - 1)(base + 4) <= (base + 4)^2` and `7/6 - (1/2) log 12 < -0.0757`. Against `X_0`, with `log tan(3 pi/8 + pi/(4 base)) >= log tan(3 pi/8) = log(1 + sqrt 2)`, `gamma + log(1 + sqrt 2) > 1.4585`, `log((base + 3)/(base + 4)) >= -1/(base + 3)` and `(base + 1)^2/(8 base) >= 0`,

```
X_0 - U >= (base/pi)(1.4585 + 0.0757 - 1/(base + 3)) - c_P base/4 > 0 ,
```

because `(1.5342 - 1/6)/pi > 0.43 > c_P/4`. So `U < X_0`, and `a U <= a X_0`, also when `U < 0`.

Assembly. `G(t) <= fill + (a + b) X_0 + (base/pi) K_2 t^2`. Since `(1 - a)(1 - b) >= 0`, `a + b <= 1 + ab = 1 + sin(pi t)/2`; since `sin(pi t) >= 2t` on `[0, 1/2]`, `t^2 <= t/2 <= sin(pi t)/4`. That is (4.2). □

The proof uses the shape of the half interval at three points: `|hat F|` is one kernel `sin(pi fill u)/sin(pi u)`, which gives the table; `fill` is the inverse of `2`, which gives Lemma 3.3; and the weights at the high odd points change sign at `delta = t`, which the signed sum `U` follows. Fact 9.1 checks on grids of shifts the bounds on the harmonic parts of `Sigma_P` and `U` and on `Sigma_Q`, the closed form of Lemma 3.3, the bound (4.2) and `T phi <= lambda phi`, and on small grids the level sums against `(3/2) lambda^i`; the elementary inequalities between those steps are not checked separately. The bound (4.2) is tight to `1.5 * 10^(-4)` of its size at `t = 0` (`lab/py/interval-digits`, verb `check`).

## The constant `2/pi`

The growth constant of Theorem 4.1 is `fill` times `(2/pi) log base + O(1)`, and no shifted-grid certificate does better in the leading term.

**Proposition 5.1 (the upper constant).** At every odd base `base >= 101`,

```
lambda/fill <= (2/pi) log base + c_inf + 3.4/base ,
c_inf = 1 + 2/pi + (2/pi)(gamma + log(1 + sqrt 2)) + c_P/4 <= 2.60043004 .
```

**Proof.** `lambda/fill = 1 + 2 X_0/(base + 1) + csc(pi/(2 base))/(base + 1)`. In `2 X_0/(base + 1)` the factor `2 base/(base + 1) < 2` multiplies `(1/pi)(log(base + 3) + gamma + log tan(3 pi/8 + pi/(4 base)))`, which is positive; `log(base + 3) <= log base + 3/base`; `log tan z` has slope `2/sin(2z)`, below `2.9` on `[3 pi/8, 3 pi/8 + pi/(4 base)]` once `base >= 101`, so `log tan(3 pi/8 + pi/(4 base)) <= log(1 + sqrt 2) + 2.9 pi/(4 base)`; and the last term of `X_0` gives `c_P (base + 1)/(4 base) = c_P/4 + c_P/(4 base)`. For the cosecant, `csc z - 1/z` is convex on `(0, pi/2]` and equals `1 - 2/pi` at the end, so `csc z <= 1/z + (2/pi)(1 - 2/pi) z` and `csc(pi/(2 base))/(base + 1) < 2/pi + (1 - 2/pi)/(base (base + 1))`. Collecting, `lambda/fill <= (2/pi) log base + c_inf + E/base` with `E = (2/pi)(3 + 2.9 pi/4) + c_P/4 + (1 - 2/pi)/(base + 1) < 3.399` at `base >= 101`. The constant `c_inf` is bounded at `120` bits, and the inequality is also checked directly, the closed form (4.1) against the right side at `120` bits, at every odd base `101..3001` and at `94939`, `200001`, `10^6 + 1` and `10^8 + 1`, the smallest gap at least `1.3332 * 10^(-7)` (`lab/py/interval-digits`, verb `wall`). □

**Proposition 5.2 (the lower constant).** At every odd base `base >= 3` and every `t`, `G(t) >= (base/pi) log((base + 2)/5) - fill/4` and `G(t) >= base`. At every odd base `base >= 9`, `min G/fill >= (2/pi) log base - 1.31 > 0`, so every grid sum, shifted or not, has `Sigma_i(s) >= (((2/pi) log base - 1.31) fill)^i` at every level `i` and shift `s`; and `Sigma_i(s) >= base^i` at every odd base.

**Proof.** Keep `t in [0, 1/2]` and the table of Theorem 4.1. Every weight is nonnegative, so drop the point `r = 0`, the high odd points and the `bQ` of the low even points. The low odd and high even points carry `b (P((m + t)/base) + P((m - t)/base)) - a (Q((m + t)/base) - Q((m - t)/base))` at `m = 2j + 1`. With `P(delta) >= 1/(pi delta)` and `1/(m + t) + 1/(m - t) >= 2/m`, the first part is at least `b (base/pi) sum_(j <= J) 2/(2j + 1) >= b (base/pi)(log(base + 1) + gamma - 2/fill)`, by `psi(x) >= log x - 1/x`, `4J + 6 >= base + 1` and `J + 3/2 >= fill/2`. `Q` runs from `1/2` to `1/sqrt 2` on `[0, 1/2]`, so each difference of `Q` lies in `[0, 1/2]`, and there are `J + 1 <= fill/2` of them, so with `a <= 1` they cost at most `fill/4`. The low even points carry at least `a (base/pi) sum_(k <= K) 1/(2k + 1/2) >= a (base/(2 pi)) log((4K + 5)/5) >= a (base/(2 pi)) log((base + 2)/5)`, by `1/(k + 1/4) >= int_k^(k+1) dx/(x + 1/4)` and `4K + 5 >= base + 2`. Since `log(base + 1) + gamma - 2/fill >= log((base + 2)/5) >= 0` and `b + a/2 >= 1` on `[0, 1/2]`, the latter because `cos z + (1/2) sin z` is concave on `[0, pi/4]` and at least `1` at both ends, `G(t) >= (base/pi) log((base + 2)/5) - fill/4`. Parseval on `Z/base` gives `sum_r |hat F((t+r)/base)|^2 = base fill` at every `t`, and each term is at most `fill`, so `G >= base`. Dividing by `fill`, `min G/fill >= (2 base/(pi (base + 1))) log((base + 2)/5) - 1/4`. That is at least `(2/pi) log base - 1.31` at every odd base `9..9999`, certified at `120` bits, and at every base from `101` because the difference is at least `1.06 - (2/pi) log 5 - (2/pi) log((base + 2)/5)/(base + 1) >= 0.016518` there, the last term decreasing from `base 11`; the factor `(2/pi) log base - 1.31` is positive from `9` and negative at `7` (`lab/py/interval-digits`, verb `wall`). Finally `T^i 1 >= (min G)^i` pointwise, because `T` is positive and `T 1 = G`, and Lemma 3.2 turns both lower bounds on `G` into lower bounds on `Sigma_i(s)`. □

So the `l^1` cost of the half interval per digit is `(2/pi) log base + O(1)` from both sides, in the normalisation `base^(alpha_1)` of the certificate: any shifted-grid certificate `(C_F, alpha_1)` has `base^(alpha_1) >= (2/pi) log base - 1.31` at every odd base from `9`, since its `C_F` is fixed while `i` grows, and the one of Theorem 4.1 has `base^(alpha_1) <= (2/pi) log base + 2.64` from `101`, by Proposition 5.1, since `2.60043004 + 3.4/101 < 2.6342`. Measured against `base`, as the Parseval floor `Sigma_i >= base^i` measures it, the cost is `(1/pi) log base` per digit. The weight of Theorem 4.1 is what reaches the constant `2/pi`: the constant weight `phi = 1` pays the one-step constant `max G`, which reads `sqrt 2` times larger in the leading term (Conjecture 9.2).

## The walls

The bar `1/5` asks `lambda/fill < base^(1/5)`, a logarithm against a power, and the power wins from a computable base on.

**Lemma 6.1 (the tails climb).** For `e` in `{1/5, 1/4}` put `h_e(base) = base^e - (2/pi) log base - c_inf - 3.4/base`. Then `h_(1/5)` increases on `base >= 327` and `h_(1/4)` increases on `base >= 43`.

**Proof.** `h_e'(base) = e base^(e-1) - 2/(pi base) + 3.4/base^2 > (e base^e - 2/pi)/base`, which is positive once `base^e > 2/(pi e)`, that is once `base > (10/pi)^5 = 326.78...` at `e = 1/5` and `base > (8/pi)^4 = 42.05...` at `e = 1/4`. □

**Proposition 6.2 (the walls).** The certificate of Theorem 4.1 has `alpha_1 < 1/5` at every odd base `base >= 94939` and `alpha_1 > 1/5` at `94937`; it has `alpha_1 < 1/4` at every odd base `base >= 3789` and `alpha_1 > 1/4` at `3787`. At the walls and the bases below them,

| bar | wall | `bar - alpha_1` at the wall | `base^bar - lambda/fill` at the wall | the same at the base below |
| --- | ---: | ---: | ---: | ---: |
| `1/5` | `94939` | `>= 1.3678 * 10^(-8)` | `>= 1.5514 * 10^(-6)` | `<= -2.6733 * 10^(-5)` |
| `1/4` | `3789` | `>= 7.9625 * 10^(-6)` | `>= 5.1473 * 10^(-4)` | `<= -1.8427 * 10^(-4)` |

and farther out `alpha_1 <= 0.1993872` at `100003`, `0.1761232` at `1000003` and `0.1331636` at `10^9 + 7`.

**Proof.** `alpha_1 < e` is `lambda/fill < base^e`, and by Proposition 5.1 `base^e - lambda/fill >= h_e(base)` at every odd `base >= 101`. In interval arithmetic at `120` bits the least integer at which the tail `h_(1/5)` is certified positive is `94946`, and for `h_(1/4)` it is `3793`; both lie above the thresholds of Lemma 6.1, so `h_e > 0` at every real `base` beyond, and the bar holds at every odd base from `94947` and from `3793`. Below the tail, the closed form (4.1) is evaluated at `120` bits at every odd base `94939..94947` and `3789..3793`, every gap `base^e - lambda/fill` certified positive, and at `94937` and `3787` the gap is certified negative, with the values in the table. Each bar is read as the exact interval `1/5` or `1/4`, and `bar - alpha_1` is printed from the interval of `log(lambda/fill)/log base`, never by differencing two rounded numbers; lower ends are truncated down and upper ends rounded up (`lab/py/interval-digits`, verb `wall`, under a second). □

The margin `1/5 - alpha_1` at the wall is thin, `1.3678 * 10^(-8)`, and it enters the proof of Theorem 1.1 only in region A: for `mu` through `eps = (1/5 - alpha_1)/4`, which fixes the implied constant of the minor-arc input, and for `Lambda`, which has no `eps`, as the power saving `y^(alpha_1 - 1/5)` of Section 7. Either way it fixes how large `x` must be before the power saving of region A takes over: it moves `C`, never `c`. The wall is the certificate's: it says where this `phi` clears the bar, not where the route stops (Conjecture 9.2).

## The reduction to the dissection

The Mobius half of the dissection is reused unchanged. This section restates it, checks each place where the digit set enters, proves the constant `kappa_F`, and proves the prime count by the same dissection with `Lambda` in place of `mu`, each step that changes written out.

**Theorem 7.1 (the dissection, restated).** Let `base >= 3` and let `F` be a set of `fill` digits with `2 <= fill < base`, so that at least one digit is avoided, with `S_F`, `A_F`, `M_F`, `D_n` and `Sigma_i` as above. Suppose `F` contains two consecutive digits and carries a shifted-grid certificate below `1/5`: constants `C_F >= 1` and `alpha_1 < 1/5` with `Sigma_i(s) <= C_F fill^i base^(i alpha_1)` at every `i >= 0` and every real `s`. Then there are `C > 0` and `c > 0`, depending on `base` and `F` alone and effectively computable, such that for every `x >= 2`

```
(i)   |M_F(x)| <= C A_F(x) exp(-c sqrt(log x)) ,
(ii)  |sum_(n in S_F, n <= x) Lambda(n) - kappa_F A_F(x)| <= C A_F(x) exp(-c sqrt(log x)) ,
```

with `kappa_F = (base/phi(base)) #{f in F : gcd(f, base) = 1}/fill`. Part (i) is proved in the earlier paper. Part (ii) is proved below, as the subsection The dissection for the von Mangoldt function of the note [The coprimality spine](../../coprime/notes/coprime.md) gives it: Lemma 7.2 carries its main term, and Lemmas 7.4 and 7.5, Proposition 7.6 and the proof after it carry the rest.

Part (i) is Theorem 1.1 of [An Unconditional Mertens Bound at Large Base](../../coprime/papers/unconditional-mertens-at-large-base.md), stated there for any set of `m >= 1` avoided digits, with the remark after its proof that no step uses `m = 1`. The proof cuts `x` into blocks `Py + D_n` with `y = base^n` within `x^(o(1))` of `x`, expands each block in additive characters modulo `y`, and cuts the frequencies into four regions by Dirichlet approximation at `Q = y^(3/5)` and `Z = exp(C_0 sqrt(log x))`: the minor arcs A, the middle denominators B, the fractions C1 whose denominator has a prime outside the base, and the fractions C2 whose denominator divides a power of the base. Part (ii) runs the same proof with `Lambda` in place of `mu` and changes three inputs: the main term, which the principal characters of region C2 carry, Lemma 7.2; the minor-arc input, Lemma 4.2 of Maynard (2022) in place of Basak, Robles and Zaharescu (2023), Lemma 7.4; and the arithmetic input, primes in progressions to moduli dividing a power of the base, where a possible exceptional zero leaves a secondary term that a finite family of characters fixed by the base bounds, so that Siegel's theorem is never used, Lemma 7.5. Two trivial bounds gain a factor `log x`, since `Lambda <= log x` where `|mu| <= 1`. Everything else is identical, the hybrid `l^1` lemma included: the earlier paper proves it under a certificate with the constant `C_F^2 (1 + pi base^2 C_F/fill)`, and part (i) uses it in that form too.

The digit set enters that proof at six places, and the half interval meets each one.

- The certificate, at shift `0` in region A, where the whole `l^1` mass is paid against the uniform `x^(4/5 + eps)`, and at every shift in the hybrid `l^1` lemma of region B, which asks only `alpha_1 < 1/4`. Theorem 4.1 and Proposition 6.2 supply it with `C_F = 3/2` from `94939`.
- Two consecutive digits, in region C1, through the contraction of `|hat F_n|` near a fraction with a denominator prime to the base, Lemma 6.1 of the earlier paper. The half interval keeps `{0, 1}`.
- `fill`, in the contraction constant `rho = 1 - (2/fill)(1 - cos(pi/(4 base)))` of that lemma, which the assembly carries into `C_0` and `c`. It depends on the base alone.
- At most `fill + 1` blocks at each scale in the block split, Lemma 3.2 of the earlier paper, at any number of avoided digits.
- The mass floor `A_F(x) >= fill^(L-1) - 1` of the same lemma, which the half interval meets directly since `0` is in `F`.
- For the prime count only, `kappa_F`, read on the last digit, Lemma 7.2 and Lemma 7.3.

No large sieve and no hybrid estimate is imported: the hybrid `l^1` bound, Lemma 5.1 of the earlier paper, is proved there from the certificate at every shift and `|hat F'| <= pi base (base - 1)`, which holds at every digit set. What the proof imports, none of it reading the digit set: for `mu`, the minor-arc bound of Basak, Robles and Zaharescu (2023), Theorem 1.4, with Theorems 23.5 and 23.6 of Koukoulopoulos (2019) inside its proof, and for the characters of region C2 Exercise 8.4 and Theorem 7.2 of Koukoulopoulos (2019) with Lemma 2.1, Lemma 2.4 and Proposition 2.3 of Chang and Martin (2019); for `Lambda`, Lemma 4.2 of Maynard (2022) and Theorems 12.3, 12.4 and 12.8 of Koukoulopoulos (2019). Effective here means computable from the base, not explicit: the implied constant of the minor-arc bound at the `eps` used, the zero-free constants of the character bounds and the constants of the progression bounds are unstated at their sources, so `C` and `c` are computable and not printed. The bound of the earlier paper on the saving, `c <= sqrt(|log rho|/24)/40` with `|log rho| <= pi^2/(15 fill base^2)`, holds unchanged.

**Lemma 7.2 (the main term).** Let `F` be any digit set, `n >= 1`, `y = base^n` and `P >= 0`, and, with `Lambda(0) = 0`, write the block sum as `sum_(u in D_n) Lambda(Py + u) = y^(-1) sum_(a mod y) hat F_n(a/y) S_P(-a/y)` with `S_P(theta) = sum_(Py <= v < (P+1)y) Lambda(v) e(v theta)`. At a residue `a` of region C2, `a/y = l/d + j/y` with `(l, d) = 1`, `d < Z` dividing a power of `base` and `|j| d = h < Z`, let the principal part of `S_P(-a/y)` be that sum with `Lambda(v)` replaced, on each class `v = r mod d`, by its mean `(d/phi(d)) 1_(gcd(r,d)=1)`. If `rad(base) < Z <= 2^n` and `rad(base) < Q`, with `Q = y^(3/5)` the Dirichlet parameter of the dissection, the principal parts vanish unless `j = 0` and `d` divides `rad(base)`, and over region C2

```
y^(-1) sum_(a in C2) hat F_n(a/y) (principal part of S_P(-a/y))
    = (base/phi(base)) #{u in D_n : gcd(u, base) = 1} = kappa_F fill^n .
```

**Proof.** Every prime power `p^e` exactly dividing `d` has `2^e <= p^e < Z <= 2^n`, so `e < n` and `d` divides `y`. On the class `r mod d`, `0 <= r < d`, the block holds the `y/d` points `v = Py + r + dw`, `0 <= w < y/d`, and `e(-va/y) = e(-rl/d) e(-rj/y) e(-wj/(y/d))`, because `Pyl/d`, `Pj` and `wl` are integers. The sum over `w` is `y/d` when `y/d` divides `j` and `0` otherwise, and `|j| < Z/d < y/d`, so only `j = 0` survives, where the principal part is `(d/phi(d))(y/d) sum_(gcd(r,d)=1) e(-rl/d) = y c_d(l)/phi(d) = y mu(d)/phi(d)`, `c_d` the Ramanujan sum, since `(l, d) = 1`. It vanishes unless `d` is squarefree, and a squarefree `d` dividing a power of `base` divides `rad(base)`. Conversely every `d` dividing `rad(base)` is below `Z`, and each residue `a = ly/d` with `(l, d) = 1` has `l/d` as its fraction, since a second fraction `l'/d'` with `d' <= Q` within `1/(d'Q)` of it would force `1/(dd') <= 1/(d'Q)`, that is `d >= Q`, against `d <= rad(base) < Q`. So the left side is `sum_(d | rad(base)) (mu(d)/phi(d)) sum_(gcd(l,d)=1) hat F_n(l/d) = sum_(u in D_n) sum_(d | rad(base)) mu(d) c_d(u)/phi(d)`. The summand is multiplicative in `d`, so the inner sum is `prod_(p | base) (1 - c_p(u)/(p - 1))`, and `c_p(u)` is `p - 1` when `p` divides `u` and `-1` otherwise, so the product is `base/phi(base)` when `gcd(u, base) = 1` and `0` otherwise. At `n >= 1`, `u` is congruent to its last digit modulo `base`, so `#{u in D_n : gcd(u, base) = 1} = fill^(n-1) #{f in F : gcd(f, base) = 1}`. □

The rest of region C2 is `Lambda` less its mean on each class, and there `Lambda` differs from `mu`: a possible exceptional zero leaves a secondary term in the primes of a class, which Lemma 7.5 below bounds and Proposition 7.6 pays.

**Lemma 7.3 (the constant on the half interval).** At every odd base, `kappa_F = (base/phi(base)) #{f in F : gcd(f, base) = 1}/fill = base/(base + 1)`; at a prime base `p` it is `p/(p + 1)`.

**Proof.** The residues prime to `base` pair as `f` and `base - f`, distinct because `base` is odd, and exactly one of each pair lies in `{1, ..., (base-1)/2}`, which is `F` less `{0}`; `0` is not prime to `base`. So `#{f in F : gcd(f, base) = 1} = phi(base)/2` and `kappa_F = base/(2 fill) = base/(base + 1)`. □

At a prime `p` this reads: `(p-1)/2` of the `(p+1)/2` last digits are prime to `p`, and each such class carries `p/(p-1)` times its share of the primes, so `kappa_F = ((p-1)/(p+1)) (p/(p-1)) = p/(p+1)`.

The rest of this section proves part (ii). As in Sections 3 to 7 of the earlier paper, a block `Py + D_n` is fixed with `x/y < base^K`, where `K = ceil(kappa sqrt(log x))` is the depth of its Proposition 3.3, `kappa` there a depth constant and not `kappa_F`; `x` is large, so `log y >= (log x)/2`; and every residue `a mod y` carries the fraction `l/d` of its Definition 4.1, with `d <= Q`, height `h = |ad - ly| <= y/Q` and `|a/y - l/d| = h/(dy)`. The block sum of Lemma 7.2 is written `Sigma_Lambda(P, n) = sum_(u in D_n) Lambda(Py + u)`, and `S_P` is as there.

**Lemma 7.4 (the minor arcs for `Lambda`).** For `X >= 2` and real `theta = l'/d' + beta` with `(l', d') = 1` and `|beta| < 1/d'^2`,

```
|sum_(v < X) Lambda(v) e(v theta)| << (X^(4/5) + X^(1/2) |d' beta|^(-1/2) + X |d' beta|^(1/2)) (log X)^4 ,
```

with an implied constant that the source does not print; its proof there involves no parameter, and this paper reads it as absolute. Consequently

```
|S_P(a/y)| << x^(4/5) base^(K/5) (log x)^4                          on region A ,
|S_P(a/y)| << (x^(4/5) base^(K/5) + x max(d, h)^(-1/2)) (log x)^4   on region B .
```

**Proof.** The first display is Lemma 4.2 of Maynard (2022), read at source with its proof in arXiv:1510.07711v1, which derives it from Vaughan's identity and its Lemma 4.1 on `sum min(M, ||theta v||^(-1))`. It says nothing at `beta = 0`, so it is read at a coarser fraction. Since `d <= Q`, `h/y <= 1/Q <= 1/d`. First, at `d >= 2`, Dirichlet's theorem at level `floor(d/2)` gives a reduced `l'/d'` with `d' <= d/2` and `|d' a/y - l'| < 2/d`. Since `d' < d` and `(l, d) = 1`, `d' l/d` is not an integer, so `|d' a/y - l'| >= |d' l/d - l'| - d' h/(dy) >= 1/d - 1/(2d)`. So `beta' = a/y - l'/d'` has `1/(2d) <= |d' beta'| < 2/d` and `|beta'| < 2/(d d') <= 1/d'^2`, and the lemma gives `(X^(4/5) + (2Xd)^(1/2) + X (2/d)^(1/2)) (log X)^4`. Second, at `l/d` itself when `h >= 1` and `d < y^(2/5)`: there `|beta| = h/(dy) < 1/d^2`, since `dh < y^(4/5)`, and `|d beta| = h/y`, so the lemma gives `(X^(4/5) + (Xy/h)^(1/2) + X (h/y)^(1/2)) (log X)^4`. `S_P` is the difference of at most two such sums with `X <= 2x`, and `y <= x < y base^K`. On region A, `y^(2/5) <= d <= y^(3/5)`, and the first reading gives `(2Xd)^(1/2) <= 2 x^(4/5)` and `X (2/d)^(1/2) <= 3 x y^(-1/5) < 3 x^(4/5) base^(K/5)`. On region B, `d < y^(2/5)` and `max(d, h) >= Z`. If `d >= h`, then `d >= Z >= 2`, and the first reading gives `(2Xd)^(1/2) <= 2 x^(7/10)` and `X (2/d)^(1/2) <= 3 x d^(-1/2)`; if `h > d`, the second gives `(Xy/h)^(1/2) <= 2 x h^(-1/2)` and `X (h/y)^(1/2) <= 2 x y^(-3/10) < 2 x^(7/10) base^(3K/10)`. In both cases every term but `x max(d, h)^(-1/2)` is at most a constant times `x^(4/5) base^(K/5)`, since `base^K < x`. □

Lemma 7.4 takes the place of Lemma 4.3 of the earlier paper and has its shape. There the bound for `mu` holds at `beta = 0` and is read at `l/d` and at a finer fraction; here it is read at a coarser fraction and at `l/d`, with `(log x)^4` for `(log x)^3` and no `eps`.

**Lemma 7.5 (`Lambda` in progressions to base-smooth moduli).** There are effective `C_2, c_2' > 0`, depending on `base` alone, such that for every `X >= 2`, every `u <= X`, every `d` dividing a power of `base` with `d <= exp(sqrt(log X)/2)` and every residue `r`,

```
|psi(u; d, r) - 1_(gcd(r,d)=1) u/phi(d)| <= C_2 X exp(-c_2' sqrt(log X)) ,
```

where `psi(u; d, r)` sums `Lambda(v)` over `v <= u` with `v = r mod d`.

**Proof.** Take `X` large; below any fixed `X` the bound holds by the choice of `C_2`. If `u <= X exp(-sqrt(log X))`, both terms are at most `u log X`, since `psi(u) <= u log u`. If `gcd(r, d) > 1`, every prime power counted shares a prime with `d`, so it is a power of a prime of `base`, and `psi(u; d, r) <= omega(base) log u`. Otherwise `gcd(r, d) = 1`, `log u >= (log X)/2` and `d <= u`. At `d >= 3`, Theorem 12.4 of Koukoulopoulos (2019), read at source with its proof, gives

```
psi(u; d, r) = (u - chi_1(r) u^(beta_1)/beta_1)/phi(d) + O(u exp(-c_2 sqrt(log u)))
```

uniformly in `u >= d >= 3` and `gcd(r, d) = 1`, with `c_2` absolute, the second term present only when the `L`-functions modulo `d` have the one real exceptional zero `beta_1` of its Theorem 12.3, of a real non-principal character `chi_1`. The display is in the form the proof of Theorem 12.4 produces, p. 122; the statement prints `u^(beta_1)` in place of `u^(beta_1)/beta_1`, and the two can differ by more than the error term, but both are at most `2 u^(1 - c_base)` below, so either serves here. That proof uses the zero-free region, the explicit formula and a zero count, and the book's ineffectivity enters only through Siegel's theorem, which it does not use. `beta_1` is a zero of the real primitive character inducing `chi_1`, whose conductor divides `d` and so divides `8 rad(base)` by Lemma 6.4 of the earlier paper: one of a finite family fixed by the base. Theorem 12.8 of Koukoulopoulos (2019) puts every real zero of that family below `1 - c_base`, with `c_base > 0` effective and fixed by the base, so the secondary term is at most `2 u^(1 - c_base)/phi(d)`. At `d <= 2` the same theorem at modulus `4` gives `psi(u; d, r) = psi(u; 4, 1) + psi(u; 4, 3) + O(log u)`, the error counting the powers of `2`, and the two secondary terms cancel, since `chi_1(1) + chi_1(3) = 0`. With `log u >= (log X)/2`, the error terms are at most a constant times `X exp(-(c_2/2) sqrt(log X))` and `2 X exp(-c_base sqrt(log X))`, so the lemma holds with `c_2' = min(1/2, c_2/2, c_base)`. □

This is Lemma 6.5 of the earlier paper for `Lambda`. There the exceptional zero is met inside a Perron integral for `1/L`; here it survives as a secondary term of the main term, and the same finite family of real characters, Lemma 6.4 there, bounds it. Siegel's theorem is never used, and this is the remark of Maynard (2022) that Siegel zeros play no role at these moduli, written out.

**Proposition 7.6 (region C2 less its main term).** Let `rad(base) < Z <= min(2^n, exp(sqrt(log x)/2))`, the second bound being `C_0 <= 1/2` in `Z = exp(C_0 sqrt(log x))`, and `rad(base) < Q`, and let `N_Z <= (1 + log_2 Z)^(omega(base))` be the number of `d < Z` dividing a power of `base`. Then for `x` large

```
|y^(-1) sum_(a in C2) hat F_n(a/y) S_P(-a/y) - kappa_F fill^n| <= 120 C_2 fill^n base^K Z^2 N_Z exp(-c_2' sqrt(log x)) .
```

**Proof.** By Lemma 7.2 the principal parts over C2 sum to `kappa_F fill^n`, so the left side is at most `y^(-1) sum_(a in C2) |hat F_n(a/y)| |R(a)|`, with `R(a)` the sum `S_P(-a/y)` less its principal part. At a residue of C2, `a = ly/d + j` with `|j| = h/d < Z/d`, as in the proof of Lemma 7.2, and `e(-va/y) = e(-vl/d) e(-vj/y)`. Partial summation against `e(-vj/y)`, whose total variation over the block is at most `2 pi |j|`, gives `|R(a)| <= (1 + 2 pi |j|) max_u |sum_(r mod d) e(-rl/d) E(u; d, r)|` over `Py <= u < (P+1) y`, with the class error

```
E(u; d, r) = sum_(Py <= v <= u, v = r mod d) Lambda(v) - (d/phi(d)) 1_(gcd(r,d)=1) #{Py <= v <= u : v = r mod d} .
```

So `|R(a)| <= (d + 2 pi Z) max |E| <= 8 Z max |E|`. The count of a class up to a point `w` is within `1` of `w/d`, so `|E(u; d, r)|` is at most twice the left side of Lemma 7.5 plus `2 d/phi(d)`, and `d/phi(d) <= base/phi(base)`. Lemma 7.5 at `X = 2x` applies, since `u < 2x` and `d < Z <= exp(sqrt(log x)/2)`, and gives `max |E| <= 4 C_2 x exp(-c_2' sqrt(log x)) + 2 base/phi(base) <= 5 C_2 x exp(-c_2' sqrt(log x))` for `x` large. C2 holds at most `3 Z N_Z` residues, `|hat F_n| <= fill^n` on them, and `x/y < base^K`. □

This is Proposition 6.6 of the earlier paper with the class error of `Lambda` in place of `M(u; d, b)`; the count of C2 and the partial summation are identical.

**Proof of Theorem 7.1(ii).** The block expansion and the split, Lemmas 3.1 and 3.2 of the earlier paper, hold with `Lambda` in place of `mu`, since neither reads the function beyond its being real. The reduction, Proposition 3.3 there, changes in two places. Its block hypothesis becomes `|Sigma_Lambda(P, n) - kappa_F fill^n| <= C' fill^n exp(-c' sqrt(log x))` for every kept block, and the main terms `kappa_F fill^n` of the kept blocks differ from `kappa_F A_F(x)` by at most `kappa_F` times the number of elements of the discarded scales, the point `x` and the element `0`. Since `Lambda <= log x` and `kappa_F <= base/phi(base)`, the discarded scales weigh at most `(log x + base/phi(base)) 6 fill^(1-K) A_F(x)`, and the point `x` and the element `0` at most `log x + 2 kappa_F`; the factor `log x` is absorbed by halving the exponent, and below any fixed `x` the bound `(log x + kappa_F) A_F(x)` is absorbed into `C`. So part (ii) follows from the block bound with `c = min(c', kappa log fill)/2`, as there.

Take

```
C_0 = min(1/2, sqrt(|log rho|/24), c_2'/4) ,   kappa = C_0/(20 log base) ,
```

the choice of the earlier paper with `c_2'` of Lemma 7.5 in place of its `c_1`, so that `base^K <= base exp((C_0/20) sqrt(log x))`, and `rad(base) < Z <= 2^n` and `rad(base) < Q` for `x` large. By the block expansion and Proposition 7.6, `|Sigma_Lambda(P, n) - kappa_F fill^n|` is at most the bound of Proposition 7.6 plus `y^(-1) sum |hat F_n(a/y)| |S_P(a/y)|` over each of A, B and C1, each `fill^n` times a saving.

- Region A: the certificate at `s = 0` gives `sum_(a mod y) |hat F_n(a/y)| <= C_F fill^n y^(alpha_1)`, and Lemma 7.4 with `x^(4/5) < y^(4/5) base^(4K/5)` gives `<< C_F fill^n y^(alpha_1 - 1/5) base^K (log x)^4`, a power saving since `log y >= (log x)/2`. This is Proposition 4.4 of the earlier paper with no `eps`.
- Region B: Proposition 5.2 of the earlier paper with Lemma 7.4 in place of its Lemma 4.3. Its classes, `D <= d < 2D` and `H <= h < 2H` or `h = 0`, number at most `(log_2 y + 2)^2`, and on each `max(d, h) >= max(D, H) > Z/2`; its Lemma 5.1, under the certificate, bounds the mass of a class by `<< fill^n max(D, H)^(2 alpha_1)`, all as there. With Lemma 7.4 a class weighs `<< y^(-1) (x^(4/5) base^(K/5) + x max(D, H)^(-1/2)) (log x)^4 fill^n max(D, H)^(2 alpha_1)`. With `max(D, H) <= y^(2/5)` and `x^(4/5) < y^(4/5) base^(4K/5)` the first term is `<< fill^n base^K y^((4/5)(alpha_1 - 1/4)) (log x)^4`; with `x/y < base^K` and `2 alpha_1 < 1/2` the second is `<< fill^n base^K (log x)^4 Z^(-(1/2 - 2 alpha_1))`. Over the classes, region B is `<< fill^n base^K (log x)^6 (y^((4/5)(alpha_1 - 1/4)) + Z^(-(1/2 - 2 alpha_1)))`. It asks `alpha_1 < 1/4`, and `1/2 - 2 alpha_1 > 1/10` under a certificate below `1/5`.
- Region C1: Proposition 6.3 of the earlier paper, through its Lemma 6.1 on the two consecutive digits, with the trivial bound `|S_P| <= y log(2x)` in place of `|S_P| <= y`, so the bound gains the factor `log(2x)`: `3 rho^(-1) fill^n log(2x) exp((2 C_0 - |log rho|/(6 C_0)) sqrt(log x))`.
- Region C2 less its main term: Proposition 7.6, a saving `exp(-(c_2' - 2 C_0 - C_0/20) sqrt(log x))` up to the factor `N_Z`, a power of `log x`, and `c_2' >= 4 C_0` makes that exponent at least `(2 - 1/20) C_0`.

As in the assembly of the earlier paper, every saving is at least `exp(-(C_0/20) sqrt(log x))` up to powers of `log x` and constants depending on `base` and `F`, so the block bound holds with `c' = C_0/21`, and part (ii) follows with `c = min(C_0/21, kappa log fill)/2`. Every constant is effective: the certificate's `C_F`, the implied constant of Lemma 7.4, the constants `c_2`, `c_base` and `C_2` of Lemma 7.5, and every threshold on `x` above. The bound `c <= sqrt(|log rho|/24)/40` of the earlier paper holds for `Lambda` too. □

**Proof of Theorem 1.1.** At an odd base `base >= 94939` the half interval contains the consecutive digits `0` and `1`, and by Theorem 4.1 and Proposition 6.2 it carries the shifted-grid certificate `C_F = 3/2`, `alpha_1 = log_base(lambda/fill) < 1/5`. Theorem 7.1 gives both bounds, with `kappa_F` as stated; `F` is determined by `base`, so `C` and `c` depend on `base` alone. □

**Proof of Corollary 1.2.** By Lemma 2.2, at the prime base `p` the set `S_F` is `{k >= 1 : p does not divide C(2k,k)}`, so `A(x) = A_F(x)`, the first sum is `M_F(x)`, and the second is the von Mangoldt sum of Theorem 1.1, where `kappa_F = p/(p+1)` by Lemma 7.3. □

**Proof of Theorem 1.3.** Theorem 5.1 of [The First Base Below a Quarter](first-base-below-a-quarter.md) proves, under the stated hypothesis and for any digit set `F` with `fill >= 2`, dimension `alpha_base` and `l^1` exponent `alpha_1`, that `|M_F(x)| <<_(base, F, eps) A_F(x)^(1 - delta + eps)` with `delta = (1/4 - alpha_1)/alpha_base`; its only Mobius input is the uniform bound `max_theta |sum_(v <= x) mu(v) e(v theta)| << x^(3/4 + eps)` of Baker and Harman (1991), and its `alpha_1` is the growth exponent of `fill^(-i) max_s Sigma_i(s)`. Theorem 4.1 bounds that exponent by `log_base(lambda/fill)`, and a larger `alpha_1` only weakens the conclusion, so the theorem holds with `alpha_1 = log_base(lambda/fill)`, which is below `1/4` at every odd base from `3789` by Proposition 6.2. As the base grows, `alpha_1 <= log((2/pi) log base + 2.64)/log base` tends to `0` by Proposition 5.1 and `alpha_base = log((base + 1)/2)/log base` tends to `1`, so `delta` tends to `1/4`. □

## What is in print

Every source below is in the references and read at source. The wider literature on missing-digit sets, with whole-text searches of the circle-method sources for a Mobius or Mertens sum, is Section 9 of [the earlier paper](../../coprime/papers/unconditional-mertens-at-large-base.md), and none of the sources read there or here carries a Mobius or Mertens sum over a missing-digit set.

Maynard (2022), Theorem 1.3, read in arXiv:1510.07711v1, proves `sum_(n < q^k) Lambda(n) 1_B(n)` asymptotic with error `O_A((q - s)^k (log q^k)^(-A))` when the `s` excluded digits are consecutive and `q - s >= q^(4/5 + eps)`, `q` sufficiently large in terms of `eps`, by a sketch in its Section 9. The half interval has `q - s = (q + 1)/2`, inside that range, so the prime asymptotic on it at sufficiently large base, with an unquantified base and a log-power error, is the content of that theorem, and at a prime base so is the prime count on `{k : q does not divide C(2k,k)}`. The constant of the sketch, `alpha_(q,s) = log((2 + 2/log q)(q/(q - s)) log q)/log q`, clears `1/5` on the half interval only from `q = 7777884825`, the same reading giving the one-missing-digit crossing `1520573` of the constant of its Lemmas 5.1 and 5.3 (`lab/py/interval-digits`, verb `wall`); in the normalisation of Section 5 it pays `4 log base` per digit where Theorem 4.1 pays `(2/pi) log base`, and the whole gain in the base lies there. The theorem prints its main term as `q(phi(q) - s')/((q - 1) phi(q)) (q - s)^k`, `s'` the excluded digits prime to `q`, which at the half interval of a prime `q = p` is `p/(2(p - 1))` times the count, against the `p/(p + 1)` of Lemma 7.3. With `q - s` in place of `q - 1`, the form of the constant `kappa_B = q(phi(q) - t)/(phi(q)(q - s))` of Maynard (2019), `t` the excluded digits prime to `q`, and of the instruction in Section 9 of Maynard (2022) to replace `q - 1` by `q - s`, it is `p/(p + 1)`; so the printed form reads as a misprint against the 2019 constant, and Lemma 7.2 gives the `q - s` form. After its Theorem 1.1, Maynard (2022) remarks that its estimates are used only at highly composite moduli, where Siegel zeros play no role, so its error terms could be replaced by effective ones of size `O((q - 1)^k exp(-c k^(1/2)))`, and its Section 9 extends the arguments to Theorem 1.3 without restating the remark: for `Lambda` on the half interval the effective shape is read from that remark, not stated at source.

Maynard (2019), Theorem 1.2, read in arXiv:1604.01041v2, gives the order of magnitude `X^(log(q - s)/log q)/log X` of the primes avoiding the top block `{q - s, ..., q - 1}` at large `q` and `s <= q - q^(57/80)`, a range that holds the half interval, and remarks the asymptotic for the primes avoiding the bottom block `{0, ..., s - 1}` at `s <= q - q^(3/4 + delta)`, with an `o(1)` error. Neither paper carries a Mobius sum.

So for `Lambda` the asymptotic on the half interval is in print at sufficiently large base, and its effective shape is a reading of a remark made there at one missing digit; what this paper adds is the written proof of that shape, Section 7 from Lemma 7.2 on, the Mobius bound, which no source read here carries, and the explicit base `94939` for both sums.

## The check and the limits

The certificate is a chain of elementary inequalities, and the links named at the end of Section 4 are recomputed on grids small enough to hold, where a wrong sign or a wrong count would show.

**Fact 9.1 (the falsification).** At every odd base `3..401` and at `1001`, `10001` and `100001`, on `201` shifts of `[0, 1/2]`, the three sums of the proof of Theorem 4.1 stay under their bounds: the harmonic parts of `Sigma_P` and `U`, in units of `base/pi`, under `log(base + 3) + gamma + K_2 t^2` and `log(base + 4) - 0.0757`, the first by at least `3.9999 * 10^(-5)`, and `Sigma_Q` under `(base/pi) log tan(3 pi/8 + pi/(4 base))`. At every odd base `3..401` on `801` shifts and at `1001`, `4001`, `10001`, `30001` and `100001` on `401`, the closed form of Lemma 3.3 meets the direct sum to relative `10^(-9)`, `G` stays under (4.2), at most `0.999857` of it and tight at `t = 0` by design, `T phi` stays under `lambda phi`, at most `0.999866` of it, and `G` stays over the lower bound of Proposition 5.2, at least `1.534638` times it. The grid sums `Sigma_i(s)` at `40` digits, at the nine bases `3, 5, 7, 9, 11, 13, 15, 17, 21`, the four shifts `0`, `1/2`, `1/(2 base)` and `(sqrt 5 - 1)/2`, and every level with `base^i <= 9261`, from level `8` at base `3` to level `3` at base `21`, stay under `(3/2) lambda^i`, at most `0.557408` of it, and the identity of Lemma 3.2 holds at level `2` at bases `5`, `7` and `9` to relative `10^(-15)`. Every ratio is printed with directed rounding, upper bounds rounded up and lower bounds down (`lab/py/interval-digits`, verb `check`, 13 seconds).

The wall of Proposition 6.2 is where this certificate clears the bar, and the readings say the route itself stops not far below.

**Conjecture 9.2 (the true rate).** The growth rate `rho_T` of `T` on the half interval, the spectral radius that `max_s Sigma_i(s)` grows at, reads `rho_T/fill = (2/pi) log base + 2.26` near base `7 * 10^4` and meets `base^(1/5)` between `70001` and `80001`; the one-step constant reads `max G/fill = (2 sqrt 2/pi) log base + 1.19`, at `t = 1/2`, and clears `1/4` from `7075` and `1/5` from `317063`.

Evidence. Power iteration of `T` with linear interpolation reads `rho_T/fill - (2/pi) log base` as `2.00766`, `2.14698` and `2.21757` at `101`, `1001` and `10001` on `1000` cells, and `2.25718`, `2.26007` and `2.26252` at `60001`, `70001` and `80001` on `300` cells, where `base^(1/5) - rho_T/fill` reads `-0.2325`, `-0.0508` and `0.1137`, and the iterated weight stays within a factor `0.81` of its maximum. The one-step reading `G(1/2)/fill - (2 sqrt 2/pi) log base` is `1.19055` at `7075` and `1.19171` at `317063`, and the maximum of `G` over `401` shifts sits at `t = 1/2` at every base checked from `1001` to `100001` (`lab/py/interval-digits`, verbs `rate` and `check`). Failure modes. The iteration interpolates and rounds, and bounds nothing: its low and high ratios agreeing to five digits says it converged on its grid, not that the grid resolves `T`. The constant beside `(2/pi) log base` still rises at `80001` and may not settle at `2.26`, and the crossing is read at three bases. The one-step scan reads `G` at one point. If the conjecture holds, then since `Sigma_i(0) = (T^i 1)(0)` and the weight is bounded below, the unshifted mass of region A grows at `rho_T` too, so no shifted-grid certificate of any kind carries this route below the crossing.

**Proposition 9.3 (uniform square-root cancellation implies RH).** If for every `eps > 0` some `C_eps` gives `|M_F(x)| <= C_eps A_F(x)^(1/2 + eps)` on the half interval at every odd base and every `x >= 1`, the Riemann hypothesis holds.

**Proof.** `S_F` contains every integer below `fill`, so at the odd base `2 ceil(x) + 1` the sum `M_F(x)` is the Mertens function `M(x)` and `A_F(x) = floor(x)`. So `M(x) = O(x^(1/2 + eps))` for every `eps > 0`, which is the Riemann hypothesis by Theorem 14.25 (C) of Titchmarsh (1986). □

The implication is trivial and the digits play no part in it: the range `x < fill` carries all of it, and there the statement is the Riemann hypothesis itself. No converse is claimed.

What the theorem does not say. The saving is `exp(-c sqrt(log x))` with a tiny `c`, never a power, and `C` and `c` are computable but not printed. Below `94939` nothing is proved without a hypothesis, and below `3789` nothing under one. The theorem says nothing about the true size of `M_F(x)`, for which square-root cancellation against `A_F(x)` is the natural guess, and Proposition 9.3 shows that guess, asked uniformly over the bases, is at least as strong as the Riemann hypothesis. It says nothing about the even-digit design `2 S_F` of Remark 2.3.

## Open problems

Three things are left undone. The even-digit design `2 S_F` of Remark 2.3 has as its Mobius sum `-sum mu(k)` over the odd `k` of `S_F` up to `x/2`, and at an odd base the parity of `k` is that of its digit sum, so that sum is half of `M_F(x/2)` less its twist by `e(k/2)`; the twist reads the certificate at the shift `1/2`, which Theorem 4.1 covers, but near a fraction `l/d` with `d` dividing a power of the base it sees the transform at `l/d + 1/2`, so regions C1 and C2 must be rerun with `2 base` in place of `base`, and that is not written. No proof is written below `94939`: a sharper weight than `1 + |sin(pi t)|/2` moves the wall toward the crossing of the true rate, which Conjecture 9.2 puts between `70001` and `80001`, and if the readings are right no certificate passes that crossing, leaving the exponent `4/5` of the minor-arc input as the other lever, which neither this paper nor Maynard (2022) moves. And the growth rate of `T` itself is read, not proved: the constant beside `(2/pi) log base` in `rho_T/fill` lies between `-1.31` and `2.60043004 + 3.4/base` by Section 5 and reads near `2.26`, and whether it settles is open.

## Reproducibility

One study prints every number of Sections 4, 5, 6, 8 and 9, run from the repository root with one verb and raising if any check fails. `uv run python research/mobius/lab/py/interval-digits/interval.py wall`, under a second, prints Propositions 5.1, 5.2 and 6.2, the constant `c_inf` and the crossing `7777884825` of the Section 9 constant of Maynard (2022) with its calibration `1520573`; `check`, in about 13 seconds, prints Fact 9.1 and the position of the maximum of `G`; `rate`, in about 20 seconds, prints the readings of Conjecture 9.2. The study evaluates every closed form at `120` bits in `mpmath` interval arithmetic, lower bounds truncated down and upper bounds rounded up; `check` runs in float64 with margins far above float error and prints its ratios with directed rounding, except at `t = 0`, where (4.2) is tight to `1.5 * 10^(-4)` of its size by design, and its level sums at `40` digits; `rate` prints readings that bound nothing. The figure is `figures/paper-half-interval-mobius.ts`, drawn in the browser in under half a second a theme; its census writer `research/shared/census/src/bin/paper-half-interval-mobius.rs` recomputes the closed form (4.1) in floating point at about `1600` odd bases from `101` to `10^7 + 1` and at the four bases beside the walls, and asserts the walls `94939` and `3789`, the sign of both margins at every base it draws, the margins `1.3678 * 10^(-8)` and `7.9625 * 10^(-6)` to five digits and `alpha_1` at `100003`. Theorem 7.1 prints no number. Part (i) is reused, and the numbers of its dissection are printed by the studies named in the earlier paper. The proof of part (ii) computes no constant; its main-term identity and its four regions for `Lambda` are checked on small designs far below any wall by `uv run python research/coprime/lab/py/prime-dissection/primes.py series` and `regions`, each in about a second, a check of the bookkeeping and not of the theorem.

## References

- Maynard 2022, Primes and polynomials with restricted digits, Int. Math. Res. Not. 2022, 10626-10648, read in arXiv:1510.07711v1. [doi.org/10.1093/imrn/rnab002](https://doi.org/10.1093/imrn/rnab002)
- Maynard 2019, Primes with restricted digits, Invent. Math. 217, 127-218, read in arXiv:1604.01041v2. [link.springer.com](https://link.springer.com/article/10.1007/s00222-019-00865-6)
- Basak, Robles and Zaharescu 2023, Exponential sums over Mobius convolutions with applications to partitions. [arxiv.org/abs/2312.17435](https://arxiv.org/abs/2312.17435)
- Koukoulopoulos 2019, The Distribution of Prime Numbers, Graduate Studies in Mathematics 203, American Mathematical Society, doi:10.1090/gsm/203, read in the author's preliminary version. [dms.umontreal.ca](https://dms.umontreal.ca/~koukoulo/documents/publications/primes.pdf)
- Chang and Martin 2019, The smallest invariant factor of the multiplicative group. [arxiv.org/abs/1908.00035](https://arxiv.org/abs/1908.00035)
- Baker and Harman 1991, Exponential sums formed with the Mobius function, J. London Math. Soc. (2) 43, 193-198. [doi.org/10.1112/jlms/s2-43.2.193](https://doi.org/10.1112/jlms/s2-43.2.193)
- Titchmarsh 1986, The Theory of the Riemann Zeta-Function, second edition revised by D. R. Heath-Brown, Clarendon Press, Oxford, ISBN 978-0-19-853369-6, Theorem 14.25 (C); the link is the publisher's record. [global.oup.com](https://global.oup.com/academic/product/the-theory-of-the-riemann-zeta-function-9780198533696)
- Kummer 1852, Uber die Erganzungssatze zu den allgemeinen Reciprocitatsgesetzen, J. reine angew. Math. 44, 93-146, doi:10.1515/crll.1852.44.93; the link is the record with the public-domain scan. [eudml.org/doc/147500](https://eudml.org/doc/147500)
