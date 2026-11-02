---
title: The Franel Threshold Is a Mean Square
lead: On the Farey fractions whose denominators lie in a digit design, Franel's square-root threshold is equivalent to square-root cancellation, in mean square over the dilates of the design, of their Mertens sums, with explicit logarithms in the sandwich; the full set recovers Franel.
date: 2026-09-21
figure: paper-franel-converse
---

Take the reduced fractions `a/b` in `(0, 1]` with `b <= Q` and `b` written, in base `q`, with digits from a fixed set `F` of at least two digits, and record how far the `j`-th of them sits from `j/m`, `m` their number. On the full set Franel proved in 1924 that the sum of the squares of those distances is `O(Q^(-1+eps))` for every `eps > 0` exactly when the Riemann hypothesis holds, through the Mertens function `M(x)`. On a design the identity behind that theorem splits: its Fourier side is a gcd-weighted quadratic form in one Mertens sum for every dilate `{c : dc in S_F}` of the design, sets that collapse to the whole line only when every digit is allowed. Jordan's totient writes the form as a sum of squares, Schur's test bounds the change of variables by `(1 + ln Q)^2` both ways, and the threshold `sum_j delta_j^2 = O(Q^(-1+eps))` becomes equivalent to `sum_(d <= Q) M_F(Q/d; d)^2 = O(Q^(alpha+eps))`, `alpha = log |F| / log q`; a primorial witness shows that no constant replaces the logarithms. The masses of the dilates sum to `O(Q^(alpha+eps))` with no Mobius input, so square-root cancellation of `mu` in each dilate against its own mass gives the threshold on every design; a hypothesis with an explicit rate in `d` also suffices and is not shown to be needed, and the equidistribution law behind its rate fails at the repunits of base 3 `{0,1}`. Every identity is checked as an exact rational at base 3 `{0,1}` to `Q = 3^6`, in floating point to `Q = 3^8`, and on the full set at `Q = 40, 81, 243`.

## Introduction

![The Farey fractions with denominators in the base 3 digit design {0,1} at Q = 81, 241 ticks on the unit line, over the 2020 ticks of the full Farey sequence of order 81, each tick as tall as its denominator is small; beneath them the sawtooth delta_j = rho_j - j/m of the design, dots joined, against the full set's sawtooth on the same scale.](paper-franel-converse)

Write every whole number in base `3` and keep the ones that use only the digits `0` and `1`: `1, 3, 4, 9, 10, 12, 13, 27, ...`, sixteen of them up to `81`. Those are the denominators. The fractions `a/b` with such a `b`, `1 <= a <= b` and `gcd(a, b) = 1`, sorted, are the top row of the figure, 241 of them at `Q = 81`; the middle row is the full [Farey sequence](https://en.wikipedia.org/wiki/Farey_sequence) of order `81`, 2020 fractions. Number the design's fractions `rho_1 < ... < rho_m` and set `delta_j = rho_j - j/m`; the orange sawtooth of the bottom row is `delta_j` against `rho_j`, the full set's sawtooth behind it in grey on the same scale. Both wander a little above and below zero and both come back. How fast they come back, in the square mean, is the question.

For the full set the answer is Franel's theorem (1924): `sum_j delta_j^2 = O(Q^(-1+eps))` for every `eps > 0` if and only if the Riemann hypothesis holds; Landau (1924), in the note printed after Franel's, gave the same for `sum_j |delta_j| = O(Q^(1/2+eps))`, and Edwards reproduces both in section 12.2 through the Mertens function `M(x) = sum_(n <= x) mu(n)`. The mechanism is an identity: the sum of squares is, up to a constant and a factor `12 m`, the form `sum_(d, e) (gcd(d, e)^2/(d e)) M(Q/d) M(Q/e)`, which is `O(Q^(1+eps))` under `M(x) = O(x^(1/2+eps))` and bounds `M(Q)^2` from above. On a denominator restricted to a digit design [the Franel note](../notes/franel.md) of this tree finds, and this paper proves, that the two directions no longer meet at one function: the form becomes `sum_(d, e) (gcd(d, e)^2/(d e)) M_F(Q/d; d) M_F(Q/e; e)` with `M_F(x; d) = sum_(c <= x, dc in S_F) mu(c)` the Mertens sum of the dilate `{c : dc in S_F}`, a different set for every `d`. The `d = 1` term is the design's own `M_F(Q)`, so the threshold still forces `M_F(Q) = O(Q^(alpha/2+eps))`; the converse needs every dilate.

The theorem of this paper is that the converse needs them in mean square and in nothing finer: for every digit design with at least two digits, `sum_j delta_j^2 = O(Q^(-1+eps))` on the denominator-restricted set if and only if `sum_(d <= Q) M_F(Q/d; d)^2 = O(Q^(alpha+eps))`, the form lying between `(6/pi^2)(1 + ln Q)^(-2)` and `(1 + ln Q)^2` times the mean square for every design and every `Q >= 1` (Theorem 4.3), with no constant in place of the logarithms (Proposition 4.4). On the full set every dilate is the whole line and the statement is Franel's with `M(x) = O(x^(1/2+eps))` standing for the hypothesis (Corollary 4.6). Section 2 sets up the objects; Section 3 proves the three identities; Section 4 the sandwich and the equivalence; Section 5 shows that the masses of the dilates sum to `O(Q^(alpha+eps))` with no Mobius input, so that square-root cancellation of `mu` in each dilate against its own mass gives the threshold, and that a hypothesis with an explicit rate in `d` suffices as well; Section 6 proves that the equidistribution law behind that rate fails at the repunits; Section 7 prints the readings; Section 8 says what is open.

The literature is positioned as read here, four texts: Franel's identity is stated over the rationals in its Fourier form by Huxley (1971, Lemma (10)) and with weights by Huxley (2012, Theorem 1), extended by Kanemitsu and Yoshimoto (1996, Theorem 3) and surveyed by Cobeli and Zaharescu (2003); none restricts the denominators to a digit design, the older literature on Farey fractions with restricted denominators is not surveyed, and this paper claims the theorems it proves and no priority beyond them. The divisor form of Ramanujan's sum is equation (2.7) of Ramanujan (1918); the factorization of a gcd matrix through a totient is Smith's (1875), and Lemma 4.1 is it with Jordan's totient; the operator bound is Schur's (1911); the equivalence of `M(x) = O(x^(1/2+eps))` with the Riemann hypothesis is Theorem 14.25 (C) of Titchmarsh (1986).

## Definitions

**Definition 2.1 (design).** Fix a base `q >= 2` and a set `F` of digits, `F` a subset of `{0, ..., q-1}` with `|F| >= 2`. The digit design `S_F` is the set of whole numbers whose every base-`q` digit lies in `F`, with `0` in `S_F`. Its dimension is `alpha = log |F| / log q`, in `(0, 1]`, and `A_F(Q) = #{n in S_F : 1 <= n <= Q}`. The full set is `F = {0, ..., q-1}`, where `S_F` is every whole number and `alpha = 1`.

**Lemma 2.2.** For `Q >= q`, `Q^alpha / |F|^2 < A_F(Q) < 2 |F| Q^alpha`.

**Proof.** Let `f = |F|` and `L = floor(log_q Q)`, so `q^L <= Q < q^(L+1)` and `L >= 1`. Every member of `S_F` in `[1, Q]` has between `1` and `L + 1` digits, all in `F`, and there are at most `f^k` such strings of length `k`, so `A_F(Q) <= sum_(k = 1)^(L+1) f^k = (f^(L+2) - f)/(f - 1) < (f^2/(f - 1)) f^L <= 2 f q^(L alpha) <= 2 f Q^alpha`, using `f/(f - 1) <= 2`. Some digit of `F` is nonzero, so the strings of exactly `L` digits with that leading digit and the other `L - 1` digits free in `F` are at least `f^(L-1)` members of `S_F`, each in `[q^(L-1), q^L)`, hence at most `Q`; and `f^(L-1) = q^((L+1) alpha) / f^2 > Q^alpha / f^2`. □

The upper bound holds with `|F|` in place of `2|F|` when `0 in F`, where every member is a padded string of `L + 1` digits, and not otherwise: at base 4 with `F = {1, 2}` and `Q = 41`, `A_F(41) = 13` against `|F| Q^alpha = 12.81`.

**Definition 2.3 (denominator set and sawtooth).** For `Q >= 1` the denominator-restricted Farey set is `F_Q^d(S_F) = {a/b : b in S_F, 1 <= b <= Q, 1 <= a <= b, gcd(a, b) = 1}`, the `d` marking the denominator as in the note, where `F_Q(S_F)` is the strict set with the numerator in the design too, of size `m_F(Q) = sum_(b in S_F, b <= Q) phi(b)` with `phi` [Euler's totient](https://en.wikipedia.org/wiki/Euler's_totient_function). Its nodes ascending are `rho_1 < ... < rho_m`, `m = m_F(Q)`, and `delta_j = rho_j - j/m`. Its exponential sum at the integer frequency `k` is `S_F(k, Q) = sum_(j) e(k rho_j)` with `e(x) = exp(2 pi i x)`.

**Lemma 2.4.** `phi(n) > n / (e^2 (ln n)^2)` for `n >= 16`, and for `Q >= 16 q^2`, `Q^(1+alpha) / (e^2 q^2 |F|^2 (ln Q)^2) < m_F(Q) < 2 |F| Q^(1+alpha)`.

**Proof.** Put `y = (ln n)^2 >= 2`. Over the primes `p <= y` dividing `n`, `p/(p-1) <= exp(1/(p-1))` and the values `p - 1` are distinct integers in `[1, y]`, so their product is at most `exp(sum_(k <= y) 1/k) <= exp(1 + ln y) = e y`. The primes `p > y` dividing `n` number at most `ln n / ln y`, since their product is at most `n`, and each contributes `p/(p-1) <= exp(1/(p-1)) <= exp(2/y)`, so together at most `exp(2 ln n / (y ln y)) = exp(1/(ln n ln ln n)) < e`. Hence `n/phi(n) = prod_(p | n) p/(p-1) < e^2 (ln n)^2`. For the node count, `m_F(Q) <= sum_(b in S_F, b <= Q) b <= Q A_F(Q) < 2 |F| Q^(1+alpha)` by Lemma 2.2; and the `|F|^(L-1) > Q^alpha / |F|^2` members of exactly `L = floor(log_q Q)` digits found there each exceed `q^(L-1) >= Q/q^2 >= 16`, so each has `phi(b) > (Q/q^2)/(e^2 (ln Q)^2)`. □

**Definition 2.5 (dilated Mertens sums).** For `d >= 1` the dilate of the design is `d^(-1) S_F = {c >= 1 : dc in S_F}`, its counting function `N_F(Q; d) = #{c : 1 <= c <= Q/d, dc in S_F} = #{n in S_F : 1 <= n <= Q, d | n}`, and its Mertens sum `M_F(x; d) = sum_(1 <= c <= x, dc in S_F) mu(c)`, with `mu` the [Mobius function](https://en.wikipedia.org/wiki/M%C3%B6bius_function). At `d = 1` it is the design's own Mertens sum `M_F(x) = sum_(n in S_F, n <= x) mu(n)` of the [Mobius page](../../mobius/notes/mobius.md). Throughout, `x_d = M_F(Q/d; d)` for `1 <= d <= Q`, a vector of integers with `|x_d| <= N_F(Q; d)` termwise, and `x_d = 0` for `d > Q`. On the full set every dilate is the whole line and `x_d = M(floor(Q/d))`.

For a design carrying `0` a dilate is a regular language, and the Franel note proves it: long multiplication by `d`, read from the least significant digit, is a deterministic automaton on the `d` carries, from carry `r` the digit `e` writing the output digit `(de + r) mod q`, which must lie in `F`, and moving to the carry `floor((de + r)/q)`, with acceptance when the terminal carry lies in `Acc_d = {0} union (S_F intersect [1, d))`; every column of its transfer matrix sums to `|F|`, so a dilate carries the design's own dimension. Nothing below uses the automaton beyond one consequence, which is elementary:

**Lemma 2.6 (base powers are free).** If `0 in F` then `qc in S_F` if and only if `c in S_F`, so `M_F(x; q^j d) = M_F(x; d)` and `N_F(q^j Q; q^j d) = N_F(Q; d)` for every `j >= 0`.

**Proof.** Appending or removing a zero digit neither enters nor leaves `S_F`. □

## The identities

**Lemma 3.1 (the divisor form of Ramanujan's sum).** For `b, m >= 1`, `c_b(m) = sum_(a mod b, gcd(a, b) = 1) e(ma/b) = sum_(d | gcd(m, b)) d mu(b/d)`.

**Proof.** Insert `sum_(g | gcd(a, b)) mu(g) = [gcd(a, b) = 1]` and write `a = g a'`: `c_b(m) = sum_(g | b) mu(g) sum_(a' mod b/g) e(m a' / (b/g))`. The inner sum is `b/g` when `(b/g) | m` and `0` otherwise; substitute `d = b/g`. □

This is equation (2.7) of Ramanujan (1918), credited there to Landau's Handbuch, page 577; the proof is included so that nothing rests on the citation.

**Theorem 3.2 (the exponential sum at frequency `m`).** For every `m >= 1` and `Q >= 1`, `S_F(m, Q) = sum_(d | m) d M_F(Q/d; d)`. In particular `S_F(1, Q) = M_F(Q)`.

**Proof.** Group the nodes by denominator: `S_F(m, Q) = sum_(b in S_F, b <= Q) c_b(m)`, and `c_b(m)` depends on `b` alone. By Lemma 3.1 this is `sum_(b in S_F, b <= Q) sum_(d | m, d | b) d mu(b/d) = sum_(d | m) d sum_(c <= Q/d, dc in S_F) mu(c)`, writing `b = dc`. □

**Theorem 3.3 (the Fourier form).** With `x_d = M_F(Q/d; d)` and the form `G_F(Q) = sum_(d, e >= 1) (gcd(d, e)^2/(d e)) x_d x_e`, a finite sum of rationals, `sum_(k != 0) |S_F(k, Q)|^2 / k^2 = (pi^2/3) G_F(Q)`.

**Proof.** The node set is carried to itself by `r -> 1 - r` away from the node `1`, which `e(k r)` sends to `1`, so `S_F(k, Q)` is real and `S_F(-k, Q) = S_F(k, Q)`. By Theorem 3.2, for `k >= 1`, `S_F(k, Q)^2 = sum_(d | k) sum_(e | k) d e x_d x_e`, and summing `k^(-2)` over the `k >= 1` divisible by both `d` and `e`, that is by `lcm(d, e)`, gives `zeta(2)/lcm(d, e)^2`; the exchange is over finitely many `d, e <= Q` and an absolutely convergent series in `k`. So `sum_(k != 0) |S_F(k, Q)|^2/k^2 = 2 zeta(2) sum_(d, e) (d e/lcm(d, e)^2) x_d x_e`, and `d e / lcm(d, e)^2 = gcd(d, e)^2/(d e)` with `2 zeta(2) = pi^2/3`. □

**Theorem 3.4 (the rank form).** Let `rho_1 < ... < rho_m` be any `m >= 1` points of `(0, 1]`, `S(k) = sum_j e(k rho_j)`, `delta_j = rho_j - j/m` and `s = sum_j rho_j - m/2`. Then `sum_(k != 0) |S(k)|^2 / k^2 = 4 pi^2 (m sum_j delta_j^2 + s(1 - s) - 1/6)`.

**Proof.** Let `B(x) = {x} - 1/2` be the sawtooth, with the Fourier series `B(x) = -sum_(k != 0) e(kx)/(2 pi i k)` in `L^2`. For `u` in `[0, 1)`, `{u + rho_j} = u + rho_j - [rho_j >= 1 - u]`, so with `A(v) = #{j : rho_j < v}` and `D(v) = A(v) - m v`, `sum_j B(u + rho_j) = m u + sum_j rho_j - m/2 - (m - A(1 - u)) = D(1 - u) + s`. Parseval on `[0, 1]` reads `sum_(k != 0) |S(k)|^2/(4 pi^2 k^2) = int_0^1 (D(v) + s)^2 dv = int_0^1 D^2 - s^2`, because `int_0^1 A(v) dv = sum_j (1 - rho_j)` gives `int_0^1 D = -s`. Put `rho_0 = 0`, `rho_(m+1) = 1` and `u_j = m rho_j - j`, so `u_0 = 0` and `u_(m+1) = -1`. On `(rho_j, rho_(j+1)]` the function `A` is `j`, and `int_(rho_j)^(rho_(j+1)) (j - mv)^2 dv = ((j - m rho_j)^3 - (j - m rho_(j+1))^3)/(3m) = (-u_j^3 + (u_(j+1) + 1)^3)/(3m)`. Summing over `j = 0, ..., m`, the term `(u_(m+1) + 1)^3` vanishes, `u_0^3` vanishes, and what is left is `(1/(3m)) sum_(j = 1)^m (3 u_j^2 + 3 u_j + 1) = m sum_j delta_j^2 + sum_j delta_j + 1/3`, using `u_j = m delta_j`. Finally `sum_j delta_j = sum_j rho_j - (m+1)/2 = s - 1/2`, so `int_0^1 D^2 - s^2 = m sum_j delta_j^2 + s(1 - s) - 1/6`. □

**Corollary 3.5.** On the denominator set with `m_F(Q) >= 1`, `G_F(Q) = 12 m_F(Q) sum_j delta_j^2 + 1` when `1 in F` and `G_F(Q) = 12 m_F(Q) sum_j delta_j^2 - 2` when `1` is not in `F`. In either case `G_F(Q) >= 0` and `G_F(Q) = 12 m_F(Q) sum_j delta_j^2 + O(1)`.

**Proof.** The node `1 = 1/1` is present exactly when `1 in S_F`, that is when `1 in F`; for `b >= 2` the nodes `a/b` with `a < b` pair under `a -> b - a` and sum to `phi(b)/2`. So `s = 1/2` when `1 in F` and `s = 0` otherwise, and Theorems 3.3 and 3.4 give `(pi^2/3) G_F(Q) = 4 pi^2 (m sum_j delta_j^2 + 1/12)` or `4 pi^2 (m sum_j delta_j^2 - 1/6)`. The left side of Theorem 3.3 is a sum of nonnegative terms. □

The constant is read as an exact rational at every `Q <= 40` at which the node set changes: `1` on base 3 `{0,1}` and base 10 without the digit 9, `-2` on base 3 `{0,2}`, base 4 `{0,2,3}`, base 5 `{0,2,4}` and base 10 `{0,2,5,7}` (`lab/py/restricted-franel`, verb `sandwich`). On the full set the identity is the one Edwards proves in section 12.2.

## The sandwich and the equivalence

**Lemma 4.1 (the Jordan form).** Let `J_2(f) = f^2 prod_(p | f) (1 - p^(-2))` be Jordan's totient and, for `f <= Q`, `y_f = sum_(m <= Q/f) x_(fm)/m`. Then `G_F(Q) = sum_(f <= Q) (J_2(f)/f^2) y_f^2`, with every weight `J_2(f)/f^2` in `[6/pi^2, 1]`. Hence `(6/pi^2) sum_f y_f^2 <= G_F(Q) <= sum_f y_f^2` for every real vector `x`.

**Proof.** `J_2 = mu * id^2` as Dirichlet convolution, so `sum_(f | n) J_2(f) = n^2` by Mobius inversion, and `gcd(d, e)^2 = sum_(f | d, f | e) J_2(f)`. Then `G_F(Q) = sum_f J_2(f) (sum_(f | d) x_d/d)^2 = sum_f (J_2(f)/f^2) (sum_m x_(fm)/m)^2`, and `prod_(p | f) (1 - p^(-2))` lies between `prod_p (1 - p^(-2)) = 1/zeta(2)` and `1`. □

The vector `y` is itself a family of dilated sums: writing `n = mc`, `y_f = sum_(n <= Q/f, fn in S_F) w(n)` with `w = mu * (1/id)`, `|w(n)| = phi(rad n)/n <= 1`.

**Lemma 4.2 (Schur's test).** For a real matrix `A = (a_(ij))` with largest absolute row sum `R` and largest absolute column sum `C`, `||A x||^2 <= R C ||x||^2` in `l^2`.

**Proof.** By Cauchy and Schwarz, `(sum_j a_(ij) x_j)^2 <= (sum_j |a_(ij)|)(sum_j |a_(ij)| x_j^2) <= R sum_j |a_(ij)| x_j^2`; summing over `i` and exchanging, `||A x||^2 <= R sum_j x_j^2 sum_i |a_(ij)| <= R C ||x||^2`. □

**Theorem 4.3 (the sandwich).** Let `H_Q = sum_(k <= Q) 1/k` and `N_Q = H_Q max_(d <= Q) sigma(d)/d`, with `sigma` the sum of divisors. Then for every real vector `x` indexed by `d <= Q`, `sum_f y_f^2 <= N_Q sum_d x_d^2` and `sum_d x_d^2 <= N_Q sum_f y_f^2`, and `N_Q <= H_Q^2 <= (1 + ln Q)^2`. Consequently, for every digit design and every `Q >= 1`,

```
(6/pi^2) (1 + ln Q)^(-2) sum_(d <= Q) M_F(Q/d; d)^2  <=  G_F(Q)  <=  (1 + ln Q)^2 sum_(d <= Q) M_F(Q/d; d)^2 .
```

**Proof.** The map `x -> y` is the upper triangular matrix with entry `1/m` at `(f, fm)`, and it is inverted by the matrix with entry `mu(m)/m` at `(f, fm)`: `sum_m (mu(m)/m) sum_n x_(fmn)/n = sum_k (x_(fk)/k) sum_(m | k) mu(m) = x_f`. Row `f` of either matrix has absolute sum `sum_(m <= Q/f) 1/m <= H_Q`. Column `d` of the first has absolute sum `sum_(f | d) f/d = sigma(d)/d`; column `d` of the inverse has absolute sum `sum_(m | d) |mu(m)|/m = prod_(p | d) (1 + 1/p) <= sigma(d)/d`. And `sigma(d)/d = sum_(k | d) 1/k <= H_d <= H_Q`. Lemma 4.2 bounds both squared operator norms by `N_Q <= H_Q^2`, and `H_Q <= 1 + ln Q`. The display combines the two norm bounds with Lemma 4.1. □

Robin's unconditional bound, `sigma(d)/d < e^gamma ln ln d + 0.6483/ln ln d` for `d >= 3`, quoted here from the Franel note as Theorem 2 of Robin (1984) and not read at source, sharpens `N_Q` to `O(ln Q ln ln Q)`.

**Proposition 4.4 (no constant).** The ratio `sum_f y_f^2 / sum_d x_d^2` is unbounded over the real vectors `x` indexed by `d <= Q` as `Q` grows: for `x_d = 1` at the divisors of a primorial `N = prod_(p <= z) p <= Q` and `0` elsewhere, the ratio is `prod_(p <= z) (1 + 1/p + 1/(2 p^2)) > sum_(p <= z) 1/p`.

**Proof.** For `f | N`, `y_f = sum_(m | N/f) 1/m = sigma(N/f)/(N/f)`, and `y_f = 0` otherwise. So `sum_f y_f^2 = sum_(g | N) (sigma(g)/g)^2 = prod_(p <= z) (1 + (1 + 1/p)^2)` by multiplicativity over the squarefree `N`, while `sum_d x_d^2 = 2^(pi(z))`; the ratio is `prod_(p <= z) (1 + 1/p + 1/(2p^2))`, which exceeds `1 + sum_(p <= z) 1/p`. By the elementary bound `prod_(p <= z) p < 4^z`, the primorial of `z = ln Q / ln 4` is at most `Q`, and `sum_(p <= z) 1/p` diverges with `z`. □

By Mertens' theorem on the prime product (1874) the ratio is of order `ln z`, so with `z` of order `ln Q` it grows like `ln ln Q`; the witness rules out a constant and says nothing about whether `N_Q` is attained.

**Theorem 4.5 (the mean-square Franel equivalence).** For every digit design with at least two digits the following are equivalent, each read for every `eps > 0`: (i) `sum_j delta_j^2 = O_eps(Q^(-1+eps))` on the denominator-restricted set; (ii) `G_F(Q) = O_eps(Q^(alpha+eps))`; (iii) `sum_(f <= Q) y_f^2 = O_eps(Q^(alpha+eps))`; (iv) `sum_(d <= Q) M_F(Q/d; d)^2 = O_eps(Q^(alpha+eps))`.

**Proof.** For `Q >= 16 q^2` the node count is at least `1`, and by Corollary 3.5 `G_F(Q) = 12 m_F(Q) sum_j delta_j^2 + kappa` with `kappa` in `{1, -2}`, while Lemma 2.4 gives `Q^(1+alpha)/(C (ln Q)^2) < m_F(Q) < 2 |F| Q^(1+alpha)`. Under (i), `G_F(Q) <= 24 |F| Q^(1+alpha) C_eps Q^(-1+eps) + 1`, which is (ii); under (ii), `sum_j delta_j^2 = (G_F(Q) - kappa)/(12 m_F(Q)) <= (G_F(Q) + 2)/(12 m_F(Q)) << Q^(alpha+eps) (ln Q)^2 Q^(-1-alpha)`, which is (i) with `2 eps`. Lemma 4.1 gives (ii) if and only if (iii) with constants `6/pi^2` and `1`. Theorem 4.3 gives (iii) if and only if (iv) at the cost of `(1 + ln Q)^2`, absorbed by the `eps`. □

**Corollary 4.6 (the full set is Franel).** On the full set the four statements of Theorem 4.5 are each equivalent to `M(x) = O_eps(x^(1/2+eps))`, hence to the Riemann hypothesis.

**Proof.** Every dilate is the whole line, so (iv) reads `sum_(d <= Q) M(Q/d)^2 = O(Q^(1+eps))`. Its `d = 1` term is `M(Q)^2`, so (iv) gives `M(Q) = O(Q^(1/2+eps))`. Conversely `|M(x)| <= C_eps x^(1/2+eps)` gives `sum_(d <= Q) M(Q/d)^2 <= C_eps^2 Q^(1+2eps) sum_(d >= 1) d^(-1-2eps) = C_eps^2 zeta(1+2eps) Q^(1+2eps)`. That `M(x) = O_eps(x^(1/2+eps))` for every `eps` is equivalent to the Riemann hypothesis is Theorem 14.25 (C) of Titchmarsh (1986). □

Landau's criterion is recovered one way only: by Cauchy and Schwarz `sum_j |delta_j| <= (m sum_j delta_j^2)^(1/2)`, so (i) gives `sum_j |delta_j| = O(Q^((alpha+eps)/2))`, Landau's `O(Q^(1/2+eps))` on the full set; the converse is Landau's own argument and is not derived here.

**Corollary 4.7 (the ceiling).** On every digit design, (i) implies `M_F(Q) = O_eps(Q^(alpha/2+eps))`, since `M_F(Q)^2 = x_1^2 <= sum_d x_d^2`.

## The surrogate and one hypothesis

**Lemma 5.1 (the divisor bound).** For every `eps > 0` there is `C_eps` with `tau(m) <= C_eps m^eps` for all `m >= 1`, `tau` the number of divisors.

**Proof.** `tau(m) m^(-eps) = prod (k + 1) p^(-k eps)` over the prime powers `p^k` exactly dividing `m`. A factor with `p >= 2^(1/eps)` is at most `(k + 1) 2^(-k) <= 1`; a factor with `p < 2^(1/eps)` is at most `c_eps = max_(k >= 0) (k + 1) 2^(-k eps)`, and there are fewer than `2^(1/eps)` such primes. So `C_eps = c_eps^(2^(1/eps))` serves. □

**Theorem 5.2 (the surrogate is a theorem).** Let `B(Q) = sum_(d, e <= Q) (gcd(d, e)^2/(d e)) sqrt(N_F(Q; d) N_F(Q; e))`. For every digit design and every `Q >= 1`,

```
sum_(n in S_F, n <= Q) tau(n)  <=  B(Q)  <=  (1 + ln Q) sum_(d <= Q) (sigma(d)/d) N_F(Q; d)  <=  (1 + ln Q)^2 sum_(n in S_F, n <= Q) tau(n) ,
```

hence `A_F(Q) <= B(Q) <= C_eps (1 + ln Q)^2 Q^eps A_F(Q)` and `B(Q) = O_eps(Q^(alpha+eps))`, with no Mobius function anywhere.

**Proof.** Every term of `B(Q)` is nonnegative and the diagonal `d = e` has kernel `1`, so `B(Q) >= sum_(d <= Q) N_F(Q; d) = sum_(n in S_F, n <= Q) tau(n)`, exchanging the count of pairs `(d, n)` with `d | n`. For the upper bound, `sqrt(N_F(Q; d) N_F(Q; e)) <= (N_F(Q; d) + N_F(Q; e))/2` and the kernel is symmetric, so `B(Q) <= sum_(d <= Q) N_F(Q; d) sum_(e <= Q) gcd(d, e)^2/(d e)`. Group the inner sum by `f = gcd(d, e)`, a divisor of `d`, and write `e = f e'`: the kernel is `f/(d e')` and the sum over `e' <= Q/f` is at most `H_Q`, so the inner sum is at most `(1/d) sum_(f | d) f H_Q = (sigma(d)/d)(1 + ln Q)`. Exchanging again, `sum_(d <= Q) (sigma(d)/d) N_F(Q; d) = sum_(n in S_F, n <= Q) sum_(d | n) sigma(d)/d`, and since `sigma(d)/d = sum_(k | d) 1/k` grows along divisibility, `sum_(d | n) sigma(d)/d <= tau(n) sigma(n)/n <= tau(n)(1 + ln n)`. Lemma 5.1 and Lemma 2.2 finish. □

**Theorem 5.3 (one hypothesis).** Call (SR) the statement that for every `eps > 0` there is `C_eps` with `|M_F(Q/d; d)| <= C_eps N_F(Q; d)^(1/2 + eps)` for every `1 <= d <= Q` and every `Q >= 1`: square-root cancellation of `mu` in each dilate against that dilate's own mass. Under (SR), `sum_j delta_j^2 = O_eps(Q^(-1+eps))` on the denominator-restricted set of every digit design with at least two digits.

**Proof.** Square (SR) and sum: `N_F(Q; d)^(1+2eps) <= A_F(Q)^(2eps) N_F(Q; d)`, so `sum_(d <= Q) M_F(Q/d; d)^2 <= C_eps^2 A_F(Q)^(2eps) sum_(n in S_F, n <= Q) tau(n) = O(Q^(alpha + eps(1 + 2 alpha)))` by Theorem 5.2 and Lemma 2.2, and Theorem 4.5 turns the mean square into the threshold. Through the kernel the same line reads `G_F(Q) <= C_eps^2 A_F(Q)^(2eps) B(Q)`. □

So the converse rests on cancellation of `mu` alone, pointwise in each dilate against its own mass, and on no statement about how `S_F` sits in residue classes; at `d = 1` (SR) is the square-root conjecture for the design's Mertens sum, at `d = q^j` the same by Lemma 2.6, at every other `d` new. The threshold is not shown to return (SR): the one route tried, the termwise bound `|S_F(k, Q)| <= k (pi^2 G_F(Q)/6)^(1/2)` from Theorem 3.3 with the Mobius inversion `d M_F(Q/d; d) = sum_(c | d) mu(d/c) S_F(c, Q)` of Theorem 3.2, yields only `|M_F(Q/d; d)| <= (sigma(d)/d)(pi^2 G_F(Q)/6)^(1/2)`, which grows in `d` where (SR) asks for the size of the dilate; no counterexample to the implication is known, and none can be exhibited while both sides are open.

**Proposition 5.4 (the hypothesis with a rate).** Write `d_co` for the part of `d` coprime to `q` and call (U') the statement that for every `eps > 0` there is `C_eps` with `|M_F(x; d)| <= C_eps (1 + d_co^((alpha-1)/2) x^(alpha/2 + eps))` for every `d >= 1` and `x >= 1`. Then (U') implies (iv) of Theorem 4.5, hence the threshold. On a design carrying both `0` and `1` with `alpha < 1`, (U') fails at every `eps > 0` without its constant term, and so does the bound (U), `|M_F(x; d)| <= C_eps d^((alpha-1)/2) x^(alpha/2 + eps)`.

**Proof.** Write `d = a d_co` with every prime of `a` dividing `q`. The sum `M_F(Q/d; d)` is empty unless `N_F(Q; d) >= 1`, and under (U') it is otherwise at most `C_eps (1 + d_co^((alpha-1)/2) (Q/d)^(alpha/2+eps))` in size, so `sum_(d <= Q) M_F(Q/d; d)^2 <= 2 C_eps^2 (#{d <= Q : N_F(Q; d) >= 1} + Q^(alpha + 2eps) sum_d d_co^(alpha-1) d^(-alpha-2eps))`. The count is at most `sum_(d <= Q) N_F(Q; d) = sum_(n in S_F, n <= Q) tau(n)`, which is `O(Q^(alpha+eps))` by Lemmas 5.1 and 2.2. The sum is `sum_(a, d_co) a^(-alpha-2eps) d_co^(-1-2eps)`, at most `prod_(p | q) (1 - p^(-alpha-2eps))^(-1) zeta(1 + 2eps)`, finite because `alpha > 0`. Without the constant term, take `x = 1` and `d = q^k + 1`, which lies in `S_F` and is coprime to `q`: `M_F(1; d) = mu(1) = 1`, while `d^((alpha-1)/2)` tends to `0` as `k` grows. For (U), Lemma 2.6 gives `M_F(x; q^j) = M_F(x)`, so (U) at `d = q^j` would read `|M_F(x)| <= C_eps q^(j(alpha-1)/2) x^(alpha/2+eps)` for every `j`, and since `alpha < 1` the right side tends to `0` at fixed `x`, forcing `M_F(x) = 0` against `M_F(1) = mu(1) = 1`. □

(U') asks, at `x = Q/d`, that `|M_F(Q/d; d)|` be at most a constant plus `(Q^alpha/d_co)^(1/2)`, up to `Q^(O(eps))` and the base-smooth factor: it is (SR) with the mass `N_F(Q; d)` replaced by the main term `A_F(Q)/d` that equidistribution of `S_F` in the residue class `0 mod d` would give. The next section shows that the main term is not a uniform upper bound for the mass, so (SR) and (U') ask for different things there: at the repunits of Theorem 6.1 the mass exceeds `A_F(Q)/d` by a factor growing like `(3/2)^t`, and (U') asks there for more cancellation than (SR) does.

## What fails: the uniform law at the repunits

The accepting-set law is the statement `N_F(Q; d) ~ A_F(Q)/d` at fixed `d` coprime to the base, the design equidistributed in the residue class `0 mod d`; the Franel note derives it from the automaton's stationary law under `gcd(d, q Delta_F) = 1`, `Delta_F` the gcd of the differences of the digits, states it as a conjecture and reads it to three decimals on eleven coprime moduli of base 3 `{0,1}`. This section proves that it is not a uniform upper bound.

**Theorem 6.1 (the repunits).** At base 3 with `F = {0, 1}` let `R_t = (3^t - 1)/2`, the number written with `t` ones, `t >= 2`. Then `R_t` lies in `S_F`, is coprime to `3`, and

```
N_F(3^(2t); R_t) = 2^t + 1 ,    N_F(3^(3t); R_t) = 2 * 3^t + 1 ,    A_F(3^(2t)) = 4^t ,    A_F(3^(3t)) = 8^t .
```

Hence `N_F(Q; R_t) R_t / A_F(Q)` equals `(2^t + 1)(3^t - 1)/(2 * 4^t)`, of order `(3/2)^t`, at `Q = 3^(2t)` with `R_t < Q^(1/2)`, and `(2 * 3^t + 1)(3^t - 1)/(2 * 8^t)`, of order `(9/8)^t`, at `Q = 3^(3t)` with `R_t < Q^(1/3)`. No bound `N_F(Q; d) <= C A_F(Q)/d` holds uniformly over the `d <= Q^(1/2)` coprime to the base, nor over the `d <= Q^(1/3)`.

**Proof.** `R_t = 1 mod 3`, so `gcd(R_t, 3) = 1`, and its digits are ones. `3^t = 1 mod R_t` because `3^t - 1 = 2 R_t`. A member `n` of `S_F` with `1 <= n <= 3^(kt)` is either `3^(kt)` itself, which is `1 mod R_t` and not counted, or `n = n_0 + n_1 3^t + ... + n_(k-1) 3^((k-1)t)` with each `n_i` a padded string of `t` digits from `{0, 1}` read as an integer, and then `n = n_0 + ... + n_(k-1) mod R_t`. Write that sum as `sum_(i < t) c_i 3^i` with `c_i` the number of strings carrying a `1` at position `i`, `0 <= c_i <= k`; it lies in `[0, k R_t]`, so `R_t | n` exactly when it equals `j R_t = sum_(i < t) j 3^i` for some `0 <= j <= k`, and `j = 0` is `n = 0`, excluded. For `k = 2`: reducing mod `3` gives `c_0 = j` and then, dividing by `3` and inducting on `t`, every `c_i = j`. So `j = 1` is the `2^t` pairs with exactly one `1` in every position and `j = 2` the single pair of two all-ones strings: `N_F(3^(2t); R_t) = 2^t + 1`. For `k = 3` the same reduction gives `c_i = j` for `j = 1, 2`, which are `3^t` triples each; for `j = 3`, `c_0` is `0` or `3`, and `c_0 = 0` would put `sum_(i >= 1) c_i 3^(i-1)`, a number of `t - 1` digits at most `3` each and so at most `(3^t - 3)/2`, equal to `R_t`, impossible, so `c_0 = 3` and induction gives every `c_i = 3`, the single triple of all-ones strings: `N_F(3^(3t); R_t) = 2 * 3^t + 1`. The counts `A_F(3^(kt)) = 2^(kt)` are the `2^(kt) - 1` nonzero strings of at most `kt` digits together with `3^(kt)`. The ratios follow, and `R_t < 3^t = Q^(1/2)` at `k = 2`, `R_t < Q^(1/3)` at `k = 3`. □

The Franel note carries the same count at `Q = 3^(2rt)` for every `r` and proves there, with a Perron root certified exactly at `r <= 5`, that the failure reaches every level `Q^(1/(2r))` down to `Q^(1/10)`; that argument is not reproduced here. What Theorem 6.1 costs is the step from (SR) to (U'): the uniform law `N_F(Q; d) << A_F(Q)/d_co` would make (SR) read `|M_F(Q/d; d)| << d_co^(-1/2) Q^(alpha/2+eps)`, which is (U') up to the base-smooth factor.

**Fact 6.2.** The counts read `N_F(3^(2t); R_t) = 5, 9, 17, 33, 65, 129, 257` for `t = 2` to `8` and `N_F(3^(3t); R_t) = 19, 55, 163, 487` for `t = 2` to `5`; over every `d <= Q^(1/2)` at `Q = 3^8` to `3^16` the ratio `N_F(Q; d) d_co/A_F(Q)` is largest at the repunit of length `ceil(L/2)` at every level `L`, from `2.6562` at `d = 40` to `12.8625` at `d = 3280`, and at those cells the Mobius values of the dilate stay at most `0.8839` of `N_F^(1/2)` in absolute partial sum, the largest at `L = 11`, `d = 364`, where `N_F = 32` and the peak is `5`. Domain: base 3 `{0,1}`, levels `8` to `16`. Generator: `lab/py/restricted-franel`, verb `accepting`.

## The readings

All numbers in this section are printed by `uv run python research/franel/lab/py/restricted-franel/restricted_franel.py sandwich` and `... accepting`, run from the repository root.

**Fact 7.1 (the sandwich read).** At base 3 `{0,1}`, `Q = 3^4` to `3^8`, and on the full set at `Q = 40, 81, 243`, the Jordan form of Lemma 4.1 equals the gcd double sum as an exact rational at every `Q <= 729` and to `2.8e-14` and `5.7e-14` in floating point at `Q = 2187` and `6561`. With `X = sum_(d <= Q) x_d^2` and `Y = sum_(f <= Q) y_f^2`:

| set | `Q` | `G_F(Q)` | `X` | `Y` | `G_F/Y` | `G_F/X` | `x_1^2/X` | `N_Q` | `x -> y` | `y -> x` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| base 3 `{0,1}` | `81` | `12.574087` | `19` | `14.918` | `0.842880` | `0.661794` | `0.210526` | `13.94` | `2.258774` | `1.940642` |
| base 3 `{0,1}` | `243` | `32.163349` | `63` | `38.754` | `0.829928` | `0.510529` | `0.253968` | `18.82` | `2.490385` | `2.080130` |
| base 3 `{0,1}` | `729` | `65.900707` | `113` | `80.499` | `0.818649` | `0.583192` | `0.035398` | `24.08` | `2.693333` | `2.196254` |
| base 3 `{0,1}` | `2187` | `152.742675` | `261` | `181.387` | `0.842080` | `0.585221` | `0.187739` | `29.29` | `2.874809` | `2.299018` |
| base 3 `{0,1}` | `6561` | `399.401957` | `673` | `480.406` | `0.831385` | `0.593465` | `0.005944` | `35.95` | `3.040083` | `2.391003` |
| full set | `40` | `62.310829` | `57` | `72.949` | `0.854165` | `1.093172` | `0.000000` | `10.82` | `2.101228` | `1.847909` |
| full set | `81` | `159.157609` | `119` | `183.073` | `0.869367` | `1.337459` | `0.134454` | `13.94` | `2.258774` | `1.940642` |
| full set | `243` | `539.496506` | `341` | `634.112` | `0.850790` | `1.582101` | `0.011730` | `18.82` | `2.490385` | `2.080130` |

**Table 1.** `G_F/Y` sits inside `[6/pi^2, 1] = [0.607927, 1]` at all eight values of `Q`, as Lemma 4.1 says it must; `G_F/X` sits inside the corridor `[(6/pi^2)/N_Q, N_Q]` of Theorem 4.3, `[0.0436, 13.94]` at `Q = 81` and `[0.0169, 35.95]` at `Q = 6561`; no constant is claimed for either. The `d = 1` share `x_1^2/X` is under a quarter throughout, so the dilates at `d > 1` carry most of the mean square. The last two columns are the `l^2` norms of the two maps of Theorem 4.3 on the `Q x Q` matrices, Rayleigh quotients after `3000` power steps and so lower bounds, against `N_Q^(1/2) = 3.733351` to `5.995693` at `Q = 81` to `6561`: both grow and both sit under the bound. Generator: verb `sandwich`.

**Fact 7.2 (the surrogate read).** At base 3 `{0,1}`, `Q = 3^4` to `3^12`, `B(Q)` reads `200.233, 571.647, 1583.160, 4043.971, 10001.178, 23553.362, 54650.025, 125375.511, 285237.429`, the Jordan form agreeing with the gcd double sum to `2.4e-11` at `Q <= 3^8`; `B(Q)/sum tau(n)` climbs from `2.9019` to `4.4174`, the middle bound of Theorem 5.2 over `B(Q)` from `2.728` to `5.402`, and `B(Q)/(Q^alpha (ln Q)^2)` falls from `0.6480` to `0.4007`, no exponent claimed. Generator: verb `accepting`.

The figure is the same computation at `Q = 81` run by its own census writer, `paper-franel-converse`: the design's sixteen denominators from `mrly_math::design::elements`, the full set's nodes from the mediant walk `mrly_num::lattice::farey` and again from a gcd sweep, the Mobius values from `mrly_num::factor::mobius_sieve`. It asserts `m_F(81) = 241`, `m(81) = 2020`, `sum_j delta_j^2 = 0.004002105` and `0.006524654` to `1e-9`, `M_F(81) = -2`, `M(81) = -4`, `sum_d x_d^2 = 19`, `G_F(81) = 12.574087` and `G(81) = 159.157609` to `1e-5` against Table 1, and the rank form on both sets to `1e-8`.

## Open problems

Three things are open. First, (SR) itself: square-root cancellation of `mu` in every dilate against its own mass, which at `d = 1` is the square-root conjecture for the design's Mertens sum and at every `d` off the base-power ladder a statement about a regular language that is not a digit design; Theorem 5.3 says it is enough and Theorem 4.5 says only its mean square is needed, and nothing here proves either; it is a conjecture. Second, the true order of the operator norms of the two maps: Proposition 4.4 puts the ratio of the mean squares above `ln ln Q` on some vectors and Theorem 4.3 with Robin's bound puts every ratio below `O(ln Q ln ln Q)`, and the gap is untouched. Third, the accepting-set law at fixed coprime `d`, the Franel note's conjecture, which Section 6 refutes only as a uniform bound; on the strict set, both numerator and denominator in the design, there is no Franel-type statement to make until its limit measure is named, since at base 3 `{0,1}` the interval `[1/2, 2/3]` holds none of its fractions at any `Q`.

## Reproducibility

One study, `lab/py/restricted-franel`, prints every number of Sections 6 and 7 and the constants of Section 3: `uv run python research/franel/lab/py/restricted-franel/restricted_franel.py sandwich` from the repository root, one core, between `2.5` and `4.7` seconds on the two machines it was timed on, prints Table 1, the rank form's constant on six designs and the norms; `... accepting`, about `12` seconds, prints Fact 7.2 and the counts of Fact 6.2; neither takes a further argument, reads a file or writes one, and the study's README names the witness of every printed line. Its other verbs check the frequency-`m` formula against the literal sum of roots of unity and the Fourier form against the truncated Fourier side with a printed tail bound; their readings are on the Franel note. The figure is `figures/paper-franel-converse.ts`, drawn in the browser in under half a second a theme, with every drawn quantity asserted inside its census writer `research/shared/census/src/bin/paper-franel-converse.rs` as listed in Section 7.

## References

- Franel 1924, Les suites de Farey et le probleme des nombres premiers, Gott. Nachr., 198-201. [eudml.org/doc/59156](https://eudml.org/doc/59156)
- Landau 1924, Bemerkungen zu der obenstehenden Abhandlung von J. Franel, Gott. Nachr., 202-206. [eudml.org/doc/59157](https://eudml.org/doc/59157)
- Edwards 1974, Riemann's Zeta Function, Academic Press, section 12.2. [archive.org](https://archive.org/details/riemannszetafunc00edwa_0)
- Ramanujan 1918, On certain trigonometrical sums and their applications in the theory of numbers, Trans. Cambridge Philos. Soc. 22, no. 13, 259-276, equation (2.7). [ramanujan.sirinudi.org](https://ramanujan.sirinudi.org/Volumes/published/ram21.pdf)
- Smith 1875, On the value of a certain arithmetical determinant, Proc. London Math. Soc. 7, 208-212. [doi.org/10.1112/plms/s1-7.1.208](https://doi.org/10.1112/plms/s1-7.1.208)
- Schur 1911, Bemerkungen zur Theorie der beschrankten Bilinearformen mit unendlich vielen Veranderlichen, J. reine angew. Math. 140, 1-28. [doi.org/10.1515/crll.1911.140.1](https://doi.org/10.1515/crll.1911.140.1)
- Robin 1984, Grandes valeurs de la fonction somme des diviseurs et hypothese de Riemann, J. Math. Pures Appl. 63, 187-213, Theorem 2, not read at source. [zbmath.org/0516.10036](https://zbmath.org/0516.10036)
- Mertens 1874, Ein Beitrag zur analytischen Zahlentheorie, J. reine angew. Math. 78, 46-62. [doi.org/10.1515/crll.1874.78.46](https://doi.org/10.1515/crll.1874.78.46)
- Titchmarsh 1986, The Theory of the Riemann Zeta-Function, second edition revised by D. R. Heath-Brown, Clarendon Press, Theorem 14.25 (C). [sites.math.rutgers.edu](https://sites.math.rutgers.edu/~zeilberg/EM18/TitchmarshZeta.pdf)
- Huxley 1971, The distribution of Farey points I, Acta Arith. 18, 281-287, Lemma (10). [web.archive.org](https://web.archive.org/web/2020id_/http://matwbn.icm.edu.pl/ksiazki/aa/aa18/aa18130.pdf)
- Huxley 2012, Identities involving Farey fractions, Trudy Mat. Inst. Steklova 276, 131-145, Theorem 1. [mathnet.ru](https://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=tm&paperid=3374&option_lang=eng)
- Kanemitsu and Yoshimoto 1996, Farey series and the Riemann hypothesis, Acta Arith. 75, Theorem 3. [eudml.org/doc/206882](https://eudml.org/doc/206882)
- Cobeli and Zaharescu 2003, The Haros-Farey sequence at two hundred years, Acta Univ. Apulensis 5, 1-38. [emis.dsd.sztaki.hu](https://emis.dsd.sztaki.hu/journals/AUA/acta5/survey3.ps_pages1-20.pdf)
