---
title: The First Base Below a Quarter
lead: Base 21 missing the digit 0 is the least base at which a one-missing-digit set has the `l^1` exponent of its digit transform below the bar `1/4`, all 108 sets of the bases 3 to 20 certified above it; 34 is the least base at which every digit clears and every base from 34 clears at every digit, from 125 by a proof; under the generalized Riemann hypothesis each clearance is a power saving for the Mobius sum over the set.
date: 2026-09-21
figure: paper-first-base-below-a-quarter
shelf: first-base-below-a-quarter
---

Fix a base and cross out one digit; the integers that can still be written form a thin set with no multiplicative structure, and everything analytic known about it passes through the `l^1` exponent `alpha_1` of its digit transform, how much of the transform survives when the whole grid of frequencies is added up in absolute value. A quarter is a different bar: under the generalized Riemann hypothesis, expanding the Mobius sum over the set in additive characters and paying the uniform bound `x^(3/4+eps)` at each frequency beats the set's own counting function exactly when `alpha_1 < 1/4`, and Theorem 5.1 proves that transfer, `|M_F(x)| << A_F(x)^(1 - delta + eps)` with `delta = (1/4 - alpha_1)/alpha`. Base `21` with the digit `0` removed is the least base at which a one-missing-digit set falls below the bar, `alpha_1 in [0.2499765, 0.2499771]`, since every one of the `108` distinct sets of the bases `3` to `20` is certified above it. Base `34` is the least base at which every excluded digit clears, and every base from `34` on does: its `17` sets and the `3663` sets of the bases `35` to `125` are certified, every base from `21` to `33` keeps a digit certified above the bar, the last being `33` at the digit `15`, and an elementary majorant that forgets the digit proves `alpha_1 < 1/4` from base `125` on. Every certified number is interval arithmetic with directed rounding, its generator named.

## Introduction

![Every one-missing-digit set of the bases 10 to 40, one column per base and one dash per set at a float upper bound on its l^1 exponent from a window of three to five digits, 395 dashes in all; the quarter line crosses the picture, dashes above it dim, dashes below it yellow, the first yellow dash at base 21 and every dash of every column yellow from base 34 on.](paper-first-base-below-a-quarter)

Write every whole number in base `q` and throw away the ones that use the digit `a_0`: base `10` missing `7` leaves `9^L` of them below `10^L`, a set of dimension `0.954`. Maynard (2019) proved that infinitely many of them are prime in base `10`, reading the set through its transform `hat F(t) = (1/(q-1)) sum_(a != a_0) e(at)`, `e(x) = exp(2 pi i x)`, whose products over `t, qt, ..., q^(L-1) t` are the transforms of the `L`-digit strings. A sieve reads how much of that product survives when the `q^L` frequencies of the grid are added up in absolute value, and the growth rate of that sum, the `l^1` exponent `alpha_1` of Definition 2.5, is the number this paper is about; the published values sit near a third.

The figure is the census: one column per base from `10` to `40`, one dash per distinct set at an upper bound on its exponent, the line at `1/4`. The columns sink left to right, the spread inside a column is the excluded digit, cheapest at the ends and dearest near the middle, one dash crosses the line at base `21`, the digit `0`, and from base `34` on whole columns are under it. The bar is worth crossing because of what it buys: on [the Mobius page](../notes/grh-power-saving.md) the Mobius sum `M_F(x) = sum_(n in S_F, n <= x) mu(n)` is expanded in additive characters modulo `q^L`, at the cost of the `l^1` mass of the transform, each frequency carries an exponential sum which under the generalized Riemann hypothesis is `O(x^(3/4+eps))` uniformly by Baker and Harman (1991), and the two costs beat `A_F(x)` exactly when `alpha_1 < 1/4`, Theorem 5.1: every dash under the line is, under the hypothesis, a power saving for a Mobius sum.

Three results say where the dashes cross. The first is a theorem: `alpha_1 < 1/4` for every base `q >= 125` and every excluded digit, with the explicit bound of Theorem 4.7, uniform in the digit because its first step throws the digit away: the triangle inequality leaves a majorant built from the Dirichlet kernel, its level product expands over subsets of the digit positions, every maximal run telescopes into one Dirichlet kernel at a higher modulus, and what remains is a sum of Lebesgue sums, each bounded by Lemma 4.4. The other two are computations by a transfer matrix on windows of digits and the test of Collatz and Wielandt: base `21` missing `0` below the bar and every base below `21` above it at every digit (Facts 6.1 and 6.2), base `34` below at every digit and every base from `21` to `33` above it at some digit (Fact 6.3), every set of the bases `35` to `125` below (Fact 6.4).

The literature, as read here. The base-`10` value behind Maynard (2019) is `alpha_1 <= 27/77 = 0.35065`, the eigenvalue bound of its display (10.5), a finite-window upper bound by Fact 6.6; Karwatowski (2022), not read at source and carried from Granville (2024) and the author's base-`9` paper, carries the same shape of bound to every base `q >= 10`. Erdos, Mauduit and Sarkozy (1998) and Konyagin (2001) distribute such sets in residue classes, Maynard (2022), Nath (2024) and Leng and Sawhney (2025) consume the exponent for primes, and Chow, Varju and Yu (2024) read the complementary exponent of the unshifted sum against `1/2`, their Proposition 2.4 giving `alpha_1^0 < 1/2` for every base `q >= 5` with one digit removed. None fixes the least base at which a bar is met, and none carries a Mobius sum over a missing-digit set, the nearest being the divisor function of Kim (2024).

## Definitions

**Definition 2.1 (design).** Fix a base `q >= 3` and a digit set `F`, a subset of `{0, ..., q-1}` with `k = |F| >= 2`; `q` is the base and `k` the fill of this tree's vocabulary, and the two letters stand for them throughout. The design `S_F` is the set of positive integers whose every base-`q` digit lies in `F`, `alpha = log k/log q` its dimension, `A_F(x) = #{n in S_F : n <= x}`. For `L >= 1`, `D_L` is the set of the `k^L` integers `0 <= n < q^L` whose `L` padded digits lie in `F`: when `0 in F` it is `S_F` below `q^L` with `0` added, and when `0` is not in `F` no string over `F` starts with a zero and `S_F` below `q^L` is the disjoint union of `D_1, ..., D_L`. The one-missing-digit set at the digit `a_0` is `F = {0, ..., q-1}` less `{a_0}`, `k = q - 1`.

**Lemma 2.2 (the mass).** For `x >= q`, `x^alpha/k^2 < A_F(x) < 2k x^alpha`.

**Proof.** Let `L = floor(log_q x) >= 1`. At most `k^j` members of `S_F` have `j` digits, so `A_F(x) <= sum_(j <= L+1) k^j < (k^2/(k-1)) k^L <= 2k q^(L alpha) <= 2k x^alpha`. Some digit of `F` is nonzero; the `k^(L-1)` strings of exactly `L` digits with that leading digit lie in `S_F` below `q^L <= x`, and `k^(L-1) = q^((L+1) alpha)/k^2 > x^alpha/k^2`. □

**Definition 2.3 (the transform).** `hat F(t) = (1/k) sum_(a in F) e(at)` and `hat F_L(t) = prod_(j < L) hat F(q^j t) = k^(-L) sum_(n in D_L) e(nt)`, the normalised transform of `D_L`; both are `1`-periodic, `|hat F| <= 1`, and `|hat F(1 - t)| = |hat F(t)|`.

**Lemma 2.4 (the closed form and the mirror).** Let `D_q(t) = sin(pi q t)/sin(pi t)`, `q` at integers. For the one-missing-digit set at `a_0`, `(q-1) |hat F(t)| = |D_q(t) - e((a_0 - (q-1)/2) t)|`, and `|hat F|` is unchanged by `a_0 -> q - 1 - a_0`, so a base carries `ceil(q/2)` distinct sets, the pair coinciding only at odd `q` and `a_0 = (q-1)/2`.

**Proof.** The geometric sum `sum_(a < q) e(at) = e((q-1)t/2) D_q(t)` less the excluded term, times the unimodular `e(-(q-1)t/2)`; and `a -> q - 1 - a` carries the set missing `a_0` to the set missing `q - 1 - a_0` and the sum to `e((q-1)t)` times its conjugate. □

**Definition 2.5 (the grid sums and the exponent).** For `N >= 1` and real `x`,

```
Sigma_N(x) = sum_(i < q^N) |hat F_N(x + i/q^N)| ,   S_N = max_x Sigma_N(x) ,   alpha_1 = lim_N (1/N) log_q S_N = inf_N (1/N) log_q S_N .
```

The limit exists by Lemma 2.6; the exponent `alpha_1^0` of `q^L int_0^1 |hat F_L|`, the average of `Sigma_L` over one cell, is at most `alpha_1`.

**Lemma 2.6 (subadditivity).** `Sigma_N` is `q^(-N)`-periodic, `S_(N+M) <= S_N S_M`, and `S_N` is nondecreasing in `N`; hence `(1/N) log_q S_N` converges to its infimum.

**Proof.** Replacing `x` by `x + q^(-N)` permutes the summands. Write `i = b + q^N a` with `b < q^N`, `a < q^M`, and `t = x + i q^(-(N+M))`. Since `hat F` is `1`-periodic, `q^M t = q^M x + b q^(-N) (mod 1)`, so the last `N` factors of `|hat F_(N+M)(t)|` are the summand of `Sigma_N(q^M x)` indexed by `b` and the first `M` factors the summand of `Sigma_M(x + b q^(-N-M))` indexed by `a`; summing over `a` first, `Sigma_(N+M)(x) = sum_(b < q^N) |hat F_N(q^M x + b q^(-N))| Sigma_M(x + b q^(-N-M))`, at most `S_M Sigma_N(q^M x) <= S_N S_M` and at least `Sigma_N(q^M x) min_y Sigma_M(y)`. At any shift `y`, Parseval on `Z/q^M` gives `sum_i |hat F_M(y + i/q^M)|^2 = q^M k^(-M)` with every term at most `1`, so `Sigma_M(y) >= q^M k^(-M) >= 1` and `S_(N+M) >= S_N`. Fekete's lemma gives the limit. □

**Lemma 2.7 (the `l^1` floor).** For every base and every digit set, `Sigma_N(0) >= q^(N(1 - alpha))`, hence `alpha_1 >= 1 - alpha`, which is `0.01602` at base `21` missing one digit.

**Proof.** The Parseval identity of Lemma 2.6 at `y = 0`, with `sum_i |z_i| >= (sum_i |z_i|^2)/max_i |z_i|` and `max_i |z_i| = |hat F_N(0)| = 1`. □

## The window machine

**Proposition 3.1 (the window machine).** Fix `n >= 2` and split the circle into the `q^n` cells `I_w = [w q^(-n), (w+1) q^(-n))`. Put `G^+(w) = sup_(I_w) |hat F|`, `G^-(w) = inf_(I_w) |hat F|`, and let `M^(+-)` be the nonnegative `q^(n-1) x q^(n-1)` matrices `(M^(+-) y)(v) = sum_(c < q) G^(+-)(vq + c) y((vq + c) mod q^(n-1))`. Then `alpha_1 <= log_q rho(M^+)`, and `alpha_1 >= log_q rho(M'')` for the restriction `M''` of `M^-` to any subset of states; and for any `y > 0`, `M^+ y <= mu y` componentwise gives `rho(M^+) <= mu`, while `M'' y >= mu y` gives `rho(M'') >= mu`.

**Proof.** Write the digits of a cell index `w < q^l` as `w = (w_1 ... w_l)`, most significant first, so that `I_w` is the set of `t` whose expansion starts `0.w_1 ... w_l`; for `t in I_w` and `j + n <= l`, `q^j t` reduced mod `1` lies in the level-`n` cell `v_j(w) = (w_(j+1) ... w_(j+n))`.

The upper bound. Take `l = L`. For `t in I_w` and `j <= L - n` the factor `|hat F(q^j t)|` is at most `G^+(v_j(w))`, the remaining `n - 1` factors at most `1`, so `U_L = sum_(w < q^L) sup_(I_w) |hat F_L| <= sum_(w < q^L) prod_(j <= L-n) G^+(v_j(w)) = 1^T (M^+)^(L-n+1) 1`, the identity being the definition of `M^+`: a state is a string of `n - 1` digits, appending `c` forms the window `vq + c` with weight `G^+(vq + c)` and moves to `(vq + c) mod q^(n-1)`. Since `x + i q^(-L)` lies in `I_i` for `x in [0, q^(-L))`, `S_L <= U_L`, and Gelfand's formula gives `limsup U_L^(1/L) <= rho(M^+)`.

The lower bound. Take `l = L + n - 1`, so every `j < L` has `j + n <= l` and on `I_w` every factor of `|hat F_L|` is at least `G^-(v_j(w))`. Hence `int_0^1 |hat F_L| >= sum_(w < q^l) q^(-l) inf_(I_w) |hat F_L| >= q^(-(L+n-1)) 1^T (M^-)^L 1`, and `q^L int_0^1 |hat F_L| <= S_L` gives `S_L >= q^(-(n-1)) 1^T (M^-)^L 1`. The entries are nonnegative, so `1^T (M^-)^L 1 >= 1^T (M'')^L 1`, and if `y > 0` on `S'` with `M'' y >= mu y` then `1^T (M'')^L 1 >= mu^L (sum_(S') y)/(max_(S') y)`, whence `alpha_1 >= log_q mu`; the best `mu` is `rho(M'')`, by the Perron vector of an irreducible matrix or a further restriction of a reducible one.

The tests. If `y > 0` and `M^+ y <= mu y` then every entry of `(M^+)^L` is `O(mu^L)`, so `rho(M^+) <= mu`; the lower test is the display above. □

The last sentence is the test of Collatz (1942) and Wielandt (1950): a float power iteration supplies `y`, one pass in interval arithmetic turns it into a rigorous bound, and a poor vector weakens the bound and never breaks it. A cell containing a zero of `|hat F|` has `G^- = 0`, and near a zero a whole row of `M^-` can vanish, so the lower test runs on the restriction `M''` to the states where the float Perron vector exceeds `10^(-9)` of its maximum, which is the plain test wherever every state is live. Cell suprema and infima are enclosed by a sub-scan of `m` points plus the Lipschitz correction `Lambda/(2 m q^n)`, `Lambda = (2 pi/k) sum_(a in F) a`.

## The digit-uniform bound

Throughout this section `F` is the one-missing-digit set at the digit `a_0`, and `c_1 = 2/pi`.

**Lemma 4.1 (the digit-blind majorant).** For every `q >= 3`, every excluded digit and every real `t`, `|hat F(t)| <= u_q(t) := min(1, (|D_q(t)| + 1)/(q-1))`, and the right side does not mention `a_0`.

**Proof.** The triangle inequality on Lemma 2.4, and `|hat F| <= 1`. □

The discarded phase makes the majorant strictly larger than `|hat F|` almost everywhere, so the uniform bound is weaker than the per-digit one by construction; Fact 6.7 measures the cost. Write `v_q(t) = (|D_q(t)| + 1)/(q-1) >= u_q(t)` and `Sigma^v_N(x) = sum_(i < q^N) prod_(j < N) v_q(q^j (x + i/q^N)) >= Sigma_N(x)`. A subset `E` of `{0, ..., N-1}` splits uniquely into maximal runs of consecutive positions, `(s, l)` the run occupying `s, ..., s+l-1`.

**Lemma 4.2 (the run identity).** For every `x`, `(q-1)^N Sigma^v_N(x) = sum_(E) sum_(i < q^N) prod_((s,l) run of E) |D_(q^l)(q^s (x + i/q^N))|`.

**Proof.** Expand `prod_(j < N) (|D_q(q^j t)| + 1)` over the subsets `E` of positions at which the kernel is taken; inside a run, with `t' = q^s t`, `prod_(r < l) sin(pi q^(r+1) t')/sin(pi q^r t') = D_(q^l)(t')`, and distinct runs occupy disjoint blocks. □

Every term is a Lebesgue sum. Let `L_M = max_x sum_(j < M) |D_M(x + j/M)|` and `lambda_l = L_(q^l)/q^l`.

**Lemma 4.3 (the peel).** Let `(s_1, l_1), ..., (s_r, l_r)` be the maximal runs of `E`, `s_1 < ... < s_r`, and `t_i = x + i q^(-N)`. Then `sum_(i < q^N) prod_k |D_(q^(l_k))(q^(s_k) t_i)| <= q^N prod_k lambda_(l_k)` for every real `x`.

**Proof.** First, for `M >= l >= 1` and real `y`, `sum_(c < q^M) |D_(q^l)(y + c/q^M)| <= q^M lambda_l`: writing `c = rho + q^(M-l) c'` with `rho < q^(M-l)`, `c' < q^l`, the inner sum over `c'` runs over `q^l` points spaced `q^(-l)` apart and is at most `L_(q^l)`. Set `a_s = i mod q^(N-s)`; by periodicity the `k`-th factor depends on `i` through `a_(s_k)` alone, `a_(s_1)` determines every `a_(s_k)`, and the top `s_1` digits of `i` enter no factor, so it suffices to prove `sum_(a_(s_1) < q^(N-s_1)) prod_k |D_(q^(l_k))(q^(s_k) x + a_(s_k)/q^(N-s_k))| <= q^(N-s_1) prod_k lambda_(l_k)`. Induct on `r`: for `r = 1` it is the first display with `M = N - s_1 >= l_1`; for `r >= 2` write `a_(s_1) = a_(s_2) + q^(N-s_2) w` with `w < q^(s_2-s_1)`, so the factors `k >= 2` depend on `a_(s_2)` alone and the first is `|D_(q^(l_1))(z + w/q^(s_2-s_1))|` at a fixed `z`; runs are maximal, so `s_2 - s_1 > l_1`, the first display bounds the inner sum by `q^(s_2-s_1) lambda_(l_1)`, and the induction hypothesis bounds the outer sum. □

**Lemma 4.4 (the Lebesgue constant).** Let `M >= 2` and `gamma' = (2/pi)(gamma + log(8/pi)) = 0.9625228...`, `gamma` Euler's constant. Then (i) `L_M` is attained at the half offset, `L_M = 2 sum_(m < ceil(M/2)) csc((2m+1) pi/(2M)) - (M mod 2)`; (ii) `L_M <= M((2/pi) log M + gamma') + 2/pi`; (iii) `lambda_l <= c_1 l log q + gamma' + c_1 q^(-l)` for every `l >= 1` and `q >= 2`.

**Proof.** (i) Write `x = theta/M`. All `M` points `x + j/M` share `|sin(pi M t)| = |sin(pi theta)| =: s`, so with `d_j` the distance from `x + j/M` to `Z` the sum is `s sum_j csc(pi d_j)`. Replacing `theta` by `1 - theta` reflects the point set, so take `theta in [0, 1/2]`; `theta = 0` gives `M`, below the half-offset value, since the pair at `m = 0` alone contributes `2 csc(pi/(2M)) > 4M/pi > M`. Let `theta in (0, 1/2]`. The distances are the pairs `{(m + theta)/M, (m + 1 - theta)/M}` for `m < ceil(M/2) - [M odd]`, with the single distance `((M-1)/2 + theta)/M` when `M` is odd. Put `mu_m = pi (m + 1/2)/M` and `delta = pi (1/2 - theta)/M in [0, pi/(2M))`, so the pair sits at `mu_m -+ delta`, `s = cos(M delta)`, and `csc(mu - delta) + csc(mu + delta) = 2 sin mu cos delta/(sin^2 mu - sin^2 delta)` writes its contribution as `2 csc(mu_m) X_m` with `X_m = cos(M delta) cos delta sin^2 mu_m/(sin^2 mu_m - sin^2 delta)`. The claim `X_m <= 1` reads `sin^2 mu_m (1 - cos(M delta) cos delta) >= sin^2 delta`, whose left side is nondecreasing in `mu_m` on `(0, pi/2]`, so it is enough at `m = 0`: put `u = mu_0 = pi/(2M) <= pi/4` and `delta = ru` with `0 <= r < 1`, so `M delta = pi r/2`. Since `1 - cos(pi r/2) cos(ru) = 2 sin^2(ru/2) + 2 cos(ru) sin^2(pi r/4)` and `sin^2(ru) = 4 sin^2(ru/2) cos^2(ru/2)`, the claim is

```
sin^2 u cos(ru) sin^2(pi r/4) >= sin^2(ru/2) (cos^2 u + cos(ru)) .
```

The function `x cot x` decreases on `(0, pi)`, its derivative being `(sin 2x - 2x)/(2 sin^2 x)`, so `f(r) = sin(pi r/4)/sin(ru/2)`, whose derivative has the sign of `(pi/4) cot(pi r/4) - (u/2) cot(ru/2) = (1/r)[(pi r/4) cot(pi r/4) - (ru/2) cot(ru/2)]` with `ru/2 < pi r/4` since `u < pi/2`, decreases on `(0, 1]` and `f(r) >= f(1) = 1/(sqrt 2 sin(u/2))`; at `r = 0` the claim is `0 >= 0`. Hence `sin^2 u cos(ru) sin^2(pi r/4) >= sin^2(ru/2) sin^2 u cos(ru)/(2 sin^2(u/2)) = sin^2(ru/2) (1 + cos u) cos(ru)`, and `(1 + cos u) cos(ru) - cos^2 u - cos(ru) = cos u (cos(ru) - cos u) >= 0`. So every pair is at most `2 csc(mu_m)`, with equality at `delta = 0`, and the odd leftover at `pi/2 - delta` contributes `cos(M delta)/cos delta <= 1` since `delta <= M delta < pi/2`. At the half offset `s = 1`, the pairs are `2 csc((2m+1) pi/(2M))` and the leftover is `csc(pi/2) = 1`, which the formula counts once.

(ii) Split `csc x = 1/x + g(x)`. The expansion `g(x) = sum_(k >= 1) 2 (2^(2k-1) - 1) |B_(2k)| x^(2k-1)/(2k)!` has nonnegative coefficients, so `g` is nonnegative, increasing and convex on `(0, pi)`, `int_0^(pi/2) g = [log tan(x/2) - log x]_0^(pi/2) = log(4/pi)`, and the midpoint rule underestimates a convex integral; two harmonic inequalities come from the convexity of `1/x`: the trapezoid rule on `[N, N+1]` gives `log(1 + 1/N) <= (1/2)(1/N + 1/(N+1))`, so `H_N - log N - 1/(2N)` is nondecreasing, and the midpoint rule on `[N + 1/2, N + 3/2]` gives `log((N + 3/2)/(N + 1/2)) >= 1/(N+1)`, so `H_N - log(N + 1/2)` is nonincreasing; both tend to `gamma`, whence `H_N <= log N + gamma + 1/(2N)` and `H_N >= log(N + 1/2) + gamma` for every `N >= 1`.

Let `M = 2n` be even. The points `x_m = (2m+1) pi/(2M)`, `m < n`, are the midpoints of the `n` intervals of length `pi/M` partitioning `(0, pi/2]`, so `sum_(m < n) g(x_m) <= (M/pi) log(4/pi)`, while `sum_(m < n) 1/x_m = (2M/pi)(H_(2n) - H_n/2) <= (2M/pi)((1/2)(log M + log 2 + gamma) + 1/(2M))`; doubling, `L_M <= (2M/pi)(log M + gamma + log(8/pi)) + 2/pi`.

Let `M = 2n + 1` be odd, `n >= 1`. By (i), `L_M = 2 sum_(m < n) csc(x_m) + 1`, the `x_m` now the midpoints of `n` intervals of length `pi/M` inside `(0, pi/2)`, so `sum_(m < n) g(x_m) <= (M/pi) log(4/pi)` as before. With `2n = M - 1` and `n + 1/2 = M/2` the harmonic inequalities give `H_(2n) - H_n/2 <= log(M-1) + gamma + 1/(2(M-1)) - (1/2)(log(M/2) + gamma) <= (1/2)(log M + log 2 + gamma) - 1/M + 1/(2(M-1))`, using `log(M-1) <= log M - 1/M`. Multiplying by `4M/pi`, the two small terms give `-4/pi + 2M/(pi(M-1)) = -2/pi + 2/(pi(M-1))`, so with the `g` part and the leftover `1`, `L_M <= (2M/pi)(log M + gamma + log(8/pi)) + 1 - 2/pi + 2/(pi(M-1))`. For `M >= 5` the constant is at most `1 - 3/(2 pi) < 0.53 < 2/pi`; at `M = 3`, `L_3 = 2 csc(pi/6) + 1 = 5` against `3((2/pi) log 3 + gamma') + 2/pi > 5.62`.

(iii) is (ii) at `M = q^l` divided by `q^l`. □

**Lemma 4.5 (the count).** Put `A_0 = 1` and `A_N = sum_(E) prod_(runs) lambda_l`. Then `A_N = A_(N-1) + sum_(l < N) lambda_l A_(N-1-l) + lambda_N` for `N >= 1`, and if `z > 1` is the unique root of `z = 1 + sum_(l >= 1) lambda_l z^(-l)` then `limsup A_N^(1/N) <= z`; moreover `z` is nondecreasing in the sequence `(lambda_l)`.

**Proof.** Split by the first position: if `0` is not in `E` the rest contributes `A_(N-1)`; if the initial run has length `l < N`, position `l` is absent and the rest contributes `lambda_l A_(N-1-l)`; `l = N` contributes `lambda_N`. Let `Lambda(x) = sum_l lambda_l x^l`, of radius `1` by Lemma 4.4 and infinite at `1`; the recursion gives `A(x)(1 - x - x Lambda(x)) = 1 + Lambda(x)` for the generating function. `phi(z) = 1 + Lambda(1/z)` decreases strictly from infinity to `1` on `(1, infinity)`, so `z - phi(z)` has a unique zero, and `phi` increases pointwise with the `lambda_l`. For `|x| < 1/z`, `|x + x Lambda(x)| <= |x|(1 + Lambda(|x|)) < z^(-1) phi(z) = 1`, so `A` is analytic there, and its coefficients being nonnegative its radius is at least `1/z`. □

**Corollary 4.6.** `alpha_1 <= log_q(z q/(q-1))` with `z` as in Lemma 4.5, since Lemmas 4.1 to 4.3 give `(q-1)^N Sigma_N(x) <= q^N A_N` for every `x`.

**Theorem 4.7 (the digit-uniform bound).** For every base `q >= 3` and every excluded digit `a_0`,

```
alpha_1(q, a_0) <= log_q (z_q q/(q-1)) ,   where z_q is the unique root z > 1 of   (z-1)^3 = c_1 (log q) z + gamma' (z-1) + c_1 (z-1)^2/(qz - 1) .
```

In particular `alpha_1 < 1/4` for every `q >= 125` and every excluded digit; the right side reads `0.249980` at `q = 125`, `0.224161` at `250`, `0.182916` at `1250` and `0.108912` at `10^6`.

**Proof.** By Lemma 4.4(iii), `lambda_l <= Lambda_l := c_1 l log q + gamma' + c_1 q^(-l)`, and by the monotonicity of Lemma 4.5 it is enough to bound the root of `z = 1 + sum_l Lambda_l z^(-l)`; summing the three series, `sum_l Lambda_l z^(-l) = c_1 (log q) z/(z-1)^2 + gamma'/(z-1) + c_1/(qz - 1)`, and multiplying by `(z-1)^2` gives the cubic, so Corollary 4.6 gives the first claim. For the second, `alpha_1 < 1/4` follows from `z_q < w_q := q^(1/4)(1 - 1/q)`. Since `phi(z) = 1 + sum_l Lambda_l z^(-l)` is strictly decreasing and meets the diagonal once, `z_q < w` if and only if `phi(w) < w`, which after multiplying by `(w-1)^2` is

```
(w-1)^3 > c_1 (log q) w + gamma' (w-1) + c_1 (w-1)^2/(qw - 1) .   (4.1)
```

Fact 4.8 verifies (4.1) at `w = w_q` for every `125 <= q < 3000`. For `q >= 3000`: at the root, `(z_q - 1)^3 >= c_1 (log q) z_q > 1` once `q >= 5`, so `z_q >= 2` and `z_q <= 2(z_q - 1)`; with `(z_q - 1)/(q z_q - 1) < 1/q` the cubic gives `(z_q - 1)^2 <= 2 c_1 log q + gamma' + c_1/q <= 2 c_1 log q + 0.97` for `q >= 86`, so `z_q <= g(q) := 1 + sqrt(2 c_1 log q + 0.97)`. At `q = 3000`, `log q < 8.01` gives `2 c_1 log q + 0.97 < 11.18 < 3.35^2` and `g(3000) < 4.35`, while `7.4^4 < 3000` gives `w_3000 > 7.39`; and `w_q - g(q)` increases from there, since `d w_q/dq > (1/4) q^(-3/4)` and `d g/dq = c_1/(q sqrt(2 c_1 log q + 0.97))`, so it suffices that `(1/4) q^(1/4) sqrt(2 c_1 log q + 0.97) > c_1`, which reads `6.1 > 0.64` at `q = 3000` with an increasing left side. □

**Fact 4.8 (the certificate).** With `gamma'` rounded up to `0.9625229`, (4.1) holds at `w = w_q` for every integer `125 <= q < 3000` in interval arithmetic at `120` bits, tightest margin `3.985 * 10^(-3)` at `q = 125`, and fails at `q = 124` with margin below `-6.16 * 10^(-2)`; a float root search finds no failure on `[125, 10^5)`, the four values in Theorem 4.7 are rounded up in the sixth place, and the input rounded to `c_1 l log q + 0.97` closes at `126`. The scan of `L_M` at ten moduli from `4` to `10^5`, three of them odd, finds the maximum at the half offset every time and the excess `L_M/M - (2/pi) log M` falling from `0.965217` to `0.962523`, under `gamma' + (2/pi)/M` throughout. Generator: `coprime/lab/py/digit-uniform-bound`, verbs `threshold` and `lebesgue`.

## The Mertens bar

Write `M_F(x) = sum_(n in S_F, n <= x) mu(n)`, with `mu` the [Mobius function](https://en.wikipedia.org/wiki/M%C3%B6bius_function), and `S_mu(x, theta) = sum_(n <= x) mu(n) e(n theta)`.

**Theorem 5.1 (the Mertens bar).** Assume the generalized Riemann hypothesis, that `L(s, chi)` has no zero in `Re s > 1/2` for every Dirichlet character `chi`. Let `F` be any digit set at base `q` with `k >= 2`, dimension `alpha` and `l^1` exponent `alpha_1`. Then for every `eps > 0` and all `x >= 2`,

```
|M_F(x)| <<_(q,F,eps) x^(alpha + alpha_1 - 1/4 + eps) ,   equivalently   |M_F(x)| <<_(q,F,eps) A_F(x)^(1 - delta + eps)   with   delta = (1/4 - alpha_1)/alpha .
```

The constant depends on `q`, `F` and `eps` alone, and the saving is a power exactly when `alpha_1 < 1/4`.

**Proof.** Step 1, orthogonality. For `0 <= n < q^L`, completeness of the characters modulo `q^L` gives `1_(D_L)(n) = q^(-L) sum_(a < q^L) k^L hat F_L(a/q^L) e(-na/q^L)`, so for `y < q^L`, `|sum_(n in D_L, n <= y) mu(n)| <= q^(-L) k^L Sigma_L(0) max_theta |S_mu(y, theta)|`. Step 2, the `l^1` mass. For every `e > alpha_1` there is `N` with `S_N <= q^(N e)`; writing `L = aN + b` with `b < N`, Lemma 2.6 gives `S_L <= S_N^a S_b <= S_N^(a+1) <= q^(e(L+N))`, so `Sigma_L(0) <<_e q^(L e)`. Step 3, the Mobius input, quoted. Under the hypothesis, `max_theta |S_mu(y, theta)| <<_eps y^(3/4 + eps)`, the case `a = 1/2` of Baker and Harman (1991), whose hypothesis is that `L(s, chi)` is zero-free in `Re s > a` for every character and whose maximum is over all real `theta`, restated at source in Zhang (2024) and Porritt (2018); the one step not derived here. Step 4, assembly. Let `L = ceil(log_q(x+1))`, so `x < q^L <= 2qx`. If `0 in F`, `S_F` below `x` is `D_L` less `{0}` in `[1, x]`, and steps 1 to 3 give `|M_F(x)| << q^(L(alpha - 1 + e)) x^(3/4 + eps) << x^(alpha + e - 1/4 + eps)`, using `q^L <= 2qx` when `alpha - 1 + e > 0` and `q^L > x` otherwise. If `0` is not in `F`, `S_F` below `x` is the disjoint union over `l <= L` of `D_l` in `[1, x]`, each block bounded the same way with `y = min(x, q^l - 1)`; by Lemma 2.7, `alpha + e > 1`, so the sum over `l` is geometric and dominated by its top term. Take `e = alpha_1 + eps`. Step 5, the yardstick. By Lemma 2.2, `x^alpha` and `A_F(x)` agree up to the constants `k^2` and `2k` for `x >= q`, so `x^(alpha(1 - delta) + eps) <<_(q,F,eps) A_F(x)^(1 - delta + eps/alpha)`. □

The Mobius page proves the same theorem with a digit-free kernel bound in place of `alpha_1`, closing at one missing digit only from base `1499`, or `1032` at an end digit; the certified exponent moves that wall to `34`. Under a zero-free half plane `Re s > a` with `a > 1/2` the bar falls below `1/4`, by the tables of Baker and Harman (1991) and Zhang (2024), and no set here is certified under it.

**Corollary 5.2 (where the bar is cleared).** Under the generalized Riemann hypothesis, `|M_F(x)| <<_(q,F,eps) A_F(x)^(1 - delta + eps)` with `delta > 0` holds (i) at base `21` missing `0`, with `delta > 2.3 * 10^(-5)`; (ii) at every base `q >= 34` and every excluded digit, with `delta > 6.3 * 10^(-4)` at `q = 34`; (iii) at every base `q >= 125` and every excluded digit with `delta >= (1/4 - u_q)/alpha`, `u_q` the bound of Theorem 4.7, above `2.0 * 10^(-5)` at `125`.

**Proof.** Theorem 5.1 with Facts 6.1, 6.3 and 6.4 and Theorem 4.7: at base `21`, `1/4 - alpha_1 > 2.29 * 10^(-5)` and `alpha < 0.984`; at base `34`, the digit `16` has `1/4 - alpha_1 > 6.299 * 10^(-4)` and `alpha < 0.9916`; at `125`, `1/4 - u_q >= 2.0 * 10^(-5)`. □

**Proposition 5.3 (nothing unconditional).** For any `B` with `max_theta |S_mu(y, theta)| <= B(y)`, steps 1, 2 and 4 return nothing better than `|M_F(x)| << B(x)`, since step 1 pays `q^(-L) k^L Sigma_L(0) >= 1` against `B` by Lemma 2.7; the unconditional `B(y) = y (log y)^(-A)` of Davenport, read at source in Porritt (2018), returns `x (log x)^(-A)`, above `A_F(x)` by the power `x^(1 - alpha)`.

## The census

Every number below is interval arithmetic with directed rounding: the transform evaluated at multiples of `2 pi/Q`, `Q = 8 m q^n`, from enclosures of the cosines at `96` bits, every operation rounded outward, the Perron root bounded by the test of Proposition 3.1 with the float Perron vector, lower bounds truncated down and upper bounds rounded up. `(n, m)` is the window length in digits and the sub-scan.

**Fact 6.1 (the first base).** Base `21` with the digit `0` removed satisfies `alpha_1 in [0.2499765, 0.2499771]` at `(6, 8)`, clearing `1/4` by more than `2.29 * 10^(-5)`, which is `1/4 - 0.2499771`; the same cell reads `[0.2499715, 0.2499821]` at `(5, 8)` and `[0.2498658, 0.2500871]` at `(4, 8)`, which straddles the bar, and a second implementation returns `[0.2499715, 0.2499822]` at `(5, 8)`. Generators: `coprime/lab/py/digit-transform-norms`, verbs `six` and `least`; `lab/py/mobius-region`, verb `params 21 123456789abcdefghijk`.

**Fact 6.2 (nothing below it).** Every one of the `108` distinct sets of the bases `3 <= q <= 20` satisfies `alpha_1 > 1/4` at `(4, 8)`, the lower test run on the restricted matrix of Proposition 3.1; the closest is base `20` missing `0` at `[0.2528608, 0.2531118]`, and the only sets with a lower bound under `0.27` are the digit `0` of the bases `16` to `20`, at `0.2678855`, `0.2635713`, `0.2596685`, `0.2561148` and `0.2528608`, each rerun at its mirror digit, whose transform is the same and whose Lipschitz correction is not, with intersecting brackets. The constant-phase digits, whose cells can empty a row, are all bracketed, base `5` missing `2` at `[0.4366508, 0.4496785]` and base `19` missing `9` at `[0.2732272, 0.2735182]`. So `21` is the least base carrying a one-missing-digit set below `1/4`; Fact 6.9 certifies `106` of the `108` sets a second time by a route that builds no window matrix. Two cross-checks come from the pair route of Section 7: at base `10` its threshold `beta` is certified above `1/4` at every digit, lower bounds `0.2502716` to `0.2541480`, and `beta <= alpha_1` since the `t = 1` term of the infimum defining `beta` is the exponent of one shift of the grid sum; its region census refutes the bases `3` and `4` and the digits `0` and `1` of base `5`, where `alpha/2 > 1/4` makes every refutation force `alpha_1 > 1/4`. Generators: `coprime/lab/py/digit-transform-norms`, verb `floor`; `lab/py/mobius-region`, verbs `threshold` and `criterion`.

**Fact 6.3 (the first family).** All `17` distinct sets of base `34` satisfy `alpha_1 < 0.2493701` at `(4, 8)`, the worst at the digit `16` with `[0.2493107, 0.2493701]` and the best at the digit `0` with `[0.2246400, 0.2247052]`. Every base `21 <= q <= 33` carries a digit with `alpha_1 > 1/4` at `(4, 8)`, the digit of largest upper bound at three window digits, the closest being base `33` missing `15` at `[0.2506145, 0.2506783]` and the next base `32` missing `10` at `[0.2528988, 0.2529666]`. Generators: `coprime/lab/py/digit-transform-norms`, verbs `least` and `floor`.

**Fact 6.4 (the family closed).** All `3663` distinct sets of the `91` bases `35 <= q <= 125` satisfy `alpha_1 < 1/4`, each base at the shortest window in `{2, 3, 4}` that clears at sub-scan `8`: bases `35` to `57` need three digits and every `q >= 58` clears at two, the tightest row being base `58` missing `28` at `alpha_1 < 0.2499305`. Generator: `coprime/lab/py/digit-transform-norms`, verb `family`.

**Corollary 6.5 (the family floor).** Every base `q >= 34` satisfies `alpha_1(q, a_0) < 1/4` at every excluded digit, and `34` is the least base at which every digit clears: every base below it has a digit above the bar.

**Proof.** Facts 6.3 and 6.4 and Theorem 4.7 from `125` for the first claim; Fact 6.2 below `21` and the witnesses of Fact 6.3 from `21` to `33` for the second. □

| `q` | digit `0` | middle digit `floor(q/2)` |
| ---: | --- | --- |
| `10` | `[0.3090500, 0.3106411]` | `[0.3498936, 0.3512569]` |
| `14` | `[0.2780846, 0.2787259]` | `[0.3129308, 0.3134889]` |
| `18` | `[0.2596685, 0.2599987]` | `[0.2909812, 0.2912724]` |
| `20` | `[0.2528608, 0.2531118]` | `[0.2828721, 0.2830946]` |
| `21` | `[0.2498658, 0.2500871]` | `[0.2664296, 0.2666524]` |
| `22` | `[0.2470967, 0.2472930]` | `[0.2760085, 0.2761833]` |
| `26` | `[0.2377912, 0.2379193]` | `[0.2649357, 0.2650507]` |
| `30` | `[0.2305272, 0.2306165]` | `[0.2563012, 0.2563819]` |
| `33` | `[0.2260074, 0.2260778]` | `[0.2401948, 0.2404222]` |
| `34` | `[0.2246400, 0.2247052]` | `[0.2493107, 0.2493701]` |

**Table 1.** Certified brackets at `(4, 8)` at the end digit and the middle digit of ten bases. At odd `q` the middle digit is the constant-phase digit, where `(q-1)|hat F(t)| = |D_q(t) - 1|` has real zeros and rows of `M^-` empty, and its lower bound is the restricted test's. The middle digit is not the dearest: base `20` missing `6` reads `[0.2831944, 0.2834134]` and base `26` missing `8` reads `[0.2652263, 0.2653395]`, each above its middle bracket. Generator: `coprime/lab/py/digit-transform-norms`, verbs `least` and `floor`.

**Fact 6.6 (the published exponent is an upper bound).** Base `10` missing `5` satisfies `alpha_1 in [0.3505775, 0.3505797]` at `(7, 16)` with the wider box `2/10^n` of the sieve literature, strictly below the published `27/77 = 0.3506494`. Generator: `coprime/lab/py/digit-transform-norms`, verb `six`.

**Fact 6.7 (the ceiling of the majorant).** Run the window machine on the majorant `u_q` of Lemma 4.1 in place of `|hat F|`, cell suprema in closed form and rounded outward. The digit-blind bound clears `1/4` first at `q = 75` at three window digits and at every base from `75` to `89`; it reads `0.339085` at base `21` and `0.297296` at base `34` at four digits, against the per-digit `0.250088` and `0.249371` at the digits `0` and `16`. So `75` is where the majorant stops working numerically and `125` where Theorem 4.7 proves it works; the exact crossing depends on the window and no value is claimed for it. Generator: `coprime/lab/py/digit-uniform-bound`, verbs `least`, `window` and `compare`, the per-digit column carried into `compare` from `coprime/lab/py/digit-transform-norms`.

**Fact 6.8 (the floor of the minorant).** `D_q` is real, so `|D_q(t) - e(ct)| >= ||D_q(t)| - 1|` at every phase and `ell_q = ||D_q| - 1|/(q-1)` is a minorant of `|hat F|` at every excluded digit. Through `M^-` at `(4, 8)` it certifies `alpha_1 > 1/4` at every digit at once for every `3 <= q <= 8`, the lowest `0.2522652` at `q = 8`; through `M^+` it puts the `l^1` exponent of the level products of `ell_q` under `1/4` at every `9 <= q <= 20`, `0.2463803` at `9` and `0.2041911` at `20`, so no bound through this minorant reaches the bar there, and Fact 6.2 stays a census of one certificate per set. Generator: `coprime/lab/py/digit-transform-norms`, verb `blind`.

**Fact 6.9 (the shift floor).** By the display in the proof of Lemma 2.6, `Sigma_(N+M)(x) >= Sigma_N(q^M x) min_y Sigma_M(y)`, so `min Sigma_N` is supermultiplicative and `alpha_1 >= (1/N) log_q min_x Sigma_N(x)` at every `N`. In interval arithmetic over one period, each factor at digit `j` lowered by `q^j h/2` times the Lipschitz constant `(2 pi/k) sum_(a in F) |a - (q-1)/2|`, `h` the spacing of the shifts, `min_x Sigma_1 > q^(1/4)` at every digit of every `3 <= q <= 15`, and `min_x Sigma_2 > q^(1/2)` at every digit of `16 <= q <= 18` and at every digit of `19` and `20` except the end digit, the tightest being base `18` missing `0` at `alpha_1 > 0.2501422`. One digit cannot go further: `Sigma_1(0) = 2` exactly, and near the grid `Sigma_1(x) = 2 + (pi^2 C/(q-1)) (qx)^2 + O(x^4)` with `C > 0` exactly at the interior digits other than the constant-phase digit and `C < 0` at the end and constant-phase digits, which therefore dip below `2`, the proof and `C` being on [the coprimality page](../../coprime/notes/coprime.md). Generator: `coprime/lab/py/digit-transform-norms`, verbs `shift` and `grid1`.

The figure is the window machine on the supremum side run by its own census writer, `paper-first-base-below-a-quarter`, over the `395` distinct sets of the bases `10` to `40` at windows of three to five digits and sub-scan `4`: a dash per set at the float Collatz-Wielandt bound, `163` under the line. The writer asserts those counts, that no dash of the bases `10` to `20` is under the line, that base `21` has one, that each of `21` to `33` keeps one above and that every base from `34` has all below; float readings, not certificates.

## What fails or is conditional

Four things. First, Theorem 5.1 is conditional and Proposition 5.3 says the same expansion returns nothing without the hypothesis; the saving at base `21` is more than `2.29 * 10^(-5)` in the exponent, and the Mobius census on [the Mobius page](../notes/meter-transfer-census.md) neither confirms nor threatens it. Second, the unconditional route is a program:

**Conjecture 7.1 (the pair route).** Let `F` be a digit set at base `q` such that every prime dividing the gcd of the differences of its digits divides `q`, with `alpha_1 < 1/4`. Then `sum_(n in S_F, n <= x, gcd(n, q) = 1) mu(n) = O_B(A_F(x) (log x)^(-B))` for every `B > 0`.

The route follows the proof of Maynard (2019) with `mu` in place of the primes and is laid out on the Mobius page: the major arcs for `mu` at base-smooth moduli, free of exceptional zeros since every real primitive character of base-smooth modulus has conductor dividing `8 rad(q)`, and the level of distribution on an initial segment are proved there, the bilinear half never sees the coefficients, the region asked for is `alpha_1 < alpha/2` and `beta <= min(1/4, (2/5)(1 - alpha_1))`, of which `alpha_1 < 1/4` is the `t = 1` proxy, and owed is the level of distribution at moduli sharing a factor with the base, which is why `gcd(n, q) = 1` stays in the statement. The gate `k >= q^(3/4)`, the Parseval floor `beta >= 1 - alpha` against `beta <= 1/4`, is met by every one-missing-digit set of base `q >= 4`, the conclusion would save a power of `log x` only, and the route is refuted at base `10` at every digit, where `beta > 1/4` is certified. Third, the floor has no digit-blind proof: the minorant of Fact 6.8 dies at base `9`, and the shift floor of Fact 6.9 is itself one certified inequality per set, which leaves the end digits of `19` and `20` to the window census. Fourth, the uniform bound is never the sharper number, Fact 6.7.

## Open problems

Four things. The gap between `75` and `125`: the largest loss in Section 4 is the peel, and treating the two-run terms directly should move `125` toward `75`. A proof of the floor: if `min_x Sigma_1 = 2` at every interior digit other than the constant-phase digit, as a float scan to base `40` says, then `alpha_1 >= log_q 2 > 1/4` there for `q <= 15` by one inequality; the end and constant-phase digits dip below `2`, and past `15` any such proof needs two digits of the phase. The exponent itself: every statement is a bracket, and whether `alpha_1` is ever rational or algebraic is untouched. And whether the end digit is always the cheapest is not claimed; the middle digit is not the dearest (Table 1), and where the dearest digit sits is not claimed either.

## Reproducibility

Three studies print every number of Sections 4 and 6, each run from the repository root with one verb that asserts its self-checks, each derived apart from the sweep it checks. `uv run python research/coprime/lab/py/digit-uniform-bound/ubound.py check` in under half a minute covers the domination, the run identity, the peel, the Lebesgue scan, the threshold with its certificate and the consistency rows; its verbs `threshold`, `lebesgue`, `least`, `window` and `compare`, each under ten seconds, print Facts 4.8 and 6.7. `uv run python research/coprime/lab/py/digit-transform-norms/norms.py floor` in `65` seconds prints Fact 6.2 and the witnesses of Fact 6.3, `blind` in `3` seconds Fact 6.8, `shift` in `15` seconds and `grid1` in one second Fact 6.9, `least` in `93` seconds base `34` of Fact 6.3, the shorter windows of Facts 6.1 and 6.6 and Table 1, `six` in `200` seconds Fact 6.6 and the six-digit bracket of Fact 6.1, `family` in `151` seconds Fact 6.4. `uv run python research/mobius/lab/py/mobius-region/mobius_region.py params 21 123456789abcdefghijk` in `12` seconds prints the second bracket of Fact 6.1, `threshold` in `41` seconds the base-`10` refutation, `criterion` in `155` seconds the census of `49` designs. The figure is `figures/paper-first-base-below-a-quarter.ts`, drawn in the browser in under a second a theme, from the numbers its census writer `research/shared/census/src/bin/paper-first-base-below-a-quarter.rs` commits.

## References

- Maynard 2019, Primes with restricted digits, Invent. Math. 217, 127-218. [link.springer.com](https://link.springer.com/article/10.1007/s00222-019-00865-6)
- Maynard 2022, Primes and polynomials with restricted digits, Int. Math. Res. Not. 2022, 10626-10648. [doi.org/10.1093/imrn/rnab002](https://doi.org/10.1093/imrn/rnab002)
- Karwatowski 2022, Primes with one excluded digit, Acta Arith. 202, 105-121, not read at source. [doi.org/10.4064/aa191002-26-8](https://doi.org/10.4064/aa191002-26-8)
- Karwatowski, base 9, Digits of primes in base b = 9, preprint. [math.hhu.de](https://www.math.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Mathematik/20_Institut-Lehrstuehle/5_Algebra_und_Zahlentheorie/Karwatowski/Digits_of_primes_in_base_b_9.pdf)
- Granville 2024, Missing digits, and good approximations, Bull. Amer. Math. Soc. 61. [arxiv.org/abs/2308.03126](https://arxiv.org/abs/2308.03126)
- Baker and Harman 1991, Exponential sums formed with the Mobius function, J. London Math. Soc. 43, 193-198. [doi.org/10.1112/jlms/s2-43.2.193](https://doi.org/10.1112/jlms/s2-43.2.193)
- Zhang 2024, On an exponential sum related to the Mobius function, Proc. Amer. Math. Soc. 152, 1373-1376. [arxiv.org/abs/2204.04613](https://arxiv.org/abs/2204.04613)
- Porritt 2018, A note on exponential-Mobius sums over F_q[t], Finite Fields Appl. 51, 298-305. [doi.org/10.1016/j.ffa.2018.02.005](https://doi.org/10.1016/j.ffa.2018.02.005)
- Chow, Varju and Yu 2024, Counting rationals and Diophantine approximation in missing-digit Cantor sets. [arxiv.org/abs/2402.18395](https://arxiv.org/abs/2402.18395)
- Erdos, Mauduit and Sarkozy 1998, On arithmetic properties of integers with missing digits I, J. Number Theory 70, 99-120. [doi.org/10.1006/jnth.1998.2229](https://doi.org/10.1006/jnth.1998.2229)
- Konyagin 2001, Arithmetic properties of integers with missing digits: distribution in residue classes, Period. Math. Hungar. 42, 145-162. [link.springer.com](https://link.springer.com/article/10.1023/A:1015256809636)
- Nath 2024, Primes with a missing digit: distribution in arithmetic progressions and an application in sieve theory, J. London Math. Soc. 109, e12837. [arxiv.org/abs/2108.09212](https://arxiv.org/abs/2108.09212)
- Leng and Sawhney 2025, Vinogradov's theorem for primes with restricted digits. [arxiv.org/abs/2409.06894](https://arxiv.org/abs/2409.06894)
- Kim 2024, The divisor function over integers with a missing digit. [arxiv.org/abs/2411.09076](https://arxiv.org/abs/2411.09076)
- Collatz 1942, Einschliessungssatz fuer die charakteristischen Zahlen von Matrizen, Math. Z. 48, 221-226. [doi.org/10.1007/BF01180013](https://doi.org/10.1007/BF01180013)
- Wielandt 1950, Unzerlegbare, nicht negative Matrizen, Math. Z. 52, 642-648. [doi.org/10.1007/BF02230720](https://doi.org/10.1007/BF02230720)
