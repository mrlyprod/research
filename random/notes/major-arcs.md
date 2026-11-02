---
title: The major arcs
lead: The circle method's major arcs, summed along the trees of a design's periodic orbits, predict the bias of every digit-parity design along the squares. The coefficient of an orbit point c/q is a product, the Gauss sum G(-c, q)/q times a B-adic integral, the mean of the design's Riesz product over the B-adic squares; for Thue-Morse it is (4/3) e(-1/3) exactly. Against the census to N = 2^32, descriptively, the observed bias has the model's sign in all 180 cells where the model exceeds 2 in z, from 21 designs. The mirror pair 6:5, 6:12345 computes to coefficients 3.666 against 0.438, so on the evidence the mirror split is a factor and not a zero; the twins 5:34 and 5:234 lose every arc of period 1 or 2 below 1/4 to a Gauss sum that vanishes on the whole tree. Along the cubes the cubic Gauss sum kills Thue-Morse's orbit, no lag passes 3.91 in z to 2^33, the two designs with the slowest arcs, 6:135 and 6:345, keep the model's sign at all twelve checkpoints, and 3:12 holds a bias the model misses.
figure: research-random-major-arcs
---

[Along the squares](squares.md) found that a lit statistic of a design `b:D` decays at the rate `N^(-2 eta)` its Fourier exponent sets, and that whether it lights at all is a coefficient no spectral invariant fixes: the mirror pair `6:5`, `6:12345` shares one spectral measure and reads `z_B = 21.28` against `-1.02` at `N = 2^30`. That chapter wrote the heuristic behind the rate, the circle method run on the design, and left its coefficient as a door. This chapter computes the coefficient. The circle method is summed exactly over every frequency whose orbit under multiplication by `b` ends on a given periodic point, which turns the conditionally convergent singular series into one product per point: a Gauss sum at the point, and an integral over the `B`-adic integers that reads the design's digit product on the `B`-adic squares. The model it gives is then tested against a fresh census of all 57 designs along `n^2`, `n < 2^32`, and, out of sample, along `n^3`, `n < 2^33`.

The model gets the sign of the bias right wherever it predicts more than 2 in `z`, 180 cells of 180, and its size to within 20 to 30 percent where the bias is strongest. The evidence goes against one hypothesis: that the mirror pair splits because the coefficient of `6:12345` vanishes. That coefficient computes to `0.4377`, against `3.6656` for `6:5`, and the gap opens at the last digit of a square, which is never 5 in base 6. The vanishing that does happen is the Gauss sum's: at odd `b` an orbit point with denominator `q = 2 (mod 4)` carries a zero quadratic Gauss sum, and the zero propagates to every frequency of its tree, which removes every arc of period 1 or 2 below `1/4` of the twins `5:34` and `5:234`. The generator is the study `lab/rs/major-arcs`, verbs `coefficients`, `model`, `terms` and `cubes`. The notation follows [along the primes](primes.md): `S_f(alpha; M) = sum_(m<M) f(m) e(m alpha)`, `e(x) = exp(2 pi i x)`.

## The model

A design has the sign `x(m) = prod_i v(m_i)` over the base-`b` digits of `m`, with `v(d) = (-1)^[d in D]` and `v(0) = 1`, and the reading along the `k`-th powers is `S(N) = sum_(n<N) x(n^k)`, `k = 2` or `3`. The `k`-th powers near `m` have density `w_k(m) = m^(1/k - 1)/k`, and the degree-`k` Gauss sum is `G_k(a, q) = sum_(r mod q) e(a r^k/q)`. The circle method replaces the indicator of the `k`-th powers by its singular series, and so `S(N)` by

`S(N) ~ 1 + sum_(alpha = a/q) (G_k(-a, q)/q) U(alpha; N^k)`, with `U(alpha; M) = sum_(1<=m<M) x(m) e(m alpha) w_k(m)`.

The `1` is `n = 0`. Each `U` is a smoothed twisted sum of the design, of the size of `S_x(alpha; M)` divided by `M^(1 - 1/k)`, and the sum over `alpha` runs over all reduced fractions. As a series in `q` it does not converge absolutely, and the order of summation is the content of the model.

**The frequencies of denominator dividing `Q` sum to the local density of the `k`-th powers modulo `Q`: `sum_(c mod Q) (G_k(-c, Q)/Q) e(cm/Q) = rho_Q(m) = #{r mod Q: r^k = m (mod Q)}` for every integer `m`. Proved.** Expand the count: `rho_Q(m) = Q^-1 sum_(c mod Q) sum_(r mod Q) e(c(r^k - m)/Q)`, and for `c/Q = a/q` in lowest terms the inner sum over `r` is `(Q/q) G_k(a, q)`. Replacing `c` by `-c` gives the identity, and grouping the `c` by their reduced denominator gives `sum_(q | Q) sum_(a mod q, (a,q)=1) (G_k(-a,q)/q) e(am/q) = rho_Q(m)`. The test `local_sums_count_the_square_roots` checks it at `Q = 36` for every residue.

So the model cut to the denominators dividing `Q` is `sum_(1<=m<N^k) x(m) rho_Q(m) w_k(m)`: the `k`-th powers spread over the residue classes modulo `Q` with their local densities and over the reals with their archimedean density, and the design summed against that measure. Every `Q` gives a well-defined model; as `Q` grows through a divisibility chain the measure sharpens, and the question is how the design's sum against it behaves.

### Trees

Multiplication by `b` acts on the frequencies, and a rational `alpha` is eventually periodic: `b^j alpha` lands on a periodic orbit whose denominator is prime to `b`. Fix `B = b^t` and a point `gamma = c/q` with `q | B - 1`, so that `B gamma = gamma (mod 1)`. The tree of `gamma` is the set of `alpha` with `B^J alpha = gamma (mod 1)` for some `J`; cut at depth `J` it is the coset `{alpha: B^J alpha = gamma}`, with `B^J` elements.

**The tree of `gamma` sums to a Gauss sum at the point times the `B`-adic local density: `sum_(alpha: B^J alpha = gamma) (G_k(-a_alpha, q_alpha)/q_alpha) e(m alpha) = (G_k(-c, q)/q) e(m gamma) rho_(B^J)(m)`. Proved.** On `(q B^J)^-1 Z/Z` the indicator of `B^J alpha = c/q` is `q^-1 sum_(u mod q) e(u (B^J alpha - c/q))`, so by the local identity the left side is `q^-1 sum_u e(-uc/q) rho_(q B^J)(m + u B^J)`. Since `q` is prime to `B`, the Chinese remainder theorem splits `rho_(q B^J)(m + u B^J) = rho_(B^J)(m) rho_q(m + u)`, using `B^J = 1 (mod q)`, and `sum_u e(-uc/q) rho_q(m + u) = sum_(r mod q) e(-c(r^k - m)/q) = e(cm/q) G_k(-c, q)`.

The twist `y(m) = x(m) e(m gamma)` is a design in base `B`: `B = 1 (mod q)` makes `e(m gamma)` a product over the base-`B` digits of `m`, so `y(m) = prod_i y(m_i)` with the tile `y(t) = x(t) e(t gamma)`, `t < B`. Write `sigma = sum_(t<B) y(t)` and `phi(0) = sigma/B`. The tree of `gamma` cut at depth `J` therefore contributes `(G_k(-c, q)/q) sum_m y(m) rho_(B^J)(m) w_k(m)`. Splitting `m = B^J m' + r` with `r < B^J`, `y(m) = y(m') y(r)`, and `sum_(r<B^J) y(r) rho_(B^J)(r) = sum_(n mod B^J) y(n^k mod B^J)`, while `sum_(r<B^J) y(r) = sigma^J`. The model replaces this sum by `(G_k(-c, q)/q) K_J U(gamma; N^k)` with

`K_J = sigma^(-J) sum_(n mod B^J) y(n^k mod B^J) = E_(n mod B^J) prod_(i<J) y(digit_i(n^k))/phi(0)`.

`K_J` is the mean of the Riesz product `R_J(m) = prod_(i<J) y(m_i)/phi(0)`, whose mean over uniform `m` is 1, taken over the `B`-adic `k`-th powers instead. The replacement is not exact. It moves the blocks below `B^J` and the change of `w_k` across each block, an amount bounded in `N` for fixed `J` that grows with `J`: `coefficients` prints the gap for Thue-Morse at `1/3` as `0.386`, `2.331` and `4.230` at `J = 1, 4, 8`, nearly constant in `M`, which at `M = 2^24` is `2.0e-3`, `7.2e-3` and `1.1e-2` of the tree, against `2.5e-2`, `8.1e-2` and `9.9e-2` at `M = 2^12`. The arc term is therefore the model's definition, not a theorem; the gap itself is Verified. With `K = lim_J K_J`, a `B`-adic integral where the limit exists, the coefficient of the arc at `gamma` is

`A(gamma) = (G_k(-c, q)/q) K`.

The model computed here, called the period-2 model below, is `S(N) ~ 1 + Re sum_gamma A(gamma) U(gamma; N^k)` over the points `c/(b^2 - 1)`, which are the orbits of period 1 and 2 under `b`, with `B = b^2`. Every periodic orbit has a tree of its own; the longer ones are left out, and the next sections say what that costs.

### The size of an arc

`U(gamma; BM) = sigma B^(1/k - 1) U(gamma; M) + O(1)`: split `m = B m' + d`, use the base-`B` multiplicativity, and compare `w_k(B m' + d)` with `B^(1/k - 1) w_k(m')`, whose difference is summable. Since `abs(sigma) = B^(1 - eta(gamma))`, with `eta(gamma) = 1 - log_B abs(sigma)` the exponent of [the squares chapter](squares.md), `abs(U(gamma; N^k))` grows as `N^(1 - k eta(gamma))` times a function periodic in `log_B N^k`. The arc stands above the noise `sqrt(N)` exactly when `eta(gamma) < 1/(2k)`: `1/4` along the squares and `1/6` along the cubes. For `B = p^t` the coefficient `K` exists exactly when `eta(gamma) < 1/k`, as the closed form below shows.

`U` is summed exactly up to a Taylor tail: `[1, M)` splits into base-`b` blocks `[c b^j, (c + 1) b^j)`, `w_k` is expanded on each block about its midpoint, and the 44 centred moments of `x(r) e(r gamma)` over `r < b^j` follow a digit recursion from level `j - 1` to `j`. The test `weighted_sums_match_the_direct_sum` checks `U` against the direct sum for four designs, both powers and six values of `M` up to `300001`.

## The B-adic integral

### Prime-power bases

For `b <= 5` the base `B = b^2` is a prime power `p^t`, and `K` is a finite sum. Write `n = p^j u` with `u` a unit, so that `j` has the probability `(1 - 1/p) p^(-j)` and `n^k = p^(kj) w` with `w` uniform on the `k`-th power units: for `k = 2`, `w mod p` a nonzero square at odd `p` and `w = 1 (mod 8)` at `p = 2`; for `k = 3`, `w = +-1 (mod 9)` at `p = 3` and any unit at `p = 2` and `5`. The digits of `w` beyond those congruences are uniform, and a uniform digit has `E y(t)/phi(0) = 1`, so only the first two base-`B` digits of `p^e w` matter. With `r_0 = B/sigma` the weight of a zero digit:

**For `B = p^t`, `K = (1 - 1/p) sum_(j<t) p^(-j) r_0^floor(kj/t) E[R(p^(kj mod t) w)] / (1 - r_0^k/B)`, and the series converges exactly when `eta(gamma) < 1/k`. Proved.** The lowest `floor(kj/t)` base-`B` digits of `n^k` are zero, each contributing `r_0`, and `j -> j + t` multiplies the weight by `p^(-t) r_0^k`. The test `closed_forms_match_the_local_sums` checks the formula against the brute-force `K_J` for five designs.

**For Thue-Morse at `gamma = 1/3`, `K = (4/3)(1 - e(1/3))` and `A(1/3) = (4/3) e(-1/3)`; the period-2 model of `S(N) = sum_(n<N) t(n^2)` is `1 + (4/3) sum_(1<=m<N^2) t(m) cos(2 pi (m - 1)/3)/sqrt(m)`. Proved.** In base 4 the tile of `t` is `(1, -1, -1, 1)`, so at `1/3` it is `y = (1, -omega, -omega^2, 1)` with `omega = e(1/3)`, `sigma = 3` and `r_0 = 4/3`. Here `t = k = 2`, so every `j` has `n^2 = 4^j w` and `K = (1/2) E[R(w)] sum_j (r_0/2)^j = (1/2) E[R(w)]/(1 - 2/3)`. A unit square `w = 1 (mod 8)` has the base-4 digits `1`, then `0` or `2` with equal chance, then uniform ones, so `E[R(w)] = r(1)(r(0) + r(2))/2 = (8/9)(1 - omega)` and `K = (4/3)(1 - omega)`. With `G(-1, 3)/3 = -i/sqrt(3)` and `1 - omega = sqrt(3) e(-1/12)`, `A = (4/3) e(-1/4 - 1/12) = (4/3) e(-1/3)`. The point `2/3` carries the conjugate, and `2 Re[(4/3) e(-1/3) t(m) e(m/3)] / (2 sqrt(m))` is the summand displayed. `coefficients` prints `K = 2.000000000000 - 1.154700538379i` against the closed form and the brute-force `K_J` at `J = 2, 4, ..., 10`, from `1.33333 - 0.76980i` to `1.97399 - 1.13968i`, closing at the rate `sqrt(B)/sigma = 2/3` per digit.

Up to period 14, the end of the search in `coefficients`, Thue-Morse has no other orbit below `1/4`. Past it the set is infinite: an independent recomputation finds orbits below `1/4` from period 19 on, `eta = 0.2460`, `0.2423` and `0.2393` at periods 19, 21 and 23, falling with the period as they shadow `1/3` longer, each with a tiny coefficient, `abs(A) = 0.0026` at `349610/524287`. The displayed sum is Thue-Morse's period-2 model, and every lit design has such a tail of shadowing orbits.

### Base 6

`B = 36` is not a prime power, and the 2-adic and 3-adic parts of `n` interact through the base-6 digits. Write `n = 2^(j2) 3^(j3) u`; then `j = min(j2, j3)` shifts `n^2` by `j` zero base-36 digits, and the rest is `4^e u^2` or `9^e u^2`, `e = abs(j2 - j3)`, uniform on a congruence class modulo `2^(2e+3) 3` or `8 . 3^(2e+1)`. Summing over `j`,

`K = (1/3) [E_0 + sum_(e>=1) (2^(-e) E_e^(2) + 3^(-e) E_e^(3))] / (1 - r_0/6)`,

with `E_0`, `E_e^(2)`, `E_e^(3)` the means of the Riesz product over the classes `{s = 1 (mod 24)}`, `{s = 4^e (mod 2^(2e+3)), s = 1 (mod 3)}` and `{s = 9^e (mod 3^(2e+1)), s = 1 (mod 8)}`. Each class mean is a digit recursion: a base-36 digit is fixed modulo 4 by the 2-adic constraint and modulo 9 by the 3-adic one, and the constraint passes to the next digit through `s -> (s - d)/36`. The cubes take the same form with `j` split by parity and the classes of `8^e u^3` and `27^e u^3`.

**The valuation formula is exact at every finite depth: cut at `J` base-36 digits it reproduces the brute-force `K_J` for `6:5` and `6:12345` at `J = 1..4` and for `6:135`, `6:14` and `6:23` at `J = 1..3`, to `1e-6`. Verified.** The test `six_valuations_match_the_local_sums` runs both routes, which share only the tile; `six_classes_match_the_enumeration` checks the class means against an enumeration of the classes. The limit is the same formula with the series in `e` cut at `e = 14`, where its terms have fallen to `6.1e-5` for `6:5`, `1.2e-6` for `6:12345` and `4.7e-5` for `6:23`, geometrically at about one half per step. No bound on the tail is proved: the trivial bound on a class mean grows as `(36/abs(sigma))^e`, faster than `2^e`, so the base-6 limits in the table below, and the statement that `K` of `6:12345` is not zero, are a Conjecture resting on the computation.

| design | `gamma` | `eta` | `abs(G)/q` | `abs(K)` | `abs(A)` |
|---|---|---|---|---|---|
| `2:1` | `1/3` | 0.2075 | 0.5774 | 2.3094 | 1.3333 |
| `4:13` | `2/5` | 0.1891 | 0.4472 | 1.7488 | 0.7821 |
| `4:23` | `1/5` | 0.1891 | 0.4472 | 0.6000 | 0.2683 |
| `5:12`, `5:23` | `1/4` | 0.2032 | 0.7071 | 0.8498 | 0.6009 |
| `5:14`, `5:124` | `1/3` | 0.2032 | 0.5774 | 1.8028 | 1.0408 |
| `5:34`, `5:234` | `1/6` | 0.2032 | 0 | 0 | 0 |
| `6:5` | `0/1` | 0.2263 | 1 | 3.6656 | 3.6656 |
| `6:3` | `0/1` | 0.2263 | 1 | 2.4761 | 2.4761 |
| `6:1` | `0/1` | 0.2263 | 1 | 1.2000 | 1.2000 |
| `6:4` | `0/1` | 0.2263 | 1 | 1.1183 | 1.1183 |
| `6:12345` | `0/1` | 0.2263 | 1 | 0.4377 | 0.4377 |
| `6:2` | `0/1` | 0.2263 | 1 | 0.2330 | 0.2330 |
| `6:14` | `2/5` | 0.1943 | 0.4472 | 1.8198 | 0.8138 |
| `6:23` | `1/5` | 0.1943 | 0.4472 | 0.5186 | 0.2319 |
| `6:135` | `3/7` | 0.1755 | 0.3780 | 1.6053 | 0.6067 |
| `6:345` | `1/7` | 0.1755 | 0.3780 | 1.4110 | 0.5333 |
| `6:125` | `2/7` | 0.1755 | 0.3780 | 1.1930 | 0.4509 |
| `6:45` | `1/7` | 0.2335 | 0.3780 | 3.0123 | 1.1386 |
| `6:13` | `3/7` | 0.2335 | 0.3780 | 2.3992 | 0.9068 |

The table is the best arc of each design listed, from `coefficients`, which prints every arc below `1/4` of all 57 designs. `K` at the fixed point `0` is real; `6:1` has `K = -1.2000`, negative, so its bias runs against those of `6:3`, `6:4` and `6:5`.

## Zero Gauss sums

**If `G_k(-c, q) = 0` at the orbit point, the frequencies of its tree sum to zero against every `e(m alpha)`, at every depth: the tree contributes nothing to the model. Proved**, by the tree identity, whose right side has the factor `G_k(-c, q)`. The quadratic Gauss sum vanishes exactly when `q = 2 (mod 4)`: by the Chinese remainder theorem it factors over the odd part and the factor 2, and `G(a, 2) = 1 + e(a/2) = 0` for odd `a`. A point `c/q` with `q | B - 1` can have `q = 2 (mod 4)` only at odd `b`.

`coefficients` lists every arc of period 1 or 2 below `1/4` with a zero Gauss sum, enumerating the points in floating point. Along `n^2` they belong to three designs: `1/6` and `5/6` of `5:34` and of `5:234`, at `eta = 0.2032`, and `1/2` of the folded class `3:1`. **The twins `5:34` and `5:234` have no arc of period 1 or 2 below `1/4` with a nonzero coefficient, so their period-2 model has no growing bias, while their twins `5:14` and `5:124` have `abs(A) = 1.0408` at `1/3`. Proved** for that model; the census agrees, the largest `abs(z_B)` of `5:34` over the eleven checkpoints `2^12..2^32` being `1.80`. Past period 2 they are not empty: their period-5 orbits on `q = 3124`, at `eta = 0.2424`, carry Gauss weights at most `0.0253`, and an independent recomputation gives them `abs(A) = 0.040` and `0.046`, about `0.1` in `z` at `2^32`, growing as `N^(0.015)`. The tree statement is stronger than the orbit statement of the squares chapter: the twin `5:34` has the same `eta` as `5:14` at every orbit, and it is the Gauss sum at `1/6`, not any property of the spectrum, that removes the whole family of frequencies `(1/6 + r)/25^J`.

The cubic Gauss sum vanishes on more. At a prime `p` with `3` not dividing `p - 1`, `r -> r^3` permutes the residues and `G_3(a, p) = sum_r e(ar/p) = 0`; at `p = 3`, `r^3 = r (mod 3)` does the same. By the Chinese remainder theorem `G_3(a, q)` is then zero whenever `q` is one of `2, 3, 5, 6`, and along `n^3` the generator lists 27 arcs below `1/4` with a zero cubic Gauss sum: the folded class at `1/2`, Thue-Morse at `1/3`, every arc of `4:13` and `4:23`, of `5:14`, `5:124`, `5:34` and `5:234`, and every fifth of base 6.

## The mirror pair

A design and its mirror `v*(d) = v(b - 1) v(b - 1 - d)` have one spectral measure, and the squares chapter's census reads them far apart: `6:5` lights all seven lags and its bias, `6:12345` neither its bias nor four of its lags. The working hypothesis was that the model's coefficient of `6:12345` vanishes. On the computed coefficients it does not.

**The coefficient of `6:12345` at its best arc, the fixed point `0`, is not zero: `K = 0.4377`, against `3.6656` for `6:5`, at the same `eta = 0.2263`. Conjecture**, resting on the valuation series without a tail bound. The brute-force partial sums over `n mod 36^J`, `J = 1..5`, are exact and Verified: `1.25`, `0.59375`, `0.53516`, `0.50098`, `0.46280`, still falling at `J = 5`; those of `6:5` are `2.25`, `3.125`, `3.45898`, `3.58997`, `3.63857` against `3.66557`.

The mirror changes the integral, not the integrand. **For a mirror pair at the fixed point `0`, `K_J* = E_(n mod b^J) R_J(b^J - 1 - (n^2 mod b^J))`: the mirror's coefficient is the design's own Riesz product read on `-1 - n^2` instead of `n^2`. Proved.** The digits of `b^J - 1 - m` are `b - 1 - m_i`, so `x*(m) = v(b - 1)^J x(b^J - 1 - m)` below `b^J`, and `phi*(0) = v(b - 1) phi(0)`; the signs `v(b - 1)^J` cancel against the normalization. So `6:12345` reads the Riesz product of `6:5`, which penalizes the digit 5, on the law of `-1 - n^2` in the 6-adic integers. A square ends in `0, 1, 3, 4` in base 6 with the chances `1/6, 1/3, 1/6, 1/3`, never in 5, so `6:5` gains the full factor `3/2` on the last digit; `-1 - n^2` ends in `5, 4, 2, 1` with the same chances, and in 5 one time in six, the uniform rate, so `6:12345` gains nothing there. Deeper digits repeat the story with smaller margins, and the product is `3.6656` against `0.4377`.

The first term of the series shows it from the Fourier side. The tree of `0` starts with `q = 1`, the point itself, and continues with every `q` dividing a power of 6, the frequencies `a/q` with `6^J a/q` an integer; its largest terms sit at `q = 3`, `9`, `12` and `27`. At `q = 3` the Gauss sum `G(1, 3) = i sqrt(3)` is purely imaginary, and the term is `(2/sqrt(3)) Im psi(1/3)` times the `q = 1` term, `psi = phi/phi(0)`. For `6:5`, `psi(1/3) = (1 + i sqrt(3))/4`, so the term is half the first; the mirror multiplies `psi(1/3)` by `e(2/3)` and conjugates it, which makes it the real number `-1/2`, and the term is zero. `terms` prints at `N = 2^32` `T_1 = 5.517` and `T_3 = 2.759` for `6:5`, and `T_1 = -0.687` and `T_3 = -0.000` for `6:12345`.

The archimedean factor splits the pair a second time. `U(0; 6M) = (sigma/sqrt(6)) U(0; M) + O(1)` with `sigma = 4` for `6:5` and `-4` for `6:12345`: a factor `1.63` per factor 6 for one, so its smoothed sum grows with a fixed sign, and `-1.63` for the other, whose smoothed sum changes sign at every factor 6 and passes through zero in between. In the census the model of `6:12345` stays within `+-0.62` in `z` over `2^12..2^32`, and the observed bias within `+-1.58`; at the rate `N^(1/2 - 2 eta) = N^(0.047)` the model's amplitude would reach 4 near `N = 2^89`. On the computed coefficients the mirror split is a factor of `8.4` in `K`, compounded by the alternation, and not a zero.

The other mirror pairs show the same freedom. Along `n^2` the largest `abs(A)` below `1/4` of each side is `1.2000` against `1.1183` for `6:1`, `6:4`; `0.2330` against `2.4761` for `6:2`, `6:3`, a factor `10.6`; `0.1753` against `0.7077` for `6:1245`, `6:25`; while at base 5 the pairs `5:12`, `5:23` and `5:14`, `5:124` have equal `abs(A)`, `0.6009` and `1.0408`, and differ only in the phase of `K`.

`6:23`, the quiet design of the squares chapter, has `abs(A) = 0.2319` at `1/5`, against `0.8138` for `6:14` at the same `eta = 0.1943`. Its computed coefficient is small, not zero, and its model stays under 2 in `z` to `2^32`, where it reads `1.63` against an observed `1.16`. The door of the squares chapter, whether that coefficient vanishes exactly or only nearly, points to nearly, under the same Conjecture as `6:12345`.

## Against the census

`model 32` reads all 57 designs along `n^2`, `n < 2^32`, in 79 s on 8 threads, walking `n^2` by its differences in limbs of `3^12`, `5^8` and `6^7` and in `u128` for bases 2 and 4, and compares `z_B = S(N)/sqrt(N)` with the model at the eleven checkpoints `N = 2^12, 2^14, ..., 2^32`. The model keeps every arc below `1/4`. An arc with `eta` between `1/4` and `1/2` has a convergent tree, but a slowly convergent one, whose limit coefficient grows without bound as `eta -> 1/2`; read at small `N` through its limit it produces a spurious bias, `52.54` in `z` at `2^12` for `6:15` in `model 32 0.5`, which keeps the arcs below `1/2`, where the census reads `1.06`. Those arcs fall under the noise at large `N` and are left out. The cut and the sign threshold below were chosen after that first run, which kept every arc below `1/2`: there, of the 226 cells with `abs(z_model) > 2`, 211 have the observed sign, and the model lights six designs the census does not. An independent recomputation finds the same model for every cut between `0.2335` and `0.2675`, the gap between the arcs.

| design | `eta` | `z_B` at `2^32` | model | same sign |
|---|---|---|---|---|
| `2:1` | 0.2075 | -16.47 | -14.01 | 11 of 11 |
| `4:13` | 0.1891 | -32.99 | -28.89 | 11 of 11 |
| `6:135` | 0.1755 | 33.79 | 29.97 | 11 of 11 |
| `6:14` | 0.1943 | -27.19 | -22.29 | 11 of 11 |
| `6:5` | 0.2263 | 23.83 | 20.22 | 11 of 11 |
| `5:124` | 0.2032 | -15.81 | -15.20 | 11 of 11 |
| `6:345` | 0.1755 | 10.68 | 8.08 | 10 of 11 |
| `5:14` | 0.2032 | -7.31 | -9.46 | 11 of 11 |
| `6:3` | 0.2263 | 7.76 | 9.62 | 11 of 11 |
| `6:4` | 0.2263 | 7.63 | 5.76 | 11 of 11 |
| `6:125` | 0.1755 | -7.15 | -8.17 | 9 of 11 |
| `6:13` | 0.2335 | 6.86 | 4.59 | 11 of 11 |
| `6:123` | 0.2142 | 6.27 | 5.77 | 10 of 11 |
| `5:23` | 0.2032 | -5.78 | -5.93 | 10 of 11 |
| `5:12` | 0.2032 | -4.16 | -4.84 | 9 of 11 |
| `6:1` | 0.2263 | -3.18 | -2.59 | 11 of 11 |
| `6:12345` | 0.2263 | -0.26 | -0.30 | 6 of 11 |
| `5:34` | 0.2032 | 0.56 | 0.00 | |

**Descriptively, along `n^2`, `n < 2^32`, the observed bias has the model's sign in `241` of the `253` cells where the model exceeds `0.5` in `abs(z)`, `215` of `219` above `1`, `196` of `197` above `1.5` and `180` of `180` above `2`, from `26`, `24`, `23` and `21` designs; above `2` it lies within a factor 2 in `169`, and the rms of `obs/model - 1` is `0.313`. Over the 594 cells of the 54 non-folded designs at the eleven checkpoints the residual has rms `1.856`, `2.116` over the 198 cells of the designs the model lights and `1.712` over the other 396. Verified.** The threshold 2 is the lowest of the four with no sign failure, and the agreement is not an independent test cell by cell: the 180 cells come from 21 designs whose residual is roughly proportional to the model. An independent recomputation on the odd checkpoints `2^13..2^31` reads `160` of `160` above `2`, `154` within a factor 2. The coin reads `z = 1.36` at `2^32`. The model tracks more than a rate: `5:12` and `5:23`, whose `K` is complex, change sign along the checkpoints, and the model changes sign with them, `5:12` reading `-0.34, -0.86, -2.57, -3.42, -2.34, 0.77, 4.72, 2.38, -2.82, -5.80, -4.16` against `1.87, 1.21, -1.62, -2.97, -1.74, 1.93, 4.90, 2.27, -2.54, -6.83, -4.84`; and `6:125` jumps from `-0.53` at `2^30` to `-7.15` at `2^32`, where the model goes from `-1.35` to `-8.17`.

With lit meaning `abs(z) > 4` at some checkpoint, 18 designs are lit by both, 32 by neither, none by the model alone, and four by the census alone: `4:1`, `3:2` and `4:3`, whose best arcs sit at `eta = 0.2982`, `0.2675` and `0.3398`, above the cut, and `6:25`, which reaches `4.00` against a model of `3.24`. The first three are the transients of the squares chapter, near-quarter arcs whose structured part decays against the noise only as `N^(1/2 - 2 eta)`, `N^(-0.096)` for `4:1`; they are the largest residuals, `7.83` for `4:1`. The split rides on the cut: `6:25` reads `4.0044` and `6:1245` `3.9881` at their largest, and `5:2` `3.9062`.

The model falls short on the strongest biases. The ratio of observed to model at `2^32` is `1.18` for Thue-Morse, `1.14` for `4:13`, `1.22` for `6:14`, `1.18` for `6:5`, `1.13` for `6:135` and `1.32` for `6:345`, and `1.04`, `0.77`, `0.81` and `0.88` for `5:124`, `5:14`, `6:3` and `6:125`. For Thue-Morse the ratio runs `1.50` at `2^24`, `1.26` at `2^30` and `1.18` at `2^32`, a residual that grows more slowly than the main term but stays between `2.2` and `4.4` times `sqrt(N)` from `2^20` to `2^32`. Within the search, period 14 for Thue-Morse and 5 for base 6, Thue-Morse, `6:5` and `6:12345` have no other orbit below `1/4`, and the shadowing orbits past it carry coefficients like Thue-Morse's `0.0026`, so a missing arc is an unlikely source, though not an excluded one; the minor arcs are the other. For other designs further orbits below `1/4` exist, of period 3 to 7, which the model omits: `coefficients` finds them at `q = 1023` and `16383` for `4:13` and `4:23` (`eta = 0.2374`, `0.2236`), at `q = 3124` for the four base-5 designs of the thirds and sixths (`0.2424`), at `q = 215`, `37` and `7775` for the sevenths designs (`0.2225`, `0.2415`, `0.2034` for `6:125`, `6:135`, `6:345`), and at `q = 7775` for `6:14`, `6:23` and four others (`0.2447`, `0.2459`).

**For every design with an arc of period 1 or 2 below `1/4` and a nonzero coefficient, the period-2 model is the leading term of its bias along `n^2`: `S(N) - 1 - Re sum A(gamma) U(gamma; N^2) = o(N^(1 - 2 eta*))`, the sum over those arcs and `eta*` the least of their `eta`. Conjecture.** The full statement would sum over every periodic point below `1/4`, a set that is infinite for every lit design, and would make a bias lit, eventually, exactly when some arc below `1/4` has `A != 0`; for `5:34` and `5:234` that comes only from their period-5 orbits, which at the rate `N^(0.015)` would reach 4 in `z` only near `N = 2^380`. The minor arcs are the obstruction: no bound on them beats the existential saving of the normality proofs (Drmota, Mauduit and Rivat 2019; Mullner 2018, [refs](../refs.md)), and the census can only say that the residual stays a few `sqrt(N)` to `2^32`.

## Term by term

The singular series converges only in the right order. `terms 32 36` prints the model as written, `T_q(N) = sum_(a mod q, (a,q)=1) (G(-a,q)/q) U(a/q; N^2)` for each `q <= 36`, in units of `sqrt(N)`. For Thue-Morse at `N = 2^32` the terms that matter are `T_3 = -1.488`, `T_12 = -5.518` and `T_24 = 1.343`, the first three levels of the tree of `1/3`, whose denominators are `3 . 2^j`; `T_6 = 0`, its Gauss sum vanishing. Their running sum to `q = 36` is `-5.687`, while the trees give `-14.011` and the census `-16.471`: the levels `q = 48, 96, ...` carry the rest, since a cut at `q <= 36` keeps the denominators `3 . 2^j` with `j <= 3` of a tree that converges only at `2/3` per factor 4. For `6:5` the running sum to 36 is `13.234` against the trees' `20.225`, through `T_1 = 5.517`, `T_3 = 2.759`, `T_9 = 1.379`, `T_12 = 2.759` and `T_27 = 0.517`. For `6:12345` it is `-0.779` against `-0.301` and an observed `-0.257`, with `T_3` zero.

The terms off the trees fade: Thue-Morse's `T_17` reads `-0.465` at `2^12` and `-0.020` at `2^32`, the orbit of `1/17` having `eta` above `1/2`. The tree organization is what makes the series a sum of finitely many products, one per periodic point.

## Along the cubes

The cube test was posed as a falsification: Thue-Morse's orbit at `1/3` meets the cubic Gauss sum `sum_(r mod 3) e(a r^3/3) = 0`, and every design of base at most 6 has `eta > 1/6`, so the model predicts that every design reads as noise along `n^3`. The trees sharpen the prediction. Along `n^3` the weight is `w_3(m) = m^(-2/3)/3`, an arc's structured part grows as `N^(1 - 3 eta)`, and the surviving arcs below `1/4` are few: the fixed point `0` of the five base-6 singletons and `6:12345`, at `eta = 0.2263`; the quarters of `5:12` and `5:23`, at `0.2032`, `abs(A) = 0.3214`; and the sevenths, with `abs(A) = 0.9416` for `6:135` and `1.0548` for `6:345` at `eta = 0.1755`, and `0.0377` for `6:125`, whose cubic Gauss sum at `2/7` is small. All fall against `sqrt(N)`, but the sevenths at `0.1755` only as `N^(-0.026)`, so slowly that the prediction at accessible `N` is not noise: it is a bias of order one to two in `z`, of a definite sign.

`cubes 33` reads all 57 designs along `n^3`, `n < 2^33`, so `n^3 < 2^99`, in 239 s, with the bias and the lags `h <= 7` of every design and the 8-block chi-square of nine streams.

**Along `n^3`, `n < 2^33`, no lag passes `abs(z) = 3.91`, the largest over the non-folded designs and the coin, at `6:2345`, `h = 7`; the overlapping 8-block chi-squares of `2:1`, `4:13`, `5:124`, `6:5`, `6:12345`, `6:14`, `6:135`, `6:345` and the coin read `312.3`, `195.9`, `303.7`, `265.8`, `275.8`, `226.6`, `265.5`, `324.5` and `319.6`, all within `1.8` sd of the iid mean 255, the sd `38.7` being that of overlapping windows, `2 sum m^2 = 1498`, as in the squares chapter. Verified.** The lags are not plainly noise: `6:2345` at `h = 7` reads `2.97, 2.15, 2.89, 1.91, 2.88, 3.91` over `2^23..2^33`, one sign for six checkpoints, and `6:23` at `h = 2` reads `-2.98, -3.19, -2.69` over `2^29..2^33`; `cubes` prints every lag row that passes 3 somewhere. Thue-Morse, whose model along the cubes is zero, reads `z_B` between `-4.09` and `0.80` over the twelve checkpoints `2^11..2^33`, ending at `-0.46`.

**Along `n^3` the bias has the model's sign in all 16 cells where the model exceeds 2, and within a factor 2 of it in 10; 13 of the cells belong to `6:135` and `6:345`, which keep the model's sign at all twelve checkpoints `2^11..2^33`, reading between `-5.30` and `-2.84` against a model between `-2.45` and `-1.70`, and between `2.17` and `5.57` against `1.77` to `2.48`, and three are single cells at `2^11`, of `6:13`, `6:5` and `6:4`. Verified.** The threshold cuts through the slow decay of the two sevenths designs, several of whose cells sit within `0.05` of 2. The residual over the 648 cells of the non-folded designs has rms `1.285`. Six designs pass `abs(z) = 4` at some checkpoint: Thue-Morse and `5:2`, whose models are zero, at `4.09` once each, and `6:13`, `6:5`, `6:135` and `6:345`, at `5.44`, `5.21`, `5.30` and `5.57`, each with the model's sign at that checkpoint. **The model misses `3:12`, which reads `+2.35` to `+3.71` at the nine checkpoints `2^17..2^33` against a model of at most `0.02`. Verified.** Its best cubic arc, `1/4` at `eta = 0.2675`, lies above the cut, and its structured part would fall as `N^(-0.30)`, too fast for what the census shows.

So the falsification half-holds. Every chi-square of the nine streams reads as noise to `2^33`, no lag passes 4, and Thue-Morse's bias is noise, its orbit removed by the cubic Gauss sum, consistent with its simple normality along the cubes (Spiegelhofer 2023, [refs](../refs.md)). But "every design reads as noise" is false at the box's reach: the sevenths carry the model's bias, of its sign and at about twice its size, which the model predicts fades only as `N^(-0.026)`, by a factor 2 over some 38 doublings; `3:12` carries a bias the model does not explain; and two lags hold one sign over several checkpoints.

## The figure

The figure draws the object of the coefficient: the 3288 squares modulo `6^6`, placed 6-adically on a `216 x 216` grid, the low base-6 digits choosing the coarse blocks, `col = 36 m_0 + 6 m_2 + m_4` and `row = 36 m_1 + 6 m_3 + m_5`. The eight coarse blocks that hold squares are the eight squares modulo 36. Above the diagonal a square is blue where `6:5` reads `+1` on its six digits and orange where it reads `-1`; below the diagonal the same for `6:12345`. Weighted by their square roots the two sums are `14168` and `2192`, which the figure checks; divided by `4^6` they are the brute-force `K_3 = 3.45898` and `0.53516` of the two designs, and the two halves are the mirror split at depth 3.

## What stays open

- **The minor arcs.** The model is a heuristic until the minor arcs are bounded by more than the existential saving of the normality proofs. The census residual, two to four `sqrt(N)` on the strongest biases, is the size such a bound must reach.
- **The base-6 tail.** The valuation series for `B = 36` converges geometrically in the computation, at about one half per step, but no proved bound controls a class mean better than `(36/abs(sigma))^e`. A sub-action, as for Thue-Morse's `eta`, or an `L^2` bound on the class means, would certify `K` and with it that `6:12345`'s coefficient is nonzero.
- **The longer orbits.** The model keeps the orbits of period 1 and 2. The orbits below `1/4` of longer period that the generator lists, and the infinite family past its search that shadows each best orbit, need `B = b^t` with `t` growing: the prime-power formula handles them in principle, base 6 needs a base-`6^t` valuation recursion, and the sum over the family needs a convergence argument.
- **The near-quarter arcs.** `4:1`, `3:2` and `4:3` light at small `N` through arcs above `1/4`. Their trees, cut at the depth the circle method allows at `N`, should give a finite-`N` model for them; the limit coefficient, used here, does not.
- **The lags.** The lag `h` pairs `n^2` with `(n + h)^2`, whose archimedean density lives on a curve, so the bivariate major arcs do not factor through single twisted sums. The mirror pair shares the linear part of every lag, by its shared spectral measure, so the split of its lags, all seven for `6:5` and three for `6:12345`, lives in the quadratic terms.
- **The cubes.** `3:12` holds `+2.35..+3.71` over nine checkpoints with no arc below `1/4` to explain it, and `6:2345` at `h = 7` holds one sign for six. A census past `2^33` would say whether they fade.
