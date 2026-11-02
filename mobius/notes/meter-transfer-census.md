---
title: The meter, the transfer and the census
lead: The Mobius meter against its design's exact count, the carry-free transfer between scaled digit sets, and the 47-column census of exact Mertens values.
figure: research-mobius
slug: mobius
---

## The meter and its yardstick

- `A_F(x)` counts `S_F` up to `x`. The count is exact at every checkpoint: `A_F(base^level) = fill^level - 1` when `0 in F` (plus 1 when `1 in F` too, for the boundary element `base^level` itself), and `A_F(base^level) = (fill^(level+1) - fill)/(fill - 1)` when `0` is not in `F`, by counting digit strings of each length. **Proved**; the lane's tests pin it against direct enumeration. Between checkpoints `A_F(x)/x^(log(fill) / log(base))` carries the log-periodic ripple every design in this tree carries - the classical fluctuation of digital sums ([Flajolet, Grabner, Kirschenhofer, Prodinger and Tichy 1994](https://doi.org/10.1016/0304-3975(92)00065-Y)) - so a checkpoint value is a grid value, never a constant.
- The meter is `M_F(x) = sum of mu(n)` over `n in S_F`, `n <= x`, and the exponent is `theta(F) = limsup of log|M_F(x)| / log A_F(x)`. A single cut of `|M_F|` is a bad estimator - the meter crosses zero freely - so the census prints two readings per level: `M_F(base^level)` itself, and the running maximum `max of |M_F(x)|` over `x <= base^level`, whose exponent `thetamax` is monotone in the numerator and is the estimator the slope tables use.
- The yardstick matters. `A_F(x)` grows like `x^(log(fill) / log(base))`, so `S_F` is sparse, and a bound of shape `o(x)` is weaker than the trivial `|M_F(x)| <= A_F(x)`. The indicator of `S_F` is automatic in that base, so [Mullner 2017](https://doi.org/10.1215/00127094-2017-0024) (automatic sequences fulfill the Sarnak conjecture) gives `M_F(x) = o(x)` for every `F`: orthogonality holds and the question is well-posed, but against the set's own mass that bound says nothing at all. The same shape repeats in base 2 through circuits: the indicator is computable in bounded depth from the binary digits, so [Green 2012](https://arxiv.org/abs/1103.4991) also gives `o(x)`, again below the trivial bound. **Verified** against the literature. The honest question is `theta`, and it is open at every `2 <= fill <= base - 1`.
- The Dirichlet series over `S_F` is built territory, and this page claims nothing about it: the abscissa is `log(fill) / log(base)` ([Kohler and Spilker 2009](https://doi.org/10.1007/s00591-009-0059-5), with position-varying digit rules in [Nathanson 2021](https://arxiv.org/abs/2010.06295)); the series continues meromorphically to `C` with simple poles among `s = log(fill) / log(base) - m + 2 pi i j / log base` (the automatic-series mechanism of [Allouche, Mendes France and Peyriere 2000](https://doi.org/10.1006/jnth.1999.2487), carried out for missing digits in [Burnol 2026](https://arxiv.org/abs/2602.19727) and unified in [Allouche, Shallit and Stipulanti 2025](https://arxiv.org/abs/2401.13524)); a pole lattice of period `2 pi i / log base` reads as log-periodic oscillation through the Mellin dictionary of [Flajolet, Gourdon and Dumas 1994](https://inria.hal.science/inria-00074307), and the oscillation is visible in the series' own numerical moments ([Burnol 2026 oscillations](https://arxiv.org/abs/2604.24754)); the Mobius function itself is automatic in no base, so the Mobius-weighted series inherits none of that continuation ([Coons 2010](https://doi.org/10.5802/jtnb.718)); and no Mobius or Mertens sum appears anywhere in that literature. **Verified** against the sources in refs.md. The series does not carry the meter the way [zeta](https://en.wikipedia.org/wiki/Riemann_zeta_function) carries Mertens: `S_F` is not multiplicatively closed - at `base 3`, `F = {0,1}`, both `4 = 11` and `13 = 111` lie in `S_F` while `4 x 13 = 52 = 1221` does not - so there is no Euler product and `M_F` is not the coefficient sum of an inverse series. **Proved** by that witness.
- The full digit set is the classical boundary. `S_F` is then every integer, `M_F` is the Mertens function of [Mertens 1897](https://www.zobodat.at/pdf/SBAWW_106_2a_0761-0830.pdf), and `M(x) = O(x^(1/2 + eps))` for every `eps > 0` is equivalent to the Riemann hypothesis ([Titchmarsh 1986](https://sites.math.rutgers.edu/~zeilberg/EM18/TitchmarshZeta.pdf), Theorem 14.25 (C)), while `limsup |M(x)|/sqrt(x) >= 1.06` unconditionally by [Odlyzko and te Riele 1985](https://doi.org/10.1515/crll.1985.357.138), so the exponent over all `x` equals `1/2` exactly when RH holds. **Verified** against the literature. This page claims nothing about RH: the full-set column below is a control rendered for scale, and the tree's own claims live in the restricted columns.
- **The even moments of the digit transform are additive energies. Proved.** `sum_{a mod base^level} |hat F_level(a/base^level)|^(2r) = base^level E_r(level)` with `E_r(level)` the number of `2r`-tuples of digit strings of length `level` with `n_1 + ... + n_r = n_(r+1) + ... + n_(2r) mod base^level`, by orthogonality, and `E_r(level)` is counted by a carry DP on the carry pairs of the two sides, so each moment is C-finite in `level` of order at most `r(r+1)/2` and its growth constant `Lambda(2r) = base rho` is an algebraic number, `rho` the Perron root of [the transfer matrix](https://en.wikipedia.org/wiki/Transfer-matrix_method_(statistical_mechanics)), certified in exact rationals (`lab/rs/rho-decoupling`, the `riesz` module).
- **The fourth-moment constants. Verified.** `Lambda(4) = 18` at `{0,1}`, `{0,2}` and `{1,2}` in base 3 (`rho = 6`), `2(23 + sqrt 353) = 83.5766` at `{0,1,2}` in base 4 (`x^2 - 23x + 44`), `(275 + 5 sqrt 2369)/2 = 259.1809` at `{0,1,2,3}` in base 5 (`x^2 - 55x + 164`), `95` at `{0,2,4}` in base 5, and `6566.412` to `6567.410` over the reflection classes of one excluded digit in base 10; every value sits strictly inside `[max(fill^4, base fill^2), base fill^3]` and a hair above `fill^4` at the dense families (`log_base(Lambda(4)/fill^4)` is `0.107` at base 3, `0.0004` at base 10), and the sixth, eighth and tenth moments at `{0,1}` base 3 are `39 + 3 sqrt 79`, `3(99 + sqrt 5265)/2` and a cubic (`lab/rs/rho-decoupling`).
- **What a moment buys the bilinear sum. Proved.** Holder with the `2r`-th moment on the digit side and Parseval on the bilinear side bounds the Type II sum over `m` up to `M` and `l` up to `N`, `4MN <= x`, by `x^(theta_p/p + 1/2 - 1/p)` with `theta_p = log Lambda(p)/log base`, which is at least `x^(alpha + 1/4)` for every even `p >= 4` and every digit set, above the trivial `x^alpha`; so no moment of the digit transform alone beats the trivial bound, and the route needs the bilinear sum on the minor arcs below its own root mean square, which random-sign coefficients defeat on the census (**Verified**, `lab/rs/rho-decoupling` the `arcs` lines).
- **The multiplicative energy of a digit column has no exponent of its own. Proved.** With `E_x(level) = #{(n_1, n_2, n_3, n_4) in D_level^4 : n_1 n_2 = n_3 n_4}` and `K = fill^level`, the two diagonals give `2K^2 - K <= E_x(level)`, and `E_x(level) = sum_m r(m)^2 <= K^2 max_m r(m)` with `r(m) <= d(m)` gives `E_x(level) = fill^(2 level) x^(o(1))` for every base and digit set (the census reads `58760487` at `{0,1}` base 3, `level 12`, the exponent `1.356938` falling toward `2 alpha = 1.261860`); so a Type II sum estimated through the energy obeys `|Sigma| <= (2MN)^(1/2) x^(alpha/2 + o(1))` and misses the trivial bound by `(1 - alpha)/2`; the excess over the diagonal is structure, not arithmetic: the shift family `(base^i u, base^j v, base^(i') u, base^(j') v)` with `i + j = i' + j'`, counted in closed form when `0` is a digit, is `0.44` of it at `{0,1}` base 3, `level 12` (`lab/rs/rho-decoupling`, the `menergy` module).
- **Above the sup the `L^p` norms of the transform buy nothing. Proved.** The sandwich `max(fill^p, base fill^(p/2)) <= Lambda(p) <= base fill^(p-1)` forces `Lambda(p)^(1/p)/fill` down to `1`, and it reads `1.224744, 1.029883, 1.004288, 1.000653, 1.000107` at `p = 2, 4, 6, 8, 10` for `{0,1}` at `base 3`, so every higher norm is the supremum up to a factor tending to `1` and no ladder of moments reaches past the bullet above (`lab/rs/rho-decoupling`, `riesz higher moments`).
- **The unbalanced kernel carries no Type II estimate uniform over bounded coefficients at any digit set containing `0`. Proved.** At `a_m = b_l = 1` the Type II sum is the box representation count and some admissible box carries `R >= x^(alpha - o(1))`, so the trivial bound is attained and the only target left is the balanced sum; there the route returns the box's own trivial bound, the ratio of bound to trivial rising through `1` (`1.0134` at `level 12` and `1.0730` at `level 14` at `{0,1}` base 3) while the margin of `alpha` over the achieved exponent falls from `0.035675` to `0.027009`, and the digit column is worse for the method than a random column of the same density at every cell of the arc regime (`lab/rs/rho-decoupling`, `menergy type II`).
- **Two box witnesses floor every coefficient-free route at `x^alpha`. Verified at the dense cells.** At `a_m = b_l = 1` the Type II sum reads `0.19` to `0.41` of `fill^level` over the nine dense cells, and at `a_m = 1_(base | m)`, `b_l = 1` the balanced sum still reads `0.0024` to `0.104` of `fill^level` there, carried by the frequencies `a'/base^j` at bounded `j`, which are exactly the major arcs the pair route below removes and computes; the statement for every digit set rests on `R >= fill^level (log x)^(-C)` and stays **Conjecture**, and the witness is void at `alpha = 0.15`, where the box is empty (`lab/rs/rho-decoupling`).
- **The Mobius signs cancel the column's correlation no better than random signs. Verified.** The digit column carries a real off-diagonal multiplicative correlation, zero in the mean for a random column of the same density, and the Mobius and Liouville signs cancel it no better than an unstructured sign vector on the same support does, `|Sigma_mu|` sitting at `0.0913` to `0.7178` of the random-sign root mean square against `0.0359` to `1.5048` for the support-matched controls over sixteen boxes, with the split against those controls `3, 9, 4` at chi-square `0.375` against the uniform-rank null; a sign vector built by greedy flips against a known column drives the same Cauchy-Schwarz bound to `0.0265` of its diagonal floor, so the census measures the arithmetic of the coefficients and not a limit of the method (`lab/rs/rho-decoupling`, `menergy signed` and `menergy signed engineered`).
- **The coefficient the method is given is not the coefficient it would need. Verified.** At seven of the eight swept boxes with both sides above `x^(2/5)`, the boxes a Vaughan decomposition actually produces, the coefficient sequence it hands the bilinear sum takes values in `{-1, 0, 1}` and needs no normalisation, and its full quadratic form sits between `0.69` and `1.21` of its own diagonal, where a sign vector engineered against the column reads `0.13` to `0.21` on the same boxes, so the sequence the method is given and the sequence the method would need are different objects (`lab/rs/rho-decoupling`, `menergy signed vaughan`).
- **The large-values refinement is the moment route itself. Verified.** The large-values refinement of the moment route is costed out and is the `l^2` route itself, exponent `(1 + alpha)/2` at every threshold (`lab/rs/rho-decoupling`, `riesz large values chain`, 66 cells over six families); the large frequencies are adjacent grid points (`407` in `331` runs at `{0,1}` base 3, `level 12`, `eta = eta_4`), so the grid offers no spacing gain, and at the dense families the bilinear sum at `a_m = b_l = 1` equals the box representation count, `0.38 fill^level` there, so no bound uniform over bounded coefficients holds at those families.

## The exact transfer between designs

The census columns are tied together by one carry-free mechanism. **Proved:**

- **Scaling.** If every digit of `F` is `a` times a digit of `F'`, so `F = aF'` inside `{0..base-1}`, then `m -> am` maps `S_F'` bijectively onto `S_F` preserving digit length: `am = sum (a d_j) base^j` and each `a d_j <= base - 1`, so no carry occurs and the digit string scales digitwise. Hence `A_F(base^level)` equals the string count of `F'` at the same depth, and `M_F(base^level) = sum of mu(am)` over `m in S_F'` with at most `level` digits.
- **Vanishing.** If `a` has a square factor then `mu(am) = 0` for every `m`, so `M_F` is identically zero: at `base 5`, `F = {0,4} = 4 x {0,1}`, the meter reads 0 at all 21 levels. A census that reads cancellation without factoring out the digit gcd reads this as infinite cancellation; the digit gcd must be squarefree before `theta` means anything.
- **Prime twist.** If `a = p` is prime then `mu(pm)` is `-mu(m)` on `p`-free `m` and `0` otherwise, so `M_(pF')(base^level) = -sum of mu(m)` over the `m in S_F'` not divisible by `p`. At `base 3`, `F' = {0,1}`: an element `m = sum of 3^j` is odd exactly when its count of 1-digits is odd, and reading the digit string as a binary index that [parity](https://en.wikipedia.org/wiki/Parity_(mathematics)) is the [Thue-Morse](https://en.wikipedia.org/wiki/Thue%E2%80%93Morse_sequence) sign, so the `{0,2}` column is the Thue-Morse-twisted `{0,1}` column.
- **Base-4 anti-symmetry.** At `base 4`, `M_{0,2}(4^level) = -M_{0,1}(4^level)` exactly: since `4 | base`, an element of `S_{0,1}` is `0` or `1 mod 4` by its unit digit, so every even element is divisible by 4 and carries `mu = 0`, and the odd-part twist above is minus the whole meter. Stronger, `M_{0,2}(x) = -M_{0,1}(x/2)` at every real `x`, and since `S_{0,1}` has no element strictly between `(4^level - 1)/3` and `4^level` the running maxima agree level by level as well. The census confirms both at all 22 levels, e.g. meters `-110/110` at `level 15`, `-342/342` at `level 17`, `34/-34` at `level 22`, and `Mmax = 1553` for both at `level 22`.

**Verified:** the generator recomputes all eight scaled census families (`{0,2}` at `base 3`; `{0,2}, {0,3}` at `base 4`; `{0,2}, {0,3}, {0,4}, {2,4}, {0,2,4}` at `base 5`) from their primitive families through `mu(am)` and asserts equality at every level. The mechanism needs a common digit factor, so it partitions the census into primitive columns and their twists and says nothing across primitive columns.

## The census

Every `M_F(base^level)` below is an exact integer: restricted families enumerated in ascending order with `mu` from deterministic factorization (trial division, Miller-Rabin on the twelve witnesses `2..37`, Pollard rho), controls by a linear Mobius sieve; one family (`base 3`, `F = {1,2}`, `level 16`) is computed by both methods and asserted equal at every level. The base-10 control reproduces [A084237](https://oeis.org/A084237) (`-1, 1, 2, -23, -48, 212, 1037, 1928` at `10^1..10^8`). Every table below is extracted by script from the generator's printed rows, never assembled by hand. **Verified.**

The three base-3 columns, checkpoint meter and running maximum `Mmax = max of |M_F(x)|` over `x <= 3^level` per row:

| `level` | `M_{0,1}` | max | `M_{0,2}` | max | `M_{1,2}` | max |
|---|---|---|---|---|---|---|
| 4 | -2 | 3 | 2 | 3 | -8 | 8 |
| 6 | 2 | 5 | 0 | 3 | -8 | 11 |
| 8 | 2 | 8 | 3 | 7 | -31 | 33 |
| 10 | 5 | 13 | 0 | 11 | -14 | 38 |
| 12 | 56 | 61 | -37 | 40 | -35 | 88 |
| 14 | 11 | 105 | -10 | 67 | -205 | 230 |
| 16 | 149 | 173 | -124 | 152 | 4 | 281 |
| 18 | -30 | 312 | 67 | 249 | -1461 | 1582 |
| 20 | 496 | 539 | -382 | 485 | -3175 | 3255 |
| 21 | 533 | 866 | -194 | 617 | -2005 | 3855 |
| 22 | 1009 | 1089 | -1205 | 1324 | -690 | 3855 |
| 23 | 1824 | 2848 | -2242 | 2942 | -3214 | 3855 |
| 24 | -1886 | 3296 | -133 | 3843 | -3248 | 4113 |

The final checkpoint of every family, with `thetamax = log(Mmax)/log A` and its drift (max minus min) over the last five levels:

| `base` | `F` | `level` | `A_F(base^level)` | `M_F(base^level)` | `Mmax` | `thetamax` | drift |
|---|---|---|---|---|---|---|---|
| 3 | 01 | 24 | 16777216 | -1886 | 3296 | 0.4869 | 0.0452 |
| 3 | 02 | 24 | 16777215 | -133 | 3843 | 0.4962 | 0.0596 |
| 3 | 12 | 24 | 33554430 | -3248 | 4113 | 0.4802 | 0.0754 |
| 4 | 01 | 22 | 4194304 | 34 | 1553 | 0.4819 | 0.0391 |
| 4 | 02 | 22 | 4194303 | -34 | 1553 | 0.4819 | 0.0391 |
| 4 | 03 | 22 | 4194303 | -541 | 1180 | 0.4638 | 0.0488 |
| 4 | 12 | 22 | 8388606 | -855 | 3965 | 0.5197 | 0.1056 |
| 4 | 13 | 22 | 8388606 | -712 | 2631 | 0.4940 | 0.0727 |
| 4 | 23 | 22 | 8388606 | -3255 | 3258 | 0.5074 | 0.0157 |
| 4 | 012 | 14 | 4782969 | -503 | 1057 | 0.4527 | 0.0475 |
| 4 | 013 | 14 | 4782969 | 2313 | 2899 | 0.5183 | 0.0487 |
| 4 | 023 | 14 | 4782968 | -753 | 1166 | 0.4591 | 0.0862 |
| 4 | 123 | 14 | 7174452 | -592 | 1644 | 0.4691 | 0.0729 |
| 5 | 01 | 21 | 2097152 | 153 | 849 | 0.4633 | 0.0485 |
| 5 | 02 | 21 | 2097151 | 250 | 889 | 0.4665 | 0.0268 |
| 5 | 03 | 21 | 2097151 | -116 | 700 | 0.4501 | 0.0311 |
| 5 | 04 | 21 | 2097151 | 0 | 0 | - | - |
| 5 | 12 | 21 | 4194302 | -128 | 1643 | 0.4856 | 0.0732 |
| 5 | 13 | 21 | 4194302 | -2875 | 3533 | 0.5358 | 0.0456 |
| 5 | 14 | 21 | 4194302 | -1511 | 2750 | 0.5193 | 0.0640 |
| 5 | 23 | 21 | 4194302 | 405 | 914 | 0.4471 | 0.0581 |
| 5 | 24 | 21 | 4194302 | 1065 | 2287 | 0.5072 | 0.0540 |
| 5 | 34 | 21 | 4194302 | -2137 | 2538 | 0.5141 | 0.0401 |
| 5 | 012 | 13 | 1594323 | -1016 | 1416 | 0.5080 | 0.0846 |
| 5 | 013 | 13 | 1594323 | -137 | 768 | 0.4652 | 0.0528 |
| 5 | 014 | 13 | 1594323 | 213 | 1005 | 0.4840 | 0.0821 |
| 5 | 023 | 13 | 1594322 | 759 | 858 | 0.4729 | 0.0455 |
| 5 | 024 | 13 | 1594322 | 686 | 959 | 0.4807 | 0.0240 |
| 5 | 034 | 13 | 1594322 | 501 | 1000 | 0.4837 | 0.0847 |
| 5 | 123 | 13 | 2391483 | 88 | 816 | 0.4565 | 0.0489 |
| 5 | 124 | 13 | 2391483 | 1036 | 1613 | 0.5029 | 0.0604 |
| 5 | 134 | 13 | 2391483 | -725 | 981 | 0.4690 | 0.0519 |
| 5 | 234 | 13 | 2391483 | -1926 | 2021 | 0.5182 | 0.1056 |
| 5 | 0123 | 11 | 4194304 | -474 | 1725 | 0.4887 | 0.0732 |
| 5 | 0124 | 11 | 4194304 | 426 | 1494 | 0.4793 | 0.0673 |
| 5 | 0134 | 11 | 4194304 | -644 | 1633 | 0.4852 | 0.0222 |
| 5 | 0234 | 11 | 4194303 | -362 | 2179 | 0.5041 | 0.0794 |
| 5 | 1234 | 11 | 5592404 | 145 | 1101 | 0.4508 | 0.0996 |

Base 10 with one digit excluded, the Kempner designs (`fill = 9`, the sets behind the convergent harmonic series of [Kempner 1914](https://doi.org/10.2307/2972074), revisited at `s = 1` in [Allouche, Hu and Morin 2024](https://arxiv.org/abs/2403.05678)), at `x = 10^8`:

| excluded | `A_F(10^8)` | `M_F(10^8)` | `Mmax` | `thetamax` |
|---|---|---|---|---|
| 0 | 48427560 | 6410 | 8177 | 0.5091 |
| 1 | 43046720 | 4108 | 6069 | 0.4956 |
| 2 | 43046721 | -183 | 3357 | 0.4619 |
| 3 | 43046721 | 455 | 3512 | 0.4644 |
| 4 | 43046721 | 56 | 4957 | 0.4841 |
| 5 | 43046721 | -7614 | 10601 | 0.5273 |
| 6 | 43046721 | -693 | 2564 | 0.4465 |
| 7 | 43046721 | -1411 | 6494 | 0.4994 |
| 8 | 43046721 | 2131 | 4495 | 0.4785 |
| 9 | 43046721 | 2181 | 5234 | 0.4871 |

The full-set controls at comparable depth: `M(3^17) = -1423` with `Mmax = 4610` (`thetamax` 0.4517), `M(4^13) = 329` with `2845` (0.4413), `M(5^11) = 617` with `2573` (0.4436), `M(10^8) = 1928` with `3448` (0.4422). The Mertens function itself - limiting exponent exactly `1/2` if and only if RH, and at least `1/2` unconditionally - reads `0.4413..0.4517` at these depths, which calibrates every reading above: at census mass even the classical meter sits a few hundredths under `1/2`.

The distribution of the apparent exponent across designs at fixed base, sorted by the generator: at `base 3` the three columns read `0.4802, 0.4869, 0.4962`; at `base 4` the ten run `0.4527` to `0.5197`; at `base 5` the twenty-four with nonzero meter run `0.4471` to `0.5358`; at base 10 the ten Kempner columns run `0.4465` to `0.5273`. All 47 readings sit within `0.054` of `1/2`, against cut readings (`theta` at the checkpoint alone) that scatter over `0.22..0.53` for the same data - the single-cut estimator is noise, the running maximum is the meter.
