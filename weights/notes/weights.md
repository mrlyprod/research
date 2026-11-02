---
title: Weighted designs
lead: Weighted designs: what weights move on the mass side, what they never move on the length side, and the pressure that closes both.
figure: research-weights
slug: weights
---

A design is a set of cells; a weighted design is the same cells carrying a probability vector, and that is a refinement equation on the base grid. Weights move the mass side of the object and only the mass side: the contraction ratios stay `1/base`, so every length-indexed observable keeps its `log base` ripple at every weight, while the mass-stopping count - the first mass-indexed observable this tree has - is log-periodic or smooth according to the arithmetic of the `log w_f` alone. The multifractal pressure of the same object closes in one line, `tau(s) = log_base sum_f w_f^s`.

`lab/py/weighted-designs` is the one generator behind every number below: exact rationals in, safe-rounded floats out, every number asserted before it prints. The lattice and nonlattice vocabulary is [the dimensions page](../../dimensions/notes/dimensions.md)'s and is used here rather than restated.

## The object

- A weighted design is the refinement equation `phi(x) = sum_(f in F) c_f phi(base x - f)`, `x` in `R^dim`, `F` the design's filled digits, inside `{0,...,base-1}^dim`, `c_f = base^dim w_f`.
- The measure form is the same object: `mu = sum_f w_f mu . S_f^(-1)` with `S_f(x) = (x + f)/base` and `sum_f w_f = 1`. `phi` is the density that `mu` does not have.
- Equal weights `w_f = 1/|F|` is the 0/1 design: the support is the design itself, the geometry is untouched, and it stays lattice ([dimensions](../../dimensions/notes/dimensions.md)).
- [Daubechies 1988](https://doi.org/10.1002/cpa.3160410705) is the same equation at base 2, dim 1, digits `{0,...,N}`: her scaling functions are weighted digit sets that **overlap**, `N >= base`, which no design has, and her remark after (4.29) is the reason regularity is bought only by widening the mask past one residue box.
- Values are exact rationals over a common denominator, `w_f = n_f/d` with `sum_f n_f = d`, so every mass at level `level` is an integer over `d^level` and floats appear only at display. The generator's object is base 3, dim 2, cells `(0,0) (2,0) (0,2)`, weights `3/8, 3/8, 1/4`, integer masses `[3, 3, 2]` over `8`, level 10.
- A symmetry of the cells acts on a weighted design by permuting the `(f, w_f)` pairs, so two weighted designs are one object exactly when one symmetry matches cell and weight together: the unweighted equivalence of [core](../../core/notes/core.md) with the weight vector carried along.

## What weights never move

The contraction ratios are `1/base` at every weight, so the geometry is the unweighted design's own and every length-indexed observable keeps its `log base` ripple forever.

**The corner lemma (Proved).** At a filled corner digit of weight `w_0`, a ball of radius `r < min_(f != 0) |f|/base` about that digit's fixed point meets its own sub-block and no other, so `mu(B(r)) = w_0 mu(B(base r))`, and `ln mu(B(r)) - alpha ln r` with `alpha = -log_base w_0` is exactly `log base`-periodic at every weight. The same one-map argument at any filled digit `d`, whose map `x -> (x + d)/base` fixes `p_d = d/(base-1)`, gives `mu(B(p_d, r)) = w_d mu(B(p_d, base r))` for `r` below the distance from `p_d` to the other filled cells; at equal weights `w_d = 1/fill` this is the identity the [spin](../../spin/notes/spin.md) and [crop](../../crop/notes/crop.md) pages use. The range is not decoration: at base 3, dim 1, `F = {0,1,2}` and `w = (1/2, 1/4, 1/4)` the radius `r = 1/2` is past `min_(f != 0) |f|/base = 1/3`, and `X = mu(B(1/2))` satisfies `X = 1/2 + X/4`, so `X = 2/3` against `w_0 mu(B(3/2)) = 1/2` and the identity fails. On the generator's object `min_(f != 0) |f|/base = 2/3`.

**The ripple survives every weight (Verified).** `M(r)`, the mass within radius `r` of the corner fixed point, is read as an exact integer shell histogram, detrended by `alpha` and folded into 24 bins of `log_3 r` over whole periods of `R` in `3^4..3^8`. Every set keeps the ripple, and the cross-level identity `M_level(r) = n_0 M_(level-1)(r)` for `r < 2 base^(level-1)` - the identity above, checked - is exact integer equality for all four. `M(r)` is the spin mass of [spin](../../spin/notes/spin.md) read at weight, and being a length observable it is the wrong detector for the arithmetic class of the `log w_f`: it cannot smooth.

| weights | `alpha` | swing | drift bar | whole-period bar | gap against equal weights |
|---|---|---|---|---|---|
| `1/3,1/3,1/3` | 1.000000000 | 0.652441 | 0.022937 | 0.005857 | 0.000000 |
| `2/5,2/5,1/5` | 0.834043767 | 0.536956 | 0.019130 | 0.004779 | 0.066973 |
| `3/8,3/8,1/4` | 0.892789261 | 0.580858 | 0.020478 | 0.005133 | 0.044604 |
| `1/2,1/4,1/4` | 0.630929754 | 0.427372 | 0.014471 | 0.003019 | 0.153880 |

Equal weights give every cell mass 1, so `M(r)` is the cell count and its gap against the 0/1 ripple is `0.000000` by identity; the other three differ from it by `0.066973`, `0.044604` and `0.153880` on drift bars of `0.014471` to `0.020478`, so weights reshape the ripple and never remove it. The `log 3` periodicity residual falls by decade of `R` - `0.0208, 0.0094, 0.0047, 0.0021` at equal weights, worst `0.0299` down to `0.0039` - which is discretisation, not drift.

**Non-Rajchman at every weight (Proved).** `hat mu(base t) = P_w(t) hat mu(t)` with `P_w(t) = sum_f w_f e(-f . t)` - the digit transform `hat F` of [mobius](../../mobius/notes/mobius.md), weighted, normalised and read at `-t`, written `P_w` on this page - and `P_w = 1` at every integer `t` because `sum_f w_f = 1`, so `hat mu(base^m t) = hat mu(t)` on `Z^dim` and the Fourier dimension is zero whatever the weights. Weights do not rescue the Mobius door ([mobius](../../mobius/notes/mobius.md)).

## What weights move

**The Dirichlet root is not a weight observable (Refuted).** That `delta`, the root of `sum_f w_f^s = 1`, reads the weighting is false: `s = 1` solves that equation for every probability vector, and the generator prints `delta = 1.0000000000` at all four of them. What weights own on the mass side is therefore the arithmetic class of the `log w_f`, never the root.

**The mass-stopping count (Proved).** `N(t) = #{words of mass >= t}`, written `N(a)` at `t = e^(-a)`, is the tree's first mass-indexed observable: it is log-periodic exactly when the group generated by the `log w_f` is cyclic, and smooth otherwise. The mechanism is [Lalley 1989](https://doi.org/10.1007/BF02392732), whose dichotomy is a property of the renewal function alone, fed here the mass variable `f = -log w`: Proposition 2.1 gives the unique `delta > 0`, Theorem 1 the nonlattice `N(a) ~ C e^(a delta)`, Theorem 2 the integer-valued `N(a) ~ C e^([a] delta)`. The class is decided exactly and never fitted: each `w_f` becomes its vector of prime exponents, and the group is cyclic exactly when that integer matrix has rank 1. This holds for positive rational weights only: there a weight carries a prime-exponent vector, and the `Q`-linear independence of the logs of the primes gives the exponent lattice and the group generated by the `log w_f` the same rank.

Four of the seven sets below sum to one and so are weighted designs proper - `1/3,1/3,1/3` is the equal-weight case, the 0/1 design itself, and `1/2,1/4,1/4`, `2/5,2/5,1/5` and `3/8,3/8,1/4` are the three proper weightings; the other three are renewal weight sets, admissible because the dichotomy reads `-log w` and never asks it to normalise. (Verified, the rank and the oscillation alike.)

| weights | primes | rank | class | `delta` | oscillation of `ln N(a) - delta a` over a window of `ln 3` at `a = 20, 40, 60, 80, 100` |
|---|---|---|---|---|---|
| `1/2,1/2,1/2` | 2 | 1 | lattice, span `ln 2` | 1.5849625007 | `1.09812  1.09852  1.09812  1.09812  1.09803` flat |
| `1/2,1/2,1/4` | 2 | 1 | lattice, span `ln 2` | 1.2715533032 | `0.88098  0.88130  0.88098  0.88098  0.88091` flat |
| `1/3,1/3,1/3` | 3 | 1 | lattice, span `ln 3` | 1.0000000000 | `1.09825  1.09825  1.09825  1.09825  1.09825` flat |
| `1/2,1/4,1/4` | 2 | 1 | lattice, span `ln 2` | 1.0000000000 | `0.69284  0.69309  0.69284  0.69284  0.69278` flat |
| `1/2,1/2,1/3` | 2, 3 | 2 | nonlattice | 1.3646005647 | `0.20259  0.14175  0.11707  0.10183  0.09030` decaying |
| `3/8,3/8,1/4` | 2, 3 | 2 | nonlattice | 1.0000000000 | `0.23133  0.16019  0.13671  0.11972  0.10734` decaying |
| `2/5,2/5,1/5` | 2, 5 | 2 | nonlattice | 1.0000000000 | `0.23297  0.19792  0.18623  0.17623  0.16776` decaying |

**Normalisation changes the class.** `1/2, 1/2, 1/4` has rank 1 and is lattice at span `ln 2`; the same weights normalised, `2/5, 2/5, 1/5`, have rank 2 and are nonlattice. Fix the probability vector before the word lattice means anything. The length side does not see the difference, since unnormalised weights differ from their normalisation by one constant per level and the detrending absorbs it; the mass side does, which is the whole trap.

**The controls, pre-registered (Verified).** Equal weights must reproduce the 0/1 design exactly, log-commensurable weights must keep every ripple, and only an irrational log ratio may move a mass observable. `1/3,1/3,1/3` lattice, length kept and mass flat; `1/2,1/4,1/4` lattice, length kept and mass flat; `3/8,3/8,1/4` nonlattice, length kept and mass moved. The length observable moves under none of the three and the mass observable under exactly one, which is the sign the controls predicted; a mass observable moving under all three would have failed here.

## The pressure

**Verified.** At equal contraction `1/base` under the open set condition the pressure equation `sum_f w_f^s (1/base)^tau(s) = 1` closes, so `tau(s) = log_base sum_f w_f^s` and `f(alpha) = inf_s (alpha s + tau(s))` is explicit, with `tau(0) = log_base |F| = 1.000000000` and `tau(1) = 0.000000000`. This is the multifractal formalism for self-similar measures under the open set condition, [Cawley and Mauldin 1992](https://doi.org/10.1016/0001-8708(92)90064-R), read in [Olsen's restatement](https://arxiv.org/abs/1307.5223), Section 3, the original being behind a publisher wall and not opened here.

| `s` | `sum_f w_f^s` | `tau(s)` | `alpha(s)` | `f(alpha(s))` |
|---|---|---|---|---|
| -2 | `272/9` | 3.102620937 | 1.088179391 | 0.926262155 |
| -1 | `28/3` | 2.033103256 | 1.050962223 | 0.982141033 |
| 0 | `3` | 1.000000000 | 1.015812676 | 1.000000000 |
| 1 | `1` | 0.000000000 | 0.985056822 | 0.985056822 |
| 2 | `11/32` | -0.971990429 | 0.959892942 | 0.947795455 |
| 3 | `31/256` | -1.921688171 | 0.940411228 | 0.899545513 |

The box partition at level `level` carries the moments exactly, `sum_i mu_i^s = (sum_f w_f^s)^level` as a rational identity at `s = -2..3` over 7, 9 and 11 distinct box masses in 729, 6561 and 59049 boxes at levels 6, 8, 10. The coarse-grained band `f_level(alpha) = log_base N(alpha)/level` sits under the transform at every level and every achievable `alpha`, in one line from `N_i mu_i^s <= base^(level tau(s))`. At the band's middle `alpha = 1.077324384` the readings are `0.769937048, 0.798858255, 0.818775202` against `f = 0.946394630`, deficits `0.176458, 0.147536, 0.127619` at levels 6, 8, 10; the gap is the Stirling volume term of the multinomial count, matched to `0.000077, 0.000058, 0.000046`, and it falls with the level. Both endpoints are exact at every level, `f_level(0.892789261) = f = 0.630929754 = log_3 2` and `f_level(1.261859507) = f = 0`, and the band's top at level 10 is `0.877427106` rising to `tau(0) = 1`. The weights demo puts the weights on sliders: the support never moves while the mass does, and `tau(s)` with its Legendre transform redraws as a closed form in the weights alone.

## The ladder

**Rung 0 is closed form (Verified).** `F` inside `{0,...,base-1}^dim` puts `supp phi` inside the unit cell, so exactly one integer translate meets it, the [Daubechies and Lagarias 1992](https://doi.org/10.1137/0523084) matrices are `1 x 1` with `T_f = [base^dim w_f]`, and the joint spectral radius - norm-independent by [Rota and Strang 1960](https://doi.org/10.1016/S1385-7258(60)50046-1) - is `base^dim max_f w_f`. Hence `alpha_Holder(phi) = -dim - log_base(max_f w_f) = alpha_min(mu) - dim`, negative for every proper design: the tree's objects are measures, never functions, and digit independence is exactly what empties the regularity question. The same word solver, run on rung 0, collapses onto the closed form.

| weights | `T_f` | `JSR` | solver | `alpha_min(mu)` | `alpha_max(mu)` | `alpha_Holder(phi)` |
|---|---|---|---|---|---|---|
| `1/3,1/3,1/3` | `[3, 3, 3]` | `3` | `[3.0000000, 3.0000000]` | 1.000000000 | 1.000000000 | -1.000000000 |
| `2/5,2/5,1/5` | `[18/5, 18/5, 9/5]` | `18/5` | `[3.6000000, 3.6000000]` | 0.834043767 | 1.464973521 | -1.165956233 |
| `3/8,3/8,1/4` | `[27/8, 27/8, 9/4]` | `27/8` | `[3.3750000, 3.3750000]` | 0.892789261 | 1.261859507 | -1.107210739 |
| `1/2,1/4,1/4` | `[9/2, 9/4, 9/4]` | `9/2` | `[4.5000000, 4.5000000]` | 0.630929754 | 1.261859507 | -1.369070246 |

**Rung 1 is the first overlap (Verified).** At base 2, mask with `sum c_even = sum c_odd = 1`, the same solver runs `T_0 = (c_(2i-j-1))` and `T_1 = (c_(2i-j))`, `1 <= i, j <= N`, on `{v : sum_i v_i = 0}`, and `alpha = -log_2 JSR`. The hat mask `(1/2, 1, 1/2)` returns `JSR` in `[0.5000000, 0.5000000]` and `alpha` in `[1.0000000, 1.0000000]`, exactly. D4 returns Gripenberg lower `0.6830127` and norm upper `0.7105812` at length 14, so `alpha` lies in `[0.4929285, 0.5500157]`, which contains the closed form `2 - log_2(1 + sqrt 3) = 0.5500157`.

The caveat travels with the bracket. D4's lower end is attained at the one-letter word `T_0`, whose [spectral radius](https://en.wikipedia.org/wiki/Spectral_radius) is `(1 + sqrt 3)/4`, so `alpha`'s **upper** end is that closed form by construction and its agreement is an identity rather than a test. The hat mask, where the bracket collapses to a single point, is the real validation of the solver, and D4's live content is the norm upper end `0.7105812`. **Refuted:** that a norm upper bound may print truncated. The bracket is a certificate - all matrix arithmetic in `Q(sqrt m)`, spectral radii and norms enclosed through trace and determinant, Gripenberg words to length 8 truncating down and norm upper bounds to length 14 rounding up - so the safe print of `alpha`'s upper end at four digits is `0.5501`, and `0.5500` would sit strictly below the closed form it is supposed to contain.

## Where the numbers live

`lab/py/weighted-designs` is the one pass behind this page: one exact generator taking `(base, dim, F, w)` and printing the arithmetic class, `delta`, the mass-stopping count, the `M(r)` ripple, the local-dimension range, the multifractal band and the JSR bracket at both rungs, over the pre-registered controls. It runs in about four seconds, prints only, writes nothing, holds one level-10 cell list at a time, asserts every number before it prints and exits nonzero on any failure; the regression is the object of the first section. The lattice geometry these weights leave alone is [dimensions](../../dimensions/notes/dimensions.md), the Fourier door they do not open is [mobius](../../mobius/notes/mobius.md), and the unweighted object underneath is [core](../../core/notes/core.md).
