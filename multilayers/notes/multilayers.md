---
title: The fractal mirror
lead: Light through a stack of glass laid on a design: the level recursion that makes `b^k` layers cost `(b - 1)k` multiplications, the Born drift `m(P_D)` of the digit polynomial, the fixed-contrast step law that reads the two-variable lift `Q_D` instead, the cyclotomic designs that still go dark, the deep drift that misses `m(Q_D)` at weak contrast, and the conjectured split between critical stacks, whose passband shrinks as `1/sqrt(k)`, and the rest, which darken geometrically.
figure: research-multilayers
slug: multilayers
---

A design at base `b` with digit set `D = {d_0 < ... < d_(m-1)}` is the set of integers whose `k` base-`b` digits all lie in `D`, the one-dimensional missing-digit set and at `b = 3`, `D = {0,2}` the [Cantor set](https://en.wikipedia.org/wiki/Cantor_set) at level `k`. Lay glass on it: cell `p` of a row of `b^k` cells is layer A (index `nA`) when `p` is in the design and host B (index `nB = 1.45`) otherwise, every cell of the same optical thickness, so light at normal incidence picks up the same phase `delta` crossing any cell. The stack sits in B on both sides; `r_k` and `t_k` are its reflection and transmission amplitudes, `T_k = abs(t_k)^2`, and the stack is open at `delta` when `abs(r_k/t_k) < 1`, that is `T_k > 1/2`. The passband fraction `f_k` is the share of `delta` in `(0, pi)` where the level-`k` stack is open.

Two polynomials carry the design. The digit polynomial `P_D(z) = sum_(d in D) z^d`, and the block-spacer lift `Q_D(x, y) = sum_i x^i y^(d_i - i)`, where `x` stands for one level-`k` block and `y` for one empty spacer of the same length; `Q_D(z, z) = P_D(z)`. Their [Mahler measures](https://en.wikipedia.org/wiki/Mahler_measure) are `m(P) = integral ln abs(P)` over the unit circle and the torus. The generator for every number below is `lab/py/design-multilayer`; its verbs are named where they print.

## The level recursion

**The level-`(k+1)` characteristic matrix is `M_(k+1)(delta) = W_D(M_k(delta), M_B(b^k delta))`, the ordered product over `d = 0 .. b-1` of `M_k` when `d` is in `D` and of the plain host slab `M_B(b^k delta)` otherwise. Proved.** The level-`(k+1)` word is `b` blocks of length `b^k` whose leading digit is `d`: the block is the level-`k` word when `d` is in `D` and `b^k` host cells otherwise, the [transfer matrix](https://en.wikipedia.org/wiki/Transfer-matrix_method_(statistical_mechanics)) of a concatenation is the product of the parts, and `b^k` host cells are one host slab of phase `b^k delta`. It is the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product) read as a substitution, and it costs `b - 1` matrix multiplications a level, so `7^12` layers cost `72`. Exact self-similar algorithms of this kind go back to [Jaggard and Sun 1990](../refs.md) for Cantor-bar multilayers.

**The recursion agrees with the layer-by-layer product to `1e-12` relative. Verified** (`check`): the brute-force product runs in 40-digit arithmetic over `3^6`, `4^5`, `5^4`, `7^4` and `9^3` cells at `nA = 1.6, 2.3, 3.5`, eight frequencies each, and the largest relative gap is `8.4e-13`.

## The Born drift

**At zero contrast the stack's reflection has the design's generating function as its derivative: `abs(d(r_k/t_k)/d nA)` at `nA = nB` equals `abs(d(r_0/t_0)/d nA)` times `abs(prod_(j<k) P_D(e^(2i b^j delta)))`, so the log of its modulus, averaged over `delta` uniform in `(0, pi)`, grows by exactly `m(P_D)` a level. Proved.** At `nA = nB` the stack is homogeneous and `r = 0`; the first derivative of the product in `nA` is a sum over A cells of one perturbed factor, each A cell at optical position `p delta` contributing the single-cell derivative times `e^(2i p delta)` up to one common phase, and the positions are the digit sums, whose generating function is `prod_(j<k) P_D(z^(b^j))`. Each factor's mean log modulus is `m(P_D)` because `b^j delta` is uniform whenever `delta` is. `check` confirms the identity at `nA = nB(1 + 1e-8)` to `5e-7` in `abs(gap)/max(1, abs(prod))` on three designs through level 4.

This is a statement about the zero-contrast derivative and nothing more; the true stack follows it only while `nA - nB` times the product stays small. **At `nA = 1.455` the mean per-level drift of `ln abs(r_k/t_k)` over six to eight levels is within `0.001` of `m(P_D)`. Verified** (`born`): `0.3832` against `m(P_D) = 0.38225` for `{0,1,3}` base 4, `0.4425` against `0.44214` for `{0,1,2,4}` base 5, `0.2818` against `0.28120` for `{0,1,5}` base 7, and within `0.004` of `0` for the three designs with cyclotomic `P_D`, standard errors below `0.002`. The agreement is to leading order, not to noise: `0.3832` sits about three standard errors off. The Born picture is the diffraction picture: for binary constant-length substitutions the Mahler measure of a Borwein polynomial is a maximal Lyapunov exponent of the Fourier matrix cocycle ([Baake, Coons and Manibo 2020](../refs.md)).

## The block-spacer lift

At fixed contrast the Born picture breaks within a few levels, and the reason is phase. A level-`k` block that is nearly transparent acts on light as a phase `u = t_k/abs(t_k)`, and at zero contrast `u` equals the spacer's phase `v = e^(i b^k delta)`; at finite contrast the two part. To first order in `r_k`, the next level's reflection is the sum of the `m` block echoes, the `i`-th reached through `i` blocks and `d_i - i` spacers, so `r_(k+1) = r_k Q_D(u^2, v^2)` up to a phase and a relative error of order `abs(r_k)`. The digit polynomial is the diagonal `u = v` of that law.

**The step law `ln abs(r_(k+1)/t_(k+1)) - ln abs(r_k/t_k) = ln abs(Q_D(u^2, v^2)) + O(abs(r_k))`. Verified** (`deep`): at every level where `abs(r_k/t_k) < 1e-3` on at least 50 frequencies, the median gap between the two sides is below `3e-6`, on the nine `deep` designs (the seven of the table, `{0,1,3}` base 4 and `{0,1,2,4}` base 5) at `nA = 1.6, 2.3, 3.5`. That covers 25 of the 27 design-contrast cells; `{0,1,3}` base 4 and `{0,1,2,4}` base 5 at `nA = 3.5` never reach 50 such frequencies. The first-order expansion above is the reason; it is not written out here as a proof.

For a digit set in arithmetic progression, `Q_D` depends on one monomial and the level step is exact at any contrast: one period of the progression is a block and its spacer, and the Abeles formula for a periodic stack gives `abs(r_(k+1)/t_(k+1)) = abs(U_(m-1)(chi_k)) abs(r_k/t_k)` with `chi_k` the half-trace of one period and `U` the Chebyshev polynomial of the second kind. It is the textbook periodic-stack formula; Chebyshev polynomials of the second kind also carry the exact transmission of generalised Cantor-like potentials in [Ogawana and Sakaguchi 2018](../refs.md).

## Cyclotomic is not enough

The Born drift vanishes exactly when `P_D` is a product of cyclotomic polynomials (Kronecker's theorem, the one-variable case of Boyd's), and the natural guess is that such a stack is critical: no drift, passband lost only slowly. The step law says the guess reads the wrong polynomial.

**"`P_D` cyclotomic implies a critical stack." Refuted** (`refute`) by `{0,2,3,4,6}` base 7, `{0,1,4,7,8}` base 9 and `{0,2,3,4,5,7}` base 8. For each, `P_D` factors into cyclotomic polynomials over the integers, while `Q_D` is not a monomial times cyclotomic polynomials in monomials, the exact test run by factoring in two variables, so `m(Q_D) > 0` by Boyd's theorem ([Boyd 1981](../refs.md)): `m(Q_D) = 0.2513, 0.2513, 0.3181`. At `nA = 2.3` their `f_k sqrt(k)` falls without pause from level 2 to level 12, `0.878` to `0.548`, `0.883` to `0.538` and `0.888` to `0.406`, while the controls `{0,2,4,6}` base 7 and `{0,1,3,4}` base 7, whose `Q_D` is Boyd-cyclotomic, hold at `1.097` and `0.901`. Their deep drift, below, is at least `0.26` a level at all three contrasts.

## The deep drift

If the pair `(u, v)` were spread evenly over the torus on the passband, the step law would make the drift `m(Q_D)`, the two-variable Mahler measure. It is close, and at weak contrast it is not equal.

The table is `deep`: the deep drift is the mean step over frequencies with `abs(r_k/t_k) < 0.02`, averaged over the last four measured levels, and `se` is the mean of their four standard errors, which bounds the standard error of that average; each cell lists `nA = 1.6, 2.3, 3.5`, on `2e5` frequencies.

| design | `m(P_D)` | `m(Q_D)` | deep drift | se | torus mean |
|---|---|---|---|---|---|
| `{0,2,3,4,6}` base 7 | `0` | `0.2513` | `0.272, 0.277, 0.288` | `0.005, 0.008, 0.012` | `0.249, 0.227, 0.183` |
| `{0,1,4,7,8}` base 9 | `0` | `0.2513` | `0.264, 0.275, 0.283` | `0.005, 0.008, 0.012` | `0.251, 0.225, 0.163` |
| `{0,2,3,4,5,7}` base 8 | `0` | `0.3181` | `0.350, 0.358, 0.355` | `0.005, 0.009, 0.011` | `0.318, 0.292, 0.232` |
| `{0,1,5}` base 7 | `0.2812` | `0.3231` | `0.308, 0.314, 0.316` | `0.005, 0.012, 0.018` | `0.321, 0.323, 0.325` |
| `{0,2}` base 3 | `0` | `0` | `-0.001, 0.000, -0.001` | `0.004, 0.005, 0.007` | `-0.005, -0.003, 0.002` |
| `{0,1,3,4}` base 7 | `0` | `0` | `0.005, 0.002, 0.009` | `0.005, 0.007, 0.009` | `-0.001, 0.001, 0.004` |
| `{0,2,4,6}` base 7 | `0` | `0` | `-0.001, 0.000, -0.002` | `0.004, 0.006, 0.008` | `0.001, 0.003, 0.001` |

**"At `nA = 1.6` the drift of `ln abs(r/t)` deep in the passband tends to `m(Q_D)`." Refuted** (`deep`): the deep drift is `0.350` against `m(Q_D) = 0.3181` for `{0,2,3,4,5,7}` base 8 and `0.272` against `0.2513` for `{0,2,3,4,6}` base 7, `6.3` and `4.3` standard errors above, each flat over its last four levels. At `nA = 2.3` and `3.5` the gaps of these three designs are `2.6` to `4.6` standard errors on series that still move (`{0,2,3,4,6}` base 7 falls from `0.327` to `0.285` at `nA = 3.5`), so there the limit is open.

The torus mean is `ln abs(Q_D(u^2, v^2))` averaged over every frequency at the last level, dark or open; it starts at `m(P_D)`, where `u = v`, and moves toward `m(Q_D)` as the phases part. At `nA = 1.6` it ends on `m(Q_D)` to two digits in all four positive rows, consistent with `(u, v)` covering the torus evenly over all frequencies; at `nA = 2.3` and `3.5` it holds flat below `m(Q_D)` on the three refutation designs. The deep drift is the same kind of average taken only where the stack is still transparent, and there the mean lift `ln abs(Q_D(u^2, v^2))` tracks the mean step to `0.003` over the last four levels (`deep`), so the gap from `m(Q_D)` is the law of `(u, v)` on the passband, not the step. What survives is the zero set. **The deep drift is at most `0.01` in absolute value, standard error at most `0.01`, at all three contrasts on the three designs whose `Q_D` is Boyd-cyclotomic, and at least `0.26`, standard error at most `0.02`, on the four whose `Q_D` is not. Verified** (`deep`).

For `{0,1,3}` base 4, `Q_D = 1 + x + x^2 y`, and the unimodular change `y -> x^(-2) y` makes it `1 + x + y`, so `m(Q_D) = 3 sqrt(3) L(chi_(-3), 2)/(4 pi) = 0.3230659472`, Smyth's constant ([Smyth 1981](../refs.md)). **`m(Q_(0,1,3))` computed on the torus equals Smyth's constant to `1e-8`. Verified** (`smyth`, against PARI's `lfun(-3, 2)`). **"At `nA = 1.6` the `{0,1,3}` base 4 stack drifts at Smyth's constant." Refuted** (`deep`): the deep drift over levels 8 to 11 is `0.51` to `0.60`, standard error at most `0.026`, while the torus mean over all frequencies ends at `0.323`. At `nA = 2.3` and `3.5` the deep set falls below 300 frequencies by level 3, too soon for a drift.

## Two regimes

**Conjecture: when `m(Q_D) = 0` the passband fraction falls as `f_k ~ C/sqrt(k) = C/sqrt(log_b L)` in the number of cells `L = b^k`; when `m(Q_D) > 0` it falls geometrically in `k`, a power of `L`.** By Boyd's theorem `m(Q_D) = 0` exactly when `Q_D` is a monomial times cyclotomic polynomials in monomials, which covers the progressions and sum sets such as `{0,1,3,4}`, where `Q_D = (1 + x)(1 + x^2 y)`. The heuristic is the step law with a driftless walk: `ln abs(r/t)` performs a walk of mean `0` and the stack stays open while the walk stays below `0`, which a centred walk does for `k` steps with probability of order `1/sqrt(k)`.

**Every design at bases 3 to 6 sorts by the Boyd class. Verified** (`census`): the 39 designs with `0` in `D`, `2 <= abs(D) < b`, up to shift and mirror, at `nA = 2.3` and `1e5` frequencies. The flatness `f_12 sqrt(12)/(f_6 sqrt(6))` is `0.921` to `1.036` on the 26 Boyd-cyclotomic designs and `0.057` to `0.315` on the other 13; the verdict is critical above `0.9`, below the predicted value `1`, and it matches the class on all 39. At these bases no design has `P_D` cyclotomic with `Q_D` not, so the census cannot tell `m(P_D)` from `m(Q_D)`; the three designs of the refutation do. At level 12, `f_k sqrt(k)` is `1.097` for `{0,2,4,6}` base 7, still falling slowly from `1.108` at level 6, and `0.901` for `{0,1,3,4}` base 7 (`refute`); all of this is at the one contrast `nA = 2.3`.

## Open

- The joint law of the block phase `u` and the spacer phase `v`: the step law is exact to first order, and every gap above sits in how `(u, v)` spreads on the torus over the passband; a proof of either regime needs that law.
- The deep drift at `nA = 2.3` and `3.5` on the three refutation designs: it still moves at the last measured level, so whether it tends to `m(Q_D)` there is not settled.
- The critical walk is a dependent walk driven by `delta -> b delta`; a Sparre Andersen bound for it would prove `1/sqrt(k)` for the progressions, where the step is exact.
- The question the page hands to optics: which design with `m(Q_D) > 0` darkens slowest. The smallest `m(Q_D) > 0` in the census is Smyth's `0.3231` and the refutation designs reach `0.2513`; whether designs can drive it toward `0`, or toward Lehmer's constant, is Lehmer's problem restricted to staircase `0/1` polynomials, read as the slowest-localising fractal mirror; the one-variable problem already has a spectral reading in [Baake, Coons and Manibo 2020](../refs.md).
- Related pages: [crop](../../crop/notes/crop.md) and [weights](../../weights/notes/weights.md) for the digit transform whose Riesz product is the Born drift, and the Rajchman property and the [random walk](https://en.wikipedia.org/wiki/Random_walk) for the two regimes.
