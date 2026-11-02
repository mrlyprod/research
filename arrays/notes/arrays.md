---
title: Sparse arrays
lead: A fractal sparse array is a digit design in the balanced base of its generator's coarray, so its coarray weight is a product over digits and, from order 2 on, its fragility is exact: a sensor is essential exactly when every digit has a unique-difference partner, the fragility is `(u/L)^r`, the power bound of Yang and coauthors is tight exactly when the generator's essential sensors all have such partners, and the bound of Cohen and Eldar exactly when every sensor has one; at a compressed base the essential count can leave the product law.
figure: research-arrays
slug: arrays
---

## The objects

A sparse array is a finite set `S` of integer sensor positions; its difference coarray is `D(S) = S - S`, and its weight `w_S(t)` counts the ordered pairs `(p, q)` in `S^2` with `p - q = t`, so `w_S(0) = card S`. A sensor `s` is essential when removing it changes the coarray, `D(S - {s}) != D(S)`, the fragility is the share of essential sensors, and `S` is maximally economic when every sensor is essential. The words are [Cohen and Eldar 2020](https://arxiv.org/abs/2001.01217), Definitions 9 to 11, who take them from Liu and Vaidyanathan; the arrays serve direction finding, where a lost lag breaks the coarray estimator.

- The generator `G` is a set of `L >= 2` sensors, shifted so that `min G = 0`, whose coarray is hole-free: `D(G) = [-a, a]` with `a = max G`. The base is `M = 2a + 1 = card D(G)`, the translation factor of Cohen and Eldar.
- The fractal array of order `r` is `F_r = {sum_(i<r) g_i M^i : g_i in G}`, the dim 1 digit design with base `M` and digit set `G` at level `r`; Cohen and Eldar grow it as `F_(r+1) = union over g in G of (F_r + g M^r)`, which is the same set. The [Cantor set](https://en.wikipedia.org/wiki/Cantor_set) array is `G = {0, 1}`, `M = 3`.
- `U(G)` is the set of paired sensors, `g in G` with some `h in G` such that `w_G(g - h) = 1`, and `u = card U(G)`. Cohen and Eldar's condition C1 says every sensor is paired, `U(G) = G`; their Lemma 1 (from Liu and Vaidyanathan) says a paired sensor is essential. `E(G)` is the set of essential sensors of `G`, so `U(G)` sits inside `E(G)`.
- Every number below is printed by the study `lab/py/fractal-array-fragility`, whose sweep is every hole-free generator with `L <= 6` and span `a <= 13`, one of each mirror pair: 119 generators. Its test of essentialness is the definition itself, remove the sensor and recompute the coarray.

## The weight is a product over digits

**Every lag `t` of `F_r` has one expansion `t = sum_(i<r) t_i M^i` with balanced digits `t_i in [-a, a]`, and `w_(F_r)(t) = prod_(i<r) w_G(t_i)`. Proved.** A difference of two sensors is `sum (g_i - h_i) M^i` with every digit difference in `[-a, a]`, and `M = 2a + 1` is exactly the number of balanced digits, so the expansion is the balanced base-`M` expansion and is unique; an ordered pair of sensors with difference `t` is therefore a tuple of `r` ordered digit pairs, the `i`-th with difference `t_i`, and they are counted independently. So the coarray is `[-(M^r - 1)/2, (M^r - 1)/2]`, hole-free, which is Cohen and Eldar's Theorem 1, and the weight is the digit-by-digit product that their Theorem 4 writes as a convolution; in the frequency variable it is the product `prod_i abs(P_G(M^i x))^2` of the generator's array factor at the scales `M^i`, the fractal product beampattern.

The Cantor array is the case `G = {0, 1}`: `w_G(0) = 2` and `w_G(1) = w_G(-1) = 1`, so the weight of `t` is `2` to the number of zero digits of `t` in balanced ternary, the `R_A(t) = 2^(z_3(t))` of the energy reduction in [cobham](../../cobham/notes/cobham.md), Object S. The sumset of that note and the array of this one carry the same product.

## Fragility is exact

**Theorem A. Let `G` be hole-free on `[-a, a]`, `M = 2a + 1` and `r >= 2`. A sensor `s` of `F_r` is essential exactly when every base-`M` digit of `s` lies in `U(G)`. So `F_r` has exactly `u^r` essential sensors and its fragility is `(u/L)^r`. Proved.**

- Removing `s` deletes the lag `t != 0` exactly when every ordered pair with difference `t` contains `s`; those pairs are at most two, `(s, s - t)` and `(s + t, s)`. So `s` is essential exactly when some `t != 0` has `w(t) <= 2` with all its pairs through `s`.
- If all digits of `s` are paired, pick `h_i` with `w_G(s_i - h_i) = 1` and put `h = sum h_i M^i`. The lag `s - h` has every digit weight `1`, so weight `1`, and its one pair is `(s, h)`: `s` is essential.
- If `w(t) = 1` and its pair is `(s, q)` or `(q, s)`, every digit pair is the unique pair for `t_i`, and `t_i != 0` because `w_G(0) = L >= 2`; so every digit of `s` is paired.
- If `w(t) = 2` with both pairs through `s`, the pairs are `(s, s - t)` and `(s + t, s)`, and the product puts weight `2` on one digit `j` and weight `1` on every other. Pick `i != j`, which needs `r >= 2`: the two pairs share the unique digit pair for `t_i`, so `s_i = (s + t)_i` and `(s - t)_i = s_i`, which gives `t_i = 0` and `w_G(t_i) = L >= 2`, a contradiction. So weight `2` never makes a sensor essential, and the three cases close the proof.

**Verified** against the definition: on all 119 generators at `r = 2` and `r = 3`, 238 cases up to `216` sensors, the set of sensors whose removal changes the recomputed coarray equals the set of sensors with all digits in `U(G)`, with 0 mismatches (`lab/py/fractal-array-fragility`, verb `exact`). The theorem is false at `r = 1`, and that is the whole gap the next section measures: `G = {0, 1, 2}` is maximally economic, its middle sensor killing the lags `+-1` of weight `2`, yet `U(G) = {0, 2}`, so `F_2` at `M = 5` has 4 essential sensors of 9 and `F_3` has 8 of 27.

## The two bounds

Two bounds are in print, both read at source. [Cohen and Eldar 2020](https://arxiv.org/abs/2001.01217), Theorem 5: if `G` satisfies C1 then so does `F_r`, hence `F_r` is maximally economic; Theorem 6: the fragility of `F_r` is at most that of `G`. [Yang, Shen, Liu, Eldar and Cui 2023](https://ieeexplore.ieee.org/document/10043745/), Part II, Lemma 1: an inessential sensor of `G` or of `F_r` passes its inessentialness to every sensor of `F_(r+1)` built on it; Proposition 1: the fragility of `F_r` is at most the `r`-th power of that of `G`. Both are one-way.

- **The strong converse of Theorem 5: for `r >= 2`, `F_r` is maximally economic exactly when `G` satisfies C1, and exactly then Theorem 6 is attained, `(u/L)^r = card E(G)/L`. Proved** (Theorem A: every sensor is essential exactly when `u^r = L^r`, that is `U(G) = G`; if `u < L` then `(u/L)^r < u/L <= card E(G)/L`, and `u >= 2` because the lag `a` has the one pair `(a, 0)`). The literal converse, C1 for `F_r` forcing C1 for `G`, holds at every order, a lag of weight `1` having every digit of weight `1`; the strong form needs `r >= 2` and fails at `r = 1` on `G = {0, 1, 2}`, maximally economic without C1.
- **For `r >= 2`, Proposition 1 is tight, `(card E(G)/L)^r = (u/L)^r`, exactly when every essential sensor of `G` is paired, `E(G) = U(G)`; otherwise it is strictly loose at every `r >= 2`. Proved** (Theorem A, with `U(G)` inside `E(G)`). The verb `exact` prints the criterion with 0 failures on the 238 cases, which follows from its Theorem A check and is not separate evidence.
- A sensor in `E(G)` but not in `U(G)` is the middle of a three-term progression `g - t, g, g + t` in `G` whose step has weight exactly `2`: the weight-2 case of the proof, which is live at `r = 1` and dies at every `r >= 2`.
- **Of the 119 generators, 21 are loose, and 14 of those 21 are maximally economic, so both published bounds print `1` while the fragility is `(u/L)^r` and tends to `0`. Verified** (`lab/py/fractal-array-fragility`, verb `exact`, which lists the 21). The 14 run from `{0, 1, 2}` and `{0, 1, 2, 4}` through `{0, 1, 4, 5, 7}` to generators with six sensors such as `{0, 1, 3, 5, 8, 9}`, and on every one of the 21 the excess `card E(G) - u` is `1`, as the study prints.
- The minimum redundancy array `{0, 1, 4, 6}` satisfies C1 and stays at fragility `1` at every order, where both bounds are tight. `{0, 1, 2, 3}` has `u = 2 = card E(G)`, so Proposition 1 is tight while Theorem 6 is not: `F_2` at `M = 7` has 4 essential sensors of 16, fragility `1/4` against the bound `1/2`.

## The compressed base

Cohen and Eldar set the base by the generator's coarray, `M = 2a + 1`, where the earlier fractal arrays they cite take any natural translation factor; the digit design allows any base `b` with `a < b <= 2a + 1`, and below `2a + 1` a lag has many balanced expansions, joined by carries.

- **For `a < b <= 2a + 1` the digit design with base `b` and digits `G` at level `r` has `L^r` sensors and the hole-free coarray `[-A, A]`, `A = a(b^r - 1)/(b - 1)`. Proved.** Digits below `b` make the sensors distinct; the lags at level `r + 1` are `d + b x` with `d` in `[-a, a]` and `x` in `[-A_r, A_r]`, and consecutive `x` give intervals of length `2a + 1 >= b` that abut or overlap. **Verified** on every row below (`lab/py/fractal-array-fragility`, verb `dial`, which asserts both).
- The essential count `e_r(G, b)` is no longer a product. **Verified** (`lab/py/fractal-array-fragility`, verb `dial`, up to 16384 sensors and 8 million lags, the fast test checked against the definition up to 250 sensors):

| `G` | `b` | `e_r` at `r = 1, 2, ...` |
| --- | --- | --- |
| `{0, 1, 2}` | 3, 4 | `3, 2, 2, 2, 2, 2, 2, 2` |
| `{0, 1, 2}` | 5 | `3`, then `2^r` |
| `{0, 1, 4, 6}` | 7 | `4, 7, 6, 6, 6, 6, 6` |
| `{0, 1, 4, 6}` | 8 | `4, 11, 25, 53, 109, 221, 445` |
| `{0, 1, 4, 6}` | 9 to 13 | `4^r` |
| `{0, 1, 2, 3, 7}` | 8 to 10 | `5, 7, 6, 6, 6, 6` |
| `{0, 1, 2, 3, 7}` | 11 | `5, 14, 29, 61, 125, 253` |
| `{0, 1, 2, 3, 7}` | 12 | `5, 17, 53, 161, 485, 1457` |
| `{0, 1, 2, 3, 7}` | 13 | `5, 21, 85, 341, 1365, 5461` |
| `{0, 1, 2, 3, 7}` | 14, 15 | `5^r` |
| `{0, 1, 4, 7, 9}` | 10 | `5, 9, 8, 8, 8, 8` |
| `{0, 1, 4, 7, 9}` | 11 | `5, 16, 35, 75, 155, 315` |
| `{0, 1, 4, 7, 9}` | 12, 13 | `5, 22, 90, 362, 1450, 5802` |
| `{0, 1, 4, 7, 9}` | 14 to 19 | `5^r` |
| `{0, 1, 2, 3, 7, 11}` | 12 to 14 | `6, 7, 7, 7, 7` |
| `{0, 1, 2, 3, 7, 11}` | 15, 16 | `6, 16, 36, 76, 156` |
| `{0, 1, 2, 3, 7, 11}` | 17 | `6, 26, 106, 426, 1706` |
| `{0, 1, 2, 3, 7, 11}` | 18, 19 | `6, 31, 156, 781, 3906` |
| `{0, 1, 2, 3, 7, 11}` | 20 to 23 | `6^r` |

- **Every row satisfies, from `r = 2` or `r = 3` on, an affine law `e_(r+1) = lambda e_r + c` with an integer `lambda` between `1` and `u`: bounded at `lambda = 1` with `c = 0`, geometric at `1 < lambda < u`, and the full `u^r`. Conjecture** (the verb `dial` fits `lambda >= 1` and `c` on two steps, taking `lambda = 1` on a flat tail, no upper bound, then checks the law on the next 1 to 4 steps; every fitted `lambda` is an integer in `[1, u]`, and the ratios seen are `1, 2, 3, 4, 5`). `{0, 1, 4, 6}` at base 8 is `e_r = 7 2^(r-1) - 3` and `{0, 1, 2, 3, 7, 11}` at base 15 is `5 2^r - 4` on the printed terms.
- The full law arrives before `2a + 1`: at `b = 9`, `14`, `14` and `20` for spans `6`, `7`, `9` and `11`, against `M = 13`, `15`, `19` and `23`. Theorem A's proof reaches no further than uniqueness of the balanced expansion, so these rows are data, not theorem.
- What the dial buys, said plainly: `{0, 1, 2}` at base 4 has `3^r` sensors, `2` essential ones at every `r` from `2` to `8`, the fragility `2/N` of the uniform line at base 3, `(4^(r+1) - 1)/3` hole-free lags, a product array factor `prod_i P_G(4^i x)` and a self-similar layout; at `r = 8` that is 6561 sensors and 87381 lags. That is not a record for lags: they grow as `N^(log 4/log 3)`, about `N^1.26`, while the minimum redundancy array reaches order `N^2` lags and is maximally economic ([Liu and Vaidyanathan 2019](https://ieeexplore.ieee.org/document/8695867/), Part II). Whether a closed-form array has fragility `2/N` and order `N^2` lags is not settled here; the gain is exact robustness in closed form with a product beampattern at every size.

## What is open

- The thresholds of the compressed base: the least `b` where the full `u^r` law starts, and where the geometric regime starts, as functions of `G`. A lag at base `b < 2a + 1` is read by a carry automaton on digit pairs, in the manner of the transfer matrices of [beneath](../../beneath/notes/beneath.md) and [transfer matrix](https://en.wikipedia.org/wiki/Transfer-matrix_method_(statistical_mechanics)), so `e_r(G, b)` should have a rational generating function; neither the automaton nor the function is built here.
- The `k`-essential family: the sets of `k` sensors whose joint removal changes the coarray while no smaller part does, the `k`-essential Sperner family of [Liu and Vaidyanathan 2019](https://ieeexplore.ieee.org/document/8695867/). Their Part II, Lemma 1, gives it for every maximally economic array, the Cantor arrays among them: the single sensors and nothing larger; so for a generator with C1 it follows from Theorem 5 of Cohen and Eldar. For the others, a lag of `F_r` loses all its pairs to `k` removals only if its weight is at most `2k`, which is the product law again, so the family should be a product of per-digit families at `r` large enough; it is not attempted here.
- The higher-order fractals of Yang, Shen, Liu, Eldar and Cui, where the coarray is a sum of `q` differences; their Proposition 2 is again an upper bound.

## Where the numbers live

- `lab/py/fractal-array-fragility`: the verb `exact` (Theorem A against the definition, the tightness criterion, the 21 loose generators) and the verb `dial` (the compressed base table and its affine fits).
- The sparse arrays demo builds `F_r` of any generator at any base `a < b <= 2a + 1` through the crate, knocks out one sensor to open the lags it alone carried, and sets the counted essential sensors beside `(u/L)^r` and the two bounds.
- The figure is `F_4` of `{0, 1, 2, 4}` at `M = 9`, folded into the plane by alternating digits, its 81 essential sensors apart from the other 175; the binary asserts Theorem A on all 256 sensors against the definition.
- The paper [The exact fragility of fractal sparse arrays](../papers/fractal-array-fragility.md) writes Theorem A and the two bounds out for an array engineer: the two-pair test, the product identity, the proof, the equality cases against the statements at source, Table I of Cohen and Eldar read through the law, two worked examples and the compressed base.
