---
title: The Exact Fragility of Fractal Sparse Arrays
lead: From order 2 on, a sensor of a fractal sparse array is essential exactly when every base-`M` digit of its position has a partner in the generator at a lag that no other pair makes, so the fragility is exactly `(u/L)^r`; from order 2 on, the bound of Cohen and Eldar is attained exactly when every sensor of the generator has such a partner, which is the strong converse of their Theorem 5, and the power bound of Yang and coauthors exactly when every essential sensor does.
date: 2026-10-05
figure: paper-fractal-array-fragility
---

A sparse array grown from a small generator, the way Cohen and Eldar grow it, has a hole-free difference coarray at every size, and the next question is how fragile it is: how many of its sensors are essential, so that losing one opens a hole in the coarray. The answers in print are upper bounds: the fractal array is at most as fragile as its generator (Cohen and Eldar 2020, and for every composite array Liu and Vaidyanathan 2019), and at most as fragile as the generator's fragility raised to the order (Yang, Shen, Liu, Eldar and Cui 2023). This paper gives the exact value. Let `G` be a generator of `L` sensors with hole-free coarray `[-a, a]`, `M = 2a + 1`, and `F_r` its fractal array of order `r`, whose sensors are the numbers with `r` base-`M` digits drawn from `G`. Call a sensor of `G` paired when it and some other sensor of `G` make a lag that no other pair makes, and let `u` count them. From order `2` on, a sensor of `F_r` is essential exactly when all its digits are paired, so `F_r` has exactly `u^r` essential sensors and fragility `(u/L)^r`. Read against the published bounds, from order `2` on `F_r` is maximally economic exactly when every sensor of `G` is paired, which is the strong converse of Cohen and Eldar's Theorem 5 and the one case where their Theorem 6 is attained, and the bound of Yang and coauthors is attained exactly when every essential sensor of `G` is paired. Everything rests on one identity: the coarray weight of `F_r` at a lag is the product of the generator's weights at the balanced base-`M` digits of the lag. Of the 119 hole-free generators with at most six sensors and span at most 13, one of each mirror pair, 21 leave the power bound strictly loose and 14 of those are maximally economic, so both published bounds read `1` at every order while the fragility falls geometrically; the smallest is the line `{0, 1, 2}`, with fragility `(2/3)^r`. With the generators' own fragilities at order `1`, the law reproduces the six fragilities of Table I of Cohen and Eldar, and it explains the gap that Table I of Yang and coauthors already prints for one of the two generators, fragility `0.03` against their bound `0.0729` at order `2`. Below the base `2a + 1` the coarray stays hole-free and the essential count can leave the product law.

## Introduction

![The fractal arrays of the generator {0, 1, 2, 4} at base 9, orders 1 to 4, as four stacked rows of 4, 16, 64 and 256 cells, one cell a sensor, each row splitting every cell of the row above into its four children: essential sensors orange, inessential blue. The top row is all orange, every sensor of the generator essential, yet its third cell, the sensor 2, has only blue children; from the second row on the orange cells are exactly those whose every digit is 0, 1 or 4, 9, 27 and 81 of them.](paper-fractal-array-fragility)

The figure is the theorem for one generator. `G = {0, 1, 2, 4}` has span `4`, a hole-free coarray `[-4, 4]` with weights `4, 2, 2, 1, 1` at the lags `0` to `4`, and base `M = 9`. Each of its four sensors is essential: `0` and `4` make the lag `4` alone, `1` and `4` the lag `3`, and `2` is the middle of `0, 2, 4`, so both pairs that make the lag `2` pass through it. The generator is maximally economic and the top row is orange. Each row below splits every sensor into four, one for each digit appended at the next power of `9`, and the cells that stay orange are those whose every digit is `0`, `1` or `4`: 9 of 16, 27 of 64, 81 of 256, the fraction `(3/4)^r`. The sensor `2` loses its standing at once: from order `2` on, a lag that carries the generator's lag `2` in one digit has either more than two pairs or two that share no sensor, so no single loss removes it. A blue cell has only blue children, the inheritance of inessential sensors proved by Liu and Vaidyanathan and by Yang and coauthors; what the figure adds is that below the top row the children of an orange cell are orange exactly at the paired digits, while the top row's `2` has only blue children. Both published bounds see the top row and say `1`.

The fractal arrays of [Cohen and Eldar 2020](https://arxiv.org/abs/2001.01217) answer a design need: a closed form at any size, a hole-free coarray (their Theorem 1), a weight function and a beampattern in product form (their Theorem 4). Robustness enters through the essential sensors and the fragility, words they take from Liu and Vaidyanathan. Their Theorem 5 says that condition (C1), every sensor at a lag of weight one with some other, passes from the generator to the fractal array, which is then maximally economic; their Theorem 6 that the fractal array is at most as fragile as its generator. [Yang, Shen, Liu, Eldar and Cui 2023](https://doi.org/10.1109/TSP.2023.3244667), Part II, sharpen Theorem 6 to the `r`-th power through an inheritance lemma, an inessential digit or an inessential lower part making a sensor inessential, and carry the bound to the arrays of higher-order cumulants. The inheritance itself is older: [Liu and Vaidyanathan 2019](https://doi.org/10.1109/ICASSP.2019.8683563) prove it for any composite array `alpha S_1 + S_2` with `alpha > max S_2`, of which `F_(r+1) = M F_r + G = M^r G + F_r` is a case both ways, and their bound, the composite at most as fragile as `S_1`, gives Theorem 6. These statements run one way: Theorem 5 is a sufficient condition for every sensor to be essential, the others bound the essential sensors from above, and none says which sensors are essential when (C1) fails.

Theorem A of Section 4 says which. Two corollaries in Section 5 settle the equality cases of both bounds and the strong converse of Theorem 5, and they show that one finite object decides everything: the sensors of `G` that are essential only as the middle of a three-term progression whose step has weight `2`. Section 2 fixes the words, Section 3 proves the product identity, Section 6 works two examples with numbers the reader can redo by hand, Section 7 shows what happens below the base `2a + 1`, where the proof of Theorem A stops, Section 8 states what is open and Section 9 names every generator. The sparse arrays demo builds `F_r` of any generator at any base above `a`, knocks out the sensor one clicks to open the lags it alone carried, and counts the essential sensors beside `(u/L)^r` and the two bounds. None of the three papers, read at source, states a converse of Theorem 5 or an equality case for either bound, though Table I of Yang and coauthors prints a case where theirs is not attained (Section 5). A search beyond them finds a test of the robustness of sparse fractal arrays to a faulty sensor ([Goel, Aggarwal and Kar 2022](https://arxiv.org/abs/2212.00341)) and no exact count of essential sensors; it is not exhaustive.

## The setting

**Definition 2.1 (coarray and weight).** A sparse array is a finite set `S` of integers, the sensor positions in units of the least spacing. Its difference coarray is `D(S) = {p - q : p, q in S}`, and its weight `w_S(t)` counts the ordered pairs `(p, q)` in `S^2` with `p - q = t`, so `w_S(0) = card S` and `w_S(-t) = w_S(t)`. The coarray is hole-free when it is an interval of integers.

**Definition 2.2 (essential, maximally economic, fragility).** A sensor `s` of `S` is essential when `D(S - {s}) != D(S)`; `S` is maximally economic when every sensor is essential; the fragility of `S` is `card E(S) / card S`, with `E(S)` the set of essential sensors. These are Definitions 9, 10 and 11 of Cohen and Eldar.

**Definition 2.3 (the fractal array).** A generator `G` is a set of `L >= 2` integers with `min G = 0` and hole-free coarray `[-a, a]`, `a = max G`; its base is `M = 2a + 1 = card D(G)`, the translation factor of Cohen and Eldar. The fractal array of order `r >= 1` is

```
F_r = { sum_(i<r) g_i M^i : g_i in G } ,
```

the recursion (5) of Cohen and Eldar, `F_0 = {0}` and `F_(r+1) = union over g in G of (F_r + g M^r)`, unrolled; `F_1 = G`, and `G = {0, 1}` at `M = 3` is the Cantor array. The digits `g_i` lie below `M`, so the `L^r` expansions are distinct and `card F_r = L^r`; `g_i` is the `i`-th digit of the sensor, the order is the level of the digit design. Cohen and Eldar state (5) for any generator, with `M` the size of the central uniform segment of its coarray; their Theorems 5 and 6 and everything below assume the coarray hole-free.

**Definition 2.4 (paired sensors, condition C1).** A sensor `g` of `G` is paired when some `h` in `G` has `w_G(g - h) = 1`. `U(G)` is the set of paired sensors and `u = card U(G)`. Condition (C1) of Cohen and Eldar is `U(G) = G`. Only `0` and `a` differ by `a`, so `w_G(a) = 1` and `u >= 2`.

**Lemma 2.5 (the two-pair test).** Let `S` have at least two sensors and let `s` lie in `S`. Removing `s` deletes a lag `t` from the coarray exactly when `t != 0` and every pair with difference `t` contains `s`, and those pairs are at most two, `(s, s - t)` and `(s + t, s)`. So `s` is essential exactly when some lag `t != 0` has `w_S(t) = 1` with its one pair through `s`, or `w_S(t) = 2` with its two pairs `(s, s - t)` and `(s + t, s)`, that is, with `s - t`, `s` and `s + t` all in `S`.

**Proof.** The lag `0` survives, since `S - {s}` is not empty. For `t != 0` a pair `(p, q)` with `p - q = t` contains `s` when `p = s` or `q = s`, which makes it `(s, s - t)` or `(s + t, s)`, two different pairs; the lag survives exactly when some pair avoids `s`. □

The first case is Lemma 1 of Cohen and Eldar, after Liu and Vaidyanathan: a paired sensor is essential, so `U(G)` lies inside `E(G)`. The second case is the only other way to be essential, and it names the difference: a sensor of `G` is essential but unpaired exactly when no lag of weight `1` ends at it and it is the middle of a progression `g - t, g, g + t` in `G` whose step `t` has weight exactly `2`.

## The weight is a product over digits

**Lemma 3.1 (the product identity).** Every lag `t` of `F_r` has one expansion `t = sum_(i<r) t_i M^i` with balanced digits `t_i` in `[-a, a]`. An ordered pair `(p, q)` of sensors has difference `t` exactly when its digit pairs `(p_i, q_i)` have difference `t_i` at every `i`, and

```
w_(F_r)(t) = prod_(i<r) w_G(t_i) .
```

**Proof.** A difference of two sensors is `sum (p_i - q_i) M^i` with every digit difference in `[-a, a]`. The `M = 2a + 1` balanced digits form a complete system of residues mod `M`, so the balanced base-`M` expansion of an integer is unique: the residue of `t` fixes `t_0`, then that of `(t - t_0)/M` fixes `t_1`, and so on. Hence `p - q = t` forces `p_i - q_i = t_i` at every `i`, and an ordered pair with difference `t` is a tuple of `r` ordered digit pairs, the `i`-th with difference `t_i`, chosen independently. □

So the coarray of `F_r` is the set of integers whose balanced digits all lie in `D(G) = [-a, a]`, the whole interval `[-(M^r - 1)/2, (M^r - 1)/2]`, which is Theorem 1 of Cohen and Eldar; the product is their Theorem 4, which writes `w_(F_r)` as the convolution of the copies of `w_G` stretched by `M^i`, and in the frequency variable the beampattern of `F_r` is the product of the generator's at the scales `M^i`. In `F_2` of `{0, 1, 2, 4}` at `M = 9`, the lag `23 = -4 + 3 * 9` has weight `w_G(-4) w_G(3) = 1`, and its one pair is `(36, 13)`, with digits `(0, 4)` and `(4, 1)`, lowest first, differing by `(-4, 3)`.

## Fragility is exact

**Theorem A (the fragility law).** Let `G` be a generator with hole-free coarray `[-a, a]`, `M = 2a + 1` and `r >= 2`. A sensor `s` of `F_r` is essential exactly when every base-`M` digit of `s` lies in `U(G)`. So `F_r` has exactly `u^r` essential sensors, and its fragility is `(u/L)^r`.

**Proof.** By Lemma 2.5, `s` is essential exactly when some lag `t != 0` has weight `1` with its pair through `s`, or weight `2` with both its pairs through `s`. Write `s_i` and `t_i` for the digits.

Paired digits suffice. If every `s_i` is paired, choose `h_i` in `G` with `w_G(s_i - h_i) = 1` and put `h = sum h_i M^i`, a sensor of `F_r`. By Lemma 3.1 the lag `s - h`, nonzero since `s_0 != h_0`, has weight `1` and its one pair is `(s, h)`, so `s` is essential.

A lag of weight 1 needs paired digits. If `w(t) = 1` and its pair is `(s, q)` or `(q, s)`, Lemma 3.1 gives `w_G(t_i) = 1` at every `i`, and `t_i != 0` because `w_G(0) = L >= 2`; the digit pair `(s_i, q_i)` or `(q_i, s_i)` is then the one pair with difference `t_i`, so `s_i` is paired, at every `i`.

A lag of weight 2 never decides. Suppose `w(t) = 2` with both pairs through `s`, so the pairs are `(s, s - t)` and `(s + t, s)`. By Lemma 3.1 the product is `2` only if one digit `j` has `w_G(t_j) = 2` and every other has weight `1`. Take `i != j`, which exists because `r >= 2`. The `i`-th digit pairs of both sensor pairs have difference `t_i`, and `w_G(t_i) = 1` leaves only one such digit pair, so `(s_i, (s - t)_i) = ((s + t)_i, s_i)`. Then `(s - t)_i = s_i`, so `t_i = 0` and `w_G(t_i) = L >= 2`, a contradiction.

So `s` is essential exactly when a lag of weight `1` ends at it, exactly when all its digits are paired, and the essential sensors are the `u^r` digit strings over `U(G)`. □

The proof never uses that the coarray of `G` is hole-free, and neither does that of Lemma 3.1, only `min G = 0`, `max G = a` and `M = 2a + 1`; the hypothesis makes `M` the translation factor of Cohen and Eldar.

The proof shows a little more: from order `2` on, the lags that the loss of a sensor opens are exactly the lags of weight `1` ending at it, read off the generator's weights digit by digit. In `F_2` of `{0, 1, 2, 4}` the loss of `13 = 4 + 1 * 9` opens `+-23` and `+-24`, the lags to `36 = 0 + 4 * 9` and `37 = 1 + 4 * 9`, the two sensors whose digits meet `4` and `1` at weight `1`; the loss of `11 = 2 + 1 * 9` opens nothing.

At order `1` the theorem is false, because the weight-2 case has no second digit to contradict it: the middle sensor of `{0, 1, 2}` deletes the lags `+-1`, of weight `2`, yet `U = {0, 2}`. The sensors of `E(G) - U(G)` are the whole difference between order `1` and order `2`.

**Fact 4.1 (Theorem A against the definition).** On all 119 hole-free generators with `L <= 6` and span `a <= 13`, one of each mirror pair `G`, `a - G`, at `r = 2` and `r = 3`, 238 arrays of up to 216 sensors, the set of sensors whose removal changes the recomputed coarray equals the set of sensors with every digit in `U(G)`, with 0 mismatches. Domain: those 238 arrays. Generator: `lab/py/fractal-array-fragility`, verb `exact`, whose test is the definition itself, remove the sensor and recompute the coarray.

## The published bounds at source

Cohen and Eldar, read at source in arXiv:2001.01217v2, Section IV.C. **Theorem 5:** let `F_r` be the fractal array generated from `G` whose difference coarray is hole-free; then `F_r` satisfies condition (C1) if `G` satisfies it, and the proof closes with `F_(r+1)` maximally economic. **Theorem 6:** under the same hypothesis the fragility of `F_r` is at most that of `G` for every `r >= 1`; the proof shows that each inessential sensor of `F_r` yields `L` inessential sensors of `F_(r+1)`.

Yang, Shen, Liu, Eldar and Cui, Part II, read at source in the journal text, Section IV.A. **Lemma 1:** in `F_(r+1) = {l + n M^r : l in F_r, n in G}`, if `n_1` is inessential in `G` then `l + n_1 M^r` is inessential for every `l` in `F_r`, and if `l_2` is inessential in `F_r` then `l_2 + n M^r` is inessential for every `n` in `G`. **Proposition 1:** for `G` with hole-free coarray, the fragility of `F_r` is at most the fragility of `G` to the power `r`, for every `r >= 1`, their (12); the proof counts the two families of Lemma 1 by inclusion and exclusion, gets `N_e(k+1) <= N_ek N_eG` for the essential counts, and is presented as an extension of Theorem 6 of Cohen and Eldar.

Liu and Vaidyanathan, Composite Singer Arrays, read at source. **Definition 6:** for `min S_1 = min S_2 = 0` and an integer `alpha > max S_2`, the composite array is `S_c = alpha S_1 + S_2`. **Corollary 1:** its coarray is `alpha D_1 + D_2`. **Proposition 1:** `S_c` is at most as fragile as `S_1`, by its Property 1, an inessential `n_1` of `S_1` making `alpha n_1 + n_2` inessential for every `n_2` in `S_2`. With `S_1 = G`, `S_2 = F_r` and `alpha = M^r`, and again with `S_1 = F_r`, `S_2 = G` and `alpha = M`, these give the inheritance behind Lemma 1 of Yang and coauthors, from an inessential digit and from an inessential block of digits, and the first gives Theorem 6.

**Corollary 5.1 (the strong converse of Theorem 5, and when Theorem 6 is attained).** For `r >= 2` the following are equivalent: `F_r` is maximally economic; `F_r` satisfies (C1); `G` satisfies (C1). Exactly then Theorem 6 is attained, `(u/L)^r = card E(G)/L = 1`; otherwise the fragility is strictly below the bound, `(u/L)^r < card E(G)/L`.

**Proof.** If `G` satisfies (C1), the first step of the proof of Theorem A gives every sensor of `F_r` a partner at a lag of weight `1`, which is (C1) for `F_r` and Theorem 5; (C1) for `F_r` makes it maximally economic by Lemma 2.5; and if `F_r` is maximally economic, Theorem A gives `u^r = L^r`, so `U(G) = G`. If `u < L` then `u >= 2` gives `0 < u/L < 1`, so `(u/L)^r < u/L <= card E(G)/L`. □

The literal converse, (C1) for `F_r` forcing (C1) for `G`, holds at every order; the stronger form, maximal economy of `F_r` forcing (C1) for `G`, needs `r >= 2` and fails at order `1` on `{0, 1, 2}`, maximally economic without (C1). Nor can the hypothesis of Theorem 5 be weakened from (C1) to maximal economy: `{0, 1, 2, 4}` is maximally economic, and its `F_2` has 9 essential sensors of 16.

**Corollary 5.2 (when Proposition 1 is attained).** For `r >= 2` the bound `(card E(G)/L)^r` of Yang and coauthors equals the fragility `(u/L)^r` exactly when `E(G) = U(G)`, that is, when every essential sensor of `G` is paired; otherwise it is strictly loose at every `r >= 2`, by the factor `(card E(G)/u)^r`.

**Proof.** `U(G)` lies inside `E(G)` by Lemma 2.5, so `card E(G) = u` exactly when the two sets agree, and Theorem A gives the fragility. □

The two bounds fail for different reasons. Theorem 6 is loose whenever `u < L`, because the fragility falls as a power of the order. Proposition 1 has the power, and its slack is exactly `E(G) - U(G)`: through Lemma 1 it counts as possibly essential every sensor whose digits are all essential in `G`, and Theorem A strikes from that count every sensor with a digit that is essential in `G` only as the middle of a progression with a step of weight `2`.

**Fact 5.3 (how often the power bound is loose).** Among the 119 hole-free generators with `L <= 6` and span `a <= 13`, one of each mirror pair, 21 have `E(G) != U(G)`, and every one of the 21 has `card E(G) = u + 1`. Fourteen of them are maximally economic, so both published bounds equal `1` at every order while the fragility is `(u/L)^r`: `{0, 1, 2}`, `{0, 1, 2, 4}`, `{0, 1, 2, 3, 6}`, `{0, 1, 3, 5, 6}`, `{0, 1, 4, 5, 7}`, `{0, 1, 2, 3, 4, 8}`, `{0, 1, 2, 3, 6, 8}`, `{0, 1, 3, 4, 8, 9}`, `{0, 1, 3, 5, 8, 9}`, `{0, 1, 5, 6, 7, 9}`, `{0, 1, 2, 6, 7, 10}`, `{0, 1, 3, 5, 9, 10}`, `{0, 1, 3, 6, 8, 10}` and `{0, 1, 4, 5, 8, 10}`. The other seven are `{0, 1, 2, 4, 5}`, `{0, 1, 2, 3, 6, 7}`, `{0, 1, 2, 4, 6, 7}`, `{0, 1, 2, 4, 7, 8}`, `{0, 1, 2, 5, 6, 8}`, `{0, 1, 3, 4, 6, 8}` and `{0, 1, 4, 5, 7, 9}`. Domain: those 119 generators. Generator: `lab/py/fractal-array-fragility`, verb `exact`, which lists the 21.

The published numbers obey the law. Section V of Cohen and Eldar takes two generators of aperture `20` from a design problem, their (P1): `S = {0, 1, 2, 4, 7, 10, 13, 16, 18, 19, 20}`, a solution, and their `G`, here `G' = {0, 1, 3, 5, 11, 13, 17, 18, 19, 20}`, a solution without the symmetry requirement (R1). Their Table I prints the fragilities `0.27`, `0.03`, `0.006` for `S` and its fractal arrays of orders `2` and `3`, and `0.30`, `0.09`, `0.027` for `G'` and its. Yang and coauthors take the same two generators, `G'` as their (26) and `S` as their (27), and their Table I already prints Proposition 1 strictly loose on `S`: `0.03` against the bound `0.0729` at order `2` and `0.006` against `0.0197` at order `3`, the bounds powered from the rounded `0.27`; Table 1 gives the exact bounds `9/121` and `27/1331` and the reason, `E(S) != U(S)`. Both coarrays are hole-free with `a = 20`, `M = 41`, and the weights at the lags `0` to `20` are

```
S  : 11, 4, 4, 6, 2, 2, 5, 2, 2, 4, 2, 2, 3, 2, 2, 2, 3, 2, 3, 2, 1
G' : 10, 4, 5, 2, 2, 2, 3, 2, 3, 1, 2, 1, 2, 2, 2, 2, 2, 3, 2, 2, 1
```

For `S` the only lag of weight `1` is `20`, so `U(S) = {0, 20}`, while the lag `10` has the two pairs `(10, 0)` and `(20, 10)`, so `E(S) = {0, 10, 20}`: `S` is a published generator with `E != U`. For `G'` the lags of weight `1` are `9`, `11` and `20`, made by `(20, 11)`, `(11, 0)` and `(20, 0)`, so `U(G') = E(G') = {0, 11, 20}`.

| array | sensors | Table I | exact | Theorem 6 | Proposition 1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `S` | 11 | `0.27` | `3/11 = 0.2727` | `0.2727` | `0.2727` |
| `S_2` | 121 | `0.03` | `4/121 = 0.0331` | `0.2727` | `9/121 = 0.0744` |
| `S_3` | 1331 | `0.006` | `8/1331 = 0.0060` | `0.2727` | `27/1331 = 0.0203` |
| `G'` | 10 | `0.30` | `3/10` | `0.3` | `0.3` |
| `G'_2` | 100 | `0.09` | `9/100` | `0.3` | `9/100` |
| `G'_3` | 1000 | `0.027` | `27/1000` | `0.3` | `27/1000` |

**Table 1.** The fragilities of Table I of Cohen and Eldar against the exact value, `card E/L` at order `1` and `(u/L)^r` from order `2` on, and the two bounds computed exactly. Every exact entry is also a count of essential sensors by the two-pair test of Lemma 2.5: 3 of 11, 4 of 121, 8 of 1331, 3 of 10, 9 of 100 and 27 of 1000. For `S`, Proposition 1 is loose by `9/4` at order `2` and `27/8` at order `3`. Domain: `S` and `G'` at orders `1` to `3`. Generator: `mrly_math::arrays::essential` over `mrly_math::arrays::fractal`, the counts, both weight rows and both sets `U` asserted by the crate test `cohen_and_eldar_table_one_reads_as_printed`.

## Two worked examples

**Example 6.1 (the line of three, the smallest gap).** `G = {0, 1, 2}`, `a = 2`, `M = 5`. The weights at the lags `-2` to `2` are `1, 2, 3, 2, 1`, so the only lags of weight `1` are `+-2`, made by `(2, 0)`, and `U(G) = {0, 2}`, `u = 2`. The middle sensor makes both pairs of the lag `1`, `(1, 0)` and `(2, 1)`, so removing it deletes `+-1`: `E(G) = G`, the generator is maximally economic, and Theorem 6 and Proposition 1 both give `1` at every order. The array of order `2` is

```
F_2 = {0, 1, 2, 5, 6, 7, 10, 11, 12} ,   coarray [-12, 12] ,
```

and Theorem A picks the sensors with both digits in `{0, 2}`: `0 = (0, 0)`, `2 = (2, 0)`, `10 = (0, 2)` and `12 = (2, 2)`, digits lowest first, 4 of 9. By hand: the lag `12 = 2 + 2 * 5` has weight `1 * 1 = 1` and the one pair `(12, 0)`; the lag `8 = -2 + 2 * 5` has weight `1` and the one pair `(10, 2)`; so the four corners each carry a lag alone. The middle sensor `6 = (1, 1)` carries nothing alone: a lag with a digit `+-1` somewhere has weight at least `2`, and when it is exactly `2`, as for `11 = 1 + 2 * 5` with the pairs `(11, 0)` and `(12, 1)`, the two pairs share no sensor. The fragility is `4/9`, then `8/27`, `16/81`, `32/243`, `64/729` at orders `3` to `6`, which is `(2/3)^r`, against `1` from both bounds.

**Example 6.2 (the array the demo opens on).** `G = {0, 1, 2, 4}`, `a = 4`, `M = 9`. The weights at the lags `0` to `4` are `4, 2, 2, 1, 1`; the lags `4` and `3` are made by `(4, 0)` and `(4, 1)` alone, so `U(G) = {0, 1, 4}`, and `2` is essential through the progression `0, 2, 4` with the lag `2` of weight `2`: `E(G) = G`. The array of order `2` has 16 sensors on `[0, 40]`,

```
F_2 = {0, 1, 2, 4, 9, 10, 11, 13, 18, 19, 20, 22, 36, 37, 38, 40} ,
```

its coarray is `[-40, 40]`, all 81 lags, and 16 of them have weight `1`, the products of the generator's four weight-1 lags `+-3`, `+-4` in both digits. Its essential sensors are the 9 with both digits in `{0, 1, 4}`, `0, 1, 4, 9, 10, 13, 36, 37, 40`; the sensor `2 = (2, 0)`, essential in `G`, is inessential here, its lag `2` now having weight `w_G(2) w_G(0) = 8`. The demo opens on this array with `13` knocked out and the four lags `+-23`, `+-24` open. The counts by the two-pair test of Lemma 2.5 at every order through `6`:

| `r` | sensors | essential | fragility | `(3/4)^r` | Theorem 6 | Proposition 1 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4 | 4 | `1` | | `1` | `1` |
| 2 | 16 | 9 | `0.5625` | `0.5625` | `1` | `1` |
| 3 | 64 | 27 | `0.4219` | `0.4219` | `1` | `1` |
| 4 | 256 | 81 | `0.3164` | `0.3164` | `1` | `1` |
| 5 | 1024 | 243 | `0.2373` | `0.2373` | `1` | `1` |
| 6 | 4096 | 729 | `0.1780` | `0.1780` | `1` | `1` |

**Table 2.** `F_r` of `{0, 1, 2, 4}` at `M = 9`: the essential sensors counted by `mrly_math::arrays::essential` over `mrly_math::arrays::fractal`, every coarray hole-free by `mrly_math::arrays::holes`, against Theorem A and the two published bounds; the crate test `theorem_a_counts_u_to_the_level` asserts the counts at orders `4` to `6` and the figure those at orders `1` to `4`.

Two more generators close the picture. The line of four, `{0, 1, 2, 3}` at `M = 7`, has `U = E = {0, 3}`, the interior sensors being inessential already; Proposition 1 is attained, `(1/2)^r`, and Theorem 6 is not: `F_2` has the 4 essential sensors `0, 3, 21, 24` of 16, fragility `1/4` against the bound `1/2`. The four-sensor minimum redundancy array `{0, 1, 4, 6}` at `M = 13` has every lag but `0` of weight `1`, so it satisfies (C1), and `F_r` stays maximally economic at every order, 16 of 16 and 64 of 64 at orders `2` and `3`, where both bounds are attained.

## The compressed base

Cohen and Eldar fix the base by the generator's coarray, `M = 2a + 1`, where the earlier fractal arrays they cite allow any natural translation factor. The digit design accepts any base `b` with `a < b <= 2a + 1`, and below `2a + 1` a lag has many balanced expansions, joined by carries.

**Proposition 7.1.** For `a < b <= 2a + 1` the set `{sum_(i<r) g_i b^i : g_i in G}` has `L^r` sensors and the hole-free coarray `[-A, A]`, `A = a (b^r - 1)/(b - 1)`.

**Proof.** The digits lie below `b`, so the sensors are distinct. Building the array by a new lowest digit, its sensors at level `r + 1` are `g + b f` with `g` in `G` and `f` at level `r`, so its coarray is `D(G) + b D_r = {d + b x : d in [-a, a], x in [-A_r, A_r]}`, granting the interval at level `r`. For consecutive `x` the intervals `[b x - a, b x + a]` have length `2a + 1 >= b` and abut or overlap, so the union is `[-(a + b A_r), a + b A_r]`, and `a + b A_r = a (b^(r+1) - 1)/(b - 1)`. □

**Fact 7.2 (the essential count at a compressed base).** The count `e_r(G, b)` of essential sensors is in general no longer `u^r`, though some bases below `2a + 1` keep it:

| `G` | `b` | `e_r` at `r = 1, 2, ...` |
| --- | --- | --- |
| `{0, 1, 2}` | 3, 4 | `3, 2, 2, 2, 2, 2, 2, 2` |
| `{0, 1, 2}` | 5 | `3`, then `2^r` to `r = 8` |
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

**Table 3.** Every base `a < b <= 2a + 1` of five generators, each row to the last order with at most 17000 sensors and 8 million lags: up to 16384 sensors, the `u^r` rows at `b = 2a + 1` being Theorem A. Domain: those rows. Generator: `lab/py/fractal-array-fragility`, verb `dial`, which asserts the coarray of Proposition 7.1 on every row and checks its fast test, a lag of weight `1` killing both its ends and one of weight `2` the sensor its pairs share, against the definition on every array up to 250 sensors.

Theorem A's proof uses the uniqueness of the balanced expansion and nothing else, and below `2a + 1` that uniqueness is gone, so the rows of Table 3 are data, not theorem. What a compressed base buys shows in the first row of Table 3. `{0, 1, 2}` at base `4` has `N = 3^r` sensors and `(4^(r+1) - 1)/3` lags without a hole, by Proposition 7.1, and only `2` essential sensors at every order from `2` to `8`, by the verb `dial`, the crate test `the_compressed_base_reads_as_printed` pinning orders `1` to `6`; they are the two ends `0` and `2 (4^r - 1)/3`, by `mrly_math::arrays::essential`, and the array has a product beampattern and a self-similar layout; at `r = 8` that is 6561 sensors and 87381 lags. Its fragility `2/N` is the least any array of `N >= 4` sensors can have, which Cohen and Eldar note after their Definition 11 for the uniform line and the robust minimum redundancy array. It is not a record: the lags grow like `N^(log 4 / log 3)`, about `N^1.262`, where nested arrays and minimum redundancy arrays reach order `N^2`, as Cohen and Eldar recall in their Section III.A, and those are maximally economic, the nested arrays with `N_2 >= 2`, by Theorem 1 of [Liu and Vaidyanathan 2019](https://doi.org/10.1109/TSP.2019.2912877), Part II; and the robust minimum redundancy arrays, found by integer programming, already have the least fragility `2/N` with order `N^2` lags, as the composite Singer arrays of [Liu and Vaidyanathan](https://doi.org/10.1109/ICASSP.2019.8683563) recall, and those arrays, built by a recursion on Singer difference sets, are hole-free with fragility `2/N` by their Theorem 2 and have order `N^2` lags for large parameters by their Remark 2.

## Open problems

This paper settles the fragility of fractal arrays at the base `2a + 1` and leaves four things open. First, the compressed base. **Conjecture:** for the five generators of Table 3 at every base `a < b <= 2a + 1`, from `r = 2` or `r = 3` on, the count obeys `e_(r+1) = lambda e_r + c` with an integer `lambda` between `1` and `u`, in three regimes, bounded at `lambda = 1`, geometric at `1 < lambda < u`, and the full `u^r`, which starts below `2a + 1` at `b = 9`, `14`, `14` and `20` for the spans `6`, `7`, `9` and `11`; on the printed terms `{0, 1, 4, 6}` at base `8` reads `e_r = 7 * 2^(r-1) - 3` and `{0, 1, 2, 3, 7, 11}` at base `15` reads `5 * 2^r - 4`. The evidence is every row of Table 3, and its failure modes are plain: the verb `dial` fits `lambda` and `c` on two steps and checks the law on only one to four further steps a row, the rows stop at 17000 sensors, and nothing here explains why `lambda` should be an integer, why it should not exceed `u`, or where the thresholds fall as functions of `G`. A lag at base `b < 2a + 1` is read by a finite carry automaton on digit pairs, so `e_r(G, b)` should have a rational generating function, which would settle the conjecture for each `G` and `b`; neither the automaton nor the function is built here. Second, the joint losses: the sets of `k` sensors whose removal together opens a hole while no smaller part does, the `k`-essential Sperner family of Liu and Vaidyanathan, are known for maximally economic arrays by their Part II, Lemma 1, so for generators with (C1) through Theorem 5; for the others a lag of `F_r` loses all its pairs to `k` removals only if its weight is at most `2k`, which is the product identity again, and the family should be a product of per-digit families at large enough order, which is not attempted here. Third, the fractal arrays of higher-order cumulants of Yang and coauthors, whose coarray is a sum of `q` differences and whose Proposition 2 is again an upper bound; the product identity has an analogue there, and its consequence for their fragility is not worked out. Fourth, the design question the compressed base raises: the robust minimum redundancy arrays reach the least fragility `2/N` with order `N^2` lags by integer programming, and the composite Singer arrays of Liu and Vaidyanathan reach both, for large parameters, by a recursion built on Singer difference sets; whether a closed form reaches both at every size is not settled here.

## Reproducibility

One study and one crate module print every number of this paper beyond the closed forms. The study `lab/py/fractal-array-fragility` runs as `uv run python research/arrays/lab/py/fractal-array-fragility/fragility.py exact` and the same with `dial` from the repository root, one core, about 21 and 20 seconds, needing numpy, taking no input and writing no file; `exact` prints Fact 4.1 and Fact 5.3, the list of the 21 loose generators and the rows `(3, 3), (4, 9), (8, 27)` of `{0, 1, 2}`, `(2, 4), (4, 16), (8, 64)` of `{0, 1, 2, 3}` and `(4, 4), (16, 16), (64, 64)` of `{0, 1, 4, 6}` at orders `1` to `3`, and `dial` prints Table 3 with its fits; its README names the witness of every line. The crate `mrly`, module `num::arrays`, computes the rest: `fractal` builds `F_r` at any base, `weights` and `holes` read the coarray, `paired` gives `U(G)`, `essential` runs the two-pair test of Lemma 2.5 and `lost` lists the lags a loss opens, `law` returns `(u^r, L^r)`. The weights of Sections 5 and 6, the essential sets of Examples 6.1 and 6.2, the opened lags `+-23`, `+-24`, `+-8` and `+-12`, Table 1 and Table 2 are those functions applied to the arrays named, a few milliseconds each; `cargo test -p mrly-math --lib arrays` runs the crate's tests: `the_fast_test_is_the_definition` checks the two-pair test and `lost` against remove-and-recompute, `theorem_a_counts_u_to_the_level` asserts 81, 243 and 729 essential sensors for `{0, 1, 2, 4}` at orders `4` to `6`, 4 of 9 and 8 of 27 for `{0, 1, 2}`, 4 of 16 for `{0, 1, 2, 3}` and 64 of 64 for `{0, 1, 4, 6}`, `cohen_and_eldar_table_one_reads_as_printed` asserts the weight rows, the sets `U` and the six counts of Table 1, and `the_compressed_base_reads_as_printed` the opening terms of four rows of Table 3. The note [arrays](../notes/arrays.md) carries the same results with their claims in [sparse arrays](../claims/sparse-arrays.md). The sparse arrays demo runs the same functions in the browser on any generator, base and order. The figure is `figures/paper-fractal-array-fragility.ts`, drawn in the browser in under a tenth of a second a theme: it builds `F_1` to `F_4` of `{0, 1, 2, 4}` at base 9 with `mrly_math::arrays::fractal`, takes the essential flags from `mrly_math::arrays::essential`, and asserts 4, 9, 27 and 81 essential sensors and, from order `2` on, that the flagged sensors are exactly those with every digit in `U(G) = {0, 1, 4}`.

## References

- Cohen and Eldar 2020, Sparse Array Design via Fractal Geometries, IEEE Trans. Signal Process. 68, doi 10.1109/TSP.2020.3016772; read at source in arXiv:2001.01217v2: the recursion (5), Theorems 1 and 4, Definitions 9 to 11, Lemma 1, condition (C1), Theorems 5 and 6 with their proofs, the remark on translation factors in Section III.B, the remark on order `N^2` coarrays in Section III.A, the least fragility `2/N` after Definition 11, and Section V with problem (P1), requirement (R1) and Table I. [arxiv.org/abs/2001.01217](https://arxiv.org/abs/2001.01217)
- Yang, Shen, Liu, Eldar and Cui 2023, High-Order Cumulants Based Sparse Array Design Via Fractal Geometries, Part II: Robustness and Mutual Coupling, IEEE Trans. Signal Process. 71, 343-357, doi 10.1109/TSP.2023.3244667; read at source: Section IV.A, Lemma 1 and Proposition 1 with its proof, the statement of Proposition 2, the generators (26) and (27) of Section IV.C and Table I. [doi.org/10.1109/TSP.2023.3244667](https://doi.org/10.1109/TSP.2023.3244667)
- Liu and Vaidyanathan 2019, Robustness of Difference Coarrays of Sparse Arrays to Sensor Failures, Part II: Array Geometries, IEEE Trans. Signal Process. 67, no. 12, 3227-3242, doi 10.1109/TSP.2019.2912877; read at source: Definitions 8 and 9, Lemma 1 on the `k`-essential Sperner family of maximally economic arrays, and Theorem 1. [doi.org/10.1109/TSP.2019.2912877](https://doi.org/10.1109/TSP.2019.2912877)
- Liu and Vaidyanathan 2019, Composite Singer Arrays with Hole-free Coarrays and Enhanced Robustness, ICASSP 2019, 4120-4124, doi 10.1109/ICASSP.2019.8683563; read at source: Definition 6, Corollary 1, Proposition 1 with Property 1, Definition 8, Theorem 2 and Remark 2. [doi.org/10.1109/ICASSP.2019.8683563](https://doi.org/10.1109/ICASSP.2019.8683563)
- Goel, Aggarwal and Kar 2022, A Highly Robust Sparse Fractal Array, arXiv:2212.00341; read at the abstract. [arxiv.org/abs/2212.00341](https://arxiv.org/abs/2212.00341)
