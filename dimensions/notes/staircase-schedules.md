---
title: Staircase schedules and the unequal split
lead: Builds the cheapest non-stationary object, a staircase schedule, and works the unequal split: its complex dimensions, cell count, detectors, and what survives.
figure: research-dimensions
slug: dimensions
---

## Staircase schedules, the cheapest non-stationary object

- Instead of a constant word, stack `carpet_3`, then `magic(3,5)`, then `magic(3,5,7)`, and so on; Kronecker associativity flattens that to the staircase word `3 | 3,5 | 3,5,7 | 3,5,7,9 | ...`.
- Letter `base_j` occurs `n - j + 1` times in the first `n` blocks, so the controls are immediate and non-negotiable: `side = prod_j base_j^(n-j+1)`, `fill = prod_j fill_j^(n-j+1)`, and `dimension_n = Sum_j (n-j+1)*ln(fill_j) / Sum_j (n-j+1)*ln(base_j)`.
- The staircase word is aperiodic and not eventually periodic, so block reduction does not apply to it. It is the cheapest concrete non-stationary schedule available.
- The weights `(n-j+1)` are a Cesaro profile - the earliest letter carries weight `n`, the newest carries 1 - so if the letter dimensions converge the whole dimension converges to their limit; the interesting regime is letter dimensions that oscillate.
- The generalisation is what makes this a programme rather than an example: a staircase is one weight profile, any letter-multiplicity schedule is another, and the question "which dimension functions are realisable by a schedule and which are not" quantifies over all schedules and is native to the construction.
- **Caveat that travels with every mixed number.** If each factor is rendered at its native base, `side_i = base_i`, the filled points have the mixed-radix form `x = a_1*(base_2...base_level) + ... + a_level` with `a_i` in `F_i`, and that is the correct arithmetic object. A factor rendered at a side unrelated to its residue base is still a valid tile product, but it is not a mixed-radix digit construction and inherits no digit theorem for free.
- `cuts/lab/py/slice-ladder-controls` prints the five staircase dimensions - `1.892789261`, `1.892315261`, `1.893034267`, `1.894190425`, `1.895495742` at `n = 1..5` - assuming the carpet at one base has fill `base^2 - ((base-1)/2)^2`; the run states that assumption before any number. (Verified under that assumption; confirm the definition against [the core page](../../core/notes/core.md) before quoting any number from it.)

## The unequal split

Objects. The *patch* is the attractor `K` of six similarities of the plane, `x -> x/2` and `x -> x/3 + t` for the five translations `t` in `{(2/3, 0), (2/3, 1/3), (2/3, 2/3), (0, 2/3), (1/3, 2/3)}`. Its letter splits one cell into children of two sizes, drawn on the 6-grid: the half covers a `3 x 3` block in a corner, each third a `2 x 2` block, and an L of seven cells stays empty. The same letter with `k` thirds beside the half has the Moran function `f(s) = 1 - 2^(-s) - k 3^(-s)`, and the patch is `k = 5`. A *cell* is the image of the unit square under a word of the maps; a word with `a` halves and `b` thirds draws a cell of side `2^-a 3^-b`, and `N(r)` counts the cells of side at least `r`, the unit square included. `D` is the patch's dimension, `P = 2^(-D)` and `Q = k 3^(-D)`, so `P + Q = 1`, and the natural measure gives the half the weight `P` and each third `3^(-D)`. Every number of this section is printed by `lab/py/unequal-split`, verbs `patch`, `poles`, `count` and `ripple`, or rounded from a printed one.

### The patch

**Proved.** The six open images of the open unit square are disjoint open squares inside it, so the open set condition holds, and by [Hutchinson 1981](../refs.md), 5.3(1), the dimension of `K` is the root of the Moran equation `2^(-s) + 5*3^(-s) = 1`, `D = 1.778602507`. Since `2^m = 3^n` has no solution in positive integers, `log 2/log 3` is irrational and the patch is nonlattice; by block reduction above no design and no periodic composition is. With the half in a corner, five thirds is the most: a third whose interior misses the half `[0, 1/2]^2` has its lower left corner at `x >= 1/2` or at `y >= 1/2`; those with `x >= 1/2` lie in a strip of width `1/2`, overlap pairwise in `x` and so stack in `y`, at most three; the others have `x < 1/2 <= y`, lie in a strip of height `1/2` and so sit side by side with left edges below `1/2`, at most two. So `k = 1..5` are planar patches, while `k >= 7` admits no open set condition in the plane at all, the squared ratios summing to `1/4 + k/9 > 1`, which disjoint images of one bounded open set forbid by area. As a string the equation stands at every `k`, and `k = 1` is the 2-3 nonlattice equation of [Lapidus and van Frankenhuijsen 2003](../refs.md), Section 2.2.4.

**Proved.** The level render, the union of the cells of the words of length `level`, is a picture on side `6^level`: a cell of side `2^-a 3^-b` with `a + b = level` covers `2^b 3^a` grid cells per axis at a grid corner, so the fill is `sum_a C(level, a) 9^a (4k)^b = (9 + 4k)^level`, `29^level` here. That is the fill of a base-6 design with 29 filled digits, and the dimension the tree would read off it is `log 29/log 6 = 1.879323585`, not `D`; the picture is not that design either, since level 2 differs from the Kronecker square of the level-1 tile in 336 of its 1296 cells. A level keeps a cell of side `2^(-level)` beside one of side `3^(-level)`, so no level is a cover at one scale. **Verified** (`patch`): the open set condition as exact rationals, the fill to level 4 with no cell covered twice, the 336, and the census of cells to word length 7, every size `2^-a 3^-b` met exactly `C(a+b, a) 5^b` times.

### Its complex dimensions

The sides of all the cells sum to the scaling zeta function, `sum over cells of side^s = sum_n (2^(-s) + k 3^(-s))^n = 1/f(s)` for `Re s > D`, the analogue of `1/(1 - fill*base^(-s))` above; its poles, the zeros of `f`, are the complex dimensions.

**Proved.** Every zero lies in the strip `D_l <= Re s <= D`, where `D_l` is the real root of `k 3^(-s) = 1 + 2^(-s)`: right of `D` the two terms have modulus summing below 1, and left of `D_l` the term `k 3^(-s)` outweighs `1 + abs(2^(-s))`. On the line `Re s = D` the only zero is `D` itself, since equality in the triangle inequality needs `t ln 2` and `t ln 3` both in `2 pi Z`, and an irrational ratio forces `t = 0`. A lattice string, all ratios powers of one `r`, makes `f` a polynomial in `r^s`, so its zeros lie on finitely many vertical lines, each an arithmetic progression; the patch's real parts take infinitely many values, by the approach to the line below. The strip is Theorem 2.5 of [Lapidus and van Frankenhuijsen 2003](../refs.md), read at source, with the same `D_l`.

**Verified** (`poles`). For every `k = 1..27` the box `Re [D_l - 1/4, D + 1/4]`, `Im [-60, 60]` holds exactly 21 zeros: the winding number of `f` around its boundary is 21, the smallest `abs(f)` on the boundary exceeds a Lipschitz bound on the change of `f` over one step by a factor of at least 161 where 1 suffices, and Newton from the roots of the lattice approximant `1 - x^41 - k x^65` (65/41 a convergent of `log2 3`) finds all 21, so the list is complete; 21 is the count `(ln 3/pi) 60 + O(1)` of the same theorem. The lowest real part comes within `0.0116` of `D_l` at every `k`. The patch's first complex dimensions, the upper half of the box:

| `Re` | `Im` |
|---|---|
| 1.778603 | 0 |
| 1.241037 | 5.976088 |
| 1.455984 | 11.106537 |
| 1.719303 | 17.347538 |
| 1.124789 | 22.929912 |
| 1.657088 | 28.338580 |
| 1.552742 | 34.630309 |
| 1.171657 | 39.856618 |
| 1.766880 | 45.665889 |
| 1.330529 | 51.782362 |
| 1.355854 | 56.873751 |

Eleven distinct real parts, and some imaginary part misses the multiples of `2 pi/ln 2`, `2 pi/ln 3` and `2 pi/ln 6` by at least 0.470, 0.058 and 0.461 of a step: no vertical line, no progression. For the 2-3 nonlattice equation, `k = 1`, the root that Section 3.1 of the same paper prints, `.7675115443 + 45.55415979 i`, lies `4.0e-8` from the root of its lattice approximant `1 - 2^(-s) - 2^(-485 s/306)`, the polynomial it was computed from, and `7.6e-5` from the true complex dimension `0.7674996132 + 45.55423466 i`, the approximation error of that lattice string at that height.

**Proved.** The complex dimensions come arbitrarily close to the line `Re s = D`, at a rate set by how well `log2 3` is approximated; this is Theorem 4.3, equation (4-9), of the same paper, with multiplicities 1 and `k`, restated here by the implicit function theorem. Take `t` with `t ln 2` in `2 pi Z` and `theta = t ln 3` reduced mod `2 pi`. Then `f(D + it + z) = 1 - P e^(-z ln 2) - Q e^(-i theta) e^(-z ln 3)` depends on `t` only through `theta`, is analytic in `(z, theta)`, vanishes at `(0, 0)` with `z`-derivative `f'(D) = P ln 2 + Q ln 3 > 0`, so the implicit function theorem gives an analytic zero `z(theta)`, and expanding to second order,

```
z(theta) = -i Q theta/f'(D) - P Q (ln 2)^2 theta^2/(2 f'(D)^3) + O(theta^3).
```

At `t = 2 pi q/ln 2` for a convergent `p/q` of `log2 3`, `theta = 2 pi (q log2 3 - p)` tends to 0, so `D - Re w = 2 pi^2 P Q (ln 2)^2 (q log2 3 - p)^2/f'(D)^3 + O(theta^3)`. **Verified** (`poles`, Newton at 60 digits in PARI, `abs(f(w))` below `1e-69` at every root):

| `q` | `Im w` | `D - Re w` | law | ratio |
|---|---|---|---|---|
| 2 | 17.347538 | `5.929932e-02` | `6.000528e-02` | 0.988235 |
| 5 | 45.665889 | `1.172287e-02` | `1.174806e-02` | 0.997856 |
| 12 | 108.687857 | `7.941580e-04` | `7.942712e-04` | 0.999857 |
| 41 | 371.728633 | `5.682884e-04` | `5.683464e-04` | 0.999898 |
| 53 | 480.416496 | `1.885987e-05` | `1.885993e-05` | 0.999997 |
| 306 | 2773.811103 | `4.519898e-06` | `4.519902e-06` | 0.999999 |
| 665 | 6028.038703 | `8.242791e-09` | `8.242791e-09` | 1.000000 |
| 15601 | 141418.701264 | `1.431886e-09` | `1.431886e-09` | 1.000000 |
| 31867 | 288865.441232 | `2.282902e-10` | `2.282902e-10` | 1.000000 |
| 79335 | 719149.583728 | `5.809059e-11` | `5.809059e-11` | 1.000000 |
| 111202 | 1008015.024959 | `5.606353e-11` | `5.606353e-11` | 1.000000 |
| 190537 | 1727164.608687 | `1.799884e-14` | `1.799884e-14` | 1.000000 |
| 10590737 | 96002068.502749 | `1.183140e-14` | `1.183140e-14` | 1.000000 |
| 10781274 | 97729233.111436 | `6.445263e-16` | `6.445263e-16` | 1.000000 |
| 53715833 | 486919000.948492 | `5.215918e-17` | `5.215918e-17` | 1.000000 |

The partial quotient 55 of `log2 3` puts a complex dimension within `1.8e-14` of the line at height `1.7e6`, and the last row within `5.3e-17` at height `4.9e8`. The lattice designs above have their whole column on the line; the patch has one point there and a sequence approaching it.

### The count of cells

**Proved.** The cells of side at least `r` are the unit square, when `r <= 1`, together with the cells of side at least `2r` inside the half and of side at least `3r` inside each third, so `N(r) = 1 + N(2r) + k N(3r)` with `N(r) = 0` for `r > 1`, and `N(r) = sum C(a+b, a) k^b` over `2^-a 3^-b >= r`. In `U = ln(1/r)` the Laplace transform of `N` is `1/(s f(s))`: its poles are the complex dimensions and `s = 0`, where the residue is `1/f(0) = -1/k`. The box count of [the oscillation section](string-to-composition.md#the-oscillation-in-the-box-count) rebuilt from cells, the cells of side at most `r` whose parent is larger, numbers `1 + k L` with `L` the count of cells larger than `r`, every such cell having been split into `1 + k`, so it carries the same reading.

**Proved.** `N(r) r^D` tends to `C = 1/(D f'(D))`, `0.573459971` at `k = 5`. The limit exists by Theorem 1 of [Lalley 1989](../refs.md): `N(e^(-T))` is his count `N(T, x)`, at every `x`, for the function `-log` of the first map's ratio on the full shift on six symbols, nonlattice because its periodic sums `a ln 2 + b ln 3` lie in no discrete subgroup. The value is forced: `F(U) = N(e^(-U)) e^(-DU)` has Laplace transform `1/((s + D) f(s + D))`, and if `F` tends to `C` then `s` times the transform tends to `C` as `s -> 0+`, which reads `1/(D f'(D))`. The carpet, the lattice control, does the opposite: `N(r) = (8^(m+1) - 1)/7` with `m` the integer part of `log_3(1/r)`, so `N(r) r^D` is a fixed non-constant `ln 3`-periodic function of `ln r` less `r^D/7`, and never converges.

**Verified** (`count`). The count as exact integers out to `r = e^(-300)`: 59448 sizes, sorted by exact comparison of the integers `2^a 3^b`, `N` reaching 232 digits, and the identity `N(r) = 1 + N(2r) + 5 N(3r)` exact at every one. Over windows of length 10 in `U`, with `F/C` sampled at 200001 points:

| window | `min F/C` | `max F/C` | mean `- 1` | swing | swing `sqrt(U)` | carpet swing |
|---|---|---|---|---|---|---|
| `[10, 20]` | 0.794332 | 1.264454 | -0.000375 | 0.470122 | 1.486657 | 2.062470 |
| `[20, 30]` | 0.836870 | 1.182388 | -0.000433 | 0.345518 | 1.545201 | 2.070132 |
| `[40, 50]` | 0.881974 | 1.121301 | -0.000075 | 0.239327 | 1.513635 | 2.081449 |
| `[80, 90]` | 0.916164 | 1.086485 | -0.000103 | 0.170321 | 1.523396 | 2.093786 |
| `[160, 170]` | 0.941725 | 1.062159 | +0.000036 | 0.120434 | 1.523381 | 2.089572 |
| `[290, 300]` | 0.956827 | 1.046449 | +0.000000 | 0.089622 | 1.526205 | 2.064570 |

Over these six windows the mean sits on `C` to `4.4e-4`, the swing falls without a floor, and the carpet's swing, relative to its own mean, stays between 2.06 and 2.10. **Conjecture:** the swing decays like `1.52 U^(-1/2)`, its product with `sqrt(U)` staying in `[1.486, 1.546]` on the six printed windows, read only at starts that are multiples of 10 and not at other starts; the jump of `F/C` at the size `2^-a 3^-b` is the binomial probability `C(a+b, a) P^a Q^b` over `C`, of order `U^(-1/2)`, but no bound is proved.

**Verified** (`count`). The detector of [the oscillation section](string-to-composition.md#the-oscillation-in-the-box-count), run on `g(u) = ln N(e^(-u)) - D u` over `u` in `[50, 300]` at step `0.002`. Folded into 40 bins, `g` keeps `0.003`, `0.005` and `0.004` of its variance at `ln 2`, `ln 3` and `ln 6`, against the carpet's `0.999` at `ln 3`. The Blackman periodogram does not go flat; it reads the complex dimensions. Its ten highest peaks between 2 and 600 sit within `0.001` of the imaginary parts of ten zeros of `f`, every one with `D - Re w` below `0.019`, and their heights match `2 e^(-(D - Re w) u)/(abs(w f'(w)) C)`, averaged under the same taper, the residue of `1/(s f(s))` at `w`, to `1.3%`:

| peak | height | complex dimension | `D - Re w` | residue height |
|---|---|---|---|---|
| 45.666 | `1.111e-02` | `1.766880 + 45.666 i` | `1.17e-02` | `1.112e-02` |
| 63.021 | `2.881e-03` | `1.760015 + 63.021 i` | `1.86e-02` | `2.845e-03` |
| 108.688 | `2.848e-02` | `1.777808 + 108.688 i` | `7.94e-04` | `2.849e-02` |
| 154.354 | `7.700e-03` | `1.772177 + 154.353 i` | `6.43e-03` | `7.723e-03` |
| 217.376 | `9.444e-03` | `1.775427 + 217.376 i` | `3.18e-03` | `9.456e-03` |
| 263.040 | `8.456e-03` | `1.775897 + 263.041 i` | `2.71e-03` | `8.467e-03` |
| 326.064 | `3.242e-03` | `1.771463 + 326.063 i` | `7.14e-03` | `3.251e-03` |
| 371.728 | `8.661e-03` | `1.778034 + 371.729 i` | `5.68e-04` | `8.664e-03` |
| 480.416 | `7.379e-03` | `1.778584 + 480.416 i` | `1.89e-05` | `7.380e-03` |
| 589.104 | `5.019e-03` | `1.777545 + 589.104 i` | `1.06e-03` | `5.021e-03` |

The peaks sit near `2 pi q/ln 2` for `q = 5, 7, 12, 17, 24, 29, 36, 41, 53, 65`, the `q` that bring `q log2 3` near an integer, each moved off it by the `-Q theta/f'(D)` of the expansion above; the carpet's six highest peaks sit at 1, 2, 3, 4, 5 and 6 times `2 pi/ln 3` to the printed three decimals. The expectation is the zeros and their residues, computed from `f` alone and never fitted to `g`. The fluctuation of `N(r) r^D` is not periodic: it is a sum of decaying waves at incommensurable frequencies over a sawtooth that shrinks, and it converges.

### What the detectors see

The count of cells sees the unequal split; the spin ripple of [the spin page](../../spin/notes/spin.md), which reads the mass `M(r) = mu(B(c, r))` of the natural measure `mu` within radius `r` of a point `c`, does not, one point at a time. `ripple` computes `M` from the six maps, enclosing it between the cells inside the ball and those not yet decided at relative size `1e-4`.

**Proved.** Let `c` be the fixed point of one map, of ratio `rho`, and suppose no other child comes within `r_0` of `c`. Then `M(r) = rho^D M(r/rho)` for `r < r_0`: in `mu = sum_i w_i mu(S_i^(-1) .)`, with `S_i` the six maps and `w_i` their weights, only the term of the map fixing `c` meets the ball, and that map pulls `B(c, r)` back to `B(c, r/rho)`. So `ln M(r) - D ln r` is exactly periodic in `ln r`, with period `ln(1/rho)`. At the corner `(0, 0)`, fixed by the half, `r_0 = 2/3` and the period is `ln 2`; at `(1, 0)` and `(1, 1)`, fixed by thirds, `r_0 = 1/3` and the period is `ln 3`. The ripple keeps its level periodicity about every such point.

**Verified** (`ripple`). Over 120 radii a 24th of a period apart, down from `0.999 r_0`, the enclosures of `M(r)` and `rho^D M(r/rho)` meet at all 96 shifted radii about each of the three corners. Folded into 24 bins, each ripple keeps all its variance at its own period and almost none at the other:

| centre | period | swing | enclosure | fold `ln 2` | fold `ln 3` | fold `ln 6` |
|---|---|---|---|---|---|---|
| `(0, 0)` | `ln 2` | 0.23906 | `8.2e-5` | 1.000 | 0.052 | 0.079 |
| `(1, 0)` | `ln 3` | 0.15388 | `4.2e-5` | 0.043 | 1.000 | 0.064 |
| `(1, 1)` | `ln 3` | 0.17843 | `6.7e-5` | 0.059 | 1.000 | 0.036 |

**Proved** (computer-assisted, `ripple`). Two ripples see what one cannot. A ripple `G` with periods `ln 2` and `ln 3` is constant: `M` is monotone and right-continuous, so `G` is right-continuous and continuous off a countable set; its periods contain the group `Z ln 2 + Z ln 3`, dense because the ratio is irrational, and invariance under a dense group carries the value at one continuity point to every other, then right-continuity to all points. The three ripples are not constant, their swings exceeding their double-precision enclosures at least 2661-fold, so the ripple about `(0, 0)` has period `ln 2` and not `ln 3`, the one about `(1, 0)` period `ln 3` and not `ln 2`, and the two share no period at all, their groups of periods being cyclic and `ln 2/ln 3` irrational. On a design every map of every word has ratio a power of `1/base`, so every identity of this kind has period a multiple of `ln base`, and any two such ripples share a period. Two fixed points with incommensurable ripples are a reading no design can give.

**Verified** against the source. The measure and the tube follow the same split. [Rapaport 2022](../refs.md), Corollary 1.6: a self-similar measure of an affinely irreducible system on `R^d` that is not Rajchman forces every ratio to be `theta^(-n_i)` for one algebraic integer `theta > 1`; here `2 = theta^(n_1)` and `3 = theta^(n_2)` would give `2^(n_2) = 3^(n_1)`. A homothety fixes a line only when the line passes through its fixed point, and the six fixed points `(0, 0)`, `(1, 0)`, `(1, 1/2)`, `(1, 1)`, `(0, 1)`, `(1/2, 1)` lie on no line, so the system is affinely irreducible and every self-similar measure on the patch with positive weights, the natural one included, has Fourier transform tending to `0`. On the line the same holds for the string of the 2-3 nonlattice equation by Theorem 1.2 of [Li and Sahlsten 2019](../refs.md), which needs only some `log r_i/log r_j` irrational and no separation. And by [Gatzouras 2000](../refs.md), Theorem 2.3(i) with Theorem 2.4, a self-similar set under the open set condition whose log-ratios lie in no `lambda Z` is Minkowski measurable, with no further hypothesis: the patch is Minkowski measurable, where the carpet, proved above, is not.

### What survives and what dies

The patch sits beside the designs as their nearest neighbour outside the family: one letter, substituted into itself, with the Kronecker product replaced by the plain composition of maps.

- Survives, the dimension: the weighted child matrix of the letter is the `1 x 1` matrix `2^(-s) + k 3^(-s)`, equal to 1 exactly at `s = D`, Hutchinson's equation above; `log(fill)/log(base)` is the case of one size.
- Survives, the box dimension: the count of cells is `C r^(-D) (1 + o(1))`, so `ln N(r)/ln(1/r)` reads `D`, as the level count reads `log(fill)/log(base)` on a design.
- Survives, the string's zeta: `1/f(s)`, the Moran function in the denominator where `1 - fill*base^(-s)` stood; the complex dimensions are its poles, `D` is simple and alone on its line, and the residues drive the count term by term.
- Survives, the corner lemma: about every fixed point with a clear window, with that map's ratio in place of `1/base`.
- Dies, the level: a level mixes the sides `2^(-level)` and `3^(-level)`, its render reads `log 29/log 6` and not `D`, and the count by size replaces the count by level.
- Dies, the code: the render is no Kronecker power of its tile, so nothing that takes a design's code, the spin census, the design zeta or the digit transform, has an input.
- Dies, the carry automaton: the cells ordered by size interleave words of every length along `a ln 2 + b ln 3`, and no base reads their addresses digit by digit.
- Dies, the lattice: the count converges instead of oscillating, the tube has a Minkowski limit, and the measures with positive weights are Rajchman.
- The one reading that keeps a period is local: the ripple about a fixed point. Its period changes from point to point, and that change is the nonlattice signature in the spin's own language.
