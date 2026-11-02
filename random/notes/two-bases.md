---
title: The moire limit
lead: The moire XOR of Thue-Morse and the ternary parity has correlations that tend to zero. Thue-Morse correlations decay in binary runs, abs C_t(h) <= 1.0093 . 0.7177^r with r the number of runs, the ternary layer decays in balanced-ternary weight, and by the S-unit theorem only finitely many integers are sparse in both bases, so two spectral measures that are not Rajchman, at the bases 2 and 3, convolve to one that is. The S-unit theorem is ineffective, but the digit-change theorem of Bugeaud, Cipu and Mignotte gives an effective rate, abs C_z(h) <= 1.0093 . 3^(-G(h)/4) with G(h) = log log h/(log log log h + C) - 1, whose constant C is effectively computable and not explicit.
figure: research-random-two-bases
---

The [design ladder](design-ladder.md) left one row open. The moire XOR `z = t . t3`, Thue-Morse times the parity of the ternary digits equal to `2`, has shift correlation `C_z(h) = C_t(h) C3(h)` at every lag; each factor stays away from `0` along the powers of its own base, and the product was conjectured to tend to `0`. This chapter proves it in three steps. The binary layer decays in runs: `abs C_t(h) <= 1.0093 . 0.7177^(r_2(h))`, `r_2(h)` the number of runs of the binary word of `h`, because one quadratic form is contracted by every run of digits; the best constant per run lies between `0.7118` and `0.7177`, and decay in the number of ones is false, since `C_t(2^k - 1) -> 1/9`. The ternary layer decays in balanced-ternary weight, as the ladder proved. An integer with few binary runs and few balanced-ternary digits is a short signed sum of powers of `2` equal to a short signed sum of powers of `3`, and the S-unit theorem of van der Poorten and Schlickewei and of Evertse allows that only finitely often. So `C_z(h) -> 0`, and the spectral measure of `z`, a convolution of two measures that are not Rajchman, is Rajchman. The S-unit theorem rests on the subspace theorem and is ineffective; the rate comes from a theorem of Bugeaud, Cipu and Mignotte on digit changes in two bases, which follows Stewart's Baker argument. A signed sum of `l` powers of `b` changes digit at most `2l - 1` times, so `abs C_z(h) <= 1.0093 . 3^(-G(h)/4)` for `h > 25`, `G(h) = log log h/(log log log h + C) - 1` in natural logarithms (the source leaves the base unstated), with `C` effectively computable from the bases and not explicit; no number can be read off the rate until `C` is computed. The census finds `15196` integers below `10^40` with at most `6` binary runs and at most `6` balanced-ternary digits; the largest is `283462598656`, and nothing lies between it and the cut.

The generator is the study `lab/rs/two-bases`, verbs `lemma`, `decay`, `windows`, `changes` and `sparse`. The figure lays the lags `h < 4096` on a `64 x 64` grid, row by row: blue where `abs C_t(h) >= 1/16`, orange where `abs C3(h) >= 1/27`, and the foreground ink where both hold, `63` cells, the last at `h = 3888`.

## The objects

Thue-Morse is `t(n) = (-1)^(s_2(n))`, `s_2` the binary digit sum, and `t3(n) = (-1)^(number of ternary digits of n equal to 2)`. Both are signed designs in the sense of the ladder, so their correlations `C_x(h) = lim (1/N) sum_(n < N) x(n) x(n + h)` exist at every lag, are even in `h` and equal `1` at `h = 0`. The ladder proved the product rule `C_z(h) = C_t(h) C3(h)` for a signed design at base `2` times one at the coprime base `3`, by the Chinese remainder theorem on periodic approximations.

The block recursion of the ladder, `C(bq + s) = (al_s C(q) + al_(b-s) C(q + 1))/b` with `al_s = sum_(r < b - s) v(r) v(r + s)`, reads at the two tiles:

- Thue-Morse, `v = (1, -1)`, `al_0 = 2`, `al_1 = -1`: `C_t(2h) = C_t(h)` and `C_t(2h + 1) = -(C_t(h) + C_t(h + 1))/2`, so `C_t(1) = -1/3` from `C_t(1) = -(1 + C_t(1))/2`.
- The ternary parity, `v = (1, 1, -1)`, `al_0 = 3`, `al_1 = 0`, `al_2 = -1`: `C3(3h) = C3(h)`, `C3(3h + 1) = -C3(h + 1)/3` and `C3(3h + 2) = -C3(h)/3`, so `C3(1) = 0` and `C3(2) = -1/3`.

Three counts of sparsity:

- `r_2(h)`, for `h >= 1`, is the number of maximal blocks of equal digits in the binary word of `h`, leading zeros excluded: `r_2(2^k - 1) = 1`, `r_2(2^k) = 2`, `r_2(1049760) = 8`.
- `w_3(h)` is the number of nonzero digits of `h` in balanced ternary, digits `-1, 0, 1`.
- `l_b(h)` is the least number of terms `+-b^a`, `a >= 0`, with sum `h`, repetitions allowed.

**`l_2(h) <= r_2(h) + 1` and `l_3(h) <= w_3(h)`. Proved.** A run of ones on the bits `i` to `j - 1` is `2^j - 2^i`, and a word with `r` runs has at most `ceil(r/2)` runs of ones, so `h` is a signed sum of at most `2 ceil(r/2) <= r + 1` powers of `2`; balanced ternary writes `h` with `w_3(h)` signed powers of `3`. Neither bound reverses for free, and only these two directions are used.

**`abs C3(h) <= 3^(-w_3(h)/2)` at every `h >= 1`. Proved in the ladder**, by induction on the recursion above: a nonzero digit costs a factor `1/3` and moves the weight by at most two. `decay` rechecks it at every `h < 2^25` with `0` failures.

## The run lemma

Write `v(h) = (C_t(h), C_t(h + 1))`. The recursion is linear in this pair: `v(2h) = M_0 v(h)` and `v(2h + 1) = M_1 v(h)` with

- `M_0 = [[1, 0], [-1/2, -1/2]]`, eigenvalue `1` on `e_0 = (1, -1/3)` and `-1/2` on `(0, 1)`;
- `M_1 = [[-1/2, -1/2], [0, 1]]`, eigenvalue `1` on `e_1 = (-1/3, 1)` and `-1/2` on `(1, 0)`.

The start is `v(0) = (1, -1/3) = e_0`, fixed by `M_0`, so leading zeros cost nothing, and reading the binary digits of `h` from the top gives `v(h) = M_(h_0) M_(h_1) ... M_(h_(L-1)) e_0`. The two matrices share their eigenvalues on different eigenvectors. Inside a run of `k` equal digits `d`, `M_d^k = P_d + (-1/2)^k (I - P_d)`, `P_d` the projection on `e_d` along the other eigenvector, so a long run forgets everything but its `e_d` component. At a switch the projections meet, `P_1 e_0 = -e_1/3` and `P_0 e_1 = -e_0/3`, so long runs lose a factor `1/3` per switch; short runs keep part of the `-1/2` direction and lose less. Nothing is lost inside a run: `C_t(2h) = C_t(h)` and `C_t(2^k - 1) -> 1/9`. So the decay must be counted in runs.

The swap `J(x, y) = (y, x)` conjugates the two digits, `J M_0 J = M_1`. With `t = (-1/2)^b`, a run of `b` zeros is `M_0^b = [[1, 0], [(t - 1)/3, t]]`, and the run operator is

`A(t) = M_0^b J = [[0, 1], [t, (t - 1)/3]]`.

**If `h` has runs of lengths `b_1, ..., b_r` from the top, then `v(h) = J^(r mod 2) A(t_r) ... A(t_1) e_0` with `t_i = (-1/2)^(b_i)`. Proved.** The top run is ones, and `M_1^(b_1) e_0 = J M_0^(b_1) J e_0 = J A(t_1) e_0`; a run of zeros after it multiplies by `M_0^(b_2)` and absorbs the `J`, giving `A(t_2) A(t_1) e_0`; a run of ones multiplies by `J M_0^(b_3) J` and leaves a `J` in front; induction on `r`. So `C_t(h)` is the first coordinate of `u = A(t_r) ... A(t_1) e_0` when `r` is even and its second when `r` is odd. At `h = 1`, `A(-1/2) e_0 = (-1/3, -1/3)` and its second coordinate is `C_t(1) = -1/3`.

Every `t_i` lies in `[-1/2, 1/4]`: `b = 1` gives `-1/2`, `b = 2` gives `1/4`, and the rest lie between them. One norm contracts the whole interval.

**Let `Q(x, y) = 20x^2 + 20xy + 41y^2`. Then `Q(A(t) u) <= (103/200) Q(u)` for every `u` in `R^2` and every `t` in `[-1/2, 1/4]`. Proved.** Write `(103/200) Q(u) - Q(A(t) u) = a(t) x^2 + B(t) xy + d(t) y^2`. Expanding `A(t) u = (y, tx + (t - 1) y/3)`:

- `a(t) = 103/10 - 41 t^2`;
- `B(t) = 103/10 + (22/3) t - (82/3) t^2`;
- `d(t) = 5807/1800 + (22/9) t - (41/9) t^2`;
- `Delta(t) = 4 a(t) d(t) - B(t)^2 = 30179/1125 - (2266/45) t - (93367/450) t^2`, the quartic and cubic terms cancelling.

The binary form is positive semidefinite when `a >= 0`, `d >= 0` and `Delta >= 0`. All three are concave quadratics in `t`, so on an interval each takes its minimum at an endpoint, and the endpoint values are positive: `a = 1/20` and `619/80`, `d = 173/200` and `1421/400`, `Delta = 133/1000` and `5077/4000` at `t = -1/2` and `t = 1/4`. `lemma` recomputes the three polynomials in exact rationals, checks their concavity and the six values, and counts their real roots on the closed interval by an exact Sturm sequence: none. The form was found by a numerical scan, and nothing in the proof depends on how it was found.

**`abs C_t(h) <= K rho^(r_2(h))` at every `h >= 1`, with `rho = (103/200)^(1/2) = 0.717635...` and `K = (41 . 161/720)^(1/2)/3 = 1.009293...`; rounded safe, `abs C_t(h) <= 1.0093 . 0.7177^(r_2(h))`. Proved.** Let `N(u) = Q(u)^(1/2)`. By the lemma each run operator shrinks `N` by `rho`, so `N(u) <= rho^r N(e_0)`, and `N(e_0)^2 = Q(1, -1/3) = 20 - 20/3 + 41/9 = 161/9`. A coordinate is bounded by the form through its inverse matrix: `Q` has matrix `[[20, 10], [10, 41]]`, determinant `720`, and the diagonal of the inverse is `41/720` and `20/720`, so `abs x <= (41/720)^(1/2) N(u)` and `abs y <= (20/720)^(1/2) N(u)`. That gives `K` for even `r` and `0.704921` for odd `r`. Since `C_t(2h) = C_t(h)`, the bound also holds with the runs of the odd part of `h`, one fewer when `h` is even.

The constant per run is nearly sharp. **No `rho < 0.711872` satisfies `abs C_t(h) <= K' rho^(r_2(h))` for any `K'`. Proved by computation.** For the run pattern `(1, 1, 2)` repeated, the operator over three runs is `A(1/4) A(-1/2) A(-1/2) = [[1/4, -1/4], [-3/16, -1/16]]`, with characteristic polynomial `x^2 - (3/16) x - 1/16` and spectral radius `(3 + sqrt 73)/32 = 0.360750`, whose cube root is `0.711872`. The eigenvectors are `(1, (5 - sqrt 73)/8)` for the larger root and `(1, (5 + sqrt 73)/8)` for the smaller, so `e_0 = (1, -1/3)` has a nonzero component on the first, whose two coordinates are both nonzero, and along `r = 3j` either coordinate of the `j`-th power applied to `e_0` is of exact order `0.360750^j`. `decay` reads the family exactly to `30` runs: at `h = 776125854900` the value is `3.5415e-5`, and the three-run ratios alternate `0.933753` and `0.542694`, as the coordinate read alternates with the parity of `r`, so the root over six runs is `(0.933753 . 0.542694)^(1/2) = 0.711858`. The best constant per run therefore lies in `[0.7118, 0.7177]`.

`decay` checks the theorem directly at every `h < 2^25`: `0` failures, the largest value of `abs C_t(h)/(1.0093 . 0.717636^(r_2(h)))` being `0.9067`, at `h = 90`. The maximum of `abs C_t` at a run count, over the same window, runs `1/3, 1/3, 1/6, 1/6, 1/8, 1/8, 1/16, 1/16, 0.02865, 0.02865, 0.01611, 0.01611` for `r = 1..12`. It repeats in pairs because a trailing zero adds a run and keeps the value, and it falls by about `1/2` every two runs while the window holds the worst words. At `r = 24` it is `2.449e-4`, below the `2.722e-4` the family reaches at `h = 3031741620`, which lies outside the window.

## Runs, not ones

**`C_t(2^k - 1) = (-1/2)^k + (1 - (-1/2)^k)/9` at every `k >= 1`, so `C_t(2^k - 1) -> 1/9`. Proved.** The word `1^k` is one run, `v(2^k - 1) = M_1^k e_0 = J M_0^k J e_0`, and the first coordinate of `M_1^k e_0` is `t + (1 - t)/9` with `t = (-1/2)^k`. `decay` checks the closed form at `k = 1, 2, 3, 4, 6, 8, 12, 16, 24, 40, 60` against the exact recursion and reads `0.111111111` from `k = 40` on.

This is the obstruction that fixes the shape of the argument. The word `1^k` has `k` ones and one run, and `C_t` stays near `1/9` along it, so no bound that decays in the binary digit sum or in the number of nonzero binary digits can hold. So the set where `abs C_t >= 1/9 - epsilon` contains infinitely many integers whose binary digit sums are unbounded. The finiteness that closes the argument must therefore be stated in signed lengths, where `2^k - 1 = 2^k - 2^0` has length `2`, and not in digit sums. Digit sums are the corner where the two-base results already in the tree live: Senge and Straus 1973 and Stewart 1980, both quoted in [two bases](../../cobham/notes/the-law-and-the-budgets.md), give finiteness when the digit sums in both bases are bounded, and they do not reach `2^k - 1`, whose binary digit sum grows without bound. Counting blocks of equal digits does reach it, since `2^k - 1` is one block: Bugeaud, Cipu and Mignotte 2013, Theorem 3.1, give an effective bound on the blocks in two independent bases for every integer, after Blecksmith, Filaseta and Nicol 1993, Theorem 1, for the powers `a^n` (the number of blocks of `a^n` in base `b` tends to infinity when `log a/log b` is irrational); both write the integer as a signed sum over its blocks ([refs](../refs.md)). The S-unit theorem below reaches the same integers for signed lengths.

The ternary layer has no such obstruction: each step of its recursion through a nonzero residue costs a factor `1/3` outright, because `al_1 = 0` for the tile `(1, 1, -1)`. The binary tile has `al_1 = -1`, and that keeps the average `-(C_t(h) + C_t(h + 1))/2` from contracting inside a run.

## Cesaro decay

A correlation is positive definite. For lags `a_1, ..., a_m` and complex `c_i`, `sum_(i,j) c_i conj(c_j) C_x(a_j - a_i) = lim (1/N) sum_(n < N) abs(sum_i c_i x(n + a_i))^2 >= 0`, since shifting the window moves each average by `O(max a/N)`. So `C_x` is positive definite with `C_x(0) = 1`, and Herglotz's theorem ([refs](../refs.md)) gives a probability measure `sigma_x` on the circle `T = R/Z` with `hat sigma_x(h) = int e(-hx) d sigma_x(x) = C_x(h)`, the spectral measure of `x`. A measure is fixed by its Fourier coefficients, so the product rule reads `sigma_z = sigma_t * sigma_3`, the convolution of the two layer measures.

**Wiener's identity: for a finite positive measure `sigma` on `T`, `(1/(2H + 1)) sum_(abs h <= H) abs(hat sigma(h))^2 -> sum_x sigma({x})^2`. Proved.** Here `abs(hat sigma(h))^2 = int int e(-h(x - y)) d sigma(x) d sigma(y)`, and the average of `e(-hu)` over `abs h <= H` is a function `D_H(u)` with `abs D_H <= 1`, `D_H(0) = 1` and `D_H(u) -> 0` at every `u != 0` in `T`. Dominated convergence on `sigma x sigma` gives the mass of the diagonal, which is `sum_x sigma({x})^2`. So `sigma` has no atoms exactly when the Cesaro mean of `abs(hat sigma)^2` tends to `0`. The proof above is complete; no source is used.

**`(1/H) sum_(1 <= h < H) C_t(h)^2 <= 2 K^2 H^(-theta)` with `theta = log2(400/303) = 0.40068...`, and `(1/H) sum_(1 <= h <= H) C3(h)^2 <= (25/9) H^(-log3(9/5))` with `log3(9/5) = 0.53502...`, at every `H >= 2`. Proved.** The words of length `m` with `r` runs number `C(m - 1, r - 1)`, one choice of the `r - 1` change points among `m - 1` gaps, so `sum_(2^(m-1) <= h < 2^m) rho^(2 r_2(h)) = rho^2 (1 + rho^2)^(m-1)` and `sum_(1 <= h < 2^L) rho^(2 r_2(h)) = (1 + rho^2)^L - 1`. With `2^(L-1) < H <= 2^L`, the run lemma gives `(1/H) sum C_t^2 <= (K^2/H)(1 + rho^2)^L <= 2 K^2 ((1 + rho^2)/2)^L <= 2 K^2 H^(-theta)`, `(1 + rho^2)/2 = 303/400`. On the ternary side, the integers `abs h <= (3^k - 1)/2` are exactly those with `k` balanced digits, so `sum 3^(-w_3(h)) = (1 + 2/3)^k = (5/3)^k` over them; take `k` least with `H <= (3^k - 1)/2`, so `3^k < 9H`, and the ladder's bound gives `(1/H) sum C3^2 <= (5/3)^k/H < (25/9) H^(log3(5/3) - 1)`.

So `sigma_t` and `sigma_3` have no atoms, and neither has `sigma_z`, whose Cesaro mean is at most that of `sigma_t` because `abs C3 <= 1`. Both layer measures are continuous and neither is Rajchman: `C_t(2^k) = -1/3` and `C3(2 . 3^k) = -1/3` for every `k`. `decay` reads the means over `1 <= h < 2^L`:

| `L` | mean `C_t^2` | bound `K^2 (303/400)^L` | mean `C3^2` | mean `C_z^2` | `sum_(1 <= h < 2^L) C_z(h)^2` |
|---|---|---|---|---|---|
| `5` | `3.993e-2` | `2.541e-1` | `1.295e-2` | `1.086e-3` | `0.034751` |
| `9` | `8.579e-3` | `8.365e-2` | `2.181e-3` | `8.042e-5` | `0.041174` |
| `13` | `1.558e-3` | `2.754e-2` | `2.435e-4` | `5.192e-6` | `0.042536` |
| `17` | `2.692e-4` | `9.069e-3` | `2.995e-5` | `3.283e-7` | `0.043031` |
| `21` | `4.573e-5` | `2.986e-3` | `3.014e-6` | `2.059e-8` | `0.043187` |
| `25` | `7.718e-6` | `9.831e-4` | `3.485e-7` | `1.289e-9` | `0.043240` |

The measured slopes of `log2` of the means from `L = 9` to `25` are `-0.6324` per bit for `C_t^2` and `-0.7882` for `C3^2`, against the proved `-0.4006` and `-0.5350`, rounded safe from `theta = 0.400682` and `log3(9/5) = 0.535026`, which `decay` prints. The proofs lose because they bound every word by its worst case.

The last column is the surprise. If the two layers were independent, the mean of `C_z^2` at height `2^L` would be the product of the two means, about `2^(-1.42 L)`, and the sum of `C_z^2` over a dyadic block would fall as `2^(-0.42 L)`. The partial sums do settle: `0.043031` at `2^17`, `0.043187` at `2^21` and `0.043240` at `2^25`, with block increments shrinking by a factor near `3` every four bits, as `2^(-0.42 L)` predicts. **That `sum_h C_z(h)^2 < infinity`, so that `sigma_z` is absolutely continuous with a square-integrable density, is a Conjecture.** It would make the convolution absolutely continuous although neither factor is Rajchman, a far stronger conclusion than the Rajchman property this chapter proves. The run lemma and the ternary bound are not enough even heuristically: they bound `C_z^2` by `K^2 rho^(2 r_2) 3^(-w_3)`, and if those two bounds decorrelated as the measured means seem to, they would bound the block sums only by a multiple of `2^(0.065 L)`, which grows; the exponent `1 - theta - log3(9/5)` is `0.064291` to six places, rounded up to `0.065`. A proof needs the two digit systems to decorrelate on average, which is a statement about joint binary and ternary digit statistics.

## Sparse in two bases

**Theorem (van der Poorten and Schlickewei 1982, Evertse 1984). Read at source.** Let `K` be a number field, `Gamma` a finitely generated subgroup of `K^*` and `a_1, ..., a_n` nonzero elements of `K`. Then `a_1 u_1 + ... + a_n u_n = 1` has only finitely many solutions `u_1, ..., u_n` in `Gamma` with no vanishing subsum on the left. The statement is read in the survey of Evertse, Gyory and Stewart ([refs](../refs.md)), which credits it to van der Poorten and Schlickewei and, independently, to Evertse. It is checked against Evertse's own paper in its Numdam scan, whose Theorem 1 gives finiteness for projective points with `S`-unit coordinates on `x_0 + ... + x_n = 0` with no vanishing sum over a proper nonempty subset of the indices; the scan's text layer drops the displayed formulas, so that equation is read from the surrounding text. The van der Poorten and Schlickewei report is unread. Evertse's introduction names the method, Schlickewei's `p`-adic form of the Thue-Siegel-Roth-Schmidt method, and the finiteness it gives is not effective.

**For every `n` only finitely many `h >= 1` have `l_2(h) + l_3(h) <= n`. Proved, from the theorem above.** Fix `h` and write it twice, `h = sum_(i in I) x_i = sum_(j in J) y_j` with each `x_i = +-2^(a_i)` and each `y_j = +-3^(b_j)`, choosing the pair with `abs I + abs J` least; then `abs I + abs J <= l_2(h) + l_3(h) <= n`, and both sides are nonempty because `h != 0`.

- Split the zero sum `sum x_i - sum y_j = 0` into minimal vanishing blocks. A block inside `I` alone (or `J` alone) could be deleted and would shorten a representation with `abs I = l_2(h)` (or `abs J = l_3(h)`), against the choice; if it were all of `I`, then `h = 0`. So every block meets both bases.
- A block `x = y` forces `+-2^a = +-3^b`, so `a = b = 0`.
- Otherwise divide a block by one of its `y`-terms `y_(j_0)` and apply the theorem to that block, which has no vanishing subsum by minimality. This gives elements of `Gamma = <-1, 2, 3>` summing to `1`; for each length there are finitely many such tuples, by the theorem with `K = Q` and every `a_i = 1`.
- In a block's tuple, `x_i/y_(j_0) = +-2^(a_i) 3^(-b_(j_0))`, and by unique factorisation the value fixes `a_i` and `b_(j_0)`. So every exponent ranges over a finite set as `h` does, and `h <= n 2^(max a_i)` is bounded.

The same proof runs for any two multiplicatively independent bases `a, b >= 2` with `Gamma = <-1, a, b>`: `a^i b^(-j) = a^(i') b^(-j')` forces `a^(i - i') = b^(j - j')`, and independence then forces equal exponents. So every set `{h : l_a(h) <= P, l_b(h) <= Q}` is finite. With the counts of the objects, every set `{h : r_2(h) <= R, w_3(h) <= W}` is finite.

This contains the digit-sum corner. An integer whose digit sums in bases `a` and `b` are at most `c` has `l_a, l_b <= c`, so the finiteness half of Senge and Straus 1973 follows for every independent pair. The converse fails: `2^k - 1 = 2^k - 2^0` has `l_2 = 2` and binary digit sum `k`, so bounded signed lengths are a strictly wider condition than bounded digit sums, and the words where `C_t` stays large, such as `2^k - 1`, lie only in the wider set. So the sentence in [two bases](../../cobham/notes/the-law-and-the-budgets.md) that finiteness arrives only at the corner where every digit sum is bounded is too narrow, as the digit-change theorem below already shows: a bounded signed length bounds the number of blocks by the lemma below, so finiteness arrives wherever both signed lengths are bounded, a set that contains the digit-sum corner. The S-unit proof reaches these sets directly.

The S-unit proof gives no bound on the largest element, but the finiteness is effective, through digit changes. For `m >= 1` with base-`b` digits `d_(L-1) ... d_0`, `d_(L-1) != 0`, let `BC(m, b)` be the number of `j < L - 1` with `d_j != d_(j+1)`, the number of times a digit different from the previous one is read; so `BC(h, 2) = r_2(h) - 1`.

**Theorem (Bugeaud, Cipu and Mignotte, Theorem 3.1). Read at source.** Let `a` and `b` be multiplicatively independent integers. Then `BC(m, a) + BC(m, b) >= log log m/(log log log m + C) - 1` for `m > 25`, where `C` is a positive number effectively computable in terms of `a` and `b` only. The survey ([refs](../refs.md)) writes `m` in each base as a signed combination of the powers at its digit changes and follows Stewart's proof for bounded digit sums, which runs on Baker's linear forms in logarithms; it gives no value of `C`. It also notes that `k` nonzero digits allow at most `2k` changes. The lemma below extends that to signed sums.

**`BC(h, b) <= 2 l_b(h) - 1` at every `h >= 1` and every base `b >= 2`. Proved.** Take a representation of `h` by `l = l_b(h)` terms `+-b^a`. It never uses one power with both signs, and never `b` copies of one signed power, since `b . b^a = b^(a+1)` would shorten it; so it merges into `h = sum_i c_i b^(a_i)` with `a_1 < ... < a_k`, `1 <= abs c_i <= b - 1` and `k <= l`. Put `e_j = c_i` at `j = a_i` and `0` elsewhere, and read digits from the bottom with a borrow: `g_0 = 0`, `d_j = (e_j + g_j) mod b`, `g_(j+1) = (e_j + g_j - d_j)/b`.

- If `g_j` is `0` or `-1`, then `e_j + g_j` lies in `[-b, b - 1]` and `g_(j+1)` is `0` or `-1`; and `sum_(j < J) d_j b^j + g_J b^J = sum_(j < J) e_j b^j` at every `J`.
- Past `a_k` the borrow is constant, and it is `0`: `g_J = -1` with `J > a_k` would give `h = sum_(j < J) d_j b^j - b^J < 0`. So the `d_j` are the standard digits of `h`, and its top digit sits at `L - 1 <= a_k`.
- Where `e_j = 0` the borrow passes unchanged and `d_j` is `0` or `b - 1` as the borrow is `0` or `-1`; so `d_j = d_(j+1)` unless `j` or `j + 1` is some `a_i`.
- Each `a_i` lies in at most two of the pairs `(j, j + 1)` with `j + 1 <= L - 1`, and `a_k` in at most one, since `a_k + 1 > L - 1`. Hence `BC(h, b) <= 2k - 1 <= 2l - 1`.

`changes` checks the lemma at every `h < 3^12` for `b = 2` and `b = 3`: `0` failures, with equality at `7141` and `232` lags, the first `2` and `3`. It computes `l_b` by the digit recursion `l(bq) = l(q)`, `l(bq + r) = min(r + l(q), b - r + l(q + 1))`, and a breadth-first search over `[-3^13, 3^13]` agrees at every lag; it also reads `l_3 = w_3` at every lag below `3^12`.

So `l_a(h) + l_b(h) <= n` gives `log log h/(log log log h + C) <= 2n - 1`, an effective bound on `h`. The constant `C` is not computed here, so the census below still cannot certify that its list is complete past its cut.

## The limit

**`C_z(h) -> 0` as `abs h -> infinity`. Proved.** `C_z` is even, so take `h >= 1`. Fix `epsilon > 0` and suppose `abs C_z(h) >= epsilon`. Both factors are at most `1` in absolute value, so `abs C_t(h) >= epsilon` and `abs C3(h) >= epsilon`. The run lemma gives `r_2(h) <= R = log(K/epsilon)/log(1/rho)`, and the ternary bound gives `w_3(h) <= W = 2 log3(1/epsilon)`. The counts of the objects then give `l_2(h) + l_3(h) <= R + 1 + W`, and the S-unit lemma leaves finitely many such `h`. So for every `epsilon` the set where `abs C_z >= epsilon` is finite.

This closes the open row of the ladder. At `epsilon = 10^(-3)` the proof confines the exceptions to `r_2 <= 20` and `w_3 <= 12`, a finite set of unknown size. `windows` finds `364` lags `1 <= h < 2^25` with `abs C_z >= 10^(-3)`; the last is `h = 3188646`, with `r_2 = 12`, `w_3 = 2` and `C_z = -2933/(2^18 3^2) = -1.2432e-3`, below the maximum `1.374e-3` of its window `[2^21, 2^22)` at `3148800`, and the sup over `[2^22, 2^25)` is `7.227e-4`; that is a reading of a finite range, not a bound.

**The rate: `abs C_z(h) <= K 3^(-G(h)/4) < 1.0093 . 3^(-G(h)/4)` at every `h > 25`, with `G(h) = log log h/(log log log h + C) - 1` and `C` the constant of Theorem 3.1 at the bases `2` and `3`. Proved, from that theorem; effective, with `C` not explicit.** Two counts feed the theorem: `r_2(h) = BC(h, 2) + 1`, and the lemma at `b = 3` with `l_3 <= w_3` gives `2 w_3(h) >= BC(h, 3) + 1`. So `r_2(h) + 2 w_3(h) >= BC(h, 2) + BC(h, 3) >= G(h)`, an inequality that holds even if `BC` counts blocks, one more than changes. The two layer bounds multiply, `abs C_z(h) <= K rho^(r_2(h)) 3^(-w_3(h)/2)`, and `rho^4 = 10609/40000 < 1/3` gives `rho < 3^(-1/4)`, so `abs C_z(h) <= K 3^(-(r_2(h) + 2 w_3(h))/4) <= K 3^(-G(h)/4)`. Where `log log log h >= C`, `G(h) >= log log h/(2 log log log h) - 1`, and the bound reads `abs C_z(h) <= 1.3284 (log h)^(-0.1373/log log log h)`, from `K 3^(1/4) = 1.328305` and `log(3)/8 = 0.137327`; `changes` prints the constants. Bounding each layer alone and splitting `r_2 + 2 w_3` in halves would give only `3^(-G(h)/8)`.

The proved bound decays slower than every power of `log h`, and on any window a computer can reach it is indistinguishable from a constant; the tail is carried by the integers that are sparse in both bases, which the census lists. What is not proved is a number. `C` is effectively computable from `2` and `3`, but the survey gives no value and none is computed here, so for no `epsilon` is a lag known past which `abs C_z <= epsilon`.

## In measure language

A finite measure on `T` is Rajchman when its Fourier coefficients tend to `0`. The two layer measures are not: `hat sigma_t(2h) = hat sigma_t(h)` and `hat sigma_3(3h) = hat sigma_3(h)`, so `hat sigma_t(2^k) = -1/3` and `hat sigma_3(2 . 3^k) = -1/3` for all `k`. That is the dilation invariance the ladder reads as the place slot, one integer ratio, and the one the [dilated stacks](../../dilations/notes/dilations.md) chapter proves for digit measures, `mu^(bt) = mu^(t)`, and [weights](../../weights/notes/weights.md) for every weighted design. Both layers are continuous, by the Cesaro section. What the limit proves is:

**`sigma_z = sigma_t * sigma_3` is a continuous Rajchman measure, the convolution of two continuous measures that are not Rajchman, invariant under multiplication by `2` and by `3` respectively. Proved.**

The proof uses one property of each factor, decay in the signed length at its own base, and it gives:

**Let `a, b >= 2` be multiplicatively independent, and let `mu`, `nu` be finite measures on `T` with `abs hat mu(h) <= f(l_a(h))` and `abs hat nu(h) <= g(l_b(h))` at every `h != 0`, for functions `f`, `g` that tend to `0` at infinity. Then `mu * nu` is Rajchman. Proved.** If `abs(hat mu(h) hat nu(h)) >= epsilon`, then `f(l_a(h)) >= epsilon/norm(nu)` and `g(l_b(h)) >= epsilon/norm(mu)`, so both signed lengths are bounded, and the S-unit lemma at the pair `(a, b)` leaves finitely many `h`; `l_a(-h) = l_a(h)` covers negative lags. Neither factor needs to be Rajchman, because `l_a(a^k) = 1` lets `hat mu(a^k)` stay where it is.

So the statement that the convolution of two non-Rajchman measures at multiplicatively independent bases is Rajchman needs a hypothesis. Invariance alone is not enough: `delta_0` is invariant under multiplication by every integer, and `delta_0 * delta_0 = delta_0`. With decay in signed lengths it holds, and Thue-Morse needed the run lemma precisely because its coefficients do not decay in digit sums.

Against the tree's own theorems:

- One base never escapes. The dilated stacks and weighted designs chapters prove that a measure invariant under one integer ratio has `hat mu(b^k t) = hat mu(t)`, so it is never Rajchman. That holds for `sigma_t` and `sigma_3`, and it fails for `sigma_z`, which is invariant at neither base: `C_z(2h) = C_t(h) C3(2h)` and `C_z(3h) = C_t(3h) C3(h)`.
- Two ratios inside one self-similar measure is a different theorem. The [staircase schedules](../../dimensions/notes/staircase-schedules.md) chapter quotes, for the unequal split, Rapaport 2022, Corollary 1.6: a self-similar measure of an affinely irreducible system that is not Rajchman has every ratio a power of one algebraic integer, so an affinely irreducible system with ratios `1/2` and `1/3` gives Rajchman measures. There the two ratios act inside one measure. Here each measure keeps one ratio, and the two meet only through the convolution, which is not self-similar. Neither statement implies the other.
- Rudin-Shapiro is the trivial case: `C_y(h) = C_r(h) C3(h) = 0` at every `h != 0`, so `sigma_y` is Lebesgue measure, and the ladder already printed that `y` is white.
- The place slot is one integer ratio. The moire XOR breaks it in the spectrum by pairing two ratios that are multiplicatively independent, and the price is the proof's rate: effective, slower than every power of `log h`, and with a constant not yet computed.

In print, from the sources read ([refs](../refs.md)):

- Salem asked about the Fourier transform of a measure carried by the sum of two Cantor sets, and Straus conjectured the digit-sum finiteness as an answer in part (Besineau 1972, page 401).
- Senge and Straus characterised the ratios `a, b` for which the Fourier transform of `mu_a * mu_b`, the convolution of the two central Cantor measures, does not tend to `0`: only when `1/a` and `1/b` are Pisot and `log b/log a` is rational (Nazarov, Peres and Shmerkin 2012, introduction and Section 4; the original is unread). At `a = 1/3` and `b = 1/4` neither factor has a Fourier transform tending to `0`, and the convolution does. So the phenomenon and its arithmetic mechanism are old; what this chapter adds is its form for the spectral measures of two digital sequences, whose coefficients decay in signed length and not in digit sums.
- For digit sums in coprime bases the product rule `C_(f_1 f_2) = C_(f_1) C_(f_2)` is Besineau's (seminar Lemma 4) and the mutual singularity of the two spectral measures is Kamae's (Theorem B, irrational weights). For one `q`-multiplicative sequence the spectral measure is atomic or singular (Coquet, Kamae and Mendes France 1977, Theorems 5 to 7; Queffelec 1979, Theorem 2), and its correlation never tends to `0` (Coquet 1979, Section II.5), the one-base statement above.
- At `h = 2^k` the limit is a corollary: `C_z(2^k) = -C3(2^k)/3`, `w_3(2^k) >= (BC(2^k, 3) + 1)/2` by the lemma above, and Blecksmith, Filaseta and Nicol's Theorem 1, which Bugeaud, Cipu and Mignotte 2013 also state in their Section 3, sends `BC(2^k, 3)` to infinity.
- The Cesaro decay of `C_t^2` is known with its sharp exponent: the correlation dimension of the Thue-Morse spectral measure is `3 - log2(1 + sqrt 17) = 0.64298` (Zaks, Pikovsky and Kurths 1997; Zaks 2001, Section 2.2), against the proved `0.4006` and the measured `0.6324`; the same section states `3 - log3 11 = 0.8173` for the 3-tupling sequence, whose correlation obeys the recursion of `C3`.
- Not found in the sources read: a Rajchman statement for the product of two `q`-multiplicative sequences in independent bases, a bound on the Thue-Morse correlation in the number of runs, the signed-length form of the convolution theorem, and `BC <= 2 l - 1` for signed lengths.

## The census

The test has two halves: every lag in a window, and every lag that is sparse in both bases up to a far cut. A counterexample to the limit would have to live among the sparse lags, since everywhere else one of the two layers is already small.

**On every window `[2^k, 2^(k+1))`, `k <= 24`, the maximum of `abs C_z` is the value below. Verified.** `windows` computes `C_z` at every lag below `2^25` with the exact recursions and prints each maximum as an exact fraction:

| `k` | max `abs C_z` | at `h` | `r_2` | `w_3` | exact |
|---|---|---|---|---|---|
| `12` | `4.123e-3` | `4374` | `8` | `2` | `19/(2^9 3^2)` |
| `13` | `8.138e-3` | `13122` | `8` | `2` | `25/(2^10 3)` |
| `14` | `2.713e-3` | `26244` | `8` | `2` | `25/(2^10 3^2)` |
| `15` | `3.328e-3` | `38880` | `6` | `4` | `23/(2^8 3^3)` |
| `16` | `2.713e-3` | `104976` | `8` | `4` | `25/(2^10 3^2)` |
| `17` | `1.551e-3` | `131220` | `8` | `4` | `343/(2^13 3^3)` |
| `18` | `1.949e-3` | `354312` | `10` | `4` | `431/(2^13 3^3)` |
| `19` | `1.374e-3` | `1048320` | `2` | `8` | `19/(2^9 3^3)` |
| `20` | `1.551e-3` | `1049760` | `8` | `4` | `343/(2^13 3^3)` |
| `21` | `1.374e-3` | `3148800` | `4` | `8` | `19/(2^9 3^3)` |
| `22` | `5.430e-4` | `6259194` | `8` | `4` | `23059/(2^19 3^4)` |
| `23` | `7.227e-4` | `9570312` | `10` | `4` | `1705/(2^18 3^2)` |
| `24` | `4.144e-4` | `25509168` | `12` | `4` | `2933/(2^18 3^3)` |

The window `k = 20` reproduces the ladder's `1.55e-3` at `1049760`. Every maximum sits at a lag with few runs and few balanced digits, such as `13122 = 2 . 3^8` and `1048320 = 2^20 - 2^8`, and the decay across windows is slow and uneven, as a limit carried by sparse integers should be. `13122` shows the mechanism at its simplest: its ternary layer is `C3(2 . 3^8) = -1/3` exactly, and its binary word has `8` runs.

**Below `10^40` exactly `15196` integers have `r_2(h) <= 6` and `w_3(h) <= 6`; the largest is `283462598656`, and none lies in `[2^39, 10^40)`. Verified.** `sparse 6 6 40` walks all `14536738777` balanced-ternary words with at most `6` nonzero digits and top position at most `84`, which covers every positive integer below `10^40` of weight at most `6`, each once since balanced ternary is unique. It keeps those whose low `64` bits have at most `6` digit changes, a necessary condition, and checks the run count, the weight and the cut exactly in `192`-bit arithmetic. All `15196` filter survivors are genuine. The generator prints the whole list. By blocks of eight bits from `2^0` they number `246, 6189, 7700, 1025, 36`, then `0` in each of the twelve blocks up to the cut, and nothing lies within `0.02` of the cut. The largest, `283462598656`, has `C_z = -7/(2^16 3^18) = -2.757e-13`. The largest `abs C_z` over the members at least `2^20` is `1.345e-3`, at `1049598`; at least `2^28` it is `1.606e-4`, at `545258464`; at least `2^32` it is `1.584e-5`, at `34896543768`. `2205` members have `C_z = 0`, through a zero of either layer, as at `C_t(5) = 0` or `C3(1) = 0`.

**Below `10^20` exactly `528556` integers have `r_2(h) <= 8` and `w_3(h) <= 8`; the largest is `4035260175885545472`. Verified.** `sparse 8 8 20` walks `20834959549` words, and two of its `528558` filter survivors fail the exact check. Its eight-bit blocks hold `255, 32706, 257842, 199689, 35657, 2295, 103, 9`, then `0` in the last block below the cut, and nothing lies within `0.02` of it. The largest member has `C_z = 6.99e-16`, and the largest `abs C_z` over the members at least `2^32` is `1.199e-4`, at `6972506136`.

Both lists thin out long before their cuts. The S-unit theorem says each list is finite even without a cut, and Theorem 3.1 bounds where it ends only through the uncomputed `C`, so the census is evidence about where the list ends and not a proof of it. What it does falsify is a cheap failure of the limit. No integer below `10^40` with `r_2 <= 6` and `w_3 <= 6` carries a correlation above `1.6e-5` once it passes `2^32`.

## Open

- **The rate.** Proved and effective, not explicit: `abs C_z(h) <= 1.0093 . 3^(-G(h)/4)`, `G(h) = log log h/(log log log h + C) - 1`, through Bugeaud, Cipu and Mignotte Theorem 3.1 and `BC(h, b) <= 2 l_b(h) - 1`. Not proved: a value of `C`, which needs Stewart's Baker argument run with explicit linear-form bounds and would also certify the census lists complete; and the true order of decay, the window maxima falling about tenfold from `2^12` to `2^24`.
- **The square-summable spectrum.** The partial sums of `C_z(h)^2` read `0.043240` at `2^25` and settle at a geometric rate. A proof that they converge would make `sigma_z` absolutely continuous with a square-integrable density. It needs the binary and ternary digit statistics of `h` to decorrelate on average, which neither layer bound sees.
- **The constant per run.** It lies in `[0.711872, 0.717636]`, the lower end coming from the pattern `(1, 1, 2)` and the upper from one quadratic form. The exact value is the joint spectral radius of `{A(t) : t = (-1/2)^b, b >= 1}`, and the guess is that the pattern `(1, 1, 2)` attains it.
- **Other layers.** The proof needs only decay in signed length at two coprime bases, which are then multiplicatively independent. Every base-`b` signed design has a block recursion of the same kind, and a run lemma for each pair of layers would make their moire XOR Rajchman in the spectrum.
