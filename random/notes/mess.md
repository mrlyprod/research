---
title: The mess
lead: The machine's nonlinear step is a symmetric Boolean function of a neighbour count; Lucas reads its degree, Krawtchouk its spectrum, and past five neighbours every Source of the tree sits below the nonlinearity ceiling, held there by its own thinness.
figure: research-random-mess
---

A round of the machine writes one side of the grid with a nonlinear function of the other. The step this page studies is the totalistic one: each written cell XORs `g(count)`, where `count` is the number of live cells in its design neighbourhood, `k` cells wide, and `g` is a Source of the tree, the indicator of a set of integers. As a function of the `k` neighbours that is a symmetric Boolean function, one that depends only on the Hamming weight of its input, so the questions a cipher designer asks of a round function have exact answers here. The page answers them for the 22 fixed Sources at every `k` from 1 to 64, which covers the design fills: `k = 8` for the carpet, `2..9` for the base-2 and base-3 designs at level 1, `16, 25, 36, 49, 64` for their level-2 squares.

Two generators print every number below and share no code: `lab/rs/mess` builds the Sources from the library and runs the census in about eight seconds, and `lab/py/mess` rebuilds them from definitions and rechecks the census three ways in about seven. The headline is a no. No Source comes near the best symmetric function, the reason is one line, and the obstruction the question assumed, that a symmetric function cannot be bent, is false.

## The object

A Source is a set `S` of integers, and `g(j) = 1` exactly when `j` is in `S`. The 22 fixed Sources are the evens, odds, primes, powers of two (binary), Fibonacci numbers, the odd squares, and the fill and void counts of the eight plane designs carpet, net, tree, void, point, dust, line and star at the odd sides `1, 3, 5, ...` (the library's `mrly_life::source::sequence`). They fall into 14 distinct sets. The sets keep 0 and 1 as `sequence()` returns them: 10 of the 22 Sources contain 0 (the evens, Fibonacci, carpet voids, net fills, tree voids, void voids and the twins of the last four, point fills, dust voids, line fills and star fills) and 14 contain 1 (the odds, binary, Fibonacci, grid squares, carpet fills, carpet voids, net voids, tree fills, void fills and the twins of the last five). The machine's `counts()` drops both by default, and the counts below then change. Dropping 0 turns the carpet voids into the net voids, so 14 sets become 13, and gives 576 of 684 tight cells (tight is defined below), 362 of 744 full-degree cells and a best of 0.7862, with the odds still affine; dropping 0 and 1 gives 13 sets, 578 of 741 tight cells, 368 of 806 full-degree cells and a best of 0.7803, with 13 non-affine classes since the odds stops being affine. The `DROPPED LOW TERMS` block of `lab/rs/mess` prints these numbers and `lab/py/mess` recomputes them from definitions. A design and its complement have `fill(F, n) + fill(complement of F, n) = n^2`, so the fills of a classic are the voids of its anti and the 16 design Sources pair into 8 sets; the pairs are carpet and point, net and dust, tree and line, void and star. The design sets are the value sets of quadratics in `m` at side `2m + 1`: the carpet fills `(m+1)(3m+1)` are `1, 8, 21, 40`, the carpet voids are all squares, the tree fills `n(m+1)` are `1, 6, 15, 28`.

The function is `f(x) = g(|x|)` for `x` in `F_2^k`, and `g(0..k)` is its simplified value vector. Five quantities are read off it: the ANF degree, the number of ones, the nonlinearity `NL = 2^(k-1) - max|W|/2` with `W` the Walsh transform, the differential uniformity `max over a != 0 and b of #{x : f(x+a) + f(x) = b}`, and the order of correlation immunity. The ceiling `C(k)` is the most nonlinearity any symmetric function can have, `2^(k-1) - 2^(k/2-1)` for even `k` and `2^(k-1) - 2^((k-1)/2)` for odd `k`.

## The degree by Lucas

**The ANF of a totalistic mess is the binary-domination transform of its Source. Proved.** The coefficient of the monomial on a set `T` of `w` cells is the XOR of `f` over the subsets of `T`, and `j` of those subsets have weight `j` in `C(w, j)` ways, so it is the XOR over `j` of `C(w, j) g(j)`. By Lucas's theorem `C(w, j)` is odd exactly when every binary digit of `j` is at most the same digit of `w`, written `j <= w` bitwise. Hence the coefficient depends only on `w`:

```
c(w) = XOR of g(j) over j with j <= w bitwise
degree = max { w <= k : c(w) = 1 }
```

The matrix of the transform has a 1 at `(j, w)` when `j <= w`; it is the level-`L` Kronecker power of the base-2 design of code 11, a [gasket](https://en.wikipedia.org/wiki/Sierpi%C5%84ski_triangle), and it is its own inverse modulo 2, since for `i <= w` the number of `j` with `i <= j <= w` is `2^(popcount(w) - popcount(i))`, which is odd only at `i = w`. So `g` and `c` determine each other, `c` does not depend on `k` except by truncation, and a Source's whole family of messes is one infinite ANF vector. This is Proposition 2 of [Canteaut and Videau 2005](../refs.md), and their Theorem 1 adds that the degree is below `2^t` exactly when `g` has period `2^t`; nothing here is new before the Sources. The figure draws the transform at `L = 6`, the 64 by 64 gasket read through the library as code 11 and asserted equal to `j <= w` entry by entry, with the prime rows in orange; the strip beneath is the column parity, which is `c` for the primes.

**The odds are the parity and the binary Source is Thue-Morse. Proved.** For `g(j) = j mod 2`, a nonempty set of digits of `w` contains digit 0 in half its subsets, so `c(w) = [w = 1]` and the mess is the XOR of the neighbours: degree 1, nonlinearity 0, resilient to order `k - 1`. The evens are its complement. For the powers of two the dominated `j` are the single digits of `w`, so `c(w) = popcount(w) mod 2`, and the first 105 terms of `c` are the Thue-Morse sequence A010060. The degree is then the largest odious number at most `k`, and it is `k`, `k - 1` or `k - 2`: among three consecutive integers two form a pair `2n, 2n + 1`, which differ in parity of popcount. The degree is `k` at 33 of the 64 values of `k`.

**The primes have full degree at every `k = 2 mod 4`. Proved.** An even `k` dominates only the prime 2, and only when its digit 1 is set, so `c(k) = [k = 2 mod 4]`: degree `k` there and below `k` at every multiple of 4. For odd `k` no closed form is claimed; the census prints the degree sequence `-1, 2, 2, 2, 5, 6, 6, 6, 6, 10, 11, 11, 11, 14, ...` to `k = 64`, which ends `58, 58, 58, 58, 62, 62, 62` and is not in the dump of the OEIS at windows of 30, 20, 14 and 10 terms from any of its first six starting terms. Across the 12 non-affine classes, 334 of 744 cells `(class, k)` with `k = 3..64` have full degree.

## The spectrum by Krawtchouk

The Walsh coefficient of a symmetric function at an input of weight `v` is `W(v) = sum over w of (-1)^g(w) K_w(v)`, where the Krawtchouk polynomial `K_w(v)` is the coefficient of `z^w` in `(1 - z)^v (1 + z)^(k-v)`: [Savicky](../refs.md) writes the generating function and Proposition 3 of Canteaut and Videau gives the formula. The autocorrelation at a shift of weight `u` has two exact routes, the direct count `sum over i, j of C(u, i) C(k-u, j) (-1)^(g(i+j) + g(u-i+j))` and the Wiener-Khinchin sum `2^(-k) sum over v of W(v)^2 K_v(u)`. Everything is integer arithmetic, 128 bits in Rust and unbounded in Python.

**The formulas match a brute-force table cell by cell. Verified.** The falsification plan was to print every class and `k` where the closed form and a truth table disagree. For `k <= 20` the Python study builds the `2^k` table, takes the fast Walsh-Hadamard transform, the Moebius transform for the degree and the transform of `W^2` for the autocorrelation, and compares seven quantities of 14 classes at 20 values of `k`: 1960 cells, 0 disagreements. The Rust census carries its own brute-force check at `k <= 14` on 980 cells and the Python exact route covers all `k <= 64` on 6272 cells with a binomial-parity ANF and Parseval's identity, 0 disagreements.

## The ceiling and the frontier

The question assumed that symmetric functions are capped below bent nonlinearity. They are not.

**Symmetric bent functions exist at every even `k`, and there are four. Refuted; Verified.** The function `sigma_2 = XOR of x_i x_j` over pairs is the symmetric function with `g(w) = C(w, 2) mod 2`, the pattern `w mod 4 in {2, 3}`; its bilinear form `J - I` is invertible modulo 2 for even `k`, so it is bent, and so are `sigma_2` plus a linear `sigma_1` and plus 1. [Savicky 1994](../refs.md) proves in Theorem 3.3 that these four are all the symmetric bent functions: a symmetric function of an even number of variables is bent exactly when `g(w + 2) = g(w) + 1` throughout. For odd `k` the same quadratics reach `2^(k-1) - 2^((k-1)/2)`, the maximum for symmetric functions by Maitra and Sarkar 2002, which this page reads only as Canteaut and Videau state it in their Section VII. The real cap is on degree: the nonlinearity ceiling is reached only by degree 2.

The lab enumerates all `2^(k+1)` symmetric functions for `k = 1..26` with a Gray code, keeping the spectrum and the ANF mask incrementally, and the Python study rebuilds the same table as a matrix product for `k <= 20`: 1380 cells agree.

**The ceiling is met by exactly four functions at every `k` from 2 to 26, and full degree costs exactly one. Verified.** The best nonlinearity over all symmetric functions is `C(k)` and exactly 4 functions attain it, all of degree 2. The best at full degree `k` is `C(k) - 1`, attained by exactly 8 functions for every `k` from 3 to 26, which has the same count as the eight functions of Proposition 20 of Canteaut and Videau. The construction is one line and is **Proved**: flip the value of a bent or maximal quadratic at the single input of weight `k`; every Walsh coefficient moves by 2, so `NL >= C(k) - 1`, and the top coefficient `c(k)`, which is 0 for a quadratic, becomes 1. Degree is the cheap half of the pair; nonlinearity is the expensive half.

| k | ceiling | functions at ceiling | best at full degree | functions there | balanced functions |
|---|---|---|---|---|---|
| 4 | 6 | 4 | 5 | 8 | 2 |
| 8 | 120 | 4 | 119 | 8 | 6 |
| 12 | 2016 | 4 | 2015 | 8 | 2 |
| 16 | 32640 | 4 | 32639 | 8 | 2 |
| 20 | 523776 | 4 | 523775 | 8 | 6 |
| 24 | 8386560 | 4 | 8386559 | 8 | 50 |
| 26 | 33550336 | 4 | 33550335 | 8 | 6 |

Between those two ends the frontier is not monotone: the best cubic symmetric function has nonlinearity `2^(k-2)` for even `k` and `C(k)/2` for odd `k`, about half the ceiling, for every `k` from 4 to 26, and the best at degree `k - 1` is `C(k) - 2` with 4 functions.

**The balanced symmetric functions are counted by A200147. Proved; Verified.** A symmetric function is balanced when `sum of C(k, w) s_w = 0` with `s_w = (-1)^g(w)`. Put `x_w = (1 + (-1)^w s_w)/2`: since `sum of (-1)^w C(k, w) = 0`, this is a bijection onto the 0-1 arrays with zero `k`-th difference. The exhaustive counts for `k = 1..26` are `2, 2, 4, 2, 8, 2, 16, 6, 32, 2, 64, 2, 144, 14, 256, 2, 512, 2, 1024, 6, 2048, 2, 4096, 50, 8192, 6`, and they are the terms `A200147(2..27)`. The trivial ones are the 2 affine functions at even `k` and `2^((k+1)/2)` at odd `k`; the others occur at `k = 8, 13, 14, 20, 24, 26` and at no other `k` up to 26. By the same bijection the record gives the counts to `k = 51`.

## The Sources against the frontier

**No Source gets near the ceiling past `k = 5`. Verified.** The primes mess is exactly the quadratic pattern `w mod 4 in {2, 3}` for `k < 5`: it is `x_1 x_2` at `k = 2`, the majority at `k = 3` and the bent `sigma_2` at `k = 4`, all at the ceiling. At `k = 5` the prime 5 breaks it, and the nonlinearity is 11 of 12. From `k = 6` to 64 no non-affine class exceeds `0.7867` of the ceiling, reached by Fibonacci at `k = 12`; the leaders are carpet voids at 6, net voids at 7 and 8, carpet voids at 9, net fills at 10, Fibonacci from 11 to 27 and the primes from 28 to 64.

| class | twins | k = 8 degree | k = 8 NL over ceiling | k = 64 degree | k = 64 NL over ceiling | k = 64 share of ones |
|---|---|---|---|---|---|---|
| primes | | 6 | 0.3000 | 62 | 0.4726 | 0.2363 |
| fibonacci | | 7 | 0.3167 | 64 | 0.1801 | 0.0901 |
| binary | | 8 | 0.3083 | 64 | 0.1987 | 0.0994 |
| grid squares | | 7 | 0.0667 | 39 | 0.0435 | 0.0217 |
| carpet fills | point voids | 8 | 0.0750 | 59 | 0.0316 | 0.0158 |
| carpet voids | point fills | 8 | 0.5417 | 47 | 0.1648 | 0.0824 |
| net fills | dust voids | 8 | 0.4750 | 64 | 0.1927 | 0.0964 |
| net voids | dust fills | 6 | 0.5500 | 64 | 0.1648 | 0.0824 |
| tree fills | line voids | 6 | 0.3000 | 63 | 0.1223 | 0.0611 |
| tree voids | line fills | 8 | 0.4750 | 64 | 0.1258 | 0.0629 |
| void fills | star voids | 3 | 0.5333 | 59 | 0.0594 | 0.0297 |
| void voids | star fills | 8 | 0.4750 | 64 | 0.0544 | 0.0272 |
| odds | | 1 | 0.0000 | 1 | 0.0000 | 0.5000 |
| evens | | 1 | 0.0000 | 1 | 0.0000 | 0.5000 |

**The bound is thinness. Proved; it binds on 632 of 684 cells, on the primes only from `k = 45`. Verified.** A Boolean function is at distance `ones` from the constant 0 and `2^k - ones` from the constant 1, so `NL <= min(ones, 2^k - ones)` with `ones = sum over j in S of C(k, j)`. For the primes `ones` is the binomial transform of the prime indicator, A121497, matched at `k = 1..33`. Every Source other than the odds and evens is thin, a set with density tending to zero, and the binomial weight concentrates within a few `sqrt(k)` of `k/2`, so `ones / 2^k` is about the local density of `S` there: 0.2363 for the primes at `k = 64`, near 0.09 for Fibonacci and binary, near 0.02 for the odd squares. The bound is nearly an equality: **`NL = min(ones, 2^k - ones)` on 632 of the 684 cells `(class, k)` of the 12 non-affine classes with `k >= 8`. Verified.** The eight design classes and the odd squares meet it at 55 to 57 of the 57 values of `k` from 8 to 64, Fibonacci at 53, binary at 51, and the primes at the 20 values from `k = 45`. For those cells the best affine approximation of the mess is a constant: the mess is zero with probability `1 - ones / 2^k`, the constant 0 agrees with the mess on 76 percent of the inputs at `k = 64` for the primes and on more than 90 percent for every other non-affine class. The lowest differential uniformity over `k >= 8` is `0.7206` of `2^k`, by Fibonacci at `k = 14`, against `0.5` for a bent function; no balanced non-affine Source cell exists past `k = 3`.

**A coin-flip `g` beats every Source at `k = 16, 24, 32, 48, 64`; the tree's own random Source is level with the best of them. Verified.** Over 2000 draws a symmetric function with each `g(w)` a fair coin has median `NL / C` of `0.5417, 0.5901, 0.6216, 0.6582, 0.6844, 0.7031` at `k = 8, 16, 24, 32, 48, 64`, with full degree in about half the draws. The Source `Random(seed)`, a subset of random size, has medians `0.3167, 0.3819, 0.4050, 0.4014, 0.4339, 0.4392`, because its density is uniform on the unit interval and so is often thin. That median is below the best Source at `k = 8, 16, 48, 64` and above it at 24 and 32; the primes at `k = 64`, `0.4726`, sit just above it and far below the coin's `0.7031`. At those five `k` the median coin-flip beats every Source: a dense `g` with no 2-adic structure is better than any Source, and the frontier is better again, full degree at `C(k) - 1`.

**Rule 30's local map is at the ceiling for three variables and has a linear structure.** The map `l XOR (c OR r)` has degree 2, four ones, nonlinearity 2, which is `C(3)`, and differential uniformity 8, so its derivative along `l` is constant. The primes mess at `k = 3` is the majority and has the same four numbers: degree 2, balanced, nonlinearity 2, differential uniformity 8. Both are optimal at three variables and both leak a linear structure; they differ in kind, the mess being symmetric in all three inputs and the rule in two. **Verified.**

## What died and what is left

- Dead: a Source that is both of full degree and of high nonlinearity at every `k`. The degree is full at 33 of 64 values of `k` for the binary Source and at the residues 2 mod 4 for the primes, and the nonlinearity is thin-limited for all.
- Dead: the premise that symmetric functions are capped below bent. The cap is on degree.
- Alive: the thin-set law as a theorem. If the density of `S` near `k/2` tends to zero, `NL / 2^(k-1)` tends to zero at the rate `2 ones / 2^k`; for the primes that is a slow `1 / ln k` by the prime number theorem, for squares `k^(-1/2)`. The census shows it at 632 cells and does not prove it.
- Alive: a mess that fills the gap. The frontier gives the target, `C(k) - 1` at full degree, and the Sources give the obstacle, density. A Source built from a dense 2-adic set, such as the residues `{2, 3} mod 4` with a defect at the all-ones weight, is the first candidate.
- Alive: the ANF vectors of the squares and of Fibonacci, which are stripes of period `2^m` with defects; the dominated count of the primes is A102210; the Fibonacci one is absent from the OEIS dump at windows of 16.
- Alive: the round, not the cell. A written cell reads one neighbourhood of one function, but a round stamps and spreads; the nonlinearity of the composed round is the next question, and this page prices only the cell.
