---
title: Necessary Exits from the Automatic Class: A Second Base or an Exercised Nonlinearity
lead: A generator built from designs at one base leaks its base to a single kernel test; multiplying in a design at a multiplicatively independent base sends every shift correlation to `0`, so that Thue-Morse times the base-`3` parity has the spectral measure `sigma_t * sigma_3`, Rajchman by the S-unit theorem though each factor is invariant under its own base, while a one-base kernel test still detects it with significance growing as `N^0.194` and its entropy stays `0`; and an elementary automaton fed by a design leaves its seed's automatic class only through a base that is not a power of its characteristic or a nonlinearity exercised at unbounded depth, a necessary condition that rules `22` and `182` show is not sufficient.
date: 2026-10-26
figure: paper-design-generators
---

## Introduction

![The certificate behind the run lemma, in the plane of the pair (C_t(h), C_t(h+1)) of Thue-Morse correlations: the outer ellipse is the unit level of the form Q(x, y) = 20x^2 + 20xy + 41y^2, the blue ellipse its level 103/200, and the twelve orange curves are the images of the outer ellipse under the run operators A(t) = [[0, 1], [t, (t - 1)/3]] at t = (-1/2)^b for runs of length b = 1 to 12. Every orange curve lies inside the blue one, which is Lemma 5.4 on these twelve values; all of them pass through the same two points, the images of the outer ellipse's points on the line of e_1 = (-1/3, 1), which every A(t) sends to the line of e_0 = (1, -1/3), and as b grows they collapse onto a segment of that line, the image of the rank-one limit A(0). The curve for b = 2, the operator A(1/4), comes closest to the blue ellipse, where the certificate is tightest.](paper-design-generators)

A design is a rule on digits carried through every level by the Kronecker product, so a design read left to right is a stream of bits, and a stream of bits can be asked to pass for random. This paper asks which generators built from designs alone can pass, and answers with three no-go results and one positive one.

The setting is elementary. A *signed design* at base `b` is a stream `x(n) = prod_i v(n_i)`, the product over the base-`b` digits `n_i` of `n` of a sign tile `v : {0, ..., b-1} -> {-1, 1}`; Thue-Morse `t` is the base-`2` signed design of tile `(1, -1)`, and `t3` is the base-`3` signed design of tile `(1, 1, -1)`, the parity of the ternary digits equal to `2`. A *digit-restricted design* is the indicator that every digit of `n` lies in a fixed set, the Cantor set at base `3` with digits `{0, 2}`. A stream is tested by correlating it with a copy of itself moved by a map of the integers: a shift `n -> n + h`, a dilation `n -> Bn + r`, or a pair of powers of one base, `b^i n` against `b^j n`.

**One base.** Every `b`-automatic stream has two powers of `b` at which it coincides with itself term by term (Theorem 3.1), so one kernel test returns correlation `1` on every window where it is defined. For a design the coincidence sits at depth one: `x(bn + r) = +-x(n)` exactly (Proposition 3.2), so the argmax of a single dilation test names the base on a window of `9 . 2^15` terms (Fact 3.6). Signed designs also fail the shift test, with `C(bh) = C(h)` and some lag below `b` correlated (Theorem 3.4), but shift tests can be blind: Rudin-Shapiro has correlation `0` at every nonzero lag and dies under the dilation `n -> 2n` (Theorem 3.5).

**Two bases.** The moire XOR `z = t . t3` multiplies a base-`2` layer by a base-`3` layer. It is automatic in neither base (Theorem 4.1), its correlations factor, `C_z = C_t C3` with `C3 = C_(t3)` (Theorem 4.2), and the dilations at `B = 2, 4, 8, 16, 3, 9, 27` have mean `0` (Theorem 4.5). The main theorem is that every shift correlation tends to `0`:

```
C_z(h) -> 0   as  abs h -> infinity     (Theorem 5.11).
```

The ternary layer decays in balanced-ternary weight (Theorem 5.2). The binary layer decays in the number `r_2(h)` of runs of the binary word of `h`, `abs C_t(h) <= 1.0093 . 0.7177^(r_2(h))` (Theorem 5.5), because one quadratic form, the one drawn above, is contracted by every run of digits (Lemma 5.4); the constant per run lies in `[0.7118, 0.7177]`, and no bound in the number of binary ones can hold, since `C_t(2^k - 1) -> 1/9` (Proposition 5.6). An integer with few runs and few balanced-ternary digits is a short signed sum of powers of `2` equal to a short signed sum of powers of `3`, and the S-unit theorem allows that only finitely often (Theorem 5.10). In measure language, the spectral measure `sigma_z = sigma_t * sigma_3` is a continuous Rajchman measure although neither factor is Rajchman, each being invariant under multiplication by its own base (Corollary 5.12), and the same holds for any two measures whose coefficients decay in signed length at two multiplicatively independent bases (Proposition 5.13). That a convolution of two measures that are not Rajchman is Rajchman is elementary in itself: two Riesz products on spectra that meet only at `0` convolve to Lebesgue measure (Section 5). The content here is these two measures, whose coefficients are nonzero together at many lags and whose joint decay the S-unit theorem forces. The rate is effective but not explicit: `abs C_z(h) <= 1.0093 . 3^(-G(h)/4)` with `G(h) = log log h/(log log log h + kappa) - 1` and `kappa` the uncomputed constant of a digit-change theorem of Bugeaud, Cipu and Mignotte (Theorem 5.16).

What the second base does not remove is the base itself, or the determinism. A one-base kernel test still detects `z`: its raw z-score grows as `N^0.194` along `N = 2^k` (Theorem 4.6). A product of automatic streams at any bases satisfies an exact identity on a cube of `2^s` dilations, and for `z` the dilation square `z(n) z(2n) z(3n) z(6n)` is identically `1` (Theorem 6.1). Its factor complexity is at most quadratic, so its topological entropy is `0` (Theorem 6.2).

**Rounds.** Feed a design to an elementary cellular automaton and read a row or a column of the spacetime diagram. A bounded number of rounds of any local rule keeps every row automatic (Theorem 7.2); a rule linear over `F_p` keeps every column `p`-automatic on a `p`-automatic seed, at every number of rounds (Theorem 7.3); and a rule that agrees with one affine rule on every window its orbit meets, from some round on, keeps every row and column `2`-automatic on a `2`-automatic seed (Theorem 7.4). So every row of a design-fed elementary orbit stays in its seed's automatic class, and a column leaves it only through a base that is not a power of `2` or a nonlinearity exercised at unbounded depth (Corollary 7.5). The condition is necessary and not sufficient: rule `18` is nonlinear and draws the gasket of rule `90` from one cell, a known picture (Remark 7.9), and rules `182` and `22` exercise their nonlinearity at unbounded depth and draw `2`-automatic diagrams (Remark 7.10). Which elementary rules are linear on their reach is decided exactly, `90` of `256` from one cell and `20` from the Thue-Morse seed (Theorem 7.7).

**What is classical, and what is not.** Cobham's theorem is used as stated in the literature, as is the kernel characterisation of automatic sequences. The correlation recursion of a signed design (Lemma 2.5), the continuity of `sigma_t` and `sigma_3` (Theorem 5.8) and the white spectrum of Rudin-Shapiro (Theorem 3.5) are classical and reproved here, not claimed; their sources are Coquet, Kamae and Mendes France 1977, read at source, and Mahler 1927, Kakutani 1972, Rudin 1959 and Shapiro 1951, cited unread (Section 8). The idea of bounding the Fourier coefficients of a digit measure by a power of a constant below `1` counted over the digit patterns of the frequency, and of then using that an integer cannot be sparse in two multiplicatively independent bases, goes back to W. Schmidt 1960, read at source, and to Cassels 1959 on a problem of Steinhaus about normal numbers, cited unread: Schmidt's Lemma 2 bounds the product that controls the Fourier coefficient of his base-`s` Cantor measure at `w` by `a_2^z`, `a_2 < 1`, when `w` has at least `z` adjacent digit pairs of a certain kind, and his Theorem 1 gives continuum many numbers normal to base `r` and not simply normal to base `s` when no power of `r` is a power of `s`. Pollington 1981, read at source, carries the method to Hausdorff dimension, and Brown, Moran and Pearce 1985, unread, to Riesz products; the Thue-Morse spectral measure is a Riesz-type product, the weak-star limit of the densities `prod_(k < L) (1 - cos 2 pi 2^k x) = 2^(-L) abs(sum_(n < 2^L) t(n) e(nx))^2`, whose coefficient at `h` is `2^(-L) sum_(n, n + h < 2^L) t(n) t(n + h)`, since `prod_(k < L) (1 - e(2^k x)) = sum_(n < 2^L) t(n) e(nx)`. The S-unit theorem is van der Poorten and Schlickewei 1982 (unread) and Evertse 1984 (read at source), and the digit-change theorem is Bugeaud, Cipu and Mignotte 2013, Theorem 3.1, read at source. Linear automata over `F_p` from finite seeds are `p`-automatic by Allouche, von Haeseler, Peitgen and Skordev 1996, read at source; Theorem 7.3 runs the Frobenius identity behind their theorems on automatic seeds. The single-cell diagrams of elementary rules are catalogued by Wolfram 1983 and 2002, cited unread, and the pictures of rules `18`, `146` and `182` used in Section 7 are taken as known. What is new here, as far as the sources read show: the run lemma with its explicit constant `rho = (103/200)^(1/2)` and its quadratic-form certificate, the sharpness bracket `[0.711872, 0.717636]`, the signed-length lemma `BC(h, b) <= 2 l_b(h) - 1`, with `BC(h, b)` the number of digit changes of `h` in base `b` and `l_b(h)` the least number of signed powers of `b` with sum `h`, that turns the digit-change theorem into an effective rate for `C_z`, the three rungs as statements about generators built from designs, the reach criterion of Theorems 7.4 and 7.6, the exact closure of Theorem 7.6 and the counts `90` and `20` of Theorem 7.7. The literature on Thue-Morse correlations, from Mahler 1927 and Kakutani 1972 to Queffelec 1987 and Baake and Grimm 2013, was not searched, and the run lemma and its bracket may be there. The convolution statement itself, Proposition 5.13, is not claimed as new: no source read states it, and Brown, Moran and Pearce 1986, which could not be read, decomposes numbers into summands with prescribed normality properties at different bases, a setting in which such a convolution arises, and may contain it. Section 8 gives the details.

## Streams, designs and three tests

Throughout, a *stream* is a map `x : N -> {-1, 1}`, `N = {0, 1, 2, ...}`, and `e(u) = exp(2 pi i u)`. A stream of bits `a` is read as the stream `(-1)^a`, so the termwise product of two streams is the XOR of their bits.

**Definition 2.1 (correlations and tests).** The *correlation* of a stream at lag `h >= 0` is `C_x(h) = lim_(N -> infinity) (1/N) sum_(n < N) x(n) x(n + h)` when the limit exists, and `C_x(-h) = C_x(h)`. On a window of `N` pairs, with Pearson's correlation coefficient written `Pearson`:

- the *shift test* at lag `h` is `Pearson((x(n), x(n + h)) : n < N)`;
- the *dilation test* at `(B, r)`, `B >= 2`, `0 <= r < B`, is `Pearson((x(n), x(Bn + r)) : n < N)`;
- the *kernel test* at `(b, i, j)`, `i < j`, is `Pearson((x(b^i n), x(b^j n)) : n < N)`;
- the *dilation cube* at multipliers `P_1, ..., P_s` and modulus `M` is the mean over `n < N` of `prod_S x(M n prod_(i in S) P_i)`, the product over the `2^s` subsets `S` of `{1, ..., s}`; at `s = 2`, `P = (2, 3)`, `M = 1` it is the *dilation square* `x(n) x(2n) x(3n) x(6n)`.

For a stream of independent uniform signs each statistic has mean `0` and standard deviation `1/sqrt N` once the points it multiplies are distinct; the run of Fact 3.6 calls a value past `4/sqrt N` a failure, `0.0221` at `N = 2^15`. A test whose coordinate is constant on the window is not defined and is skipped.

**Definition 2.2 (designs).** Fix a base `b >= 2`. A *code* is a set `c` of digits and its *tile* is `v = (-1)^c`, `v(d) = -1` exactly when `d` is in `c`. The *signed design* of level `L` is `x_L(n) = prod_(i < L) v(n_i)` for `0 <= n < b^L`, the product over the `L` base-`b` digits `n_i` of `n`, leading zeros included. When `v(0) = 1` the levels agree and define one stream `x(n) = prod_i v(n_i)`; a code with `v(0) = -1` gives the stream of the tile `v(0) v` up to the global sign `v(0)^L`, so streams are taken with `v(0) = 1`. *Thue-Morse* is `t(n) = (-1)^(s_2(n))`, `s_2(n)` the number of ones in the binary word of `n`, the tile `(1, -1)` at base `2`, and `t3` is the tile `(1, 1, -1)` at base `3`, the parity of the number of ternary digits of `n` equal to `2`, with correlation `C3 = C_(t3)`. For a digit set `D` the *digit-restricted design* of level `L` is `A_L(n) = [every one of the L digits of n lies in D]`, the indicator that a Kronecker power renders, and its stream is `(-1)^(A_L(n))`; the Cantor set is `D = {0, 2}` at base `3`. *Rudin-Shapiro* is `r(n) = (-1)^(u(n))`, `u(n)` the number of adjacent pairs `11` in the binary word of `n`, so that `r(2n) = r(n)` and `r(2n + 1) = (-1)^n r(n)`.

**Definition 2.3 (kernel).** The *`b`-kernel* of a sequence `x` on `N` is the set of its subsequences `n -> x(b^k n + r)`, `k >= 0`, `0 <= r < b^k`, and `x` is *`b`-automatic* when its `b`-kernel is finite. Every signed design at base `b` is `b`-automatic, its kernel being `{x, -x}`, by Lemma 2.5 below. Rudin-Shapiro is `2`-automatic: with `r_1(n) = (-1)^n r(n)` its two rules give `r(2n) = r(n)`, `r(2n + 1) = r_1(n)`, `r_1(2n) = r(n)` and `r_1(2n + 1) = -r_1(n)`, so its kernel lies inside `{+-r, +-r_1}`.

That a sequence is `b`-automatic in the sense of Definition 2.3 exactly when a finite automaton reading base-`b` digits generates it is the kernel criterion, cited from Allouche and Shallit 2003 and not read here. It is used once, to apply Cobham's theorem in Theorem 4.1; every other statement of this paper is about kernels.

**Lemma 2.4 (closure).** Fix a base `b`.

1. A function of finitely many `b`-automatic sequences, termwise, is `b`-automatic.
2. If `x` is `b`-automatic and `x(m)` is read as `0` for `m < 0`, then `n -> x(n + j)` is `b`-automatic for every integer `j`.
3. A sequence that differs from a `b`-automatic one in finitely many terms is `b`-automatic, and every ultimately periodic sequence is `b`-automatic.
4. For `c >= 1`, a sequence is `b`-automatic exactly when it is `b^c`-automatic.

*Proof.* (1) The kernel element at `(k, r)` of `F(x_1, ..., x_m)` is `F` of the kernel elements at `(k, r)` of the `x_i`, so the kernel lies in a finite set of tuples. (2) Write `r + j = b^k i + r'` with `0 <= r' < b^k`. Then `x(b^k n + r + j) = g(n + i)` with `g = x(b^k . + r')` in the kernel of `x`, extended by `0` to negative arguments, and for `0 <= r < b^k` the integer `i` lies between `floor(j/b^k)` and `floor(j/b^k) + 1`; so every kernel element of the shift is some `n -> g(n + i)` with `g` in a finite set and `abs i <= abs j + 1`. (3) If `y = x` off `[0, n_0)`, a kernel element `y(b^k n + r)` with `b^k >= n_0` differs from `x(b^k n + r)` at `n = 0` at most, so it is one of finitely many sequences; the depths with `b^k < n_0` are finitely many. A periodic sequence of period `P` has `x(b^k n + r)` determined by `b^k mod P` and `r mod P`, finitely many choices, and an ultimately periodic one differs from a periodic one in finitely many terms. (4) The `b^c`-kernel lies inside the `b`-kernel. Conversely, a `b`-kernel element at depth `cq + s`, `0 <= s < c`, is `x(b^(cq+s) n + r) = g(b^s n + rho)` with `r = b^(cq) rho + r''`, `0 <= r'' < b^(cq)`, `0 <= rho < b^s`, and `g = x(b^(cq) . + r'')` in the `b^c`-kernel; `g`, `s` and `rho` range over finite sets. □

**Lemma 2.5 (the block recursion; classical, reproved here).** Let `x` be a signed design at base `b` with `v(0) = 1`. Then `x(bm + d) = v(d) x(m)` for `0 <= d < b`, `C_x(h)` exists at every lag, and with `al_s = sum_(r < b - s) v(r) v(r + s)`, so that `al_0 = b`,

```
C(bq + s) = (al_s C(q) + al_(b-s) C(q + 1))/b ,   0 <= s < b .
```

In particular `C(bq) = C(q)`, `C(1) = al_1/(b - al_(b-1))` and `C(s) = (al_s + al_(b-s) C(1))/b` for `1 <= s < b`.

*Proof.* The digits of `bm + d` are `d` followed by those of `m`, which gives the first identity, and iterating it, `x(u b^j + rho) = x(u) x(rho)` for `rho < b^j`. Fix `h` and `j` with `b^j > h`. If `rho = n mod b^j` satisfies `rho + h < b^j`, then `n` and `n + h` have the same quotient by `b^j` and `x(n) x(n + h) = x(rho) x(rho + h)`. So `x(n) x(n + h)` agrees with a `b^j`-periodic function outside a set of residues of density `h/b^j`, the upper and lower means of `x(n) x(n + h)` differ by at most `2h/b^j`, and letting `j` grow the mean exists. For the recursion write `S(N, h) = sum_(n < N) x(n) x(n + h)`, `n = bm + r` and `h = bq + s`. When `r + s < b`, `x(n) x(n + h) = v(r) v(r + s) x(m) x(m + q)`; when `r + s >= b`, it is `v(r) v(r + s - b) x(m) x(m + q + 1)`. Summing over `r` and `m < N`,

```
S(bN, bq + s) = al_s S(N, q) + al_(b-s) S(N, q + 1) ,
```

the second coefficient being `sum_(r' < s) v(r' + b - s) v(r') = al_(b-s)`. Divide by `bN` and let `N` grow. The two special values are the cases `q = 0`, with `C(0) = 1` and `al_(b-1) = v(0) v(b-1) = +-1 != b`. □

## One base: the leak

**Theorem 3.1 (the kernel coincidence).** Let `x` be `b`-automatic with a `b`-kernel of `m` elements. There are `0 <= e < e + d <= m` with `x(b^e n) = x(b^(e+d) n)` at every `n`. So the two coordinates of the kernel test at `(b, e, e + d)` agree term by term, and the test returns `1` on every window on which `n -> x(b^e n)` is nonconstant.

*Proof.* The `m + 1` sequences `n -> x(b^k n)`, `0 <= k <= m`, lie in the kernel, so two of them coincide. □

Nothing here is asymptotic. An independent uniform stream satisfies `x(b^e n) = x(b^(e+d) n)` on a window of `N` distinct pairs with probability `2^(-N)`. The proviso on the window is needed: for `(-1)^n` at base `2` the only coincidence is `x(2n) = x(4n) = 1`, a constant.

**Proposition 3.2 (designs leak at depth one).** For a signed design at base `b`, `x_L(bn + r) = v(0) v(r) x_L(n)` at every `n < b^(L-1)` and every `0 <= r < b`, so the dilation test at `(b, r)` returns exactly `v(0) v(r) = +-1` at every residue and on every window. The termwise product of signed designs at one base is the signed design of the product tile, so it leaks the same way.

*Proof.* The digits of `bn + r` at level `L` are `r` followed by the `L - 1` low digits of `n`, and since `n < b^(L-1)` the top digit of `n` at level `L` is a leading `0`. So `x_L(bn + r) = v(r) prod_(i < L-1) v(n_i)` and `x_L(n) = v(0) prod_(i < L-1) v(n_i)`. □

The generator checks the identity at every pair `(n, bn + r)`, `n < b^(L-1)`, for every base-`2` and base-`3` level-`1` code, `L = 20` at base `2` and `L = 13` at base `3`, `1048576` and `1594323` pairs, and finds `0` violations. For digit-restricted designs this paper states only the window readings of Fact 3.6.

**Remark 3.3 (depth zero cannot always be taken).** In Theorem 3.1 the coincidence need not start at `e = 0`. For `x_1(n) = (-1)^n t(n)` and `d >= 1`, `x_1(2^d n + r) = (-1)^r t(r) t(n)`, so `x_1(n) x_1(2^d n + r) = +-(-1)^n` and every two-point dilation of `x_1` at a power of `2` is balanced to the last term: the generator prints `max abs Pearson = 0.000000` over `d <= 3`, `r < 2^d`, `n < 2^20`. Yet `x_1(2n) = t(n) = x_1(4n)`, and the kernel pair `(2n, 4n)` prints `1.000000`. The shift test catches `x_1` anyway, since `C_(x_1)(h) = (-1)^h C_t(h)`.

**Theorem 3.4 (no signed design is white).** Let `x` be a signed design with `v(0) = 1` and `v` not constant. Then `C(bh) = C(h)` at every `h`, and some lag `1 <= s < b` has `C(s) != 0`; so `C(b^k s) = C(s) != 0` for every `k`.

*Proof.* `C(bh) = C(h)` is Lemma 2.5 at `s = 0`. At `b = 2`, `al_1 = v(0) v(1) = -1` and `C(1) = -1/(2 + 1) = -1/3`. At `b >= 3`, if `C(1) = 0` then Lemma 2.5 gives `C(b - 1) = al_(b-1)/b = v(0) v(b - 1)/b != 0`, and `b - 1 >= 2`. □

The exact values, from Lemma 2.5 in rational arithmetic, for every level-`1` code at bases `2` and `3`; a code is written as its digit set, two codes in one row having the same tile up to the sign `v(0)`:

| base | codes | class | dilation `v(0) v(r)`, `r = 0, ..., b-1` | `C(1), ..., C(b)` | `max abs C(h)`, `1 <= h <= 4096` |
|---|---|---|---|---|---|
| 2 | `{}`, `{0, 1}` | constant | `1, 1` | `1, 1` | `1` |
| 2 | `{0}`, `{1}` | Thue-Morse | `1, -1` | `-1/3, -1/3` | `1/3` at `h = 1` |
| 3 | `{}`, `{0, 1, 2}` | constant | `1, 1, 1` | `1, 1, 1` | `1` |
| 3 | `{1}`, `{0, 2}` | `(-1)^n` | `1, -1, 1` | `-1, 1, -1` | `1` |
| 3 | `{2}`, `{0, 1}` | parity of the `2`s | `1, 1, -1` | `0, -1/3, 0` | `1/3` at `h = 2` |
| 3 | `{0}`, `{1, 2}` | parity of the nonzero digits | `1, -1, -1` | `0, -1/3, 0` | `1/3` at `h = 2` |

The base-`3` code `{1}` is periodic: `3` is odd, so the parity of the ternary digits equal to `1` is the parity of `n`. The empirical `C(h)` of the level-`19` and level-`12` renders, `N = 524288` and `531441`, agrees with the recursion within `1.02e-5` and `5.65e-6` at every `h <= 8`. In measure language, the spectral measure of a signed design, the measure `sigma_x` on the circle with `hat sigma_x(h) = C_x(h)` (Section 5), satisfies `hat sigma_x(bh) = hat sigma_x(h)`: it is invariant under `u -> bu` and is never Rajchman. That invariance under one integer ratio is what the second base will break.

**Theorem 3.5 (Rudin-Shapiro is white; classical, reproved here).** `sum_(n < 2^k) r(n) r(n + h) = 0` for every `k >= 2` and `1 <= h <= 2^(k-2)`; for every `N` and `h >= 1`, `abs(sum_(n < N) r(n) r(n + h)) <= 4h (log2 N + 1)`. So `C_r(h) = 0` at every `h >= 1`, and the spectral measure of `r` is Lebesgue measure. Yet `r(2n) = r(n)`, and the dilation test at `(2, 0)` returns `1`.

*Proof.* Write `S(N, h) = sum_(n < N) r(n) r(n + h)` and `D(N, h) = sum_(n < N) (-1)^n r(n) r(n + h)`. Splitting `n` by parity with the two rules gives four exact identities:

```
S(2N, 2h)     = (1 + (-1)^h) S(N, h)
S(2N, 2h + 1) = (-1)^h D(N, h) + D(N, h + 1)
D(2N, 2h)     = (1 - (-1)^h) S(N, h)
D(2N, 2h + 1) = (-1)^h D(N, h) - D(N, h + 1)
```

If `S(N, j) = 0` for `1 <= j <= J` and `D(N, j) = 0` for `0 <= j <= J`, the identities give the same at `2N` with `2J` in place of `J`, the lag `0` entering only through `D(2N, 0)`, whose coefficient `1 - (-1)^0` vanishes. At `N = 4`, where `r(0), ..., r(4) = 1, 1, 1, -1, 1`, the sums `S(4, 1)`, `D(4, 0)` and `D(4, 1)` are `0`, so `J = 1` at `k = 2` and `J = 2^(k-2)` at `k`.

For a general `N`, cut `[0, N)` into at most `log2 N + 1` dyadic blocks `[u 2^j, (u + 1) 2^j)`. On a block, `r(u 2^j + m) = r(u) r(m) (-1)^(u_0 m_(j-1))` for `m < 2^j`, since the pairs `11` of `u 2^j + m` are those of `u`, those of `m` and possibly the pair of the last digit `u_0` of `u` with the top digit `m_(j-1)` of `m`. So the pairs `(n, n + h)` inside the block sum to the aperiodic autocorrelation `A_P(h)` of `P_j = (r(m))_(m < 2^j)` when `u` is even and to `A_Q(h)` of `Q_j = (r(m) (-1)^(m_(j-1)))_(m < 2^j)` when `u` is odd, and at most `h` pairs leave the block. The blocks obey `P_(j+1) = (P_j, Q_j)` and `Q_(j+1) = (P_j, -Q_j)`, so on the unit circle `abs(P_(j+1)(z))^2 + abs(Q_(j+1)(z))^2 = 2 (abs(P_j(z))^2 + abs(Q_j(z))^2)`, which from `P_1 = (1, 1)`, `Q_1 = (1, -1)` is the constant `2^(j+1)`; hence `A_P(h) + A_Q(h) = 0` at every `h != 0`. Since `S(2^j, h)` is `A_P(h)` plus at most `h` pairs, `abs(A_Q(h)) = abs(A_P(h)) <= h` once `h <= 2^(j-2)`, and a block with `2^j < 4h` holds fewer than `4h` pairs. Each block therefore contributes at most `4h`. Dividing by `N` gives `C_r(h) = 0`, and a measure with all nonzero coefficients `0` is Lebesgue measure. The last sentence is the rule `r(2n) = r(n)`. □

The generator finds `S(2^k, h) = 0` exactly for every `h <= 64` and `k = 10, ..., 20`; the largest `abs(S(N, h))` over all `N <= 2^20` and `h <= 64` is `118`, at `h = 59`, against `sqrt(2^20) = 1024`; and `r(2n) = r(n)` and `r(4n + 1) = r(n)` hold at all `2^20` values of `n`. One dilation reads what no fixed lag can. A lag comparable to the sample can as well: `r(n + 2^k) = r(n)` for `n < 2^(k-1)`, so the shift test at lag `2N` on a window of `N` returns `1`, a dilation in disguise.

**Fact 3.6 (the forensic reading).** At `N = 2^15` pairs, the largest dilation test over `2 <= B <= 9` and `0 <= r < B` is exactly `1.0000`, attained at the stream's own base, for Thue-Morse at `(2, 0)`, Rudin-Shapiro at `(2, 0)`, the Cantor stream at `(3, 0)` and the Sierpinski carpet read in Morton order at `(9, 0)`; for the centre column of Rule 30 from one cell and for the `splitmix64` control it is `0.0095` and `0.0093`, inside the band `0.0221`. The Morton carpet is the base-`9` digit-restricted design that drops the digit `4`, the centre pair `(1, 1)`. Generator: `lab/rs/design-ladder`, verb `battery`, whose Rule 30 column is checked against the first `36` terms of A051023.

So the argmax of one dilation test reads the base of a one-base generator off a window of `9 . 2^15` terms.

## Two bases: the leak diluted

The *moire XOR* multiplies a base-`2` layer by a base-`3` layer. Two are studied: `z = t . t3`, Thue-Morse times the ternary parity, and `y = r . t3`, Rudin-Shapiro times the same. Both factors of `z` satisfy `abs C(h) <= 1/3` at every `h >= 1`: by Lemma 2.5, `C_t(2q) = C_t(q)` and `C_t(2q + 1) = -(C_t(q) + C_t(q + 1))/2` with `C_t(1) = -1/3`, and `C3(3q) = C3(q)`, `C3(3q + 1) = -C3(q + 1)/3`, `C3(3q + 2) = -C3(q)/3` with `C3(0) = 1`, `C3(1) = 0`, and induction on `h` closes both bounds.

**Theorem 4.1 (two bases leave both automatic classes).** `z` and `y` are neither `2`-automatic nor `3`-automatic.

*Proof.* An ultimately periodic stream of eventual period `P` has `C(P) = 1`. If `z` were `2`-automatic, then `t3 = z . t` would be `2`-automatic by Lemma 2.4, and it is `3`-automatic; by Cobham's theorem, which needs the kernel criterion to be read in the sense of Definition 2.3, a stream automatic in the two multiplicatively independent bases `2` and `3` is ultimately periodic, against `abs C3(h) <= 1/3` at `h >= 1`. If `z` were `3`-automatic, `t = z . t3` would be ultimately periodic in the same way, against `abs C_t(h) <= 1/3`. For `y` the same two steps run with `C_r(h) = 0` (Theorem 3.5) in place of `C_t`. □

Whether `z` or `y` is automatic in another base is not decided here (Question 9.7).

**Theorem 4.2 (correlations factor at coprime bases).** Let `x` be a signed design at base `b` with `v(0) = 1`, or Rudin-Shapiro at `b = 2`, and let `w` be a signed design at a base `b'` with `gcd(b, b') = 1` and `v'(0) = 1`. Then `C_(xw)(h)` exists and equals `C_x(h) C_w(h)` at every lag.

*Proof.* By the proof of Lemma 2.5, `F(n) = x(n) x(n + h)` agrees with a `b^j`-periodic function `F_j` outside a set of residues of density `h/b^j`. Rudin-Shapiro does the same with period `2^(j+1)`: for `n = u 2^j + rho` with `rho + h < 2^j`, `r(n) r(n + h) = r(rho) r(rho + h) (-1)^(u_0 (rho_(j-1) + (rho + h)_(j-1)))` depends only on `n mod 2^(j+1)`. Likewise `G(n) = w(n) w(n + h)` agrees with a `b'^l`-periodic `G_l` outside density `h/b'^l`. Since the periods are coprime, the Chinese remainder theorem makes the residues of `n` modulo the two periods jointly uniform over a full period of the product, so the mean of `F_j G_l` is the product of the means of `F_j` and `G_l`. The upper and lower means of `FG` are therefore within `2h/b^j + 2h/b'^l` of the product of the means of `F_j` and `G_l`, each of which is within `2h/b^j` or `2h/b'^l` of `C_x(h)` or `C_w(h)`. Let `j` and `l` grow. □

**Corollary 4.3 (the shift test catches `z` and not `y`).** `C_z(h) = C_t(h) C3(h)`, so `C_z(2) = 1/9` and `C_z(6) = -1/9`, and `C_z(1), ..., C_z(12)` are `0, 1/9, 0, -1/27, 0, -1/9, 0, -1/27, 0, 0, 0, 1/27`. And `C_y(h) = C_r(h) C3(h) = 0` at every `h >= 1`, so the spectral measure of `y` is Lebesgue measure.

The empirical `C_z(h)` at `N = 2^22` agrees with the product within `7.97e-4` at every `h <= 256`; the empirical `abs C_y(h)` over `h <= 256` at `N = 2^22` is at most `6.96e-4`, against `1/sqrt N = 4.88e-4`. Generator: `lab/rs/design-ladder`, verb `moire`.

A dilation at one base strips that layer exactly. For `r < 2^a`, `t(2^a n + r) = t(r) t(n)`, so `z(n) z(2^a n + r) = t(r) . t3(n) t3(2^a n + r)`, a test on the ternary layer alone; symmetrically `z(n) z(3^c n + r) = t3(r) . t(n) t(3^c n + r)` for `r < 3^c`. What is left is a two-point dilation of one signed design at a foreign multiplier, and those sums obey a finite linear recursion.

**Lemma 4.4 (the state recursion).** Let `x` be a signed design at base `b` with `v(0) = 1`, `A, B >= 1`, and for `0 <= a < A`, `0 <= r < B` put `S_(a,r)(N) = sum_(n < N) x(An + a) x(Bn + r)`. There is an integer matrix `M` on the `AB` states `(a, r)` with `S(bN) = M S(N)` and `S(bN + e) = M S(N) + c` with `abs c <= b - 1` in every coordinate, `0 <= e < b`. If the spectral radius of `M` is `lambda < b`, every mean `S_(a,r)(N)/N` tends to `0`, and `S(N) = O(N^(theta + epsilon))` for every `epsilon > 0`, `theta = log_b max(lambda, 1)`.

*Proof.* Split `n = bm + d`, `0 <= d < b`, and write `Ad + a = b q_1 + s_1` and `Bd + r = b q_2 + s_2` with `0 <= s_1, s_2 < b`; then `0 <= q_1 < A` and `0 <= q_2 < B`, and by Lemma 2.5 `x(An + a) x(Bn + r) = v(s_1) v(s_2) x(A m + q_1) x(B m + q_2)`. Summing over `m < N` and `d < b` gives `S(bN) = M S(N)` with `M[(a, r), (q_1, q_2)]` the sum of `v(s_1) v(s_2)` over the digits `d` that lead there. For `bN + e` the digits `d < e` carry one extra term `m = N`, each of modulus `1`. Unrolling along the base-`b` digits of `N`, `S(N) = sum_(i <= L) M^i c_i` with `L = floor(log_b N)` and `abs c_i <= b - 1`, so `max abs S(N) <= (b - 1) sum_(i <= L) norm(M^i)`, and `norm(M^i) = O((lambda + epsilon)^i)`. □

**Theorem 4.5 (one-base dilations of `z` have mean `0`).** For `B = 2, 4, 8, 16, 3, 9, 27` and every `0 <= r < B`, `(1/N) sum_(n < N) z(n) z(Bn + r) -> 0`.

*Proof (by exact computation).* By the stripping identities, the sums are `t(r)` times the sums of Lemma 4.4 for `t3` at `A = 1`, `B = 2^a` and base `3`, or `t3(r)` times those for `t` at `A = 1`, `B = 3^c` and base `2`. The generator builds each `M`, computes `M^64` exactly in `128`-bit integers with every product and sum checked for overflow, and bounds the spectral radius by `norm(M^64)^(1/64)`, the maximum row sum, the root taken in floating point and enlarged by a relative `10^(-9)`. It prints `lambda <= 1.0789, 1.7503, 1.7871, 1.8011` at `B = 2, 4, 8, 16` on the ternary layer, each below `3`, and `lambda <= 1.4311, 1.4997, 1.4607` at `B = 3, 9, 27` on the binary layer, each below `2`. Lemma 4.4 closes. □

The exponents are `theta <= 0.0692, 0.5096, 0.5285, 0.5356` and `0.5171, 0.5847, 0.5466`. Only the bound at `B = 2` lies below `1/2`, so at the other six multipliers Lemma 4.4 does not exclude that these sums grow faster than `sqrt N`: Theorem 4.5 is a statement about means, not about detection. The battery at `N = 2^15` sees no excess: the largest dilation test of `z` over `2 <= B <= 9` is `0.0111`, at `(6, 2)`. Generator: `lab/rs/design-ladder`, verbs `moire` and `battery`.

The kernel pairs are sharper, because a pair of powers of one base strips that layer and leaves the other layer along two foreign multipliers. Since `t3(3n) = t3(n) = t3(9n)`, `z(3n) z(9n) = t(3n) t(9n)`.

**Theorem 4.6 (the kernel pair at base `3` diverges).** The sums `S_k = sum_(n < 2^k) t(3n) t(9n)` satisfy `S_(k+4) = 2 S_(k+3) - 2 S_(k+2) + S_(k+1) + 2 S_k` at every `k >= 0`, so `S_k = c_1 phi^k + c_2 psi^k + O(2^(k/2))` with `phi` the golden ratio, `psi = -1/phi` and `c_1 != 0`. So the raw kernel statistic `S/sqrt N` of `z` at `(3, 1, 2)` grows as `N^(log2 phi - 1/2) = N^(0.1942...)` along `N = 2^k`. The Pearson form also subtracts the means of `z(3n)` and `z(9n)`, which this theorem does not bound.

*Proof (by exact computation).* Lemma 4.4 at `A = 3`, `B = 9`, base `2`, gives `27` states and `S(2^k) = M^k S(1)`, so `S_k` satisfies the linear recurrence of the characteristic polynomial of `M`, of order `27`, and so does `U_k = S_(k+4) - 2 S_(k+3) + 2 S_(k+2) - S_(k+1) - 2 S_k`. The generator checks `U_k = 0` at every `k <= 60`; `27` consecutive zeros of a sequence obeying a monic recurrence of order `27` force it to vanish. The characteristic polynomial of the order-`4` law is `(x^2 - x - 1)(x^2 - x + 2)`, whose second factor has roots of modulus `sqrt 2`. Applying `E^2 - E + 2`, `E` the shift, `T_k = S_(k+2) - S_(k+1) + 2 S_k = 3 c_1 phi^k + 3 c_2 psi^k`, since `phi^2 - phi + 2 = psi^2 - psi + 2 = 3`; so `T` is an integer sequence with `T_(k+2) = T_(k+1) + T_k`. From `S_0, ..., S_3 = 1, 2, 4, 6`, `T_0 = 4` and `T_1 = 6`, so `T` is positive and increasing, while `c_1 = 0` would make `T_k = 3 c_2 psi^k -> 0`. □

The generator prints `S_60/phi^60 = 1.263038` and a growth exponent `0.694232` against `log2 phi = 0.694242`. From these two printed constants the predicted raw value `1.263038 (2^15)^(log2 phi - 1)` is `0.0526`; the exact one is `S_15/2^15 = 1662/32768 = 0.0507`, by the order-`4` law from `S_0, ..., S_3 = 1, 2, 4, 6`, and the battery reads the Pearson value `0.0508` for `z`. So a one-base kernel test still detects `z`, with significance growing as `N^0.194`: the second base dilutes this leak and does not remove it.

## Non-Rajchman by Design, Rajchman by Convolution

A correlation is positive definite: for lags `a_1, ..., a_m` and complex `c_i`, `sum_(i,j) c_i conj(c_j) C_x(a_j - a_i) = lim (1/N) sum_(n < N) abs(sum_i c_i x(n + a_i))^2 >= 0`, since moving the window by `max a_i` changes each average by `O(max a_i/N)`. By Herglotz's theorem, cited through its encyclopedia statement, there is a probability measure `sigma_x` on the circle `T = R/Z` with `hat sigma_x(h) = int e(-hu) d sigma_x(u) = C_x(h)`, the *spectral measure* of `x`. A measure is determined by its coefficients, so Theorem 4.2 reads `sigma_z = sigma_t * sigma_3`, `sigma_3` the spectral measure of `t3`. A finite measure on `T` is *Rajchman* when `hat sigma(h) -> 0` as `abs h -> infinity`. Neither factor is: `C_t(2^k) = C_t(1) = -1/3` and `C3(2 . 3^k) = C3(2) = -1/3` for every `k`, the invariance of Theorem 3.4. This section proves that their convolution is. That such a convolution exists is elementary (the Riesz pair after Corollary 5.12); the content is that these two measures, each invariant under its own base and with coefficients nonzero together at many lags, are forced to decay jointly by the S-unit theorem.

Three counts of sparsity, for `h >= 1`:

- `r_2(h)` is the number of maximal blocks of equal digits in the binary word of `h`, leading zeros excluded: `r_2(2^k - 1) = 1`, `r_2(2^k) = 2`;
- `w_3(h)` is the number of nonzero digits of `h` in balanced ternary, digits `-1, 0, 1`;
- `l_b(h)` is the least number of terms `+-b^a`, `a >= 0`, with sum `h`, repetitions allowed; `l_b(-h) = l_b(h)`.

**Lemma 5.1 (signed lengths).** `l_2(h) <= r_2(h) + 1` and `l_3(h) <= w_3(h)`.

*Proof.* A run of ones on the bits `i` to `j - 1` is `2^j - 2^i`, and a word with `r` runs has at most `ceil(r/2)` runs of ones, so `h` is a signed sum of at most `2 ceil(r/2) <= r + 1` powers of `2`. Balanced ternary writes `h` with `w_3(h)` signed powers of `3`. □

**Theorem 5.2 (the ternary layer decays in balanced weight).** `abs C3(h) <= 3^(-w_3(h)/2)` at every `h >= 1`; hence `abs C_z(h) <= 3^(-1 - w_3(h)/2)`.

*Proof.* Write `h = 3q + delta` with `delta` in `{-1, 0, 1}`. Lemma 2.5 at the tile `(1, 1, -1)` gives `C3(3q) = C3(q)` and `C3(3q + delta) = -C3(q + delta)/3` for `delta = +-1`, the case `delta = -1` being `C3(3(q - 1) + 2) = -C3(q - 1)/3`. Adding `+-1` to a balanced-ternary word changes its weight by exactly one: the carry runs through the low digits equal to the added sign, turning each into its negative, and stops at the first other digit, which moves by one toward or away from `0`. So `w_3(q + delta) >= w_3(q) - 1 = w_3(h) - 2` when `delta != 0`, and induction on `h` from `C3(0) = 1`, `C3(1) = 0` closes. The second bound adds `abs C_t(h) <= 1/3`. □

The bound is attained: the generator prints the maximum of `abs C3` over `h < 2^20` at each even weight `2, 4, 6, 8, 10` as `1/3, 1/9, 1/27, 1/81, 1/243`, and `0` at every odd weight, where `h` is odd and `C3(h) = 0`.

### The run lemma

Write `vec(h) = (C_t(h), C_t(h + 1))`. By Lemma 2.5, `vec(2h) = M_0 vec(h)` and `vec(2h + 1) = M_1 vec(h)` with

```
M_0 = [[1, 0], [-1/2, -1/2]] ,    M_1 = [[-1/2, -1/2], [0, 1]] ,
```

and `vec(0) = e_0 = (1, -1/3)` is fixed by `M_0`, so leading zeros cost nothing and reading the binary digits `h_(L-1) ... h_0` of `h` from the top gives `vec(h) = M_(h_0) M_(h_1) ... M_(h_(L-1)) e_0`. `M_0` has eigenvalue `1` on `e_0` and `-1/2` on `(0, 1)`; `M_1` has `1` on `e_1 = (-1/3, 1)` and `-1/2` on `(1, 0)`. Inside a run of `k` equal digits `d`, `M_d^k = P_d + (-1/2)^k (I - P_d)` with `P_d` the projection on `e_d` along the other eigenvector, so a long run forgets everything but its `e_d` component, and nothing is lost inside a run. Decay is counted in runs. With the swap `J(x, y) = (y, x)`, `J M_0 J = M_1`, and with `t = (-1/2)^b`, `M_0^b = [[1, 0], [(t - 1)/3, t]]`. The *run operator* is

```
A(t) = M_0^b J = [[0, 1], [t, (t - 1)/3]] .
```

**Lemma 5.3 (run factorisation).** If the binary word of `h >= 1` has runs of lengths `b_1, ..., b_r` from the top, then `vec(h) = J^(r mod 2) A(t_r) ... A(t_1) e_0` with `t_i = (-1/2)^(b_i)`. So `C_t(h)` is the first coordinate of `u = A(t_r) ... A(t_1) e_0` when `r` is even and its second when `r` is odd, and every `t_i` lies in `[-1/2, 1/4]`.

*Proof.* `M_0^b = [[1, 0], [c_b, t]]` with `c_1 = -1/2` and `c_(b+1) = -1/2 - c_b/2`, which `c_b = (t - 1)/3` solves. The top run is ones, and `M_1^(b_1) e_0 = J M_0^(b_1) J e_0 = J A(t_1) e_0`. A run of zeros after a leading `J` multiplies by `M_0^(b)` and absorbs the `J`, giving `A(t) (...)`; a run of ones multiplies by `J M_0^b J` and leaves a `J` in front. Induction on `r`. The values: `b = 1` gives `-1/2`, `b = 2` gives `1/4`, and the rest lie between. □

At `h = 1`, `A(-1/2) e_0 = (-1/3, -1/3)`, and its second coordinate is `C_t(1) = -1/3`.

**Lemma 5.4 (the certificate).** Let `Q(x, y) = 20x^2 + 20xy + 41y^2`. Then `Q(A(t) u) <= (103/200) Q(u)` for every `u` in `R^2` and every `t` in `[-1/2, 1/4]`.

*Proof.* Since `A(t) u = (y, tx + (t - 1) y/3)`, the difference `(103/200) Q(u) - Q(A(t) u)` is the binary form `a(t) x^2 + B(t) xy + d(t) y^2` with

```
a(t)     = 103/10 - 41 t^2
B(t)     = 103/10 + (22/3) t - (82/3) t^2
d(t)     = 5807/1800 + (22/9) t - (41/9) t^2
Delta(t) = 4 a(t) d(t) - B(t)^2 = 30179/1125 - (2266/45) t - (93367/450) t^2 ,
```

the quartic and cubic terms of `Delta` cancelling. The form is positive semidefinite when `a >= 0`, `d >= 0` and `Delta >= 0`. All three are concave quadratics in `t`, so on an interval each is least at an endpoint, and the endpoint values are positive: `a = 1/20` and `619/80`, `d = 173/200` and `1421/400`, `Delta = 133/1000` and `5077/4000` at `t = -1/2` and `t = 1/4`. □

The generator recomputes the three polynomials in exact rationals, checks their concavity and the six values, and counts their real roots on the closed interval by an exact Sturm sequence: none. The form was found by a numerical scan, and nothing in the proof depends on how it was found. The figure draws the lemma on the twelve run lengths `b = 1, ..., 12`.

**Theorem 5.5 (the run lemma).** `abs C_t(h) <= K rho^(r_2(h))` at every `h >= 1`, with `rho = (103/200)^(1/2) = 0.717635...` and `K = (41 . 161/720)^(1/2)/3 = 1.009293...`; rounded safe, `abs C_t(h) <= 1.0093 . 0.7177^(r_2(h))`.

*Proof.* Let `N(u) = Q(u)^(1/2)`. By Lemmas 5.3 and 5.4, `N(u) <= rho^r N(e_0)` with `N(e_0)^2 = Q(1, -1/3) = 161/9`. The matrix of `Q` is `[[20, 10], [10, 41]]`, of determinant `720`, and the diagonal of its inverse is `41/720` and `20/720`, so `abs x <= (41/720)^(1/2) N(u)` and `abs y <= (20/720)^(1/2) N(u)`. That gives `K` when `r` is even and `0.704921` when `r` is odd. □

Since `C_t(2h) = C_t(h)`, the bound also holds with the runs of the odd part of `h`, one fewer when `h` is even. The generator checks the theorem at every `h < 2^25` with `0` failures, the largest ratio of `abs C_t(h)` to the bound being `0.9067`, at `h = 90`.

**Proposition 5.6 (sharpness, and runs rather than ones).**

1. No `rho' < ((3 + sqrt 73)/32)^(1/3) = 0.711872...` satisfies `abs C_t(h) <= K' rho'^(r_2(h))` for any constant `K'`. So the best constant per run lies in `[0.7118, 0.7177]`.
2. `C_t(2^k - 1) = (-1/2)^k + (1 - (-1/2)^k)/9` at every `k >= 1`, so `C_t(2^k - 1) -> 1/9`, and no bound that decays in the number of binary ones of `h` can hold.

*Proof.* (1) For the run pattern `(1, 1, 2)` repeated, the operator over three runs is `A(1/4) A(-1/2) A(-1/2) = [[1/4, -1/4], [-3/16, -1/16]]`, with characteristic polynomial `x^2 - (3/16) x - 1/16` and eigenvalues `(3 +- sqrt 73)/32`, the larger `0.360750...`. Its eigenvectors are `(1, (5 - sqrt 73)/8)` for the larger root and `(1, (5 + sqrt 73)/8)` for the smaller; `e_0 = (1, -1/3)` is not parallel to the second, so it has a nonzero component on the first, whose coordinates are both nonzero. Along `r = 3j` runs either coordinate of the `j`-th power applied to `e_0` is therefore of exact order `0.360750^j`, which a bound `K' rho'^(3j)` with `rho'^3 < 0.360750` contradicts. (2) The word `1^k` is one run, and `vec(2^k - 1) = M_1^k e_0 = J M_0^k J e_0`; with `J e_0 = (-1/3, 1)` and `M_0^k = [[1, 0], [(t - 1)/3, t]]`, `t = (-1/2)^k`, the first coordinate is `t + (1 - t)/9`. □

The generator reads the family of (1) exactly to `30` runs, where at `h = 776125854900` the value is `3.5415e-5` and the root over six runs is `0.711858`, and checks the closed form of (2) at `k = 1, 2, 3, 4, 6, 8, 12, 16, 24, 40, 60`. Part (2) fixes the shape of everything that follows: the finiteness that closes the argument must be stated in signed lengths, where `2^k - 1 = 2^k - 2^0` has length `2`, and not in digit sums.

### Continuity

**Lemma 5.7 (Wiener's identity).** For a finite positive measure `sigma` on `T`, `(1/(2H + 1)) sum_(abs h <= H) abs(hat sigma(h))^2 -> sum_u sigma({u})^2`.

*Proof.* `abs(hat sigma(h))^2 = int int e(-h(u - u')) d sigma(u) d sigma(u')`, and the average of `e(-hs)` over `abs h <= H` is a function `D_H(s)` with `abs D_H <= 1`, `D_H(0) = 1` and `D_H(s) -> 0` at every `s != 0` in `T`. Dominated convergence on `sigma x sigma` gives the mass of the diagonal. □

**Theorem 5.8 (Cesaro decay).** At every `H >= 2`, `(1/H) sum_(1 <= h < H) C_t(h)^2 <= 2 K^2 H^(-theta)` with `theta = log2(400/303) = 0.40068...`, and `(1/H) sum_(1 <= h <= H) C3(h)^2 <= (25/9) H^(-log3(9/5))` with `log3(9/5) = 0.53502...`. So `sigma_t`, `sigma_3` and `sigma_z` have no atoms; for `sigma_t` and `sigma_3` this is classical, reproved here (Section 8).

*Proof.* The words of length `m` with `r` runs number `binom(m - 1, r - 1)`, so `sum_(2^(m-1) <= h < 2^m) rho^(2 r_2(h)) = rho^2 (1 + rho^2)^(m-1)` and `sum_(1 <= h < 2^L) rho^(2 r_2(h)) = (1 + rho^2)^L - 1`. With `2^(L-1) < H <= 2^L`, Theorem 5.5 gives `(1/H) sum C_t^2 <= (K^2/H)(1 + rho^2)^L <= 2 K^2 ((1 + rho^2)/2)^L <= 2 K^2 H^(-theta)`, since `(1 + rho^2)/2 = 303/400`. On the ternary side, the integers `abs h <= (3^k - 1)/2` are exactly those with `k` balanced digits, so `sum 3^(-w_3(h)) = (5/3)^k` over them; take `k` least with `H <= (3^k - 1)/2`, so that `3^k < 9H`, and Theorem 5.2 gives `(1/H) sum C3^2 <= (5/3)^k/H < (25/9) H^(log3(5/3) - 1)`. Lemma 5.7 then excludes atoms, for `sigma_z` too since `abs C_z <= abs C_t`. □

The measured slopes of `log2` of the means over `1 <= h < 2^L`, from `L = 9` to `L = 25`, are `-0.6324` per bit for `C_t^2` and `-0.7882` for `C3^2`, against the proved `-0.4006` and `-0.5350`, rounded safe from `theta = 0.400682` and `log3(9/5) = 0.535026`. Generator: `lab/rs/two-bases`, verb `decay`.

### Sparse in two bases

**Theorem 5.9 (van der Poorten and Schlickewei 1982, Evertse 1984).** Let `Gamma` be a finitely generated subgroup of `Q^*` and `a_1, ..., a_n` nonzero rationals. Then `a_1 u_1 + ... + a_n u_n = 1` has only finitely many solutions `(u_1, ..., u_n)` in `Gamma^n` with no vanishing subsum on the left.

This is the rational case of the statement of the survey of Evertse, Gyory and Stewart 2018, read at source, which credits it to van der Poorten and Schlickewei and, independently, to Evertse; it was checked against Theorem 1 of Evertse 1984 in its Numdam scan, and the report of van der Poorten and Schlickewei is unread. The proof runs through Schlickewei's `p`-adic subspace theorem and is not effective.

**Theorem 5.10 (sparse in two bases).** Let `a, b >= 2` be multiplicatively independent. For every `n` only finitely many `h >= 1` have `l_a(h) + l_b(h) <= n`. In particular every set `{h : r_2(h) <= R, w_3(h) <= W}` is finite.

*Proof.* Write `h = sum_(i in I) x_i = sum_(j in J) y_j` with `x_i = +-a^(alpha_i)`, `y_j = +-b^(beta_j)`, `abs I = l_a(h)` and `abs J = l_b(h)`, and split the zero sum `sum x_i - sum y_j = 0` into minimal vanishing blocks.

- A block inside `I` alone could be deleted from a representation of length `l_a(h)`, and if it were all of `I` then `h = 0`; so no block lies in `I` alone, and likewise none in `J` alone.
- A block of two terms `x = y` forces `+-a^alpha = +-b^beta`, so `alpha = beta = 0` by independence.
- A larger block, divided by `-y_(j_0)` for one of its `y`-terms, becomes a solution of `u_1 + ... + u_m = 1` in `Gamma = <-1, a, b>` with no vanishing subsum, by minimality, and `m < n`; Theorem 5.9 with every coefficient `1` leaves finitely many such tuples for each `m`.
- In a block's tuple, `x_i/y_(j_0) = +-a^(alpha_i) b^(-beta_(j_0))`, and `a^alpha b^(-beta) = a^(alpha') b^(-beta')` forces `alpha = alpha'` and `beta = beta'` by independence; `y_j/y_(j_0)` fixes `beta_j - beta_(j_0)` the same way.

So every exponent ranges over a finite set as `h` does, and `h <= n a^(max alpha_i)` is bounded. The second sentence follows by Lemma 5.1. □

This contains the finiteness half of Senge and Straus 1973, cited through Bugeaud, Cipu and Mignotte: digit sums at most `c` in both bases give `l_a, l_b <= c`. The converse fails, since `2^k - 1` has `l_2 = 2` and binary digit sum `k`, and the lags where `C_t` stays large lie only in the wider set.

### The limit

**Theorem 5.11 (the moire limit).** `C_z(h) -> 0` as `abs h -> infinity`.

*Proof.* `C_z` is even, so take `h >= 1`. Fix `epsilon > 0` and suppose `abs C_z(h) >= epsilon`. Both factors are at most `1` in absolute value, so `abs C_t(h) >= epsilon` and `abs C3(h) >= epsilon`. Theorem 5.5 gives `r_2(h) <= log(K/epsilon)/log(1/rho)`, Theorem 5.2 gives `w_3(h) <= 2 log3(1/epsilon)`, and Theorem 5.10 leaves finitely many such `h`. □

At `epsilon = 10^(-3)` the proof confines the exceptions to `r_2 <= 20` and `w_3 <= 12`, a finite set of unknown size.

**Corollary 5.12 (non-Rajchman by design, Rajchman by convolution).** `sigma_z = sigma_t * sigma_3` is a continuous Rajchman measure, and it is the convolution of two continuous measures that are not Rajchman, invariant under `u -> 2u` and `u -> 3u` respectively.

*Proof.* Theorems 5.11 and 5.8, and Theorem 3.4 for the factors. □

A convolution of two measures that are not Rajchman can be Rajchman, and even Lebesgue measure, for an elementary reason. Let `mu` and `nu` be the Riesz products `prod_(k >= 0) (1 + cos 2 pi 100^k u)` and `prod_(k >= 0) (1 + cos 2 pi 10 . 100^k u)`, the weak-star limits of their partial products. A sum `sum_j delta_j 10^j` with every `delta_j` in `{-2, ..., 2}` vanishes only when every `delta_j` does, since its lowest nonzero digit is not divisible by `10`. So expanding the product gives `hat mu(h) = 2^(-w)` when `h = sum_k epsilon_k 100^k` with `epsilon_k` in `{-1, 0, 1}` and `w` of them nonzero, and `hat mu(h) = 0` otherwise, and `nu` is the same with `10 . 100^k`. Neither is Rajchman, since `hat mu(100^k) = hat nu(10 . 100^k) = 1/2`. By the same fact the two spectra, on the even and on the odd powers of `10`, meet only at `0`, so `hat mu(h) hat nu(h) = 0` at every `h != 0` and `mu * nu` is Lebesgue measure. So the phenomenon is not new in itself. What Corollary 5.12 adds is about these two measures: each is invariant under one base, their coefficients are nonzero together at many lags (Corollary 4.3, Fact 5.17), and the decay of the product is forced by the S-unit theorem, not by disjoint spectra.

The proof uses one property of each factor, decay in signed length at its own base, and so gives more.

**Proposition 5.13 (a convolution criterion).** Let `a, b >= 2` be multiplicatively independent, and let `mu`, `nu` be finite measures on `T` with `abs hat mu(h) <= f(l_a(h))` and `abs hat nu(h) <= g(l_b(h))` at every `h != 0`, for functions `f`, `g` tending to `0` at infinity. Then `mu * nu` is Rajchman.

*Proof.* If `abs(hat mu(h) hat nu(h)) >= epsilon`, then `f(l_a(h)) >= epsilon/norm(nu)` and `g(l_b(h)) >= epsilon/norm(mu)`, so both signed lengths are bounded, and Theorem 5.10 leaves finitely many `h`; `l_a(-h) = l_a(h)` covers negative lags. □

Neither factor needs to be Rajchman, since `l_a(a^k) = 1` leaves `hat mu(a^k)` free. Invariance under the two multiplications is not enough on its own: `delta_0` is invariant under every `u -> ku`, and `delta_0 * delta_0 = delta_0`. With decay in signed lengths the statement holds, and Thue-Morse needed the run lemma precisely because its coefficients do not decay in digit sums (Proposition 5.6).

### The rate

For `m >= 1` with base-`b` digits `d_(L-1) ... d_0`, `d_(L-1) != 0`, let `BC(m, b)` be the number of `j < L - 1` with `d_j != d_(j+1)`; so `BC(h, 2) = r_2(h) - 1`.

**Theorem 5.14 (Bugeaud, Cipu and Mignotte, Theorem 3.1).** Let `a` and `b` be multiplicatively independent integers. Then `BC(m, a) + BC(m, b) >= log log m/(log log log m + kappa) - 1` for every `m > 25`, where `kappa` is a positive number effectively computable in terms of `a` and `b` only.

It is read at source in their survey, which writes `m` in each base as a signed combination of the powers at its digit changes and follows Stewart's proof for bounded digit sums, which runs on Baker's linear forms in logarithms; it gives no value of `kappa`. The survey defines `BC(m, b)` as the number of times a digit different from the previous one is read, which can be read as the digit changes above or as the blocks, one more; Theorem 5.16 holds under either reading.

**Lemma 5.15 (digit changes of a signed sum).** `BC(h, b) <= 2 l_b(h) - 1` at every `h >= 1` and every base `b >= 2`.

*Proof.* Take a representation of `h` by `l = l_b(h)` terms `+-b^a`. It never uses one power with both signs, and never `b` copies of one signed power, since `b . b^a = b^(a+1)` would shorten it; so it merges into `h = sum_(i <= k) c_i b^(a_i)` with `a_1 < ... < a_k`, `1 <= abs c_i <= b - 1` and `k <= l`. Put `e_j = c_i` at `j = a_i` and `e_j = 0` elsewhere, and read digits from the bottom with a borrow: `g_0 = 0`, `d_j = (e_j + g_j) mod b`, `g_(j+1) = (e_j + g_j - d_j)/b`.

- If `g_j` is `0` or `-1`, then `e_j + g_j` lies in `[-b, b - 1]` and `g_(j+1)` is `0` or `-1`; and `sum_(j < J) d_j b^j + g_J b^J = sum_(j < J) e_j b^j` at every `J`.
- Past `a_k` the borrow is constant, and it is `0`, since `g_J = -1` with `J > a_k` would give `h = sum_(j < J) d_j b^j - b^J < 0`. So the `d_j` are the standard digits of `h`, and its top digit sits at `L - 1 <= a_k`.
- Where `e_j = 0` the borrow passes unchanged and `d_j` is `0` or `b - 1` as the borrow is `0` or `-1`; so `d_j = d_(j+1)` unless `j` or `j + 1` is some `a_i`.
- Each `a_i` lies in at most two pairs `(j, j + 1)` with `j + 1 <= L - 1`, and `a_k` in at most one, since `a_k + 1 > L - 1`. Hence `BC(h, b) <= 2k - 1 <= 2l - 1`. □

The generator checks the lemma at every `h < 3^12` for `b = 2` and `b = 3` with `0` failures, `l_b` computed by a digit recursion and confirmed by a breadth-first search over `[-3^13, 3^13]`.

**Theorem 5.16 (an effective rate).** `abs C_z(h) <= K 3^(-G(h)/4) < 1.0093 . 3^(-G(h)/4)` at every `h > 25`, with `G(h) = log log h/(log log log h + kappa) - 1` and `kappa` the constant of Theorem 5.14 at the bases `2` and `3`. Where `log log log h >= kappa`, this reads `abs C_z(h) <= 1.3284 (log h)^(-0.1373/log log log h)`.

*Proof.* `r_2(h) = BC(h, 2) + 1`, and Lemma 5.15 at `b = 3` with Lemma 5.1 gives `2 w_3(h) >= 2 l_3(h) >= BC(h, 3) + 1`. So `r_2(h) + 2 w_3(h) >= BC(h, 2) + BC(h, 3) + 2`, the number of blocks of `h` in base `2` plus that in base `3`, and this is at least `G(h)` by Theorem 5.14, an inequality that holds even if `BC` counts blocks, one more than changes. Theorems 5.5 and 5.2 give `abs C_z(h) <= K rho^(r_2(h)) 3^(-w_3(h)/2)`, and `rho^4 = 10609/40000 < 1/3` gives `rho < 3^(-1/4)`, so `abs C_z(h) <= K 3^(-(r_2(h) + 2 w_3(h))/4) <= K 3^(-G(h)/4)`. Where `log log log h >= kappa`, `G(h) >= log log h/(2 log log log h) - 1`, and `K 3^(1/4) = 1.328305...`, `log(3)/8 = 0.137327...` give the second form. □

The bound decays more slowly than every power of `log h`, and no lag is known past which `abs C_z <= epsilon` for any `epsilon`, since `kappa` is not computed. Bounding each layer alone and splitting `r_2 + 2 w_3` in halves would give only `3^(-G(h)/8)`.

**Fact 5.17 (the window maxima).** On every window `[2^k, 2^(k+1))`, `k <= 24`, the maximum of `abs C_z` is exact, from the recursions; from `k = 12` on:

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

Generator: `lab/rs/two-bases`, verb `windows`. Every maximiser has few runs and few balanced digits, such as `13122 = 2 . 3^8`, whose ternary layer is `C3(2 . 3^8) = -1/3` exactly, and `1048320 = 2^20 - 2^8`; the decay across windows is slow and uneven, as a limit carried by sparse integers should be.

**Fact 5.18 (the sparse census).** Below `10^40` exactly `15196` integers have `r_2(h) <= 6` and `w_3(h) <= 6`; the largest is `283462598656`, none lies in `[2^39, 10^40)`, and the largest `abs C_z` among those at least `2^32` is `1.584e-5`, at `34896543768`. Below `10^20` exactly `528556` integers have `r_2(h) <= 8` and `w_3(h) <= 8`; the largest is `4035260175885545472`. Generator: `lab/rs/two-bases`, verbs `sparse 6 6 40` and `sparse 8 8 20`.

The census walks every balanced-ternary word with at most `W` nonzero digits below the cut, `14536738777` and `20834959549` words, each integer once since balanced ternary is unique; it keeps those whose low `64` bits have at most `R` digit changes, a necessary condition, and checks the run count, the weight and the cut exactly in `192`-bit arithmetic. Nothing lies within `0.02` of either cut. Theorem 5.10 says each list is finite without a cut, and Theorem 5.14 bounds where it ends only through the uncomputed `kappa`, so the census is evidence about where the list ends and not a proof of it.

## What the second base keeps

The second base sends every shift correlation to `0` (Theorem 5.11) and the dilation means at `B = 2, 4, 8, 16, 3, 9, 27` to `0` (Theorem 4.5); it does not make the stream random. A one-base kernel test still detects `z`, with significance growing as `N^0.194` (Theorem 4.6), and two certificates of determinism survive any finite number of bases.

**Theorem 6.1 (the dilation cube).** Let `x = x_1 ... x_s` be a termwise product of sign sequences, `x_i` being `b_i`-automatic, the bases arbitrary. There are `P_i = b_i^(d_i)`, `d_i >= 1`, and a modulus `M = prod_i b_i^(e_i)` with

```
prod over S of x(M n prod_(i in S) P_i) = 1   at every n ,
```

the product over the `2^s` subsets `S` of `{1, ..., s}`.

*Proof.* Theorem 3.1 gives `e_i`, `d_i` with `x_i(b_i^(e_i) u) = x_i(b_i^(e_i + d_i) u)` at every `u`. Fix `i` and pair each subset `S` without `i` with `S + {i}`. The point `M n prod_S P_j` is `b_i^(e_i) u` for an integer `u`, and its partner is `b_i^(e_i + d_i) u`, so `x_i` takes the same value at both, and the `x_i` part of the cube is a product of squares. □

When the products `prod_S P_i` are distinct, as for pairwise coprime bases, the `2^s` points are distinct at `n >= 1`, and the cube has mean `0` on an independent uniform stream; one base needs two points, `s` bases need `2^s`. For `z` and `y` every `e_i = 0` and `d_i = 1`, since `t(2n) = t(n)`, `r(2n) = r(n)` and `t3(3n) = t3(n)`, so `z(n) z(2n) z(3n) z(6n) = y(n) y(2n) y(3n) y(6n) = 1` at every `n`. The battery at `N = 2^15` prints the square mean `+1.0000` for both, `-0.0030` for the centre column of Rule 30 and `+0.0010` for the `splitmix64` control. Generator: `lab/rs/design-ladder`, verb `battery`.

**Theorem 6.2 (zero entropy).** The factor complexities `p_z(n)` and `p_y(n)`, the numbers of distinct words of length `n` in `z` and `y`, are `O(n^2)`, so the topological entropy `lim (1/n) log p(n)` of `z` and `y` is `0`.

*Proof.* Let `b^(k-1) < n <= b^k`. A factor of length `n` of a signed design `x` sits inside two consecutive level-`k` blocks, and by Lemma 2.5 the block of index `j` is `x(j) B_k`, `B_k = (x(m))_(m < b^k)`; so the factor is fixed by its offset below `b^k` and the pair of signs `(x(j), x(j + 1))`, and `p_x(n) <= 4 b^k < 4bn`. Rudin-Shapiro's level-`k` blocks are `+-P_k` and `+-Q_k` by the block identity in the proof of Theorem 3.5, at most `16` ordered pairs, and `2^k < 2n`, so `p_r(n) <= 32n`. A factor of a termwise product is the product of the two factors at the same position, so `p_(xw)(n) <= p_x(n) p_w(n)`. □

The census is exact where every pair of layer factors occurs jointly, which the generator checks in a prefix of length `2^22`; the layer counts come from a prefix `18 b^k`, which holds every offset and every sign pair, and the Thue-Morse column is A005942.

| `n` | `p_t(n)` | `p_t3(n)` | `p_z(n)` | `p_z(n)/n^2` | `p_y(n)` |
|---|---|---|---|---|---|
| 8 | 22 | 28 | 176 | 2.75 | 232 |
| 16 | 46 | 60 | 1244 | 4.86 | 3100 |
| 32 | 94 | 124 | 5564 | 5.43 | 14368 |
| 64 | 190 | 252 | 23420 | 5.72 | 61472 |

`p_z(64)` is `0.49` of `p_t(64) p_t3(64)`, as the global sign pair `(u, w) -> (-u, -w)` predicts for two layers closed under complement. Zero entropy is what a compressor sees in the limit: the LZ78 rate at `2^18` bits is `0.840` bits per bit for `z` and `0.989` for `y`, against `1.215` for the `splitmix64` control, the phrase overhead of an incompressible stream. Generator: `lab/rs/design-ladder`, verbs `moire` and `battery`.

## Rounds: leaving the automatic class

**Definition 7.1 (rounds).** A *local rule* of radius `w` over a field `F` is a map `f : F^(2w+1) -> F`; it acts on rows `a(t, .) : Z -> F` by `a(t + 1, n) = f(a(t, n - w), ..., a(t, n + w))`. The *seed* is `a(0, n) = s(n)` for `n >= 0` and `a(0, n) = 0` for `n < 0`, and a generator reads a *row* `n -> a(T, n)`, `n >= 0`, or a *column* `t -> a(t, n_0)`. An *elementary* rule is a rule of radius `1` over `F_2`, numbered by Wolfram's convention, `sum f(l, c, r) 2^(4l + 2c + r)`. The *reach* `W_t` of an orbit at round `t` is the set of windows `(l, c, r)` that occur in row `t`, the zero exterior included. A rule `g = alpha l + beta c + delta r + epsilon` over `F_2` is *affine*; sixteen of the `256` elementary rules are. A rule `f` is *linear on its reach from round `t_0` with shadow `g`* when `g` is affine and `f = g` on `W_t` for every `t >= t_0`. A rule *exercises a nonlinearity at unbounded depth* when for every affine `g` infinitely many rounds `t` have a window in `W_t` on which `f` and `g` differ.

**Theorem 7.2 (bounded rounds).** If the seed is `b`-automatic, every row `a(T, .)`, read on `n >= 0`, is `b`-automatic, for every local rule of every radius.

*Proof.* After `T` rounds `a(T, n) = G(s(n - Tw), ..., s(n + Tw))` for one function `G`, with `s` read as `0` on negative arguments. By Lemma 2.4(2) each shift `n -> s(n + j)` is `b`-automatic, and by Lemma 2.4(1) so is `G` of them. □

So a fixed number of rounds, linear or not, stays on the first rung and dies under its kernel test (Theorem 3.1), at a depth that grows with `T`. The kernel count of a stream, the number of distinct words `(s(2^k n + r))_(n < 24)` over `r < 2^k` accumulated over `k = 0, ..., 10`, stops growing at the kernel size of an automatic stream and reads `1, 3, 7, ..., 2047` on a stream with no structure. On the Thue-Morse seed, rows `1, 2, 3, 4` of rule `30` stop at `8, 11, 14, 21` words and row `4` of rule `110` at `11`, the seed itself at `2`, while the `splitmix64` control reads `2047`. Generator: `lab/rs/design-ladder`, verb `rounds`.

**Theorem 7.3 (linear rounds at a matching base).** Let `p` be prime and `f` a rule linear over `F_p`, `a(t + 1, n) = sum_(abs j <= w) phi_j a(t, n + j)`, and let the seed be `p`-automatic. Then every column `t -> a(t, n_0)` is `p`-automatic, at every number of rounds.

*Proof.* Write the row `t` as the Laurent series `R_t = sum_n a(t, n) x^n`, so that `R_t = Phi^t u` with `Phi = sum_j phi_j x^(-j)` and `u = sum_(n >= 0) s(n) x^n`, and for `0 <= eta < p` let `Lambda_eta(Y) = sum_m Y(pm + eta) x^m` be the section. Over `F_p`, `Phi(x)^p = Phi(x^p)`, and `Lambda_eta(Phi(x^p) Y) = Phi(x) Lambda_eta(Y)`. So for `0 <= e < p`,

```
a(pt + e, pm + eta) = [x^(pm + eta)] Phi(x^p)^t Phi^e u = [x^m] Phi^t Lambda_eta(Phi^e u) :
```

the cell `(pt + e, pm + eta)` of the diagram of `u` is the cell `(t, m)` of the diagram of the seed `Lambda_eta(Phi^e u)`. Let `G` be the `p`-kernel of `s`, each element extended by `0` to negative arguments, and `V` the `F_p`-span of the shifts `x^i g`, `g` in `G`, `abs i <= w + 2`; `V` is finite and contains `u`. It is closed under `v -> Lambda_eta(Phi^e v)`: multiplying by `Phi^e` moves a shift by at most `(p - 1) w`, and `Lambda_eta(x^i g) = x^(-q) Lambda_(eta')(g)` with `eta - i = pq + eta'`, `0 <= eta' < p`, where `Lambda_(eta') g` lies in `G` and `abs q <= (abs i + p - 1)/p`, at most `w + 1` when `abs i <= pw + 2`. Iterating, `a(p^k t + tau, n_0) = a_(u_k)(t, floor(n_0/p^k))` for `0 <= tau < p^k` and some seed `u_k` in `V`, `a_(u_k)` its diagram. Once `p^k > abs n_0` the column index is `0` or `-1`, so the kernel elements of the column at those depths are among the columns `0` and `-1` of the diagrams of the finitely many seeds of `V`. □

Rule `90` is `Phi = x + x^(-1)` over `F_2` and rule `150` is `x + 1 + x^(-1)`, as in Martin, Odlyzko and Wolfram 1984, cited unread; from one cell rule `90` draws Pascal's triangle modulo `2`, which the generator checks on every cell `t < 4096`. Allouche, von Haeseler, Peitgen and Skordev 1996, Theorem 2, read at source, prove the double sequence of a linear automaton over `Z/p^l Z` from a finite seed `p`-automatic; their Theorem 3 rests on the same identity `r(X)^m = r(X^m)`, and the proof above runs it for automatic seeds. The centre columns of rules `90` and `150` on the Thue-Morse seed stop at `5` and `4` kernel words at base `2`, and rule `90` from one cell at `2`, while at base `3` the same columns keep growing, `38` and `1078` words at depth `6`, as Cobham's theorem requires of a stream that is not periodic. Generator: `lab/rs/design-ladder`, verb `rounds`.

**Theorem 7.4 (linear on the reach).** Let the seed be `2`-automatic and let the elementary rule `f` be linear on its reach from round `t_0` with shadow `g`. Then every column of the orbit of `f` is `2`-automatic, and by Theorem 7.2 so is every row.

*Proof.* From round `t_0` every window of row `t` lies in `W_t`, where `f = g`, so row `t + 1` is `g` of row `t`: the orbit of `f` from round `t_0` is the orbit of `g` from row `t_0`. Write `g = L + epsilon` with `L` linear and `a(t, n) = c_t + d(t, n)`, where `c_t` is the value of row `t` far to the left; then `c_(t+1) = (alpha + beta + delta) c_t + epsilon`, an ultimately periodic sequence, and `d(t + 1, .) = L d(t, .)`. Row `t_0` of `d` vanishes at `n < -t_0`. Shifted by `t_0`, `m -> d(t_0, m - t_0)` on `m >= 0` is `2`-automatic: the formula `a(t_0, n) = G(s(n - t_0), ..., s(n + t_0))` in the proof of Theorem 7.2 holds at every integer `n`, so Lemma 2.4(2) and (1) apply to it and to its sum with the constant `c_(t_0)`. Theorem 7.3, applied to this seed and the linear rule `L`, makes every column of its orbit `2`-automatic, and from round `t_0` on a column of `d` is one of them, delayed by `t_0` rounds. So a column of `f` is, after finitely many first terms, a column of that orbit shifted in time by `t_0`, plus `c_t`; Lemma 2.4(2) for the time shift and Lemma 2.4(1) and (3) close. □

**Corollary 7.5 (the necessary condition).** Let the seed be `b`-automatic with `b` a power of `2`, and let `f` be an elementary rule.

1. Every row of the orbit is `b`-automatic (Theorem 7.2).
2. A column that is not `b`-automatic needs a nonlinearity exercised at unbounded depth (Theorem 7.4).

So a design-fed elementary stream leaves its seed's automatic class only in a column, and only through a base that is not a power of `2` or a nonlinearity the orbit exercises at unbounded depth.

*Proof.* (1) is Theorem 7.2, which holds at every base. (2) By Lemma 2.4(4), `b`-automatic and `2`-automatic agree, so the seed is `2`-automatic. If for some affine `g` only finitely many rounds meet a window where `f != g`, then `f` is linear on its reach from the round after the last of them, and Theorem 7.4 makes every column `2`-automatic, hence `b`-automatic. □

**Theorem 7.6 (the shadow criterion).** An elementary rule `f` is linear on its reach from round `0` with shadow `g` exactly when `f = g` on `M(g, s)`, the set of windows of the `g`-orbit of the seed `s` over all rounds and cells. For the one-cell seed and the Thue-Morse seed, `M(g, s)` is computed exactly, as the window set of a finite closure, for each of the `16` affine `g`.

*Proof.* If `f = g` on `M(g, s)`, the two orbits agree by induction on `t`, so the reach of `f` lies in `M(g, s)`, where `f = g`. Conversely, if `f = g` on the reach of `f` at every round, the orbits agree and `M(g, s)` is that reach. □

*The computation (by exact computation).* Write `g = L + epsilon` with `L` linear of Laurent polynomial `Phi` and `sigma = alpha + beta + delta`. As in the proof of Theorem 7.4, the `g`-diagram of `s` is its `L`-diagram plus the constant `c_T` of row `T`, which is `0` when `epsilon = 0`, `T mod 2` when `epsilon = sigma = 1`, and `[T != 0]` when `epsilon = 1` and `sigma = 0`. By the proof of Theorem 7.3 at `p = 2`, the cell `(2T + e, 2M + eta)` of the `L`-diagram of a seed `v` is the cell `(T, M)` of the `L`-diagram of `Lambda_eta(Phi^e v)`.

- A *state* is a triple `((v_j, o_j))_(j = -1, 0, 1)` and a tag bit `q`. Each `v_j` lies in the `F_2`-span `V` of the fifteen functions `n -> chi(n - i)`, `abs i <= 2`, with `chi` the indicator of `n = 0`, the Thue-Morse bit `s_2(n) mod 2` on `n >= 0` extended by `0`, or the indicator of `n >= 0`; each offset `o_j` lies in `{-1, 0, 1}`. The state stands for the window of cells `(T, M + o_j)` of the `L`-diagrams of the `v_j`, plus `q`.
- The start is `((s, -1), (s, 0), (s, 1))` with `q = 0`, the window at `(T, M)` of the diagram of `s`.
- Reading the low binary digits `e` of `T` and `eta` of `M` sends each `(v, o)` to `(Lambda_(eta'')(Phi^e v), o')`, where `eta + o = 2 o' + eta''` and `eta''` is `0` or `1`. The tag becomes `e` at the first digit when `epsilon = 1`, and then stays `q` when `sigma = 1` and becomes `q or e` when `sigma = 0`; it stays `0` when `epsilon = 0`. Since `t(2m) = t(m)` and `t(2m + 1) = 1 + t(m)` on `m >= 0` over `F_2`, a section of a shift of one of the three functions is a combination of shifts of them, and the generator asserts at every step, after `Phi^e` and after the section, that no shift leaves `abs i <= 2`.
- Halving `T` and `M` ends at `T = 0` and `M` in `{0, -1}`, where a state's window is `(v_j(M + o_j) + q)_j`.

The generator lists every state reachable from the start by breadth-first search, and `M(g, s)` is the union of the windows of the reachable states at `M = 0` and `M = -1`, the start included. □

The generator prints the exact set `M(g, s)` beside the set seen by simulation to `t <= 4096` for the `16` affine orbits of both seeds: `32` pairs, `0` differences. Generator: `lab/rs/design-ladder`, verb `shadow`.

**Theorem 7.7 (the rules linear on their reach).** Of the `256` elementary rules, exactly `90` are linear on their reach from round `0` from one cell, and exactly `20` from the Thue-Morse seed.

*Proof (by exact computation, with Theorem 7.6).* From one cell the `90` are the `16` affine rules and `74` more, grouped by shadow: rule `90` shares its shadow with `18, 26, 82, 146, 154, 210` and `218`, eight rules that draw the gasket; rules `0, 170, 204` and `240` each have fifteen non-affine twins, whose cell dies, shifts or stays; rule `255` has seven, whose row is all ones from round `1`; the other ten affine rules stand alone. From the Thue-Morse seed the `20` are the `16` affine rules and `42, 76, 112` and `128`, which differ from `170, 204, 240` and `0` only on the window `111`, which the orbits of those four rules, shifts and the zero row of the cube-free Thue-Morse word, never make. □

Rule `18` can be checked by hand: `18 = (l + r)(1 + c)` equals rule `90` on every window without two adjacent ones, and from one cell the live cells of rule `90` at round `t` sit at positions of the parity of `t`. Simulation to `t <= 4096` finds the same two sets, and no rule meets a new window after round `3` from one cell or round `1` from Thue-Morse. The generator prints `0` cells differing from rule `90` on `t < 4096` for rules `18` and `146`, against `88574` and `265720` for rules `22` and `126`, all counted on `n >= 0`. On the Cantor seed rule `18` is rule `90`, since every point of the seed has digits `0` and `2` and so sits at an even position.

**Fact 7.8 (the linear-on-reach columns).** Among the `90` and the `20` rules of Theorem 7.7, the largest depth-`10` kernel count of a column `n_0 = -1, 0, 1`, read to `t < 2^15`, is `5` from one cell and `13` from the Thue-Morse seed, against `2047` for the control, and none of the `110` pairs grows. Generator: `lab/rs/design-ladder`, verb `escape`.

The necessary condition is not sufficient, and the two refutations below are the reason Corollary 7.5 is stated as one implication.

**Remark 7.9 (nonlinearity and unbounded rounds do not force escape; refuted).** The statement that a design-fed automaton stream escapes the automatic class whenever its rounds are unbounded and its rule is nonlinear is false. Rules `18` and `146`, of algebraic degree `2` and `3`, draw Pascal's triangle modulo `2` from one cell at every round (Theorem 7.7), so their columns are `2`-automatic by Theorem 7.3; the generator finds `0` differing cells on `t < 4096`. That the two rules draw the gasket from one cell is a known picture of the single-site catalogues (Section 8); the remark records only what it refutes.

**Remark 7.10 (exercise at unbounded depth does not force escape; refuted).** The statement that exercising a nonlinearity at unbounded depth forces escape is false, by two rules that exercise in two different places. Write `P(s, n) = [s + n even, binom(s, (s + n)/2) odd]`, the diagram of rule `90` from one cell.

- *Rule `182`* draws, from one cell, `a(t, n) = [abs(n) <= t] (1 + P(t + 1, n))`, the complemented gasket inside the cone, a known picture of the single-site catalogues (Section 8). By induction on `t`: inside the cone the window is the complement of a window of a row of `P`, which lies in `{000, 001, 010, 100, 101}`, and rule `182` agrees on the complements `{111, 110, 101, 011, 010}` with rule `165 = 1 + l + r`, the Pascal recursion complemented. At the edge `f(100) = 1` makes the new cell, and the cell beside it reads `([t odd], 1, 0)`, where `f(110) = 0` and `f(010) = 1` give `[t even]`, the value of `1 + P(t + 2, t)`. Its columns are `2`-automatic by Theorem 7.3 and Lemma 2.4. The nonlinearity sits on the frontier of the light cone only: `4084` of the `4097` rounds to `4096` meet a window no affine rule fits, and the interior windows never do.
- *Rule `22`*, exactly one live cell in the window, draws for `t >= 1` `a(t, n) = P(t, n)` when `t + n` is even and `[P(t - 1, n) = 1, P(t - 1, n - 2) = P(t - 1, n + 2) = 0]` when it is odd, a block code of one row of Pascal's triangle. By induction on `t`, with a computation: the step from row `t` to `t + 1` concerns seven consecutive cells of one row of `P`; the closure of Theorem 7.6, run on windows of seven cells with a parity tag, lists the `22` windows that occur, the same `22` the rows to `4096` show, and on each the rule `22` step of the block code is the block code of the next row, `0` violations; the base `a(1, .) = 111` holds. Its columns are `2`-automatic by Theorem 7.3 and Lemma 2.4. The exercise is in the bulk: row `2^k - 1`, `k >= 2`, is the blocks `111` around the isolated live cells of row `2^k - 2`, with the windows `000, 001, 011, 111, 110, 100`, which no affine `g` fits, since `g(000) = 0`, `g(001) = g(100) = 1` and `g(011) = 0` force `g(111) = 1` while `f(111) = 0`.

The generator prints `0` differing cells for `t <= 4096` for both closed forms. Generator: `lab/rs/design-ladder`, verb `lemma`.

**Proposition 7.11 (rule 126 is rule 18 on differences).** For rule `126`, not all three equal, the difference field `d(t, n) = a(t, n) + a(t, n + 1)` obeys rule `18`: `d(t + 1, n) = (1 + d(t, n))(d(t, n - 1) + d(t, n + 1))`. From one cell `d` is the rule-`18` orbit of the pair `11` at `-1, 0`.

*Proof.* The new cell is `1` exactly when its window is not constant, that is `a(t + 1, n) = d(t, n - 1) or d(t, n)`. So `d(t + 1, n) = (d(t, n - 1) or d(t, n)) + (d(t, n) or d(t, n + 1))`, which is `0` when `d(t, n) = 1` and `d(t, n - 1) + d(t, n + 1)` when `d(t, n) = 0`. One live cell at `0` gives `d(0, -1) = d(0, 0) = 1` and `d = 0` elsewhere. □

The pair `11` is an input with two adjacent ones, so rule `18` does not reduce to rule `90` here; the generator prints `0` differing cells for `t <= 4096`, and the interior windows of rule `126` are fit by rule `90` at every round `t <= 4096` except the `11` rounds `2^k - 2`, `k = 2, ..., 12`.

**Fact 7.12 (the one-cell census).** From one cell, of the `166` elementary rules not linear on their reach, `56` exercise for the last time at round `2` or earlier, and `110` exercise to round `4094` or later, none lying between; the last exercise of a rule `f` is the least, over the sixteen affine `g`, of the last round `t <= 4096` whose reach `W_t` holds a window on which `f` and `g` differ. Reading three columns `n_0 = -1, 0, 1` to `t < 2^15`, calling a column periodic when it has a period of at most `512` after a preperiod of at most an eighth of its length, and a plateau when the kernel count of its last depth equals that three depths earlier: of the `110`, `92` have all three columns ultimately periodic, `13` have a column that does not plateau (rules `30, 45, 73, 75, 86, 89, 101, 110, 124, 135, 137, 149` and `193`, with `403` to `2047` kernel words), and `5` have an aperiodic column and plateau at `4` to `5` words: rules `22, 126, 129, 161` and `182`. Generator: `lab/rs/design-ladder`, verb `escape`.

Rules `22` and `182` are proved automatic in Remark 7.10, each a different way for a nonlinearity to be exercised without effect; for `126`, `129` and `161` the plateau is evidence only (Section 9). An automatic stream with a large kernel looks the same at depth `10` as one with none, so every kernel count of this section is evidence about escape and never a proof of it.

## What is in print

**Normal numbers in two bases.** W. Schmidt 1960, read at source in its scan, proves that normality to base `r` implies normality to base `s` exactly when `r^n = s^m` for some positive integers (Theorem 1), and builds the counterexamples as images of Lebesgue measure on a Cantor set of base-`s` expansions with digits below `t`. His Lemma 2 bounds the product that controls the Fourier coefficient of that measure at `w` by `a_2^z`, `0 < a_2 < 1`, once `w` has at least `z` adjacent digit pairs satisfying a nondegeneracy condition, and his Fundamental Lemma averages these bounds along `r^n` through the orders of `r` modulo prime powers. Pollington 1981, read at source in its scan, states the averaged form as his Lemma 1, a bound `2N^(1 - a(r, s))` on a sum of products of cosines, credits it to Schmidt's Acta Arithmetica paper of 1961-62, unread here, and proves that the numbers normal to every base of one class and to no base of the other have Hausdorff dimension `1`. Cassels 1959 is cited unread. Brown, Moran and Pearce 1985 and 1986 are cited unread: the first, by its title, takes the method to Riesz products, and the second decomposes numbers into summands with prescribed normality properties. In this lineage the decay of the coefficients is counted over the digit pairs of the frequency in one base, and the second base enters through the orbit `r^n`. Here the two bases enter symmetrically, each measure decaying in signed length at its own base, and the integers sparse in both are excluded by the S-unit theorem. The Thue-Morse layer is a Riesz-type product, as the Introduction notes, but not a product of independent digit factors: its coefficient is the matrix product of Lemma 5.3, which is why Theorem 5.5 needs a joint quadratic form rather than a factor-by-factor bound, and why its natural count is runs, not digit pairs of one kind. No source read states Theorem 5.5, Lemma 5.4 or Lemma 5.15: the sources on Thue-Morse correlations read below, Besineau, Kamae, Coquet, Kamae and Mendes France, Coquet 1979, Queffelec 1979 and Zaks 2001, hold no bound in the number of runs, no signed-length form of the convolution theorem and no `BC <= 2l - 1`, while Mahler 1927, Kakutani 1972, Queffelec 1987 and Baake and Grimm 2013 remain unread, and the run lemma and its bracket may be there.

**Rajchman measures.** Lyons 1995, read at source in the author's copy, surveys Rajchman measures from Riemann and Lebesgue through the uniqueness problem: Rajchman and Milicer-Gruzewska proved that the Rajchman measures form a band, closed under absolute continuity; F. Riesz built Riesz products to exhibit a continuous measure that is not Rajchman; and the Cantor-Lebesgue measure has `hat mu(n) = hat mu(3n)`, the one-ratio invariance of Theorem 3.4. The survey says nothing on convolutions. That the Rajchman measures are an ideal of the convolution algebra is the inequality `abs(hat(mu * nu)) <= abs(hat mu) norm(nu)`; Corollary 5.12 is about the opposite direction, a convolution of two factors outside the ideal landing in it, and that direction is elementary in itself: the two Riesz products after Corollary 5.12, on spectra that meet only at `0`, convolve to Lebesgue measure. The content of Corollary 5.12 is the specific measures: each is invariant under one base, their coefficients are nonzero together at many lags, and the decay of their product is forced by the S-unit theorem. Kamae 1977, read at source, proves that the shift systems of `e(a s_p(n))` and `e(b s_q(n))` are disjoint and have mutually singular spectral measures for coprime `p`, `q` and irrational `a`, `b`; the signs here have `a = 1/2` and a ternary digit count, outside its hypotheses. Coquet, Kamae and Mendes France 1977, on the spectral measure of `q`-multiplicative sequences, a class that contains the signed designs, is read at source, as the one-base paragraph below says; the correlation recursion of Lemma 2.5 and the continuity of `sigma_t` and `sigma_3` in Theorem 5.8 are classical, the recursion theirs, and are reproved, not claimed. The continuity of the Thue-Morse measure goes back to Mahler 1927 and Kakutani 1972, and the Lebesgue spectrum of Rudin-Shapiro, Theorem 3.5, to Rudin 1959 and Shapiro 1951, all cited unread. Host 1995, cited unread, and Hochman and Shmerkin 2015, read at the abstract, which gives a fractal-geometric condition for a measure to be supported on points normal in a given base, are statements about typical points of a measure, not about the coefficients of a convolution, and nothing here uses them.

**One ratio.** That a measure invariant under one integer ratio is never Rajchman is the identity `hat sigma(bh) = hat sigma(h)`, used in Theorem 3.4 and proved for the digit measures of weighted designs in the [dilated stacks](../../dilations/notes/dilations.md) chapter of this tree. For self-similar measures the converse direction is the Rajchman problem. Bremont 2021, read in its arXiv version, Theorem 2.3: the invariant measure of a family of affine maps of the line with no common fixed point, contracting on average, is not Rajchman exactly when an affine change of variables puts the family in Pisot form, every ratio a power of one number `theta` with `1/theta` Pisot, and the coefficients of the transformed measure along `theta^(-k)` do not tend to `0`. Varju and Yu 2022, read at the abstract, give quantitative decay for self-similar measures and characterise the self-similar sets of uniqueness. Rapaport 2022, Corollary 1.6, read at source in this tree's dimensions topic, proves for affinely irreducible systems in `R^d` that a self-similar measure that is not Rajchman has every ratio an integer power of one algebraic integer. In all three the ratios act inside one self-similar measure. In Corollary 5.12 each measure keeps one ratio and the two meet only through a convolution, which is not self-similar; neither statement implies the other.

**Convolutions of Cantor measures.** Senge and Straus 1973, unread, is read through Nazarov, Peres and Shmerkin 2012, read at source in its arXiv version: answering a question of Salem on the Fourier transform of a measure carried by the sum of two Cantor sets, Senge and Straus characterised the pairs `a, b` for which the transform of `mu_a * mu_b`, the convolution of the two central Cantor measures, does not tend to `0`, namely only when `1/a` and `1/b` are Pisot and `log b/log a` is rational; at `a = 1/3` and `b = 1/4` neither factor's transform tends to `0` and the convolution's does. Besineau 1972, page 401, read as images, records Salem's question and Straus's conjecture, the finiteness of the integers with bounded digit sums in two coprime bases, as its answer in part, and that finiteness is the Senge and Straus theorem that Theorem 5.10 contains. So the phenomenon of Corollary 5.12, two measures outside the Rajchman class whose convolution lies in it through the arithmetic of two ratios, is theirs; what this paper adds is its form for the spectral measures of two digital sequences, whose coefficients decay in signed length and not in digit sums. Nazarov, Peres and Shmerkin's own Theorem 1.1, the correlation dimension of `mu_a * (mu_b o S_lambda)` for `log b/log a` irrational, and Theorem 4.1, singular convolutions on a dense set of `lambda` for Pisot pairs with irrational `log b/log a`, are statements about Cantor measures and say nothing of digital sequences.

**Products across bases.** Besineau 1971-72, the seminar version of Besineau 1972, read at source in its Numdam scan: its Lemma 4 is the product rule of Theorem 4.2 for `f_i(n) = e(alpha_i s_(g_i)(n))` at coprime bases, by the same mechanism, progressions of ratio a power of `g_i` on which the digit-sum difference is constant, with ratios coprime across the bases; Theorem 4.2 runs it for any sign tile and for Rudin-Shapiro. His Theorem 2 gives the product over pairwise coprime bases a correlation of null mean square as soon as one `alpha_i (g_i - 1)` is not an integer, which is the continuity of the convolution in Theorem 5.8 for digit sums; nothing there says that the correlation of the product tends to `0`, the content of Theorem 5.11. Kamae 1977, as above, proves the mutual singularity of the two layer measures for irrational weights and leaves the weight `1/2` of the signs outside.

**One base.** Coquet, Kamae and Mendes France 1977 is read at source in its Numdam scan, whose text layer reads except for the displayed formulas of Theorem 3: that theorem gives a `q`-multiplicative sequence of modulus `1` a correlation and the recursion `gamma_m(q_(m+1) n + a) = A_m(a) gamma_(m+1)(n) + B_m(a) gamma_(m+1)(n + 1)`, of which Lemma 2.5 is the case of a signed design at one base; Theorem 4 gives the mass of an atom of its spectral measure as an infinite product, and Theorems 5 and 6 make the measure atomic when it is not diffuse and singular when it is, for bounded `q`; the continuity of `sigma_t` and `sigma_3` in Theorem 5.8 is their diffuse case, reproved there through Wiener's identity. Queffelec 1979, read at source in pages 385-402, writes that measure as a generalised Riesz product (Proposition 1), proves it discrete, or with a convolution power in `L^1`, or with strongly independent powers and then singular (Theorem 2), and states on page 400, after Coquet, Kamae and Mendes France, that for bounded `q` it is singular to every Rajchman measure; a search of its whole text finds no Senge, no Straus, no Stewart and no convolution across bases. Coquet 1979, read at source in its Numdam scan, calls a sequence whose correlation tends to `0` pseudo-random in the sense of Bass, which is the Rajchman property of its spectral measure (Section I.2), and proves in Section II.5 that no `q`-multiplicative sequence of modulus `1` is, through `abs gamma(0) <= q abs gamma(q^(k + 1) - q^k) + (q - 1) abs gamma(q^(k + 1))` at every `k`; Theorem 3.4 is the signed-design case, with the lag named. In those words Theorem 5.11 says that `z` is pseudo-random in the sense of Bass while neither of its layers is. Coquet's Section II.7, that a product of a `q_1`-multiplicative and a `q_2`-multiplicative sequence of non-empty Fourier-Bohr spectrum has non-empty spectrum, does not bear on `z`, whose layer measures are continuous.

**Digit changes.** Blecksmith, Filaseta and Nicol 1993, read at source in the Acta scan, prove as Theorem 1 that for `log a/log b` irrational the number of blocks of equal digits of `a^n` in base `b` tends to infinity, by writing the number as the signed sum `sum (d_(j-1) - d_j) (b^(k_j) - 1)/(b - 1)` over its blocks and applying a Baker lower bound for linear forms in logarithms; the number of blocks is one more than `BC`. Bugeaud, Cipu and Mignotte 2013 quote that theorem in their Section 3 and extend it in Theorem 3.1 from the powers `a^n` to every integer, with Stewart 1980, unread, as the digit-sum model; the same signed sum over blocks, read in the other direction, is Lemma 5.15. At `h = 2^k` the limit of Theorem 5.11 is a corollary of the 1993 theorem alone: `C_z(2^k) = -C3(2^k)/3`, Lemma 5.15 at `b = 3` with Lemma 5.1 gives `2 w_3(2^k) >= BC(2^k, 3) + 1`, and `BC(2^k, 3)` tends to infinity. Neither the Senge and Straus nor the Stewart finiteness reaches the lags `2^k - 1` of Proposition 5.6, whose binary digit sum is unbounded; blocks do, as Section 5 needs.

**The sharp Cesaro exponents.** Zaks 2001, read at source in Section 2.2 of the open copy, derives from the recursion of Lemma 2.5 at the Thue-Morse tile that the sums of `C_t^2` over `(n, 2n]` obey `S_(2j) = (3/2) S_j + (1/2) Pi_j` and `Pi_(2j) = -S_j - Pi_j`, whose ratio map `kappa -> 1/kappa + 1/2` has the attracting fixed point `(1 + sqrt 17)/4`, so that the correlation dimension of the Thue-Morse spectral measure is `D_2 = 3 - log2(1 + sqrt 17) = 0.64298`, a value it credits to Zaks, Pikovsky and Kurths 1997, unread here; the same section states, without derivation, `3 - log3 11 = 0.8173` for the 3-tupling sequence, whose correlation obeys the recursion of `C3`. Against Theorem 5.8 these are the sharp exponents of the two Cesaro means, where the proof gives `0.4006` and `0.5350` and the measured slopes are `0.6324` and `0.7882` per bit: the theorem's exponents are not sharp, and the sharp ones are in print from the same recursion. The two printed exponents sum to more than `1`, so the decorrelation heuristic of Conjecture 9.1, which with the proved exponents bounds the dyadic block sums of `C_z^2` only by a growing multiple, would with the sharp ones make them fall geometrically, as the measured partial sums do; a proof still needs the binary and ternary digits of `h` to decorrelate on average, which no source read gives.

**Automata.** Cobham's theorem, Cobham 1969, is used as stated through Bes and through Durand and Rigo, read at source in this tree's [Cobham](../../cobham/notes/cobham.md) topic; the kernel criterion is cited from Allouche and Shallit 2003, unread. Rule `90` and rule `150` as Laurent polynomials over `F_2` are Martin, Odlyzko and Wolfram 1984, cited unread. Allouche, von Haeseler, Peitgen and Skordev 1996, read at source, prove the double sequence of a linear automaton over `Z/p^l Z` from a finite seed `p`-automatic and do not treat nonlinear rules; Theorem 7.3 runs the same Frobenius identity for automatic seeds. The single-site diagrams of the elementary rules are catalogued in Wolfram 1983 and Wolfram 2002, both cited unread, and the pictures of rules `18`, `146` and `182` from one cell, the gasket of rule `90` for the first two and its complement inside the cone for the third, are taken as known. What is not in any source read is the reach criterion of Theorems 7.4 and 7.6, the exact closure of Theorem 7.6 and the counts `90` and `20` of Theorem 7.7. Rule `30` as a generator is Wolfram 1986, cited unread; nothing here depends on it.

## Open problems

Four conjectures and four questions are left.

**Conjecture 9.1 (the square-summable spectrum).** `sum_h C_z(h)^2 < infinity`, so that `sigma_z` is absolutely continuous with a square-integrable density.

The partial sums over `1 <= h < 2^L` read `0.043031`, `0.043187` and `0.043240` at `L = 17, 21, 25`, with block increments falling by a factor near `3` every four bits (`lab/rs/two-bases`, verb `decay`). It would make `sigma_z` absolutely continuous. Absolute continuity of a convolution of two measures that are not Rajchman is not the point, since the Riesz pair after Corollary 5.12 convolves to Lebesgue measure; the interest is this measure, whose factors are nonzero together at many lags and whose decay is known only through the S-unit theorem. The layer bounds are not enough even heuristically: they bound `C_z^2` by `K^2 rho^(2 r_2) 3^(-w_3)`, and if those decorrelated as the measured means seem to, they would bound the block sums only by a multiple of `2^(0.065 L)`, which grows; the exponent `1 - theta - log3(9/5)` is `0.064291` to six places, rounded up to `0.065` (verb `decay`). A proof needs the binary and ternary digits of `h` to decorrelate on average.

**Question 9.2 (the constant of the rate).** Compute the constant `kappa` of Theorem 5.14 at the bases `2` and `3`. That needs Stewart's Baker argument run with explicit linear-form bounds; it would turn Theorem 5.16 into a number and certify the lists of Fact 5.18 complete.

**Conjecture 9.3 (escape).** The centre columns of rule `90` and rule `150` on the base-`3` Cantor seed, of rule `18` on the Thue-Morse seed, and of Rule `30` from one cell are automatic in no base.

Their kernel counts at base `2`, depth `10`, reach `803`, `874`, `1074` and `2047`, the control `2047`, and the first two read `891` and `1092` at base `3`, growing at every depth (`lab/rs/design-ladder`, verb `rounds`). The first two would show a design-fed linear round escaping through the base alone, the third a nonlinearity that bites because the seed's adjacent ones break the parity pattern of Theorem 7.7. For Rule `30` not even non-periodicity is proved here.

**Conjecture 9.4 (three more automatic rules).** Rule `126` from one cell, and rules `129` and `161` with it, have `2`-automatic columns.

They read `5`, `3` and `5` kernel words at the columns `-1, 0, 1` over `t < 2^15`, and the interior windows of rule `126` are fit by rule `90` at every round `t <= 4096` except the `11` rounds `2^k - 2` (Fact 7.12, Proposition 7.11). Every round but `2^k - 2` is then linear on its reach, with a `2`-automatic schedule of shadows; partial products of an automatic schedule of Laurent polynomials need not be automatic, so Theorem 7.4 does not apply as it stands.

**Conjecture 9.5 (the density law on one cell).** From one cell, an aperiodic column of an elementary rule leaves the kernel plateau exactly when the exercise has positive density, the density being the least, over the sixteen affine `g`, of the share of interior cone cells `(t, n)`, `2048 <= t <= 4096` and `abs n <= t - 2`, whose window is one on which the rule and `g` differ.

The five rules that stay read at most `0.0092` and the thirteen that leave read `0.0778` to `0.2491`, a gap of a factor `8` over `18` rules, `11` up to mirror images and equal diagrams, which is evidence and not a law (`lab/rs/design-ladder`, verb `escape`). It fails on the Thue-Morse seed: in its bulk `44` of `102` rules exercised to round `4094` or later plateau at `9` to `39` words at densities `0.047` to `0.219`. A rule with density between `0.0092` and `0.0778` would break it.

**Question 9.6 (the constant per run).** Is the constant per run of Theorem 5.5, which lies in `[0.711872, 0.717636]`, the joint spectral radius of `{A((-1/2)^b) : b >= 1}` attained by the run pattern `(1, 1, 2)`?

**Question 9.7 (no base at all).** Is `z` or `y` automatic in some base other than `2` and `3`? At a base `k` with `k^a = 2^c` or `k^a = 3^c`, Lemma 2.4(4) reduces the question to Theorem 4.1; at bases such as `5`, `6` or `10` nothing here decides it. A superlinear lower bound on the factor complexity would settle it, since an automatic stream has linear complexity, a classical fact not read here; the census of Theorem 6.2 is exact to `n = 64` and a lower bound for every `n` needs the joint pairs complete at every length.

**Question 9.8 (Rudin-Shapiro along `3n`).** Do the sums `sum_(n < N) r(n) r(3n + c)` grow faster than `sqrt N`? At `N = 2^21` and `c = 2` the sum is `-14682`, `10.1` times `sqrt N`, and the battery flags `y` at the dilation `(3, 2)` with `0.0277`, `5.0` standard deviations (`lab/rs/design-ladder`, verbs `moire` and `battery`). An exact recursion on the carry states of `3n + c`, the analogue of Lemma 4.4 for Rudin-Shapiro, would decide it.

## Reproducibility

Every number of this paper is printed by one of two single-threaded Rust studies, each run from the repository root as `cargo run --release --manifest-path research/random/lab/rs/<study>/Cargo.toml -- <verb>`.

- `design-ladder`, verbs `codes` (Proposition 3.2, Theorem 3.4 and its table, about 2 s), `rudin` (Remark 3.3, Theorem 3.5, under 1 s), `battery` (Fact 3.6, the battery values of Sections 4 and 6, about 2 s), `moire` (Corollary 4.3, Theorems 4.5, 4.6 and 6.2, Question 9.8, about 4 s), `rounds` (Theorems 7.2, 7.3 and 7.7, Conjecture 9.3, under 1 s), `shadow` (Theorems 7.6 and 7.7, about 5 s), `escape` (Facts 7.8 and 7.12, Conjecture 9.5, about 30 s) and `lemma` (Remark 7.10, Proposition 7.11, under 1 s).
- `two-bases`, verbs `lemma` (Lemma 5.4 in exact rationals with its Sturm counts, Proposition 5.6(1), under 0.1 s), `decay` (Theorems 5.2, 5.5 and 5.8 with their exponents, Proposition 5.6, Conjecture 9.1, about 3 s), `windows` (Fact 5.17, about 1.5 s), `changes` (Lemma 5.15 and the constants of Theorem 5.16, under 0.3 s), and `sparse 6 6 40` and `sparse 8 8 20` (Fact 5.18, about 10 s and 24 s).

Correlations are exact: `C_t(h) = x/(3 . 2^L)` for an integer `x` read off the binary digits of `h`, in `192`-bit integers past `2^64`, and `C3(h)` is `0` or `+-3^(-m)`. The figure `paper-design-generators` recomputes in floating point the six endpoint values of Lemma 5.4 and the discriminant identity, checks `A(-1/2) e_0 = (-1/3, -1/3)`, and computes for each of the twelve run operators drawn, and for `97` equally spaced `t` in `[-1/2, 1/4]`, the largest value of `Q(A(t) u)/Q(u)` as a generalized eigenvalue, refusing to draw if any exceeds `103/200` or if the run of length `2` is not the largest of the twelve; the exact check of Lemma 5.4 is the `two-bases` verb `lemma`. Every Theorem, Proposition, Lemma proved here and every Fact is a dated Proved or Verified row of this tree's claim pages, [the design ladder](../claims/design-ladder.md) and [the moire limit](../claims/two-bases.md); the Conjectures are rows of the same pages, and the two Remarks of Section 7 are its Refuted rows.

## References

- Allouche and Shallit 2003, Automatic Sequences: Theory, Applications, Generalizations, Cambridge University Press; cited for the kernel criterion, not read. [doi.org/10.1017/CBO9780511546563](https://doi.org/10.1017/CBO9780511546563)
- Allouche, von Haeseler, Peitgen and Skordev 1996, Linear cellular automata, finite automata and Pascal's triangle, Discrete Appl. Math. 66, 1-22; read at source in the authors' copy, Theorems 1 to 3. [webusers.imj-prg.fr](https://webusers.imj-prg.fr/~jean-paul.allouche/ahps96.pdf)
- Baake and Grimm 2013, Aperiodic Order, Vol. 1: A Mathematical Invitation, Cambridge University Press; cited unread, as an account of the Thue-Morse correlations that was not searched.
- Besineau 1972, Independance statistique d'ensembles lies a la fonction "somme des chiffres", Acta Arith. 20, 401-416; read at source in the Acta scan as images, page 401; its seminar version, Seminaire Delange-Pisot-Poitou 13 (1971-72), exp. 23, 1-8, is read at source in its Numdam scan, Lemmas 3 and 4 and Theorems 2 and 3. [doi.org/10.4064/aa-20-4-401-416](https://doi.org/10.4064/aa-20-4-401-416)
- Blecksmith, Filaseta and Nicol 1993, A result on the digits of a^n, Acta Arith. 64, 331-339; read at source in the Acta scan, Theorem 1. [doi.org/10.4064/aa-64-4-331-339](https://doi.org/10.4064/aa-64-4-331-339)
- Bremont 2021, Self-similar measures and the Rajchman property; read at source in its arXiv version, Theorem 2.3. [arxiv.org/abs/1910.03463](https://arxiv.org/abs/1910.03463)
- Brown, Moran and Pearce 1985, Riesz products and normal numbers, J. London Math. Soc. (2) 32, 12-18; cited unread, the publisher page being paywalled.
- Brown, Moran and Pearce 1986, A decomposition theorem for numbers in which the summands have prescribed normality properties, J. Number Theory 24, 259-271; cited unread, the publisher page refusing the fetch. [doi.org/10.1016/0022-314X(86)90034-X](https://doi.org/10.1016/0022-314X(86)90034-X)
- Bugeaud, Cipu and Mignotte 2013, On the representation of Fibonacci and Lucas numbers in an integer base, Ann. Math. Quebec 37, 31-43; read at source in the authors' copy, Theorem 3.1, its definition of `BC` and the statements of Senge and Straus and of Stewart it carries. [doi.org/10.1007/s40316-013-0002-y](https://doi.org/10.1007/s40316-013-0002-y)
- Cassels 1959, On a problem of Steinhaus about normal numbers, Colloq. Math. 7, 95-101; cited unread, the open scan refusing the fetch. [doi.org/10.4064/cm-7-1-95-101](https://doi.org/10.4064/cm-7-1-95-101)
- Coquet, Kamae and Mendes France 1977, Sur la mesure spectrale de certaines suites arithmetiques, Bull. Soc. Math. France 105, 369-384; read at source in its Numdam scan, Theorems 3 to 7. [doi.org/10.24033/bsmf.1856](https://doi.org/10.24033/bsmf.1856)
- Coquet 1979, Correlation de suites arithmetiques, Seminaire Delange-Pisot-Poitou 20, exp. 15, 1-12; read at source in its Numdam scan, Sections I.2, II.5 and II.7 and Theorems 3 and 4. [numdam.org](https://www.numdam.org/item/SDPP_1978-1979__20_1_A11_0/)
- Cobham 1969, On the base-dependence of sets of numbers recognizable by finite automata, Math. Systems Theory 3, 186-192; stated through Bes, Theorem 24, and Durand and Rigo, Theorem 1.1. [doi.org/10.1007/BF01746527](https://doi.org/10.1007/BF01746527)
- Evertse 1984, On sums of S-units and linear recurrences, Compositio Math. 53, 225-244; read at source in its Numdam scan, Theorem 1. [numdam.org](https://www.numdam.org/item/CM_1984__53_2_225_0)
- Evertse, Gyory and Stewart 2018, Mahler's work on Diophantine equations and subsequent developments; read at source, the S-unit statement and its credits. [arxiv.org/abs/1806.00355](https://arxiv.org/abs/1806.00355)
- Herglotz's theorem, positive-definite functions on `Z` as Fourier-Stieltjes transforms; the 1911 paper is unread, the statement is read in the encyclopedia entry. [wikipedia.org](https://en.wikipedia.org/wiki/Bochner%27s_theorem)
- Hochman and Shmerkin 2015, Equidistribution from fractal measures, Invent. Math. 202, 427-479; read at the abstract, context only. [arxiv.org/abs/1302.5792](https://arxiv.org/abs/1302.5792)
- Host 1995, Nombres normaux, entropie, translations, Israel J. Math. 91, 419-428; cited unread, context only. [doi.org/10.1007/BF02761660](https://doi.org/10.1007/BF02761660)
- Kakutani 1972, Strictly ergodic symbolic dynamical systems, Proc. Sixth Berkeley Symp. Math. Statist. Probab. 2, 319-326; cited unread, for the continuity of the Thue-Morse measure.
- Kamae 1977, Mutual singularity of spectra of dynamical systems given by "sums of digits" to different bases, Asterisque 49, 109-114; read at source in its Numdam scan. [numdam.org](https://www.numdam.org/item/AST_1977__49__109_0)
- Lyons 1995, Seventy years of Rajchman measures, J. Fourier Anal. Appl., Kahane special issue, 363-377; read at source in the author's copy, Sections 1, 2 and 5. [rdlyons.pages.iu.edu](https://rdlyons.pages.iu.edu/pdf/seventy.pdf)
- Mahler 1927, The spectrum of an array and its application to the study of the translation properties of a simple class of arithmetical functions, Part Two, J. Math. Phys. 6, 158-163; cited unread, for the Thue-Morse spectral measure. [doi.org/10.1002/sapm192761158](https://doi.org/10.1002/sapm192761158)
- Martin, Odlyzko and Wolfram 1984, Algebraic properties of cellular automata, Comm. Math. Phys. 93, 219-258; cited unread, for rules `90` and `150` as Laurent polynomials. [doi.org/10.1007/BF01223745](https://doi.org/10.1007/BF01223745)
- Nazarov, Peres and Shmerkin 2012, Convolutions of Cantor measures without resonance, Israel J. Math. 187, 93-116; read at source in its arXiv version, the introduction, Section 4 and Theorems 1.1 and 4.1. [arxiv.org/abs/0905.3850](https://arxiv.org/abs/0905.3850)
- Pollington 1981, The Hausdorff dimension of a set of normal numbers, Pacific J. Math. 95, 193-204; read at source in its scan, Theorem 1 and Lemma 1. [doi.org/10.2140/pjm.1981.95.193](https://doi.org/10.2140/pjm.1981.95.193)
- Queffelec 1979, Mesures spectrales associees a certaines suites arithmetiques, Bull. Soc. Math. France 107, 385-421; read at source in its Numdam scan, pages 385-402. [numdam.org](https://www.numdam.org/item/BSMF_1979__107__385_0)
- Queffelec 1987, Substitution Dynamical Systems - Spectral Analysis, Lecture Notes in Math. 1294, Springer; cited unread, as an account of the Thue-Morse spectrum that was not searched.
- van der Poorten and Schlickewei 1982, The growth conditions for recurrence sequences, Macquarie Univ. Math. Rep. 82-0041; cited unread, read through Evertse, Gyory and Stewart 2018.
- Rapaport 2022, On the Rajchman property for self-similar measures on `R^d`, Adv. Math. 403, 108375; Corollary 1.6, read at source in this tree's dimensions topic. [arxiv.org/abs/2104.03955](https://arxiv.org/abs/2104.03955)
- Rudin 1959, Some theorems on Fourier coefficients, Proc. Amer. Math. Soc. 10, 855-859; cited unread, for the Rudin-Shapiro sequence. [doi.org/10.1090/S0002-9939-1959-0116184-5](https://doi.org/10.1090/S0002-9939-1959-0116184-5)
- Schmidt 1960, On normal numbers, Pacific J. Math. 10, 661-672; read at source in the text layer of its scan, Theorems 1 and 2, Lemma 2 and the Fundamental Lemma. [doi.org/10.2140/pjm.1960.10.661](https://doi.org/10.2140/pjm.1960.10.661)
- Schmidt 1961-62, Uber die Normalitat von Zahlen zu verschiedenen Basen, Acta Arith. 7, 299-309; cited unread, carried from Pollington 1981.
- Senge and Straus 1973, PV-numbers and sets of multiplicity, Period. Math. Hungar. 3, 93-100; unread, read through Bugeaud, Cipu and Mignotte 2013 and Nazarov, Peres and Shmerkin 2012. [doi.org/10.1007/BF02018464](https://doi.org/10.1007/BF02018464)
- Shapiro 1951, Extremal problems for polynomials and power series, Master's thesis, MIT; cited unread, for the Rudin-Shapiro sequence.
- Steele, Lea and Flood 2014, Fast splittable pseudorandom number generators, OOPSLA 2014; cited unread, `splitmix64` being the control. [doi.org/10.1145/2660193.2660195](https://doi.org/10.1145/2660193.2660195)
- Stewart 1980, On the representation of an integer in two different bases, J. reine angew. Math. 319, 63-72; read through Bugeaud, Cipu and Mignotte. [doi.org/10.1515/crll.1980.319.63](https://doi.org/10.1515/crll.1980.319.63)
- Varju and Yu 2022, Fourier decay of self-similar measures and self-similar sets of uniqueness, Anal. PDE 15, 843-858; read at the abstract, context only. [arxiv.org/abs/2004.09358](https://arxiv.org/abs/2004.09358)
- Wolfram 1983, Statistical mechanics of cellular automata, Rev. Mod. Phys. 55, 601-644; cited unread, for the single-site diagrams of the elementary rules. [doi.org/10.1103/RevModPhys.55.601](https://doi.org/10.1103/RevModPhys.55.601)
- Wolfram 1986, Random sequence generation by cellular automata, Adv. Appl. Math. 7, 123-169; cited unread. [doi.org/10.1016/0196-8858(86)90028-X](https://doi.org/10.1016/0196-8858(86)90028-X)
- Wolfram 2002, A New Kind of Science, Wolfram Media; cited unread, for the single-site diagrams of the elementary rules.
- Zaks, Pikovsky and Kurths 1997, On the correlation dimension of the spectral measure for the Thue-Morse sequence, J. Stat. Phys. 88, 1387-1392; cited unread, through Zaks 2001. [doi.org/10.1007/BF02732440](https://doi.org/10.1007/BF02732440)
- Zaks 2001, Fractal Fourier spectra in dynamical systems, Habilitationsschrift, University of Potsdam; read at source in the open copy, Section 2.2. [publishup.uni-potsdam.de](https://publishup.uni-potsdam.de/files/145/zaks.pdf)
- OEIS A005942, the factor complexity of Thue-Morse. [oeis.org/A005942](https://oeis.org/A005942)
- OEIS A051023, the centre column of Rule 30. [oeis.org/A051023](https://oeis.org/A051023)
