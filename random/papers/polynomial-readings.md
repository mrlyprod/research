---
title: Normal Is Not Random: Digit Designs Read Along Polynomials
lead: Thue-Morse along the squares is normal, yet `t(4n^2) = t(n^2)`, so one dilation test returns `1` on it; a design with `x(bm) = x(m)` read along a monomial keeps that identity exactly. For a polynomial `P` of degree `k >= 2` and `lambda > 0`, `lambda != 1`, `P(lambda u + mu) - lambda^k P(u)` is constant exactly when `P` is a monomial centred at the fixed point `mu/(1 - lambda)` plus a constant, so for a polynomial reading of degree at least `2` of Thue-Morse, Rudin-Shapiro or the ternary parity at base `b`, a dilation pair `(n, Bn + r)` of a battery has `P(Bn + r) = b^L (P(n) + q) + s`, `0 <= s < b^L`, exactly when `B` is a power of `b` and `P` is a monomial centred at `-r/(B - 1)` plus a constant, and there the mean of the product tends to a local constant read off the `b`-adic digits, such as `-1/2` for Thue-Morse along `n^2 + 4`. The readings of Rudin-Shapiro along `n^2 + n` and `(2n + 1)^2` are automatic in no base, as every reading of Thue-Morse of degree at least `2` is, yet the free readings of Thue-Morse along `n^2 + n` and `n^3 + n` fail the battery at `N = 2^20` to `2^26`, and the Thue-Morse correlations over the multiples of `3` below `2^k` sum to exactly `(16/45)(3/2)^k + 8/15` at even `k >= 2` and `+ 7/15` at odd `k`, a mass at frequency one third that such a pair is conjectured to leak at the rate `N^(-log2(4/3))`. Rudin-Shapiro along `n^2 + n` stays at the control's level in every test to `N = 2^26`, the evidence for the conjecture that its battery correlations decay at the rate of independent signs, and a proof needs a joint normality of `(P(n), P(Bn + r))` that no theorem read here gives.
date: 2026-10-30
figure: paper-polynomial-readings
---

## 1. Introduction

![The three sums `U_0(H)`, `U_1(H)` and `U_2(H)` of the Thue-Morse correlations `C_t(h)` over the lags `0 <= h < H` in each residue class modulo `3`, drawn against `H` for `H <= 2^12` on one vertical scale, the thin line at `0`. The blue walk, the multiples of `3`, climbs to `U_0(2^k) = (16/45)(3/2)^k + 8/15` at even `k >= 2` and `(16/45)(3/2)^k + 7/15` at odd `k` (Theorem 9.1), that is `(16/45) H^0.585` up to a bounded term at `H = 2^k`; the two orange walks, strong for the lags `h = 1 mod 3` and pale for `h = 2 mod 3`, fall and braid, and the three sums total `2/3` at every `2^k`. The figure computes the correlations by their recursion `C_t(2q) = C_t(q)`, `C_t(2q + 1) = -(C_t(q) + C_t(q + 1))/2` and refuses to draw if the closed form or the total fails at any `2^k`, `1 <= k <= 12`.](paper-polynomial-readings)

A digit design read at `n` is an automatic sequence, and an automatic sequence is easy to tell from random: the [companion paper](design-generators.md) reads the base of a design off one dilation test, since `x(bn + r) = +-x(n)` exactly, and assembles the battery of correlation tests used here. A cheap way out of the automatic class is to read the design along a polynomial. Thue-Morse `t(n) = (-1)^(s_2(n))`, `s_2(n)` the number of binary ones of `n`, is normal along the squares (Drmota, Mauduit and Rivat 2019), and so is Rudin-Shapiro (Mullner 2018); a normal sequence has every block of every length, so neither reading is automatic in any base, and the squares have been proposed as a remedy for the low complexity of automatic sequences used as pseudorandom sources (Sun and Winterhof 2019, Popoli 2020, Merai and Winterhof 2022). This paper runs the battery on such readings and finds that normal is not random. Since `(2n)^2 = 4n^2` and `t(4m) = t(m)`,

```
t((2n)^2) = t(n^2)   at every n,
```

so the dilation test, the correlation of `y(n) = t(n^2)` with `y(2n)`, returns exactly `1` on a normal sequence (Theorem 3.1). The question is which polynomial readings escape the battery, and the answer comes in four parts.

**Where the battery sees a reading.** A *reading* is `y(n) = x(P(n))` for a design `x` at base `b` and an integer polynomial `P`; the designs are Thue-Morse and Rudin-Shapiro at `b = 2` and the ternary parity `t3` at `b = 3`, and all three satisfy `x(bm) = x(m)`. A pair of affine maps `(phi, psi)` compared by the battery is a *near-head* of the reading when `P(psi(n)) = b^L (P(phi(n)) + q) + s` as polynomials, `0 <= s < b^L`, so that for a signed design `y(psi(n)) = x(s) x(P(phi(n)) + q)` reads the design a bounded distance from the point `y(phi(n))` reads. The rigidity lemma says that for `P` of degree `k >= 2` and `lambda > 0`, `lambda != 1`, the difference `P(lambda u + mu) - lambda^k P(u)` is constant exactly when `P(u) = c (u - u_0)^k + c_0` with `u_0 = mu/(1 - lambda)` the fixed point (Lemma 4.1). So a dilation pair `(n, Bn + r)` is a near-head exactly when `B` is a power of the base and `P` is a monomial plus a constant centred at `-r/(B - 1)` (Theorem 4.2), and along an integer contraction `n -> lambda n + (1 - lambda) u_0`, which exists exactly when the denominator of the centre `u_0` is prime to the base, a reading of Thue-Morse or `t3` repeats itself up to a sign at every large depth (Theorem 4.4). At a near-head the mean of the product converges to a local constant computed from the `b`-adic distribution of the values, inside an explicit enclosure (Theorem 5.1): `-1/2` for Thue-Morse along `n^2 + 4` at `(2, 0)`, `-1/3` for `t3` along `n^2 + 1` at `(3, 0)`, and the pair law of a design along the squares, the limit of the mean of `x(n^2 + c) x(n^2)`, which for Thue-Morse reads `0, -1, 1, -1/2, -1/2, 1/2` at `c = 1, ..., 6` and differs from the Thue-Morse correlations at each of the eight offsets computed. Every other pair is *free*.

**Out of the automatic class.** A `b`-automatic sequence whose kernel has `kappa` elements has at most `b k abs(A)^(2 kappa)` factors of length `k` (Proposition 7.1, a classical bound), and two explicit families of kernel elements show that `t(n^2)` is not `2`-automatic and `t3(n^2)` is not `3`-automatic (Theorem 7.2). Every reading of Thue-Morse along a polynomial of degree at least `2`, and the readings of Rudin-Shapiro along `n^2`, `n^2 + 1`, `2n^2 + 1`, `n^2 + n` and `(2n + 1)^2`, are automatic in no base (Theorem 7.3). For Thue-Morse and for `r(n^2)` this is in print; the readings of Rudin-Shapiro along `n^2 + n` and `(2n + 1)^2` follow from Mullner's normality and an identity between readings, `r((2n + 1)^2) = r(n^2 + n)` (Proposition 6.1).

**What the battery still sees.** At `N = 2^15` every exact near-head reads `1.0000` where Theorem 4.2 puts it, and the local constants reappear as the dilation maxima `0.5002`, `0.3333` and `0.3334` (Fact 8.1). The battery reads depth `1` of a contraction only at the centre `0`: Thue-Morse along `(n + 1)^2 + 1`, centred at `-1`, reads inside the band in its dilation and kernel tests while `y(2n + 1) = y(4n + 3)` at every `n` (Fact 8.2). Freedom is not enough either: Thue-Morse along `n^2 + n`, free at every battery pair, fails the dilation test at `N = 2^20`, `2^22`, `2^24` and `2^26` with `z` from `9.48` to `14.42`, and Thue-Morse and Rudin-Shapiro along `n^3 + n` fail the kernel test (Remark 8.3).

**The mass at one third.** The sums `U_r(H)` of the Thue-Morse correlations `C_t(h)` over `0 <= h < H`, `h = r mod 3`, obey a linear recursion in `H -> 2H` with eigenvalues `0`, `3/2` and `-3/2`, and

```
sum of C_t(h) over 3 | h < 2^k  =  (16/45)(3/2)^k + 8/15    (k >= 2 even)
                                =  (16/45)(3/2)^k + 7/15    (k odd)          (Theorem 9.1),
```

so the correlations at the multiples of `3` have mean `(16/15)(2^k)^(-log2(4/3)) + O(2^(-k))` against `2/(3 . 2^k)` over all lags: a mass at frequency one third, where every factor `1 - cos(2 pi 2^j theta)` of the Thue-Morse Riesz product equals `3/2`. The largest correlations of Thue-Morse along `n^2 + n` sit at the pairs `(4, 0)`, `(4, 3)` and `(8, 5)`, whose remainders `P(Bn + r) - B^2 P(n)` are `-12n`, `12n + 12` and `24n + 30`, all with slopes divisible by `3`, and in the runs behind Conjecture 11.1 they decay with slope `-0.399` in `log2 N` from `2^20` to `2^26`, where the noise has `-0.5`; the conjecture is that such a pair leaks the mass at the rate `N^(-log2(4/3)) = N^(-0.415...)`. Rudin-Shapiro, whose correlations vanish at every nonzero lag, stays at the control's level in every test along `n^2 + n` to `N = 2^26`, the evidence for Conjecture 11.2: every battery correlation of that reading is `O(N^(-1/2 + epsilon))`, the rate of independent signs. The obstruction to a proof is a joint normality of `(P(n), P(Bn + r))` that no theorem read here gives (Problem 11.3).

**What is classical, and what is not.** The correlation recursion of Thue-Morse is classical: Baake and Coons 2024, read at source, credit it to Mahler 1927, Kakutani 1972 and Baake and Grimm 2013, all three cited unread, and the companion paper proves it (Lemma 2.5 there); the vanishing correlations of Rudin-Shapiro go back to Rudin 1959 and Shapiro 1951, cited unread, and are reproved there (Theorem 3.5). The linear factor complexity of automatic sequences is classical, Corollary 10.3.2 of Allouche and Shallit 2003 as Popoli 2020 quotes it; Proposition 7.1 proves it in kernel terms because that book is cited unread. The normality theorems are used as stated, and much of Theorem 7.3 is in print: Allouche and Shallit 2003, Theorem 6.10.1, makes Thue-Morse along a polynomial non-automatic exactly when the degree is at least `2`, as Popoli 2020 quotes it in its Section 2, and for the squares Sun and Winterhof 2019 in their Section 1; Shen 2022, read at source, gives every block of `t(H(n))` for every polynomial `H` of degree at least `2`, which with Proposition 7.1 excludes every base; and Popoli 2020, in its Section 1, records `r(n^2)` as non-automatic, through Mullner's normality and Christol's theorem. Lemma 4.1 is an exercise in centring a polynomial and Theorem 5.1 a compactness argument in the `b`-adic integers. Not found in the sources read, which the References list with what was read of each: a test that compares a polynomial reading at `n` with the same reading at `Bn + r`, the rigidity criterion as a statement about such tests (Theorems 4.2 and 4.4), the local constants with their explicit enclosures (Theorem 5.1), the explicit kernel witnesses of Theorem 7.2, the witness for `t3` among them, the non-automaticity of the Rudin-Shapiro readings along `n^2 + n` and `(2n + 1)^2`, and the closed form of Theorem 9.1. That the pair-law limits exist and are rational at a prime base follows from Adamczewski, Drmota and Mullner 2022, and the literature on the Thue-Morse correlations was searched only in part, so the closed form may be in print. Section 10 gives the details.

## 2. Readings and the battery

Throughout, a *stream* is a map `x : N -> {-1, 1}`, `N = {0, 1, 2, ...}`, `e(u) = exp(2 pi i u)`, and `m_i` is the base-`b` digit of `m` at place `i`, `m_0` the last. A stream is *normal* when, for every `k`, each of the `2^k` blocks of `k` signs occurs at `n, n + 1, ..., n + k - 1` for a share `2^(-k)` of the starting points `n`.

**Definition 2.1 (designs).** Three designs are read:

- *Thue-Morse* `t(n) = (-1)^(s_2(n))`, `s_2(n)` the number of binary ones of `n`;
- *Rudin-Shapiro* `r(n) = (-1)^(s_11(n))`, `s_11(n)` the number of adjacent pairs `11` in the binary word of `n`;
- the *ternary parity* `t3(n) = (-1)^(d_2(n))`, `d_2(n)` the number of ternary digits of `n` equal to `2`.

Thue-Morse is a *signed design* at `b = 2` and `t3` one at `b = 3`: for `L >= 1`, `m >= 0` and `0 <= s < b^L`,

```
x(b^L m + s) = x(m) x(s),                                    (2.1)
```

since the digits of `b^L m + s` are those of `s`, padded with zeros to `L` digits, followed by those of `m`. Rudin-Shapiro reads one pair across the boundary,

```
r(2^L m + s) = r(m) r(s) (-1)^(m_0 s_(L-1)),                 (2.2)
```

so `r(2m) = r(m)` and `r(2m + 1) = (-1)^m r(m)`. All three designs satisfy `x(bm) = x(m)`, with `b = 2` for `t` and `r` and `b = 3` for `t3`; this `b` is the *base* of the design. The *correlation* of a stream at lag `h >= 0` is `C_x(h) = lim_(N -> infinity) (1/N) sum_(n < N) x(n) x(n + h)` when the limit exists. For Thue-Morse it exists at every lag, with `C_t(0) = 1`, `C_t(1) = -1/3` and

```
C_t(2q) = C_t(q),   C_t(2q + 1) = -(C_t(q) + C_t(q + 1))/2,  (2.3)
```

and Rudin-Shapiro has `C_r(h) = 0` at every `h >= 1` (the companion paper, Lemma 2.5 and Theorem 3.5).

**Definition 2.2 (readings).** A *polynomial reading* of a design `x` is the stream `y(n) = x(P(n))` for an integer polynomial `P` with `P(n) >= 0` at every `n >= 0`, and the reading *along the primes* is `y(n) = x(p_n)`, `p_0 = 2 < p_1 = 3 < ...` the primes in order. The readings run here are along the eight sequences `n`, `n^2`, `n^2 + 1`, `n^2 + n`, `(2n + 1)^2`, `2n^2 + 1`, `n^3` and `p_n`, and along six probes that test the edges of the statements: `n^2 + 2`, `n^2 + 3`, `n^2 + 4`, `n^2 + 5`, `(n + 1)^2 + 1` and `n^3 + n`.

**Definition 2.3 (the battery).** On a window of `N` pairs, with Pearson's correlation coefficient written `Pearson`:

- the *shift test* is the largest `abs Pearson((y(n), y(n + h)) : n < N)` over `1 <= h <= 64`;
- the *dilation test* is the largest `abs Pearson((y(n), y(Bn + r)) : n < N)` over the `44` pairs `2 <= B <= 9`, `0 <= r < B`;
- the *kernel test* is the largest `abs Pearson((y(bn), y(b^2 n)) : n < N)` over `b = 2, 3`;
- the *dilation square* is the mean of `y(n) y(2n) y(3n) y(6n)` over `n < N`;
- the *LZ78 rate* of the first `N` terms, read as bits, is `c (log2 c + 1)/N` bits per bit, `c` the number of phrases of their LZ78 parse.

For a stream of independent uniform signs each correlation and the square have mean `0` and standard deviation `1/sqrt N` once the points they multiply are distinct; a value past `4/sqrt N`, `0.0221` at `N = 2^15`, is a failure, and `z = max abs(Pearson) sqrt N` measures a test on that scale. The control is `splitmix64` (Steele, Lea and Flood 2014). This is the protocol of the companion paper (Fact 3.6 there), run on readings. A *battery pair* is one of the `46` pairs of affine maps `(phi, psi)` that the dilation and kernel tests compare: `(n, Bn + r)` and `(bn, b^2 n)`.

**Definition 2.4 (near-heads, exact and twisted pairs, contractions).** Let `y = x o P` be a polynomial reading of a design at base `b`.

- A battery pair `(phi, psi)` is a *near-head* of the reading when `P(psi(n)) = b^L (P(phi(n)) + q) + s` as polynomials in `n`, for integers `L >= 1`, `q` and `0 <= s < b^L`. A battery pair that is not a near-head is *free*, and a reading is free when no battery pair is a near-head; the reading along the primes is called free.
- A pair of maps `(phi, psi)` is *exact* when `y(psi(n)) y(phi(n))` is constant in `n`, and *twisted* when it is periodic in `n` and not constant.
- A *contraction* is an integer affine map `phi(n) = lambda n + mu` with `lambda >= 2`. The *head scan* reads `y(phi^j(n)) y(phi^(j+1)(n))` over `n < N` for `lambda` in `{2, 3, 4, 8, 9}`, `0 <= mu < lambda` and `j = 0, 1, 2`, and reports the products that are periodic in `n` with period at most `81`: plain at period `1`, twisted above. The dilation test at `(B, r)` reads the contraction `Bn + r` at depth `j = 0`, and the kernel test at `b` the contraction `bn` at depth `j = 1`.

## 3. Monomial readings

**Theorem 3.1 (the dilation identity of monomial readings).** Let `x` be a stream with `x(bm) = x(m)` at every `m >= 0`, and let `P(n) = c n^k` with integers `c >= 1` and `k >= 1`. Then the reading `y(n) = x(P(n))` satisfies `y(bn) = y(n)` at every `n >= 0`. So the dilation test at `(b, 0)` and the kernel test at `b` compare `y` with itself and return `1` on every window on which `y` is nonconstant, and the dilation square equals `1` at every `n`, at `b = 2` and at `b = 3`.

*Proof.* `P(bn) = b^k P(n)`, and `x(bm) = x(m)` applied `k` times gives `y(bn) = y(n)`, hence also `y(b^2 n) = y(bn)`. At `b = 2` the factors of the square pair as `y(2n) = y(n)` and `y(6n) = y(3n)`, and at `b = 3` as `y(3n) = y(n)` and `y(6n) = y(2n)`. □

So Thue-Morse and Rudin-Shapiro along every monomial keep `y(2n) = y(n)`, and `t3` along every monomial keeps `y(3n) = y(n)`. For `t(n^2)` this is an identity on a normal sequence whose shift correlations all vanish: Theorem 2 of the preprint of Drmota, Mauduit and Rivat bounds `sum_(n < N) e((1/2) sum_(l < k) alpha_l s_2((n + l)^2))` by `N^(1 - eta)` for every nonzero `alpha` in `{0, 1}^k`, so every shift correlation of `t(n^2)` tends to `0` at a power rate, with an exponent `eta > 0` that is never given a value. The battery at `N = 2^15` reads `t(n^2)` as mean `-0.0150`, shift `0.0280` at `h = 15`, dilation `1.0000` at `(2, 0)`, kernel `1.0000` and square `+1.0000`, and its LZ78 rate is `1.214` bits per bit at `2^18`, against the control's `1.215` (Fact 8.1). The compressor sees a normal sequence, as it must; one dilation sees the identity. The measures used where the squares are proposed as a pseudorandom source read consecutive terms and their translations only: the maximum order complexity, and the correlation measures of Mauduit and Sarkozy 1998, which compare `x(n + d_1)` with `x(n + d_2)`. A test that compares `y(n)` with `y(2n)` is not found in the sources read (Section 10).

## 4. The rigidity criterion

**Lemma 4.1 (rigidity).** Let `P` be a polynomial of degree `k >= 2`, and `chi(u) = lambda u + mu` with `lambda > 0` rational, `lambda != 1`, and fixed point `u_0 = mu/(1 - lambda)`. Then `P(chi(u)) - lambda^k P(u)` is constant if and only if `P(u) = c (u - u_0)^k + c_0` for constants `c != 0` and `c_0`, and the constant is then `c_0 (1 - lambda^k)`.

*Proof.* Write `Q(w) = P(u_0 + w) = sum_(i <= k) c_i w^i`. Since `chi(u_0 + w) = u_0 + lambda w`, putting `w = u - u_0`,

```
P(chi(u)) - lambda^k P(u) = Q(lambda w) - lambda^k Q(w) = sum_(i <= k) c_i (lambda^i - lambda^k) w^i,
```

and `lambda^i != lambda^k` for `i != k` because `lambda > 0` and `lambda != 1`. So the difference is constant exactly when `c_i = 0` for `0 < i < k`, and its value is then `c_0 (1 - lambda^k)`. □

Such a `P` is a *centred monomial plus a constant*, with *centre* `u_0`. Every polynomial of degree `1` has this form, which is why the lemma asks for `k >= 2`.

**Theorem 4.2 (near-heads of the battery).** Let `b` be prime and `y = x o P` a polynomial reading of degree `k >= 2` of a design at base `b`. A battery pair `(phi, psi)` is a near-head of `y` exactly when `P` is a centred monomial plus a constant whose centre is the fixed point of `chi = psi o phi^(-1)` and whose slope ratio `lambda`, the slope of `psi` over the slope of `phi`, has `lambda^k` a power of `b`. For the dilation pair `(n, Bn + r)` this asks for the centre `-r/(B - 1)` and `B` a power of `b`; for the kernel pair `(b'n, b'^2 n)` it asks for `b' = b` and the centre `0`.

*Proof.* Write `phi(n) = alpha n + beta` and `psi(n) = gamma n + delta`, so that `chi(u) = (gamma/alpha)(u - beta) + delta` has slope `lambda = gamma/alpha`, which is `B` or `b'`, never `1`. If `(phi, psi)` is a near-head, `P(psi(n)) - b^L P(phi(n)) = b^L q + s` is constant, its coefficient of `n^k` vanishes, and `lambda^k = b^L`. Since `phi` is an affine bijection of the rationals, the identity in `n` is the identity `P(chi(u)) - lambda^k P(u) = b^L q + s` in `u = phi(n)`, and Lemma 4.1 makes `P` a centred monomial plus a constant centred at the fixed point of `chi`. Conversely, if `P` is centred at that fixed point and `lambda^k = b^L`, Lemma 4.1 gives `P(psi(n)) - b^L P(phi(n)) = c_0 (1 - b^L)`, a constant integer polynomial in `n`, hence an integer, and writing it as `b^L q + s` with `0 <= s < b^L` gives a near-head; `L >= 1` since `lambda > 1`. For the dilation pair, `chi(u) = Bu + r` has fixed point `-r/(B - 1)` and `lambda = B`, and `B^k` is a power of the prime `b` exactly when `B` is; for the kernel pair, `chi(u) = b'u` has fixed point `0` and `lambda = b'`. □

**Corollary 4.3 (the centres the battery visits).** A polynomial reading of degree `k >= 2` has a near-head among the battery pairs exactly when it is a centred monomial plus a constant whose centre is one of `0`, `-1/3`, `-2/3`, `-1` and `-r/7`, `0 < r < 7`, at base `2`, or one of `-r/8`, `0 <= r <= 8`, at base `3`. Every other polynomial reading of degree at least `2` is free.

*Proof.* Theorem 4.2 with `B` in `{2, 4, 8}` at base `2` and in `{3, 9}` at base `3`; the kernel pair adds the centre `0`. □

The centre decides, not the degree. `n^2 + n = (n + 1/2)^2 - 1/4` is centred at `-1/2 = -4/8`: a near-head of `t3` at `(3, 1)` and `(9, 4)` and free for `t` and `r`, and so is `(2n + 1)^2 = 4 (n + 1/2)^2`. `(n + 1)^2 + 1` is centred at `-1`, a near-head at `(2, 1)`, `(4, 3)` and `(8, 7)` for `t` and `r` and at `(3, 2)` and `(9, 8)` for `t3`. `n^3 + n` is not a centred monomial plus a constant, since its coefficients would need `3 u_0 = 0` and `3 u_0^2 + 1 = 0`, so it is free for all three designs.

**Theorem 4.4 (identities along a contraction).** Let `b` be prime and `P` an integer polynomial of degree `k >= 2` with `P(n) >= 0` for `n >= 0` that is a centred monomial plus a constant, `P(u) = c (u - u_0)^k + c_0`, with centre `u_0 = p/d` in lowest terms. An integer contraction `phi(n) = lambda n + (1 - lambda) u_0` with `lambda` a power of `b` exists exactly when `d` is prime to `b`. Along such a `phi`, the reading `y` of `t` at `b = 2`, or of `t3` at `b = 3`, satisfies for every large `j`

```
y(phi^(j+1)(n)) = eps_j y(phi^j(n))   at every integer n >= max(0, u_0),
```

with a sign `eps_j` that does not depend on `n`.

*Proof.* The map `phi` has integer coefficients exactly when `(lambda - 1) p/d` is an integer, that is when `d` divides `lambda - 1`. With `lambda = b^i`, every multiple `i >= 1` of the order of `b` modulo `d` achieves this when `d` is prime to `b`, and no `i` does when `b` divides `d`, since `b^i - 1` is then prime to `b`. Fix such a `phi` and put `beta = lambda^k = b^(ik)`. Since `phi(u) - u_0 = lambda (u - u_0)`, `P(phi^j(u)) = beta^j (P(u) - c_0) + c_0`, and for `n >= max(0, u_0)`, `phi^j(n) = u_0 + lambda^j (n - u_0) >= n >= 0`, so `m_j = P(phi^j(n))` is a value of the reading, `m_j >= 0`. Here `m_j = beta^j P(n) + c_0 (1 - beta^j)`, and `c_0 (1 - beta^j) = P(phi^j(n)) - beta^j P(n)` is an integer, so `m_j` is congruent modulo `beta^j` to the residue `c_j` of that integer, which does not depend on `n`. Also `m_(j+1) = beta m_j + c_0 (1 - beta)`, and writing the integer `c_0 (1 - beta)` as `beta q + s` with `0 <= s < beta`,

```
m_(j+1) = beta (m_j + q) + s,
```

with `q` and `s` independent of `n` and `j`. The number `c_0` is a `b`-adic integer, since `c_0 (1 - beta)` is an integer and `1 - beta` is a `b`-adic unit, and `c_j` is its residue modulo `beta^j`. If `c_0` is a nonnegative integer, `c_j + q = floor(c_0/beta)` once `beta^j > c_0`; if it is a negative integer, `c_j + q = beta^j + floor(c_0/beta)` once `beta^j > -c_0`; otherwise its `b`-adic digits are eventually periodic and neither eventually `0` nor eventually `b - 1`, so `c_j/beta^j` stays in a closed subinterval of `(0, 1)` for large `j`. In each case `0 <= c_j + q < beta^j` for every large `j`. Then `m_j = beta^j A + c_j` and `m_j + q = beta^j A + (c_j + q)` with `A >= 0`, and (2.1), applied at `L = ik` and at `L = ikj`, gives

```
y(phi^(j+1)(n)) y(phi^j(n)) = x(s) x(m_j + q) x(m_j) = x(s) x(c_j + q) x(c_j) = eps_j.  □
```

**Remark 4.5 (what the kernel test sees).** The kernel test at `b` is the case `u_0 = 0`, `lambda = b`, `j = 1`. For Thue-Morse along `n^2 + 1`, `c_0 = 1`, `beta = 4`, `q = -1`, `s = 1` and `c_1 + q = 0`, so the identity holds at `j = 1` and the kernel test reads `1`; along `n^2 + 5`, `q = -4`, `s = 1`, `c_1 + q = -3` and `c_2 + q = 1`, so the identity starts at `j = 2` and the kernel test misses it; and `(n + 1)^2 + 1`, centred at `-1`, needs the contraction `2n + 1`, along which `y(4n + 3) = y(2n + 1)` already at `j = 1` (Fact 8.2). The head scan finds all three (Fact 8.1, last column).

## 5. Local constants

At a near-head the reading is not free, and the mean of its product has a limit that the low digits decide.

**Theorem 5.1 (local constants).** Let `x` be `t` or `t3` at its base `b`, or `r` at `b = 2`, let `F` be a nonconstant integer polynomial with `F(n) >= 0` for `n >= 0`, and let `q` be an integer; `x` is read as `0` at the finitely many negative arguments `F(n) + q`. Then the mean of `x(F(n) + q) x(F(n))` over `n < N` converges as `N -> infinity`. For `J >= 1`, call a class of `n` modulo `b^J`, modulo `2^(J+1)` for `r`, *good* when `a = F(n) mod b^J` has `0 <= a + q < b^J`. On a good class the product is the constant `x(a + q) x(a)`, times `(-1)^(h ((a + q)_(J-1) + a_(J-1)))` for `r`, `h` the binary digit of `F(n)` at place `J`. Let `kappa_J` be the sum of these constants over the good classes divided by the number of classes, and `eps_J` the share of the classes that are not good. Then the limit lies in `[kappa_J - eps_J, kappa_J + eps_J]` for every `J`, and `eps_J -> 0`.

*Proof.* On a good class `F(n) = b^J A + a` with `A >= 0` and `F(n) + q = b^J A + (a + q)` with `0 <= a + q < b^J`, so the addition does not reach the digits above the window: (2.1) gives `x(F(n) + q) x(F(n)) = x(a + q) x(a)`, and (2.2) adds for `r` the factor `(-1)^(A_0 ((a + q)_(J-1) + a_(J-1)))` with `A_0 = h`. Since `F` has integer coefficients, `F(n) mod b^J` and `h` depend only on `n` modulo the period, `b^J` or `2^(J+1)`. Over a whole number of periods the good classes contribute exactly `kappa_J` times the number of terms and each other term has absolute value at most `1`, so the mean lies in `[kappa_J - eps_J, kappa_J + eps_J]`; over any `N` it lies within `2 b^(J+1)/N` of that interval, so every limit point of the means lies in it, for every `J`. A class that is not good has `F(n) + e = 0` modulo `b^J` for one of `abs(q)` integers `e`. In the `b`-adic integers the sets `{n : F(n) + e = 0 mod b^J}` are open and closed and decrease with `J`, and, `b` being prime, their intersection is the finite set of roots of the nonzero polynomial `F + e`, of Haar measure `0`. So their measures, the shares of bad classes, tend to `0`, the intervals shrink to one point, and the means converge to it. □

**Corollary 5.2 (near-head limits).** Let `(phi, psi)` be a near-head of a polynomial reading `y = x o P`, `P(psi(n)) = b^L (m + q) + s` with `m = P(phi(n))`. For `x = t` or `t3` the raw mean of `y(psi(n)) y(phi(n))` over `n < N` tends to `x(s) kappa`, `kappa` the limit of Theorem 5.1 for `F = P o phi` and this `q`. For `x = r` it tends to the limit of the same class sums with the boundary factor `(-1)^((m + q)_0 s_(L-1))` of (2.2) inside them, which the same proof bounds; that factor varies with `n` when `s_(L-1) = 1`. When `eps_J = 0` from some `J` on, every class is good and the product is periodic in `n`, so the pair is exact or twisted. The Pearson coefficient of the battery has the same limit when the means of `y(phi(n))` and `y(psi(n))` both tend to `0`. That is known for `t(n^2)` and `r(n^2)` by their normality and for `t(n^3)` by Spiegelhofer 2023, Theorem 2.1, and is not proved here for `t3` along the squares, for `t` along `n^3 + c` or for the probes along `n^2 + c`, `c != 0`, whose limits are statements about raw means.

*Proof.* By (2.1), `y(psi(n)) = x(s) x(m + q)` for a signed design, and by (2.2), `r(2^L (m + q) + s) = r(m + q) r(s) (-1)^((m + q)_0 s_(L-1))`; here `m + q >= 0`, since `b^L (m + q) = P(psi(n)) - s > -b^L`. On a good class modulo `b^J`, `J >= 1`, the digit `(m + q)_0 = (a + q)_0` is constant, so the factor joins the class constant. The Pearson coefficient of two sign sequences is the mean of their product minus the product of their means, divided by the square roots of `1 - mean^2`. □

The boundary factor matters: Rudin-Shapiro along `n` has the near-head `(2, 1)` with `s = 1`, and `r(2n + 1) r(n) = (-1)^n` has mean `0`, not `r(1) = 1`. The study `readings`, verb `local`, takes every near-head battery pair of the thirteen polynomial readings of Definition 2.2 under the three designs, `143` pairs, and computes `kappa_J` and `eps_J` at `J = 20` (`J = 13` at base `3`) beside the raw mean at `N = 2^22` (`3^14` at base `3`), a whole number of periods, at which the raw mean must lie inside the enclosure; every one does. Of the `143` pairs, `91` are exact, `+-1` with `eps_J = 0`, `23` have the limit `0` with `eps_J = 0`, and `29` have a limit of absolute value strictly between `0` and `1`. The census counts polynomial readings, not streams: by Proposition 6.1, Thue-Morse along `n^2`, `n^2 + 2`, `n^2 + 3` and `2n^2 + 1` is one stream up to sign, and along `n^2 + 1` the same stream twisted by `(-1)^n`, yet each polynomial's pairs are counted. For the readings of Definition 2.2:

| reading | centre | `t` | `r` | `t3` |
|---|---|---|---|---|
| `n^2`, `n^3` | `0` | exact at `(2, 0)`, `(4, 0)`, `(8, 0)` and the kernel pair | the same | exact at `(3, 0)`, `(9, 0)` and the kernel pair |
| `n^2 + 1` | `0` | `0` at the dilations, exact at the kernel pair | exact at all four | `-1/3` at the dilations, exact at the kernel pair |
| `2n^2 + 1` | `0` | exact at all four | `0` at the dilations, exact at the kernel pair | `0.333333` at the dilations, exact at the kernel pair |
| `n^2 + n`, `(2n + 1)^2` | `-1/2` | free | free | exact at `(3, 1)` and `(9, 4)` |
| `p_n` | none | free | free | free |

The entry `0.333333` is an enclosure, `eps_J <= 1.26e-6`, not a proved rational; every other entry has `eps_J = 0`. The linear reading `n` of `t` and `t3` is exact at every dilation by a power of its base, and `r` along `n` is exact at `8` of its `15` near-head pairs and twisted with mean `0` at the other `7`, since `r(2m + 1) = (-1)^m r(m)`. Two probes show the constants away from `+-1` and `0`:

- Thue-Morse along `n^2 + 4`: `P(2n) = 4 (P(n) - 3)`, so the dilation `(2, 0)` compares `t(n^2 + 4)` with `t(4 (n^2 + 1)) = t(n^2 + 1)`. For odd `n`, `n^2 = 1 mod 8`, and `n^2 + 1` and `n^2 + 4` end in `010` and `101` below the same higher digits, so the product is `-1`; for `n = 2a`, `t(4a^2 + 1) t(4a^2 + 4) = -t(a^2) t(a^2 + 1) = (-1)^a` by Proposition 6.1. The product is `+1` at `n = 0 mod 4` and `-1` otherwise, a twisted pair of period `4` with mean `-1/2`, and the battery prints `0.5002` at `(2, 0)` (Fact 8.1).
- `t3` along `n^2 + 5` reads `0.555556` at `(3, 0)` with `eps_J <= 2.51e-6` at `J = 13`: the carries leave the window on a set of classes whose share shrinks with `J` and has not vanished at `J = 13`, and the head scan finds no periodic identity at depth `0`.

Theorem 5.1 with `F(n) = n^2`, run by `local` at `J = 26` digits (`16` at base `3`), gives the *pair law* of a design along the squares, the limit of the mean of `x(n^2 + c) x(n^2)`:

| `c` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| `t` | `0` | `-1` | `1` | `-1/2` | `-1/2` | `1/2` | `1/6` * | `-2/3` * |
| `t3` | `-1/3` | `-1/3` * | `1/3` * | `1/3` * | `-5/9` * | `-5/9` * | `-1/9` * | `1/3` * |
| `C_t(c)` | `-1/3` | `-1/3` | `1/3` | `-1/3` | `0` | `1/3` | `0` | `-1/3` |

The cells without a star are exact, `eps_J = 0`, and `local` prints them to nine decimals. A starred cell is not proved equal to the rational `p/q` shown: `local` prints the bound `abs(kappa_J - p/q) + eps_J` on the distance from the limit to it, at most `1.40e-7`, and `7.95e-8` for `t` at `c = 7` and `8`. The last row is Thue-Morse's own correlation, from (2.3). Along `floor(n^c)`, `1 < c < 7/5`, the pairs of Thue-Morse keep the Thue-Morse law, `-1/3` at offset `1` (Deshouillers, Drmota and Morgenbesser 2012). Along the squares they do not, at any of the eight offsets: a square is `0` or `1` modulo `4`, so its binary digit at place `1` is always `0` and adding `2` always adds one binary one. The squares are not equidistributed modulo powers of `2`, and the pair law is the Thue-Morse law read through their `2`-adic distribution.

## 6. Identities between readings

Several readings are one stream up to a sign or a periodic twist.

**Proposition 6.1 (identities between readings).** At every `n >= 0`:

1. `t(n^2 + 1) = -(-1)^n t(n^2)`, `t(n^2 + 2) = -t(n^2)`, `t(n^2 + 3) = t(n^2)` and `t(2n^2 + 1) = -t(n^2)`;
2. `t((2n + 1)^2) = -t(n^2 + n)`;
3. `r(n^2 + 1) = r(n^2)`, `r(2n^2 + 1) = (-1)^n r(n^2)` and `r((2n + 1)^2) = r(n^2 + n)`;
4. `t3(n^2 + 1) = -t3(n^2)` when `3` does not divide `n`, and `t3(n^2 + 1) = t3(n^2)` when it does;
5. `t3` along `n^2 + n` has `y(3n + 1) = -y(n)`, along `(2n + 1)^2` it has `y(3n + 1) = y(n)`, and along `2n^2 + 1` it has `y(3n) = y(9n) = t3(2n^2)`.

*Proof.* Each is one look at the low digits. (1) A square is `0` modulo `4` or `1` modulo `8`. For even `n`, adding `1` to `n^2` sets the digit at place `0`; for odd `n`, `...001 + 1 = ...010` keeps the number of ones. The digit of a square at place `1` is `0`, so adding `2` sets it. Adding `3` sets two ones when `n` is even, and turns `8a + 1` into `8a + 4`, one one for one, when `n` is odd. The binary word of `2n^2 + 1` is that of `n^2` followed by a `1`. (2) `(2n + 1)^2 = 4 (n^2 + n) + 1`. (3) `...00 + 1` and `...001 + 1 = ...010` create and destroy no pair `11`; `r(2m + 1) = (-1)^m r(m)` at `m = n^2`, with `n^2 = n mod 2`; and `r(4m + 1) = r(m)` by (2.2). (4) The last ternary digit of a square is `1` when `3` does not divide `n` and `0` when it does, and adding `1` turns a `1` into a `2` and a `0` into a `1`. (5) `(3n + 1)(3n + 2) = 9 (n^2 + n) + 2`, `(6n + 3)^2 = 9 (2n + 1)^2`, `2 (3n)^2 + 1 = 9 (2n^2) + 1` and `2 (9n)^2 + 1 = 81 (2n^2) + 1`, with (2.1). □

The study `readings`, verb `identities`, checks every identity at every `n < 2^24` and finds `0` violations. So the eight readings of Thue-Morse, with the probes `n^2 + 2` and `n^2 + 3`, are five streams up to sign and twist: along `n`, along the squares (with `n^2 + 1`, `n^2 + 2`, `n^2 + 3` and `2n^2 + 1`), along `n^2 + n` (with `(2n + 1)^2`), along `n^3` and along the primes. The squares stream is the `n^2 + n` stream spread over the odd parts, `t((2^a (2m + 1))^2) = -t(m^2 + m)`, and for Rudin-Shapiro `n^2 + 1` collapses onto `n^2` exactly.

## 7. Out of the automatic class

The `b`-kernel of a sequence `x` is the set of its subsequences `n -> x(b^k n + r)`, `k >= 0`, `0 <= r < b^k`, and `x` is `b`-automatic when its `b`-kernel is finite, as in the companion paper (Definition 2.3 there).

**Proposition 7.1 (linear complexity).** A `b`-automatic sequence over an alphabet `A` whose `b`-kernel has `kappa` elements has at most `b k abs(A)^(2 kappa)` factors of length `k`, for every `k >= 1`.

*Proof.* Take `j` with `b^(j-1) < k <= b^j`. A factor of length `k` starting at `b^j q + rho`, `0 <= rho < b^j`, lies in `[b^j q, b^j (q + 2))`, and `x(b^j q + s)` for `0 <= s < b^j` is the kernel element `n -> x(b^j n + s)` read at `q`. So the factor is fixed by `rho` and by the two vectors of the values of all kernel elements at `q` and at `q + 1`, which leaves at most `b^j abs(A)^(2 kappa) < b k abs(A)^(2 kappa)` factors. □

**Theorem 7.2 (kernel witnesses).** `t(n^2)` is not `2`-automatic and `t3(n^2)` is not `3`-automatic.

*Proof.* Since `(2^k n + 1)^2 = 2^(k+1) (2^(k-1) n^2 + n) + 1`, (2.1) makes the kernel element `n -> t((2^k n + 1)^2)` equal to `g_k(n) = -t(2^(k-1) n^2 + n)`. Fix `k >= 3` and `n = 2^(k-1) + 1`. Then `2^(k-1) n^2 + n = 2^(3k-3) + 2^(2k-1) + 2^k + 1` has four binary ones, so `g_k(n) = -1`. For `k' > k`, `n < 2^(k'-1)`, so the binary word of `2^(k'-1) n^2 + n` is that of `n^2 = 2^(2k-2) + 2^k + 1` above that of `n`, `3 + 2 = 5` ones, and `g_(k')(n) = 1`. So the kernel elements `g_k`, `k >= 3`, are pairwise distinct, and the `2`-kernel is infinite. In base `3`, (2.1) and `t3(1) = 1` make `n -> t3((3^k n + 1)^2)` equal to `e_k(n) = t3(3^k n^2 + 2n)`. At `n = (3^k + 1)/2`, `2n = 3^k + 1` and `e_k(n) = t3(3^k (n^2 + 1) + 1) = t3(n^2 + 1)`, while for `k' > k`, `2n < 3^(k')` and `e_(k')(n) = t3(n^2) t3(3^k + 1) = t3(n^2)`. Since `n = 2 mod 3`, `3` does not divide `n` and `t3(n^2 + 1) = -t3(n^2)` by Proposition 6.1. So the `e_k`, `k >= 1`, are pairwise distinct. □

The verb `identities` checks the two kernel links at every `k <= 20` and `n < 1024`, and the witnesses split all `703` pairs `3 <= k < k' <= 40` and all `276` pairs `1 <= k < k' <= 24`.

**Theorem 7.3 (automatic in no base).** Every reading of `t` along a polynomial of degree at least `2`, and the readings of `r` along `n^2`, `n^2 + 1`, `2n^2 + 1`, `n^2 + n` and `(2n + 1)^2`, are automatic in no base.

*Proof.* A sequence over `{-1, 1}` with all `2^k` factors of length `k` for every `k` exceeds the bound of Proposition 7.1 for large `k` in every base, so it is automatic in no base. For `t` along a polynomial `H` of degree at least `2` with `H(N)` in `N` this is Theorem 1 of Shen 2022. For `r`, `r(n^2)` is normal by Theorem 1.4 of Mullner with `q = 2` and `m' = 2`, where `gcd(q - 1, m') = 1` and the values of `s_11` have greatest common divisor `1`, and a normal sequence has all `2^k` factors. A reading equal to the squares stream times a periodic sign is automatic in a base exactly when the squares stream is, since periodic sequences are automatic in every base and termwise products keep automaticity (the companion paper, Lemma 2.4); by Proposition 6.1 this covers `r` along `n^2 + 1` and `2n^2 + 1`. For `w(n) = r(n^2 + n) = r((2n + 1)^2)`, the squares stream `y(n) = r(n^2)` satisfies `y(2m) = y(m)` and `y(2m + 1) = w(m)`, so a factor of `y` of length `2k` is fixed by the parity of its start, a factor of `y` of length `k` and a factor of `w` of length `k`. Hence `4^k <= 2 . 2^k p_w(k)`, with `p_w(k)` the number of factors of `w` of length `k`, so `p_w(k) >= 2^(k-1)`, more than Proposition 7.1 allows in any base. □

The same interleaving, from the normality of `t(n^2)` (Theorem 1 of the preprint of Drmota, Mauduit and Rivat), recovers the readings of `t` along these five polynomials without Shen's theorem. For Thue-Morse, Theorem 7.3 is in print: Allouche and Shallit 2003, Theorem 6.10.1, unread here, makes Thue-Morse along a polynomial non-automatic exactly when the degree is at least `2`, as Popoli 2020 quotes it in its Section 2, and Shen's theorem excludes every base. Popoli 2020, in its Section 1, also records `r(n^2)` as non-automatic, through Mullner's normality and Christol's theorem. What Theorem 7.3 adds is Rudin-Shapiro along `n^2 + n` and `(2n + 1)^2`. For the ternary parity Mullner's theorem does not apply, since its hypothesis `gcd(q - 1, m') = 1` fails at `q = 3`, `m' = 2`, and only the `3`-kernel witness of Theorem 7.2 is proved here.

## 8. The battery on the readings

**Fact 8.1 (the battery at `N = 2^15`).** At `N = 2^15` pairs per test, with the band `4/sqrt N = 0.0221`, the battery of Definition 2.3 reads the three designs along the eight readings and three of the probes of Definition 2.2, in five rows, and the control, as follows; the last column counts the identities of the head scan, plain and twisted.

| design | reading | mean | shift | dilation | kernel | square | identities |
|---|---|---|---|---|---|---|---|
| `t` | `n` | `+0.0000` | `0.3340` at 32 | `1.0000` at `(2, 0)` | `1.0000` | `+1.0000` | 42, 0 |
| `t` | `n^2` | `-0.0150` | `0.0280` at 15 | `1.0000` at `(2, 0)` | `1.0000` | `+1.0000` | 9, 0 |
| `t` | `n^2 + 1` | `+0.0112` | `0.0284` at 15 | `0.0135` at `(9, 5)` | `1.0000` | `+1.0000` | 6, 3 |
| `t` | `n^2 + n` | `+0.0145` | `0.0144` at 28 | `0.0293` at `(4, 3)` | `0.0045` | `+0.0014` | 0, 0 |
| `t` | `(2n + 1)^2` | `-0.0145` | `0.0144` at 28 | `0.0293` at `(4, 3)` | `0.0045` | `+0.0014` | 0, 0 |
| `t` | `2n^2 + 1` | `+0.0150` | `0.0280` at 15 | `1.0000` at `(2, 0)` | `1.0000` | `+1.0000` | 9, 0 |
| `t` | `n^3` | `-0.0008` | `0.0130` at 15 | `1.0000` at `(2, 0)` | `1.0000` | `+1.0000` | 9, 0 |
| `t` | `p_n` | `-0.0433` | `0.0568` at 2 | `0.0111` at `(6, 5)` | `0.0029` | `-0.0037` | 0, 0 |
| `r` | `n` | `+0.0078` | `0.0001` at 1 | `1.0000` at `(2, 0)` | `1.0000` | `+1.0000` | 35, 7 |
| `r` | `n^2`, `n^2 + 1` | `+0.0082` | `0.0104` at 35 | `1.0000` at `(2, 0)` | `1.0000` | `+1.0000` | 9, 0 |
| `r` | `n^2 + n`, `(2n + 1)^2` | `+0.0089` | `0.0153` at 8 | `0.0169` at `(4, 1)` | `0.0156` | `-0.0034` | 0, 0 |
| `r` | `2n^2 + 1` | `+0.0010` | `0.0105` at 35 | `0.0223` at `(3, 0)` | `1.0000` | `+1.0000` | 6, 3 |
| `r` | `n^3` | `+0.0142` | `0.0197` at 4 | `1.0000` at `(2, 0)` | `1.0000` | `+1.0000` | 9, 0 |
| `r` | `p_n` | `-0.0007` | `0.0115` at 8 | `0.0143` at `(4, 3)` | `0.0012` | `+0.0004` | 0, 0 |
| `t3` | `n` | `+0.0001` | `0.3342` at 54 | `1.0000` at `(3, 0)` | `1.0000` | `+1.0000` | 36, 0 |
| `t3` | `n^2` | `+0.0185` | `0.0129` at 61 | `1.0000` at `(3, 0)` | `1.0000` | `+1.0000` | 6, 0 |
| `t3` | `n^2 + 1` | `-0.0024` | `0.0115` at 18 | `0.3333` at `(3, 0)` | `1.0000` | `+1.0000` | 4, 2 |
| `t3` | `n^2 + n` | `-0.0047` | `0.0128` at 2 | `1.0000` at `(3, 1)` | `0.0108` | `+0.0051` | 6, 0 |
| `t3` | `(2n + 1)^2` | `+0.0159` | `0.0131` at 8 | `1.0000` at `(3, 1)` | `0.0024` | `+0.0032` | 6, 0 |
| `t3` | `2n^2 + 1` | `+0.0024` | `0.0172` at 54 | `0.3334` at `(3, 0)` | `1.0000` | `+0.1112` | 4, 0 |
| `t3` | `n^3` | `+0.0060` | `0.0165` at 41 | `1.0000` at `(3, 0)` | `1.0000` | `+1.0000` | 6, 0 |
| `t3` | `p_n` | `+0.0027` | `0.0694` at 1 | `0.0162` at `(7, 2)` | `0.0131` | `-0.0035` | 0, 0 |
| `t` | `n^2 + 4` | `+0.0108` | `0.0254` at 60 | `0.5002` at `(2, 0)` | `0.0036` | `+1.0000` | 5, 4 |
| `t` | `(n + 1)^2 + 1` | `+0.0112` | `0.0284` at 15 | `0.0171` at `(5, 0)` | `0.0098` | `+0.0014` | 6, 3 |
| `t` | `n^3 + n` | `+0.0170` | `0.0138` at 33 | `0.0524` at `(4, 0)` | `0.0658` | `+0.0038` | 0, 0 |
| `r` | `n^3 + n` | `+0.0016` | `0.0124` at 25 | `0.0171` at `(4, 0)` | `0.0231` | `+0.0253` | 0, 0 |
| `t3` | `n^3 + n` | `+0.0061` | `0.0130` at 64 | `0.0238` at `(9, 0)` | `0.0415` | `+0.0178` | 0, 0 |
| `splitmix64` | seed `1` | `+0.0122` | `0.0114` at 4 | `0.0093` at `(9, 4)` | `0.0061` | `+0.0010` | 0 |

Every exact near-head reads `1.0000` where Theorem 4.2 puts it, and the local constants of Section 5 reappear as the dilation maxima `0.5002`, `0.3333` and `0.3334`. A twisted pair of mean `0` hides from its own dilation test, Thue-Morse along `n^2 + 1` reading `0.0135`, and the kernel test catches the reading at depth `1`. The free readings split: Thue-Morse fails the dilation test along `n^2 + n`, `0.0293` at `(4, 3)`, and the kernel test along `n^3 + n`, `0.0658`, while Rudin-Shapiro along `n^2 + n` reads at most `0.0169` and along the primes at most `0.0143`. The primes fail the shift test for Thue-Morse, `0.0568` at `h = 2`, and for `t3`, `0.0694` at `h = 1`. Every reading but `n` has the control's LZ78 rate, `1.208` to `1.216` bits per bit at `2^18` against `1.215`: LZ78 does not separate them.

**Fact 8.2 (an identity at a depth the battery does not read).** In the run of Fact 8.1, Thue-Morse along `(n + 1)^2 + 1` reads dilation `0.0171`, kernel `0.0098` and square `+0.0014`, all inside the band, while `y(2n + 1) = -t((n + 1)^2) = y(4n + 3)` at every `n`. Its centre `-1` is visited at depth `0`, where `(2, 1)`, `(4, 3)` and `(8, 7)` are near-heads (Corollary 4.3) whose limit `local` prints as `0`, but its constant `1` puts the identity at depth `1` of the contraction `2n + 1` (Theorem 4.4 and Remark 4.5), which no battery pair reads; the head scan finds it with five more. The reading does fail the shift test, `0.0284` at `h = 15`, as the squares do with `0.0280`, since `y(n) = (-1)^n t((n + 1)^2)` by Proposition 6.1.

The battery reads depth `0` at the centres of Corollary 4.3 and depth `1` at the centre `0` alone: an identity of Theorem 4.4 at another centre or at a greater depth escapes its two-point tests, and a centre `-p/d` needs a contraction whose ratio `b^i` has `d` dividing `b^i - 1`.

**Remark 8.3 (a free reading can fail).** Freedom does not imply passing: a polynomial reading with no near-head among the battery pairs need not pass the battery. At `N = 2^20`, `2^22`, `2^24` and `2^26`, Thue-Morse along `n^2 + n` fails the dilation test at `(4, 0)` or `(4, 3)` with `z = 9.48`, `10.51`, `11.86` and `14.42`, against the control's `3.83`, `3.34`, `2.33` and `2.81`; Thue-Morse along `n^3 + n` fails its kernel test with `z = 29.60`, `36.75`, `42.69` and `65.52`; and Rudin-Shapiro along `n^3 + n`, though its design is white, fails its kernel test with `z = 9.41`, `16.75`, `17.43` and `17.10`. Generator: the verb `falsify 26`.

## 9. The mass at one third

The free readings are not all alike, and the difference is in the design.

**Theorem 9.1 (the mass at one third).** For `H >= 1` and `r` in `{0, 1, 2}` let `U_r(H)` be the sum of `C_t(h)` over `0 <= h < H` with `h = r mod 3`, and `U(H) = (U_0(H), U_1(H), U_2(H))`. Then, with `[H = j]` the indicator of `H = j mod 3`,

```
U(2H) = M U(H) + f(H),    M = [[1, -1/2, -1/2], [-1/2, -1/2, 1], [-1/2, 1, -1/2]],
                          f(H) = (-[H = 2] C_t(H)/2, -[H = 1] C_t(H)/2, 1/2 - [H = 0] C_t(H)/2),
```

and for every `k >= 1`

```
U_0(2^k) = (16/45)(3/2)^k + 8/15   (k even),    U_0(2^k) = (16/45)(3/2)^k + 7/15   (k odd).
```

So the mean of `C_t(h)` over the multiples of `3` below `2^k` is `(16/15)(2^k)^(-log2(4/3)) + O(2^(-k))`, `log2(4/3) = 0.415...`, against `2/(3 . 2^k)` over all lags below `2^k`.

*Proof.* Split the lags `h < 2H` by parity and apply (2.3), indices of `U` modulo `3`. The even lags `h = 2q = r mod 3` are the `q < H` with `q = 2r mod 3`, and contribute `U_(2r)(H)`. The odd lags `h = 2q + 1 = r mod 3` are the `q < H` with `q = 2r - 2 mod 3`, and contribute `-1/2` times the sum of `C_t(q) + C_t(q + 1)` over them, which is `U_(2r-2)(H) + U_(2r-1)(H) - [2r - 1 = 0] C_t(0) + [H = 2r - 1] C_t(H)`. With `C_t(0) = 1`, the three cases `r = 0, 1, 2` are the three rows of `U(2H) = M U(H) + f(H)`. The columns of `M` sum to `0` and the entries of `f(H)` to `1/2 - C_t(H)/2`, so the sum `S(H)` of all `C_t(h)`, `h < H`, obeys `S(2H) = 1/2 - C_t(H)/2`; as `C_t(2^j) = C_t(1) = -1/3` for every `j >= 0`, `S(2^k) = 2/3` for `k >= 1`.

`M` is symmetric with eigenvalues `0`, `3/2` and `-3/2` on the orthogonal eigenvectors `(1, 1, 1)`, `(2, -1, -1)` and `(0, 1, -1)`. At `H = 2^k`, `C_t(H) = -1/3`, and `H = 1 mod 3` for even `k` and `2 mod 3` for odd `k`, so `f(2^k) = (0, 1/6, 1/2)` at even `k` and `(1/6, 0, 1/2)` at odd `k`. Write `U(2^k) = a_k (2, -1, -1) + c_k (0, 1, -1) + e_k (1, 1, 1)`. Then `e_k = S(2^k)/3 = 2/9` for `k >= 1`, the vector `(0, 1, -1)` has first coordinate `0`, and `U_0(2^k) = 2 a_k + 2/9`. Projecting the recursion on `(2, -1, -1)` gives `a_(k+1) = (3/2) a_k + f(2^k) . (2, -1, -1)/6`, that is

```
a_(k+1) = (3/2) a_k - 1/9   (k even),    a_(k+1) = (3/2) a_k - 1/36   (k odd),    a_0 = 1/3,
```

`a_0` from `U(1) = (1, 0, 0)`. The formulas `a_k = (8/45)(3/2)^k + 7/45` at even `k` and `a_k = (8/45)(3/2)^k + 11/90` at odd `k` give `a_0 = 1/3` and satisfy both steps, since `(3/2)(7/45) - 1/9 = 11/90` and `(3/2)(11/90) - 1/36 = 7/45`, so they hold at every `k` by induction, and `U_0(2^k) = 2 a_k + 2/9` is the closed form. There are `2^k/3 + O(1)` multiples of `3` below `2^k`, which gives the mean. □

The study `readings`, verb `leak`, checks the recursion against direct sums of `C_t` within `2.27e-12` at every `k <= 24`, prints `U_0(2^k) (2/3)^k = 0.355556` at `k = 40`, `50`, `60` and `61`, and checks the closed form in exact integers, as `45 V_k[0] = 96 . 3^k + 144 . 2^k` at even `k` and `96 . 3^k + 126 . 2^k` at odd `k` for the integer vector `V_k = 6 . 2^k U(2^k)`, at every `k = 1, ..., 70`. The growth `(3/2)^k` sits in the residue classes modulo `3`, not in the total, which stays `2/3`: it is a growth at frequency one third, since the real part of `sum_(h < 2^k) C_t(h) omega^h`, `omega` a primitive cube root of unity, is `U_0 - (U_1 + U_2)/2 = 3 a_k`. In spectral language, the `C_t(h)` are the Fourier coefficients of the spectral measure of Thue-Morse, the Riesz product of the factors `1 - cos(2 pi 2^j theta)` (Baake and Coons 2024, Remark 3.2), and at `theta = 1/3` every factor is `3/2`.

## 10. What is in print

**Normality along polynomials.** Drmota, Mauduit and Rivat 2019 is read in the authors' 2013 preprint, whose numbering is used throughout; the journal text is unread. Theorem 1 there makes `t(n^2)` normal, and Theorem 2 bounds each sum `sum_(n < N) e((1/2) sum_l alpha_l s_2((n + l)^2))`, `alpha` nonzero in `{0, 1}^k`, by `N^(1 - eta)` with `eta > 0` never given a value; the power saving belongs to Theorem 2 and Lemma 1, not to Theorem 1. Their Conjecture 1 on polynomials is read only through quotations that disagree on its degree, and is not used. Mullner 2018, read at source, makes `b(n^2) mod m'` normal for a strongly block-additive `b` in base `q` under `gcd(q - 1, m') = 1` and `gcd(m', gcd{b(n)}) = 1` (Theorem 1.4), which covers Rudin-Shapiro at `q = 2`, `m' = 2` and not the ternary parity. Mauduit and Rivat 2018, on Rudin-Shapiro along the squares, is unread and not used. Before normality, every block was known to occur in `t(n^2)` (Moshe 2007, cited unread, as the preprint of Drmota, Mauduit and Rivat and the Mullner 2018 talk state it), and Shen 2022, read at source, proves that `t(H(n))` has all `2^m` blocks of every length `m` for every `H` in `Q[T]` of degree at least `2` with `H(N)` in `N` (Theorem 1). With Proposition 7.1 this puts every polynomial reading of Thue-Morse of degree at least `2` outside every automatic class, so the Thue-Morse half of Theorem 7.3 is a case of Shen's theorem; that these readings are not automatic is older still, Allouche and Shallit 2003, Theorem 6.10.1, unread here and quoted by Popoli 2020 in its Section 2 and, for the squares, by Sun and Winterhof 2019 in their Section 1. Popoli 2020, in its Section 1, also records `r(n^2)` as non-automatic, through Mullner's normality and Christol's theorem. Along the cubes, Spiegelhofer 2023, read at source, proves that `t(n^3)` takes each value with density `1/2`, with a power saving (Theorem 2.1), and the normality of `t(n^3)` and of `t(n^3 + n)` is open; Merai and Winterhof 2022 record that it is not known whether `t(f(p))`, along the primes, is normal for any nonconstant polynomial `f`, in their Section 7. Drmota, Mauduit and Rivat 2011, read in the authors' submitted copy, treat the digit sums of polynomial values in prime bases `q >= q_0(d)`, small bases staying open there (Remark 2).

**Pseudorandomness and the tests used.** Sun and Winterhof 2019, read at source, bound the maximum order complexity of Thue-Morse along the squares from below by `sqrt(2N/5)` (Theorem 1) and offer the squares as a cryptographic candidate; Popoli 2020, read at source, proves the lower bound `>> N^(1/d)` along monic polynomials of degree `d` (Theorems 3 and 4) and calls the squares a potentially better candidate. Merai and Winterhof 2022, read at source in Sections 5, 7 and 8, survey these sequences through the well-distribution and correlation measures of Mauduit and Sarkozy 1998, cited unread through them and Konieczny 2019: the correlation measure of order `k` is the largest `abs(sum_n s_(n + d_1) ... s_(n + d_k))`, a statistic of translations. They leave open whether the measures of order `k` of `t(n^2)` and `r(n^2)` are `o(N)` (Problem 3), and read their plots as suggesting square-root size (Figures 12 and 13, a numerical remark). No test that compares a reading at `n` with the same reading at `Bn + r` was found in these sources or in the others read, and by Theorem 3.1 the readings they propose, along the squares, keep an exact dilation identity.

**Dilations in the literature.** Mauduit and Rivat 2015, read at source, use the dilation `b(2^kappa n) = b(n)` to read off their Fourier property for Rudin-Shapiro-like sequences along the primes; here the same identity is the obstruction. Drmota, Mauduit, Rivat and Spiegelhofer 2020, read at source, prove `sum mu(n) g(n^2) = o(N)` for every strongly `b`-multiplicative `g` of modulus `1` (Theorems 1 and 2), and their Section 1.4 needs `sum_n g(p^2 n^2) conj(g(q^2 n^2)) = o(N)` for distinct primes `p`, `q` prime to `b`: a correlation of the reading `g(n^2)` at the pair `(pn, qn)`, which is not a near-head in the sense of Definition 2.4, since `(q/p)^2` is not a power of `b`. It is the one family of free dilation pairs of a polynomial reading proved to decorrelate in the sources read; the pairs whose ratio is a power of `2`, where the leak of Conjecture 11.1 lives, are not among them.

**Densities along the squares.** Adamczewski, Drmota and Mullner 2022, read at source, prove that the logarithmic densities of every automatic sequence along the squares exist and are computable (Theorems 1.1 to 1.4), and give the densities of primitive prolongable sequences along the squares (Theorem 8.1), letter densities only. A pair `x(n^2 + c) x(n^2)` is a letter of the automatic sequence `m -> x(m + c) x(m)` read along the squares, so the existence of the pair law of Section 5 in the logarithmic sense is a case of their theorems, and their Theorems 1.3 and 9.5 make it rational and computable at a prime base; Theorem 5.1 gives the natural density, equal to the logarithmic one, along every nonconstant integer polynomial, with an enclosure at every `J`. Drmota and Morgenbesser 2012, read in the authors' preprint, give the distribution of `T(n^2)` for generalized Thue-Morse sequences with values in a compact group, one letter at a time (Theorem 1). Deshouillers, Drmota and Morgenbesser 2012, read at source, prove that along `floor(n^c)`, `1 < c < 7/5`, the pairs `(t(floor(n^c)), t(floor(n^c) + k))` keep the Thue-Morse pair law, `-1/3` at `k = 1`, and Mullner and Spiegelhofer 2017, read in its abstract, that `t(floor(n^c))` is normal for `1 < c < 3/2`; along the squares the pair law of Section 5 differs from the Thue-Morse law at each of the eight offsets computed.

**The correlations of Thue-Morse.** Baake and Coons 2024, read at source in Sections 1 to 3, give the recursion (2.3) in matrix form, credit it to Mahler, Kakutani and Baake and Grimm, and write the spectral measure as the Riesz product `prod (1 - cos(2^(l + 1) pi x))` (Remark 3.2). Coons, Mazac, Pincus-Kazmar and Stout 2025, read at source, record `C_t(2^n) = -1/3` and `max abs C_t(m) = 1/3` and bound `sum_(m <= x) abs C_t(m)` between constant multiples of `x^0.6274882485` and `x^0.6464616661` (Theorem 1.1); Zaks 2001, read at source, obtains a linear recursion for the sums of `C_t(h)^2` over dyadic blocks. None of them sums `C_t` over a residue class, and the closed form of Theorem 9.1 is not found in the sources read; Queffelec 1987 and Baake and Grimm 2013, cited unread, were not searched and may contain it.

**The primes.** Aloui, Mauduit and Mkaouar 2021, cited unread through Popoli 2023, treat the blocks `t(p), t(p + 1), ..., t(p + k)` at one prime, not the consecutive primes of Conjecture 11.4. Popoli 2023, read at source in Sections 1 and 2, records in its Section 1.4 that the normality of Thue-Morse along the primes is open, and proves that `r(p) r(p + 1) ... r(p + k)` has mean `0` along the primes (Theorem 2.1).

## 11. Open problems

Three conjectures, one problem and three questions are left.

**Conjecture 11.1 (the leak at one third).** A free battery pair `(phi, psi)` of a polynomial reading of Thue-Morse whose remainder `P(psi(n)) - 2^L P(phi(n))`, `2^L` the ratio of the leading coefficients, is linear in `n` with slope divisible by `3` has a correlation of order `N^(-log2(4/3)) = N^(-0.415...)`; a reading with such a pair then fails the battery at every large `N`.

The verb `falsify 26` runs the free readings at `N = 2^20`, `2^22`, `2^24` and `2^26`. For Thue-Morse along `n^2 + n` the largest Pearson dilation correlation is `0.009261`, `0.005130`, `0.002895` and `0.001760`, always at `(4, 0)` or `(4, 3)` with `(8, 5)` third, three pairs whose remainders are `-12n`, `12n + 12` and `24n + 30`; the fourth battery pair at a power of `2` whose remainder has slope divisible by `3`, `(8, 2)` with remainder `-24n + 6`, is not among these peaks; the slopes in `log2 N` are `-0.426`, `-0.413` and `-0.359` step by step and `-0.399` overall, where the noise has `-0.5`, and `z` is `9.48`, `10.51`, `11.86` and `14.42` against the control's `3.83`, `3.34`, `2.33` and `2.81`. Thue-Morse along `n^3 + n`, whose kernel pair leaves `P(4n) - 8 P(2n) = -12n`, fails harder, with kernel `z = 29.60`, `36.75`, `42.69` and `65.52` (Remark 8.3). And `t(n^2)` itself, whose shift pairs leave `(n + h)^2 - n^2 = h (2n + h)`, has its largest shift correlation at `h = 15` at every size run, `2^15` and `2^20` to `2^26`, from `0.0280` down to `0.0011`, with slope `-0.384` from `2^20` to `2^26`, slower than the noise. The evidence is two free polynomials of Thue-Morse, and the step slopes move from `-0.426` to `-0.359`, so these sizes do not pin the rate down.

A heuristic links the failing pairs to Theorem 9.1. For `n^2 + n` at the dilation `(4, r)`, `P(4n + r) - 16 P(n) = 4 (2r - 3) n + r (r + 1)`, so `y(4n + r) = t(16 A + R(n))` with `A = n^2 + n`, `R` that remainder and `t(16 A) = t(A) = y(n)`. Writing `R(n) = 16 a + c` with `0 <= c < 16`, `t(16 A + R) t(16 A) = t(c) t(A + a) t(A)`. If the high digits of `A` were uniform and independent of `n` modulo `4`, the mean would be a weighted mean of `C_t(a)` with weights `t(c)`, where `a` is near `(2r - 3) n/4` and `c` repeats with period `4` in `n`. As functions of the lag the weights then repeat with period `abs(2r - 3)`, and when `3` divides `2r - 3`, at `r = 0` and `r = 3`, they have a component at frequency one third, where Theorem 9.1 puts the mass. The same holds at `(8, r)` with slope `2r - 7`, divisible by `3` at `r = 2` and `r = 5`. A proof needs the joint distribution of the high digits of `n^2 + n` with `n`: a carry lemma for two polynomials.

**Conjecture 11.2 (a fourth rung).** Every battery correlation of Rudin-Shapiro along `n^2 + n`, equivalently along `(2n + 1)^2`, is `O(N^(-1/2 + epsilon))` for every `epsilon > 0`, the rate of independent signs.

The reading leaves the automatic class in every base (Theorem 7.3), every battery pair is free for it (Corollary 4.3), and its design has no correlation at any nonzero lag to leak. At `N = 2^20`, `2^22`, `2^24` and `2^26` its largest `z` over the `64` shift, `44` dilation and `2` kernel correlations and the square is `2.14`, `2.56`, `2.68` and `3.41`, against the control's `3.83`, `3.34`, `2.76` and `2.81`, and Rudin-Shapiro along the primes reads at most `3.32`. The ladder of the companion paper climbs three rungs of generators built from designs, one base, two bases and rounds of a cellular automaton; this reading would be a fourth, a design read along a polynomial that leaves the automatic class and whose battery correlations decay as those of independent signs. Whiteness alone is not enough, since Rudin-Shapiro along `n^3 + n` fails (Remark 8.3).

**Problem 11.3 (joint normality).** Prove that every dilation correlation of Rudin-Shapiro along `n^2 + n` tends to `0`. More generally, decide for which pairs of integer polynomials `P_1`, `P_2` that are not related as in Definition 2.4 the mean of `x(P_1(n)) x(P_2(n))` tends to `0`.

The dilation correlation of a free pair `(n, Bn + r)` of a reading along `P` is the mean of `x(P(n)) x(P(Bn + r))`, a correlation of the design at two polynomial values; its decay needs the joint distribution of the base-`b` digits of `P(n)` and `P(Bn + r)`, a joint normality of the pair, which no theorem read here gives. The normality theorems control the design at the consecutive values `P(n + l)`, through exponential sums over the shifts (Theorem 2 of the preprint of Drmota, Mauduit and Rivat), and the only free dilation pairs of a polynomial reading proved to decorrelate in the sources read are the pairs `(pn, qn)` of the squares of Drmota, Mauduit, Rivat and Spiegelhofer 2020, `p` and `q` distinct primes prime to the base (Section 10). Theorem 5.1 settles the near-heads; the free pairs are open.

**Conjecture 11.4 (the primes).** Along the consecutive primes, Thue-Morse and the ternary parity fail the shift test at every large `N`.

At `N = 2^20`, `2^22`, `2^24` and `2^26`, Thue-Morse reads `0.0545`, `0.0522`, `0.0503` and `0.0491` at `h = 2`, and `t3` reads `0.0552`, `0.0524`, `0.0492` and `0.0467` at `h = 1`, with `z` up to `402`, decaying more slowly than any power of `N` seen. A heuristic: consecutive primes past `2` are odd, and `t(2a + 1 + g) t(2a + 1) = t(a + g/2) t(a)` for even `g`, so a gap `g` contributes the Thue-Morse correlation at `g/2`; the gaps of consecutive primes are, heuristically, biased toward multiples of `6`, and a weighted mean of `C_t` over lags biased modulo `3` feels the mass of Theorem 9.1, spread over a range of order `log p_N` and so decaying in `log N`, not in `N`. For `t3`, a gap `g` divisible by `6` keeps the last ternary digit and moves the digits above by `g/3`, and a gap of `6` contributes `C_(t3)(2) = -1/3`. Neither heuristic is a proof.

**Question 11.5 (Rudin-Shapiro along the cubic).** Which statistic carries the failure of Rudin-Shapiro along `n^3 + n` (Remark 8.3)? Its kernel pair leaves `P(4n) - 8 P(2n) = -12n`, a slope divisible by `3`, but Rudin-Shapiro has no mass at one third; a correlation of `r` restricted to a class of the low digits is a candidate.

**Question 11.6 (the pair law along the squares).** Compute the starred cells of the pair law exactly. They are rational: the natural densities of Theorem 5.1 equal the logarithmic ones, which Adamczewski, Drmota and Mullner 2022 prove rational and computable along the squares at a prime base (Theorems 1.3 and 9.5), and `local` puts them within `1.40e-7` of `1/6`, `-2/3`, `-1/3`, `1/3`, `-5/9` and `-1/9`. Are they these rationals, and is there a closed form as an integral over the `b`-adic distribution of the squares?

**Question 11.7 (the ternary squares).** Is `t3(n^2)` normal? Mullner's theorem does not apply, since `gcd(3 - 1, 2) = 2`; the `3`-kernel is infinite (Theorem 7.2); and the battery cannot tell, since the identity `t3(9n^2) = t3(n^2)` of Theorem 3.1 catches the reading first.

## Reproducibility

Every number of this paper that the text does not derive by hand is printed by one single-threaded Rust study without dependencies, run from the repository root as `cargo run --release --manifest-path research/random/lab/rs/readings/Cargo.toml -- <verb>`:

- `identities`: Proposition 6.1, Theorem 7.2, the twisted pair of Thue-Morse along `n^2 + 4` and the identity of Fact 8.2, at every `n < 2^24`; under `1` s.
- `battery`: the readings of Section 3, Facts 8.1 and 8.2 and the head scan; about `13` s.
- `local`: Theorem 5.1 and Corollary 5.2, the census of the `143` near-head pairs at `J = 20` (`13` at base `3`) and the pair law at `J = 26` (`16` at base `3`); about `25` s.
- `falsify 26`: Remark 8.3 and the evidence of Conjectures 11.1, 11.2 and 11.4; about `130` s, of which `23` s is a segmented sieve of `603979904` primes.
- `leak`: Theorem 9.1, its recursion against direct sums and its closed form in exact integers; under `1` s.

The figure `paper-polynomial-readings` recomputes `C_t(h)` for `h <= 2^12` by the recursion (2.3) in floating point, and refuses to draw if `U_0(2^k)` misses the closed form of Theorem 9.1, or the sum of `C_t(h)` over `h < 2^k` misses `2/3`, by more than `1e-9` at any `1 <= k <= 12`; the exact check is the verb `leak`. Every Theorem, Proposition, Lemma and Corollary proved here and every Fact is a dated Proved or Verified row of [the readings](../claims/readings.md) claim page; the Conjectures are rows of the same page, and Remark 8.3 is its Refuted row.

## References

- Adamczewski, Drmota and Mullner 2022, (Logarithmic) densities for automatic sequences along primes and squares, Trans. Amer. Math. Soc. 375, 455-499, arXiv:2009.14773; read at source in the arXiv v2, Theorems 1.1 to 1.4 and 8.1. [arxiv.org/abs/2009.14773](https://arxiv.org/abs/2009.14773)
- Allouche and Shallit 2003, Automatic Sequences: Theory, Applications, Generalizations, Cambridge University Press; not read, cited for the kernel criterion and, as Popoli 2020 and Sun and Winterhof 2019 quote them, for Theorem 6.10.1 and Corollary 10.3.2. [doi.org/10.1017/CBO9780511546563](https://doi.org/10.1017/CBO9780511546563)
- Aloui, Mauduit and Mkaouar 2021, On the correlation of the sum of digits along prime numbers, Acta Arith. 199(3), 275-301; unread, its theorem read as Popoli 2023 states it. [doi.org/10.4064/aa200514-22-12](https://doi.org/10.4064/aa200514-22-12)
- Baake and Coons 2024, Correlations of the Thue-Morse sequence, Indag. Math. 35, 914-930, arXiv:2209.07102; read at source in the arXiv v2, Sections 1 to 3. [arxiv.org/abs/2209.07102](https://arxiv.org/abs/2209.07102)
- Baake and Grimm 2013, Aperiodic Order, Vol. 1: A Mathematical Invitation, Cambridge University Press; cited unread, as an account of the Thue-Morse correlations.
- Coons, Mazac, Pincus-Kazmar and Stout 2025, On the absolute value of the autocorrelations of the Thue-Morse sequence, arXiv:2511.06386; read at source in the arXiv v1, Section 1, Theorem 1.1 and Section 3. [arxiv.org/abs/2511.06386](https://arxiv.org/abs/2511.06386)
- Deshouillers, Drmota and Morgenbesser 2012, Subsequences of automatic sequences indexed by floor(n^c) and correlations, J. Number Theory 132(9), 1837-1866; read at source in the authors' copy, the abstract and Theorems 1 to 3. [dmg.tuwien.ac.at](https://www.dmg.tuwien.ac.at/drmota/ddm.pdf)
- Drmota and Morgenbesser 2012, Generalized Thue-Morse sequences of squares, Israel J. Math. 190, 157-193; read at source in the authors' 2010 preprint, Theorem 1 and its remarks; the published text is unread. [dmg.tuwien.ac.at](https://www.dmg.tuwien.ac.at/drmota/generalized_thue_morse.pdf)
- Drmota, Mauduit and Rivat 2011, The sum-of-digits function of polynomial sequences, J. London Math. Soc. 84(1), 81-102; read at source in the authors' submitted copy, Theorem 1 and Remarks 1 and 2. [dmg.tuwien.ac.at](https://www.dmg.tuwien.ac.at/drmota/dmr6.pdf)
- Drmota, Mauduit and Rivat 2019, Normality along squares, J. Eur. Math. Soc. 21(2), 507-548, doi:10.4171/JEMS/843. The text read is the authors' 2013 preprint, The Thue-Morse sequence along squares is normal, and its numbering is used here: Theorem 1, Theorem 2 and Lemma 1. The journal text is unread, and its Conjecture 1 is read only through quotations. [dmg.tuwien.ac.at](https://www.dmg.tuwien.ac.at/drmota/alongsquares.pdf)
- Drmota, Mauduit, Rivat and Spiegelhofer 2020, Moebius orthogonality of sequences with maximal entropy, arXiv:2001.07898; read at source in the arXiv v1, the abstract, Theorems 1 and 2 and Section 1.4. [arxiv.org/abs/2001.07898](https://arxiv.org/abs/2001.07898)
- Kakutani 1972, Strictly ergodic symbolic dynamical systems, Proc. Sixth Berkeley Symp. Math. Statist. Probab. 2, 319-326; cited unread, for the Thue-Morse correlations as Baake and Coons 2024 credit them.
- Konieczny 2019, Gowers norms for the Thue-Morse and Rudin-Shapiro sequences, Ann. Inst. Fourier 69(4), 1897-1913, arXiv:1611.09985; read at source in the arXiv v2, here for its account of Mauduit and Sarkozy 1998. [arxiv.org/abs/1611.09985](https://arxiv.org/abs/1611.09985)
- Mahler 1927, The spectrum of an array and its application to the study of the translation properties of a simple class of arithmetical functions, Part Two, J. Math. Phys. 6, 158-163; cited unread, for the Thue-Morse correlations as Baake and Coons 2024 credit them. [doi.org/10.1002/sapm192761158](https://doi.org/10.1002/sapm192761158)
- Mauduit and Rivat 2015, Prime numbers along Rudin-Shapiro sequences, J. Eur. Math. Soc. 17, 2595-2642; read at source in the publisher file, page 2597, Definitions 1 and 2 and Theorem 1. [doi.org/10.4171/JEMS/566](https://doi.org/10.4171/JEMS/566)
- Mauduit and Rivat 2018, Rudin-Shapiro sequences along squares, Trans. Amer. Math. Soc. 370(11), 7899-7921; unread and not used. [doi.org/10.1090/tran/7210](https://doi.org/10.1090/tran/7210)
- Mauduit and Sarkozy 1998, On finite pseudorandom binary sequences. II. The Champernowne, Rudin-Shapiro, and Thue-Morse sequences, a further construction, J. Number Theory 73(2), 256-276; unread, its measures read as Merai and Winterhof 2022 and Konieczny 2019 state them.
- Merai and Winterhof 2022, Pseudorandom sequences derived from automatic sequences, Cryptogr. Commun. 14(4), 783-815, arXiv:2105.03086; read at source in the arXiv v1, Sections 5, 7 and 8. [arxiv.org/abs/2105.03086](https://arxiv.org/abs/2105.03086)
- Moshe 2007, On the subword complexity of Thue-Morse polynomial extractions, Theoret. Comput. Sci. 389(1-2), 318-329; unread, read as the preprint of Drmota, Mauduit and Rivat and the Mullner 2018 talk state it.
- Mullner 2018, The Rudin-Shapiro sequence and similar sequences are normal along squares, Canad. J. Math. 70(5), 1096-1129, arXiv:1704.06472; read at source in the arXiv v2, Definition 1.1 and Theorem 1.4. [arxiv.org/abs/1704.06472](https://arxiv.org/abs/1704.06472)
- Mullner 2018 talk, The Rudin-Shapiro sequence and similar sequences are normal along squares, talk slides; read at source. [dmg.tuwien.ac.at](https://dmg.tuwien.ac.at/cmuellner/talks/Marseille_10_2018.pdf)
- Mullner and Spiegelhofer 2017, Normality of the Thue-Morse sequence along Piatetski-Shapiro sequences, II, Israel J. Math. 220, 691-738, arXiv:1511.01671; read in its abstract. [arxiv.org/abs/1511.01671](https://arxiv.org/abs/1511.01671)
- Popoli 2020, On the maximum order complexity of Thue-Morse and Rudin-Shapiro sequences along polynomial values, arXiv:2011.03457; read at source in the arXiv v1, Sections 1 and 2 and Theorems 3 and 4. [arxiv.org/abs/2011.03457](https://arxiv.org/abs/2011.03457)
- Popoli 2023, Correlation of the Rudin-Shapiro sequence along prime numbers, arXiv:2306.10629; read at source in Sections 1 and 2. [arxiv.org/abs/2306.10629](https://arxiv.org/abs/2306.10629)
- Queffelec 1987, Substitution Dynamical Systems - Spectral Analysis, Lecture Notes in Math. 1294, Springer; cited unread, as an account of the Thue-Morse spectrum that was not searched.
- Rudin 1959, Some theorems on Fourier coefficients, Proc. Amer. Math. Soc. 10, 855-859; cited unread, for the Rudin-Shapiro sequence. [doi.org/10.1090/S0002-9939-1959-0116184-5](https://doi.org/10.1090/S0002-9939-1959-0116184-5)
- Shapiro 1951, Extremal problems for polynomials and power series, Master's thesis, MIT; cited unread, for the Rudin-Shapiro sequence.
- Shen 2022, The subword complexity of polynomial subsequences of the Thue-Morse sequence, C. R. Math. Acad. Sci. Paris 360, 503-511; read at source in the open copy, Theorem 1 and its proof. [doi.org/10.5802/crmath.321](https://doi.org/10.5802/crmath.321)
- Spiegelhofer 2023, Thue-Morse along the sequence of cubes, arXiv:2308.09498; read at source in the arXiv v2, the abstract, Section 1 and Theorem 2.1. [arxiv.org/abs/2308.09498](https://arxiv.org/abs/2308.09498)
- Steele, Lea and Flood 2014, Fast splittable pseudorandom number generators, OOPSLA 2014; cited unread, `splitmix64` being the control. [doi.org/10.1145/2660193.2660195](https://doi.org/10.1145/2660193.2660195)
- Sun and Winterhof 2019, On the maximum order complexity of subsequences of the Thue-Morse and Rudin-Shapiro sequence along squares, Int. J. Comput. Math. Comput. Syst. Theory 4(1), 30-36, arXiv:1910.13763; read at source in the arXiv v1, the abstract, Section 1, Theorems 1 and 2 and Problems 1 to 6. [arxiv.org/abs/1910.13763](https://arxiv.org/abs/1910.13763)
- Zaks 2001, Fractal Fourier spectra in dynamical systems, Habilitationsschrift, University of Potsdam; read at source in the open copy, Chapter 2, Section 2.2. [publishup.uni-potsdam.de](https://publishup.uni-potsdam.de/files/145/zaks.pdf)
