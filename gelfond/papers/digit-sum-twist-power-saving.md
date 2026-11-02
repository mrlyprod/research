---
title: A Power Saving for the Digit Sum on Missing Digits
lead: On the integers whose base-`b` digits avoid a fixed set, keeping two consecutive digits and the one-step certificate (W), the Mobius and von Mangoldt functions modulated by the digit sum, `e(alpha s_b(n))`, save a power of `x` against the count of the set whenever `(b - 1) alpha` is not an integer, with an explicit exponent, an effective constant and no L-function; (W) holds at every one-missing-digit set from base `92317`, and the primes of the set then equidistribute by digit sum modulo every `q` coprime to `b - 1`. The power is tiny, `delta_0 >= 4.678 * 10^-9` at `92317` for every `alpha` with `norm((b - 1) alpha) >= 0.2`: the point is the mechanism.
date: 2026-10-05
figure: paper-digit-sum-twist-power-saving
---

Write the integers in a large base `b` and keep only those whose digits avoid a fixed set: they form a thin set `S_F`, `F` the `fill` digits kept, with `A_F(x) = x^(log(fill)/log(b) + o(1))` members up to `x` and no multiplicative structure. How the sum of the digits of a prime is distributed in residue classes is a question of Gelfond, answered for all primes by Mauduit and Rivat (2010): equidistributed modulo every `q` coprime to `b - 1`. This paper proves the same for the primes of `S_F` at large base. The main theorem bounds the Mobius and von Mangoldt sums over `S_F` modulated by the digit sum `s_b(n)`: if `F` keeps two consecutive digits and satisfies the one-step certificate (W) of [An Unconditional Mertens Bound at Large Base](../../coprime/papers/unconditional-mertens-at-large-base.md), and `(b - 1) alpha` is not an integer, then `sum_(n in S_F, n <= x) mu(n) e(alpha s_b(n))` and the same sum with `Lambda` are at most `C A_F(x) x^(-delta)`, with `delta` explicit and `C` effective, both uniform in `alpha` once `norm((b - 1) alpha)` is bounded below; (W) holds at every one-missing-digit set from base `92317`. Without the modulation that paper saves only `exp(-c sqrt(log x))`: near the fractions of small denominator it pays a contraction of the transform and, where the denominator divides a power of the base, `mu` in progressions through a zero-free region, both at root-log strength. The modulation removes that ceiling. A two-position lemma shows that on every pair of adjacent digits the modulated digit transform loses a fixed share of its mass at every frequency, exactly when `(b - 1) alpha` is not an integer, so every frequency near a fraction of small denominator is paid by the supremum of the transform times the number of such frequencies, and no L-function enters. The `l^1` bounds are imported unchanged, once a peel shows that the one-step certificate does not see the modulation, and the minor-arc bounds are re-read with the ratio `r = x/y` in place of `x^(o(1))`. The power is tiny, `delta_0 >= 4.678 * 10^-9` at base `92317` for every `alpha` with `norm((b - 1) alpha) >= 0.2`, and the size is not the point; the mechanism is, with its corollary for the primes. The mechanism is in print for the full digit set, in the survey of Maynard (2019), and the digit-sum classes of the integers of a missing-digit set are in Aloui (2015); what is written here is the mechanism carried onto the primes and `mu` of a missing-digit set.

## Introduction

Write every integer in base `b` and throw away those that use a forbidden digit. What is left is sparse, with no multiplicative structure, yet its primes still behave: Maynard (2019) proved that there are infinitely many of them in base `10`, and Maynard (2022) counted them at large base. Add up the digits of each such prime and sort the primes by that sum modulo some `q`. This paper proves that at large base every class gets its fair share, for every `q` coprime to `b - 1`. The tool is a sum that weighs each member `n` of the set by `mu(n)` or `Lambda(n)` and by the unit complex number `e(alpha s_b(n))`, the modulation, with `e(t) = exp(2 pi i t)` and `s_b(n)` the digit sum, which turns once around the circle each time the digit sum grows by `1/alpha`. That sum cancels by a power of `x` unless `(b - 1) alpha` is an integer, and the reason fits in one picture.

![The torus of phase pairs (v_1, v_2), each phase read in \[-1/2, 1/2) with v_1 across and v_2 up, painted at base 10 missing the digit 7 by the share |g_F(v_1) g_F(v_2)|/81 of the full mass on a ramp from the ground through blue to yellow: one yellow peak at the centre, where both phases are integers, a cross of fading blue lobes along the two axes, where one phase is an integer, and a faint plaid elsewhere. Eleven thin dim segments of slope 10 cross the square, the pairs (u, 10 u) of the lattice frequency alpha = 0, one of them through the middle of the peak. Ten thicker orange segments of the same slope sit halfway between them, the pairs (u + 1/2, 10 u + 1/2) of the frequency alpha = 1/2: the two nearest the centre cross the horizontal axis at v_1 = -1/20 and 1/20 on the flanks of the peak and meet the vertical axis only at the top and bottom edges, and no point of any orange segment has both phases within 1/22 of the integers.](paper-digit-sum-twist-power-saving)

The figure is the two-position lemma at base `10` missing the digit `7`. A frequency `t` meets two adjacent digit positions `j` and `j + 1` through the phase pair `(v_1, v_2) = (b^j t + alpha, b^(j+1) t + alpha)`. The square is the torus of all such pairs, each phase read in `[-1/2, 1/2)`, painted by `|g_F(v_1) g_F(v_2)|/fill^2`, the share of the full mass that the two positions keep, with `g_F(v) = sum_(d in F) e(d v)` the digit transform: it is `1` only at the centre, where both phases are integers, and the cross through the centre holds the pairs where one of them is. As `t` runs, the pair runs along the line `v_2 = b v_1 - (b - 1) alpha`, ten segments of slope `10` on the torus. On the lattice, at `alpha = 0`, the dim segments pass through the centre, and the frequency `t = 0` keeps the full mass at every pair of positions. At `alpha = 1/2` the orange segments are moved by `(b - 1) alpha = 9/2`, half a turn, and pass beside the centre: at every point of them one phase is at least `beta/(b + 1) = 1/22` from the integers, `beta = norm((b - 1) alpha) = 1/2` the distance from `(b - 1) alpha` to the nearest integer, so the pair keeps at most `1 - c` of the mass, `c >= 0.1500` here (Fact 3.5). That fixed loss per pair of positions, at every frequency, is what pays the major arcs. Base `10` is far below the range of the theorem; the picture is the lemma, not the bound.

Fix a base `b >= 3`, a set `E` of `m >= 1` digits and the rest `F = {0, ..., b - 1}` less `E`, with `fill = b - m`. `S_F` is the set of positive integers whose base-`b` digits all lie in `F`, `A_F(x)` its count up to `x`, `s_b(n)` the sum of the base-`b` digits of `n`, `e(t) = exp(2 pi i t)` and `norm(t)` the distance from `t` to the nearest integer. The certificate (W) reads two explicit one-step constants, `P_b(m)` and, at one missing digit `e_0`, `P'_b(e_0)`, recalled in Definition 2.2 from Definition 2.4 of [An Unconditional Mertens Bound at Large Base](../../coprime/papers/unconditional-mertens-at-large-base.md), the large-base paper below.

**Theorem 1.1.** Suppose

```
(P)   F contains two consecutive digits ,
(W)   P_b(m) < fill b^(-4/5) ,   or   m = 1, E = {e_0}, b >= 36 and P'_b(e_0) < (b - 1) b^(-4/5) ,
```

and let `alpha_1 = log_b(b P_W/fill)`, `P_W` either constant (W) reads, so that `alpha_1 < 1/5`. Fix `beta_0` with `0 < beta_0 <= 1/2`, let `c = c(b, F, beta_0)` be the constant (1.1), `gamma = log(1/(1 - c))/(2 log b)`, and put

```
sigma = gamma (1 - 4 alpha_1)/(5 - 4 alpha_1) ,   eta = min(1/5 - alpha_1, sigma) ,   delta_0 = (eta/4) log(fill)/log(b) .
```

Then for every `delta < delta_0` there is `C > 0`, depending on `b`, `F`, `beta_0` and `delta` alone and effectively computable, such that for every `x >= 2` and every real `alpha` with `norm((b - 1) alpha) >= beta_0`

```
|sum_(n in S_F, n <= x) mu(n) e(alpha s_b(n))| <= C A_F(x) x^(-delta) ,
|sum_(n in S_F, n <= x) Lambda(n) e(alpha s_b(n))| <= C A_F(x) x^(-delta) .
```

The constant is explicit. With `x_0 = beta_0/(b + 1)`, `x_f = 1/(2(b - 1))`, `K_b(x) = max(sin(pi b x)/sin(pi x), 1/sin(pi/b))` and `V_F = fill^(-1) sum_(d in F) (d - mean(F))^2` the variance of the digits of `F`,

```
c = 1 - min(G_1, G_2, G_3)/fill ,     G_1 = fill - 4 sin^2(pi x_0/2) ,     G_2 = m + K_b(x_0) ,
G_3 = max(fill sqrt(1 - 16 x_0^2 V_F), min(m + K_b(x_f), fill - 4 sin^2(pi x_f/2))) ,          (1.1)
```

the pair, kernel and variance forms of Lemma 3.1, and `0 < c <= 1 - 1/pi`. The hypothesis `norm((b - 1) alpha) >= beta_0` keeps `alpha` off the lattice `(b - 1)^(-1) Z`, where the method gives nothing: on the lattice the modulation is a character of `n mod b - 1` and the modulated transform keeps its full mass at some frequency at every depth (Proposition 3.3). (W) forces (P) (Lemma 2.3), which is listed because the proof uses it in its own right.

**Corollary 1.2 (where it holds).** At one missing digit, (P) and (W) hold at every base `b >= 92317` and every missing digit, so Theorem 1.1 holds on every one-missing-digit set from base `92317`. At `m >= 2` missing digits (W) forces `m < b^(2/5)`, and at `b = 10^7` it holds exactly for `m <= 176`.

**Corollary 1.3 (the digit sums of primes).** Let `F` satisfy (P) and (W), let `q >= 2` be coprime to `b - 1` and let `delta` be below the `delta_0` of Theorem 1.1 at `beta_0 = 1/q`. Write `psi_F(x) = sum_(n in S_F, n <= x) Lambda(n)` and `pi_F(x)` for the number of primes of `S_F` up to `x`. Then for every residue `a mod q`

```
sum_(n in S_F, n <= x, s_b(n) = a mod q) Lambda(n) = psi_F(x)/q + O(A_F(x) x^(-delta)) ,
#{p in S_F : p <= x, s_b(p) = a mod q} = (1 + o(1)) pi_F(x)/q .
```

The saving is tiny. At base `92317` and `beta_0 = 0.2` or `1/2` the exponent is `delta_0 >= 4.678 * 10^-9`, set by the margin `1/5 - alpha_1 >= 1.871 * 10^-8` by which (W) holds at its wall; it is `2.819 * 10^-4` at base `10^6`, and it is always below `log(fill)/(20 log b)` (Fact 6.1, `lab/py/digit-sum-twist`, verb `lemma`). No computation at any reachable height can see such a power. What the theorem shows is the mechanism, that the major arcs need no L-function once the modulation is present, and its corollary for the primes.

The proof, after blocks carry `x` off the powers of the base, expands each block sum in additive characters modulo `y = b^k` and cuts the frequencies `a/y` by Dirichlet approximation, as the large-base paper does, but into three regions. Region A, the minor arcs, pays the whole `l^1` mass of the modulated transform against the minor-arc bound `x^(4/5 + eps)` of Basak, Robles and Zaharescu (2023), or Lemma 4.2 of Maynard (2022) for `Lambda`. Region B, the middle denominators, pays a hybrid `l^1` mass against the decay `max(d, h)^(-1/2)` of the same bounds. Region C holds the frequencies within `Z/y` of a fraction of denominator below `Z`, the major arcs. A and B read the transform only through the one-step constant `B_b(F)`, a supremum over every shift, and the modulation only moves the shift (Lemma 4.1), so both are imported, the minor-arc bounds re-read with the ratio `r = x/y`. In the unmodulated dissection C splits into C1, where the denominator has a prime outside the base and the transform is small, and C2, where it divides a power of the base, the transform is full and the block sum is `mu` in progressions paid by a zero-free region; both are paid at root-log strength, and they cap the saving at `exp(-c sqrt(log x))`. With the modulation all of C is paid by the supremum of the transform, `fill^k y^(-gamma)` up to a constant (Corollary 3.2), times its `3 Z^2` frequencies times the trivial bound on the block sum; so `Z` can be a power of `y`, the blocks can be cut at a depth that is a power of `x`, and the saving is a power.

Mauduit and Rivat (2010) prove `sum_(n <= x) Lambda(n) e(alpha s_q(n)) = O(x^(1 - sigma_q(alpha)))` whenever `(q - 1) alpha` is not an integer, for unrestricted digits, by Type I and Type II sums with a carry truncation, and Martin, Mauduit and Rivat (2014) extend it to digital functions with real coefficients; neither restricts the digits. The survey of Maynard (2019) proves a weak form of the Mauduit-Rivat theorem at large base by the circle method, with the full-digit two-position bound and the major arcs paid by supremum times count, so the mechanism of this paper is in print for the full digit set (Section 9 says exactly how). What is written here is its object: the modulation carried onto a missing-digit set, for `mu` and `Lambda`, with an effective constant, from an explicit base, and with a two-position constant that does not decay with `b`. For the integers of a missing-digit set, Aloui (2015), read in its abstract only, carries the digit-sum classes with power-saving errors (Section 9). A search for a digit-sum theorem for the primes or `mu` of a missing-digit set found none; the search is stated in Section 9.

Section 2 fixes the objects and the imports. Section 3 proves the two-position lemma, its sharpness and the lattice dichotomy; Section 4 the modulated peel and the modulated hybrid bound; Section 5 the dissection, each changed step written out and each imported step stated with its pointer, and Theorem 1.1. Section 6 is the exponent and which term binds where, with Corollary 1.2; Section 7 is Corollary 1.3 and why it needs (W) itself; Section 8 the falsification; Section 9 the literature; Section 10 what is left open.

## Definitions and imports

**Definition 2.1 (the modulated transform).** For `k >= 0`, `D_k` is the set of the `fill^k` integers `0 <= n < b^k` whose `k` padded digits lie in `F`. The digit transform is `g_F(v) = sum_(d in F) e(d v)`, written `hat F` in the large-base paper. For real `alpha` the modulated transform at depth `k` is

```
hat w_k(t) = sum_(n in D_k) e(alpha s_b(n) + n t) = prod_(j < k) g_F(b^j t + alpha) ,
```

since with `n = sum_(j < k) d_j b^j` the phase is `sum_j d_j (b^j t + alpha)` and the digits `d_j` run over `F` independently; it has period `1`, and at `alpha = 0` it is the level transform `hat F_k` of the large-base paper. Put `beta = norm((b - 1) alpha)`; the lattice is `(b - 1)^(-1) Z`, where `beta = 0`. The modulated grid sums are `Sigma_i(s) = sum_(a < b^i) |hat w_i(s + a/b^i)|` at real `s`, `Sigma_0 = 1`. `V_F`, `K_b` and `x_f` are as in (1.1).

**Definition 2.2 (the one-step constant and (W)).** `B_b(F) = sup_t sum_(r mod b) |g_F((t + r)/b)|`, a supremum over every shift `t`. With `H(n) = log n + gamma_E + 1/(2n)`, `gamma_E` Euler's constant,

```
Phi_b = (4/pi) b + (2b/pi) H(ceil((b - 2)/2)) + (1 - 2/pi)(b - 2) + 0.727 ,     P_b(m) = sqrt(m) + Phi_b/b ,
```

and at one missing digit `e_0` and `b >= 36`, with `u_0 = e_0 - (b - 1)/2` and `p = floor(b/2)`, the chord constant of Definition 2.4 of the large-base paper is

```
Psi'_b = (b/pi)(2 H(p - 1) - 1 + 1/p) + (1 - 2/pi) b/2                       at even b ,
Psi'_b = (b/pi)(2 H(p - 1) - 1 + 2/p) + (1 - 2/pi)(b/2 + 1/(2b))              at odd b ,
P'_b(e_0) = ((4/pi) b + Psi'_b + b/2 - sec(pi u_0/b)/2)/b ,
```

below `P_b(1)` at every `b >= 36` and every `e_0`, and smallest at the end digits. (W) is the condition displayed in Theorem 1.1, `P_W` either constant it reads and `alpha_1 = log_b(b P_W/fill)`; (W) holds exactly when `alpha_1 < 1/5`.

**Lemma 2.3 (imported).** (i) `B_b(F) >= b` for every digit set. (ii) Under (W), `B_b(F) <= b P_W = fill b^(alpha_1)`. (iii) (W) forces (P), and at `m >= 2` it forces `b^(1/5) > 2` and `m < b^(2/5)`.

**Proof.** (i) is Lemma 2.3(iii) of the large-base paper: the digits of `F` are distinct modulo `b`, so Parseval on `Z/b` gives `sum_(r mod b) |g_F((t + r)/b)|^2 = b fill`, and each term is at most `fill`. (ii) is its Lemmas 2.5 and 2.6. (iii) is its Lemma 6.2: at `m = 1` two consecutive digits remain; at `m >= 2`, `P_b(m) >= sqrt(m) + 4/pi > 2` while (W) gives `P_b(m) < b^(1/5)`, so `b^(1/5) > 2` and `m < b^(2/5) < (b - 1)/2`, and a set with no two consecutive digits misses at least `floor(b/2)` digits. □

## The two-position lemma

**Lemma 3.1 (two positions).** Let `F` contain two consecutive digits, let `alpha` be real with `beta = norm((b - 1) alpha)`, put `x_0 = beta/(b + 1)` and let `G_F(x_0)` be the supremum of `|g_F(v)|` over `norm(v) >= x_0`. Then for every real `u`

```
|g_F(u + alpha) g_F(b u + alpha)| <= fill G_F(x_0) ,     G_F(x_0) <= min(G_1, G_2, G_3) ,
```

with `G_1`, `G_2` and `G_3` the forms of (1.1) at this `x_0`.

**Proof.** Locating. Put `v_1 = u + alpha` and `v_2 = b u + alpha`. Then `b v_1 - v_2 = (b - 1) alpha`, so `beta = norm(b v_1 - v_2) <= b norm(v_1) + norm(v_2) <= (b + 1) max(norm(v_1), norm(v_2))`, and one of the two phases has `norm(v_i) >= x_0`. Its factor is at most `G_F(x_0)` and the other at most `fill`.

The pair form. With `d` and `d + 1` in `F`, `|g_F(v)| <= fill - 2 + |e(d v) + e((d + 1) v)| = fill - 2 + 2 |cos(pi v)| = fill - 4 sin^2(pi norm(v)/2)`, which falls as `norm(v)` rises; so `G_F(x_0) <= G_1`.

The kernel form. `g_F = D_b - g_E` with `D_b(v) = sum_(n < b) e(n v)` and `|g_E| <= m`. `|D_b|` is even and has period `1`, so put `x = norm(v)`. On `(0, 1/b]`, `e(-(b - 1) x/2) D_b(x) = sum_(n < b) cos(2 pi (n - (b - 1)/2) x)`, the sines cancelling in the pairs `n` and `b - 1 - n`; each term falls in `x` because `2 pi |n - (b - 1)/2| x <= pi (b - 1)/b < pi`, and the sum is `sin(pi b x)/sin(pi x) >= 0`, so `|D_b(x)| = sin(pi b x)/sin(pi x)` falls on `(0, 1/b]`. On `[1/b, 1/2]`, `|D_b(x)| <= 1/sin(pi x) <= 1/sin(pi/b)`. Since `x_0 <= 1/(2(b + 1)) < 1/b`, the supremum of `|D_b|` over `x >= x_0` is at most `K_b(x_0)`, and `G_F(x_0) <= m + K_b(x_0) = G_2`.

The variance form. `|g_F(v)|^2 = sum_(d, d' in F) cos(2 pi (d - d') x) = fill^2 - 2 sum_(d, d') sin^2(pi (d - d') x)`. At `x <= x_f` every `|d - d'| x` is at most `1/2`, where `sin^2(pi y) >= 4 y^2` by the concavity of `sin` on `[0, pi/2]`, and `sum_(d, d') (d - d')^2 = 2 fill^2 V_F`; so `|g_F(v)| <= fill sqrt(1 - 16 x^2 V_F)` on `[0, x_f]`, the radicand nonnegative since `V_F <= (b - 1)^2/4`, and the bound falls in `x`. As `x_0 < x_f`, the supremum over `x_0 <= x <= x_f` is at most `fill sqrt(1 - 16 x_0^2 V_F)`, and over `x >= x_f` the pair and kernel forms at `x_f` apply, since `x_f <= 1/b`. So `G_F(x_0) <= G_3`. □

**Corollary 3.2 (the supremum).** Under the hypotheses of Lemma 3.1 put `c = 1 - min(G_1, G_2, G_3)/fill` and `gamma = log(1/(1 - c))/(2 log b)`. If `beta > 0`, then `0 < c <= 1 - 1/pi`, `gamma < 0.53`, and at every `k >= 0` and every real `t`

```
|hat w_k(t)| <= fill^k (1 - c)^(floor(k/2)) <= (1 - c)^(-1/2) fill^k b^(-gamma k) .
```

`c` does not fall as `beta` rises, so the bound at `beta_0` holds at every `alpha` with `norm((b - 1) alpha) >= beta_0`.

**Proof.** `G_1 < fill` because `x_0 > 0`. Each form is at least `min(fill - 1, m + b/pi)`: `G_1 >= fill - 4 sin^2(pi/16) > fill - 1`, `G_2 >= m + 1/sin(pi/b) >= m + b/pi`, and `G_3` is at least its second entry. That minimum is at least `fill/pi`, since `fill >= 2` and `b >= fill`; so `1 - c >= 1/pi` and `gamma <= log(pi)/(2 log 3) < 0.53`. Pair the positions `2i` and `2i + 1` of the product of Definition 2.1 with `u = b^(2i) t`, bound an odd last factor by `fill`, and use `floor(k/2) >= (k - 1)/2`. Each form is nonincreasing in `x_0`, `G_2` because `sin(pi b x)/sin(pi x)` falls on `(0, 1/b]`, so `c` is nondecreasing in `beta`. □

**Proposition 3.3 (the lattice).** If `(b - 1) alpha` is an integer, then `e(alpha s_b(n)) = e(alpha n)` for every `n >= 0`, `|hat w_k(-alpha)| = fill^k` at every `k`, and

```
sum_(n in S_F, n <= x) mu(n) e(alpha s_b(n)) = sum_(a mod b - 1) e(alpha a) M_F(x; b - 1, a) ,
```

with `M_F(x; b - 1, a)` the sum of `mu(n)` over the `n <= x` of `S_F` with `n = a mod b - 1`. So with two consecutive digits `sup_t |hat w_k(t)|/fill^k` is `1` at every depth when `alpha` is on the lattice, and at most `(1 - c)^(floor(k/2))` with `c > 0` when it is off.

**Proof.** `n - s_b(n) = sum_j d_j (b^j - 1)`, and each `b^j - 1` is a multiple of `b - 1`. At `t = -alpha` every phase `b^j t + alpha = -(b^j - 1) alpha` is an integer, where `|g_F| = fill`. The identity follows from the first statement, and the decay off the lattice is Corollary 3.2. □

The mechanism is the phase map. The phases `v_j = b^j t + alpha` obey `v_(j+1) = b v_j - (b - 1) alpha`. With two consecutive digits `|g_F(v)| = fill` exactly when `v` is an integer, so full mass at depth `k` needs `v_0, ..., v_(k-1)` all integral, and the affine map carries the integers into themselves exactly when `(b - 1) alpha` is one. At `alpha = 0` the frequencies of full mass are the `b`-adic preimages `l/b^j` of `0`: the region C2 of the large-base paper, the one place its zero-free region enters (its Lemma 6.5 and Proposition 6.6). Off the lattice no frequency keeps the mass. On the lattice the method returns nothing: the sum is a combination of `mu` in progressions modulo `b - 1` over `S_F`, whose major arcs sit at `-alpha + l/b^j` and would need zero-free regions for the characters modulo `(b - 1) b^j`. That case is carried at `exp(-c sqrt(log x))` in the section Progressions and every modulus of [the digit-sum twist note](../notes/progressions-every-modulus.md).

**Example 3.4 (the consecutive pair is needed).** At `b = 5`, `F = {0, 2}` and `alpha = 1/8`, off the lattice with `beta = 1/2`, `|hat w_k(3/8)| = fill^k` at every `k`.

**Proof.** The phases are `5^j (3/8) + 1/8 = (3 * 5^j + 1)/8`, and `5^j` is `1` or `5` modulo `8`, so `3 * 5^j + 1` is `4` or `0` modulo `8` and every phase lies in `Z/2`, where `g_F(v) = 1 + e(2v) = 2 = fill`. The generator reads `|hat w_k(3/8)|/fill^k = 1` to twelve digits at `k = 4, 8, 12` (`lab/py/digit-sum-twist`, verb `lemma`). □

Without the pair, `|g_F|` is full on a larger set than the integers, here `Z/2`, and the locating step finds one phase far from the integers without finding it far from that set. The pair form is the step that rules this out.

**Fact 3.5 (the constants).** At base `10` missing `7`, `c >= 0.1500` and `gamma >= 0.0352` at `alpha = 1/2`, and `c >= 0.0463` at `1/7`, `0.0353` at `1/4`, `0.0224` at `1/5`, `0.0055` at `1/10` and `5.55 * 10^-5` at `1/900`, where `beta = 0.01`, each from the variance form, which beats the pair form by a factor `66` at `1/2` and `61` at `1/900`. At base `92317` missing its extreme digit, whose variance is the least, the kernel form gives `c >= 0.3633` at `beta = 1/2` and `0.06448` at `0.2`, and the variance form `c >= 6.666 * 10^-7` at `beta = 10^-3`. Generator: `lab/py/digit-sum-twist`, verb `lemma`.

To leading order the variance form gives `c = 8 x_0^2 V_F`; since `V_F` is about `b^2/12` and `x_0` about `beta/b`, that is `c` about `(2/3) beta^2` at small `beta` and large base, independent of `b`, where the pair form alone decays like `beta^2 b^(-3)`. The full-digit bound in print (Section 9) has the shape of the pair form.

## The modulated peel

**Lemma 4.1 (the modulated peel).** For every real `alpha` and `s` and every `i >= 1`, `Sigma_i(s) <= B_b(F) Sigma_(i-1)(b s)`, so `Sigma_i(s) <= B_b(F)^i`. The same grid sum with the factor of any one position replaced by `1` is at most `b B_b(F)^(i-1)`, and `int_0^1 |hat w_i| <= b^(-i) B_b(F)^i`. Under (W), `Sigma_i(s) <= fill^i b^(i alpha_1)`: the shifted-grid certificate of the large-base paper with `C_F = 1`, at every shift, for the modulated transform.

**Proof.** `hat w_i(t) = g_F(t + alpha) hat w_(i-1)(b t)`, and `hat w_(i-1)` has period `1`. Write `a = a' + r b^(i-1)` with `a' < b^(i-1)` and `r < b`. Then `b (s + a/b^i) = b s + a'/b^(i-1) + r`, so the second factor depends on `a'` alone, and `s + a/b^i + alpha = (T + r)/b` with `T = b (s + alpha) + a'/b^(i-1)`. Summing over `r` first gives at most `sum_(r mod b) |g_F((T + r)/b)| <= B_b(F)` times `|hat w_(i-1)(b s + a'/b^(i-1))|`, and summing over `a'` gives `B_b(F) Sigma_(i-1)(b s)`; induct from `Sigma_0 = 1`. If the factor at position `0` is replaced by `1`, the sum over `r` is exactly `b`; if it sits at position `j > 0`, peel the `j` positions below it as above, then that one at the cost `b`, then the rest; this gives `b B_b(F)^(i-1)`. The integral is `int_0^(b^(-i)) Sigma_i(s) ds`. Under (W) apply Lemma 2.3(ii). The modulation moves only the shift `T`, and `B_b(F)` is a supremum over every shift. □

**Lemma 4.2 (the hybrid `l^1` bound, modulated).** Let `y = b^k` and give each residue `a mod y` a reduced fraction `l/d` with `1 <= d <= Q = y^(3/5)` and `|a/y - l/d| <= 1/(d Q)`, by Dirichlet's theorem, and its height `h = |a d - l y| <= y^(2/5)`. Let `D >= 1` and `H >= 0` with `16 b^2 D (D + H) <= y`, let `R(D, H)` be the residues with `D <= d < 2D` and `h < 2H`, or `h = 0` when `H = 0`, and let `V_1 = b^(i_1)` be the least power of `b` at least `4 D^2` and `V_2 = b^(i_2)` the least at least `4H/D + 1`. Then under (W), at every real `alpha`,

```
sum_(a in R(D, H)) |hat w_k(a/y)| <= (1 + pi b) fill^k (V_1 V_2)^(alpha_1) .
```

**Proof.** This is Lemma 5.1 of the large-base paper, whose proof reads the transform through three things only: a split into three products, the certificate at a shift, and a bound on the derivative of one digit factor. The split holds for the modulated transform: grouping the positions `j < i_1`, the middle ones and the top `i_2` gives `hat w_k(t) = hat w_(i_1)(t) hat w_(k - i_1 - i_2)(V_1 t) hat w_(i_2)(y t/V_2)`, since `b^(k - i_2) = y/V_2`, and the middle factor is at most `fill^(k - i_1 - i_2)`. At `t = a/y` the top factor is `hat w_(i_2)(a/V_2)`. The residues attached to one fraction are at most `V_2` consecutive integers, distinct modulo `V_2`, so they pay at most `Sigma_(i_2)(0) <= B_b(F)^(i_2)` by Lemma 4.1. The bottom factor `f = hat w_(i_1)` at each residue of `l/d` is at most its supremum over `J_(l,d) = [l/d - 1/(8 D^2), l/d + 1/(8 D^2)]`, and by `|f(t)| <= |J|^(-1) int_J |f| + int_J |f'|` on the disjoint intervals `J_(l,d)` the suprema sum to at most `4 D^2 ||f||_1 + ||f'||_1`, norms on `[0, 1]`. Lemma 4.1 gives `||f||_1 <= V_1^(-1) B_b(F)^(i_1)`. Differentiating one factor at a time, `|d/dt g_F(b^j t + alpha)| <= b^j 2 pi sum_(d in F) d <= b^j pi b (b - 1)`, and the integral over `[0, 1]` of the product with the factor at position `j` removed is `int_0^(1/V_1)` of its grid sum on `V_1` points, at most `V_1^(-1) b B_b(F)^(i_1 - 1)` by Lemma 4.1; with `sum_(j < i_1) b^j < V_1/(b - 1)` that is `||f'||_1 <= pi b^2 B_b(F)^(i_1 - 1) <= pi b B_b(F)^(i_1)`, by Lemma 2.3(i). So the fractions pay at most `(1 + pi b) B_b(F)^(i_1)`, since `4 D^2 <= V_1`, and `B_b(F)^(i_1 + i_2) <= fill^(i_1 + i_2) (V_1 V_2)^(alpha_1)` by Lemma 2.3(ii). The count of the residues of one fraction, the containment of those residues in `J_(l,d)` and the disjointness of the intervals use `16 D H <= y` and the fractions alone, and are as in Lemma 5.1 there. □

**Remark 4.3 (what does not transfer).** The large-base paper proves its certificate from base `584` by a digit-uniform chain (its Proposition 8.5), which telescopes a run of positions into one Dirichlet kernel, `prod_(j < l) D_b(b^j u) = D_(b^l)(u)`; that needs consecutive arguments in ratio `b`, and the modulated arguments `b^j t + alpha` are not, since `b (b^j t + alpha) = b^(j+1) t + b alpha` differs from `b^(j+1) t + alpha` modulo `1` unless `beta = 0`. Its window certificates (Facts 8.6 and 8.7) bound grid sums of `prod_j |g_F(b^j t)|` over the one-parameter family of shifts, while a modulated grid sees the shift vector `(b^j s + alpha)_j`. Neither is refuted; neither is carried. Only the one-step bound of Lemma 4.1 survives the modulation, and so the base of Theorem 1.1 is the base of (W).

## The dissection

The dissection is that of the large-base paper for `mu` and of Section 7 of [An Unconditional Mertens Bound on the Half Interval](../../mobius/papers/half-interval-mobius.md), the half-interval paper below, for `Lambda`. Three things change: the blocks are cut at a depth that is a power of `x`, the major-arc parameter `Z` is a power of the block length, and the major arcs, one region C in place of C1 and C2, are paid by the two-position lemma. Imported by statement: the block split of Lemma 3.2 of the large-base paper, the bound of Basak, Robles and Zaharescu (its Lemma 4.2), the first display of Lemma 7.4 of the half-interval paper, and the classes of region B of Proposition 5.2 of the large-base paper. Imported by proof, with the ratio `r` in place of `x^(o(1))`: the two approximations of Lemma 4.3 there and the second display of Lemma 7.4 (Lemma 5.4). Not used: the contraction near a fraction with a prime outside the base (Lemma 6.1 there), `mu` and `Lambda` in progressions (Lemma 6.5 there, Lemma 7.5 of the half-interval paper), the main term of Lemma 7.2 of the half-interval paper, and every character, zero-free region and exceptional zero those lemmas carry.

**Lemma 5.1 (the modulated blocks).** Let `f` be `mu` or `Lambda`, with `f(0) = 0`. For `k >= 0`, `y = b^k` and an integer `P >= 0` put `S_P(theta) = sum_(P y <= n < (P+1) y) f(n) e(n theta)`. Then

```
sum_(t in D_k) e(alpha s_b(P y + t)) f(P y + t) = e(alpha s_b(P)) Sigma(P, k) ,     Sigma(P, k) = y^(-1) sum_(a mod y) hat w_k(a/y) S_P(-a/y) ,
```

so `|Sigma(P, k)| <= y^(-1) sum_(a mod y) |hat w_k(a/y)| |S_P(a/y)|`.

**Proof.** For `0 <= t < y` the digits of `P y + t` are the `k` padded digits of `t` below those of `P`, so `s_b(P y + t) = s_b(P) + s_b(t)`. Completeness of the additive characters modulo `y`, applied to `t -> 1_(D_k)(t) e(alpha s_b(t))`, whose transform at `a/y` is `hat w_k(a/y)`, gives `1_(D_k)(t) e(alpha s_b(t)) = y^(-1) sum_(a mod y) hat w_k(a/y) e(-a t/y)`; put `n = P y + t`, so `e(-a t/y) = e(-a n/y)`. `f` is real, so `|S_P(-theta)| = |S_P(theta)|`. □

**Lemma 5.2 (the split at power depth).** Let `x >= b` have `L` digits and `0 < tau < 1`. By Lemma 3.2 of the large-base paper, `S_F` below `x` is the point `x`, when its digits lie in `F`, together with sets `P b^k + D_k` at scales `k < L`, the element `0` discarded, at most `fill + 1` of them at each scale, and `A_F(x) >= fill^(L-1) - 1`. The sets at scales with `b^k < x^(1 - tau)` hold at most `2 (fill + 1) fill (A_F(x) + 1) x^(-tau log(fill)/log(b))` elements together, and every other set has `y = b^k >= x^(1 - tau)`, that is `r = x/y <= x^tau`.

**Proof.** With `k_0` the largest such scale, those sets hold at most `(fill + 1) sum_(k <= k_0) fill^k <= 2 (fill + 1) fill^(k_0)` elements, and `fill^(k_0) = (b^(k_0))^(log(fill)/log(b)) < x^((1 - tau) log(fill)/log(b))`. Since `b^(L-1) > x/b`, the mass floor gives `A_F(x) + 1 >= fill^(L-1) > (x/b)^(log(fill)/log(b)) = x^(log(fill)/log(b))/fill`. □

**Definition 5.3 (three regions).** On a set `P y + D_k` with `y >= x^(1 - tau)`, give each residue `a mod y` the fraction `l/d` and height `h` of Lemma 4.2, at `Q = y^(3/5)`, and put `Z = y^zeta` with `zeta = 2 gamma/(5 - 4 alpha_1)`. Region A is `d >= y^(2/5)`; region B is `d < y^(2/5)` with `max(d, h) >= Z`; region C is `d < Z` and `h < Z`. By Corollary 3.2, `gamma < 0.53` and `5 - 4 alpha_1 > 4.2`, so `zeta < 0.26 < 2/5` and the three regions partition the residues. By the choice of `zeta`, `gamma - 2 zeta = sigma` and `zeta (1/2 - 2 alpha_1) = sigma`.

**Lemma 5.4 (the minor arcs, with the ratio).** On a set with `r = x/y`, for `f = mu` and every fixed `eps > 0`,

```
|S_P(a/y)| <<_eps x^(4/5 + eps) r^(1/5)                                   on region A ,
|S_P(a/y)| <<_eps x^(4/5 + eps) r^(1/5) + x (log x)^3 max(d, h)^(-1/2)      on region B ,
```

the implied constant effective and depending on `eps` alone, and for `f = Lambda` the same with `(log x)^4` in place of both `x^eps` and `(log x)^3` and an absolute constant.

**Proof.** For `mu` this is Lemma 4.3 of the large-base paper, whose proof applies the bound of Basak, Robles and Zaharescu (2023), Theorem 1.4, its Lemma 4.2, `|sum_(n <= X) mu(n) e(n theta)| <<_eps X^(4/5 + eps) + X d^(-1/2) (log X)^3 + (X d)^(1/2) (log X)^3` at `|theta - l/d| <= 1/d^2`, `(l, d) = 1`, to `S_P`, a difference of two such sums with `X <= 2x`, at `l/d` and, on region B with `h >= 1`, at a second fraction `l'/d'` with `y/(2h) <= d' <= 2y/h`. There `x/y` was at most `b^K = x^(o(1))`; here it is `r`. On region A, `X d^(-1/2) <= 2 x y^(-1/5) = 2 x^(4/5) r^(1/5)` and `(X d)^(1/2) <= (2 x y^(3/5))^(1/2) <= 2 x^(4/5)`. On region B, `l/d` leaves `x (log x)^3 d^(-1/2)` and `(2 x y^(2/5))^(1/2) (log x)^3 << x^(7/10) (log x)^3`, and `l'/d'` leaves `2 x h^(-1/2)` and `X d'^(-1/2) <= 2 x (2h/y)^(1/2) <= 3 x y^(-3/10) = 3 x^(7/10) r^(3/10)`, as `h <= y^(2/5)`; and `x^(7/10) r^(3/10) <= x^(4/5) r^(1/5)` because `r <= x`. For `Lambda` this is Lemma 7.4 of the half-interval paper, Lemma 4.2 of Maynard (2022) read at a coarser fraction and at `l/d`, whose proof gives the terms `X^(4/5)`, `3 x y^(-1/5)`, `2 x^(4/5)`, `2 x^(7/10)`, `3 x d^(-1/2)`, `2 x h^(-1/2)` and `2 x y^(-3/10)` times `(log x)^4`, with `r` in place of its `b^K`. □

**Proposition 5.5 (region A).** `y^(-1) sum_(a in A) |hat w_k(a/y)| |S_P(a/y)| <<_eps fill^k r^(1 + eps) y^(alpha_1 - 1/5 + eps)` for `mu`, and `<< fill^k r y^(alpha_1 - 1/5) (log x)^4` for `Lambda`.

**Proof.** By Lemma 4.1 at `s = 0`, `sum_(a mod y) |hat w_k(a/y)| = Sigma_k(0) <= fill^k y^(alpha_1)`, and by Lemma 5.4 every residue of A carries `|S_P| <<_eps (r y)^(4/5 + eps) r^(1/5)`. □

**Proposition 5.6 (region B).** Once `y >= (32 b^2)^5`, `y^(-1) sum_(a in B) |hat w_k(a/y)| |S_P(a/y)| <<_eps fill^k (log y)^2 (r^(1 + eps) y^((4/5) alpha_1 - 1/5 + eps) + r (log x)^4 y^(-sigma))`, and the same for `Lambda` with `(log x)^4` in place of `r^eps y^eps`.

**Proof.** Group the residues of B into the classes of Proposition 5.2 of the large-base paper, `D <= d < 2D` over powers of two `D < y^(2/5)`, and `H <= h < 2H` over powers of two `H <= y^(2/5)` or `h = 0`: at most `(log_2 y + 2)^2` classes, each inside `R(D, H)`, with `16 b^2 D (D + H) <= 32 b^2 y^(4/5) <= y`, and on each `max(d, h) >= max(D, H) > Z/2`. By Lemma 4.2, with `V_1 V_2 < 16 b^2 D (D + H) <= 32 b^2 max(D, H)^2`, and Lemma 5.4, a class weighs `<<_eps y^(-1) (x^(4/5 + eps) r^(1/5) + x (log x)^4 max(D, H)^(-1/2)) fill^k max(D, H)^(2 alpha_1)`. With `max(D, H) <= y^(2/5)` the first term is `<<_eps fill^k r^(1 + eps) y^((4/5) alpha_1 - 1/5 + eps)`. With `2 alpha_1 - 1/2 < 0` and `max(D, H) > Z/2` the second is `<< fill^k r (log x)^4 Z^(-(1/2 - 2 alpha_1)) = fill^k r (log x)^4 y^(-sigma)`. □

**Proposition 5.7 (region C, with no L-function).** For `f = mu` or `Lambda` and every `alpha` with `norm((b - 1) alpha) >= beta_0`,

```
y^(-1) sum_(a in C) |hat w_k(a/y)| |S_P(a/y)| <= 3 (1 - c)^(-1/2) fill^k y^(-sigma) log(2x) .
```

**Proof.** The residues attached to a fraction `l/d`, read modulo `1`, have `|a - l y/d| < Z/d`, so region C holds at most `sum_(d < Z) phi(d) (2Z/d + 1) <= 3 Z^2` residues. On each, `|S_P(a/y)| <= y log(2x)`, the block holding `y` terms each at most `log(2x)`, and `|hat w_k(a/y)| <= (1 - c)^(-1/2) fill^k y^(-gamma)` by Corollary 3.2. So region C is at most `3 (1 - c)^(-1/2) fill^k Z^2 y^(-gamma) log(2x) = 3 (1 - c)^(-1/2) fill^k y^(2 zeta - gamma) log(2x)`, and `gamma - 2 zeta = sigma`. □

No contraction lemma, no progression, no Perron integral and no exceptional zero enters region C: the regions C1 and C2 of the large-base paper are paid together by the supremum, and the modulated sums carry no main term.

**Proof of Theorem 1.1.** Fix `beta_0`, `delta < delta_0` and `alpha` with `norm((b - 1) alpha) >= beta_0`; `c` and `gamma` are those of Corollary 3.2 at `beta_0`, which serve at `alpha`. Put `tau = eta/4`, take `eps > 0` small in terms of `delta_0 - delta`, and let `x >= b`; below `b` the bound is trivial. By Lemmas 5.1 and 5.2 the modulated sum is at most `log x`, for the point `x`, plus the sum of `|Sigma(P, k)|` over the sets of the split. The sets at scales with `b^k < x^(1 - tau)` weigh `<< A_F(x) x^(-tau log(fill)/log(b)) log x = A_F(x) x^(-delta_0) log x`, by Lemma 5.2 and `f <= log x`. Every other set has `y >= x^(1 - tau)`, large once `x` is, and `r <= x^tau`, and by Propositions 5.5, 5.6 and 5.7, with `alpha_1 - 1/5 <= -eta`, `(4/5) alpha_1 - 1/5 <= alpha_1 - 1/5`, as `alpha_1 > 0` by Lemma 2.3(i) and (ii), and `sigma >= eta`,

```
|Sigma(P, k)| <<_eps fill^k r^(1 + eps) y^(-eta + eps) (log x)^6 .
```

As `r^(1 + eps) y^(-eta + eps) = x^(1 + eps) y^(-1 - eta)` falls in `y`, it is at most `x^(tau (1 + eps) - (1 - tau)(eta - eps))`, whose exponent is at most `-3 eta/4 + eta^2/4 + 2 eps <= -eta/2 + 2 eps` because `eta <= 1/5`. At most `fill + 1` sets sit at each scale, so these sets carry `(fill + 1) sum_(k < L) fill^k <= (fill + 1) fill (A_F(x) + 1)` times that bound. So the modulated sum is `<<_eps A_F(x) (x^(-eta/2 + 2 eps) + x^(-delta_0)) (log x)^6`, and `delta_0 < eta/2` gives `C A_F(x) x^(-delta)` once `eps` is small. Below the effective point past which every step that asked `y` large holds, the trivial bound `A_F(x) log x` goes into `C`. Every constant is effective: the implied constants of the two minor-arc inputs at this `eps`, `(1 - c)^(-1/2)`, `1 + pi b`, `(32 b^2)^(alpha_1)` and the thresholds on `x`, each explicit in `b`, `F`, `beta_0` and `delta`. □

Which region needs what:

| region | paid by | asks |
| --- | --- | --- |
| A: `d >= y^(2/5)` | the peel at `s = 0` against the minor-arc exponent `4/5` | `alpha_1 < 1/5`, the wall |
| B: `d < y^(2/5)`, `max(d, h) >= Z` | the modulated hybrid bound against `max(d, h)^(-1/2)` | `alpha_1 < 1/4` |
| C: `d < Z`, `h < Z` | the two-position lemma times `3 Z^2` | two consecutive digits, `beta > 0` |
| the blocks | depth `tau log_b x`, a power of `x` | nothing |

In the unmodulated dissection the blocks sit at depth `kappa sqrt(log x)` and `Z = exp(C_0 sqrt(log x))`, because regions C1 and C2 are paid at root-log strength and cap both. Here nothing caps `Z` below a power of `y`. The modulation trades the base for the shape: unmodulated, the bound is proved from base `584` at `exp(-c sqrt(log x))`; modulated, from base `92317` at `x^(-delta)`, because only the one-step certificate survives the modulation (Remark 4.3).

## The exponent

The exponent has three sources. Region A saves `y^(alpha_1 - 1/5)`, as small as the margin by which (W) holds. Regions B and C save `y^(-sigma)` each, and `zeta` is chosen to balance them: B pays `Z^(2 alpha_1 - 1/2)`, rising with `Z`, and C pays `Z^2 y^(-gamma)`, falling, and they meet at `Z = y^zeta`. So the kept sets save `eta = min(1/5 - alpha_1, sigma)`, less what the ratio `r <= x^tau` costs. The discarded scales save `x^(-tau log(fill)/log(b))`, and at `tau = eta/4` that term sets `delta_0 = (eta/4) log(fill)/log(b)`; the depth `tau = eta/4` is a convenience and is not optimised. Since `eta <= 1/5 - alpha_1 < 1/5`, `delta_0 < log(fill)/(20 log b)` at every base.

**Fact 6.1 (the exponent).** Each entry is a lower bound, computed at `40` digits and floored at four; `c` and `gamma` are read at the extreme missing digit, whose variance is the least, so each row bounds every one-missing-digit set of its base, and the chord row reads the middle digit, where the chord constant is largest.

| base | constant in (W) | `beta_0` | `1/5 - alpha_1` | `c` | `gamma` | `sigma` | `delta_0` |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `92317` | `P_b(1)` | `1/2` | `1.871 * 10^-8` | `0.3633` | `0.01974` | `9.403 * 10^-4` | `4.678 * 10^-9` |
| `92317` | `P_b(1)` | `0.2` | `1.871 * 10^-8` | `0.06448` | `0.002915` | `1.388 * 10^-4` | `4.678 * 10^-9` |
| `92317` | `P_b(1)` | `10^-3` | `1.871 * 10^-8` | `6.666 * 10^-7` | `2.915 * 10^-8` | `1.388 * 10^-9` | `3.47 * 10^-10` |
| `39363` | `P'_b(e_0)` | `1/2` | `2.696 * 10^-7` | `0.3633` | `0.02133` | `1.015 * 10^-3` | `6.741 * 10^-8` |
| `10^6` | `P_b(1)` | `1/2` | `0.02411` | `0.3633` | `0.01634` | `1.127 * 10^-3` | `2.819 * 10^-4` |
| `10^6` | `P_b(1)` | `0.2` | `0.02411` | `0.0645` | `0.002413` | `1.665 * 10^-4` | `4.163 * 10^-5` |

Region A binds just above the wall: by bisection between `92317` and `10^6`, `1/5 - alpha_1` is the smaller term up to `100083` and `sigma` from `100084` at `beta_0 = 1/2`, and the switch falls at `93409` at `beta_0 = 0.2`; at `beta_0 = 10^-3`, `sigma` binds already at `92317`, and at `10^6` it binds at all three. Generator: `lab/py/digit-sum-twist`, verb `lemma`, under a second.

So past about `10^5` the two-position constant limits `delta_0`, through regions B and C, and below that the margin of (W) does: past the switch `delta_0` is `sigma log(fill)/(4 log b)`, proportional to `gamma = log(1/(1 - c))/(2 log b)`.

**Proof of Corollary 1.2.** At `m = 1`, (P) holds because one missing digit leaves two consecutive ones at `b >= 3`. (W) holds through `P_b(1)` at every `b >= 92317`: by Lemma 8.2 of the large-base paper the gap `(b - 1) b^(-4/5) - P_b(1)` increases at every `b >= 11221`, and it is positive at `92317` and negative at `92316`, certified in interval arithmetic at `200` bits: `1/5 - alpha_1 >= 1.87123612 * 10^-8` at `92317` and `<= -4.82092732 * 10^-8` at `92316` (`lab/py/digit-sum-twist`, verb `lemma`), where Fact 8.1 there reads the same gap as a float re-read at `40` digits; this is Fact 8.3 there, with its one numerical step certified. At `m >= 2`, Lemma 2.3(iii) gives (P) and `m < b^(2/5)`, and Fact 8.4 of the large-base paper gives `P_b(m) < (b - m) b^(-4/5)` at `b = 10^7` exactly for `m <= 176` (`mobius/lab/rs/mertens-numerology`). □

**Fact 6.2 (below the proved wall).** On the float scan of Fact 8.1 of the large-base paper, the chord form of (W) holds at every base from `39363` at every missing digit, and from `28352` at the end digits `e_0 in {0, b - 1}`; no base of that scan is interval-certified, so on `39363 <= b < 92317` Theorem 1.1 holds on the scan and is Verified, never Proved, with `delta_0 >= 6.741 * 10^-8` at `39363`, `beta_0 = 1/2` (Fact 6.1). Generators: `mobius/lab/py/mobius-dissection`, verb `wall`, and `lab/py/digit-sum-twist`, verb `lemma`.

## Digit sums of primes

**Proof of Corollary 1.3.** Orthogonality modulo `q` writes the left side of the first display as `q^(-1) sum_(r mod q) e(-r a/q) psi_F(x; r/q)`, with `psi_F(x; alpha) = sum_(n in S_F, n <= x) Lambda(n) e(alpha s_b(n))`. The term `r = 0` is `psi_F(x)/q`. At `r != 0`, `(b - 1) r/q` is not an integer, because `gcd(q, b - 1) = 1` and `q` does not divide `r`; it lies in `q^(-1) Z`, so `norm((b - 1) r/q) >= 1/q`, and Theorem 1.1 at `beta_0 = 1/q` bounds each of these `q - 1` terms by `C A_F(x) x^(-delta)`. That is the first display.

For the prime count, the prime powers `p^j <= x` with `j >= 2` weigh at most `x^(1/2) (log x)^2` in any class, so the sum of `log p` over the primes of the class is `psi_F(x)/q + O(A_F(x) x^(-delta))` too, because `A_F(x) x^(-delta) >> x^(log(fill)/log(b) - delta)` and `log(fill)/log(b) > 1/2 + delta`: the dimension `log(fill)/log(b)` is at least `1 - alpha_1 > 4/5`, since `b <= B_b(F) <= fill b^(alpha_1)` by Lemma 2.3(i) and (ii), and `delta < delta_0 < 1/20`. Partial summation then gives `#{p in S_F : p <= x, s_b(p) = a mod q} - pi_F(x)/q << x^(log(fill)/log(b) - delta)`, since the error `E(t)` of the class sum is `<< A_F(t) t^(-delta)` and `A_F(t) <= fill t^(log(fill)/log(b))`, so `int_2^x |E(t)| dt/(t log^2 t)` is of that order. It remains to see that `pi_F(x) >> A_F(x)/log x`. Theorem 7.1(ii) of the half-interval paper proves, for any digit set with two consecutive digits and a shifted-grid certificate below `1/5`, which (W) supplies with `C_F = 1` by Lemma 2.3(ii) and the peel at `alpha = 0`,

```
psi_F(x) = kappa_F A_F(x) + O(A_F(x) exp(-c' sqrt(log x))) ,     kappa_F = (b/phi(b)) #{f in F : gcd(f, b) = 1}/fill ,
```

with `c' > 0` effective. And `kappa_F > 0` under (W): at `m = 1` the base has `phi(b) >= 2` digits coprime to it and one missing digit leaves one; at `m >= 2`, (W) forces `b > 32` and `m < b^(2/5)` (Lemma 2.3(iii)), while `phi(b) >= sqrt(b/2) > b^(2/5)` at `b > 32`, so some digit coprime to `b` stays in `F`. The bound `phi(n) >= sqrt(n/2)` holds because `phi(n)/sqrt(n)` is the product over `p^j` exactly dividing `n` of `p^(j/2 - 1)(p - 1)`, each factor at least `1` except `1/sqrt 2` at `p^j = 2`. □

The proof needs (W) itself, the inequality on the per-digit constants `P_b(m)` and `P'_b(e_0)`, and not only the one-step bound `B_b(F) <= fill b^(alpha_1)` that the proof of Theorem 1.1 reads, at one point: `kappa_F > 0`, which at `m >= 2` comes from `m < b^(2/5)`.

**Example 7.1 (the bare one-step bound does not suffice).** Let `b = 6^15` and let `F` be the digits `d < b` with `gcd(d, 6) > 1`, so `fill = b - phi(b) = 2b/3`. Then `F` holds the consecutive digits `2` and `3` and `B_b(F) < 90 fill < fill b^(1/5)`, so the proof of Theorem 1.1 runs and its conclusion holds; but no digit of `F` is coprime to `b`, every member of `S_F` shares a factor with `b`, the set holds no prime beyond `3`, and Corollary 1.3 has no main term. (W) fails: `F` misses `m = b/3` digits and `P_b(m) >= sqrt(m)`.

**Proof.** `1_F(d) = 1_(2 | d) + 1_(3 | d) - 1_(6 | d)`, so `g_F(v) = D_(b/2)(2v) + D_(b/3)(3v) - D_(b/6)(6v)` with `D_N(v) = sum_(n < N) e(n v)`. For `e` in `{2, 3, 6}` and `N = b/e`, `D_N` has period `1`, so `sum_(r mod b) |D_N(e (t + r)/b)| = e sum_(r mod N) |D_N((t + r)/N)| <= e Phi_N` by the kernel bound in the proof of Lemma 2.5 of the large-base paper, and `Phi_N <= N ((2/pi) log N + 2)` at `N >= 10`, from `H(n) <= log(N/2) + gamma_E + 1/(N - 2)`. So `B_b(F) <= 3 b ((2/pi) log b + 2)` and `B_b(F)/fill <= (9/pi) log b + 9 < 3 log b + 9 < 90`, since `log b = 15 log 6 < 27`, while `b^(1/5) = 6^3 = 216`. A member `n` of `S_F` is congruent modulo `b` to its last digit, which shares a factor with `6`, so `n` does. □

**Remark 7.2 (the other moduli).** When `g = gcd(q, b - 1) > 1`, the class of `s_b(p)` modulo `g` is the class of `p` modulo `g`, by Proposition 3.3, so the count asks for the primes of `S_F` in progressions modulo `g`: the lattice case, where the modulated transform keeps its full mass. Mauduit and Rivat carry that case for all primes in their Theorem 3, with main term `(gcd(m, q - 1)/m) pi(x; a, gcd(m, q - 1))` in their notation; on `S_F` the section Progressions and every modulus of [the digit-sum twist note](../notes/progressions-every-modulus.md) carries it at `exp(-c sqrt(log x))`.

## The falsification

Every reading below sits at base `10` or `37`, far below the wall of (W), so it tests the shape of the statements and not the theorem. Each check raises if it fails.

**Fact 8.1 (the supremum, read).** At base `10` missing `7` the supremum `sup_t |hat w_k(t)|/fill^k` is bracketed at `k = 1..8` and sixteen classes of `alpha`, five on the lattice and eleven off it: above by a max-product pass over the windows of six base-`10` digits of `t`, each window's factor bounded by its centre value plus the Lipschitz margin `pi h sum(F)` on a cell of width `h`, a float reading and not a certificate; below by the value attained at the best path refined on a local grid. At every off-lattice class and every `k >= 2` the upper bracket sits below the proved `(1 - c)^(floor(k/2))`, and at `k = 8`:

| `alpha` | `beta` | attained | upper | proved |
| --- | ---: | ---: | ---: | ---: |
| `1/2` | `0.5` | `0.071152` | `0.071165` | `0.521973` |
| `1/7` | `0.2857` | `0.310870` | `0.310915` | `0.826973` |
| `1/4` | `0.25` | `0.411592` | `0.411649` | `0.866052` |
| `1/10` | `0.1` | `0.870069` | `0.870167` | `0.977921` |
| `1/90` | `0.1` | `0.870069` | `0.870172` | `0.977921` |
| `1/900` | `0.01` | `0.998609` | `0.998716` | `0.999778` |

At the lattice classes `0`, `1/9`, `1/3`, `2/3` and `4/9` the value at `t = -alpha` is `1` to twelve digits at `k = 8`, as Proposition 3.3 demands; `1/3` is on the lattice at base `10`, `9/3 = 3`, and does not decay. The measured rate is far better than the proved one: the upper bracket's eighth root is `0.718676` at `alpha = 1/2` against the proved `(1 - c)^(1/2) = 0.921947`, since the pairing books one deficit per two positions where the readings show about one per position. At base `37` missing `18`, `k = 1..4`, the same holds at its nine off-lattice classes, and the classes `1/2` and `1/4`, on the lattice there because `36` is divisible by `4`, read `1`. Generator: `lab/py/digit-sum-twist`, verb `sup`, `10` seconds.

**Fact 8.2 (the peel).** At base `10` missing `7` the grid maximum of `sum_(r mod b) |g_F((T + r)/b)|` over `2 * 10^5` shifts reads `24.702183`, so `B_10(F) <= 24.702781` with the Lipschitz margin, in floating point, against the floor `2(b - 1) = 18` and the step-3 bound `b P_b(1) = 39.6623`. Over `16` classes of `alpha`, `12` shifts and depths `1..5`, `960` checks, the one-step identity of Lemma 4.1 holds to `1.40 * 10^-11`, and `Sigma_i(s)/B^i` and `Sigma_i(s)/(B Sigma_(i-1)(b s))` stay below `0.999914`. The modulation raises the unshifted mass, `Sigma_5(0)^(1/5)` reading `21.7151` at `alpha = 1/2` and `22.3176` at `1/7` against `19.5656` unmodulated: it moves mass toward the shifted supremum, which (W) already pays. Generator: `lab/py/digit-sum-twist`, verb `peel`, `9` seconds.

**Fact 8.3 (the modulated meters below `10^7`).** On base `10` missing `7` and missing `0`, `|sum_(n in S_F, n <= x) mu(n) e(alpha s_b(n))|/A_F(x)` is at most `0.0005` at `x = 10^7` over all sixteen classes, and its running maximum at most `0.0006` of `A_F(10^7)`, with no visible difference between the lattice and the rest: at this height `mu` cancels at every `alpha`, and a power `x^(-delta)` with `delta` near `10^-8` is invisible. The von Mangoldt meter separates them. On the lattice `psi_F(10^7; alpha)/psi_F(10^7)` sits within `0.0013` of the value the primes carry, `-0.5000` at `alpha = 1/3` and `2/3`, where the modulation is `e(p/3)` and `mu(3)/phi(3) = -1/2`, and near `0` at `1/9` and `4/9`, where `mu(9) = 0`. Off the lattice it sits within `0.0076` of the average of `e(alpha s_b(n))` over the members of `S_F` coprime to `b (b - 1) = 90`, which is near `0` at most classes and still large at `1/27`, `1/90` and `1/900`, `0.17`, `0.86` and `0.998` in modulus on missing `7`: those classes carry a small frequency `alpha`, and at `10^7` the frequency-zero factor `(|g_F(alpha)|/fill)^7` has barely begun to decay. Generator: `lab/py/digit-sum-twist`, verb `meter`.

**Fact 8.4 (the digit sums of the primes below `10^7`).** The `266823` primes of base `10` missing `7` below `10^7` deviate from `1/q` in their digit-sum classes by at most `0.00825`, `0.00880` and `0.00503` of the class size at `q = 2, 4, 5`, and by `0.05348`, `0.10263` and `0.24272` at `q = 7, 8, 10`. The slow moduli are the design's and not the primes': the members of `S_F` coprime to `90`, which carry the structure of the primes modulo `3` because `s_b(n) = n mod 3`, deviate by `0.04846`, `0.09911` and `0.23939` on the same moduli, and the primes sit within `0.02213` of that model at every `q` coprime to `9` tested. The `397866` primes of base `10` missing `0` read the same, `0.08164`, `0.14497` and `0.29072` at `q = 7, 8, 10` against the model's `0.07811`, `0.13695` and `0.29155`, within `0.02584` of the model. The slow class at `q = 10` is `r = 3`, at `alpha = 3/10 = 1/3 - 1/30`: far from the lattice in `beta = 0.3` but a small frequency `-1/30` away from the lattice point `1/3`, so it mixes only at the rate `|g_F(1/30)|/fill` per digit. At `q = 3`, which divides `b - 1`, the class `0` holds one prime, `3` itself, the lattice case of Remark 7.2. Generator: `lab/py/digit-sum-twist`, verb `gelfond`.

## What is in print

Mauduit and Rivat (2010), read in the published text, prove in Theorem 1 that `sum_(n <= x) Lambda(n) e(alpha s_q(n)) = O(x^(1 - sigma_q(alpha)))` with `sigma_q(alpha) > 0` whenever `(q - 1) alpha` is not an integer, and in Theorem 3 the count of primes with `s_q(p) = a mod m`, with main term `(gcd(m, q - 1)/m) pi(x; a, gcd(m, q - 1))`. The digits are unrestricted, `mu` is in none of their theorems, and their proof runs through Type I and Type II sums with a carry truncation, not through major arcs. Martin, Mauduit and Rivat (2014), read at source in the authors' copy, extend Theorem 1 to `e(f(n) + beta n)` for every digital function `f(n) = sum_k alpha_k |n|_k` with real `alpha_k`, with the explicit exponent `c_q norm((q - 1) theta)^2` when the `alpha_k` form a progression, all in their notation, which they credit to Drmota, Mauduit and Rivat, and remark that `mu` follows by the same method. A real `f` makes `e(f(n))` unimodular, so the indicator of a missing-digit set is outside that class.

The survey of Maynard (2019), read at source, proves a weak form of the Mauduit-Rivat theorem at large base by the circle method. Its Lemma 8.4 (2) is the two-position lemma for the full digit set, `|hat g(theta)| << b^(k (1 - delta_(alpha, b)))` with `delta_(alpha, b) = norm((b - 1) alpha)^2/(4 b^4)`, by the same locating step `norm(b^i theta + alpha) + norm(b^(i+1) theta + alpha) >= norm((b - 1) alpha)/b`. Its Lemma 10.1, as printed, bounds by that supremum times their count the major arcs whose denominator has a prime outside `b`, while its proof of Theorem 2.4 applies the same supremum times count to every small denominator: the mechanism of region C here. It treats the missing-digit set and the digit sum as two separate theorems and never their product. So the shape of this paper's argument is in print. What is written here is its object, the modulation carried onto a missing-digit set, for `mu` and `Lambda`, with an effective constant, from an explicit base, and with a two-position constant that does not decay with `b` (Fact 3.5), where the survey's per-digit loss `1 - b^(-delta_(alpha, b))`, about `beta^2 log b/(4 b^4)`, does.

On the missing-digit side, Maynard (2019) and Maynard (2022) carry the primes of a missing-digit set and no digit-sum modulation; the dissection of the large-base paper came from them.

Aloui (2015), read at the publisher page in its abstract only; its definition of `T_N(alpha, beta)` and its Theorems 2.1, 2.5 and 2.6 are known here through a machine summary of that page and through the citation in Saavedra-Araya, the full text unread, carries Gelfond's theorem onto the integers of a missing-digit set, jointly with residue classes, and onto their `z`-free members, with power-saving error terms, through the sums `sum_(n in S_F, n < N) e(n theta + alpha s_b(n))`, a sum of the shape of the modulated transform of Definition 2.1; no prime and no Mobius sum enters. It is the integer neighbour of Corollary 1.3. Saavedra-Araya (2024), read at source in its abstract and Theorem A, gives the necessary and sufficient form for the integers: for a missing-digit set and `gcd(g, a) = gcd(a, a') = 1`, the classes of `n mod a` and `s_g(n) mod a'` are jointly uniform exactly when `gcd(a a', d_2 - d_1, ..., d_t - d_1) = 1`, by Markov chains, with no primes. Aloui, Mauduit and Mkaouar (2017), read in its abstract, carry the digit-sum classes onto palindromes with missing digits, which are not primes. A search of the literature on the digit sum of primes with missing, restricted or ellipsephic digits, of the publication lists of Martin and of Rivat, of the reference lists of Mauduit and Rivat (2010) and of the survey, and of the prior-work list of Saavedra-Araya (2024) found no theorem on the digit sum of primes, or of `mu`, restricted to a missing-digit set. As far as that search reaches, no digit-sum theorem for the primes of a missing-digit set is in print, and the shape of the argument is.

## Open problems

Five things are left undone. The pairing sets the exponent past about `10^5`: Lemma 3.1 books one deficit per two positions, while Fact 8.1 reads close to one per position, `0.718676` per digit at `alpha = 1/2` against the proved `0.921947`, and a deficit per position would raise `sigma`, which binds from `100084` at `beta_0 = 1/2`. The base is that of (W), `92317` by proof and `39363` on the scan, against `584` and `115` for the unmodulated bound; a modulated chain would need, in place of the Lebesgue sums of `D_(b^l)`, the grid `l^1` norms of the full-digit modulated kernel `prod_(j < l) D_b(b^j t + alpha)`, the discrete transform Mauduit and Rivat estimate, and a modulated window would need a supremum over shift vectors, and neither is written. The lattice is carried only at `exp(-c sqrt(log x))`, in the section Progressions and every modulus of [the digit-sum twist note](../notes/progressions-every-modulus.md): `mu` in progressions modulo `(b - 1) b^j` on the set gives that bound at every lattice `alpha` and Corollary 1.3 at every `q`, the missing-digit form of Mauduit-Rivat Theorem 3, and a power saving there would ask a zero-free strip for the characters modulo `(b - 1) b^j`. In the unmodulated dissection the region C1 also caps the saving at `exp(-c sqrt(log x))`, through the rate of its contraction against the `3 Z^2` residues of the region and with no L-function, so the root-log ceiling there has two sources, and whether C1 can be priced at a power is open. And the exponent stays below `log(fill)/(20 log b)` at every base, nothing here touches sets of bounded fill, and the true size of the modulated sums is not examined.

## Reproducibility

One study prints every number of Sections 3, 6 and 8, run from the repository root with one verb and raising if any check fails: [lab/py/digit-sum-twist](../lab/py/digit-sum-twist/). `uv run python research/gelfond/lab/py/digit-sum-twist/twist.py lemma`, under a second, prints Fact 3.5, Example 3.4, Fact 6.1 with the bisection of the binding switch on `92317..10^6`, the margin at `92316`, and the interval enclosures of the wall used in Corollary 1.2; `sup`, in `10` seconds, prints Fact 8.1 over `k = 1..8` at base `10` and `k = 1..4` at base `37`; `peel`, in `9` seconds, prints Fact 8.2; `meter gelfond`, in `7` seconds with a peak of `1.7` GB, prints Facts 8.3 and 8.4 by a sieve to `10^7`. The large-base constants are computed in `mpmath` at `40` digits and floored at four significant digits, so every printed `c`, `gamma`, `sigma`, `delta_0` and `1/5 - alpha_1` is a lower bound; the base-`10` constants are floats; the upper bracket of Fact 8.1 is a float reading with an explicit Lipschitz margin and never a certificate. The chord walls of Fact 6.2 and the budget of Corollary 1.2 are printed by the studies the large-base paper names, `mobius/lab/py/mobius-dissection`, verb `wall`, and `mobius/lab/rs/mertens-numerology`. Theorem 1.1 and Corollary 1.3 compute no constant beyond these. The figure is `figures/paper-digit-sum-twist-power-saving.ts`, drawn in the browser in under a second a theme; it evaluates `g_F` in closed form and asserts the full mass at the centre, a lattice segment through it, the counts of `11` and `10` segments, and at `2 * 10^5` points of the line at `alpha = 1/2` a phase at least `1/22` from the integers and a share at most `1 - c` with `c >= 0.15`.

## References

- Mauduit and Rivat 2010, Sur un probleme de Gelfond: la somme des chiffres des nombres premiers, Ann. of Math. 171(3), 1591-1646. [doi.org/10.4007/annals.2010.171.1591](https://doi.org/10.4007/annals.2010.171.1591)
- Martin, Mauduit and Rivat 2014, Theoreme des nombres premiers pour les fonctions digitales, Acta Arith. 165(1), 11-45. [doi.org/10.4064/aa165-1-2](https://doi.org/10.4064/aa165-1-2)
- Maynard 2019 survey, Digits of primes, read at source. [arxiv.org/abs/1910.13402](https://arxiv.org/abs/1910.13402)
- Aloui 2015, Sur les entiers ellipsephiques: somme des chiffres et repartition dans les classes de congruence, Period. Math. Hungar. 70, 171-208, abstract read at the publisher page, full text unread. [doi.org/10.1007/s10998-014-0066-8](https://doi.org/10.1007/s10998-014-0066-8)
- Saavedra-Araya 2024, Distribution of integers with digit restrictions via Markov chains. [arxiv.org/abs/2411.07418](https://arxiv.org/abs/2411.07418)
- Aloui, Mauduit and Mkaouar 2017, Somme des chiffres et repartition dans les classes de congruence pour les palindromes ellipsephiques, Acta Math. Hungar. 151, 409-455, read in its abstract. [doi.org/10.1007/s10474-017-0688-4](https://doi.org/10.1007/s10474-017-0688-4)
- Basak, Robles and Zaharescu 2023, Exponential sums over Mobius convolutions with applications to partitions. [arxiv.org/abs/2312.17435](https://arxiv.org/abs/2312.17435)
- Maynard 2022, Primes and polynomials with restricted digits, Int. Math. Res. Not. 2022, 10626-10648, read in arXiv:1510.07711v1. [doi.org/10.1093/imrn/rnab002](https://doi.org/10.1093/imrn/rnab002)
- Maynard 2019, Primes with restricted digits, Invent. Math. 217, 127-218. [link.springer.com](https://link.springer.com/article/10.1007/s00222-019-00865-6)
