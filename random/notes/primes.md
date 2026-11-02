---
title: Along the primes
lead: A digit-parity design read along the primes is biased by its major-arc sum, to the noise floor. Replace the von Mangoldt function by its Ramanujan expansion truncated at q <= 1024, sum the design's twisted sums exactly, and the result meets a census of the primes below 2^32 within 2.35 sigma for all 55 distinct designs of base at most 6, 47 independent streams, at each of 49 checkpoints from 2^20, sigma^2 the sum of Lambda^2. The pooled rms is 0.715 against a random-sign floor of 0.710; at 2^32 alone it is 0.936 against 0.781, an excess that higher q do not remove. For Thue-Morse the leading part is Shevelev's closed form -3^floor(k/2) at 2^k, and the full sum meets the census within 0.54 sigma at every power of 2 from 2^20, where Q <= sqrt x, to 2^32, where the bias reaches -146.8 sigma; the q = 3 term alone misses by 31.7 sigma at 2^31, and the integers prime to 6 miss by 4.1 sigma at 2^32. Conjecturally the Gelfond sum is -G(log_2 x) x^(log_4 3)/log x, with G a limit of period 2, 1 at even log_2 x and down to 0.4834. The exponent, the closed form and the model by the integers prime to 6 are Shevelev's; the terms past q = 6, the test at the noise floor and the law across designs were not found in the sources read.
figure: research-random-primes
---

[Along the squares](squares.md) read every design `f_{b,D}`, `2 <= b <= 6`, along the first `2^26` primes and found it biased far beyond the noise. Thue-Morse `t(n) = (-1)^(s_2(n))` read `z_B = -89.7`, that is, a mean of about `-0.011` over `2^26` primes. At base 6 the units decided: `6:1` and `6:5` read about 0 while `6:2`, `6:3` and `6:4` read `63.6`, `95.0` and `152.6`. And `6:23`, noise along the squares, was slow along the primes. Mauduit and Rivat 2010 prove that the Thue-Morse bias is `o(pi(x))` with a power saving ([refs](../refs.md)); they do not say how large it is.

This chapter gives the size. **The bias of a design along the primes is its major-arc sum, to the noise floor.** Replace `Lambda(n)` by its Ramanujan expansion `sum_q (mu(q)/phi(q)) c_q(n)` truncated at `q <= Q`. Then `sum_(n<x) x(n) Lambda(n)` becomes a sum of the design's twisted sums at the rationals `a/q`, each an exact finite product. With `Q = 1024` this sum meets a census of the primes below `2^32` for every design of base at most 6, at every checkpoint, within the scatter that random signs would show. For Thue-Morse the leading part has Shevelev's closed form, `-3^floor(k/2)` at `x = 2^k`, with a coefficient that is periodic in `log_4 x` in the limit. That is the true order the heuristic predicts for Gelfond's sum. The heuristic is the circle method's, and its leading part for Thue-Morse, the integers prime to 6 with the exponent `log 3/log 4`, is Shevelev's 2007 model. What the census adds is the match at the noise floor for every design and the terms past the leading one.

The generator is the study `lab/rs/primes`, verbs `census`, `thue`, `law`, `ladder` and `exact`. The census sieves to `2^32` in seconds on 8 threads and reproduces the squares chapter's numbers to the last digit, `-89.7` for Thue-Morse along the first `2^26` primes and `-46.7` for `6:23`. The figure plots `-sum_(n<x) t(n) Lambda(n)/x^(log_4 3)` against `log_2 x` from 14 to 22, the primes as dots. The blue curve is the `q | 30` part of the major-arc sum, which the dots settle onto as the noise falls like `x^(-0.29)`; the grey curve is the `q = 3` term alone, which the dots leave at every odd power of 2. The figure checks that `pi(2^22) = 295947` and that the `q | 6` part and `T_3` at `2^22` are both `-3^11`.

## The objects

A design `b:D` has the sign `x(n) = prod_i v(n_i)` over the base-`b` digits of `n`, with `v(d) = (-1)^[d in D]` and `v(0) = 1`, as in [the squares](squares.md). Beside the 57 designs run Rudin-Shapiro and a coin, the top bit of splitmix64 at `n`. Two sums measure the bias along the primes:

- the weighted sum `Psi(x) = sum_(n<x) x(n) Lambda(n)`, in units of `sigma(x) = sqrt(sum_(n<x) Lambda(n)^2)`, which is its standard deviation for random signs;
- the plain sum `sum_(p<x) x(p)`, in units of `sqrt(pi(x))`.

The weighted sum is the clean object, and the plain one follows from it by partial summation. Below `2^32`, `sigma = 301611.6` and `sqrt(pi) = 14257.6`.

The model is the major-arc sum. Write `e(t) = exp(2 pi i t)`, `S(alpha; x) = sum_(n<x) x(n) e(n alpha)` for the twisted sum, and `c_q(n) = sum_((a,q)=1) e(a n/q)` for Ramanujan's sum. Then

`M_Q(x) = sum_(q<=Q) T_q(x)`, `T_q(x) = (mu(q)/phi(q)) sum_((a,q)=1) S(a/q; x) = (mu(q)/phi(q)) sum_(n<x) x(n) c_q(n)`.

It is `sum_(n<x) x(n) Lambda_Q(n)` with `Lambda_Q(n) = sum_(q<=Q) (mu(q)/phi(q)) c_q(n)`, the truncated Ramanujan expansion of the von Mangoldt function. In the language of the circle method, the primes carry the mass `mu(q)/phi(q)` at every reduced `a/q`, the singular-series weight, and `M_Q` collects what the design's spectrum meets there. Only squarefree `q` count, and `T_q` sums over the units mod `q`, so the expansion knows that a prime avoids `0 mod q`. The term `q = 1` is exact: it is `sum_(n<x) x(n)`, the integers at density one. Here `Q = 1024` is fixed and the checkpoints start at `2^20`, so that `Q <= sqrt(x)` throughout.

Every term is a finite product. **For a design of base `b`, `S(alpha; b y) = P(e(alpha)) S(b alpha; y)` for every integer `y >= 0`, `P(z) = sum_(d<b) v(d) z^d`, so `S(alpha; b^k) = prod_(j<k) P(e(b^j alpha))`; if `b^p alpha = alpha (mod 1)`, then `S(alpha; b^p y) = S(alpha; b^p) S(alpha; y)` at every `y`. Proved.** The integers below `b y` are `n = b m + d` with `m < y` and `d < b`, and `x(b m + d) = x(m) v(d)` because the digits of `b m + d` are `d` followed by those of `m`; the sum factors. Iterating `p` times around the orbit gives the second identity. So along a periodic orbit `abs(S(alpha; x))` is exactly `x^g` times a function periodic in `log_(b^p) x`, where `g = (1/p) sum_(j<p) log_b abs(P(e(b^j alpha)))` is the orbit's growth exponent. The twisted sum itself turns by the phase of `S(alpha; b^p)` every `p` levels, so a term is periodic only when that phase is a root of unity. The orbit of `1/7` under doubling has product `-i sqrt 7`, and Thue-Morse's `T_7` at `8^m`, `m = 1..6`, runs `0, 7, 0, -49, 0, 343`. The growth exponent is `1 - eta(alpha)` in the squares chapter's notation, and a term eventually rises above the noise exactly when `g > 1/2`. Each `S(a/q; x)` is summed by a digit walk: the product above gives `S(a/q; b^j)`, and the digits of `x`, read from the top, add the partial blocks. Rudin-Shapiro, which is not a product, walks its two-state automaton the same way.

The noise has a scale the census can predict. For signs drawn at random, `Psi - M_Q = sum s(n) (Lambda(n) - Lambda_Q(n))` has variance `sum (Lambda - Lambda_Q)^2`. Orthogonality of Ramanujan's sums gives `sum_(n<x) Lambda_Q(n)^2 ~ x H`, and the prime number theorem in progressions gives `sum_(n<x) Lambda Lambda_Q ~ theta(x) H`, with `H = sum_(q<=Q) mu^2(q)/phi(q) = 8.266` at `Q = 1024`. So random signs leave a residual of sd `sqrt(1 - theta H/sum Lambda^2)` in units of `sigma`, read off the census at each checkpoint: `0.598` at `2^20` and `0.781` at `2^32`, and `0.710` pooled over the 49 quarter-octave checkpoints from `2^20`. This is the floor. A residual of a few units of `sigma` is noise; a residual of tens is a missing term.

## Thue-Morse to 2^32

At base 2 the tile polynomial is `P(z) = 1 - z`, and the twisted sums are the products `prod_(j<k) (1 - e(2^j alpha))`. The point `1/3` sits on the 2-cycle `1/3, 2/3`, where `(1 - omega)(1 - omega^2) = 3`, `omega = e(1/3)`. That is the spectral peak of Thue-Morse, `g = log_4 3 = 0.7925`, the squares chapter's `eta = 0.2075` at the orbit `1/3`. Three families of terms carry the bias: the denominators dividing 6, those dividing 30 with a factor 5, and the rest.

**For Thue-Morse at `x = 2^k`, `k >= 2`: the `q | 6` part `T_1 + T_2 + T_3 + T_6 = (6/phi(6)) sum_(n<x, (n,6)=1) t(n)` equals `-3^floor(k/2)`; `T_5 + T_10 = -5^floor(k/4) (1, 1/2, 0, 5/2)[k mod 4]`; `T_15` and `T_30` are bounded. For every integer `y >= 0`, `S(1/3; 4y) = 3 S(1/3; y)`. Proved.** The last identity is the product formula with `P(omega) P(omega^2) = 3`. For the first, `T_1 = sum_(n<2^k) t(n) = 0` and `T_2 = -S(1/2; 2^k) = -2 S(0; 2^(k-1)) = 0`. Then `T_3 = -Re S(1/3; 2^k)` and `T_6 = Re S(1/6; 2^k)`, with `S(1/3; 4^m) = 3^m`, `S(1/3; 2 . 4^m) = (1 - omega) 3^m` and `S(1/6; 2^k) = e(-1/6) S(1/3; 2^(k-1))`, since `1 - e(1/6) = e(-1/6)`. At even `k = 2m` this gives `T_3 = -3^m` and `T_6 = 0`; at odd `k = 2m + 1`, `T_3 = -(3/2) 3^m` and `T_6 = 3^m/2`, so the part is `-3^m` either way. This is Shevelev's 2007 Theorem 2(1), `Delta^odd_3([0, 2^n)) = 3^(floor(n/2) - 1)` for the evil-minus-odious count of the odd multiples of 3, since the `q | 6` part at `2^k` is `-3 Delta^odd_3([0, 2^k))`. For the fifths write `k = 4m + r` and `zeta = e(1/5)`. The orbit of `1/5` has period 4 and product `Phi_5(1) = 5`, so `T_5 = -(5^m/4) s_r` with `s_r = sum_(a=1..4) prod_(j<r) (1 - zeta^(a 2^j))`, which is `4, 5, 5, 10` for `r = 0, 1, 2, 3`; the last because the product is `5/(1 - zeta^(3a))` and `sum_c 1/(1 - zeta^c) = 2`. With `w` over the primitive tenth roots, `S(a/10; 2^k) = (1 - w) S(a/5; 2^(k-1))` gives `T_10 = (5^m/4) sum_w (1 - w) prod_(j<r-1) (1 - w^(2^(j+1)))` for `r >= 1`, the sums being `3, 5, 0` for `r = 1, 2, 3`; at `r = 0` it is `(5^m/4) sum_w (1 - w)/(1 + w) = 0`, since `w^6 = -w` and `(1 - w)/(1 + w)` is imaginary on the unit circle. So `T_5 + T_10` is `-5^m`, `-(5/4) 5^m + (3/4) 5^m = -5^m/2`, `-(5/4) 5^m + (5/4) 5^m = 0` and `-(5/2) 5^m + 0` for `r = 0, 1, 2, 3`. The two orbits mod 15 are negatives of each other and their products multiply to `Phi_15(1) = 1`, so each has modulus 1. The `q | 6` part is the multiples of 3 seen from the units: since `t(2m + 1) = -t(m)` and `t(2m) = t(m)`, it equals `3 (N(x/2) - N(x)) + O(1)` with Newman's sum `N(y) = sum_(n<y, 3 | n) t(n)`, whose excess of evil multiples of 3 is the Moser-Newman phenomenon, proved by Newman 1969 and given an exact formula by Coquet 1983 ([refs](../refs.md)). The primes avoid the multiples of 3, and they inherit the opposite excess. The generator checks both closed forms at `k = 12..32`, to `4.5e-15` and `8.2e-15`, and its test counts the first by brute force to `2^20`. Over the same range it prints `abs(T_15 + T_30) <= 5.8e-3 sigma` and at most `0.5215 sigma` for the eight terms `q | 210` with `7 | q`: the orbit of `1/7` has product `(1 - zeta_7)(1 - zeta_7^2)(1 - zeta_7^4) = -i sqrt 7` over its period 3, so `g = log_2 7/6 = 0.468`, under the noise.

The census reads `Psi` at every power of 2. The table gives it with the terms in units of `sigma`, the rest of `M_Q`, the residual `r = (Psi - X)/sigma` for four truncations `X`, the normalised coefficient `G = -Psi/x^(log_4 3)` observed and modelled, and the residual of the plain sum `sum_(p<x) t(p)` against the model weighted by `1/log n`, in units of `sqrt(pi)`.

| `k` | `Psi` | `z` | `T_3` | `T_6` | `T_5` | `T_10` | rest | `r` for `T_3` | for `q \| 6` | for `q \| 30` | for `M_Q` | `G` observed | `G` model | plain `r` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 | -868 | -5.03 | -4.22 | -0.00 | -0.72 | -0.00 | -0.09 | -0.80 | -0.80 | -0.07 | 0.01 | 1.1905 | 1.1934 | 1.01 |
| 13 | -1049 | -4.11 | -4.29 | 1.43 | -0.61 | 0.37 | -1.04 | 0.17 | -1.26 | -1.01 | 0.02 | 0.8309 | 0.8357 | 0.77 |
| 14 | -2376 | -6.30 | -5.80 | -0.00 | -0.41 | 0.41 | -0.45 | -0.50 | -0.50 | -0.50 | -0.05 | 1.0863 | 1.0769 | 0.55 |
| 15 | -3006 | -5.43 | -5.92 | 1.97 | -0.56 | -0.00 | -0.84 | 0.50 | -1.48 | -0.92 | -0.08 | 0.7936 | 0.7815 | 0.42 |
| 16 | -8322 | -10.25 | -8.08 | -0.00 | -0.77 | -0.00 | -1.20 | -2.17 | -2.17 | -1.40 | -0.20 | 1.2683 | 1.2435 | 0.17 |
| 17 | -7073 | -5.95 | -8.28 | 2.76 | -0.66 | 0.39 | -0.06 | 2.33 | -0.43 | -0.17 | -0.11 | 0.6224 | 0.6110 | 0.14 |
| 18 | -19659 | -11.34 | -11.36 | -0.00 | -0.45 | 0.45 | -0.78 | 0.01 | 0.01 | 0.01 | 0.80 | 0.9988 | 1.0689 | 0.94 |
| 19 | -21185 | -8.39 | -11.69 | 3.90 | -0.62 | -0.00 | -0.39 | 3.30 | -0.59 | 0.02 | 0.41 | 0.6214 | 0.6521 | 0.58 |
| 20 | -64566 | -17.59 | -16.08 | -0.00 | -0.85 | -0.00 | -0.42 | -1.50 | -1.50 | -0.65 | -0.23 | 1.0934 | 1.0790 | -0.09 |
| 21 | -67263 | -12.62 | -16.61 | 5.54 | -0.73 | 0.44 | -1.13 | 4.00 | -1.54 | -1.25 | -0.12 | 0.6577 | 0.6514 | 0.00 |
| 22 | -191212 | -24.74 | -22.92 | -0.00 | -0.51 | 0.51 | -1.76 | -1.82 | -1.82 | -1.82 | -0.06 | 1.0794 | 1.0768 | 0.02 |
| 23 | -214276 | -19.14 | -23.74 | 7.91 | -0.70 | -0.00 | -2.33 | 4.60 | -3.32 | -2.62 | -0.29 | 0.6984 | 0.6879 | -0.19 |
| 24 | -569970 | -35.19 | -32.82 | -0.00 | -0.96 | -0.00 | -1.92 | -2.38 | -2.38 | -1.41 | 0.50 | 1.0725 | 1.0878 | 0.52 |
| 25 | -550530 | -23.52 | -34.06 | 11.35 | -0.83 | 0.50 | -1.02 | 10.54 | -0.82 | -0.48 | 0.54 | 0.5981 | 0.6118 | 0.56 |
| 26 | -1618514 | -47.89 | -47.17 | -0.00 | -0.58 | 0.58 | -0.62 | -0.72 | -0.72 | -0.72 | -0.10 | 1.0152 | 1.0131 | -0.06 |
| 27 | -1653767 | -33.92 | -49.05 | 16.35 | -0.80 | -0.00 | -0.48 | 15.13 | -1.22 | -0.42 | 0.06 | 0.5989 | 0.6000 | 0.11 |
| 28 | -4902049 | -69.74 | -68.04 | -0.00 | -1.11 | -0.00 | -0.88 | -1.69 | -1.69 | -0.58 | 0.30 | 1.0249 | 1.0293 | 0.32 |
| 29 | -4878842 | -48.18 | -70.85 | 23.62 | -0.96 | 0.58 | -0.33 | 22.67 | -0.95 | -0.56 | -0.23 | 0.5889 | 0.5861 | -0.15 |
| 30 | -14417833 | -98.90 | -98.42 | -0.00 | -0.67 | 0.67 | -0.85 | -0.47 | -0.47 | -0.47 | 0.38 | 1.0048 | 1.0086 | 0.36 |
| 31 | -14878589 | -70.93 | -102.61 | 34.20 | -0.93 | -0.00 | -1.46 | 31.68 | -2.53 | -1.59 | -0.14 | 0.5987 | 0.5975 | -0.06 |
| 32 | -44276255 | -146.80 | -142.72 | -0.00 | -1.30 | -0.00 | -2.51 | -4.08 | -4.08 | -2.78 | -0.27 | 1.0286 | 1.0267 | -0.25 |

**The major-arc sum `M_Q`, `Q = 1024`, meets `sum_(n<2^k) t(n) Lambda(n)` within `0.54 sigma` at every `k` from 20 to 32, where `Q <= sqrt(x)`, and the census reads `z = -146.80` against `-146.53` at `2^32`. Below `2^20` the same `Q` exceeds `sqrt(x)`, outside the rule, and the residual stays within `0.80 sigma` at `k = 12..19`; the plain sum stays within `1.01 sqrt(pi)` at every `k` from 12. Verified.** At `2^32` the coefficient reads `1.0286` against the model's `1.0267`, and at `2^31` `0.5987` against `0.5975`.

**The `q = 3` term alone is the main term of the Thue-Morse bias. Refuted.** It misses by `22.67 sigma` at `2^29` and `31.68 sigma` at `2^31`. At odd `k` the term `q = 6` cancels a third of `T_3`: `1/6` is a preimage of `1/3` under doubling, and `mu(6)/phi(6) = 1/2` does not vanish. At even `k` the term `q = 6` is zero, and the `q = 3` term still misses by `4.08 sigma` at `2^32`. The `q | 6` part, the integers prime to 6 reweighted by 3, misses there by the same `4.08 sigma`, and the `q | 30` part by `2.78 sigma`. What closes the gap is not one term but many. In units of `sigma` at `2^32`, after `T_3 = -142.72`, the largest are `T_17 = -2.09`, `T_257 = -1.39`, `T_5 = -1.30`, `T_31 = 1.20`, `T_65 = 1.13`, `T_255 = -1.07`, `T_127 = -1.01`, `T_43 = -0.87` and `T_341 = 0.72`. Their best orbits grow at `g = 0.633, 0.702, 0.580, 0.636, 0.675, 0.659, 0.686, 0.634` and `0.685`, all above `1/2`, each point of an orbit weighted by `1/phi(q)`. The denominators divide `2^p +- 1` for small `p`: the orbit of `3/17` has period 8, that of `5/31` period 5. No single power describes the bias past its leading term: orbits that follow the cycle of `1/3` for a long stretch grow at rates up to `log_4 3`, under weights that fall like `2^(-p)`.

The leading coefficient is the `q | 6` part `A(x)` divided by `x^(log_4 3)`. Only `T_3` scales exactly. **`A(4x) = 3 A(x) + O(1)`: for `x = 1..2^22`, `A(4x) - 3 A(x)` lies in `[-9, 9]` and is nonzero at 2796203 of the 4194304 values. Verified.** So `G_6(u) = lim_m -A(2^(u + 2m))/2^((u + 2m) log_4 3)` exists along each residue of `u` and has period 2. **Sampled at 2048 points of `u = log_2 x` in `[30, 32)`, where the `O(1)` costs under `1e-6`, `G_6` is `1` at even `u` and `1/sqrt 3 = 0.5774` at odd `u`, with mean `0.7481`; its least sample is `0.4835` at `u = 30.585`, and 4096 points within `2^-10` of it bring the minimum to `0.48346`, so `G_6` reaches below `0.4835` and the chapter writes its lower end as `0.4834`. Verified.** The minimum sits near `x = 3 . 2^(2m - 1)`.

**Conjecture: `sum_(n<x) t(n) Lambda(n) = -G_6(log_2 x) x^(log_4 3) (1 + o(1))` and `sum_(p<x) (-1)^(s_2(p)) = -G_6(log_2 x) x^(log_4 3)/log x (1 + o(1))`, with `G_6` the limit of period 2 defined by the `q | 6` part.** The plain form follows from the weighted one by partial summation, the correction being a relative `O(1/log x)`. The exponent `log 3/log 4` is Shevelev's 2007 Conjecture 2 for the race between odious and evil primes, the closed form at `2^k` is his Theorem 2(1), and his heuristic, three times the race among the odd multiples of 3 over `log n`, is the `q | 6` part here with Coquet's periodic function inside it. The terms past `q = 6` and a test at the noise floor were not found in the sources read. The model also fits the plain sums along the first `2^26` primes that the squares chapter measured: Thue-Morse `-89.7` against `-90.0`, `6:2` `63.6` against `63.8`, `6:3` `95.0` against `94.9`, `6:4` `152.6` against `153.2`, `6:23` `-46.7` against `-45.9`, `6:1` `0.2` against `-0.5` and `6:5` `-0.9` against `-0.3`. Since the model meets `z` itself there, the squares chapter's slope of `0.270` per doubling against `1/2 - eta = 0.292` is the factor `1/sqrt(log x)` in `z` and the period-2 wave of `G_6`, not a missing term.

## Every design

`law` runs the same model for every design and for Rudin-Shapiro, at `2^32` and at the 49 quarter-octave checkpoints from `2^20`. The table gives `z = Psi/sigma`, the model `M_Q/sigma`, the residual `r` at `2^32` and its largest absolute value over the checkpoints, the largest term `q:T_q/sigma`, and the part from `q | rad(b)`, the restriction to the residues prime to `b`. Its last two columns are the plain `sum_(p<x) x(p)/sqrt(pi)` and its model.

| design | `z` | model | `r` | max | top `q:T_q` | `q \| rad b` | plain `z` | model |
|---|---|---|---|---|---|---|---|---|
| `2:1` | -146.80 | -146.53 | -0.27 | 0.88 | `3:-142.72` | 0.00 | -146.91 | -146.66 |
| `3:1` | -14240.03 | -14240.06 | 0.04 | 0.13 | `2:-14240.06` | 0.00 | -14257.64 | -14257.91 |
| `3:2` | -0.23 | 0.10 | -0.33 | 1.42 | `10:0.58` | 0.00 | -0.12 | 0.14 |
| `3:12` | 0.23 | -0.37 | 0.60 | 1.93 | `5:-0.58` | -0.00 | 0.12 | -0.40 |
| `4:1` | -10.84 | -10.68 | -0.15 | 1.54 | `3:-7.86` | 0.00 | -10.34 | -10.20 |
| `4:2` | -17.24 | -17.91 | 0.67 | 1.62 | `3:-7.86` | 0.43 | -16.79 | -17.41 |
| `4:12` | -146.80 | -146.53 | -0.27 | 0.88 | `3:-142.72` | 0.00 | -146.91 | -146.66 |
| `4:3` | -8.47 | -7.69 | -0.78 | 2.02 | `5:-3.79` | 0.00 | -8.62 | -7.83 |
| `4:13` | -116.68 | -116.68 | -0.00 | 1.48 | `5:-107.43` | 0.00 | -116.43 | -116.52 |
| `4:23` | -115.09 | -115.45 | 0.35 | 2.35 | `5:-107.43` | 0.00 | -115.66 | -116.02 |
| `4:123` | -3.36 | -3.23 | -0.13 | 1.51 | `5:-3.79` | 0.43 | -3.43 | -3.25 |
| `5:1` | -4.86 | -3.84 | -1.02 | 2.16 | `2:-10.63` | 4.84 | -5.20 | -4.28 |
| `5:2` | -2.84 | -2.84 | -0.00 | 1.26 | `3:-12.01` | 6.80 | -2.31 | -2.32 |
| `5:12` | 1.87 | 2.17 | -0.30 | 1.33 | `3:-2.57` | 0.00 | 1.77 | 2.08 |
| `5:3` | 4.86 | 3.56 | 1.31 | 1.50 | `1:10.63` | 8.86 | 5.20 | 3.97 |
| `5:13` | -14240.03 | -14240.06 | 0.04 | 0.13 | `2:-14240.06` | 0.00 | -14257.64 | -14257.91 |
| `5:23` | -1.87 | -2.25 | 0.38 | 1.48 | `6:2.57` | 0.00 | -1.77 | -2.17 |
| `5:123` | 2.84 | 2.47 | 0.38 | 1.02 | `6:12.01` | 0.00 | 2.31 | 1.93 |
| `5:4` | 17.79 | 18.62 | -0.83 | 1.10 | `1:21.09` | 17.57 | 18.06 | 18.85 |
| `5:14` | -136.69 | -137.00 | 0.32 | 1.88 | `3:-133.34` | -0.00 | -137.31 | -137.63 |
| `5:24` | -11.62 | -11.55 | -0.07 | 1.38 | `2:-5.81` | 0.00 | -10.93 | -10.91 |
| `5:124` | -121.78 | -120.81 | -0.97 | 1.46 | `3:-141.87` | -0.00 | -123.93 | -122.94 |
| `5:34` | 136.69 | 136.51 | 0.18 | 1.14 | `6:133.34` | 0.00 | 137.31 | 137.10 |
| `5:134` | -17.79 | -19.74 | 1.96 | 1.96 | `2:-21.09` | -0.00 | -18.06 | -20.02 |
| `5:234` | 121.78 | 120.92 | 0.87 | 1.85 | `6:141.87` | -0.00 | 123.93 | 123.02 |
| `5:1234` | 11.62 | 11.27 | 0.35 | 1.26 | `1:5.81` | 9.68 | 10.93 | 10.60 |
| `6:1` | -0.16 | -1.68 | 1.52 | 1.56 | `1:3.29` | 0.00 | -0.12 | -1.65 |
| `6:2` | 163.83 | 163.59 | 0.23 | 1.33 | `1:108.02` | 162.03 | 164.54 | 164.37 |
| `6:12` | -10.94 | -11.81 | 0.87 | 1.42 | `5:-10.48` | 0.00 | -11.05 | -11.89 |
| `6:3` | 159.98 | 159.96 | 0.02 | 1.40 | `1:108.02` | 162.03 | 162.11 | 162.09 |
| `6:13` | -54.40 | -53.31 | -1.09 | 1.49 | `7:-29.33` | 0.00 | -54.08 | -53.00 |
| `6:23` | -31.23 | -31.32 | 0.09 | 1.74 | `5:-71.52` | 0.07 | -33.34 | -33.36 |
| `6:123` | -20.39 | -21.34 | 0.94 | 1.78 | `5:-36.02` | 0.00 | -22.28 | -23.23 |
| `6:4` | 161.41 | 161.87 | -0.47 | 1.29 | `1:108.04` | 162.07 | 164.71 | 165.20 |
| `6:14` | -124.90 | -123.22 | -1.68 | 1.76 | `5:-97.77` | -0.00 | -126.17 | -124.51 |
| `6:24` | -2.22 | -2.51 | 0.29 | 1.08 | `5:-4.22` | 0.07 | -3.20 | -3.48 |
| `6:124` | -12.41 | -10.50 | -1.91 | 1.96 | `5:-49.22` | -0.00 | -12.38 | -10.47 |
| `6:34` | -20.75 | -22.47 | 1.73 | 2.13 | `5:-15.81` | 0.07 | -20.29 | -22.04 |
| `6:134` | -108.00 | -108.79 | 0.79 | 1.50 | `5:-66.22` | -0.00 | -108.55 | -109.27 |
| `6:234` | -69.86 | -69.56 | -0.29 | 1.35 | `5:-55.62` | -0.00 | -70.74 | -70.43 |
| `6:1234` | 2.14 | 0.62 | 1.52 | 1.52 | `5:-2.16` | -0.00 | 1.80 | 0.30 |
| `6:5` | -1.56 | -0.84 | -0.72 | 1.41 | `1:108.05` | 0.00 | -1.59 | -0.84 |
| `6:15` | -6.19 | -6.31 | 0.12 | 1.53 | `35:-1.30` | -0.01 | -6.33 | -6.47 |
| `6:25` | -2.32 | -3.41 | 1.09 | 1.93 | `7:-12.55` | 0.00 | -2.97 | -4.10 |
| `6:125` | -26.03 | -24.79 | -1.24 | 1.83 | `7:-66.96` | 0.00 | -29.49 | -28.26 |
| `6:35` | -4.20 | -4.85 | 0.66 | 1.81 | `35:-0.91` | 0.00 | -4.18 | -4.77 |
| `6:135` | -84.81 | -83.36 | -1.46 | 2.01 | `7:-102.64` | 0.00 | -84.93 | -83.40 |
| `6:235` | -3.19 | -3.86 | 0.67 | 0.98 | `14:4.58` | 0.00 | -3.40 | -4.07 |
| `6:1235` | -0.67 | -0.84 | 0.18 | 1.12 | `35:-1.39` | 0.01 | -0.81 | -1.03 |
| `6:45` | -39.56 | -39.89 | 0.34 | 1.28 | `7:-25.17` | -0.00 | -39.81 | -40.17 |
| `6:145` | -14.13 | -14.10 | -0.03 | 1.06 | `7:-6.74` | -0.00 | -14.33 | -14.25 |
| `6:245` | 0.27 | 1.34 | -1.07 | 1.42 | `10:0.77` | -0.00 | 0.16 | 1.20 |
| `6:1245` | -23.97 | -22.87 | -1.11 | 1.11 | `7:-19.19` | 0.01 | -24.74 | -23.66 |
| `6:345` | -84.07 | -84.69 | 0.62 | 1.09 | `7:-88.05` | -0.00 | -85.86 | -86.42 |
| `6:1345` | -7.15 | -5.09 | -2.06 | 2.06 | `7:-1.50` | 0.01 | -7.23 | -5.21 |
| `6:2345` | -9.53 | -8.94 | -0.59 | 1.30 | `7:-25.18` | -0.00 | -10.31 | -9.73 |
| `6:12345` | 7.48 | 5.24 | 2.24 | 2.24 | `1:3.23` | 4.85 | 9.99 | 7.86 |
| Rudin-Shapiro | -0.05 | -0.51 | 0.46 | 1.30 | `3:0.27` | 0.00 | -0.14 | -0.56 |

**At each of the 49 quarter-octave checkpoints from `2^20` to `2^32`, all 55 distinct designs read within `2.35 sigma` of their major-arc sum, the largest at `4:23`; 48 stay within `2 sigma`, and the plain sums all stay within `2.5 sqrt(pi)`. The pooled rms of the residual is `0.715` against the pooled floor `0.710`, a ratio of `1.008`; at `2^32` alone it is `0.936` against `0.781`, `chi^2/55 = 1.44`. Verified.** The 55 hold 8 sign-twin pairs along the odd primes, `3:2`/`3:12`, `5:1`/`5:3`, `5:2`/`5:123`, `5:12`/`5:23`, `5:4`/`5:134`, `5:14`/`5:34`, `5:24`/`5:1234` and `5:124`/`5:234`, so they are 47 independent streams; over the 47 the pooled rms is again `0.715`. The readings run from noise, `0.23` at `3:12`, to `163.83` at `6:2` and `-14240.03` at the folded class `3:1`. The folded class is `n mod 2`, `-1` at every odd prime, and its model is `-x` up to `O(Q)`, from `q = 2`; its residual `0.04 sigma` is the error in the prime number theorem. Rudin-Shapiro, whose spectrum is flat, has a model of noise size whose residual has rms `0.517` and largest value `1.30`. The coin, which has no model, reads `z` with rms `0.912` over the checkpoints, where random signs give 1.

The truncation converges. `ladder` sums the expansion to `Q = 4096` at three points, 3.6 million fractions in 5 s, and reports the rms and the largest residual over the 54 distinct non-folded designs and Rudin-Shapiro. Beside them it prints the floor for each `Q`, which falls as `Q` grows because `H` does, the ratio of the rms to it, and the excess variance `rms^2 - floor^2` in units of `sigma^2`.

| `x` | | `Q = 6` | 30 | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096 |
|---|---|---|---|---|---|---|---|---|---|---|
| `2^28` | rms | 8.668 | 2.380 | 1.073 | 0.986 | 0.950 | 0.869 | 0.825 | 0.764 | 0.816 |
| `2^28` | max | 32.281 | 8.054 | 2.560 | 2.155 | 2.020 | 2.079 | 2.055 | 1.687 | 1.769 |
| `2^28` | floor | 0.907 | 0.862 | 0.839 | 0.815 | 0.791 | 0.767 | 0.742 | 0.717 | 0.690 |
| `2^28` | ratio | 9.552 | 2.761 | 1.279 | 1.210 | 1.201 | 1.132 | 1.111 | 1.065 | 1.183 |
| `2^28` | excess | 74.312 | 4.922 | 0.447 | 0.308 | 0.277 | 0.166 | 0.129 | 0.069 | 0.190 |
| `2^30` | rms | 16.832 | 2.818 | 1.284 | 1.399 | 0.984 | 0.810 | 0.796 | 0.777 | 0.762 |
| `2^30` | max | 75.696 | 9.091 | 3.625 | 4.590 | 3.684 | 1.989 | 2.071 | 2.309 | 1.850 |
| `2^30` | floor | 0.914 | 0.872 | 0.851 | 0.830 | 0.808 | 0.786 | 0.763 | 0.740 | 0.716 |
| `2^30` | ratio | 18.411 | 3.231 | 1.508 | 1.686 | 1.219 | 1.031 | 1.043 | 1.051 | 1.065 |
| `2^30` | excess | 282.475 | 7.182 | 0.924 | 1.268 | 0.317 | 0.039 | 0.052 | 0.057 | 0.069 |
| `2^32` | rms | 21.399 | 3.764 | 1.793 | 1.444 | 1.190 | 0.942 | 0.938 | 0.909 | 0.902 |
| `2^32` | max | 82.532 | 11.990 | 5.977 | 4.943 | 3.971 | 2.395 | 2.243 | 2.236 | 2.456 |
| `2^32` | floor | 0.920 | 0.881 | 0.862 | 0.842 | 0.822 | 0.801 | 0.781 | 0.760 | 0.738 |
| `2^32` | ratio | 23.258 | 4.271 | 2.081 | 1.715 | 1.449 | 1.175 | 1.201 | 1.196 | 1.223 |
| `2^32` | excess | 457.079 | 13.388 | 2.474 | 1.376 | 0.742 | 0.245 | 0.269 | 0.249 | 0.269 |

**At `2^32` the rms falls from `21.399 sigma` at `Q = 6` to `0.942` at `Q = 512` and `0.902` at `Q = 4096`, while the floor falls from `0.920` to `0.801` and `0.738`; for `Q = 512..4096` the rms stays `1.175` to `1.223` times the floor, an excess variance of `0.245` to `0.269 sigma^2` that does not shrink with `Q`. At `2^30` the ratio is `1.031` to `1.065` over the same `Q`, and at `2^28` `1.065` to `1.183`. Verified.** The excess at `2^32` is a finding, not noise. Denominators beyond 4096 are candidates: Thue-Morse's best orbit on `q = 255` grows at `g = 0.659` and on `q = 65535` at `g = 0.703`, and `T_65535` is already `0.20 sigma`, or `0.93 sqrt(x)`, at `2^32`. The rms first falls below `1.1` at `Q = 64` at `2^28` but only at `Q = 512` at `2^32`, as orbits with larger `q` and larger `g` rise above the noise. The heaviest terms name the orbits. Base 4 is led by the fifths, a 2-cycle under `x 4`, with `T_5 = -107.43` for `4:13` and `4:23`, then `T_17 = -5.18` and `T_257 = -2.11`. Base 6 is led by the fifths, fixed by `x 6`, for `6:14`, `6:23`, `6:123`, `6:124`, `6:134` and `6:234`. The sevenths, a 2-cycle under `x 6`, lead `6:125`, `6:135` and `6:345`, where `T_7 = -102.64` for `6:135`. The fixed point 0 leads the singletons. The preimages follow each orbit: `T_30 = -16.18` beside `T_5 = -97.77` for `6:14`, and `T_42 = 12.84` beside `T_7` for `6:135`.

**Conjecture: for every design of base at most 6, `sum_(n<x) x(n) Lambda(n) = M_(Q(x))(x) + O(x^(1/2 + epsilon))` with `Q(x) = x^(1/2)`, so the bias along the primes has the order `x^g` of the best orbit through a squarefree denominator, `g > 1/2`, or is noise.** A slowly growing `Q(x)` would make it false. An orbit that follows the cycle of `1/3` for most of its period has `g` close to `log_4 3`, and its denominator `2^m - 1` can be squarefree: `T_65535`, with `g = 0.703`, is `0.93 sqrt(x)` at `2^32` and grows past `x^(1/2 + epsilon)`, while it lies beyond any `Q(x)` of order `log x`. The census supports the conjecture with `Q = 1024` at every checkpoint above, up to the excess at `2^32`; a proof is out of reach.

## What the weights see

Along the primes the weight of `a/q` is `mu(q)/phi(q)`; along the squares it is the quadratic Gauss sum `G(a, q)/q` ([the squares](squares.md)). The two readings differ exactly where the two weights do.

**The unit digits.** The `q | rad(b)` part of the expansion is `(rad(b)/phi(rad(b))) sum_(n<x, (n,b)=1) x(n)`, the integers whose last digit is a unit, reweighted. **At `x = 6^k` the `q | 6` part of a base-6 design is `3 (v(1) + v(5)) P(1)^(k-1)`, `P(1) = sum_(d<6) v(d)`. Proved.** The integers below `6^k` prime to 6 are `6m + r` with `r` in `{1, 5}` and `m < 6^(k-1)`, and `x(6m + r) = x(m) v(r)`. So the fixed point 0, the spectral peak of the five singletons `6:d` with `abs(P(1)) = 4` and `g = log_6 4 = 0.774`, lights along the primes only if `v(1) + v(5) != 0`. **At `2^32` the `q | 6` part reads `162.03 sigma` for `6:2` and `6:3` and `162.07 sigma` for `6:4`, and `0.00` for `6:1` and `6:5`, which give the units opposite signs; the full models are `163.59`, `159.96` and `161.87` against the census's `163.83`, `159.98` and `161.41`, and `-1.68`, `-0.84` against `-0.16`, `-1.56`. Verified.** This is the squares chapter's rule of the units, now with its coefficient. It reverses a mirror pair. `6:5` and its mirror `6:12345` share one spectral measure, and along the squares they read `X2 = 13152` against `714`. Along the primes `6:5` is dark, `-1.56 sigma`, its units cancelling, while `6:12345`, with `P(1) = -4` and `v(1) + v(5) = -2`, has a `q | 6` part `-6 (-4)^(k-1)` at `6^k` that alternates in sign; it reads `4.85 sigma` at `2^32` and the census `7.48`.

**The fifths of base 6.** Since `6 = 1 (mod 5)`, every `a/5` is fixed by `x 6`, and the fifths enter with their preimages `a/10`, `a/15`, `a/30`. **At `x = 6^k` the `q in {5, 10, 15, 30}` part of a base-6 design is `-(3/4) sum_(a=1..4) (v(1) zeta^a + v(5)) P(zeta^a)^(k-1)`, `zeta = e(1/5)`. For `6:23` it is `-(3/4) L_(3k-1)` and for `6:14` it is `-(15/4) F_(3k-2)`, Lucas and Fibonacci numbers. Proved.** Since `c_(5q) = c_5 c_q` for `q` prime to 5, the four terms sum to `-(3/4) sum_(n<x, (n,6)=1) x(n) c_5(n)`. Write `n = 6m + r`; then `c_5(n) = sum_a zeta^(a m) zeta^(a r)` and `sum_(m<6^(k-1)) x(m) zeta^(a m) = P(zeta^a)^(k-1)`. For `6:23`, `P(zeta^(+-1)) = 2 + sqrt 5 = phi^3` and `P(zeta^(+-2)) = 2 - sqrt 5 = psi^3`, with `phi` and `psi` the roots of `z^2 = z + 1`. Then `2(1 + cos 72 deg) = phi^2` and `2(1 + cos 144 deg) = psi^2` give `L_(3k-1)`. For `6:14` the two values trade places and `v(1) = -1`, which gives `5 F_(3k-2)`. The generator checks the general form for all 31 base-6 designs at `k <= 12` to `6.0e-14`, and the two closed forms to `1.5e-15`: `-2.25`, `-8.25`, `-35.25`, `-149.25` for `6:23` and `-3.75`, `-11.25`, `-48.75`, `-206.25` for `6:14` at `k = 1..4`. The growth exponent is `log_6 phi^3 = 0.8057 = 1 - 0.1943`, the squares chapter's `eta` for both designs.

So `6:23` is slow along the primes for a reason that holds at every scale: its fifths carry a Lucas number with coefficient `-3/4`, and nothing cancels it. Along the squares the same orbit is weighted by Gauss sums, and its coefficient there cancels, as the squares chapter found. The census reads `6:23` at `-31.23 sigma` against a model of `-31.32` at `2^32`. Along the first `2^26` primes the plain sum reads `-46.7` against `-45.9`, and `g - 1/2 = 0.306` is the slope the squares chapter measured, `0.307`.

**The quarter.** `3:2` has its spectral peak on the orbit `1/4`, with `eta = 0.2675`, so `g = 0.7325`. That would make it slow along the primes, but `mu(4) = 0`: the primes are equidistributed mod 4 and carry no mass at `1/4`. **`3:2` reads `-0.23 sigma` along the primes at `2^32`, against a model of `0.10`, while its best orbit has `g = 0.7325 > 1/2`. Verified.** Along the squares the same orbit is visible, since `G(1, 4) != 0`. The two weights fail in opposite places: the Gauss sum vanishes at `q = 2 (mod 4)`, where `mu` does not, and `mu` vanishes at `4 | q`, where the Gauss sum does not. The twin `5:34` shows the other side: its peak `1/6` is invisible to the bias along the squares, `G(1, 6) = 0`, and visible along the primes, `mu(6) = 1`. It reads `136.69 sigma`, the exact negative of `5:14` at every odd prime.

**The floor.** Along the squares the twisted sums run to `N^2` and a term is of size `N^(1 - 2 eta)` against the noise `sqrt N`, so the floor sits at `eta = 1/4`. Along the primes they run to `x` against the noise `sqrt(x log x)`, so the floor sits at `eta = 1/2`, or `g = 1/2`. Every design of base at most 6 has `eta < 0.342` ([the squares](squares.md)), so whether it lights along the primes depends on which denominators carry its good orbits and on the weights there.

## Two truncations

The expansion can be cut in two ways. Cutting at `q <= Q` keeps the term `q = 1` exact. Cutting at the divisors of `P(z) = prod_(p<=z) p` gives the sieve model, since `(P/phi(P)) 1[(n, P) = 1] = sum_(q | P) (mu(q)/phi(q)) c_q(n)` for squarefree `P`: the primes as the integers free of small factors, reweighted. The census counts the sieve model directly, binning every `n` prime to 210 by its least prime factor up to 1021. **The sieve model at `z = 7` equals the `q | 210` part of the expansion, computed by the digit walk, to `1.6e-12 sigma` at each of the 98 quarter-octave checkpoints from `2^8`. Verified.** The two codes share nothing but the designs.

**At `z = 1021` the sieve model misses the folded class `3:1` by `45.32 sigma` at `2^32`, while the cut at `q <= 1024` meets it within `0.04 sigma`. Verified.** The sieve keeps every `q | P(1021)`, most of them far above `x`, and its density `(P/phi(P)) #{n < x: (n, P) = 1}/x` is still `1.00318` at `2^32`, where `u = log x/log z = 3.2`: the fundamental lemma of the sieve converges only as `u` grows. That is invisible for small biases and decisive for the folded class, whose bias is `-x`. On the designs the sieve is otherwise close: Thue-Morse reads `0.93 sigma` from it and the coin `-0.05`, with sd `0.644` for random signs. The cut at `q <= Q` is the cleaner model.

## What a proof needs

The obstruction is the minor arcs. The main term of Thue-Morse alone needs `sum_(n<x) t(n) Lambda(n) + G_6 x^(log_4 3) = o(x^(log_4 3))`, a saving of `x^(-0.2075)` over the trivial bound together with an asymptotic. Mauduit and Rivat 2010, Theorem 1, give `O(x^(1 - sigma))` with `sigma` only asserted to exist, and their explicit constants place it below a 28th of a type II quantity ([refs](../refs.md)). The law at the noise floor asks for more: every major arc to square-root accuracy, as from the generalised Riemann hypothesis, and minor arcs bounded by `x^(1/2 + epsilon)`. Every term of the main sum is exact; the difficulty is to show that nothing else is there.

What is known. The integers prime to 6 as a model for the primes, with the exponent `log 3/log 4`, are Shevelev's: his 2007 Conjecture 2, restated in 2008 as Conjecture 5.2, puts `log(pi^o(n) - pi^e(n))/log n -> log 3/log 4` for the odious and evil primes, and his 2007 Theorem 2(1) is the closed form `-3^floor(k/2)` above. In 2008, Theorem 5.4, he proves the analogue for the integers whose least prime factor is a fixed `p >= 5`. Mauduit and Moreira 2015 prove the analogue for the squarefree integers, `sum_(n<N) mu^2(n) t(n) = -(2/pi^2)(1 + o(1)) sum_(n<N, 3 | n) t(n)`, and record Shevelev's conjecture that the sum over the first `N` primes is negative from `N = 11`. OEIS A130911 counts the race between odious and evil primes: the odious lead holds to `10^13` and, in a later check, to `10^19`, where the race has a fractal profile like the race among the integers prime to 3; A127977 gives its minimum lead on dyadic intervals to `2^63`. Not found in the sources read: the terms past `q = 6`, the major-arc sum for any other design, a test of any of it at the noise floor, and any Omega-result for the prime sum.

## Doors

- **An Omega-result.** Prove that `sum_(p<x) t(p)` is not `o(x^(log_4 3 - epsilon))`, or that it is negative for `x >= 32`, under the generalised Riemann hypothesis and a minor-arc hypothesis the census can test.
- **The second exponent.** Past `q = 6` the Thue-Morse bias is a sum of orbits; among its largest terms at `2^32` the growth reaches `0.702`, at `q = 257`. Find the exponent of the remainder, the supremum of `g` over the orbits weighted by `1/phi(q)`, and whether it is a power at all.
- **The excess at `2^32`.** For `Q = 512..4096` the residual variance sits `0.245` to `0.269 sigma^2` above the floor at `2^32` and not at `2^30`. Sum the expansion to `Q = 2^16` at `2^32` and find whether the large orbits, `T_65535` among them, carry it.
- **The coefficients of the sevenths.** `6:125`, `6:135` and `6:345` carry their bias on the 2-cycle of sevenths under `x 6`, where `P(zeta_7^a) P(zeta_7^(-a))` is an algebraic integer of the cubic field. Their main terms at `6^k` should be linear recurrences like the Lucas numbers of `6:23`.
- **Other readings.** Along `p + 2`, along the primes in a progression, or along `floor(p^c)`, the weights change from `mu(q)/phi(q)` to the singular series of the reading. The census and the digit walk carry over unchanged.
- **The lags.** The squares chapter's gap model gives the lags along the consecutive primes and the major-arc sum gives the bias; together they could predict the 8-block chi-square of a design along the primes.
