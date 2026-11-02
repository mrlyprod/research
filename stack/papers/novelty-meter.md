---
title: The Stack Hears the Zeros
lead: The novelty of the stack, read through a smooth window, has an error that is a sum of waves at the frequencies of the zeros of zeta; its decay exponent is a Riemann hypothesis meter, and a sharp window hears the primes instead.
date: 2026-09-21
figure: paper-novelty-meter
---

Lay a grid of `n` equal cells on the unit interval for every `n` up to `N` and stack the layers. A reduced fraction `a/b` is drawn first by layer `b`, so layer `n` lights `phi(n)` nodes for the first time: the novelty of scale `n` is Euler's totient, and its Dirichlet series is `zeta(s-1)/zeta(s)`. Read the novelty of all scales near `1/y` through a window `f` supported on `[1, 2]`, as `S_f(y) = sum_n phi(n) f(ny)`; the main term is `(6/pi^2) F(2) y^-2` with `F` the Mellin transform of `f`, and the scaled error `E_f(y) = y^2 S_f(y) - (6/pi^2) F(2)` is the meter. A theorem of Verjovsky, read here at source, says that the meter decays like `y^(3/2 - eps)` for every `eps` and every smooth `f` if and only if the Riemann hypothesis holds, while for the indicator of `[1, 2]` the exponent is `1` and no better. This paper assembles that equivalence from the Mellin identity with both inputs named, proves the sharp exponent from the jump of the totient sum at a prime, writes the smoothed error under the hypothesis as a sum over the zeros, and reads the meter: on `y = 2^-j` for `j` from `8` to `23.5` with the totients sieved to `3 * 10^7`, the fitted exponents in `q = y^2` are `0.5023` for the indicator and `0.7498` and `0.7471` for a `C^2` and a `C^infinity` bump, against `1/2` and `3/4`, and the residue sum over the first 138 zeros reproduces the smooth bump's error over the whole grid to a relative `2.0e-6`. The stack hears the zeros. It cannot tell where they are: on any finite range the smoothed error is a rendering of the first zeros, a zero off the line at a height already checked would be invisible by the same decay that makes the sum converge, and the meter is the hypothesis restated, not a way in.

## Introduction

![The smoothed novelty error, scaled by its own `y^(3/2)`, as beads on the wave the first 29 zeros of zeta predict from their heights, two values of `zeta` and the window's transform at each; below it the sharp window's error, scaled by `y`, which no zero sum follows.](paper-novelty-meter)

Draw the lines `x = a/n` for `a = 0, ..., n` on the unit interval, once for every `n` from `1` to `N`, with the ink thin enough that lines add. That is the [stack](../notes/stack.md), and the [Farey page](../../farey/notes/farey.md) says what it draws: the [Farey fractions](https://en.wikipedia.org/wiki/Farey_sequence) of denominator at most `N`. Layer `n` adds the fractions of denominator exactly `n`, `phi(n)` of them in `(0, 1]`, so a prime scale adds `n - 1` new nodes and a composite one fewer. That count is the novelty of the scale, and it is the only arithmetic this paper needs.

The question is how regularly the novelty grows. Summed sharply, `sum_(n <= x) phi(n)` is `3x^2/pi^2` up to an error that jumps by about `0.39 p` at every prime `p`, so no bound smaller than `x` is possible. Summed through a smooth window, the jumps average out and the error becomes a sum of waves, one per zero of the [Riemann zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function), the wave of `rho = 1/2 + i gamma` oscillating like `cos(gamma log y)` with an amplitude falling like `y^(3/2)`. The figure above is that sum: the dots are the measured error of one smooth window scaled by `y^(3/2)` at 201 values of `y` between `2^-8` and `2^-20.5`, the line through them is the prediction of the first 29 zeros, computed from their heights, two values of `zeta` at each and the window's Mellin transform there, with nothing fitted, and below it the sharp window's error at the same scales, scaled by `y`, is a cloud.

The decay exponent of the smoothed error is a meter for the Riemann hypothesis. If the hypothesis holds, every wave has amplitude `y^(3/2)` and the error is `O(y^(3/2 - eps))`; if it fails, a zero off the line contributes a wave that decays more slowly, and the rate of decay over all smooth windows locates the rightmost zero. That is a theorem of Verjovsky (1994), stated in Section 3 as at source and translated into this paper's variable. Section 4 assembles the equivalence from the Mellin transform and derives the wave sum; from the rate to the zeros needs nothing beyond absolute convergence and Weierstrass, the other way needs one input from the literature, the bound on `1/zeta` right of the critical line under the hypothesis with the horizontal lines the contour closes on, cited at the page that proves it. Section 5 proves the sharp exponent from the prime jump, Section 6 reads the meter, every number a Fact with its finite domain and the study that prints it, and Section 7 says what the reading is worth.

## Definitions

**Definition 2.1 (novelty).** The novelty of scale `n >= 1` is `phi(n)`, the number of reduced fractions in `(0, 1]` with denominator `n`. The novelty series is `sum_n phi(n) n^-s`.

**Lemma 2.2.** For `Re s > 2`, `sum_n phi(n) n^-s = zeta(s-1)/zeta(s)`.

**Proof.** `sum_(d | n) phi(d) = n`, since every `a/n` with `1 <= a <= n` reduces to exactly one fraction of denominator `d | n`. Multiplying the two absolutely convergent series, `zeta(s) sum_n phi(n) n^-s = sum_n (sum_(d | n) phi(d)) n^-s = sum_n n^(1-s) = zeta(s-1)`. □

**Definition 2.3 (window and meter).** A window is a bounded measurable `f` on `(0, infinity)` supported in `[1, 2]`. Its Mellin transform `F(s) = int_0^infinity f(u) u^(s-1) du` is entire. The smoothed novelty is `S_f(y) = sum_n phi(n) f(ny)`, a finite sum over the scales `n` between `1/y` and `2/y`, and the meter is

```
E_f(y) = y^2 S_f(y) - (6/pi^2) F(2)
```

read as a power of `y` as `y -> 0`, or of `q = y^2`, which is the variable of the source. An exponent `a` in `q` is `2a` in `y`.

**Definition 2.4 (the three windows).** The sharp window is the indicator of `[1, 2]`, `F(s) = (2^s - 1)/s`, `F(2) = 3/2`. The `C^2` window is `64 (u-1)^3 (2-u)^3` on `(1, 2)`, `F(2) = 24/35`. The `C^infinity` window is `exp(4 - 1/((u-1)(2-u)))` on `(1, 2)`, `F(2) = 0.575725895994` by a 2000-node Gauss-Legendre rule (`lab/py/smoothed-novelty`, `main`). The constants `64` and `4` make both bumps peak at `1` at `u = 3/2`.

## The theorem at source

Verjovsky, *Discrete measures and the Riemann hypothesis*, Kodai Math. J. 17 (1994), defines on page 596, for `y > 0`, the measure `m_y(f) = sum_(n in N) y phi(n) f(y^(1/2) n)` on functions of compact support in the positive reals, and on page 597 `m_0(f) = int_0^infinity (6/pi^2) u f(u) du`. Its `y` is this paper's `q`, so `m_q(f) = y^2 S_f(y)` at `q = y^2`, `m_0(f) = (6/pi^2) F(2)`, and `E_f(y) = m_q(f) - m_0(f)`. Page 597 states, and the paper proves:

**Theorem A.** For every `f` in `C_c^0(R*)`, `m_y(f) = (6/pi^2) int_0^infinity u f(u) du + O(y^(1/2) log y)` as `y -> 0`.

**Theorem B.** (1) The Riemann hypothesis is true if and only if for every `f` in `C_c^r(R*)`, `2 <= r <= infinity`, `m_y(f) = m_0(f) + o(y^(3/4 - eps))` as `y -> 0`, for all `eps > 0`. (2) For `alpha` in `(1/2, 3/4)`, `m_y(f) = m_0(f) + o(y^(alpha - eps))` for all `f` in `C_c^2(R*)` and all `eps > 0` holds if and only if `zeta` has no zeroes in the half-plane `Re s > 2(1 - alpha)`. (3) If `f` is the characteristic function of an interval then `limsup_(y -> 0) y^-alpha |m_y(f) - m_0(f)| = infinity` if `alpha > 1/2`; hence `1/2` is the best possible exponent of the error for some nonsmooth functions. Part (4), on one specific window and zeros near `Re s = 1`, is not used here.

The same page defines `m_y` of the characteristic function of `[a, b]` as the sum over `a y^(-1/2) <= n <= b y^(-1/2)`, both ends included, the convention of Definition 2.4, and page 598 derives Theorem A for it from `sum_(n <= x) phi(n) = 3x^2/pi^2 + O(x log x)`, attributed there to Mertens (1874). The 2017 text of the author's 1993 Pitman notes, arXiv:1711.03593, restates the smooth case as its Theorem 5.1 in the variable `y = q^(1/2)`: for all `g` in `C_0^infinity(R*)` the error is `o(y)`, and `o(y^(3/2 - eps))` for all `eps > 0` if and only if the Riemann hypothesis holds, the proof referred to the Kodai article. Both texts were read at source. In this paper's variable the smooth meter sits at `3/2 - eps` in `y` exactly when the hypothesis holds, and the sharp meter at `1` in `y` unconditionally and sharply.

## The meter is a Riemann hypothesis meter

**Theorem 4.1 (the Mellin identity).** Let `f` be a `C^2` window and `c > 2`. Then

```
S_f(y) = (1/(2 pi i)) int_(Re s = c) F(s) (zeta(s-1)/zeta(s)) y^-s ds
```

and the residue of the integrand at `s = 2` is `(6/pi^2) F(2) y^-2`.

**Proof.** Two integrations by parts, the boundary terms vanishing because `f` and `f'` vanish at `1` and `2`, give `F(s) = (1/(s(s+1))) int f''(u) u^(s+1) du`, so `|F(s)| <= 2^(Re s + 1) ||f''||_1 / (|s| |s+1|)` uniformly on vertical strips. Mellin inversion therefore holds pointwise with an absolutely convergent integral, `f(u) = (1/(2 pi i)) int_(Re s = c) F(s) u^-s ds`. Put `u = ny`, multiply by `phi(n)` and sum; `sum_n phi(n) n^-c` converges for `c > 2` and `int |F(s)| |ds|` is finite, so the sum and the integral exchange, and Lemma 2.2 gives the display. At `s = 2` the factor `zeta(s-1)` has a simple pole of residue `1` and `zeta(2) = pi^2/6`. □

**Theorem 4.2 (the hypothesis gives the rate).** Assume the Riemann hypothesis. Then for every `C^2` window `f` and every `eps > 0`, `E_f(y) = O(y^(3/2 - eps))`.

**Proof.** Move the line of Theorem 4.1 from `Re s = c` to `Re s = 1/2 + eps`. Two facts about `zeta` under the hypothesis are needed, and both are taken from Hu, Kaneko, Martin and Schildkraut (2023), where they are proved for any Dedekind zeta and, at `K = Q`, are Corollary 13.16 and Theorem 13.22 of Montgomery and Vaughan (2007) by that paper's own attribution. Its Lemma 5.4, with `tau = |t| + 4` and `|t| >= 1`, bounds `|log zeta(sigma + it)|` by `log(1/(sigma - 1)) + O(sigma - 1)` on `1 + 1/log log tau <= sigma <= 3/2`, by `log log log tau + O(1)` on `1 - 1/log log tau <= sigma <= 1 + 1/log log tau`, and by `log(1/(1 - sigma)) + O((log tau)^(2 - 2 sigma) / ((1 - sigma) log log tau))` on `1/2 + 1/log log tau <= sigma <= 1 - 1/log log tau`. The third range covers `1/2 + eps <= sigma <= 1 - 1/log log tau` once `tau` is large, the second and first cover `sigma` up to `3/2`, and beyond `3/2` the Euler product gives `|1/zeta(s)| <= zeta(3/2)`; together they exponentiate to `1/zeta(s) << tau^delta` for every `delta > 0`, uniformly on `Re s >= 1/2 + eps`. Its Lemma 2.4 gives, for each integer `n >= 4`, a height `T_n` in `[n, n + 1)` with `|zeta(sigma + i T_n)| >= exp(-C log n / log log n)` for `-1 <= sigma <= 2`. On the new line `zeta(s - 1)` is `O(|t|^(1 - eps))` by the functional equation, `F(s) = O(|s|^-2)` by the bound in the proof of Theorem 4.1, and `1/zeta(s) = O(|t|^delta)`, so the integrand is `O(|t|^(-1 - eps + delta))` and the integral converges absolutely for `delta < eps`, giving `O(y^(-1/2 - eps))`. On the horizontal segments at height `T_n` between the two lines the integrand is `O(T_n^(-1 + o(1)) y^-c)` and vanishes as `n -> infinity`. Under the hypothesis no zero of `zeta` lies in `Re s > 1/2 + eps`, so the only residue crossed is the one at `s = 2`, and `S_f(y) = (6/pi^2) F(2) y^-2 + O(y^(-1/2 - eps))`. Multiply by `y^2`. □

**Theorem 4.3 (the rate gives the hypothesis).** Fix `eps > 0`. If `E_f(y) = O(y^(3/2 - eps))` for every `C^2` window `f`, then `zeta` has no zero in `Re s > 1/2 + eps`. If it holds for every `eps > 0`, the Riemann hypothesis holds.

**Proof.** `S_f(y) = 0` for `y > 2`, since `f(ny) = 0` once `ny > 2`. For `Re s > 2`, by absolute convergence, `int_0^2 S_f(y) y^(s-1) dy = sum_n phi(n) int_0^infinity f(ny) y^(s-1) dy = F(s) zeta(s-1)/zeta(s)`. Writing `S_f(y) = ((6/pi^2) F(2) + E_f(y)) y^-2`, the left side is

```
(6/pi^2) F(2) 2^(s-2)/(s-2) + int_0^2 E_f(y) y^(s-3) dy .
```

Under the rate the integrand is `O(y^(Re s - 3/2 - eps))` near `0` and bounded near `2`, so the integral converges absolutely and locally uniformly on `Re s > 1/2 + eps` and the display is holomorphic there except for the simple pole at `s = 2`. So `F(s) zeta(s-1)/zeta(s)` continues holomorphically to `Re s > 1/2 + eps` away from `s = 2`. Let `rho` be a zero of `zeta` there. Since `zeta` has no zeros on `Re s >= 1`, `Re rho < 1`, and `Re(rho - 1)` lies in `(-1/2, 0)`, where `zeta` has no zeros; so `zeta(rho - 1) != 0`, and holomorphy at `rho` forces `F(rho) = 0`. Now `u^k f(u)` is again a `C^2` window for every integer `k >= 0`, with Mellin transform `F(s + k)`, so the same argument gives `F(rho + k) = 0` for every `k`: the continuous function `g(u) = f(u) u^(rho - 1)` on `[1, 2]` is orthogonal to every polynomial. By Weierstrass, polynomials are dense in `C[1, 2]`, so `g` is orthogonal to its own conjugate and `g = 0`, whence `f = 0`. Taking any nonzero window gives the contradiction. The second sentence follows by letting `eps -> 0`. □

Theorems 4.2 and 4.3 together are Theorem B part (1) at `r = 2` for windows on `[1, 2]`, and the two halves do not cost the same: the backward half is elementary, the forward half rests on a bound on `1/zeta` that is itself a consequence of the hypothesis.

**Theorem 4.4 (the wave sum).** Assume the Riemann hypothesis and that the nontrivial zeros of `zeta` are simple, and let `f` be a `C^infinity` window, or a `C^2` window whose transform is `O(|s|^-4)`, as the `C^2` bump of Definition 2.4 is. Then for every `0 < delta < 1/2`,

```
E_f(y) = sum_rho F(rho) (zeta(rho - 1)/zeta'(rho)) y^(2 - rho) + O(y^(3/2 + delta)) ,
```

the sum over the nontrivial zeros taken as residues in the order of the heights `T_n` of Theorem 4.2. Each term has modulus `|c_rho| y^(2 - Re rho)` with `c_rho = F(rho) zeta(rho - 1)/zeta'(rho)`, which is `|c_rho| y^(3/2)` on the line; a zero of real part `beta` would contribute a term of exponent `1 - beta/2` in `q`, which is `3/4` on the line and smaller at the member of a pair `beta, 1 - beta` that sits right of it.

**Proof.** Integrating by parts `k` times, `F(s) = O(|s|^-k)` for every `k` when `f` is `C^infinity`, since every derivative vanishes at the endpoints; for the `C^2` bump three integrations by parts and one more boundary term give `O(|s|^-4)`, `f'''` being its first derivative not vanishing at the endpoints. Shift the line of Theorem 4.1 to `Re s = 1/2 - delta`, crossing the pole at `s = 2` and the zeros. On the new line `zeta(s - 1) = O(|t|^(1 + delta))` and `1/zeta(s) = O(|t|^(delta'))` for every `delta' > 0`, the latter by the functional equation and the bound of Theorem 4.2 at `Re s = 1/2 + delta`, so against `F = O(|s|^-4)` the integrand is `O(|t|^(-3 + delta + delta'))`, the integral converges absolutely and is `O(y^(-1/2 + delta))`. The horizontal segments at the heights `T_n` vanish as before, and the residue at a simple zero `rho` is `F(rho) zeta(rho - 1) y^-rho / zeta'(rho)`. Multiply by `y^2`. The exponent count is `|y^(2 - rho)| = y^(2 - Re rho)` and `q = y^2`. □

The wave of `rho = 1/2 + i gamma` is `2 Re(c_rho y^(3/2 - i gamma)) = 2 |c_rho| y^(3/2) cos(gamma log y - arg c_rho)`, one cosine in `log y` per zero at the zero's height as frequency: the top panel of the figure. What smoothness buys is not the shift, which the `C^2` bump already admits, but the speed at which the coefficients die: for the `C^infinity` window `F(rho)` falls faster than any power of `gamma`, so a few dozen zeros are the whole error, while for the `C^2` bump they fall like a power and the sum converges slowly (Fact 6.2). For the indicator `F(s) = (2^s - 1)/s` is only `O(1/|t|)`, the contour cannot cross the line at all, and the sharp exponent comes from somewhere else.

## The sharp cutoff

**Theorem 5.1.** Let `R(x) = sum_(n <= x) phi(n) - 3x^2/pi^2`. Then `R(x) = O(x log x)`, and `R(x) = Omega(x)`: across every prime `p` it jumps by `phi(p) - 3(2p - 1)/pi^2 = (1 - 6/pi^2) p - 1 + 3/pi^2 > 0.39 p - 1`, so one of `|R(p - 1)|` and `|R(p)|` exceeds `0.19 p - 1`.

**Proof.** The identity `sum_(d | n) phi(d) = n` of Lemma 2.2 inverts to `phi = mu * id`, so `sum_(n <= x) phi(n) = sum_(d <= x) mu(d) T(floor(x/d))` with `T(m) = m(m+1)/2`. Since `T(floor(x/d)) = x^2/(2d^2) + O(x/d)` and `sum_(d <= x) mu(d)/d^2 = 6/pi^2 + O(1/x)`, the sum is `3x^2/pi^2 + O(x) + O(x sum_(d <= x) 1/d) = 3x^2/pi^2 + O(x log x)`. The jump is the definition of `R` at `p - 1` and at `p`, with `phi(p) = p - 1`. □

The upper bound is Mertens' (1874), recalled on page 598 of the source; the lower bound is the only part the meter feels.

**Corollary 5.2 (the sharp meter).** For the indicator `f` of `[1, 2]`, `E_f(y) = O(y log(1/y))` and `E_f(y) = Omega(y)`: at `y = 2/p` for an odd prime `p` the meter jumps by `y^2 phi(p) = (1 - 1/p) 2y`, so the exponent of the sharp meter is `1` in `y` and `1/2` in `q`, unconditionally, and it is attained.

**Proof.** `S_f(y) = sum_(1/y <= n <= 2/y) phi(n)`, so `E_f(y) = y^2 (R(floor(2/y)) - R(ceil(1/y) - 1)) + y^2 (3/pi^2)(floor(2/y)^2 - (ceil(1/y) - 1)^2) - 9/pi^2`, and the last two terms cancel to `O(y)`, which with Theorem 5.1 gives the upper bound. As `y` decreases through `2/p` the upper end `floor(2/y)` steps from `p - 1` to `p` while `ceil(1/y) = (p+1)/2` does not move, so `S_f` jumps by `phi(p)` and `E_f` by `y^2 (p - 1)` with `y = 2/p`; one side of the jump has `|E_f| >= (1 - 1/p) y`. □

This is Theorem B part (3) at `[1, 2]`, its proof reduced to the prime jump. The sharp meter sits at `1` in `y` because the primes' jumps are louder than the zeros' waves, which live at `3/2`: its exponent says nothing about the zeros. Smoothing is the whole difference between the two readings, and the difference is the quarter.

## The meter read

All numbers in this section are printed by `uv run python research/stack/lab/py/smoothed-novelty/smoothed_novelty.py`. The grid is `y = 2^-j` for `j` from `8` to `23.5` in steps of `1/16`, 249 samples, the largest window summing `1.19 * 10^7` scales; the totients are sieved to `3 * 10^7` and agree with brute-force gcd counts to `2000` and with the Mobius route of Theorem 5.1 at `x = 10^6`, both exactly (`totients`, `totient_sum_mobius`). The summation noise at `j = 23.5`, pairwise against compensated summation, is `5.6e-17` against `|E_f| = 5.4e-12`.

**Fact 6.1 (the slopes).** An octave is the set of samples with `j` in `[k, k + 1)`, 16 of them for `k` from `8` to `22` and the 9 of `[23, 23.5]` for the last, each placed at `j = k + 1/2`. Fit `log` of the root mean square of `E_f` over each octave against `log q` over the 16 octaves, and over the lower and upper eight. The slopes in `q` read:

| window | slope | lower eight | upper eight | residual |
| --- | ---: | ---: | ---: | ---: |
| indicator of `[1, 2]` | `0.5023` | `0.5296` | `0.4883` | `0.179` |
| `C^2` bump `64 (u-1)^3 (2-u)^3` | `0.7498` | `0.7553` | `0.7487` | `0.092` |
| `C^infinity` bump `exp(4 - 1/((u-1)(2-u)))` | `0.7471` | `0.7474` | `0.7460` | `0.115` |

**Table 1.** Exponents of the meter in `q = y^2`, against `1/2` for the sharp window and `3/4` for the smooth ones. Residuals are root mean square in `log` units; the eight-octave windows carry residuals between `0.065` and `0.169`, and no two windows of one row differ by more than `0.05`. Domain: the 249-sample grid above. Script: `slopes`. No exponent is claimed beyond that window.

The smoothed slopes do not read `1/2`: the smoothing gains the quarter that Theorem B states, and the `C^2` bump gains all of it, as part (1) says it must at `r = 2`.

**Fact 6.2 (the wave sum against the first 138 zeros).** The 138 zeros of `zeta` to height `300`, with `zeta(rho - 1)` and `zeta'(rho)` at each, are read from PARI in one call (`zeros_from_pari`; `lfunzeros`, `zeta`, `lfun` at derivative order `1`, 30 digits). For the `C^infinity` bump the residue sum of Theorem 4.4, `2 Re sum_rho c_rho y^(2 - rho)`, reproduces the measured `E_f` over the whole grid to a relative `2.0e-6` in the maximum norm; one zero alone gives `0.50`, ten `2.7e-2`, thirty `2.0e-3`, and the coefficients `|c_rho|` fall from `1.879e-1` at the first zero to `1.475e-7` at the 138th. The scaled amplitude `|E_f|/y^(3/2)` lies in `[1.3e-4, 0.558]` over the grid against `2 sum |c_rho| = 0.755` over the same zeros. For the `C^2` bump the same sum stops at `2.3e-4`, its coefficients falling only like a power, the printed least-squares power against `gamma` being `-3.22` on the 138. Domain: the 249-sample grid, zeros to height `300`. Script: `main`, `mellin_cinf`, `mellin_c2`.

**Fact 6.3 (the prime jump).** `R(p - 1) = 108941.6` and `R(p) = 501014.9` at `p = 1000003`, a jump of `392073.4`; `R(p - 1) = 1011363.8` and `R(p) = 4932099.6` at `p = 10000019`, a jump of `3920735.7`; both equal `phi(p) - 3(2p - 1)/pi^2` to the printed digit. Script: `main`, SHARP CUTOFF.

The figure is the same computation at a smaller height, run by its own figure, `paper-novelty-meter`: `phi` sieved to `3 * 10^6` with `mrly_num::factor::totients`, `j` from `8` to `20.5` in steps of `1/16`, 201 samples, so that every window fits under the sieve, the 29 zeros below height `100` from `mrly_num::zeta::Line::zeros`, `zeta(rho - 1)` and `zeta'(rho)` by `mrly_num::zeta::Line::pair`, checked against the crate's value at `1/2 + 30i` and against `zeta(-1/2)` and `zeta'(1/2)`, and `F(rho)` by a 4096-node quadrature. It asserts `F(2)` against `0.575725895994`, `|c_rho|` at the first and tenth zeros against `1.879e-1` and `4.286e-3`, and the 29-zero sum against the drawn samples within one percent of their peak. The top panel is `E_f(y)/y^(3/2)` for the `C^infinity` bump, dots, with the 29-zero sum as the line; the bottom panel is `E_f(y)/y` for the indicator, dots alone; both run from `j = 8` on the left to `20.5` on the right.

## What the meter does not do

Three things, said plainly. First, the equivalence of Section 4 is exactly as hard as the Riemann hypothesis: it is the hypothesis rewritten in `y`, and the same wall is reached by the Farey discrepancy on [the Farey page](../../farey/notes/farey.md) and by the Gaussian Franel identity on [the stack page](../notes/franel-and-eisenstein.md). Second, a finite reading cannot certify or refute. Fact 6.2 says the smoothed error on the whole grid is the first 138 zeros to six digits, and the reason is Theorem 4.4: the Mellin transform of a smooth window kills the coefficients of high zeros faster than any power, so a zero off the line at a height where the hypothesis has already been checked would contribute a wave too small to see, by the same decay that makes the sum converge. The smoothed novelty meter is not a route to the hypothesis. Third, the slopes of Table 1 are fits on sixteen octaves with residuals of a tenth in `log`, and the residual is not noise but the waves themselves, whose octave averages beat against one another; no exponent is claimed beyond the window read.

## Open problems

Theorem 4.4 is written under the hypothesis and for simple zeros. A multiple zero contributes the residue of the same integrand and changes nothing else. Dropping the hypothesis is not so cheap: the horizontal lines of the contour come from Lemma 2.4 of Hu, Kaneko, Martin and Schildkraut, which assumes it, and the unconditional substitute, `1/zeta << exp(C log^2 T)` on chosen heights (Titchmarsh 1986, chapter 9), beats every polynomial decay of `F`, so a general `C^infinity` window does not obviously survive; a window whose transform decays like `exp(-c sqrt(|t|))` would, and whether the bump used here does is not checked. Neither version is written here. Two questions are left open beyond that: which window of a given support makes a given zero's coefficient `|F(rho) zeta(rho - 1)/zeta'(rho)|` largest relative to the rest, so that the meter hears that zero best, nothing here optimising `f`; and whether the octave-to-octave residual of the smooth fits in Table 1, `0.092` and `0.115`, is the beat of the first few waves alone, which it should be by Fact 6.2 and which has not been checked. Neither is claimed.

## Reproducibility

One study, `lab/py/smoothed-novelty`, prints every number of Section 6: `uv run python research/stack/lab/py/smoothed-novelty/smoothed_novelty.py` from the repository root, one core, about fifteen seconds, needing numpy and `gp` on the path, taking no arguments, reading and writing no file; its README names the witness of every printed line. The figure is `figures/paper-novelty-meter.ts`, drawn in the browser in under a second a theme, with its sieve and `F(2)` in the census writer `research/shared/census/src/bin/paper-novelty-meter.rs`, and every quantity it draws is asserted by one or the other as listed in Section 6, so a wrong sieve, zero, `zeta` or Mellin transform stops the writer or the press.

## References

- Verjovsky 1994, Discrete measures and the Riemann hypothesis, Kodai Math. J. 17, no. 3, 596-608. [doi.org/10.2996/kmj/1138040054](https://doi.org/10.2996/kmj/1138040054)
- Verjovsky 2017, Arithmetic, geometry and dynamics in the unit tangent bundle of the modular orbifold, the updated text of the 1993 Pitman Research Notes article, Theorem 5.1. [arxiv.org/abs/1711.03593](https://arxiv.org/abs/1711.03593)
- Hu, Kaneko, Martin and Schildkraut 2023, On a Mertens-type conjecture for number fields, Math. Proc. Cambridge Philos. Soc., Lemma 5.4 and Lemma 2.4. [arxiv.org/abs/2109.06665](https://arxiv.org/abs/2109.06665)
- Montgomery and Vaughan 2007, Multiplicative Number Theory I: Classical Theory, Cambridge Studies in Advanced Mathematics 97, Corollary 13.16 and Theorem 13.22. [cambridge.org](https://www.cambridge.org/core/books/multiplicative-number-theory-i/4E45519B26115AEEA4839C6C38206ACD)
- Titchmarsh 1986, The Theory of the Riemann Zeta-Function, second edition revised by D. R. Heath-Brown, Clarendon Press, chapter 9. [sites.math.rutgers.edu](https://sites.math.rutgers.edu/~zeilberg/EM18/TitchmarshZeta.pdf)
- Mertens 1874, Ueber einige asymptotische Gesetze der Zahlentheorie, J. reine angew. Math. 77, 289-338. [doi.org/10.1515/crll.1874.77.289](https://doi.org/10.1515/crll.1874.77.289)
