---
title: Along the squares
lead: A digit-parity design read along the squares turns random at a rate its spectrum fixes and with a strength it does not. The 30 statistics that light most fit N^(-2 eta) within 0.05 over 2^18..2^30, eta the Fourier sup exponent at their orbit, certified to 7.8e-5 for all 57 designs of base at most 6. How strongly a statistic lights is a coefficient that no spectral invariant decides: a design and its mirror have one spectral measure, yet along n^2 the mirror pair 6:5, 6:12345 reads a chi-square of 13152 against 714 at N = 2^30. At N = 2^30 no design with eta > 1/4 sits above the noise floor while 23 of the 30 below it do, the lit lags pass through orbits with eta < 1/4, 43 of 43, and the bias sees only orbits with a nonzero quadratic Gauss sum. Along the consecutive primes the lags match each design's own autocorrelation at the prime gaps, to the noise floor.
figure: research-random-squares
---

A design is a digit rule carried through every level by the Kronecker product, and read in order it is automatic ([the design ladder](design-ladder.md)). Read along the squares it turns into noise: Thue-Morse `t(n)` shows 22 of the 256 blocks of 8 bits in order and all 256 along `n^2`, and `t(n^2)` is normal (Drmota, Mauduit and Rivat 2019, [refs](../refs.md)). This chapter asks how fast, for every design of the simplest kind, `f_{b,D}(n)`, the parity of the number of base-`b` digits of `n` lying in `D`, with `2 <= b <= 6` and `D` a nonempty subset of `{1..b-1}`. There are 57 of them. Each is read along `n^2` for `n < 2^30` and along the first `2^26` primes, through its bias, its lags `h <= 7` and the chi-square of its 8-blocks, beside Rudin-Shapiro and a coin.

The literature settles which of these designs are normal along the squares, and it settles no speed. Every even-base design is normal along `n^2` by Mullner 2018, Theorem 1.4. The odd-base designs fall outside its hypothesis, and for them only the letters are proved uniform (Drmota and Morgenbesser 2012). The proofs give `N^(1 - eta')` with `eta' > 0` existential (Drmota, Mauduit and Rivat 2019, Theorem 2; Mullner 2018, Theorem 1.6), and they run on one parameter of the design, the exponent `gamma` of its Fourier property, with a carry constant (Mauduit and Rivat 2015, Definitions 1 and 2; Adamczewski, Drmota and Mullner 2022, Definitions 7.1 and 7.2). A ranking of designs by their speed along the squares was not found in the sources read. The natural guess is that the Fourier exponent ranks it.

Half of that guess holds up. `eta` is certified here to within `7.8e-5` for all 57 designs, and the 30 statistics that light most strongly at `N = 2^30` fit `N^(-2 eta)` within `0.05` over `2^18..2^30`, Thue-Morse's bias at `-0.403` against `-log2(4/3) = -0.415`, the local dimension of its spectral measure at `1/3`, the same peak that [the readings](readings.md) find leaking into the free readings of Thue-Morse. The other half fails, and no spectral invariant can rescue it: the level of a statistic, and whether it lights at all, is a coefficient the spectrum does not fix. A design and its mirror, the tile read backwards, have the same autocorrelation at every lag, hence one spectral measure, yet along `n^2` the mirror pair `6:5`, `6:12345` reads `X2 = 13152` against `714`. What `eta` says about the coefficient is a floor, seen at one checkpoint: at `N = 2^30` no design with `eta > 1/4` is above the noise floor, 24 of 24, while 23 of the 30 below `1/4` are. Every lag that lights up passes through an orbit with `eta < 1/4`, 43 of 43, and the bias rises only on orbits whose quadratic Gauss sum is nonzero, which splits the twins `5:14` and `5:34` at every checkpoint. Along the consecutive primes the picture is different again: the lags there match each design's own autocorrelation, sampled at the prime gaps and conditioned on the last digit, to the noise floor over 171 cells.

The generator is the study `lab/rs/squares`, verbs `census`, `eta`, `diagonal`, `speed` and `thue`. The figure is a `64 x 64` square read row by row, `n = 64 i + j`. Above the diagonal it draws Thue-Morse in order, a plaid since `t(64 i + j) = t(i) t(j)`; below it draws Thue-Morse along the squares, `t(n^2)`, which is noise. Blue and orange mark the ones. The figure checks that the 4096 cells meet 22 and 256 distinct 8-blocks, and that `sum_(n<4096) (-1)^(t(n^2)) = -260`, which the census prints as `z_B = -4.1` at `2^12`.

## The objects

A design `b:D` has the sign `x(n) = (-1)^(f_{b,D}(n)) = prod_i v(n_i)` over the base-`b` digits `n_i` of `n`, with the tile `v(d) = (-1)^[d in D]` and `v(0) = 1`. So `2:1` is Thue-Morse, `3:2` is the ternary parity `t3` of [the moire limit](two-bases.md), and `6:135` counts the odd digits in base 6. Every such `x` is a signed design in the sense of the ladder, strongly `b`-multiplicative: `x(a b^k + r) = x(a) x(r)` for `r < b^k`.

A reading is a sequence of integers fed to the design: the squares give `y(n) = x(n^2)`, the primes `y(i) = x(p_i)` with `p_0 = 2`. A reading of length `N` is cut into its `N` overlapping 8-blocks `(y(n), ..., y(n + 7))`, `n < N`, and the 256 block counts are Walsh transformed. That gives every sum `T_S = sum_(n<N) prod_(i in S) y(n + i)`, `S` a nonempty subset of `{0..7}`, and three statistics:

- the bias `z_B = T_{0}/sqrt(N)`, so `z_B = S(N)/sqrt(N)` with `S(N) = sum_(n<N) y(n)`;
- the lags `z_Ch = T_{0,h}/sqrt(N)` for `h = 1..7`;
- the chi-square `X2 = sum_(S) T_S^2/N`, which is exactly Pearson's statistic of the 256 counts. It splits by the order `abs(S)` into `X2(1)`, `X2(2)` and `X2(3+)`.

**For an iid coin `X2` has mean exactly 255 and standard deviation tending to 38.7. Proved.** Each `T_S/sqrt(N)` has mean 0 and variance 1, since two products over distinct translates of `S` are uncorrelated. The 255 sets fall into 128 classes under translation, the class with span `s` holding `8 - s` sets whose sums agree up to boundary terms, so `X2` tends to `sum_c (8 - s_c) Z_c^2` with independent standard normals `Z_c`. Its variance is `2 sum_c (8 - s_c)^2 = 2 . 749`, and its square root is `38.7`. The census calls a value above the mean plus 4 sd, `409.8`, above the floor. A `z` beyond about 4 is out of the noise for an iid stream; for a deterministic stream it is a structure, not a probability.

The Fourier side of a design is the trigonometric polynomial of its tile, `phi(alpha) = b^-1 sum_(d<b) v(d) e(d alpha)`, `e(x) = exp(2 pi i x)`. Through the product formula `sum_(m < b^L) x(m) e(m alpha) = b^L prod_(j<L) phi(b^j alpha)`, `abs(phi(alpha))^2 = b^-2 (b + 2 sum_(s=1)^(b-1) a_s cos(2 pi s alpha))` depends only on the tile autocorrelation `a_s = sum_r v(r) v(r + s)`. Two invariants follow.

- **The Fourier sup exponent** `eta = -lim (1/L) log_b sup_alpha prod_(j<L) abs(phi(b^j alpha))`, the limit existing by submultiplicativity. It is the best exponent `gamma` of the Fourier property in the sources above. Parseval gives `eta <= 1/2`, and Rudin-Shapiro reaches `1/2`. For Thue-Morse `eta = log2(2/sqrt 3) = 0.2075`, at `alpha = 1/3`, proved below. The two-step bound of Mauduit and Rivat 2010, Lemmas 20 and 21, `q^gamma = max_t sqrt(phi_q(alpha + t) phi_q(alpha + q t))`, does not reach it: at `q = 2`, `alpha = 1/2` its maximum is `1.75477`, at `t = 0.3041` and not at `1/3`, which gives only `eta >= 1 - gamma = 0.18872`. A periodic point `alpha = a/(b^k - 1)` has its own exponent `eta(alpha) = -(1/k) log_b prod_(j<k) abs(phi(b^j alpha))`, and `eta <= eta(alpha)`.
- **The energy exponent** `D2 = 3 - log_b rho(T)`. The aperiodic autocorrelation `A_L(h)` of the first `b^L` terms obeys `A_(L+1)(bq + s) = al_s A_L(q) + al_(b-s) A_L(q + 1)`, with `al_s = sum_(r < b-s) v(r) v(r + s)`, so the pair `(sum_h A_L(h)^2, sum_h A_L(h) A_L(h + 1))` moves by an integer `2 x 2` matrix `T` from `L` to `L + 1`. Since `sum_h A_L(h)^2 = integral abs(F_L)^4` with `F_L` the polynomial of the first `b^L` terms, `D2` is the correlation dimension of the spectral measure read through its Riesz-product approximants. For Thue-Morse `T = [[6, 2], [-4, -4]]` and `D2 = 3 - log2(1 + sqrt 17) = 0.64298`. [The moire limit](two-bases.md) measured the Cesaro means of `C_t^2` falling at `0.6324` per bit over `2^9` to `2^25`.

Four relations sort the 57 designs, each printed by `census`.

- **Copies.** `4:12` is Thue-Morse, since the base-4 tile `(1, -1, -1, 1)` is `(-1)^(a + c)` on the digit `2a + c`. `5:13` is `3:1`. That leaves 55 distinct sequences.
- **The folded class.** At odd `b` the design of the odd digits is `n mod 2`, since `b^j` is odd and `n = sum_j n_j b^j = sum_j n_j (mod 2)`. Its 8-blocks are two, and it is a class of its own, `3:1` and `5:13`.
- **Twins.** **At odd `b`, `x_D(n) = (-1)^n x_{D'}(n)` with `D'` the symmetric difference of `D` and the odd digits, so `C_h(D) = (-1)^h C_h(D')` along the squares at every `N`, and `x_D(p) = -x_{D'}(p)` at every odd prime. Proved.** The tile of `D'` is `v(d) (-1)^d`, and `n` and its digit sum agree mod 2. Along the squares, `(-1)^(n^2) = (-1)^n`, so `y_D(n) y_D(n + h) = (-1)^h y_{D'}(n) y_{D'}(n + h)`. There are 8 twin pairs, at `b = 3, 5`.
- **Mirrors.** **The mirror `v*(d) = v(b-1) v(b-1-d)` satisfies `x*(m) = v(b-1)^L x(b^L - 1 - m)` for `m < b^L`, has the same tile autocorrelation, and so has the same correlation sequence and the same spectral measure. Proved.** The digits of `b^L - 1 - m` are `b - 1 - m_i`, with no borrow, which gives the first identity. A block and its reversal have the same aperiodic autocorrelation, so `a*_s = a_s`, `abs(phi*) = abs(phi)` pointwise, and `A*_L(h) = A_L(h)` at every level. Then `C*(h) = lim A*_L(h)/b^L = C(h)`, and the spectral measure is the one with these Fourier coefficients (Herglotz). There are 21 mirror pairs and 15 self-mirror designs; Thue-Morse, `4:13`, `6:14`, `6:23` and `6:135` are their own mirrors.

## What the sources prove

The census prints a status for each design, read from [refs](../refs.md).

- **Proved normal along `n^2`: the 39 even-base designs.** Mullner 2018, Theorem 1.4, makes `b(n^2) mod m'` normal for every strongly block-additive `b` in base `q` once `gcd(q - 1, m') = 1` and `gcd(m', gcd{b(n)}) = 1`. Here `q = b`, `m' = 2` and `b(n) = f_{b,D}(n)` with blocks of length 1; `b - 1` is odd exactly when `b` is even, and `f` takes the value 1 on any digit of `D`. Thue-Morse is Drmota, Mauduit and Rivat 2019, and Rudin-Shapiro is Mullner 2018 and Mauduit and Rivat 2018.
- **Open: the 16 non-folded odd-base designs.** At odd `b`, `gcd(b - 1, 2) = 2` and Theorem 1.4 does not apply; the folded class shows the hypothesis is needed. For these designs only the letters are proved uniform along `n^2`: Drmota and Morgenbesser 2012, Theorem 1, distributes `T(n^2)` on cosets of a subgroup of index `m` dividing `q - 1` with the weights `Q(u, m) = #{n mod m: n^2 = u}`, uniform for `m <= 2` (their Remark 1). The census verdict on their blocks is evidence for a conjecture, not a proof. The same paper asserts that digit sums mod `m > 2` are not uniform along the squares when `m` and `q - 1` share a factor beyond 2, as `s_q(n^2) = n^2 (mod q - 1)` makes plain; this does not bear on the digit counts here, whose integer `m` in Theorem 1 is at most 2.
- **Periodic: the folded class**, `n mod 2` along every reading.
- **No rate.** The exponents in Drmota, Mauduit and Rivat 2019, Theorem 2, and in Mullner 2018, Theorem 1.6, are existential. The proofs need the Fourier property, a bound `abs(b^-L sum_(u < b^L) x(u) e(-u alpha)) <= c b^(-gamma L)` uniform in `alpha`, and a carry property (Mauduit and Rivat 2015). The best `gamma` is the `eta` of this chapter, known in closed form in the sources only for digit sums and digit counts.
- **Along the primes.** The letters of every non-folded design are equidistributed along the primes. Martin, Mauduit and Rivat 2014, Theorem 2, gives a power saving for a digit count whose coefficients are not an arithmetic progression mod 1, which for a parity design excludes only the odd digits, and their Theorem 1 covers the odd digits at even `b`; at odd `b` the odd digits are the folded class, constant on the odd primes. Mauduit and Rivat 2010 is the digit sum. Normality of Thue-Morse along the primes is open (Popoli 2023, Section 1.4; Merai and Winterhof 2022, Section 7), so the primes column below is a measurement and no comparison with a theorem.
- **The block length.** Hanna 2016, Proposition 11.1, shows the Fourier property fails for blocks of length `k >= A log2 N`, `A > 1`, a threshold on the block-length axis. This chapter reads blocks of 8 only and does not test it.

## The census

`census 30 26` reads every design along `n^2`, `n < 2^30`, and along the first `2^26` primes, up to `p = 1339484197`, in about 15 s on 8 threads. The table gives the status, the certified bracket on `eta` from `eta`, printed outward, with its orbit, the energy exponent `D2`, the chi-square and bias along the squares at `N = 2^30`, the lags `h <= 7` with `abs(z_Ch) > 4`, and the bias along the primes at `N = 2^26`.

| design | class | normal | `eta` | orbit | `D2` | `X2` | `z_B` | lit | `z_B` primes |
|---|---|---|---|---|---|---|---|---|---|
| `2:1` |  | proved | `[0.207510, 0.207519]` | `1/3` | 0.64298 | 3545 | -15.63 | 3,6 | -89.7 |
| `3:1` | folded | periodic | `[0.000000, 0.000001]` | `1/2` | 0.00000 | 1.36e+11 | 0.00 | 1,2,3,4,5,6,7 | -8192.0 |
| `3:2` | twin | open | `[0.267493, 0.267514]` | `1/4` | 0.81734 | 268 | -0.14 | - | 0.3 |
| `3:12` | twin | open | `[0.267493, 0.267514]` | `1/4` | 0.81734 | 282 | -1.41 | - | -0.3 |
| `4:1` |  | proved | `[0.298140, 0.298162]` | `1/3` | 0.76226 | 372 | -2.37 | - | 4.9 |
| `4:2` |  | proved | `[0.298140, 0.298162]` | `1/3` | 0.76226 | 236 | 0.42 | - | -1.8 |
| `4:12` | copy | proved | `[0.207491, 0.207519]` | `1/3` | 0.64298 | 3545 | -15.63 | 3,6 | -89.7 |
| `4:3` |  | proved | `[0.339827, 0.339841]` | `1/5` | 0.88343 | 232 | 0.50 | - | -5.7 |
| `4:13` |  | proved | `[0.189034, 0.189078]` | `2/5` | 0.57430 | 9001 | -28.93 | 5 | -46.2 |
| `4:23` |  | proved | `[0.189049, 0.189078]` | `1/5` | 0.57430 | 1722 | -4.84 | 5 | -74.5 |
| `4:123` |  | proved | `[0.339827, 0.339841]` | `1/5` | 0.88343 | 297 | 1.31 | - | -4.0 |
| `5:1` | twin | open | `[0.317377, 0.317394]` | `0/1` | 0.90778 | 327 | -2.23 | - | -2.5 |
| `5:2` | twin | open | `[0.315418, 0.315447]` | `7/26` | 0.75641 | 272 | 2.30 | - | 10.1 |
| `5:12` | twin | open | `[0.203134, 0.203154]` | `1/4` | 0.63479 | 744 | -5.80 | 2 | -0.8 |
| `5:3` | twin | open | `[0.317377, 0.317394]` | `0/1` | 0.90778 | 264 | -0.88 | - | 2.5 |
| `5:13` | copy | periodic | `[0.000000, 0.000001]` | `1/2` | 0.00000 | 1.36e+11 | 0.00 | 1,2,3,4,5,6,7 | -8192.0 |
| `5:23` | twin | open | `[0.203134, 0.203154]` | `1/4` | 0.63479 | 485 | 0.65 | 2 | 0.8 |
| `5:123` | twin | open | `[0.315405, 0.315447]` | `2/13` | 0.75641 | 272 | 0.84 | - | -10.1 |
| `5:4` | twin | open | `[0.317371, 0.317394]` | `0/1` | 0.83515 | 253 | -1.00 | - | 10.2 |
| `5:14` | twin | open | `[0.203128, 0.203154]` | `1/3` | 0.63479 | 2399 | -6.68 | 3,6 | -66.5 |
| `5:24` | twin | open | `[0.317371, 0.317394]` | `1/2` | 0.83515 | 331 | -0.57 | - | 14.7 |
| `5:124` | twin | open | `[0.203128, 0.203154]` | `1/3` | 0.63479 | 3237 | -15.55 | 3,6 | -97.6 |
| `5:34` | twin | open | `[0.203122, 0.203154]` | `1/6` | 0.63479 | 2015 | 0.86 | 3,6 | 66.5 |
| `5:134` | twin | open | `[0.317371, 0.317394]` | `1/2` | 0.83515 | 301 | -0.52 | - | -10.2 |
| `5:234` | twin | open | `[0.203122, 0.203154]` | `1/6` | 0.63479 | 1398 | -2.35 | 3,6 | 97.6 |
| `5:1234` | twin | open | `[0.317371, 0.317394]` | `0/1` | 0.83515 | 273 | 0.14 | - | -14.7 |
| `6:1` |  | proved | `[0.226263, 0.226295]` | `0/1` | 0.77812 | 1007 | -2.65 | 1,2 | 0.2 |
| `6:2` |  | proved | `[0.226254, 0.226295]` | `0/1` | 0.77812 | 762 | 1.91 | 1 | 63.6 |
| `6:12` |  | proved | `[0.274246, 0.274283]` | `4/5` | 0.79683 | 298 | -0.83 | - | -9.1 |
| `6:3` |  | proved | `[0.226254, 0.226295]` | `0/1` | 0.77812 | 1766 | 8.92 | 1,2,4,5,6 | 95.0 |
| `6:13` |  | proved | `[0.233468, 0.233507]` | `3/7` | 0.59529 | 622 | 5.58 | 7 | -29.9 |
| `6:23` |  | proved | `[0.194227, 0.194292]` | `4/5` | 0.57411 | 271 | 0.14 | - | -46.7 |
| `6:123` |  | proved | `[0.214170, 0.214209]` | `1/5` | 0.59689 | 381 | 2.41 | - | -41.8 |
| `6:4` |  | proved | `[0.226263, 0.226295]` | `0/1` | 0.77812 | 1164 | 6.61 | 1,2 | 152.6 |
| `6:14` |  | proved | `[0.194244, 0.194292]` | `2/5` | 0.57411 | 6184 | -24.05 | 5 | -73.4 |
| `6:24` |  | proved | `[0.233468, 0.233507]` | `3/7` | 0.59529 | 390 | -1.23 | - | -20.4 |
| `6:124` |  | proved | `[0.214168, 0.214209]` | `2/5` | 0.59689 | 479 | 4.13 | 5 | -20.4 |
| `6:34` |  | proved | `[0.274246, 0.274283]` | `4/5` | 0.79683 | 298 | -1.80 | - | -8.3 |
| `6:134` |  | proved | `[0.214168, 0.214209]` | `2/5` | 0.59689 | 399 | 2.68 | 5 | -72.9 |
| `6:234` |  | proved | `[0.214170, 0.214209]` | `1/5` | 0.59689 | 438 | -1.29 | 5 | -53.0 |
| `6:1234` |  | proved | `[0.301727, 0.301784]` | `1/7` | 0.79683 | 238 | -1.13 | - | -5.5 |
| `6:5` |  | proved | `[0.226259, 0.226295]` | `0/1` | 0.72771 | 13152 | 21.28 | 1,2,3,4,5,6,7 | -0.9 |
| `6:15` |  | proved | `[0.333283, 0.333346]` | `9/35` | 0.83664 | 293 | -1.08 | - | -4.6 |
| `6:25` |  | proved | `[0.233451, 0.233507]` | `2/7` | 0.70993 | 431 | 3.84 | - | -11.0 |
| `6:125` |  | proved | `[0.175434, 0.175480]` | `2/7` | 0.51913 | 870 | -0.53 | 7 | -76.0 |
| `6:35` |  | proved | `[0.324756, 0.324807]` | `3/35` | 0.83664 | 248 | -0.03 | - | -1.1 |
| `6:135` |  | proved | `[0.175402, 0.175480]` | `3/7` | 0.51913 | 7255 | 27.25 | 7 | -52.7 |
| `6:235` |  | proved | `[0.282103, 0.282146]` | `2/7` | 0.83003 | 256 | -1.13 | - | -5.8 |
| `6:1235` |  | proved | `[0.333283, 0.333346]` | `9/35` | 0.83664 | 270 | -1.83 | - | -1.9 |
| `6:45` |  | proved | `[0.233476, 0.233507]` | `1/7` | 0.70993 | 364 | 3.08 | 7 | -27.2 |
| `6:145` |  | proved | `[0.282103, 0.282146]` | `2/7` | 0.83003 | 327 | 1.64 | - | -9.8 |
| `6:245` |  | proved | `[0.341251, 0.341313]` | `4/37` | 0.83003 | 241 | -0.12 | - | -1.7 |
| `6:1245` |  | proved | `[0.233451, 0.233507]` | `2/7` | 0.70993 | 319 | 3.13 | - | -26.2 |
| `6:345` |  | proved | `[0.175440, 0.175480]` | `1/7` | 0.51913 | 2354 | 7.04 | 7 | -74.3 |
| `6:1345` |  | proved | `[0.324756, 0.324807]` | `3/35` | 0.83664 | 208 | -0.89 | - | -5.5 |
| `6:2345` |  | proved | `[0.233476, 0.233507]` | `1/7` | 0.70993 | 231 | -0.95 | - | -12.1 |
| `6:12345` |  | proved | `[0.226259, 0.226295]` | `0/1` | 0.72771 | 714 | -1.02 | 1,2,3 | 52.6 |
| Rudin-Shapiro | | proved | `0.5` | | | 343 | 3.37 | - | -2.4 |
| coin | |  |  | | | 289 | 0.56 | - | 1.5 |

**Along the squares every non-folded design meets all 256 8-blocks below `2^30`, against 18 to 28 along `n < 2^20`. Verified.** Thue-Morse meets 22 in order, Rudin-Shapiro 56 and the coin 256; the folded class meets 2 either way. The coin reads `X2 = 289`, `0.9` sd above the mean, and Rudin-Shapiro `343`, with a bias `z_B = 3.37`.

**The slow designs grow, and the fast ones stay near the floor without staying under it. Verified.** The growth table of `census` prints `z_B` and `X2` at every second power of 2. Thue-Morse climbs from `X2 = 502` at `2^12` to `3545` at `2^30`, and `z_B` from `-4.1` to `-15.6`; `4:13` climbs from `774` to `9001`, `6:5` from `2714` to `13152`. The coin stays between `206` and `373`. A design just above `eta = 1/4` can sit over the floor for a long while: `4:1`, `eta = 0.2982`, reads above `409.8` at 18 of the 19 checkpoints from `2^12` to `2^30`, peaks at `847.7` at `2^20`, and is down to `372.3` at `2^30`; `3:2` peaks at `522.4` at `2^25` and ends at `267.8`.

## The exponent, certified

**For each of the 57 designs `eta` lies in a bracket of width at most `7.80e-5`. Proved.** The computation is finite and runs in interval arithmetic, every float rounded outward after each operation, with cosines by a Taylor series on `[0, pi/4]` and logarithms by the `atanh` series. `eta` prints the brackets in about 6 s, outward: the lower end floored and the upper end ceiled at 6 places.

- **The lower end.** Cut the circle into the `b^r` cylinders `I_i = [i b^-r, (i + 1) b^-r)` with `b^r <= 2^17`, so `r = 17, 10, 8, 7, 6` for `b = 2..6`. Multiplication by `b` maps `I_i` onto the `b` cylinders `bi + d mod b^r`. Bound `abs(phi)^2` on `I_i` by its value at the centre plus `2 pi sum_s s abs(a_s)/(b^2 b^r)`, a Lipschitz term, and call the bound `G_i`. For any positive vector `w`, let `c = max_i G_i max_j w_j/w_i`, `j` over the successors of `i`. Then along any orbit, with `i_t` the cylinder of `b^t alpha`, `prod_(t<L) abs(phi(b^t alpha))^2 <= prod_t G_(i_t) <= c^L w_(i_0)/w_(i_L)`, so the growth rate is at most `c` and `eta >= -(1/2) log_b c`. Howard's policy iteration for the maximal cycle mean chooses `w`; it only needs to be positive, and the bound is exact whatever it is.
- **The upper end.** A periodic point `alpha = a/(b^k - 1)` gives `sup_alpha prod_(j<L) abs(phi(b^j alpha)) >= prod_(j<L) abs(phi(b^j a/(b^k - 1)))`, so `eta <= eta(alpha)`, with the product taken in intervals. The candidates are every orbit up to period `log_b 2^20` and the optimal cycle of the certificate.

The ends meet on short orbits. Thue-Morse sits in `[0.207510, 0.207519]`, which holds `log2(2/sqrt 3) = 0.2075187` at `1/3`. `3:2` sits in `[0.267493, 0.267514]`, which holds `log3(3/sqrt 5) = 0.2675132` at `1/4`. The five singletons `6:d` end at `0.226295`, just above `log6(3/2) = 0.2262944`, the fixed point `0`, where `abs(phi(0)) = 4/6`. `5:12` and `5:14` end at `0.203154`, above `log5(5/sqrt 13) = 0.2031537`, on `1/4` and on `1/3`. `6:14` and `6:23` end at `0.194292`, above `log6(6/(2 + sqrt 5)) = 0.1942917`, on the fixed points `2/5` and `4/5`, where both polynomials take the value `2 + (zeta + zeta^4) - (zeta^2 + zeta^3) = 2 + sqrt 5`, `zeta = e(1/5)`; their lower ends differ, `0.194244` and `0.194227`. The least exponents are `6:125`, `6:135` and `6:345`, which end at `0.175480` on period-2 orbits of sevenths, and the largest is `6:245`, in `[0.341251, 0.341313]` on a period-4 orbit of `4/37`.

**For `2:1` and its copy `4:12`, `3:2`, `3:12`, `5:12`, `5:23`, `5:1` and `5:3`, `eta` is exactly the value on the printed orbit. Proved.**

- Thue-Morse. With `s(a) = abs(sin pi a)`, `abs(phi(a)) = s(a)` and `s(2a) = 2 s(a) abs(cos pi a)`, so `s(a)^2 s(2a) = 2 abs(sin^3 x cos x)` with `x = pi a`. Its maximum is at `tan^2 x = 3`, `x = pi/3`, where it is `(sqrt 3/2)^3`: `s(a)^(2/3) s(2a)^(1/3) <= sqrt 3/2`, with equality at `a = 1/3`. Multiplying over `j = 0..L-2` with `u_j = s(2^j a)` gives `u_0^(2/3) u_1 ... u_(L-2) u_(L-1)^(1/3) <= (sqrt 3/2)^(L-1)`, and since `u_0, u_(L-1) <= 1`, `prod_(j<L) u_j <= (sqrt 3/2)^(L-1)`. The orbit `{1/3, 2/3}` reaches `(sqrt 3/2)^L`, so `eta = log2(2/sqrt 3)`. This is a sub-action: a weight that turns the `L`-step product into a one-step bound.
- `3:2` and `3:12`. Here `a_1 = 0` and `a_2 = -1`, so `abs(phi)^2 = (3 - 2 cos 4 pi alpha)/9 <= 5/9`, with equality on `{1/4, 3/4}`, a cycle of multiplication by 3; the twin shifts `alpha` by `1/2`, which leaves `cos 4 pi alpha`. So `eta = log3(3/sqrt 5)`.
- `5:12` and `5:23`. Here `a_1 = a_3 = 0`, `a_2 = -3`, `a_4 = 1`, so with `c = cos 4 pi alpha`, `abs(phi)^2 = (3 - 6c + 4c^2)/25`, largest at `c = -1`, `13/25`, on the fixed points `1/4` and `3/4`. So `eta = log5(5/sqrt 13)`.
- `5:1` and `5:3`. Here `a_1 = a_3 = 0`, `a_2 = a_4 = 1`, so `abs(phi)^2 = (3 + 2c + 4c^2)/25`, largest at `c = 1`, `9/25`, on the fixed points `0` and `1/2`. So `eta = log5(5/3)`.

In each one-step case the bound `sup prod_(j<L) abs(phi(b^j alpha)) <= (max abs(phi))^L` is reached on an invariant set. **That the upper end is the exact `eta` for the other designs, on the printed orbit, is a Conjecture.**

The folded class has `eta = 0`, at the fixed point `1/2` of odd `b`, where `abs(phi(1/2)) = 1`. No other design reaches `0`: `abs(phi) = 1` forces the tile to be `(-1)^d`, and at even `b` the point `1/2` falls to `0`, where `phi` vanishes. Of the 54 distinct non-folded designs, 24 have `eta > 1/4` and 30 have `eta < 1/4`. Thue-Morse is not the least: `4:13` and `4:23` at `0.189078`, `6:14` and `6:23` at `0.194292`, and the three sevenths designs at `0.175480` sit below its `0.2075`.

## The spectrum does not fix the level

**One spectral invariant of a design ranks how far from random it reads along the squares. Refuted.** By the mirror theorem, `6:5` and its mirror `6:12345` have the same autocorrelation at every lag and the same spectral measure, so they agree on every invariant of the measure: `eta`, `D2`, the decay of the correlations, the local dimensions. Along `n^2`, `n < 2^30`, they read

| | `X2` | `z_B` | `z_C1` | `z_C2` | lit lags |
|---|---|---|---|---|---|
| `6:5` | 13152 | 21.28 | 24.42 | 21.80 | 1 to 7 |
| `6:12345` | 714 | -1.02 | 5.67 | -5.37 | 1, 2, 3 |

so one reads 18.4 times the chi-square of the other. The pair refutes a shared level, not a shared rate: the bias of `6:5` decays at a fitted `-0.447` against `-2 eta = -0.453`, while that of `6:12345` never lights, so there is no second rate to compare. The mirror maps `m` to `b^L - 1 - m`, which keeps the integers and their spectrum but not the squares. The squares see the phases of the Fourier coefficients, which the spectral measure forgets.

The pair `6:14`, `6:23` makes the same point without a symmetry. Neither is the other's mirror, both being palindromes, and their autocorrelations differ (`C(1) = -3/5` against `1/5`). Yet they share the upper end of their `eta` brackets, `0.1942917` on both orbits, with lower ends `0.194244` and `0.194227`, and `D2 = 0.57411`, their energy matrices `[[82, 6], [-52, 12]]` and `[[82, -26], [12, 12]]` both having the characteristic polynomial `x^2 - 94x + 1296`. Along the squares they read `X2 = 6184` against `271`, `z_B = -24.05` against `0.14`, and `z_C5 = 19.78` against `-0.52`.

As a rank statistic over the 54 distinct non-folded designs, Spearman's correlation of `X2` at `2^30` is `-0.799` with `eta` and `-0.656` with `D2`, and `eta` against `D2` is `0.900`. Of the 855 pairs whose `eta` differ by more than `0.01`, with one of them above the floor, 50 are ranked backwards by a factor 2 in `X2`. `eta` explains most of the order and decides none of it.

## What the exponent decides

### The rate

**The 30 statistics with `abs(z) >= 8` at `N = 2^30`, in 15 designs, have fitted decay exponents within `0.05` of the prediction over `2^18..2^30`, mean gap `+0.011` and mean absolute gap `0.017`. Verified.** `speed 30` fits the least-squares slope of `log2 abs(T_S/N)` over every `N = 2^k`, `k = 18..30`, for the bias and for each lit lag, against `-2 eta` for the bias and `-2 e_h` for the lag `h`, `e_h` the best orbit through the points `k/(2h)`. The noise alone has `-0.5`. Over the late window `2^24..2^30` the mean absolute gap is `0.023`.

| design | `eta` | statistic | `z` at `2^30` | fit | predicted | gap |
|---|---|---|---|---|---|---|
| `2:1` | 0.2075 | `B` | -15.63 | -0.403 | -0.415 | 0.012 |
| `2:1` | 0.2075 | `C3` | 14.17 | -0.377 | -0.415 | 0.038 |
| `2:1` | 0.2075 | `C6` | 10.82 | -0.371 | -0.415 | 0.044 |
| `4:13` | 0.1891 | `B` | -28.93 | -0.403 | -0.378 | -0.025 |
| `4:13` | 0.1891 | `C5` | 25.70 | -0.378 | -0.378 | 0.000 |
| `4:23` | 0.1891 | `C5` | -20.04 | -0.360 | -0.378 | 0.018 |
| `5:14` | 0.2032 | `C3` | 17.81 | -0.397 | -0.406 | 0.010 |
| `5:14` | 0.2032 | `C6` | 9.11 | -0.401 | -0.406 | 0.005 |
| `5:124` | 0.2032 | `B` | -15.55 | -0.390 | -0.406 | 0.016 |
| `5:124` | 0.2032 | `C3` | 13.81 | -0.387 | -0.406 | 0.020 |
| `5:34` | 0.2032 | `C3` | -17.81 | -0.397 | -0.406 | 0.010 |
| `5:34` | 0.2032 | `C6` | 9.11 | -0.401 | -0.406 | 0.005 |
| `5:234` | 0.2032 | `C3` | -13.81 | -0.387 | -0.406 | 0.020 |
| `6:1` | 0.2263 | `C1` | -8.88 | -0.448 | -0.453 | 0.004 |
| `6:3` | 0.2263 | `B` | 8.92 | -0.465 | -0.453 | -0.013 |
| `6:13` | 0.2335 | `C7` | 8.14 | -0.489 | -0.467 | -0.022 |
| `6:4` | 0.2263 | `C2` | -8.01 | -0.434 | -0.453 | 0.018 |
| `6:14` | 0.1943 | `B` | -24.05 | -0.414 | -0.389 | -0.025 |
| `6:14` | 0.1943 | `C5` | 19.78 | -0.383 | -0.389 | 0.006 |
| `6:5` | 0.2263 | `B` | 21.28 | -0.447 | -0.453 | 0.005 |
| `6:5` | 0.2263 | `C1` | 24.42 | -0.437 | -0.453 | 0.016 |
| `6:5` | 0.2263 | `C2` | 21.80 | -0.409 | -0.453 | 0.044 |
| `6:5` | 0.2263 | `C3` | 13.48 | -0.428 | -0.453 | 0.024 |
| `6:5` | 0.2263 | `C4` | 12.27 | -0.427 | -0.453 | 0.026 |
| `6:5` | 0.2263 | `C5` | 11.62 | -0.440 | -0.453 | 0.013 |
| `6:5` | 0.2263 | `C6` | 11.72 | -0.432 | -0.453 | 0.020 |
| `6:5` | 0.2263 | `C7` | 9.68 | -0.430 | -0.453 | 0.023 |
| `6:135` | 0.1755 | `B` | 27.25 | -0.337 | -0.351 | 0.014 |
| `6:135` | 0.1755 | `C7` | 30.21 | -0.367 | -0.351 | -0.016 |
| `6:345` | 0.1755 | `C7` | 40.35 | -0.338 | -0.351 | 0.013 |

The exponent is the local dimension of the spectral measure at the orbit, `r^(2 eta)` for the mass within `r`. For Thue-Morse that is `log2(4/3)` at `1/3`. [The readings](readings.md) prove the mean of `C_t` over the multiples of 3 below `2^k` to be `(16/15)(2^k)^(-log2(4/3)) + O(2^(-k))`. They measure the dilation correlation of `t` along `n^2 + n` at `(4, 0)` falling with slope `-0.399`, and the shift correlation of `t(n^2)` at `h = 15` with `-0.384`; this chapter's `C3` and `C6` read `-0.377` and `-0.371`, and the bias `-0.403`. Each is the one third peak read through a different window. The fits are finite-window slopes and they move with the window. `thue 40` reads Thue-Morse alone, in 128-bit squares, to `N = 2^40`: its bias decays at `-0.4330` over `2^18..2^40` and `-0.4407` over `2^30..2^40`, against `-0.4150`, and its lags `C3`, `C6` at `-0.3855`, `-0.3799` over `2^18..2^40` and `-0.4022`, `-0.4006` over `2^30..2^40`. An independent recomputation reads `6:5` at `-0.4547` over `2^18..2^38` against `-0.4526`. The slopes are evidence that the exponent of a lit statistic is set by `eta`, which stays a Conjecture, and not a measurement of it; the mirror pair shows that its level is not set by the spectrum.

### A floor at one quarter

**At the single checkpoint `N = 2^30` every distinct non-folded design with `eta > 1/4` reads under the floor, 24 of 24, the largest `X2 = 372.3` at `4:1`, while 23 of the 30 with `eta < 1/4` read above it. Verified.** The split belongs to that checkpoint. Over the 19 checkpoints from `2^12` to `2^30` five of the 24 cross the floor, all close to `1/4`: `4:1` at 18 of them, with a peak of `847.7` at `2^20` and the last at `2^29`; `3:2` at 6, peak `522.4` at `2^25`; `3:12`, `4:123` and `6:12` once each. At `2^30` six designs sit within 1 sd of the cut: `4:1` at `372.3`, `6:123` at `380.7`, `6:24` at `390.2`, `6:134` at `399.0`, `6:25` at `430.7` and `6:234` at `438.0`.

The heuristic behind the quarter is the circle method run on the design rather than on the squares. A square near `m` falls in the class `c mod q` with relative weight `#{n mod q: n^2 = c}`, and expanding those weights over the residues gives the heuristic `1[m is a square] ~ (1/(2 sqrt m)) sum_(a/q) (G(a, q)/q) e(-am/q)`, the sum over reduced fractions with the quadratic Gauss sum `G(a, q) = sum_(n mod q) e(a n^2/q)`. So

`S(N) = sum_(n<N) x(n^2) ~ sum_(a/q) (G(a, q)/q) sum_(m < N^2) x(m) e(-am/q)/(2 sqrt m)`.

On a periodic orbit of multiplication by `b` the twisted sums of the design are exact products. **For `alpha = a/(b^k - 1)`, `abs(sum_(m < b^(kL)) x(m) e(m alpha)) = b^(kL (1 - eta(alpha)))`. Proved**, by the product formula: the `kL` factors run `L` times around the orbit. With `m` up to `N^2`, the term of `alpha` is of size `N^(2 (1 - eta(alpha)))/N = N^(1 - 2 eta(alpha))`. It stands above the noise `sqrt(N)` exactly when `eta(alpha) < 1/4`, and only when `G(a, q) != 0`. **That a lit statistic is `c N^(-2 eta)` up to a bounded log-periodic factor, with a coefficient `c` no spectral invariant fixes, so that a design with `eta > 1/4` is fast at every scale, is a Conjecture**: the circle method here is a heuristic, and the minor arcs are not bounded.

`4:1`, `eta = 0.2982`, shows the other side of the quarter: its structured part falls against the noise as `N^(1/2 - 2 eta) = N^(-0.096)`, a slow decay, and its chi-square falls from `847.7` at `2^20` to `372.3` at `2^30`.

### Which lags light

The lags see the spectral measure where the bias sees Gauss sums. For `n` large, `(n + h)^2 - n^2 = 2hn + h^2` is small beside `n^2`, and averaging `x(n^2) x((n + h)^2)` over `n` samples the design's own correlation along the progression `2hn + h^2`. That average is the mass of the spectral measure within `1/(hN)` of the points `k/(2h)`, and near a point whose orbit is periodic the mass scales as `r^(2 eta(orbit))`. A lag `h` can light only through an orbit with `eta < 1/4` that passes through some `k/(2h)`.

**At `N = 2^30` each of the 43 lags `h <= 7` with `abs(z_Ch) > 4` passes through an orbit with `eta < 1/4`, and none lights elsewhere; 32 lags that pass through one stay quiet. Verified.** The best orbit names the lag:

- Thue-Morse, orbit `1/3`: `3 | 2h` gives `h = 3, 6`, which read `14.17` and `10.82`, while every other lag reads under `2.3`.
- `4:13`, `4:23` and `6:14`, orbits of fifths: `h = 5`, which reads `25.70`, `-20.04` and `19.78`.
- `6:135` and `6:345`, orbits of sevenths: `h = 7`, which reads `30.21` and `40.35`.
- `5:14`, `5:124`, `5:34` and `5:234`, orbits of thirds and sixths: `h = 3, 6`, with `z_C3 = 17.81, 13.81, -17.81, -13.81`, the twin signs `(-1)^h` exact.
- `6:5` and the other singletons, the fixed point `0`, which every `k/(2h)` reaches: `6:5` lights all seven lags.

**The linear part of the lag is not the lag. Verified.** `diagonal 24 22` sums `sum_(n<N) C_x(2hn + h^2)` exactly. `C_x` comes from the `2 x 2` digit matrices `M_s = b^-1 [[al_s, al_(b-s)], [al_(s+1), al_(b-s-1)]]`, `(C(bq + s), C(bq + s + 1)) = M_s (C(q), C(q + 1))`, and a digit walk over `n < N` carries `2hn + h^2`. The model is nonzero where the lags light, but at `N = 2^24` it reads `1.37` against the observed `7.74` for Thue-Morse at `h = 3`, `3.37` against `15.55` for `4:13` at `h = 5`, `5.99` against `22.59` for `6:345` at `h = 7` and `0.92` against `19.61` for `6:5` at `h = 1`; the rms difference over the 399 cells is `2.828`. The rest of each lag comes from the quadratic terms, the Gauss sums again.

### What the bias sees

`G(a, q) = 0` exactly when `q = 2 (mod 4)`. At odd `b` an orbit can carry such a denominator, and then the bias cannot see it while the lags can. The twins at `b = 5` make the test, since twins share `eta` and their lags up to sign. The best orbit of `5:14` is `1/3`, with `G != 0`; that of its twin `5:34` is the shifted orbit `1/6`, with `G(1, 6) = 0`. **Over the 19 checkpoints from `2^12` to `2^30` the largest `abs(z_B)` is `1.80` for `5:34` and `3.18` for `5:234`, whose orbit is `1/6`, against `11.45` for `5:14` and `15.55` for `5:124`, orbit `1/3`, and `5.80` for `5:12` and `7.15` for `5:23`, orbit `1/4` with `G(1, 4) != 0`; the lags of each pair agree, `C_h(5:34) = (-1)^h C_h(5:14)` exactly. Verified.** One checkpoint is no witness: `5:12` and `5:23`, both visible, read `-5.80` and `0.65` at `2^30`, and `5:23` reads its largest, `7.15`, at `2^28`. An independent recomputation to `2^38` keeps `5:34` at `1.80` and `5:234` at `3.18` while `5:14` reaches `20.3`, `5:124` `24.8`, `5:12` `11.0` and `5:23` `7.3`. The folded class is the extreme case: its orbit is `1/2`, with `G(1, 2) = 0`, and its bias is exactly 0 while its lags are `+-1`.

### The quiet seven

Seven designs with `eta < 1/4` read under the floor at the checkpoint `2^30`: `6:23` at `0.1943`, `6:123` and `6:134` at `0.2142`, and `6:24`, `6:45`, `6:1245` and `6:2345` at `0.2335`. For the last four, `1/2 - 2 eta = 0.033`, a growth of `2^0.033` per doubling that `2^30` cannot show. `6:23` is different: its prediction, `0.111`, is that of `6:14`, which reads `X2 = 6184`. Its coefficient on the fifths nearly cancels along the squares, and only there: an independent recomputation finds it above the floor at `2^34..2^37`, `X2` from 420 to 462, all in the orders 3 and up, with its bias and lags still under `abs(z) = 3`. Along the primes `6:23` reads `z_B = -46.7`, with a slope of `0.307` per doubling against the predicted `1/2 - eta = 0.306`. So the cancellation belongs to the squares, not to the design. A design with `eta < 1/4` is slow along the squares unless the arithmetic of the squares cancels its orbit, and `eta` cannot tell which.

## Along the primes

The same census runs along the first `2^26` primes. Two mechanisms show, and both are far from the noise.

**The bias follows `eta` and the last digit.** The primes see the rationals `a/q` with `q` squarefree, through the weight `mu(q)/varphi(q)` of the singular series, `varphi` Euler's function. The twisted sums now run to `m ~ N log N`, not to `N^2`, so the term of an orbit is of size `N^(1 - eta)` against the noise `sqrt(N)`, up to logarithms, and since `eta < 1/2` for every design, a design is slow along the primes whenever its coefficient survives. Summed over the divisors of `b`, the singular series restricts `m` to the residues prime to `b`. The coefficient of the fixed point `0` is then `sum_(d prime to b) v(d)`, which can vanish. **Along the first `2^26` primes Thue-Morse reads `z_B = -89.7`, with a slope of `0.270` per doubling over `2^18..2^26` against `1/2 - eta = 0.292`. At base 6 the units are 1 and 5: `6:1` and `6:5`, which give them opposite signs, read `z_B = 0.2` and `-0.9`, while `6:2`, `6:3` and `6:4` read `63.6`, `95.0` and `152.6`. Verified.** The twins at odd `b` are exact negatives along the odd primes, `5:14` reading `-66.5` and `5:34` `+66.5`. Spearman's correlation of `abs(z_B)` along the primes with `eta` is `-0.626`. Rudin-Shapiro reads `-2.4` and the coin `1.5`.

**The lags match the design's own autocorrelation at the prime gaps.** Consecutive primes share every digit above their gap. So `x(p_i) x(p_(i+h))` reads only the last digits and the carry, and, averaged over the digits above the last, it is exactly `x(r) x((r + g) mod b) C_x(floor((r + g)/b))`, with `r = p_i mod b` and `g = p_(i+h) - p_i`, by strong `b`-multiplicativity. **Summed over the primes, this model matches the observed lags in all 171 cells `h <= 3` at `N = 2^22`, with an rms difference of `0.753` and a largest of `2.17` in the units `z = sum/sqrt(N)`, while 160 of the cells exceed 4 and `6:12345` reads `712.1` at `h = 1`. Verified.** `diagonal 24 22` prints the cells. Rudin-Shapiro, whose autocorrelation vanishes at every nonzero lag, reads `z_C1 = -0.1` along `2^26` primes, where its excess sits in the orders 3 and up, `X2(3+) = 2668` of `2731`. Thue-Morse there reads `z_C2 = 403.3` and `X2 = 2481818`. [The readings](readings.md) find Thue-Morse and the ternary parity failing the shift test along the consecutive primes at about `0.05`, nearly flat in `N`, while Rudin-Shapiro passes; the gap model is that observation, cell by cell.

The lags along the consecutive primes therefore fall only as the gaps spread, and the gaps grow like `log N`. If Thue-Morse is normal along the consecutive primes, an open question, the rate at which its lags vanish is governed by the Cesaro means of `C_t` over gaps near `log N`: a power of `1/log N`, not of `1/N`. That is a heuristic, read off the verified model.

## Doors

- **The leading coefficient.** If the conjecture holds, the bias along the squares is `c(x) N^(1 - 2 eta)` with a coefficient `c(x)` that the spectral measure does not determine. It is a sum over the preimage tree of the best orbit, of `G(a, q)/q` times the partial digit products along the preperiod. Compute it in closed form, and decide whether that of `6:23` vanishes exactly or only nearly.
- **The floor as a theorem.** Prove that `eta > 1/4` makes the structured part of `S(N)` `o(sqrt N)`. That needs the minor arcs, which is the Mauduit-Rivat machinery with explicit constants, and none of the sources read gives them.
- **The open class.** The 16 non-folded odd-base designs each meet all 256 8-blocks below `2^30`, and the 10 of them with `eta > 1/4` read `X2 <= 331.0`, at `5:24`. **That every non-folded odd-base design is normal along `n^2` is a Conjecture**; it lies outside Mullner 2018, Theorem 1.4, and the twin identity ties each pair into one question about `(-1)^n`.
- **Cubes.** Along `n^3` the twisted sums run to `N^3` and the Gauss sums turn cubic, so the same heuristic predicts a structured part `N^(1 - 3 eta)` and a floor at `eta = 1/6`. Thue-Morse, at `0.2075 > 1/6`, would be fast along the cubes, where Spiegelhofer 2023 proves it simply normal.
- **The block length.** Hanna 2016 puts the failure of the Fourier property at blocks of length `A log2 N`. Read the census with growing blocks and find where each design's chi-square leaves the floor.
- **The gap model.** Prove that the lags along the consecutive primes converge to the gap model, and give the rate of normality of Thue-Morse along the consecutive primes as a power of `1/log N`.
