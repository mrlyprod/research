# gaussian-digits

- The numbers of the Gaussian design on [coprime](../../../notes/gaussian-primes-on-a-design.md), section PRIMES ON A DESIGN, subsection Gaussian primes on a design: base `b` on both coordinates of `Z[i]`, `F = {0..b-1}^2` less one vector `a`, its two-dimensional `l^1` exponent `alpha_2` in units of `b^(2k)` against the gates `1/5`, `1/6` and `1/8`, and a Gaussian prime count against the main term.
- `chain`: the square chain `z_c - 1 = sum_l (lambda_l^+)^2 z_c^(-l)`, `lambda_l^+ = (2/pi) l log b + gamma' + (2/pi) b^(-l)`, in closed form, against `z_c < b^(2 delta)(1 - b^(-2))`; located in floats, certified in `mpmath.iv` at 120 bits on interval blocks from the wall to the base where the cap `max(3, 1 + sqrt(3 c_pi^2 L_b^2 + 3 c_pi L_b gamma' + gamma'^2 + tau_b))` sits below the target with the slope test positive, a failing base below each wall; then the box chain against the Morton chain at base `b^2` at `b = 10^1..10^6`.
- `uwin`: the square uniform window, one Collatz-Wielandt certificate per base for every missing vector, from the cell suprema `(Delta[k_1] Delta[k_2] + 1)/(b^2 - 1)` of the majorant `(|D_b(t_1)| |D_b(t_2)| + 1)/(b^2 - 1)`; scanned down from the chain wall to the first failing base at `1/5` and `1/6`, and read at five bases from `300` to `3000`.
- `census`: every orbit of `m = a - ((b-1)/2, (b-1)/2)` under sign changes and swap at every base `3..44`, 4046 orbits, each printed with the one-digit shift floor and maximum of `Sigma_1` and the two-digit window matrix read with sampled cell suprema and infima; the closed form `|K(t_1) K(t_2) - e(<m, t>)|/(b^2 - 1)` is checked against the digit sum at every vector of bases 3, 4, 5 first. Pass `lo hi` for a smaller range.
- `floor`: `min_s Sigma_1` near the grid at every orbit of `3..31`, a `13 x 13` scan of `|s_i| <= 2.5/b^3` refined by a shrinking pattern search, after the witness base 7 missing `(3, 3)` at `s = (0.019396/7, 0.980604/7)`; it prints every orbit that dips below `2` and tests that the dips are exactly the orbits with a zero coordinate of `m` at odd bases. Pass `lo hi` for another range.
- `morton`: the same two-digit window reading on the box and on the Morton image at base `b^2`, at `b = 8, 12, 16, 24, 32` and four vectors each.
- `minor`: the minor-arc lemma read on the box `[1, 2048)^2`, `X = 8380418`, which fills its norm range: `|S(xi)|` over all prime-power elements, prime powers included, against `X^(7/8) + X N(q)^(-1/2) + X^(11/16) N(q)^(1/2) + X^(1/2) N(q)` with constant 1 and no logarithm, at 1000 frequencies near sampled `lambda/q` in ten bands of `N(q)` from `X^0` to `X^0.45`, half at random within `4/N(q)` and half within `1/X`; it asserts that the box mass exceeds the least value of the bound, so the reading can fail, and that every ratio stays below 1.
- `gate`: the sum `S(xi)` on the boxes `[1, y)^2`, `y = 3^6, 3^7, 3^8`, all prime-power elements, against `y^(8/5)`, at 40 frequencies near sampled `lambda/q` in each of five bands of `N(q)` near `y^0.8..y^1.2`; it asserts that the box mass exceeds ten times `y^(8/5)`, so the reading can fail, and that in every band the largest ratio falls from `729` to `6561`.
- `hcount`: the count lemma by brute force at `H = 10^6`, `#{h : N(h) <= H, ||h xi|| <= eta}` over the whole disc against `(H/N(q) + 1)(N(q) eta^2 + 1)`, at 184 frequencies within `4/N(q)` of sampled coprime `lambda/q` and eight `eta` each, and at the cell `q = 9 + 4i`, `lambda = 53 + 77i`, `xi = lambda/q`, `H = 1940`, `eta^2 = 5/97`; it asserts the ratio stays below the `2^18` of the proof.
- `bilinear`: the bilinear bound, `B` over the annuli `N(v) in [V_0, 2 V_0)`, `N(u) in [U_0, 2 U_0)` with `alpha` the extremal choice for Cauchy-Schwarz and `beta` random unimodular or `Lambda`, against `norm2(alpha) norm2(beta) (U_0 + V_0 + N(q) + U_0 V_0/N(q))^(1/2)`, at five shapes and ten bands of `N(q)`, then the operator norm, the largest singular value of `e(Re(v u xi))` on the discs `N(v), N(u) < 2 U_0` of the statement, at `U_0 = 200`, `q = 1`, `xi = 0` and at `U_0 = 100`, `q = 1 + i`, `xi = 1/(1 + i)`; it asserts every ratio stays below the `2 * 10^4` of the proof.
- `sep`: the Fourier `l^1` norm in `(log abs(z), arg z)` of the square `[1, 2]^2` smoothed at `h_s = 2^(-4)..2^(-9)` by a product of `C^infinity` steps, against a polar box smoothed at the same scale, with a grid check at `2^(-7)` and `2^(-8)`; it asserts the square's growth per halving stays above `1.3` and the polar box's below `1.2` at the last step.
- `count`: `sum_(z in S_k) Lambda(z)` over a sieve of the norms against `(4/pi) kappa_F |S_k|`, `kappa_F = prod_(pi | b) N(pi)/(N(pi) - 1)` times the share of `F` coprime to `b`, at ten designs of bases 3, 4, 5, 7.

## THE ROUNDING

- `chain`: margins, caps and slope tests are `mpmath.iv` intervals at 120 bits and pass only on the sign of the unsafe endpoint; `gamma' = 0.9625229` is asserted above `(2/pi)(gamma + log(8/pi))` first; on a block `[lo, hi]` the root equation is read with `log hi` in the terms increasing in `b` and `lo` in the terms decreasing in `b`, against the target at `lo`.
- `uwin`: sines are evaluated after reducing the argument mod `2 pi` in integers, then widened by `10^-12`, against a measured deviation of `1.2 * 10^-15` from 80-bit values at 800 sampled points; the Collatz-Wielandt ratio is inflated by `10^-9`; printed bounds round up at six decimals, failing readings round down.
- `kern` reduces its argument to `(-1/2, 1/2]` and carries the sign `(-1)^((b-1) n)` of the shift by `n`, so the kernel and the phase `e(<m, t>)` stay consistent at negative shifts and at `t = 1`.
- `minor` is a float reading of a bound that carries an unstated constant and `(log X)^4`; read with constant 1 and no logarithm it can fail, and it proves nothing.
- `gate`, `hcount` and `bilinear` are float readings: `gate` of the sum itself, `hcount` an exact integer count at float frequencies, `bilinear` of the extremal form; `sep` is a discrete Fourier sum on a grid fine enough that halving its step moves the reading by `0.0018` relative, a reading of a norm and not a bound on it.
- `census`, `floor` and `morton` are float readings with sampled cell suprema or a local search, never bounds; `count` is a float sum of exact sieve data, a reading of the main term's shape.

## RUN

- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py chain` in under a second.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py uwin` in 13 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py census` in 84 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py floor` in 10 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py morton` in under a second.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py count` in 2 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py minor` in 8 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py gate` in 18 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py hcount` in 23 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py bilinear` in 40 seconds.
- `uv run python research/coprime/lab/py/gaussian-digits/gauss.py sep` in 9 seconds.
- Every check raises if it fails.

## WITNESSES

- The square chain: walls `72`, `267`, `5054` against `1/5`, `1/6`, `1/8`, certified on `[72, 153]`, `[267, 806]`, `[5054, 30781]`; `alpha_2 < 0.199772` at `72`, `0.166614` at `267`, `0.124999` at `5054`.
- The box against the Morton chain: exponent ratios `1.1863`, `1.2696`, `1.3247`, `1.3648`, `1.3959`, `1.4211` at `b = 10^1..10^6`.
- The square uniform window: `alpha_2 < 1/5` at every base `43..71`, largest `0.199242` at `43`, `0.200340` at `42`; `alpha_2 < 1/6` at every base `107..266`, largest `0.166554` at `107`, `0.166814` at `106`; `0.125280` at `1000`, `0.117274` at `2000`.
- The census: below `1/5` first at `31`, at every orbit from `33`; orbit spread at most `0.0137` from `base 6` and `0.0034` from `base 21`.
- The floor: `Sigma_1 = 1.997526` at the base 7 witness; 135 orbits dip, exactly those with a zero coordinate of `m` at the odd bases `3..31`, least `1.948655` at base 3 missing `(1, 1)`; 1359 orbits stay at `2`, among them every orbit of every even base `4..30`.
- Morton: the box reading exceeds the Morton reading by `0.046..0.106` at all 20 rows.
- The count: ratios `0.9991..1.0002` at all ten designs.
- The minor-arc reading: 370005 prime-power elements, 369593 of them prime; box mass `1.8884` of the bound's least value; largest `|S|/bound` `0.7456` in the bands up to `X^0.10`, `0.1170` at `X^0.25`, `0.0153` at `X^0.35`.
- The gate read on full blocks: masses `17.74`, `27.58`, `42.81` times `y^(8/5)`; largest `abs(S)/y^(8/5)` in the bands `y^0.8..y^1.2` `0.2648, 0.1326, 0.0927, 0.1105, 0.1192` at `729`, `0.1437, 0.0674, 0.0601, 0.0903, 0.0616` at `2187`, `0.1000, 0.0493, 0.0290, 0.0294, 0.0432` at `6561`.
- The count lemma: largest sampled ratio `18.6834` at `N(q) = 194625`, `eta = 10^(-3)`, over `3141549` points and `N(q)` from `2` to `441650`; the cell `q = 9 + 4i` reads `10.8651`, count `1369`; readings, not the supremum.
- The bilinear bound: largest ratio on the annuli `1.4444` at `U_0 = 300`, `V_0 = 6000`, `N(q) = 1105`, below the operator norm; the operator norm reads `6.1940` times the bracket at `U_0 = 200`, `q = 1`, `sigma_1 = 1245`, and `6.1640` at `U_0 = 100`, `q = 1 + i`.
- The separation cost: square `6.711, 9.706, 13.735, 19.242, 26.791, 37.321` at `h_s = 2^(-4)..2^(-9)`, growth `1.3923..1.4462` per halving; polar box `5.747..14.209`, growth `1.2590` falling to `1.1573`.
