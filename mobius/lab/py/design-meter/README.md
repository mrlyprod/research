# design-meter

- Reads the oscillation spectrum of the digit-design Mobius meter `M_F(x) = sum of mu(n)` over `n in S_F`, `n <= x`, with `F` the kept digits, normalised by `x^(alpha/2)` with `alpha = log_base fill`, and asks which frequencies it carries.
- The full set is the control: `M(x)/sqrt(x)` to `10^8` reproduces the nontrivial zeta zeros, so any frequency a design does or does not carry is measured against a pipeline that provably sees zeros.
- Three candidate frequency families are tested per design: the zeta ordinates `gamma`, the design's pole lattice `2 pi j / log base`, the ordinates of the primitive quadratic Dirichlet `L`-function of conductor `base`, plus the `alpha`-rescaled zeta ordinates and two rigid displacements of the zeta list as null controls.
- Every series is resampled uniformly in `log x` on 32768 points by exact step evaluation, Hann-windowed, mean-removed, and read as `gamma = 2 pi f` from the real FFT power spectrum; the noise floor is a 101-bin running median, so a score is power over local floor.
- A family is scored by the mean of `log10` of its targets' scores against a null of 4000 rigid circular shifts of the same target list inside the band `4 < gamma < 60`; the null therefore keeps the family's own spacing statistics and only moves it.
- Two controls run on every design at the same depth: the counting function `A_F(x)/x^alpha` on the same elements, and a sign vector supported exactly where `mu` is nonzero. The first says whether a lattice would be seen; the second says what a meter with no arithmetic looks like.
- The split test cuts the log range in half and correlates the two disjoint halves' spectra, calibrated against 24 support-matched random-sign draws; the family verb reports, for each design it compares, the peak count reached by eight support-matched random-sign meters, which is the null the peak counts are read against.
- The mean-field split defines `R_F` by `M_F(x) = sum of mu(n) A_F(n)/n + R_F(x)` and reads the spectrum of all three terms, which separates a design's own frequencies from the classical Mertens function seen through the design's density; it also prints the echo's share of the meter at three depths against the rate `x^(-min(alpha, 1 - alpha)/2)`, and prints two designs below `alpha = 1/2` where the echo converges instead of decaying.

## RUN

- `uv run python research/mobius/lab/py/design-meter/design_meter.py sieve`
- `uv run python research/mobius/lab/py/design-meter/design_meter.py spectrum`
- `uv run python research/mobius/lab/py/design-meter/design_meter.py family`
- About fifty seconds, three minutes and two minutes; prints only, writes nothing. Peak memory is the `10^8` linear Mobius sieve.
- `sieve` enumerates each design and factors every element exactly (trial division to the cube root with a shrinking survivor set, then a perfect-square test and deterministic Miller-Rabin on the twelve witnesses `2..37`), and asserts the meters against an independent generator. `spectrum` is the full read. `family` is the equal-`alpha` comparison across bases.

## READS

- Control, base 10 full set to `10^8`: `M(10^8) = 1928`, which is [A084237](https://oeis.org/A084237); bin width `0.3544`; all ten of the strongest peaks fall within one bin of a nontrivial zeta zero, offsets `0.018` to `0.196`.
- Base 10 with the digit `9` missing, `A_F = 43046721` elements, `alpha = 0.954243`: all six of the strongest peaks fall within one bin of a zeta zero, offsets `0.068` to `0.216`.
- The same design split by the density echo: the echo carries six of six top peaks at zeta zeros, the residual `R_F` carries none of the two it has, and the echo is `0.1342` of the meter in root mean square against `0.6476`.
- The echo's share decays with depth: `0.356028, 0.242495, 0.207229` at base 10 missing `9` and `0.208549, 0.099001, 0.047902` at base 3 `{0,1}`, against the rate `x^(-min(alpha, 1 - alpha)/2)`; the printed share over prediction reads `1.0000, 0.7387, 0.6846` and `1.0000, 0.9219, 0.8663`, so each design decays at least as fast as the rate. The generator prints ratios and never an exponent, and the share is measured against the meter, so reading it as a rate against `x^(alpha/2)` assumes the square-root conjecture.
- Base 16 `{0,1}` to `10^7`, `alpha = 0.25`, 63 elements: the echo `sum mu(n) A_F(n)/n` settles near `-0.11`, reading `-0.0937, -0.1330, -0.1242, -0.1051, -0.1099` at `10^3` to `10^7` against `x^(alpha - 1/2)` reading `0.1778, 0.1000, 0.0562, 0.0316, 0.0178`, the ratio climbing `0.53, 1.33, 2.21, 3.33, 6.18`. Base 10 `{0,1}`, `alpha = 0.301030`, 128 elements: `-0.0500` at `10^7` against `0.0405`, ratio `1.24`.
- Below `alpha = 1/2` the tail `sum_(m in S_F) H(m-1)` converges absolutely, so the echo tends to a nonzero constant and `O(x^(alpha - 1/2 + eps))` is false there; the surviving rate against `x^(alpha/2)` is `x^(-alpha/2)`.
- Base 3 `F = {0,1}` at level 14, 16, 18, 20, 22: the zeta family scores `0.528, 0.621, 0.693, 0.742, 0.759` against nulls `0.167` to `0.302`, every rung at `p <= 0.008`; the same list displaced by `+2` scores `0.187` and reflected scores `0.445`.
- The pole lattice on the meter of base 3 `{0,1}`, base 3 `{0,2}` and base 5 `{0,1}` scores `-0.592`, `-0.640`, `-0.734` against nulls `0.30`, `0.27`, `0.23`, below its own null; on the counting function over the identical elements it scores `3.602`, `3.764`, `3.973` against nulls `0.72`, `0.95`, `1.15`.
- Base 3 `{0,1}` against base 9 `{0,1,2,3}`, equal `alpha = 0.630930` and equal element count `1048575`: ten peaks against none, where support-matched random-sign meters on the same two supports reach `0` to `4` and `0` to `2` peaks over eight draws each. The spectra correlation printed beside it, `-0.047`, runs against no null here and carries nothing; the peak counts are the read.
- Base 9 `{0,1,3,4}` is base 3 `{0,1}`, its digits being the base-3 pairs `00, 01, 10, 11`, so the generator asserts the elements and the `mu` values equal and the printed correlation `1.0000` is an assertion of that identity and not a measurement.

## CHECKS

- Meters asserted against the exact census of [mobius-designs](../../rs/mobius-designs/): base 3 `{0,1}` reads `(11, 105)`, `(149, 173)`, `(-30, 312)`, `(496, 539)`, `(1009, 1089)` at level 14, 16, 18, 20, 22, base 3 `{0,2}` reads `(-382, 485)` at level 20, base 3 `{1,2}` reads `(-1461, 1582)` at level 18, base 10 missing `9` reads `(2181, 5234)` at `10^8`, each as `(M_F, max abs M_F)`.
- The two `mu` paths agree: exact factorisation of a design's elements against the linear sieve.
- The `L`-function zero finder is checked by its own output, the first ordinates reading `8.0397` at conductor 3, `6.0209` at conductor 4 and `6.6485` at conductor 5, each with real part `0.5` to the printed precision.

## NOTE

- A design has far fewer samples than the full set at the same depth, so resolution comes from the log range and not the element count; the bin widths here are `0.2579` to `0.3871` against the control's `0.3544`.
- The family score against its shift null is deterministic to the printed digits; the `p` beside it is Monte Carlo over 4000 shifts and moves by a factor of two in the tail between runs, so a claim rests on the score and its null, never on the `p` alone.
- A score of 3 over the local median is not a detection: for exponential noise it fires on one bin in eight, and the printed base rate per design says so. The threshold here is 8, and no claim rests on a single peak.
- The lattice band starts at `gamma = 4`, below the first lattice line of every base read, so a lattice is inside the band it is tested in.
