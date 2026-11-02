# interval-digits

- The numbers of the half interval on [mobius](../../../notes/half-interval-even-digits.md): the digit set `F = {0..fill-1}` at the odd base `base = 2 fill - 1`, its transfer certificate, the walls it sets at the bars `1/5` and `1/4`, a falsification of every inequality behind it, readings of the true rate, and a sanity meter.
- `wall`: the constant `c_inf` of the tail bound `lambda/fill <= (2/pi) log base + c_inf + 3.4/base`, proved from base `101` and checked directly at every odd base `101..3001` and at `94939`, `200001`, `10^6 + 1`, `10^8 + 1`, the least base where the tail bound clears each bar, the closed form `lambda/fill` certified at 120 bits at every odd base from the wall to that point and failing at the base below, `bar - alpha_1` and the gap at the wall, the bars read as the exact intervals `1/5` and `1/4`, the lower bound `min G/fill >= (2/pi) log base - 1.31` at every odd base `9..9999` with its tail from `101`, the density floor `1 - alpha_base < 1/5`, `alpha_1` at three larger bases, and the base where the constant of the Section 9 sketch of Maynard 2022 clears `1/5` on the same set, calibrated against its one-missing-digit crossing `1520573`.
- `check`: the three sums of the certificate's proof, the `bP`, `bQ` and signed `aP` sums, against their bounds at every odd base `3..401` and at `1001`, `10001`, `100001` on `201` shifts; then at every odd base `3..401` on `801` shifts and at five larger bases on `401` shifts, the closed form of `S(t)` against the direct sum, the bound on `G(t)`, the inequality `T phi <= lambda phi` and the lower bound on `G(t)`, each asserted; the readings `G(0)/fill - (2/pi) log base` against `1 + gamma'` and `max G/fill - (2 sqrt2/pi) log base`; the level sums at 40 digits at nine small bases, four shifts and every level up to `base^level <= 9261` against `(3/2) lambda^level`; and the transfer identity at level 2.
- `rate`: power iteration of the transfer operator on a grid of cells with linear interpolation, readings of the true rate near the wall, and the one-step constant `G(1/2)/fill` against both bars.
- `meter`: the Mertens meter and the prime count on the half interval at the prime bases `101`, `1009` and `10007` below `10^7`, a sanity print far below every wall and never evidence for the theorem.

## THE ROUNDING

- `wall` evaluates every closed form in `mpmath` interval arithmetic at 120 bits; lower bounds print their interval's lower end truncated down, upper bounds the upper end rounded up.
- `check` prints every upper bound rounded up and every lower bound truncated down at its last printed digit; it runs in float64 and asserts strict inequalities whose margins sit far above float error except at `t = 0`, where the bound on `G` is tight by design to `1.5 * 10^-4` of its size; the level sums run at 40 digits.
- `rate` prints readings: the interpolated power iteration and the one-step scan bound nothing.

## RUN

- `uv run python research/mobius/lab/py/interval-digits/interval.py wall` in under a second.
- `uv run python research/mobius/lab/py/interval-digits/interval.py check` in 13 seconds.
- `uv run python research/mobius/lab/py/interval-digits/interval.py rate` in about 20 seconds.
- `uv run python research/mobius/lab/py/interval-digits/interval.py meter` in 2 seconds.
- Prints only, writes nothing; every check raises if it fails.

## WITNESSES

- The tail constant `c_inf <= 2.60043004`, the direct check's smallest gap `>= 1.3332 * 10^-7`; the wall at the bar `1/5` is `94939`, the closed form certified at every odd base `94939..94947` and the tail bound from `94946`, failing at `94937` with gap `<= -2.6733 * 10^-5`, `1/5 - alpha_1 >= 1.3678 * 10^-8` at the wall.
- The lower bound holds at every odd base from `9`, its tail margin `>= 0.016518` from `101`, and the factor is negative at `7`; the density floor `1 - alpha_base < 1/5` holds at every odd base from `27` and fails at `25`.
- The wall at the bar `1/4` is `3789`, certified at `3789..3793`, failing at `3787`, `1/4 - alpha_1 >= 7.9625 * 10^-6` at the wall.
- `alpha_1 <= 0.1993872` at `100003`, `0.1761232` at `1000003`, `0.1331636` at `10^9 + 7`.
- The Section 9 constant of Maynard 2022 at the consecutive half interval clears `1/5` from `7777884825`.
- The check: the `bP` sum under its bound by at least `3.99 * 10^-5`, largest `G/bound` `0.999857`, largest `T phi/(lambda phi)` `0.999866`, smallest `G/lower` `1.534638`, largest level sum over `(3/2) lambda^level` `0.557408`, transfer identity to `10^-15`.
- Readings: `G(0)/fill - (2/pi) log base` is `1.962430` at `100001` against `1 + gamma' = 1.962523`; the true rate reads `(2/pi) log base + 2.26` near `7 * 10^4` and crosses `base^(1/5)` between `70001` and `80001`; the one-step constant reads `(2 sqrt2/pi) log base + 1.19` and clears `1/4` from `7075` and `1/5` from `317063`.
