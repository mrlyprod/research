# kummer-carry

- The numbers of [kummer](../../../notes/kummer.md): the set of `k` with `p` not dividing `C(3k, k)` as a two-state carry automaton, the mass of every Kummer set `{k : p does not divide C(ak, bk)}`, the matrix transfer certificate in the row-0 basis, its walls, a falsification of every inequality behind it, readings of the true rate and a sanity meter.
- `kummer`: the four digit intervals, the row sums `(p+1)/2` and the pair `{0, 1}` at every prime `5..2999`; the witness `u, v p, u + v p` against a digit design at every one of them; the automaton against Kummer's carry count `(s_p(k) + s_p(2k) - s_p(3k))/(p - 1)` at every `k < p^5` for `p = 5, 7, 11, 13`, at every `k < 101^3` directly, and at every `k < 101^5` through `k = 101^2 q + r`, the digit sums splitting exactly over the carry classes of `r`.
- `family`: the count `((p+1)/2)^L` read off Legendre's digit sums for every `2 <= a <= 7`, `1 <= b < a`, prime `a < p <= 37`, and where it fails at `p <= a`; the row count and the pair over every state of every `(a, b)` with `a <= 12` and prime `a < p < 400`; the row deviations at `p = 10007`; the product-set test below `p^3` at every `2 <= a <= 6`, `1 <= b < a`, prime `a < p <= 31`; the row-0 limit designs and their constant against `3 - 4 gcd(a, b)/a`.
- `wall`: the constants `nu^+`, `H^+`, `nu_1^+`, `H_1^+` as grid maxima plus a Lipschitz slack, then certified: the Lipschitz sums `S` and `S'` against their closed forms, `nu_1^+` against `(sqrt 6 + sqrt 2)/3`, and `nu^+`, `H^+`, `H_1^+` by a dyadic interval cover of `[0, 1/2]`; the closed form `lambda/fill` in `mpmath` interval arithmetic at 120 bits; the tail constant `c_inf`; the walls at `1/5` and `1/4`, each certified at every prime from the wall to the point where the tail takes over and failing at the prime below; the spread by class of `p mod 12`; the crossings of `0.22` and `0.18`.
- `check`: `G*`, `G_01`, `C`, `D` summed directly at every prime `13..401` and six larger primes on `201` shifts against their pointwise bounds and closed forms, and `T psi <= lambda psi`; the product formula and the row-0 basis against the direct transform at level `3`; the shifted grid sums by FFT against `lambda^i`.
- `rate`: power iteration of the row-0 operator on a grid of cells with linear interpolation, the exact unshifted grid sums of `K` and of the row-0 design `F*` by FFT, and both offsets carried to the bar `1/5`; readings that bound nothing.
- `meter`: the Mertens meter and the prime count on `K` at `p = 101`, `1009`, `10007` below `10^7`, a sanity print far below every wall.
- `small`: the prime count with prime powers on `K` at `p = 5, 7, 11, 13` up to `10^8`, `A_K(10^8)` counted exactly through the blocks by state and checked against direct membership at `10^6`, a falsification print.

## THE ROUNDING

- `wall` takes each constant as its grid maximum plus half a grid step times a Lipschitz constant plus `10^-12`, rounded up at the sixth decimal, the Lipschitz constant the upper end of its interval value; then it proves each bound by halving `[0, 1/2]` into dyadic cells until every cell has its center value, in interval arithmetic at 120 bits from the moduli themselves, plus half-width times the Lipschitz constant under the printed bound. Every closed form runs in interval arithmetic at 120 bits, lower bounds print their lower end truncated down, upper bounds their upper end rounded up.
- `kummer` and `family` are exact integer computations.
- `check` runs in float64 and asserts strict inequalities whose margins sit far above float error.
- `rate`, `meter` and `small` print readings.

## RUN

- `uv run python research/kummer/lab/py/kummer-carry/kummer.py kummer` in under a second.
- `uv run python research/kummer/lab/py/kummer-carry/kummer.py family` in 5 seconds.
- `uv run python research/kummer/lab/py/kummer-carry/kummer.py wall` in 15 seconds.
- `uv run python research/kummer/lab/py/kummer-carry/kummer.py check` in 10 seconds.
- `uv run python research/kummer/lab/py/kummer-carry/kummer.py rate` in 22 seconds.
- `uv run python research/kummer/lab/py/kummer-carry/kummer.py meter` in 2 seconds.
- `uv run python research/kummer/lab/py/kummer-carry/kummer.py small` in 3 seconds.
- Prints only, writes nothing; every check raises if it fails.

## WITNESSES

- The automaton: `0` mismatches at every `k < p^5` for `p = 5, 7, 11, 13, 101`; `abs(K below 101^5) = 345025251 = 51^5`.
- The mass: `0` mismatches of `((p+1)/2)^L` at `p > a`; `0` of `73654` rows off `(p+1)/2`; the pair in every state at `p >= 2a - 1`; `(3, 1)` at `p = 3` holds `1` integer below `3^11`.
- The constants: `nu^+ = 1.854292`, `H^+ = 5.672774`, `nu_1^+ = 1.287903`, `H_1^+ = 3.861981`; Lipschitz constants `76.8912161.. <= 76.892` for `H` and `19.890173.. <= 19.891` for `H_1`; the cover closes in `2529`, `207` and `472` cells for `nu`, `H`, `H_1`; `2 nu^+/pi <= 1.180479`, `c_inf <= 7.482879`.
- The wall at `1/5` is `13567831`, certified at all `34446` primes `13567831..14133194`, least gap `>= 4.8500 * 10^-6` at `13567859`, failing at `13567811` with gap `<= -9.8373 * 10^-6`; `1/5 - alpha_1 >= 7.5864 * 10^-5` at the wall and `>= 1.1060 * 10^-8` at `13567859`, the least over the certified primes. The wall at `1/4` is `227189`, failing at `227177`.
- The row-0 constant `3 - 4 gcd(a, b)/a` at every `a <= 7`; row deviations at most `2a - 4` digits at `p = 10007`; over `127` cells only `(2, 1)` is a product set below `p^3`.
- The small primes: `psi_K/((p/(p+1)) A_K)` at `10^8` reads `1.0032`, `1.0041`, `1.0010`, `1.0051` at `p = 5, 7, 11, 13`.
- The check: `G*/bound <= 0.765814`, `G_01/bound <= 0.897928`, `T psi/(lambda psi) <= 0.905536`, grid sums at most `0.480589` of `lambda^i`.
- Readings: the row-0 operator `(10/(3 pi)) log p + 3.25..3.39`; the exact grid sums of `K` `(10/(3 pi)) log p + 1.40..1.49`.
