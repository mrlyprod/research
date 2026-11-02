# adjacent-digits

- The numbers of [adjacent](../../../notes/adjacent.md): the set `R_b` of integers with no two equal adjacent base-`b` digits as a `(b+1)`-state automaton, its blocks and units digit, the run expansion of its transform, the chain certificate and its walls, the two-step contraction, a falsification of every inequality behind them, readings of the true rate and the prime count.
- `count`: the automaton against brute force below `b^6` at `b = 3, 4, 5, 10`; the counts `b (b-1)^(L-1)` of strings and `(b-1)^L` of members with `L` digits; the block split against direct membership at `420` values of `x`; the units-digit identity at every `b = 3..12`, `j = 1..6` and top digit `f`; the witness `1, b, b + 1` against a digit design; `(J - I)^n` in closed form; the embedded design of alternating halves at base `b^2`, `b = 4, 6, 8, 10`, inside `R_b` with `kappa_F = 1` and the half interval's transform.
- `identity`: the run expansion and the matrix product against the direct transform at random frequencies.
- `wall`: the constant `h*` as a cell-by-cell upper bound of `sin(pi t)(-psi(t) - psi(1-t))` at 40 digits; the closed form `Q_b(y)` in `mpmath` interval arithmetic at 120 bits; the walls at `1/5` and `1/4`, each certified at every base from the wall to the point where the decreasing tail bound takes over and failing at the base below; `C_1` and `C_e` at the wall; the spread of `alpha_1` beside the Cauchy-Schwarz exponent `1 - log(b - 1)/(2 log b)`; the crossings of `0.22` and `0.18`, scanned on `300..1000` and `2500..5000` with no tail.
- `check`: the full-block and long-run grid sums against `U(N)` and `R_r` on grids of shifts; the grid sums `Sigma_i(s)` by FFT against `C_1 z^i` at the certified `h* = 3.927267` and the blocks against `(C_1 + C_e) z^i`; the two-step rows against `b abs(S_b(v)) + 2b - 1`; `c_E` by base; the C1 bound at rational points with exact phases.
- `rate`: the exact unshifted grid sums by FFT, readings that bound nothing.
- `primes`: `psi_R/A` with prime powers, six decimals, at `10^6`, `10^7`, `10^8` for `b = 10, 7, 3, 4, 5`, `A` counted exactly through the blocks; the Mobius sum over `R_10` below `10^8` by a segmented sieve, checked against `pi(10^8) = 5761455` and the Mertens value `M(10^8) = 1928`.

## THE ROUNDING

- `wall` runs every closed form in interval arithmetic at 120 bits; `h*` is the largest cell bound rounded up at the sixth decimal; upper ends print rounded up, gaps truncated down.
- `count` is exact integer arithmetic.
- `check` runs in float64, the C1 phases exact in rationals, and asserts margins far above float error.
- `rate` and `primes` print readings.

## RUN

- `uv run python research/adjacent/lab/py/adjacent-digits/adjacent.py count` in 3 seconds.
- `uv run python research/adjacent/lab/py/adjacent-digits/adjacent.py identity` in under a second.
- `uv run python research/adjacent/lab/py/adjacent-digits/adjacent.py wall` in 2 seconds.
- `uv run python research/adjacent/lab/py/adjacent-digits/adjacent.py check` in 3 seconds.
- `uv run python research/adjacent/lab/py/adjacent-digits/adjacent.py rate` in under a second.
- `uv run python research/adjacent/lab/py/adjacent-digits/adjacent.py primes` in 8 seconds.
- Prints only, writes nothing; every check raises if it fails.

## WITNESSES

- The automaton: `0` mismatches below `b^6` at `b = 3, 4, 5, 10`, `597870` members at `b = 10`; the block split exact at `420` values of `x`; the units-digit identity at all `450` cells.
- The constant: `h* <= 3.927267`, against `2 gamma + 4 log 2 = 3.927020` at `t = 1/2`.
- The wall at `1/5` is `1153`, certified at every base `1153..2400` with the tail from `2400`, failing at `1152`; `1/5 - alpha_1 >= 1.1311 * 10^-5` at `1153`, `C_1 <= 3.2416`, `C_e <= 2.0484`. The wall at `1/4` is `185`, failing at `184`.
- The band: `alpha_1 < 0.22` at every base `489..1000`, not at `488`; `alpha_1 < 0.18` at every base `3452..5000`, not at `3451`.
- The check: full blocks at most `0.996523` of `U`, long runs at most `0.981280` of `R_r`, grid sums at most `0.304639` of `C_1 z^i`, blocks at most `0.132486` of `(C_1 + C_e) z^i`, two-step rows at most `0.948152` of their bound; `c_E < 1` from `b = 19`, `0.787312` at `1153`.
- The embedded design: `abs(hat F_i) = abs(hat H_(2i))` to `3.4 * 10^-11`.
- The prime count: `psi_R/A` at `10^8` reads `0.999501` at `b = 10`, `0.999148` at `b = 7` and `0.993479` at `b = 3`; the Mobius meter at `b = 10` reads `abs M(10^8)/A^(1/2) = 0.655`.
