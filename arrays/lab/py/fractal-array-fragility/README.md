# fractal-array-fragility

- Builds the fractal sparse array `F_r = {sum_(i<r) g_i b^i : g_i in G}` of a generator `G` with hole-free coarray `[-a, a]`, at the base `M = 2a + 1` of Cohen and Eldar and at every compressed base `a < b <= 2a + 1`.
- Tests essentialness by the definition: remove the sensor, recompute the difference coarray, compare. A fast test by pair counts (a lag of weight 1 kills both ends, a lag of weight 2 kills the sensor its two pairs share) is checked equal to the definition on every array up to 250 sensors before it is used alone; `dial` builds the weight level by level as the convolution of `w_G` stretched by `b^(r-1)`, checked against the counted weight up to 250 sensors.
- `exact`: every hole-free generator with `L <= 6` and span `a <= 13`, one of each mirror pair, 119 generators, at `r = 2, 3`; checks the essential set of `F_r` against the digit rule of Theorem A (every digit in `U(G)`), checks the tightness criterion `e_r = card E(G)^r` iff `E(G) = U(G)`, asserts the coarray hole-free, and lists the loose generators with `card E(G)`, `u` and whether `G` is maximally economic.
- `dial`: the generators `{0, 1, 2}`, `{0, 1, 4, 6}`, `{0, 1, 2, 3, 7}`, `{0, 1, 4, 7, 9}`, `{0, 1, 2, 3, 7, 11}` at every base `a < b <= 2a + 1`, `e_r` for every `r` with at most 17000 sensors and 8 million lags, the hole-free coarray `[-A, A]`, `A = a(b^r - 1)/(b - 1)`, asserted on every row, and the affine law `e_(r+1) = lambda e_r + c` the row obeys from `r = 2` or `3`, or `u^r`; `lambda >= 1` and `c` are fitted on two steps, taking `lambda = 1` on a flat tail, no upper bound, and the line prints how many further steps check the law.

## RUN

- `uv run python research/arrays/lab/py/fractal-array-fragility/fragility.py exact` from the repo root; about 21 seconds.
- `uv run python research/arrays/lab/py/fractal-array-fragility/fragility.py dial`; about 21 seconds.
- Standard library and numpy; no argument runs both verbs.

## WITNESSES

- [arrays](../../../notes/arrays.md), Fragility is exact - Theorem A against the definition on 238 cases with 0 mismatches, and `{0, 1, 2}` at 4 of 9 and 8 of 27: verb `exact`.
- [arrays](../../../notes/arrays.md), The two bounds - the tightness criterion with 0 failures (a consequence of the Theorem A check), `{0, 1, 2, 3}` at 4 of 16, the 21 loose generators, 14 of them maximally economic, and the excess `card E(G) - u = 1`: verb `exact`.
- [arrays](../../../notes/arrays.md), The compressed base - the hole-free coarray, the table of `e_r`, the affine fits with their checks, and 6561 sensors with 87381 lags for `{0, 1, 2}` at base 4: verb `dial`.
