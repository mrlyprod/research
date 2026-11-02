# occupancy-decay

- Counts the occupied rays of the ternary gasket at `level` inside a height window, against the first moment of the same window, and measures the digit-congruence density available to any residue-class attack on Conjecture O.
- Conventions, printed by the census itself: the height of a ray is `max(z_1, z_2)` of its primitive direction, the window `octave <= alpha level` is read as the threshold `height <= 3^(alpha level)`, and the octave is `floor(log_3 height)`, one below the generator `lab/rs/dimension-one-ladder/src/census.rs`, which returns `floor(log_3 height) + 1` above `height = 1`; ray totals here exclude the two fibre rays that the Rust census totals include.
- `A(level, alpha)` counts distinct occupied rays in the window and `F(level, alpha) = Sum_z M_level(z)` counts the gasket points carrying them, so `F/A` is the mean multiplier count of a windowed ray and `F` is the first moment Conjecture O is allowed no access to.
- Ratio set: `R_k = {u v^(-1) mod 3^k : (u,v) in G_k, u > 0, 3 does not divide v}` and `sigma_k = |R_k|/3^k`, the exact density of the mod-`3^k` occupancy constraint, with the congruence-collinear pair count `M_2(k) = Sum_r m(r)^2` and its Cauchy-Schwarz floor `sigma_k >= (3^(k-1) - 2^(k-1))^2 / (3^k M_2(k))`.
- Prime sum: `Sum_{p > 3^(beta level)} N_level(p)` at target zero from a smallest-prime-factor sieve over the gcd histogram of `G_level`, against the first-moment bound `(F(level, 3^((1-beta) level)) + 2^(level+1)) / beta`.
- `constants` derives the printed thresholds from the standing window edges: `0.3597878`, `0.5524022`, the congruence caps `c <= 0.2618596` and `alpha <= 0.575328`, and the exponent `theta < 1.8073` that O asks of the ratio set.

## RUN

- `uv run python research/coprime/lab/py/occupancy-decay/occupancy.py check` matches `R_k` against a brute-force image for `k = 2..9`, matches `A` and `F` against direct ray enumeration for `level 4..10`, pins the regression `A(9, 3^7) = 2818`, and confirms the non-fibre count `3^level - 2^(level+1) + 1`, one coordinate of every occupied ray divisible by 3, and every weight coprime to 3.
- `... constants`, `... census 12 13 14 15 16`, `... ratios 15`, `... sieve 12 14` are the kept rows, together about ten seconds.
- `... census 18` is the deepest sweep, 386896202 non-fibre points in 27.5 s at 0.4 GB; `... census 16 --cut 1.0` costs 20.9 s at 3.0 GB and `... ratios 18` 61 s at 3.1 GB, the two memory walls; `--mem` sets the chunk exponent.

## CLAIMS

- `A(level, 0.5533)` reads `822, 1976, 3770, 8000, 16366, 34716, 65342, 139050, 258036` at `level 10..18` with `log_3 A / level` inside `[0.6109, 0.6345]` and `theta = log A / log X` inside `[1.1041, 1.1467]`; at `alpha = 1/2` the same bands are `[0.5416, 0.5798]` and `[1.0833, 1.1596]`, so the earlier reading `0.543 .. 0.557` does not reproduce.
- Occupied non-fibre ray totals `3151656, 9491964, 28545340` at `level 14, 15, 16` and `1044840` at `level 13`, two above each in the fibre-counting convention, a second builder for the census rows.
- `F(level, 0.5533)/A(level, 0.5533)` reads `5.41, 5.20, 5.52, 5.64, 5.63, 5.86, 5.79, 6.08, 5.92` at `level 10..18` while `log_3 F / level` falls `0.7645 -> 0.7201`: the first moment and the ray count share one exponent, the gap `0.0966` at `level 17` being `log(F/A)/(level log 3)`.
- At a fixed height the two part company: `A(level, 3^5) = 384 .. 474` over `level 10..18` while `F(level, 3^5) = 2728 .. 51694`, mean multiplicity `7.10 -> 109.06`, the shift rays carrying it.
- `R_k` is indexed by the modulus `3^k` and `R_1` is empty under the hypothesis `u > 0`, so the sequence starts at `k = 2` and carries that offset wherever it is quoted, `lab/py/ratio-set-saving` included.
- `|R_k| = 1, 3, 9, 23, 63, 168, 457, 1245, 3423, 9447, 26285, 73440, 206149, 580920, 1643545, 4663382, 13272515` at `k = 2..18`, `sigma_18 = 0.034259`, growth `|R_(k+1)|/|R_k|` rising monotonically `2.794 -> 2.8461`, local decay `c_k = 1 - log_3 growth` falling `0.0647562` at `k = 13` to `0.0479313` at `k = 18`, with `k c_k` inside `[0.8418, 0.8628]` over `k = 13..18`: polynomial through the measured range, no exponential floor in sight and none proved.
- `M_2(k)/4^k = 0.3914, 0.4036, 0.4061, 0.4098, 0.4077, 0.4071, 0.4029` at `k = 10..16`, still falling, so `M_2 = O(4^k)` is a measured hypothesis and not a limit; on that hypothesis the Cauchy-Schwarz floor caps every congruence-only decay at `c <= 0.2618596` and `alpha <= 0.575328`, which excludes neither `0.5533` nor `0.5524022`.
- The first-moment inequality holds with room at `level 10, 12, 14` and `beta = 0.45, 0.5, 0.6`, ratio of prime sum to bound between `0.0846` and `0.1517`, worst `0.1517` at `(14, 0.6)`.

## WITNESSES

- the-window-occupancy-and-the-band.md THE WINDOW AT DIMENSION ONE: the occupied-ray exponents `[0.5416, 0.5798]` at `c = 1/2` and `[0.6109, 0.6345]` at `c = 0.5533`, and the non-fibre totals `3151656, 9491964, 28545340` at `level 14, 15, 16`.
- the-window-occupancy-and-the-band.md THE WINDOW AT DIMENSION ONE: Conjecture O trivial below one half, and the first-moment inequality checked at `level 10, 12, 14` with `beta = 0.45, 0.5, 0.6`, worst ratio `0.1517`.
- the-window-occupancy-and-the-band.md THE WINDOW AT DIMENSION ONE: `F/A` at `alpha = 0.5533` reading `5.41 .. 5.92` over `level 10..18`, `log_3 F / level` falling `0.7645 -> 0.7201`, and the fixed-height split `A(level, 3^5) = 384 .. 474` against `F(level, 3^5) = 2728 .. 51694`.
- the-window-occupancy-and-the-band.md THE WINDOW AT DIMENSION ONE: the digit-congruence bound and its seed, `sigma_k` from `0.046063` at `k = 13` to `0.034259` at `k = 18`, `M_2(k)/4^k` over `k = 13..16`, and the caps `c <= 0.2618596`, `alpha <= 0.575328`.
- the-window-occupancy-and-the-band.md THE WINDOW AT DIMENSION ONE: `theta` inside `[1.1041, 1.1467]` at `alpha = 0.5533` and `[1.0833, 1.1596]` at `alpha = 1/2`, against the `1.8073` O asks.
- research/claims/ the `R_k` modulus-indexing and offset row, carried identically by `lab/py/ratio-set-saving`.
