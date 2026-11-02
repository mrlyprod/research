# base3-exponent

- The `l^1` exponent `alpha_1` of the digit transform of `F = {0,1}` at base 3 on [the unconditional dissection](../../../notes/unconditional-dissection.md), enclosed in interval arithmetic: here `|hat F(t)| = |cos(pi t)|`, `hat F_level(t) = prod_(j<level) hat F(3^j t)`, `alpha_1 = 1 + lim (1/level) log_3 int_0^1 |hat F_level|`, which is also the limit of `(1/N) log_3 max_x Sigma_N(x)` for `Sigma_N(x) = sum_(i<3^N) prod_(j<N) |hat F(3^j (x + i/3^N))|`.
- It runs the library half of [digit-transform-norms](../../../../coprime/lab/py/digit-transform-norms/README.md), everything above its verb dispatch: the closed form, the interval kernel, the window tables and the Collatz-Wielandt passes; nothing there is edited.
- The one-step number `log 2/log 3 = 0.630930` is the exponent of the one-digit sum `2 = max_x Sigma_1`, an upper bound that is not the value.

## WHAT IT COMPUTES

- `check`: the closed form against the character sum, the exact one-digit sum `1 + 1/2 + 1/2 = 2`, the entropy floor, the shift sandwich at `N = 1, 2`, one window row, and the order floor < lower < upper < two-step bound.
- `enclosure`: all of that at full strength, the shift sandwich at `N = 1, 2, 4, 6`, the two-step norm `||T^2 1|| = fill^2 max_x Sigma_2`, the window ladder at `6, 8, 9, 10, 11, 12` digits, the enclosure that takes the best side of each, and the float grid sums `Sigma_L(0)` for `L = 1..12` as a reading.
- The floor. Jensen's formula gives `int_0^1 log|cos(pi t)| dt = -log 2`, so `int |hat F_level| >= 2^(-level)`; the mean over the shift `x` of `Sigma_N(x)` is `3^N int |hat F_N|`, so `max_x Sigma_N >= (3/2)^N` and `alpha_1 >= 1 - log 2/log 3 = log(3/2)/log 3`. It is exact arithmetic, checked against 30-digit quadrature.
- The shift sandwich. `min_x Sigma_N` and `max_x Sigma_N` over one period `1/3^N`, a grid of `m` shifts plus the Lipschitz slack, give `(1/N) log_3 min_x Sigma_N <= alpha_1 <= (1/N) log_3 max_x Sigma_N`: the sandwich of Chow, Varju and Yu 2024, Theorem 4.2, both sides run in interval arithmetic. Both sides move slowly in `N`.
- The window bracket. `G(w)` bounds `|hat F|` above and below on the cell behind a window `w` of `n_d` digits, the transfer matrix on `3^(n_d - 1)` states carries a window to its three successors, and the Perron root of the supremum matrix bounds `alpha_1` above and that of the infimum matrix below, the lower test on the states the Perron vector keeps. This is Maynard's Markov eigenvalue bound in Karwatowski's `lambda(s, J)` form, and per digit of work it closes faster than the sandwich does.
- The masked lower side. `|hat F|` vanishes at `t = 1/2`, so cells near it have infimum `0` and the infimum matrix loses states. The tree's restriction clause takes the test vector from the Perron vector of the whole infimum matrix and its lower end falls with the digits (`0.6045801` at 9, `0.4220878` at 10, nothing at 12). Here the test vector is the Perron vector of the restricted matrix itself: power iteration with the dropped states held at zero, the dropped set iterated to a fixed point, then the same certified Collatz-Wielandt pass on the kept states. The test vector is only a test vector, and the entries are nonnegative, so words inside the kept states are a subset of the path sum that bounds `int |hat F_level|` below, and the Collatz-Wielandt pass on the kept block gives `M_K y >= mu y` with `y > 0`, so the lower end stays a bound whatever the dropped states are.

## THE ROUNDING

- Every enclosure comes from the tree's interval kernel with outward `nextafter` on every operation; lower bounds truncate down and upper bounds round up at the seventh digit, each checked as `mu > 3^lower` or `mu < 3^upper` before it prints.
- The floor is exact arithmetic, its integral read at 30 digits as a check.
- The grid sums `Sigma_L(0)` are float64 readings, bound nothing, and are printed beside the enclosure.

## RUN

- `uv run python research/mobius/lab/py/base3-exponent/exponent.py check` in 3 seconds.
- `uv run python research/mobius/lab/py/base3-exponent/exponent.py enclosure` in 80 seconds, 50 of them the two-step norm.
- Prints only, writes nothing; every asserted row raises if it fails.

## WITNESSES

- The enclosure `alpha_1(3, {0,1})` in `[0.6055183, 0.6055241]`, width `5.8 * 10^(-6)`; the upper end from the supremum matrix at 12 digits and sub-scan 16 (`mu_hi = 1.9449498`), the lower end from the masked infimum matrix at 12 digits on `167416` of `177147` states after 5 rounds; equivalently `hat kappa_1 = 1 - alpha_1` in `[0.3944759, 0.3944817]`, below `1/2` as Remark 6.1 of Chow, Varju and Yu 2024 says at base 3.
- The ladder, lower and upper at `6, 8, 9, 10, 11, 12` digits: `0.6036172, 0.6053102, 0.6054510, 0.6054978, 0.6055130, 0.6055183` and `0.6074103, 0.6057322, 0.6055918, 0.6055449, 0.6055296, 0.6055241`; the lower end rises and the upper falls at every step.
- The two-step bound: `max_x Sigma_2 < 3.879398` at `m = 500000`, so `||T^2 1|| = 4 max_x Sigma_2 < 15.517591` (the exact `4 Sigma_2(0)` reads `15.51754`) and `alpha_1 < 0.6169966`; the one-step `0.6309298` is exact because `Sigma_1(x) = 2 cos(pi x)` on `[-1/6, 1/6]`, the period `1/3` and even.
- The shift sandwich rows: `(0.5584958, 0.6169966)` at `N = 2`, `(0.5814728, 0.6115566)` at `N = 4` with `m = 4000`, `(0.5884817, 0.6104209)` at `N = 6` with `m = 1500`.
- The floor `alpha_1 > 0.3690702`, above region B's `1/4` and region A's `1/5` at every shift.
- The reading: `log_3 (Sigma_L(0)/Sigma_(L-1)(0))` is `0.605521` at every `L = 6..12`, inside the enclosure.
