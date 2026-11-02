# digit-sum-twist

- The numbers of [gelfond](../../../notes/gelfond.md): the two-position constant and the exponent it buys, a falsification of the lemma on the supremum of the twisted transform, a check of the twisted peel, and the twisted Mobius and von Mangoldt meters and the digit-sum classes of the primes of two base-10 designs below `10^7`.
- `lemma`: the constant `c = 1 - min(G_1, G_2, G_3)/fill` in its pair, kernel and variance forms at base 10 missing `7` over sixteen classes of `alpha`, with `gamma`; then `1/5 - alpha_1`, `c`, `gamma`, `sigma` and `delta_0` at `base 92317` (step 3 constant), `39363` (chord constant, middle digit) and `10^6`, at `beta = 1/2`, `0.2` and `10^-3`, the variance form read at the extreme excluded digit; the base where `sigma` overtakes `1/5 - alpha_1`, by bisection between `92317` and `10^6` at `beta = 1/2` and `0.2`; the sharpness example `F = {0, 2}` at base 5, `alpha = 1/8`, read at `k = 4, 8, 12`; and the step 3 margin at `92316`; then the step 3 wall in `mpmath` interval arithmetic at 200 bits, `1/5 - alpha_1` enclosed at `92317` and `92316`, lower end floored and upper end rounded up at nine digits.
- `sup`: `sup_t abs(hat w_k(t))/fill^k` at `k = 1..8` at base 10 missing `7` and `k = 1..4` at base 37 missing `18`, bracketed by a max-product pass over windows of `W` base-`b` digits of `t` (6 at base 10, 4 at base 37), each window factor bounded by its centre value plus the Lipschitz margin `pi h sum(F)`, and by the value attained at the best path refined on a local grid; each off-lattice class is set against `(1 - c)^(floor(k/2))`, each lattice class against the exact value `1` at `t = -alpha`.
- `peel`: the shifted-grid supremum `B_10(F)` at missing `7` on `2 * 10^5` shifts with its Lipschitz margin, then `960` checks of `Sigma_i(s) <= B^i`, of the one-step bound and of the one-step identity at depths `1..5`, and the unshifted twisted mass `Sigma_5(0)`.
- `meter`: `abs(M_F(x; alpha))/A_F(x)` at `10^4..10^7`, its running maximum, and `psi_F(x; alpha)/psi_F(x)` at `10^5..10^7` against the Ramanujan value on the lattice and against the average of `e(alpha s_b(n))` over the design's elements coprime to `90` off it, for base 10 missing `7` and missing `0`.
- `lattice`: the primes of base 10 missing `7` and missing `0` above `30` and below `10^6` and `10^8`, counted by `p mod 9`, `s_10(p) mod 2, 3, 6, 18` and `(p mod 9, s_10(p) mod 6)`, `(p mod 7, s_10(p) mod 3)`, each class against the law of the section Progressions and every modulus in units of the square root of its predicted count, and against the design's `8`-digit strings prime to `30030`, read exactly by a digit recursion on `(n mod 90090, s mod 2)`; the even share of those strings prime to `30`, `210`, `330`, `390`, `2310` and `30030`; `abs(g_F(1/(2(b + 1))))/fill`; and the bias from `11`, the even share prime to `330` less prime to `30`, at `6`, `7`, `8`, `10` and `12` digits, with its ratio over two digits at even counts.
- `gelfond`: the primes of the same two designs below `10^7` counted by `s_10(p) mod q` at `q = 2, 4, 5, 7, 8, 10, 3, 9`, against `1/q` and against the elements of the design coprime to `90`.

## THE ROUNDING

- The large-base constants are computed in `mpmath` at 40 digits and floored at four significant digits, so every printed `c`, `gamma`, `sigma`, `delta_0` and `1/5 - alpha_1` is a lower bound.
- The base-10 constants are floats printed at six digits and quoted on the page truncated at three or four.
- The supremum's upper bracket is a float computation with an explicit Lipschitz margin; it is a reading, never a certificate; the generator raises if the attained bracket exceeds the proved bound, a refutation, or if the upper bracket does at `k >= 2`, an inconclusive reading.

## RUN

- `uv run python research/gelfond/lab/py/digit-sum-twist/twist.py lemma` in under a second.
- `uv run python research/gelfond/lab/py/digit-sum-twist/twist.py sup` in 10 seconds.
- `uv run python research/gelfond/lab/py/digit-sum-twist/twist.py peel` in 10 seconds.
- `uv run python research/gelfond/lab/py/digit-sum-twist/twist.py meter gelfond` in 7 seconds, peak resident memory 1.7 GB.
- `uv run python research/gelfond/lab/py/digit-sum-twist/twist.py lattice` in 1 second, peak resident memory 0.6 GB.
- Prints only, writes nothing; every check raises if it fails.

## WITNESSES

- The constants: `c >= 0.1500`, `gamma >= 0.0352` at `alpha = 1/2` at base 10 missing `7`; `delta_0 >= 4.678 * 10^-9` at `base 92317`, `beta = 1/2`, region A binding through `1/5 - alpha_1 >= 1.871 * 10^-8`; `3.47 * 10^-10` at `beta = 10^-3`; `6.741 * 10^-8` at `39363`; `2.819 * 10^-4` at `10^6`; `sigma` binds from the switch at `100084` at `beta = 1/2` and `93409` at `beta = 0.2`; the sharpness example reads `1` at `k = 4, 8, 12`; the step 3 wall is interval-certified, `1/5 - alpha_1 >= 1.87123612 * 10^-8` at `92317` and `<= -4.82092732 * 10^-8` at `92316`.
- The supremum: at `k = 8`, `alpha = 1/2`, attained `0.071152`, upper `0.071165`, proved `0.521973`; upper below proved at all eleven off-lattice classes and every `k >= 2`; value `1` at `t = -alpha` at the five lattice classes, `1/3` among them.
- The peel: `B_10(F) <= 24.702781`; `Sigma_i(s)/B^i <= 0.999914` over `960` checks; identity to `1.40 * 10^-11`.
- The meters: `abs(M)/A <= 0.0005` and running maximum `<= 0.0006` at `10^7`; the von Mangoldt ratio within `0.0013` of the Ramanujan value on the lattice and within `0.0076` of the model off it.
- The lattice: `2079504` and `3103338` primes below `10^8`; lattice classes within `0.83` units of the law; the even share of `s_10(p)` `0.48942` and `0.49464` below `10^6` and `10^8` missing `7`, `0.48418` and `0.49144` missing `0`; the even share below `10^7` `0.49588` and `0.49576`; strings even share `0.50000` prime to `30`, `0.49583` prime to `330` and `2310`, `0.49996` prime to `390`; the bias from `11` at `6, 7, 8, 10, 12` digits `-0.00769, -0.00248, -0.00417, -0.00210, -0.00103` missing `7`, ratios `0.54, 0.50, 0.49` per two digits against `0.480`; primes within `2.37` units of the strings prime to `30030`; `abs(g_F(1/22))/fill = 0.6927` and `0.7491`.
- The classes: `266823` and `397866` primes; at `q = 10` deviation `0.24272` against the model's `0.23939` on missing `7`, primes within `0.02213` and `0.02584` of the model at every `q` coprime to `9`.
