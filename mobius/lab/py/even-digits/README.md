# even-digits

- The numbers of the even digits on [mobius](../../../notes/half-interval-even-digits.md): the digit set `E = {0, 2, .., base-1} = 2F` at the odd base `base = 2 fill - 1`, `F = {0..fill-1}` the half interval, its shifted-grid certificate read through `F`, the dissection re-cut at the denominators `2e`, the arithmetic at `2 base`, the witness at `d = 2`, and a sanity meter.
- `transfer`: at every odd base `3..41` and at `99`, `101`, `103`, `105`, `1001`, `1003`, `94939`, `94941`, every class of `(base mod 4, base mod 3)` printed, on six shifts `s` and every level `i` under the cap, the grid sum `Sigma^E_i(s)`, read from `E`'s own digits at bases `3..41` and through the closed form `|sin(2 pi fill t)/sin(2 pi t)|` above `41`, where it is `|hat F(2t)|` by construction, against `Sigma^F_i(2s)` read from `F`'s, against `(3/2) lambda^i` with `lambda` the half interval's closed form, and against the Parseval floor `base^i`; each asserted; the reading `Sigma^F_i(1/2)/Sigma^F_i(0)` at the top level.
- `regions`: at every grid point of five designs, base 3, 5, 7, 9, 11, the regions A, B, C1', C2' against the exact `sum mu` over the strings of `E`, the counts of C1' and C2', the coprime part of the halved denominator at every C1' point, the `l^1` mass of the `d = 2e` points under `E` and under `F`, and partial summation at those points; the progression identity `M(u; 2e, r)` against the sums `M(u/2^j; e, r 2^(-j))` as integers; Lemma A' through the halving at seeded draws off `d = 2e`; the witness `a = (y-1)/2` at `d = 2`, its limit over the levels, and the level where it beats the pair's bound.
- `meter`: below `10^7` at six bases `5..10007`, the identity `2 M_E(2X) = M^tw_F(X) - M_F(X)` at every `X`, the readings `max |M_E|/A_E^(1/2)` and `max |M^tw_F|/A_F^(1/2)`, the prime count `psi_E(10^7)` against `kappa_E A_E(10^7)`, and `kappa_E = base/(base + 1)` as a fraction.

## THE ROUNDING

- Every check runs in float64, except the progression identity and the meter identity, which are exact integers, and `kappa_E`, an exact fraction.
- Upper bounds print rounded up and lower bounds truncated down at the last printed digit; masses, meters and the reading `Sigma^F_i(1/2)/Sigma^F_i(0)` are readings with no direction.

## RUN

- `uv run python research/mobius/lab/py/even-digits/even.py transfer` in 30 seconds.
- `uv run python research/mobius/lab/py/even-digits/even.py regions` in 3 seconds.
- `uv run python research/mobius/lab/py/even-digits/even.py meter` in 5 seconds.
- Prints only, writes nothing; every check raises if it fails.

## WITNESSES

- The certificate transfers: `Sigma^E_i(s) = Sigma^F_i(2s)` to `2.2 * 10^-12`, at most `0.775356 (3/2) lambda^i`, at least `1.154700 base^i`, all six classes of `(base mod 4, base mod 3)` sampled.
- The dissection on the even digits: the regions sum to `-37`, `105`, `-37`, `21`, `-3` at bases 3, 5, 7, 9, 11; the `d = 2e` points carry `3.374999 fill^k` to `5.799374 fill^k` under `E` against at most `0.001722 fill^k` under `F`; partial summation at most `0.310113` of its bound.
- The arithmetic at `2 base`: `1950` integer identities.
- Lemma A' through the halving: `6000` draws, log margin at least `11.9560`; the witness limit `0.466274` at base 3, `0.527585` at 5, `0.636613` at `94939`, beating the pair's bound from level `35`, `117`, `281` at bases 3, 5, 7.
- Readings: `max |M_E|/A_E^(1/2)` `0.900` to `1.296`, the same at `1009` and `10007` because there the maximum sits below the base, where `S_E` is every even number; `psi_E(10^7)` from `4 log 2` to `18 log 2` against `kappa_E A_E(10^7)` from `49206.7` to `2501749.0`.
