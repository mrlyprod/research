# mrly-euler

- The multiplicative structure of a digit design `S_F` over the digit set `F`: which designs carry an Euler product over primes, and what stands in its place when none does.
- `wall`: for every digit set of every base `2 <= base <= 12`, whether the indicator of `S_F` is multiplicative, with the constructed coprime witness verified and an independent minimal witness searched.
- `pair`: the Dirichlet coefficients of `zeta_F(s) M_F(s)`, and the least `n > 1` where the product leaves `1`.
- `position`: the position-product identity `sum_(n in D_level) n^(-s) = int_0^1 G_level(t) Z(s,t) dt` with `G_level(t) = prod_(i<level) sum_(d in F) e(d base^i t)` and `Z` the periodic zeta, evaluated against the direct sum; the Lerch evaluator is checked against mpmath `polylog`.
- `dual`: the periodic zeta rebuilt from Hurwitz's formula, DLMF 25.13.3, solved for `F(-t, s)`.
- `fibre`: the Lerch-Mobius series `M(s, a/Q) = sum_(n>=1) mu(n) e(-n a/Q) n^(-s)` written as a finite combination of inverse Dirichlet L-functions of modulus dividing `Q`, checked as a coefficient identity.
- `word`: the free-monoid zeta `1/(1 - fill base^(-s))` expanded as an Euler product over Lyndon words in the digit alphabet.
- `beurling`: the free semigroup `N_B` on the primes lying in `S_F`, with `pi_F(x)`, `N_B(x)`, the Mertens `M_B(x)` and its exponent by residue class of `log_base x`, to `x = 10^6`.
- Domain: `base <= 12` for the wall, minimal witnesses searched to product `4000`; coefficients to `n = 4000` for the pair; `level <= 5` for the position identity; `n = 3000` for the fibre identity; `u^16` for the word product; `x = 10^6` for the Beurling census.

## RUN

- `uv run python research/zeta/lab/py/mrly-euler/euler.py wall`
- `uv run python research/zeta/lab/py/mrly-euler/euler.py pair`
- `uv run python research/zeta/lab/py/mrly-euler/euler.py position`
- `uv run python research/zeta/lab/py/mrly-euler/euler.py dual`
- `uv run python research/zeta/lab/py/mrly-euler/euler.py fibre`
- `uv run python research/zeta/lab/py/mrly-euler/euler.py word`
- `uv run python research/zeta/lab/py/mrly-euler/euler.py beurling`
- `position` runs a few minutes; every other verb is seconds.

## WITNESSES

- `8177` digit sets at `2 <= base <= 12`: `4083` with `1` outside `F`, `11` full, `4083` witnessed, split `4072` repunit, `5` odd base, `6` even base; hardest minimal witness `base 12`, `F = {1}`, pair `(5, 377)`, product `1885`
- least witnesses at small base: base 3 `{0,1}` `(2, 5)`, base 3 `{1,2}` `(2, 5)`, base 4 `{1,2,3}` `(2, 9)`, base 5 `{1,2,3,4}` `(2, 13)`, base 2 `{1}` `(3, 5)`
- `zeta_F M_F` leaves `1` at `n = 4` for base 3 `{0,1}`, `n = 10` for base 3 `{1,2}`, `n = 9` for base 10 missing `9`, `n = 10` for base 10 missing `0`; over 257 sets the least such `n` is at most `50`, and the eight full sets have none below `4000`
- position identity error `1.95e-16` and `2.04e-16` at base 10 missing `9`, `level 3` and `level 4`, `s = 3.3`; `1.67e-16` and `1.75e-16` at `s = 2.7 + 1.9i`; `2.9e-16` at base 3 `{0,1}`, `level 3, 4, 5`; Lerch evaluator against `polylog` below `1.7e-25`
- Hurwitz reflection error below `2.1e-30` at `s = 3.3`, `2.7 + 1.9i` and `0.6 + 4.1i`
- fibre identity maximum coefficient error `2.6e-12` over eleven pairs `(Q, a)` including `Q = 3, 9, 27, 100`, Euler-factor step below `7.4e-16`
- Lyndon expansion equals `1/(1 - fill u)` through `u^16` at `fill 2, 3, 4, 9, 10`; `c_2(level) = 2, 1, 2, 3, 6, 9, 18, 30, 56, 99`, [A001037](https://oeis.org/A001037)
- Beurling census: base 3 `{0,2}` has the single prime `2` and `M_B(x) = 0` for `x >= 2`; base 3 `{0,1}` has `525` primes, `N_B(920483) = 2198`, running `max abs(M_B) = 98`, exponent `0.3339` against `alpha/2 = 0.3155`; base 10 missing `9` has `35139` primes, `N_B(10^6) = 488864` against `x^alpha = 531441`, `M_B(10^6) = 1860`, running max `1866`, exponent `0.5452` against `alpha/2 = 0.4771`, climbing `0.4203, 0.4882, 0.5452` at `10^4, 10^5, 10^6`
- the full base-10 control reproduces `M(10^4) = -23`, `M(10^5) = -48`, `M(10^6) = 212`, [A084237](https://oeis.org/A084237)

## SOURCES

- [DLMF 25.13](https://dlmf.nist.gov/25.13) - the periodic zeta `F(x,s)` at 25.13.1 and Hurwitz's formula at 25.13.3.
- [DLMF 25.12](https://dlmf.nist.gov/25.12) - the polylogarithm expansion at 25.12.12, the evaluator used here.
- [Beurling 1937](https://doi.org/10.1007/BF02546666) - generalised prime systems.
- [Diamond, Montgomery and Vorhauer 2006](https://doi.org/10.1007/s00208-005-0638-2) - a Beurling system with a regular integer count whose zeta has infinitely many zeros off any fixed half plane.
