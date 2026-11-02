# prime-dissection

- The numbers of the dissection for the von Mangoldt function on [coprime](../../../notes/primes-on-a-design.md), section PRIMES ON A DESIGN: the wall where a shifted-grid `l^1` certificate below `1/5` holds at one excluded digit, the singular series as the principal characters of region C2, the four regions for `Lambda` on the whole grid of three small designs, and a prime count against the main term.
- It imports `pb_step3`, `strings`, `convergent` and `smooth` from [mobius-dissection](../../../../mobius/lab/py/mobius-dissection/) and `uniform_alpha` from [digit-uniform-bound](../digit-uniform-bound/), read only.
- `wall`: the digit-uniform chain `(z-1)^3 = (2/pi)(log base) z + gamma'(z-1) + (2/pi)(z-1)^2/(base z - 1)` against `z < base^(1/5)(1 - 1/base)`, certified in `mpmath.iv` at 120 bits from the least base where it clears to the base where the closed-form cap `1 + sqrt(2 (2/pi) log base + 0.97)` takes over, with the slope test that keeps the cap below from `base 100` on; the least base where Maynard's written constant `log((q/(q-1)) log q + 3q/(q-1))/log q` drops below `1/5`, and its value at `2000001`; the missing-digit budget at `base 10^7` and `10^8` for Maynard's `C_(q,s) = 1 + (2+s)/log q`, for the wall condition (W) through `pb_step3`, and for the chain; the chain's float wall at two and three excluded digits.
- `walls`: one wall per number `m = 1..13` of excluded digits. The chain at `m` digits, `(z - m)(z - 1)^2 = m((2/pi)(log base) z + gamma'(z - 1) + (2/pi)(z - 1)^2/(base z - 1))` against `z < base^(1/5)(1 - m/base)`, is located in floats, certified in `mpmath.iv` at 120 bits from the wall to the base where the cap `m + sqrt(m((2/pi) log base + 0.97))` takes over (point by point near the wall, on interval blocks split until they clear above it), with the slope test `base^(1/5) sqrt((2/pi) log base + 0.97) > (5/2) sqrt(m)(2/pi)` there and a failing base below the wall; the (W) wall `sqrt(m) + Phi_base/base < base^(1/5)(1 - m/base)` is certified at its points and above by a smooth upper bound on `Phi_base/base` whose gap grows from `base 327`; the margins of both certificates at and one below every wall, the thinnest over all 26; then the better certificate per `m`, and the `m >= 14` argument through the smooth (W) gap at `base m^5`.
- `window`: the digit-uniform window at two window digits against `1/5`, one outward-rounded Collatz-Wielandt certificate per base covering every excluded digit, scanned down from `583` to the first base that fails.
- `series`: `sum_(d | base) mu(d)/phi(d) sum_f c_d(f)` over `F`, divided by `fill`, against `kappa_F = (base/phi(base)) #{f in F : gcd(f, base) = 1}/fill` in exact rationals at every set missing one or two digits of every base `3..30`, and the constant printed in Maynard 2022 Theorem 1.3 of arXiv v1 beside it at six sets; the identity is checked against its own closed form, so the prime count of `count` is the test of the main term.
- `regions`: every grid point `a mod y` of base 10 missing 5 and missing 1 at level 6 and base 5 missing 2 at level 9, cut into regions A, B, C1, C2 as on [mobius](../../../../mobius/notes/mobius.md), with `Lambda` in place of `mu`; the region sums against the exact sum of `Lambda` over the strings, the principal characters `y mu(d)/phi(d)` at the C2 points with `h = 0` against `kappa_F fill^k`, and each region's share.
- `count`: `sum_(n <= x, n in S_F) Lambda(n)` over a sieve to `10^8` against `kappa_F A_F(x)`, `A_F` by digit counting checked against enumeration to `10^5`, at the powers of the base and twelve seeded `x`, on eight sets, two of them without two consecutive digits.

## THE ROUNDING

- `1/5 - alpha_1` at the wall is read from the float root's upper endpoint and floored at three digits.
- The chain margin, the cap gap and the slope test are `mpmath.iv` intervals at 120 bits, and a row passes only on the sign of the unsafe endpoint; `gamma' = 0.9625229` is asserted above the true `(2/pi)(gamma + log(8/pi))` first.
- The printed `alpha_1` of the chain is `log_base(z base/(base - 1))` rounded up at six decimals from the float root; the window rows are the rounded-up exponents of [digit-uniform-bound](../digit-uniform-bound/).
- `series` is exact; `regions` and `count` are floats and their ratios are readings of the main term's shape at bases far below any wall, never a bound.

## RUN

- `uv run python research/coprime/lab/py/prime-dissection/primes.py wall` in 2 seconds.
- `uv run python research/coprime/lab/py/prime-dissection/primes.py walls` in one second.
- `uv run python research/coprime/lab/py/prime-dissection/primes.py window` in 21 seconds.
- `uv run python research/coprime/lab/py/prime-dissection/primes.py series` in under a second.
- `uv run python research/coprime/lab/py/prime-dissection/primes.py regions` in one second.
- `uv run python research/coprime/lab/py/prime-dissection/primes.py count` in 3 seconds the first time, writing `count-100000000.json` under `data/`, then read from that cache.
- Every check raises if it fails.

## WITNESSES

- The chain: certified on `[584, 1272]`, tightest margin `6.0170 * 10^-3` at `584` in the units of the root equation cleared of denominators, `1/5 - alpha_1 >= 1.79 * 10^-5` there, margin `-8.3138 * 10^-3` at `583`; cap gap `4.6127 * 10^-4` at `1272`, slope test `3.3175` at `100`; `alpha_1 < 0.199983` at `584`.
- The walls: chain `584`, `4692`, `17596`, `46798`, `102133`, `195891` at `m = 1..6`, certified up to `1272`, `5477`, `27761`, `87303`, `212356`, `439640`; (W) `92317`, `124332`, `153738`, `182147`, `210153`, `238044` at `m = 1..6` and `265981`, `294064`, `322357`, `350906`, `379742`, `408886`, `438358` at `m = 7..13`; the chain better at `m <= 6`, (W) from `m = 7`; smooth (W) gap at least `0.2947` at `base 14^5`; over all 26 walls the thinnest pass `1/5 - alpha_1 >= 8.42 * 10^-10`, (W) at `m = 12`, and the closest failure one below `-9.835 * 10^-8`, the chain at `m = 11`; `w - z = -1.083 * 10^-7` at `4691`, `m = 2`, and `w - PB = -9.992 * 10^-8` at `265980`, `m = 7`.
- The window: `alpha_1 < 1/5` at every base `301..583`, largest `0.199923` at `301`, `0.200021` at `300`.
- Maynard's constant below `1/5` from `1520573`, `0.197309` at `2000001`; budgets `s <= 7`, `(W) m <= 176`, chain `m <= 16` at `10^7` and `19`, `703`, `29` at `10^8`; the chain's float walls `4692` and `17596` at two and three excluded digits.
- The series: equal at all `4956` sets; `5/4` at base 10 missing `{0, 5}` against the printed `10/9`.
- The regions: the four sums meet the exact total at all three designs; the principal characters return `1.000000 kappa_F fill^k`; C2 on the grid `0.999546`, `0.999021`, `1.000427`.
- The count: `0.9994..1.0001` at the largest power of the base and `0.9984..1.0009` at the seeded `x` on the six sets with two consecutive digits; below `0.0006` at base 3 missing `1` and base 5 missing `{1, 3}`.
