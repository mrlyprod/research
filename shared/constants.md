# Constants

Real numbers fall out of the fractal work beside the integer sequences: sieve limits, dimensions, pressure zeros, alignment reaches. This page is the ledger they are cited from. A row holds the digits its generator prints, truncated and never rounded, so every digit shown is a digit of the constant; a bracket row shows the digits its two ends share.

- `name` is the row's anchor, `#name` below, for a note, a claim or an OEIS link to point at.
- `generator` is what prints the digits: a crate function, a lab study and its verb, or a figure census.
- `digits` counts the significant digits shown; `OEIS` is the entry whose digits match the row at its full printed window in the OEIS stripped dump, or none.

| name | value | digits | generator | topic note | OEIS |
| --- | --- | --- | --- | --- | --- |
| [wallis-sieve-cube](#wallis-sieve-cube) | `0.948815485719679` | 15 | `mrly_num::sieve::solid_limit` | [pi](../pi/notes/pi.md) | none |
| [wallis-sieve-cube-even](#wallis-sieve-cube-even) | `1.052420668` | 10 | `mrly_num::sieve::parity_product(3, false)`, `magic/lab/py/staircase-constants` section parity | [pi](../pi/notes/pi.md) | none |
| [wallis-sieve-plane](#wallis-sieve-plane) | `0.785398163` | 9 | `mrly_num::sieve::PLANE_LIMIT`, `magic/lab/py/staircase-constants` section parity | [pi](../pi/notes/pi.md) | [A003881](https://oeis.org/A003881) |
| [wallis-sieve-dim4](#wallis-sieve-dim4) | `0.985352084` | 9 | `mrly_num::sieve::parity_product(4, true)`, `magic/lab/py/staircase-constants` section parity | [pi](../pi/notes/pi.md) | [A371219](https://oeis.org/A371219) |
| [wallis-sieve-dim6](#wallis-sieve-dim6) | `0.998553027` | 9 | `mrly_num::sieve::parity_product(6, true)`, `magic/lab/py/staircase-constants` section parity | [pi](../pi/notes/pi.md) | [A371465](https://oeis.org/A371465) |
| [dim-e2](#dim-e2) | `0.53128050627720514162446864736847178549305910` | 44 | `beneath/lab/py/periodic-determinant` verb control | [beneath](../beneath/notes/beneath.md) | [A279903](https://oeis.org/A279903) |
| [log2-phi](#log2-phi) | `0.694241913630` | 12 | `beneath/lab/py/question-mark` verb table, code 126 | [beneath](../beneath/notes/beneath.md) | [A242208](https://oeis.org/A242208) |
| [orphan-11892](#orphan-11892) | `0.41681776443307102575640207215` | 29 | `beneath/lab/py/periodic-determinant` verb orphans | [beneath](../beneath/notes/beneath.md) | none |
| [carpet-reach-7-3](#carpet-reach-7-3) | `0.79100941515` | 11 | `spirograph/lab/rs/roulette-reaches`, track `7/3` inside, marks `(-6, -9)` | [spirograph](../spirograph/notes/spirograph.md) | none |
| [slice-ladder-rate](#slice-ladder-rate) | `1.346122006764217277640827069917397630838` | 40 | `cuts/lab/py/slice-ladder-rate` verb product | [cuts](../cuts/notes/cuts.md) | none |

## ROWS

### wallis-sieve-cube

- The volume of the solid Wallis sieve: cut every surviving cube into `(2k+1)^3` at level `k` and drop the centre one, so `prod_(n odd >= 3) (1 - n^-3) = pi^(3/2) / (8 |Gamma(7/4 - i sqrt3/4)|^2)`, the odd parity design of `dim 3` at code `150` on the row word.
- `solid_limit` reads it from the log series with an Euler-Maclaurin tail and lands within four ulps of `0.9488154857196796`; the lab's Gamma form at roots of unity meets the zeta log series `exp(-sum_m ((1 - 2^(-3m)) zeta(3m) - 1)/m)` below `1e-50`, and at `150` digits in mpmath the Gamma form meets that series to `1e-150`, continuing `0.94881548571967959932844868318037787316013952414766658520425938671492997254356892986`.
- Claims: [wallis-sieve](../pi/claims/wallis-sieve.md), [staircase-law](../magic/claims/staircase-law.md). The exponents `2`, `4` and `6` are in the OEIS; this odd member is not.

### wallis-sieve-cube-even

- The even parity design of `dim 3`, code `105`: `prod_(n odd >= 3) (1 + n^-3) = pi^(3/2) / (8 |Gamma(5/4 - i sqrt3/4)|^2)`, and the cube constant times this is the exponent-`6` product A371465.
- The lab prints `1.052420668` against the zeta log series below `1e-50`; the Gamma form continues `1.05242066816787920549`.
- Claim: [staircase-law](../magic/claims/staircase-law.md).

### wallis-sieve-plane

- The plane Wallis sieve keeps the area `pi/4`, Wallis's product; the odd parity design of `dim 2`, code `6`.
- The lab's Gamma form meets `pi/4` and the zeta log series below `1e-50`; `mrly_num::sieve::ratio` on the odd word of two million letters reads `0.785398262`, the truncated product.
- Claims: [wallis-sieve](../pi/claims/wallis-sieve.md), [staircase-law](../magic/claims/staircase-law.md).

### wallis-sieve-dim4

- `prod_(n odd >= 3) (1 - n^-4) = pi cosh(pi/2) / 8`, the odd parity design of `dim 4`, code `27030`.
- The lab's Gamma form meets the closed form and the zeta log series below `1e-50`.
- Claim: [staircase-law](../magic/claims/staircase-law.md).

### wallis-sieve-dim6

- `prod_(n odd >= 3) (1 - n^-6) = pi cosh(pi sqrt3/2) / 24`, the odd parity design of `dim 6`, code `7608434000728254870`.
- The lab's Gamma form meets the closed form and the zeta log series below `1e-50`.
- Claim: [staircase-law](../magic/claims/staircase-law.md).

### dim-e2

- The Hausdorff dimension of `E_2`, the irrationals whose continued fraction quotients are `1` and `2`: the pressure zero of the transfer operator of `A = {1,2}`, which is the memory dial's run-length code `126` and the gasket's `E_2` exponent `delta_2`.
- The periodic-point determinant at periods to `18` in interval arithmetic brackets it in `[0.5312805062772051416244686473684717854930591089, 0.5312805062772051416244686473684717854930591092]`, containing the interval of Jenkinson and Pollicott; the `40`-mode collocation of `question-mark` verb table reads `0.5312805062772048`, gap `3.3e-16`.
- Claims: [question-mark](../beneath/claims/question-mark.md), [apollonian-gasket](../apollonian/claims/apollonian-gasket.md).

### log2-phi

- `log_2 phi`, the growth rate of the golden rule, code `7` at `k = 2`, whose transfer matrix has `det(I - x T) = 1 - x - x^2`, so its first pole comb sits on `Re s = log_2 phi`; also `log_2 rho` of the run-length code `126`, `A = {1,2}`.
- The table verb prints `rho = 1.618033988749` and `log_2 rho = 0.694241913630` on the code `126` row; the closed form continues `0.6942419136306173`.
- Claims: [memory-dial](../beneath/claims/memory-dial.md), [question-mark](../beneath/claims/question-mark.md).

### orphan-11892

- The pressure zero `dim_CF` of the width-`4` orphan code `11892`, the run words over `{1,2}` without the pair `22`, at the supergolden `rho`.
- The periodic-point determinant on its six-state recurrent graph brackets it in `[0.416817764433071025756402072152, 0.416817764433071025756402072154]`; the collocation reading `0.416817764433` is the bracket truncated at twelve digits.
- Claim: [question-mark](../beneath/claims/question-mark.md).

### carpet-reach-7-3

- The reach at which four of the carpet wheel's trochoids on the circle track `7/3`, inside, run through one point: the transversal alignment with marks `(-6, -9)`, corner seat modulus `0.527339610104`, ring `1.005704332357` wheel radii.
- The study bisects it to an f64 bracket of width `1.1e-16` and prints `0.791009415157`, the exact law residual `4.18e-15`; the twelve printed decimals round safe, so eleven are truncated safe.
- Claim: [spirograph-reaches](../spirograph/claims/spirograph-reaches.md).

### slice-ladder-rate

- `r_inf = 1/prod_(k>=2) cos(2 pi/3^k)`, the per-dimension decay rate on the eigenvalue scale of the slice excess `rho_dim - fill/3` of the base-3 central diagonal slice; the product itself is `0.7428747134175312701067637570301888840694`, and `2 r_inf = 2.69224401352843455` is the rate on the dimension scale.
- The verb multiplies 139 cosines at 60 and 120 digits, the two agreeing on all 40 truncated digits, the tail below `pi^2/(4 9^140) = 6.28e-134`; the verb `rate` reads the same constant off exact Perron roots, the one-step ratio deflated by the form's `dim`-dependence within `1.06e-17` of it at `dim = 120`, and retires the fit `1.34612251727283689`, seven true digits.
- Claims: [diagonal-slice-ladder](../cuts/claims/diagonal-slice-ladder.md), [slice-ladder-rate-and-dead-routes](../cuts/claims/slice-ladder-rate-and-dead-routes.md).
