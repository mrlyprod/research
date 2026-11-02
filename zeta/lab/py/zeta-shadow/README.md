# zeta-shadow

- The MECHANISM of the ordinate shadow: why the zeros of the design zeta `zeta_F(s) = sum_(n in S_F) n^(-s)`, with `F` the digit set, sit near the zeros of `zeta`, and how far.
- The shadow is a first-order perturbation, so the offset from a zeta zero `rho_0` is a computed number and not a trend: it is `abs(zeta_F(rho_0))/abs(zeta_F'(rho_0))` with no constant, and `abs(zeta_F(rho_0))/(c abs(zeta'(rho_0)))` under a split `zeta_F = c zeta + E_F`.
- The evaluator, the truncation bounds and the pole data come from `../design-zeta`; the second family and its census come from `../zeta-family` and `../zeta-locus`. This study adds the fibre weight, the per-zero prediction and two new rungs.

## THE SPLIT

- The position identity of `../mrly-euler` reads `sum_(n in D_level, n >= 1) n^(-s) = int_0^1 G_level(t) Z(s,t) dt` with `G_level(t) = prod_(i<level) sum_(d in F) e(d base^i t)` the position product, `Z(s,t)` the periodic zeta and `Z(s,0) = zeta(s)`; `D_level` is the set of strings of `level` digits over `F`.
- Its discrete form is exact and finite: `1_(D_level)(n) = base^(-level) sum_(a mod base^level) G_level(a/base^level) e(-n a/base^level)` for `0 <= n < base^level`, so `zeta_(F,level)(s) = base^(-level) sum_(a mod base^level) G_level(a/base^level) S_level(s, a/base^level)` with `S_level(s,x) = sum_(1 <= n < base^level) e(-nx) n^(-s)`. Verb `mass` reproduces `1_(D_level)` from the transform to `8.326e-40` at level 2 on every design of the ladder.
- THE PRINCIPAL FIBRE. The term `a = 0` is `base^(-level) G_level(0) S_level(s,0) = (fill/base)^level` times the partial sum of `zeta` to `base^level`, since `G_level(0) = fill^level`. That weight is exact and needs no arc, and the verb prints it against `(fill/base)^level` at every rung.
- THE LEVEL IS NOT FORCED. The identity splits the POLYNOMIAL at level `level` against a TRUNCATED `zeta`, while the object of this study is the continued `zeta_F` against the full `zeta`; `(fill/base)^level` falls to `0` with `level`, and to the right of the abscissa both series tend to `1`, so the fibre gives `fill/base` only at level 1 and that reading is a choice.
- The continuous form is the mass of `G_level` on the principal arc `abs(t) < 1/(2 base^level)`, an exact sinc sum `1/base^level + sum_(n in D_level, n > 0) sin(pi n/base^level)/(pi n)`; it equals `kappa_level(F) (fill/base)^level` with `kappa_level` running `0.6015221` to `0.96774464` over the ladder at level 1, 2, 3, a factor of `1.6`, and the full set's own readings falling `0.81830989, 0.70926038, 0.65068647` toward `Si(pi)/pi = 0.5894898722`. It is the same weight up to that shape factor and adds no constant the fibre does not already give.

## THE FIRST-ORDER SHADOW

- At a zero `rho_0 = 1/2 + i gamma` of `zeta` one has `zeta(rho_0) = 0`, hence for ANY constant `c` the split `zeta_F = c zeta + E_F` gives `E_F(rho_0) = zeta_F(rho_0)` exactly. That identity carries no information about `c`.
- A zero of `zeta_F` near `rho_0` solves `c zeta(s) = -E_F(s)`, so to first order `s = rho_0 - zeta_F(rho_0)/(c zeta'(rho_0))`.
- STEP, THE PRIMARY COLUMN. Reading `c zeta'(rho_0)` as `zeta_F'(rho_0)` removes the constant and gives `s = rho_0 - zeta_F(rho_0)/zeta_F'(rho_0)`, Newton's own first step, which is Taylor's theorem at a simple zero of `zeta_F` and needs no split at all.
- PRED, THE SECONDARY COLUMN. The same law at the level 1 fibre reading `c = fill/base`. The lab prints both and lets the sweep choose.
- Both columns are computed from `zeta_F`, `zeta` and the digit density alone and see no design zero. The offset is one complex number, so at the zeros the law pairs the ordinate offset and the real-part offset are one quantity with no preferred phase.
- The full digit set is the control and it is exact: `zeta_F = zeta`, so `E_F(rho_0) = 0` and both predicted offsets are `0`.
- **Proved**: the discrete identity, the fibre weight, `E_F(rho_0) = zeta_F(rho_0)` for every `c`, and the first-order step at a simple zero. **Verified**: every ratio and every rate, resting on the ladder's propagated bound and on a resolved and uncertified Newton step.

## THE FALSIFICATION PLAN

- Both predictions are computed BEFORE the comparison and from quantities that do not see the measured zero.
- The measured zero is reached by Newton from `rho_0` and accepted only at `abs(zeta_F) < 1e-16`, within `1.5` of `rho_0` and `0.02` clear of the pole lattice; the printed ratio is per zero, never a mean of the two sides. A miss is a zero the trust region does not reach, and the rung's medians are conditioned on that.
- A ratio far from `1` at the dense rungs kills the mechanism. The full-set control must read `E_F = 0` and offset `0`.
- The coupling `zeta_F'(rho_0)/zeta'(rho_0)` is a second reading of the same split and it decides between the candidate constants only if it resolves them; the lab prints its distance to `fill/base` and to `1` side by side.

## WHAT THE SWEEP FOUND

- Nine designs on one ladder in `alpha`, the missing digit always the top one, so base 20 missing one digit is `F = {0..18}` and base 50 missing one digit is `F = {0..48}`. Twelve zeta zeros to `Im s = 56.4462476971` on the seven rungs below `alpha = 0.96` and six to `Im s = 37.5861781588` on the two densest, so the rungs do not share one height. Largest ladder bound anywhere `9.001e-23`.
- THE CONTROL IS EXACT. At the base 2 full set the evaluator reads `abs(zeta_F(rho_0))` between `1.85e-34` and `1.329e-25` at all twelve zeros, so `E_F = 0` and both predicted offsets are `0`, with nothing fitted to make it so.
- THE CONSTANT-FREE STEP IS THE LAW. Its median ratio reads `1.3843088, 1.284225, 1.2481449, 1.2060106, 1.2042502, 0.89075541, 1.0195598, 1.005076, 0.99741809` at `alpha = 0.430676558, 0.5, 0.630929754, 0.792481250, 0.861353116, 0.903089987, 0.954242509, 0.982877878, 0.994835739`, and its band closes at the dense end: largest `abs(ratio - 1)` is `0.14041` at base 20 missing one digit and `0.01734` at base 50 missing one digit, band `[0.94875, 1.14041]` and `[0.98266, 1.01144]`.
- IT SHARPENS WITH THE OFFSET, WHICH IS THE COMB LAW'S OWN SHAPE. Pooled over the ladder, the step's largest `abs(ratio - 1)` runs `0.01734, 0.0508884, 0.193158, 0.83912, 3.32327` over the buckets `abs off < 0.05`, `< 0.1`, `< 0.2`, `< 0.4` and above, on `7, 4, 14, 18, 44` zeros.
- THE `fill/base` READING IS WEAKER AND IS NOT SELECTED. Its median ratio reads `1.4129353, 1.2842149, 1.0955991, 1.1806066, 1.1372453, 1.276577, 1.1347487, 1.0320127, 1.0507079` with largest `abs(ratio - 1)` `0.24964` and `0.0821168` at the two dense rungs, five times looser than the step at base 50. The coupling cannot decide it either: `median abs(coupling - fill/base)` is `0.24057225` and `0.08291158` there against `median abs(coupling - 1)` `0.27619434` and `0.079335871`, so the reading flips between the two rungs while the candidates differ only by `m/base = 0.05` and `0.02`. `c = fill/base` stays the level 1 fibre reading and nothing more.
- THE RATE. Median offset over `m/base = 1 - fill/base` reads `1.6463532, 1.2495026, 1.8345578, 1.5731321, 2.2102406, 1.6634381, 2.2424916, 1.8779239, 1.051349`, and over `1 - alpha` it reads `1.7350628, 1.2495026, 1.6569184, 1.8951686, 3.1883019, 3.432954, 4.9008186, 5.4839111, 4.0716338`. A least squares in the logs, a fit and not a theorem, gives `(m/base)^1.04544` with `R2 0.957842` against `(1-alpha)^0.71691` with `R2 0.944011`; the first column spans `2.13297` and the second `4.38888`, so `m/base` carries the exponent by a factor of `2.05764` inside the `4.28797` that `(1-alpha)/(m/base)` itself spans over this ladder.
- TWO NEW RUNGS BETWEEN `0.954` AND `1`. Base 20 missing one digit at `alpha = 0.9828778777` and base 50 missing one digit at `alpha = 0.9948357391`, all six zeros located at each: median `abs(E_F(rho_0))` `0.11830158` and `0.028066806`, median offset `0.093896196` and `0.021026979`. The paired offset falls fast across the interval the family sweep left empty.
- THE SHADOW STATISTIC ON THE PAIRED ZEROS. Read in the form of the family row, the mean distance from a located design ordinate to the nearest `zeta` ordinate over a quarter of the mean gap between consecutive `zeta` ordinates in the range, the ladder reads `0.81218635, 0.57141859, 0.50488757, 0.37447728, 0.20954319, 0.29197634, 0.14794812, 0.052888241, 0.011098646` and `0` at the full set. The pairing is zeta-zero-first and the family row's is design-zero-first, so this is a parallel ladder and not that row recomputed.

## WHAT DIED

- THE SHADOW DOES NOT SEPARATE ORDINATE FROM REAL PART AT THE ZEROS IT PAIRS. The offset is one complex number with no preferred phase, and per zero `abs(Im off)/abs(Re off)` spans `0.137681` to `6.11895` at base 20 missing one digit and `0.14167` to `18.7749` at base 50 missing one digit.
- Rung by rung the two medians are `0.55734029/0.43095421, 0.49670656/0.28603903, 0.48081533/0.22535435, 0.30633741/0.16748977, 0.12823995/0.36138728, 0.25093621/0.12181102, 0.11574693/0.19332005, 0.058239278/0.049215607, 0.011954894/0.010995712`: the ordinate offset is the larger on six rungs and the smaller on three, and both fall along the ladder.
- The law binds only zeros Newton reaches from a `zeta` zero inside `1.5`, so it constrains no other design zero and says nothing about a census taken design-zero-first.

## WHAT IS OWED

- The census at base 20 and base 50 missing one digit, so the located zeros are named second family or teeth and the new rungs carry the family row's own statistic.
- The count of design zeros with no `zeta` partner per rung, which this lab cannot produce: the pairing is zeta-zero-first and it enumerates no design zero.
- The nine zeros at base 5 `{0,1}`, base 4 `{0,1}` and base 3 `{0,1}` where Newton reaches nothing inside the trust region: predicted offsets there run `0.95618855` to `3.0967393`, outside the first order's own domain.
- A rigorous bound on `zeta_F` near `rho_0`, which turns the Newton step into a Rouche count and the Verified ratios into a proved enclosure.
- The two new rungs at twelve zeros, which would remove the height incomparability with the other seven.

## RUN

- `uv run python research/zeta/lab/py/zeta-shadow/zeta_shadow.py mass ladder` - the discrete identity check, the principal fibre weight `(fill/base)^level` and the exact principal-arc mass at level 1, 2, 3, with the full-set control.
- `uv run python research/zeta/lab/py/zeta-shadow/zeta_shadow.py predict ladder` - per design and per zeta zero the remainder `abs(E_F(rho_0))`, the coupling, the constant-free step, the `fill/base` reading, the located zero and both ratios.
- `uv run python research/zeta/lab/py/zeta-shadow/zeta_shadow.py rungs ladder` - the same plus the rate in `m/base` and `1 - alpha` with its log-log fit, the shadow statistic per rung, and the sharpening of both columns as the offset shrinks.
- The second argument is a family: `ladder`, `old`, `new`. Prints only, writes nothing. The full ladder is about thirteen minutes and `mass` about four minutes (217 s).
