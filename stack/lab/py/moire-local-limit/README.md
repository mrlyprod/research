# moire-local-limit

- The local limit of the overlay of two `dim 2` designs at sides `N` and `N + h`, both on the unit square, under the difference filter (a point is on when exactly one design inks it), computed exactly from window integrals.
- The window is a box `I x J`; the box mean of a product picture is the product of the 1D means, so every box mean is a bilinear form in the joint law of the digit parities on `I` and on `J`.
- `law1` computes that joint law exactly on `W` windows by merging the cut points `i/N`, `j/(N + h)` and `k/W` over the common denominator `N (N + h) W`; `law2` does the same with the four level-2 parities on the denominator `N^2 (N + 2)^2 W`.
- `check_walsh` checks the Walsh form of the limit, `(1 - sum_S w_S^2 Lambda(hu)^[1 in S] Lambda(hv)^[2 in S])/2` with `Lambda(s) = 1 - 2 dist(s, 2Z)`, against the lag-law kernel in exact rationals for all 16 codes at gaps 1, 2, 3, 4, 6 on 169 points, and prints the five classes the codes fall into.
- `check_direct` rasters the overlay at sides 21 and 23 on the full `483 x 483` grid, with no factorisation, and matches the factorised box means exactly for all 16 codes on 7 by 7 windows.
- `ladder` runs code 7 at gap 2 on 32 by 32 windows of side `1/32` at the pairs `21/23`, `101/103`, `2321/2323`, `23231/23233` and `232321/232323`, printing the worst window error against `H(u,v) = (1 - |1 - 2u| |1 - 2v|)/2`, the error scaled by `N l` and the proved bound `(36/l + 20h)/N`.
- `gaps` runs the gap-`h` limit at `2320/2322` (gcd 2), `2321/2322` (gap 1), `2321/2325` (gap 4) and `2319/2325` (gap 6, gcd 3), with the exact global covariance of the two 1D parities.
- `obstruction` prints the window correlation of the two parities against the window mean of `Lambda(2u)`, and the exact global overlay mean at `2321/2323`.
- `octagon` prints the regular chord octagon level `1 - 1/sqrt 2`, the sag of its curved sides, and two box means at `232321/232323` on windows of side `1/1024` on the diagonal, one at the chord midpoint and one on the arc.
- `level2` computes the level-2 limit law in exact rationals by integrating the lag law over the top-level phase, checks it is quadratic on each quarter, runs the direct box means at `201, 1001, 2001` and `200, 1000, 2000`, checks the closed form `H2 = 5/8 - alpha alpha - gamma gamma - 2 beta beta` at 169 points, prints a nonzero 3 by 3 minor, the parity gap `9/256`, the distance to the naive formula and the kernel rank 3.

## RUN

`uv run python research/stack/lab/py/moire-local-limit/moire_local_limit.py`

Eleven seconds, prints only, writes nothing; the level-2 runs at `2001` and `2000` take six of them.

## WITNESSES

- The local limit on the lag torus for all 16 codes and every gap: `check_walsh`, `check_direct`, `law1`, `kernel1`, `limit_law1`.
- The code 7 limit `H`, its mean `3/8` and the rate: `ladder`, `gaps`.
- Local correlation against global decorrelation: `obstruction`.
- The curved octagon: `octagon`.
- The level-2 limit, its closed form, its parity split and the death of the naive formula: `level2`, `law2`, `limit_law2`, `closed_form`.
