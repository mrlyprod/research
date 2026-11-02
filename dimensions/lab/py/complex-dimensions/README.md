# Complex Dimensions

- The complex dimensions of the four one-base designs, base 3 `{0,2}`, base 5 `{0,2,4}`, base 15 `{0,4,10,14}` and base 15 `{0,2,4,10,12,14}`: the 81 predicted poles `s = dimension + 2 pi i m/ln(base)`, `m = -40..40`, checked against `1 - fill base^(-s)`, the numerator `D(s)` at each, and `D(1) = 1 - fill/base`.
- The two-ratio control with ratios `1/3` and `1/5`: its `dimension`, its 21 Moran roots in `Re [-3,3]`, `Im [-40,40]` by Newton from a grid, the winding number of `1 - 3^(-s) - 5^(-s)` over that box, and the worst offset of the roots from each candidate spacing.
- Composition: alternating bases multiplies into the product base, checked as integers and as 256 exact intervals.
- The box count `N(eps)` of each cover, `g(u) = ln N(e^-u) - dimension u` on `u in [ln(1/0.03), ln(1e6)]` at 3000 points, its Blackman periodogram peak on a direct DFT grid, and the variance explained by folding `u` modulo `ln 3`, `ln 5`, `ln 15` in 40 bins, for the four designs, the two-ratio control and a Thue-Morse aperiodic control.
- The inner tube `V(eps)` in closed form, `M(eps) = eps^(dimension-1) V(eps)` to `u = 60`, its swing per window, the periodicity defect `M(eps) = M(eps/p)`, the decay of the two-ratio swing, and the Cantor limit profile `2^(1-dimension) (t^(dimension-1) + t^dimension)` against the measured tube.
- `carpet_tube.py`: the Sierpinski carpet's tube by digit level, `V(eps) = sum_m 8^(m-1) h(3^(-m), eps)`, the hole census cell by cell to level 5, the closed form against the level sum as exact rationals, the limit profile `G(t)` on its two branches with its seam, its extrema as quadratic roots, its swing, the level sum at `u in [50, 60]`, and a level-6 distance transform equal to the closed form as rationals.

## RUN

- `uv run python research/dimensions/lab/py/complex-dimensions/complex_dimensions.py`
- About 15 seconds; prints only, writes nothing.
- `uv run python research/dimensions/lab/py/complex-dimensions/carpet_tube.py`
- About 1 second; prints only, writes nothing, and asserts the closed profile against the direct sum.

## WITNESSES

- dimensions.md (the lead and the last section) name two scratch passes; this study is the one pass that replaces both.
- string-to-composition.md all 81 poles at `m = -40..40` kill the denominator to `5e-14`, minimum `|D(s)| = 0.500000` for the Cantor design, `D(1) = 1 - fill/base` for all four.
- string-to-composition.md `dimension` 0.630930, 0.682606, 0.511916, 0.661642 and `omega` 5.719202, 3.903963, 2.320188, 2.320188.
- string-to-composition.md every periodogram peak within 1% of `2*pi/ln(base)`.
- string-to-composition.md folding 0.192, 0.030, 0.016; 0.018, 0.509, 0.034; 0.007, 0.022, 0.667; 0.009, 0.030, 0.534; 0.082, 0.047, 0.038; 0.007, 0.012, 0.115.
- string-to-composition.md the alternation equals base 15 `{0,4,10,14}` and the 256 intervals are identical.
- string-to-composition.md `dimension = 0.518370`, 21 roots and winding number 21, real parts `-0.699926` to `0.518370`, offsets 0.43, 0.17, 0.38.
- measurability-and-sponge.md minimum `2.494975716` at `t = 0.584963`, maximum `2.583040469`, swing 3.53%, measured tube within `4.9e-9` at the minimum.
- measurability-and-sponge.md swings flat from `u = 15` to `u = 60` (`eps = 8.8e-27`), periodicity at the own base to `1e-9` or better, the two-ratio swing decaying `3.79%` to `0.42%`.
- research/claims/:360 the poles of `1/(1 - fill base^(-s))` on one vertical line of period `2 pi/ln(base)`.
- measurability-and-sponge.md the sponge's hole boundary is not in the sponge: `(1/2, 1/2, 1)` at distance `1/6`, `(2/3, 1/2, 5/6)` at distance `1/18`; read off the digit rule, no generator.
- measurability-and-sponge.md `dimension = 1.892789`; every hole ringed by filled cells and the holes exhausting the complement at levels 3, 4, 5 (`carpet_tube.py`, HOLES).
- measurability-and-sponge.md the two branches of `G(t)`, the level sum equal to the closed form at 96 exact rationals (`carpet_tube.py`, TUBE and PROFILE).
- measurability-and-sponge.md seam `(44/35) 2^(2-dimension) = 1.354123517`, ends `379/280`, maximum `1.35561708227` (rounded up) at `t = 0.429638415`, minimum `1.3506702097` (rounded down) at `t = 0.692137253`, swing `0.36625%`, the level sum within `1e-9` of `G` on `u in [50, 60]`, the raster equal to `V(21/729)` as rationals (`carpet_tube.py`, PROFILE, MEASURED, RASTER).
- measurability-and-sponge.md the area of `F_eps` inside `Gamma`, `4 eps/3 - 4 eps^2` on `(0, 1/6]`: the level-1 term of the sum, read off the closed form; the theorem numbers are the arXiv version's.
- measurability-and-sponge.md the class profile `G(t) = t^(dimension-dim) sum_j fill^(j-1) base^(-j dim) h(t base^j)`: at base 3, dim 2, fill 8 it is the two branches the script checks; no other member is run.
