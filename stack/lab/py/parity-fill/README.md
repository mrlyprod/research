# parity-fill

- Computes the inked fraction of the unit square under the parity blend of the odd carpet stack: the XOR of the layers `C_n(u, v) = chi_n(u) chi_n(v)`, `n = 1, 3, ..., N`, with `chi_n(u) = 1` iff `floor(n u)` is odd.
- The parity blend is the sum of the layers folded to its parity, `mrly_math::spin::Blend::Parity`, which on `0/1` layers is exactly their XOR.
- Layer `n = 1` is blank on `[0, 1)^2`, so the XOR is carried by the scales `3, 5, ..., N` alone and the expansion runs over subsets of those.
- The exact route is the subset expansion. With `s_n = 1 - 2 C_n` in `{+1, -1}`, `prod_n s_n = (-1)^(sum_n C_n) = 1 - 2 XOR`, so `fill = (1 - E[prod_n s_n])/2`; expanding `prod_n (1 - 2 C_n) = sum_S (-2)^|S| prod_(n in S) C_n` and splitting `prod_(n in S) C_n(u, v) = [prod_(n in S) chi_n(u)][prod_(n in S) chi_n(v)]` factorises each mean, so `E[prod_(n in S) C_n] = m_S^2` with `m_S` the measure of the set of `u` in `[0, 1)` where every `floor(n u)`, `n in S`, is odd. Hence `fill(N) = (1 - sum_S (-2)^|S| m_S^2)/2`.
- `m_S` is exact and rational: every `chi_n`, `n in S`, is constant on each cell of the grid of side `1/lcm(S)`, so `m_S` is a cell count over `lcm(S)`. All `2^L` masses, `L` the number of layers, come from one pass over the grid of `lcm(3, 5, ..., 21) = 14549535`, which records the scale set of each cell and then sums over supersets (`mass_table`).
- The independent comparison gives layer `n` the mean `p_n = ((n-1)/(2n))^2` and pretends the layers are independent Bernoulli, which would make the XOR fill `(1 - prod_n (1 - 2 p_n))/2` (`fill_independent`).

## RUN

- `uv run python research/stack/lab/py/parity-fill/parity_fill.py`
- From the repo root. One core, about three seconds, numpy and `fractions` only.
- Domain is the exact fill at `N = 3, 5, ..., 21`, the literal 2D cross-check at `N <= 9`, a `4096 x 4096` raster check at every `N`, the spun readings at `N = 55` on rasters `R = 256, 512, 1024, 2048`, and the whole-degree increment sweep `0..90` at `R = 512`.

## WHAT IT PRINTS

- `section_exact` prints `fill(N)` as an exact fraction and a decimal for each `N`, beside the `4096 x 4096` raster mean and the gap, and beside the literal 2D XOR count at `N <= 9`.
- `fill_exact` sums the `2^L` subset terms in integers over the common denominator `lcm^2`, so no rounding enters; `fill_literal` XORs the actual layer products on the `lcm x lcm` cell grid and counts inked cells, sharing no line with the expansion.
- `section_independent` prints the exact fill, the independent-Bernoulli fill, their difference, the two deviations from `1/2` and their ratio, then the successive ratios of each deviation.
- `section_spun` prints the disc fill at `N = 55` for eight schedules, the resolution table for three of them, and fill against layer count with the exact and square-raster columns beside the disc readings.
- `section_sweep` prints the disc fill of the fixed-increment stack at every whole increment from `0` to `90`, its extremes, the eyes, and two symmetry checks.

## WHAT IT FINDS

- The expansion is exact and the literal count confirms it: `fill(N)` reads `1/9`, `53/225`, `3524/11025`, `36284/99225`, `19619/51975`, `117419647/289864575`, `109067744/289864575`, `17006699344/45107387325`, `6812188030619/19244451701475` and `1114185811873/2749207385925` at `N = 3, 5, ..., 21`, matching the literal 2D XOR count at all four `N <= 9` with zero mismatches and the `4096 x 4096` raster mean at every `N` within `2.42e-04` (`fill_exact`, `fill_literal`, `raster_row`).
- The fill is not monotone in `N`: it climbs to `0.405084502` at `N = 13`, falls to `0.376271381` at `N = 15`, falls again to `0.353981924` at `N = 19`, and returns to `0.405275287` at `N = 21` (`section_exact`).
- Most subsets contribute nothing: `587` of the `1024` subsets of `{3, 5, ..., 21}` have `m_S = 0`, the smallest being `{3, 5, 7}`, since `floor(3u)`, `floor(5u)` and `floor(7u)` are never all odd at the same `u` (`mass_table`).
- That zero is the whole mechanism. At two coprime scales the layers are independent and the exact and independent fills agree to the digit, `1/9` at `N = 3` and `53/225` at `N = 5`; the triple `{3, 5, 7}` has joint mass `0` against the independent `2/35`, so the two fills part company at `N = 7` and the exact fill is below the independent one at every `N` from `7` to `21`, by `-1.31e-02` up to `-1.40e-01` (`section_independent`).
- The independent approximation decays as advertised and the exact quantity does not. The independent deviation `1/2 - fill` has successive ratios `0.680000, 0.632653, ..., 0.546485` falling toward `1/2`, while the exact deviation's ratios run `0.680000, 0.682044, 0.744755, 0.912184, 0.774630, 1.303566, 0.993894, 1.187399, 0.648719` and rise above `1` twice; the ratio of the two deviations grows from `1.000000` to `29.337331` over `L = 2` to `L = 11` (`section_independent`).
- On the square the raster fill of the unspun stack still climbs toward `1/2`: `0.320427, 0.375175, 0.424397, 0.444069` at `L = 4, 8, 14, 28` at `R = 1024`, the first two against the exact `0.319637188` and `0.376271381`, the raster reading high by `7.9e-04` and low by `1.1e-03` (`section_spun`).
- Spinning pushes the fill to `1/2`. On the disc at `N = 55` and `R = 1024` the fill reads `0.472574` unspun, `0.490510` at one degree per layer index, `0.505291` prime degrees, `0.503331` golden, and `0.500821, 0.491967, 0.500658, 0.502868` at the eyes `90/j`, `j = 2, 3, 4, 5`; the unspun stack is the outlier at `2.74e-02` from `1/2` and no spun schedule exceeds `9.49e-03` (`section_spun`).
- The raster reading is good to about `3e-03`: the same schedule over `R = 256, 512, 1024, 2048` reads `0.478511, 0.475647, 0.472574, 0.473205` unspun, `0.485506, 0.492253, 0.490510, 0.489691` at one degree, and `0.492967, 0.494657, 0.491967, 0.492404` at the eye `90/3`. The unspun offset survives that band by a factor of nine; the individual spun readings mostly do not (`section_spun`).
- The eyes do not stand out. Over the whole-degree sweep at `R = 512` the `89` nonzero increments span `0.490165` to `0.511161`, the eyes at `18`, `30` and `45` reading `0.500185`, `0.494657` and `0.500282`, inside that band and not at either end; the extremes are `0` and `90` at `0.475647` on one side and `23` and `67` at `0.511161` on the other (`section_sweep`).
- Two symmetries close the sweep. Every layer is invariant under the quarter turn about the centre, since `chi_n(1 - u) = chi_n(u)` at odd `n` and `C_n` is symmetric in its two coordinates, so increment `90` reproduces the unspun stack exactly and increments `d` and `90 - d` are mirror images: the sweep reads `fill(90) = fill(0) = 0.475647` and matches `45` of the `46` mirror pairs to the bit, the one exception differing by `9.71e-06` in `float32` (`section_sweep`).

## WITNESSES

- `franel-and-eisenstein.md` THE SPUN PICTURE: the same field this study folds by parity is the one that study averages, scales `n = 1, 3, ..., 55`, `28` layers, disc-masked raster; the mean-blend paper coverage `0.7682` there and the parity fill here are two blends of one stack.
- `stacking-and-spinning.md` SELECTING THE SCALES: the moire correlation law's coprime independence is what makes the exact and independent fills agree at `N = 3` and `N = 5`, and the joint mass `m_{3,5,7} = 0` is what makes them disagree at `N = 7`.
- `mrly_math::spin::Blend::Parity`: the blend this study measures, the sum of the copies folded to its parity.
