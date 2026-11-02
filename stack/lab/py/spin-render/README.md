# spin-render

- Renders the odd parity-carpet stack with each layer spun about the centre, and measures what survives a change of raster.
- Layer `k` is the carpet `C_n(u, v) = chi_n(u) chi_n(v)` at the `k`-th odd scale `n = 2k - 1`, `n <= 55`, rotated about `(1/2, 1/2)` by `theta_k`, sampled at pixel centres of an `R x R` raster.
- The stack is the paper coverage `P` in `[0, 1]`, the mean over the 28 layers of the not-inked indicator; the ink field is `D = 1 - P`, and `D` takes values `k/28`.
- Peaks are read on `D`. In paper coverage the value `1` is attained on wide regions near the rim and carries no node structure; the isolated features are the ink concentrations, so every peak, plateau and drift number below is a statement about `D`.
- Measurement is restricted to the inscribed disc of radius `1/2` about the centre, so every layer covers every sampled pixel at every angle; the layers themselves are the periodic extension of the parity rule, so no pixel is ever undefined.
- The centre is the fixed point of every rotation and its cell has inradius `1/(2n)` in each layer, so a disc of radius `1/110` about the centre lies inside one cell of all 28 layers whatever the schedule; the centre reads `D = 14/28` in every schedule at every `R`, the 14 scales `n = 3 mod 4` (`main`, `dark_layers`).
- Six schedules (`schedules`): `unspun` at `0`; `degrees` at `k` degrees; `golden` at `k * 137.507764` degrees; `primes` at `0, 2, 3, 5, 7, ...`, layer 1 at `0` and layer `k > 1` at `p_(k-1)`, ending at `103`; `random` at a seeded uniform angle, `default_rng(1)`; `gaussian` at `atan2(b, a)` for `a >= b >= 1` with `a^2 + b^2 = n` and `0` when no such pair exists, which rotates only the 9 scales `5, 13, 17, 25, 29, 37, 41, 45, 53`.

## RUN

- `uv run python research/stack/lab/py/spin-render/spin_render.py`
- From the repo root. One core, about 10 seconds.
- Domain is `R = 256, 512, 1024, 2048` for all six schedules, the fade table at `R = 1024` and layer count `L = 4, 8, 14, 28`, and a `0.02`-wide zoom at effective `R = 51200` on each schedule's top peak.
- Writes the five PNGs beside this file and nothing else.

## WHAT IT PRINTS

- `report` prints, per schedule and per `R`: the peak of `D`, the unit-square area of all disc pixels at that peak, the mean of `P`, the RMS contrast of `D`, the centre value, and the top three of the 20 local maxima with the pixel count of each maximum's connected plateau.
- Maxima are the pixels equal to their `3 x 3` maximum filter (`maximum_filter3`), taken greedily by value with an exclusion radius `0.01` in unit-square coordinates (`peaks`); each is reported at the centroid of its flood-filled plateau (`plateau`), not at its seed pixel.
- `drift` matches each of the 20 maxima at `R` to its nearest at `2R` and prints the median and maximum displacement, in unit coordinates and in pixels at `R`, with the count of the 20 that land within one pixel.
- `diagonal_maximum` computes the unspun maximum and its whole locus exactly, with no raster. `D(u, v) = |S(u) inter S(v)|` with `S(u) = {n : floor(nu) odd}`, so a maximum needs `S(u) = S(v)` of maximal size, and the 1D maxima are read off the `636` breakpoints `k/n` in exact rationals; it returns every interval attaining the maximum and every product cell whose two factors carry the same `S`.
- `exact_variance` sums the paper's rational covariance formula over the `L x L` scale pairs, an expectation for the fade table computed without any render.
- `dark_layers` lists the scales inked at a peak's seed pixel, independently of the accumulator that measured it.

## WHAT IT FINDS

- The unspun global maximum is `18/28`, the 18 scales `n = 3, 5 mod 6`, attained on two diagonal intervals `[1/3, 18/53)` and `[35/53, 2/3)`, each of width `1/159`, and so on four square cells: `S(1 - x) = S(x)` exactly for odd scales, since `floor(n(1 - x)) = n - 1 - floor(nx)` with `n - 1` even, so both intervals carry the same scale set and all four products are maxima. The locus is `4` cells of total area `4/25281 = 0.000158222`, one cell `1/25281 = 3.95554e-05` (`diagonal_maximum`).
- The raster finds all four: peak `18/28` at `(0.33667, 0.33667)` and its three mirrors at every `R`, plateaus of `4, 9, 49, 169` pixels at `R = 256, 512, 1024, 2048`, measured locus area `0.000244141, 0.000137329, 0.000186920, 0.000161171` against the exact `0.000158222`, and the zoom measures the single cell at the top peak as `3.95523e-05` against the exact `3.95554e-05` (`report`, `main`).
- Median drift of the unspun top 20 is `0.71, 0.35, 0.35` pixels over `256 -> 512 -> 1024 -> 2048`, under one pixel at every step, with `20/20` and `19/20` inside one pixel on the last two steps (`drift`).
- Every whole-degree schedule loses the unspun maximum: `degrees` peaks at `14/28`, exactly the centre's value, at all four resolutions; `primes` peaks at `16/28`; neither reaches `18/28`, and the zoom at effective `R = 51200` finds nothing above the raster peak in any schedule (`report`).
- The spun peaks sit on far smaller cells. The area of the single cell at the top peak, from the zoom, is `3.95523e-05` unspun, `1.65596e-05` primes, `9.80453e-06` degrees, `2.27585e-06` random, `1.30005e-06` gaussian, `1.0128e-06` golden; the golden peak is one pixel wide even at `R = 2048` (`report`).
- The discriminator is the single cell, never the whole locus: a schedule whose maximum is lower spreads it over more cells, so the total locus area at the peak runs to `0.000415802` under `degrees` at `R = 2048`, larger than the unspun `0.000158222` at a value four layers higher (`report`).
- Off-centre maxima do stay put for whole-degree schedules: the `primes` peak `16/28` at `(0.19469, 0.15501)` drifts `0.29` pixels over `1024 -> 2048`. A raster cannot read this as an exact coincidence. The stack is piecewise constant on cells of positive area, so any cell wider than a pixel is found at every resolution; plateau stability measures cell area, not lattice coincidence, and the only honest statements a raster supports are resolution-drift statements.
- Top-20 lists reshuffle under refinement wherever the peak value is degenerate: `8/20` to `15/20` of the spun maxima land within one pixel of a maximum at `2R`, against `19/20` unspun, and the maximum displacement runs to `470` pixels where a small cell drops out of the top 20 (`drift`).
- The centre is never the strict global maximum: `18/28` beats it unspun, `16/28` or `17/28` beats it in five schedules, and only `degrees` leaves it tied at the top with `14/28` (`report`).
- The fade constant `c = rms * sqrt(L)` at `R = 1024` over the full square runs `0.309907, 0.390837, 0.427756, 0.460395` at `L = 4, 8, 14, 28` unspun, against the exact `0.309477, 0.389754, 0.426869, 0.458411` from `exact_variance`; the tree's `c = 0.522` is the `L -> infinity` limit, not the value at `L = 28`, and `0.4604` against `0.4584` is the raster's `0.4%` discretisation bias (`main`).
- Spinning does not change the fade law: `c` at `L = 28` reads `0.401417` degrees, `0.427986` golden, `0.425189` primes, `0.423811` random, `0.440940` gaussian, against `0.460395` unspun, all within `13%` and all falling as `1/sqrt(L)` (`main`).
- The mean of `P` over the disc is `0.7682` in every schedule at every resolution, `0.7689` unspun at `R = 2048`; rotation moves no ink (`report`).
- `float32` against `float64` at `R = 512` on the golden schedule: `8` pixels of `262144` differ, by one layer (`main`).

## FIGURES

- `stack-unspun.png`, `stack-degrees.png`, `stack-primes.png`, `stack-gaussian.png`: the `R = 1024` render of each schedule, block-averaged to `512` and quantised to the 29 stack levels, `45820, 61142, 63814, 55405` bytes; `stack-sheet.png` is the four at `256`, `69865` bytes (`save_png`, `contact_sheet`).
- Unspun is a lattice of straight moire rays through the diagonals; the whole-degree fan pleats them into curved bands; the prime and golden schedules wind them into a vortex of arms about the centre, with the centre node visible as a single dark speck.

## WITNESSES

- `farey.md` THE STACK IS AN ADDRESS: the unspun render this study spins is the same field `farey/lab/py/carpet-stack-address` renders at `512 x 512`, and its global maximum `18/28`, on four cells of side `1/159` against the two intervals `[1/3, 18/53)` and `[35/53, 2/3)`, is the closed form read at a denominator-3 node and its mirror.
- The stack's `L^2` fade: `c = 0.522` is the limit of `sqrt(L Var(G_L))`, whose exact values `0.309477 .. 0.458411` at `L = 4 .. 28` this study reproduces from the raster to `0.4%`.
