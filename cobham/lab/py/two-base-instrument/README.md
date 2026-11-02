# two-base-instrument

- The falsification test that must pass before any two-base instrument is built: does the count of a two-base intersection, read on a `log N` grid, carry both lattice frequencies at once.
- The cell is `dim 1`, base `3` digits `{0, 1}` against base `5` digits `{0, 1, 2}`, a multiplicatively independent pair of positive budget `log_3 2 + log_5 3 - 1`.
- `C(N)` is the number of `1 <= n < N` whose base-3 digits all lie in `{0, 1}` and whose base-5 digits all lie in `{0, 1, 2}`. The enumeration walks the base-3 digits from the top down and prunes a branch as soon as a base-5 digit that the branch has already settled leaves `{0, 1, 2}`: after `j` base-3 digits the value lies in an interval of width `3^(level-j)`, and every base-5 digit position on which the two ends of that interval agree is settled for the whole branch.
- The mechanics of the reading are inherited from `dimensions/lab/py/complex-dimensions`, function `periodogram_peak`: `g(u) = ln C(e^u)` detrended in `u`, a Blackman window, the mean removed, the power on a direct grid over `[0.5, 14]`, each local maximum refined by golden section.
- The decision rule is this study's own and is weaker than the one `dimensions.md` states. That page takes the global argmax of the band and asks it to land within `1%` of the prediction; a two-base count has to be read at frequencies that are not its loudest, so here a prediction counts present when the nearest local maximum lies within `1%` of it and carries at least `10x` the median power of the band.
- A power ratio is a property of its window. Every reading prints the window span, the number of maxima above `10x` the median and the fraction of the band their `1%` neighbourhoods cover, which is the chance a prediction is met by position alone, and ratios are comparable only across windows of the same span and the same crowding.
- Every reading is printed twice, detrended by a straight line and by a cubic, because the local exponent of the two-base count is still falling at the height reached and a linear detrend leaves curvature at the low-frequency end.
- Each strongest maximum is tested against the whole lattice `m 2 pi / ln 3 + n 2 pi / ln 5` with `|m|, |n| <= 8`, so a maximum matching no prediction is a maximum matching no small combination of the two.
- Three one-base controls and one block model calibrate the rule: the base-3 factor alone, the base-5 factor alone, the multiplicatively dependent pair base `3` digits `{0, 1}` against base `9` digits `{0, 1, 3}`, which the collapse theorem turns into the single base-9 design `{0, 1, 3}` of exponent exactly `log_9 3 = 1/2`, and the block model `C3(N) C5(N) / N`, which multiplies the two band structures with no joint arithmetic.
- The detector reads a count that grows, and this one does not grow everywhere. A member of the base-3 design with `k+1` digits lies in `[3^k, (3^(k+1)-1)/2]` and a member of the base-5 design with `i+1` digits in `[5^i, (5^(i+1)-1)/2]`, so the designs occupy `0.369070` and `0.569323` of their own decades and `C` is exactly constant on every stretch where the two bands miss.
- Exact integer arithmetic in the enumeration and in the band test; floating point only in the detector. Pure Python, no dependency outside the standard library.

## VERBS

- `control M` enumerates the cell to `3^level` for `level 0..M` twice, once by the pruned walk and once by testing every integer below `3^level` against both digit rules, prints whether the two lists agree, and prints the first members of the cell. Runtime `0.03` s at `M = 9`.
- `cell level` prints the height `3^level`, the budget rounded up, the node count, the hit count, the enumeration time, and then the reading on the intersection over the full window and over its upper half, followed by the base-3 and base-5 controls. Runtime `27` s at `level 44` and `82` s at `level 47`, of which the walk prints `78` s, one core.
- `collapse M` runs the dependent pair to `3^M` through the same pipeline. Runtime `6` s at `M = 28`, where the hit count `4782968` is within a factor `1.2` of the hit count of the independent cell at `level 44`.
- `ladder A B` reads the full window at every height from `3^A` to `3^B` under both detrends and prints one row per height and detrend: the hit count, the hits new in that decade, the window span, the crowding of the band, and the error, power ratio and verdict at both frequencies. Runtime `4` min `19` s at `A = 38`, `B = 47`.
- `blocks level` prints the decades of `3` below `3^level` on which the two design bands do not meet, so that the cell has no member there and `C` is exactly constant across the decade, and then runs the same reading on the block model `C3(N) C5(N) / N`. Runtime `0.8` s at `level 47`.

## RUN

```
uv run python research/cobham/lab/py/two-base-instrument/detector.py control 9
uv run python research/cobham/lab/py/two-base-instrument/detector.py cell 44
uv run python research/cobham/lab/py/two-base-instrument/detector.py cell 47
uv run python research/cobham/lab/py/two-base-instrument/detector.py collapse 28
uv run python research/cobham/lab/py/two-base-instrument/detector.py ladder 38 47
uv run python research/cobham/lab/py/two-base-instrument/detector.py blocks 44
uv run python research/cobham/lab/py/two-base-instrument/detector.py blocks 47
```

## WITNESSES

- Each row is a number this study prints and stands behind; the page section they belong to is `gasket-three-base-and-two-frequencies.md`.
- The enumerator agrees with direct digit filtering at every level to `3^9`, where both routes return the same `42` values, and the first members of the cell are `1, 10, 12, 27, 30, 31, 36, 37, 252, 255, 256, 280`.
- The height `3^44 = 984770902183611232881`, the hit count `5667470`, and at `3^47` the hit count `19042219`.
- The budget `0.313536` rounded up, against a fitted exponent `0.323301` over the whole window and `0.295978` over its upper half at `3^44`, `0.318728` and `0.292184` at `3^47`.
- The base-3 control at `3^44`, window span `38.05`, `11` maxima above `10x` the median covering `0.113` of the band: `2 pi / ln 3 = 5.719202` met at `5.719220`, error `0.000%`, power `1.7e7` times the median under the linear detrend and `1.9e7` under the cubic; `2 pi / ln 5` at error `1.115%` and power `35.9` under the linear detrend, and at error `0.966%` and power `39.2` under the cubic, which passes the rule at a frequency absent by construction.
- The same control at `3^47`, window span `41.34`, `10` loud maxima covering `0.104`: `2 pi / ln 5` at error `3.151%` and power `2.54` under the linear detrend and `2.791%` and `2.6` under the cubic, so the pass at `3^44` does not recur.
- The base-5 control at `3^44`, window span `38.49`, `4` loud maxima covering `0.048`: `2 pi / ln 5 = 3.903963` met at `3.903959`, error `0.000%`, power `4.0e6`; `2 pi / ln 3` at error `1.060%` and power `2.21` under the linear detrend and at `0.963%` and `1.95` under the cubic, inside the tolerance and held out by the power gate alone.
- The same control at `3^47`, window span `41.79`, `17` loud maxima covering `0.150`: `2 pi / ln 3` at error `1.012%` and power `1.02` under the linear detrend and `1.000%` and `0.941` under the cubic.
- The dependent control at `3^28 = 9^14`, window span `17.47`, `4` loud maxima covering `0.042`: hit count exactly `3^14 - 1 = 4782968`, fitted exponent `0.497809` against the exact `1/2`, `2 pi / ln 9 = 2.859601` met at `2.859889`, error `0.010%`, power `886` times the median, and `2 pi / ln 3` met at error `0.003%`, power `2.5e3`.
- The block model at `3^44`, window span `28.26`, `4` loud maxima covering `0.043`: both frequencies met at errors `0.002%` and `0.003%`, and the nearest maximum to the difference frequency `1.815239` at `1.895843`, error `4.440%`, power `0.588` times the median, so block structure puts nothing loud where the cell's strongest maximum sits.
- The block model at `3^47`, window span `31.55`, `4` loud maxima covering `0.043`: `2 pi / ln 3` met at error `0.001%` and power `1.33e6` times the median, `2 pi / ln 5` met at `0.004%` and `9.79e5`, so the two band structures alone carry both frequencies at five to six orders of magnitude above the floor.
- The independent cell at `3^44`, window span `27.41`, `4` loud maxima covering `0.020`: `2 pi / ln 3` not met, nearest maximum `5.638828`, error `1.405%`, power `24`; `2 pi / ln 5` met at `3.921970`, error `0.461%`, power `16.5`.
- The independent cell at `3^47`, window span `30.70`, `4` loud maxima covering `0.020`: `2 pi / ln 3` not met, nearest maximum `5.634120`, error `1.488%`, power `22.5`; `2 pi / ln 5` met at `3.927166`, error `0.594%`, power `15.7`.
- The independent cell on the upper half of its window at `3^44`, counts `149469` to `5667470`, span `13.71`: neither frequency met, but `2 pi / ln 3` sits at error `0.825%` under the linear detrend and `0.791%` under the cubic, inside the tolerance and held out by the power gate alone at `5.55` and `6.1` times the median.
- The ladder from `3^38` to `3^47`: both frequencies are met on the full window under both detrends at `level 38, 39, 40, 41`, and `2 pi / ln 3` is not met from `level 42` up, its error rising from `1.034%` to `1.574%` while its power ratio stays between `22` and `25`.
- The decades `[3^j, 3^(j+1))` below `3^47` whose two design bands do not meet, so that the cell has no member in them, are `j = 1, 4, 17, 20, 23, 26, 36, 39, 42, 45`, which is `10` of the `47`; the ladder confirms each of the three inside its range with a height whose new hits are `0`, at `level 40, 43, 46`.
- The strongest maximum of the independent cell is `1.702087` at `3^44` and `1.697925` at `3^47`, at `57` and `68` times the median, and the nearest `m 2 pi / ln 3 + n 2 pi / ln 5` with `|m|, |n| <= 8` is `(1, -1) = 1.815239` at errors `6.233%` and `6.463%`.
