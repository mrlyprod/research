# dilated-receptive-field

- The effective receptive field of a dilated convolution stack read as a digit count, the numbers of [dilations](../../../notes/dilations.md).
- `r_Q(n) = [z^n] prod_(j < L) Q(z^(b^j))` is the weighted base-`b` digit count of the block polynomial `Q`, computed in exact integers by `product`.
- `count` checks `product` against the digit recursion `r_L(n) = sum_(c = n mod b) q_c r_(L-1)((n - c)/b)` on every lag, for `Q = 1 + z + z^2`, `b = 2`, `L = 1..exact`, and for `1 + z`, `2 + z`, `1 + z + z^2 + z^3` at `b = 3` and `1 + (1 + z + z^2)^2` at `L = exact - 4`; checks Stern's recursion against the first `32` terms of A002487, copied from OEIS; then per level `L = 1..top` asserts `r(n) = s(n + 1)` below `2^L`, the mirror, the maximum `F_(L+1)` at the four closed-form lags, the single-path lags `2^k - 1` and their mirrors numbering `2L + 1`, the fold `s(2^L + j) = r(j - 1) + r(j - 1 + 2^L)`, and prints peak over mean `F_(L+1)(2^(L+1) - 1)/3^L` and its step ratio against `2 phi/3`.
- `limit` decides absolute continuity of the depth limit by the exact cyclotomic check in rational arithmetic, `Q` vanishing at `e^(2 pi i k/b^i)` for some `i <= depth` for every `k < b^depth` not divisible by `b`, on `35` block polynomials: uniform `K = 2..8` at `b = 2, 3, 4`, `1 + c(1 + z)` and `1 + c(1 + z + z^2)^2` over a grid of `c`, and `(1 + z + z^2)^2`; asserts the verdict against the Fourier product `prod_i Q(e^(2 pi i t/b^i))/Q(1)` over `t < b^3` not divisible by `b`; checks the real-gain `Q = (z^2 - sqrt2 z + 1)(z^4 + sqrt2 z^2 + 1)(1 + z + z^2)^2` at `b = 2`: its least coefficient, the level at which each odd `k < 2^depth` is covered, `|Q|` at the primitive `2^i`-th roots for `i = 1..4`, and the Fourier product over odd `t < 256`; prints the share of lags carrying half the mass at `L = 8, 12, 16, 20` for four polynomials.
- `gradient` draws random linear stacks, `channels` channels, `levels` levels, `draws` times, for the pure `K = 3` stack at `C sigma^2 = 1` and the residual block `1 + c(1 + z)` at `c = 1/2`, and compares the mean squared gradient over input channels with `r_Q`, printing the largest z-score, the share past `3`, the median relative error and the correlation.
- `copy` trains a linear `K = 3`, `b = 2` stack by full-batch gradient descent on the lag-`n` copy loss `sum_m (f(m) - [m = n])^2`, `f` the end-to-end filter, from `seeds` initialisations at `scale` times the variance-preserving scale, counting steps until `f(n) >= hit`, a run past `cap` read as above the cap in its median, at one lag in `2^(L-1)..2^L - 1` for each `r` in `1, 2, 3, 5, 8, 13, 21, 34, 55` present, and fits log median steps against `log 1/r`.

## RUN

- `uv run python research/dilations/lab/py/dilated-receptive-field/drf.py count` from the repo root, `1.8` seconds.
- `uv run python research/dilations/lab/py/dilated-receptive-field/drf.py limit`, `8.5` seconds.
- `uv run python research/dilations/lab/py/dilated-receptive-field/drf.py gradient`, `1.0` seconds.
- `uv run python research/dilations/lab/py/dilated-receptive-field/drf.py copy`, `4.2` seconds; `copy --scale 0.5`, `112` seconds.
- `uv run python research/dilations/lab/py/dilated-receptive-field/drf.py all`, `16` seconds.
- `--top 20`, `--exact 12`, `--depth 6`, `--channels 8`, `--levels 8`, `--draws 4000`, `--seeds 5`, `--lr 0.005`, `--scale 1`, `--cap 20000` and `--hit 0.5` are the dials.
- Prints only; reads and writes nothing.

## WITNESSES

- The lines of [dilations](../../../notes/dilations.md), every section.
- `count`: `0` mismatches between the product and the digit recursion at `L = 1..12` and on the four other polynomials; the first `32` A002487 terms agree; every assertion holds at `L = 1..20`; maximum `89` at `682, 852, 1194, 1364` and `21` single-path lags at `L = 10`; fold `0` mismatches at `L = 1..20`; peak over mean `3.0853` at `L = 10`, `6.5835` at `L = 20`, step ratio `1.07143` at `L = 3`, `1.14815` at `L = 4`, `1.07869` from `L = 17`.
- `limit`: `35` rows, verdict and Fourier product agreeing; uniform taps absolutely continuous exactly when `b` divides `K`; `1 + c(1 + z)` singular at every `c` on the grid; `1 + c(1 + z + z^2)^2` absolutely continuous only at `c = 1`; the real-gain `Q` has least coefficient `0.5858`, cover levels `3, 4`, `|Q|` up to `11.6569, 0.8284, 0.6863, 9.9446` at the primitive `2, 4, 8, 16`-th roots, Fourier product at most `2.1e-18`; half the mass on `0.2740, 0.2561, 0.2415, 0.2280` of the lags for `K = 3`, `b = 2`, and on `0.0010` at `L = 20` for `1 + (1 + z)/4`.
- `gradient`: median relative error `0.0202`, correlation `0.9989` over `511` lags, share of `|z| > 3` `0.0196`, pure `K = 3`; `0.0224`, `0.9999` over `256` lags, residual.
- `copy`: lags `255, 191, 223, 207, 215, 192, 213, 212` at `r = 1, 2, 3, 5, 8, 13, 21, 34`, median steps `615, 505, 376, 270, 250, 198, 108, 94`, slope `0.551`, correlation `0.984`, no run capped; at `--scale 0.5` median steps `15760, 12768, 8779, 7081, 6233, 6085, 4270, 3788`, one run at `r = 1` capped, slope `0.404`, correlation `0.983` over the `8` lags.
