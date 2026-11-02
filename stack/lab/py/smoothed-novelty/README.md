# smoothed-novelty

- The novelty meter: the totient sum read through a test function, and the exponent the error carries against the same sum read through a sharp cutoff.
- The object is `S_f(y) = sum_n phi(n) f(ny)` for `f` supported on `[1, 2]`, main term `(6/pi^2) F(2) y^-2` with `F(s) = int f(u) u^(s-1) du`, and the meter is `E_f(y) = y^2 S_f(y) - (6/pi^2) F(2)` read as a power of `q = y^2`.
- Three test functions: the indicator of `[1, 2]`, the `C^2` bump `64 (u-1)^3 (2-u)^3` whose Mellin transform is the closed form `64 sum_k a_k (2^(s+k) - 1)/(s+k)`, and the `C^infinity` bump `exp(4 - 1/((u-1)(2-u)))` whose Mellin transform is a 2000-node Gauss-Legendre quadrature.
- `totients` sieves `phi(n)` to `3 * 10^7` in numpy, dividing each prime out of every multiple and finishing with the one prime factor above the square root; checked against brute-force gcd counts to 2000 and against the Mobius route `sum_d mu(d) T(floor(x/d))` at `x = 10^6`, both exact.
- The grid is `y = 2^-j` for `j` from `8` to `23.5` in steps of `1/16`, 249 samples, the largest window summing the `1.19 * 10^7` scales between `2^23.5` and `2^24.5`.
- `slopes` fits `log` of the per-octave root mean square of `E_f` against `log q` over all 16 octaves and over the lower and upper 8, an octave being the samples with `j` in `[k, k + 1)`, 16 of them except the last, `[23, 23.5]`, which holds 9, each placed at `j = k + 1/2`; it prints each slope with its window and residual, and the least-squares slope over every sample beside it.
- `zeros_from_pari` calls `gp -q` once for the 138 zeros of `zeta` to height 300 with `zeta(rho - 1)` and `zeta'(rho)` at each, and the explicit formula `E_f(y) = 2 Re sum_rho F(rho) zeta(rho - 1)/zeta'(rho) y^(2 - rho)` is compared with the measured `E_f` over the whole grid at 1, 10, 30 and 138 zeros; the least-squares power of `|c_rho|` against `gamma` on the 138 zeros is printed for each bump.
- The summation noise is printed at the smallest `y` as the gap between numpy's pairwise sum and `math.fsum` of the same products.

## WHAT IT PRINTS

- The per-octave slopes in `q`: indicator `0.5023` (residual `0.179`; windows `0.5296` and `0.4883`, residuals `0.161` and `0.169`), `C^2` bump `0.7498` (residual `0.092`; windows `0.7553` and `0.7487`), `C^infinity` bump `0.7471` (residual `0.115`; windows `0.7474` and `0.7460`).
- The explicit formula: with 138 zeros the `C^infinity` bump's `E_f` is reproduced to a relative `2.0e-06` over the grid, the `C^2` bump's to `2.3e-04`, the gap being the tail of its coefficients, whose printed least-squares power against `gamma` is `-3.22` on the 138 zeros; one zero alone gives `0.50` and `0.33`.
- `|E_f|/y^(3/2)` stays between `1.3e-04` and `0.558` for the `C^infinity` bump against `2 sum |c_rho| = 0.755` over the same 138 zeros.
- The sharp cutoff: the totient error `sum_{n <= x} phi(n) - 3x^2/pi^2` jumps by `phi(p) - 3(2p - 1)/pi^2` across the prime `p`, `392073.4` at `p = 1000003` and `3920735.7` at `p = 10000019`.
- Summation noise at `j = 23.5`: `5.6e-17` against `|E_f| = 5.4e-12`.

## RUN

- `uv run python research/stack/lab/py/smoothed-novelty/smoothed_novelty.py`
- From the repository root, one core, 15 seconds: 2 for the sieve, 8 for the 498 smoothed sums, under 1 for PARI. Needs `gp` on the path.
- Nothing is written to disk.

## WITNESSES

- `lean-and-meter.md`, "The novelty meter": the three slopes with their windows, the explicit formula and its residuals, the jump at a prime.
- The Mellin transform of the `C^2` bump at `s = 2` is `24/35`, printed from the closed form beside the exact value.
