# major-arcs

- The generator of the chapter The major arcs: one binary, four verbs, no dependencies, up to 8 threads.
- A design is `f_{b,D}(n)`, the parity of the number of base-`b` digits of `n` in `D`, `2 <= b <= 6`: 57 designs named `b:D`, sign `x(n) = prod_i v(n_i)`, `v(d) = (-1)^[d in D]`.
- The model of `S(N) = sum_(n<N) x(n^k)`, `k = 2` or `3`, is the circle method's major arcs summed by trees: `1 + Re sum_gamma A(gamma) U(gamma; N^k)` with `U(gamma; M) = sum_(1<=m<M) x(m) e(m gamma) m^(1/k - 1)/k`. The points `gamma = c/q` are the `b^2 - 1` fractions `c/(b^2 - 1)`, each fixed by `B = b^2`; the arcs are computed below `eta(gamma) = 1 - log_B abs(sigma) = 1/2` (squares) or `1/4` (cubes), `sigma = sum_(t<B) x(t) e(t gamma)`, and the model keeps those below a cut, by default `1/4`.
- `A = (G(-c, q)/q) K`, `G` the degree-`k` Gauss sum, `K = lim_J E_(n in Z_B) prod_(i<J) y(digit_i(n^k))/phi(0)` with `y(t) = x(t) e(t gamma)` and `phi(0) = sigma/B`. For `B = p^t` (`b <= 5`) `K` is a finite sum over the `p`-adic valuation of `n`; for `B = 36` it is a series over the 2-adic and 3-adic valuations, each term the average of a digit product over one congruence class, computed by a memoized digit recursion; the series is cut at `e = 14` (squares) and `e = 10` (cubes).
- `U` is summed exactly up to the truncation of a Taylor series: the interval `[1, M)` splits into base-`b` blocks `[c b^j, (c + 1) b^j)`, `c >= 1`, and on each block `m^s` is expanded about the midpoint in `44` centred moments of `x(r) e(r gamma)`, which a digit recursion carries from level `j - 1` to `j`.
- `model 32` takes about 80 s cold and `cubes 33` about 240 s, on 8 threads. The census walks `n^k` by forward differences, in `u128` for bases 2 and 4 and in limbs of `3^12`, `5^8` and `6^7` for bases 3, 5 and 6, read through digit tables; per dyadic block it keeps the bias, the lags `h <= 7` by bit-sliced counters, and, along the cubes, the 8-blocks of nine streams. A coin rides in bit 57.

## VERBS

- `coefficients`: the arcs of every design along `n^2` and `n^3`, with `eta`, `abs(G)/q`, `K` and `abs(A)`; the mirror pairs; the orbits below `1/4` of period 3 and more, which the model omits; the checks, the tree cut at depth `J` by brute force against the limit `K` and the `e`-terms of the series; the arcs killed by a zero Gauss sum. The arcs are cached in the data folder; about 55 s cold.
- `model K [cut]`: the census along `n^2`, `n < 2^K`, against the model at `N = 2^12, 2^14, ..., 2^K`, design by design, with the residual, the signs and the lit split; the model keeps the arcs with `eta < cut`, default `1/4`. The census is cached.
- `terms K Q`: for `2:1`, `6:5` and `6:12345`, the model term by term, `T_q(N) = sum_(a mod q, (a,q)=1) (G(-a,q)/q) U(a/q; N^2)` for every `q <= Q`, its running sum, the trees and the census.
- `cubes K [cut]`: the same along `n^3`, `n < 2^K`, from `N = 2^11`, with the lags and the chi-square of nine streams.

## RUN

```
cargo run --release --manifest-path research/random/lab/rs/major-arcs/Cargo.toml -- coefficients
cargo run --release --manifest-path research/random/lab/rs/major-arcs/Cargo.toml -- model 32
cargo run --release --manifest-path research/random/lab/rs/major-arcs/Cargo.toml -- terms 32 36
cargo run --release --manifest-path research/random/lab/rs/major-arcs/Cargo.toml -- cubes 33
cargo test --release --manifest-path research/random/lab/rs/major-arcs/Cargo.toml
```
