# squares

- The generator of the chapter Along the squares: one binary, five verbs, no dependencies, up to 8 threads.
- A design is `f_{b,D}(n)`, the parity of the number of base-`b` digits of `n` in `D`, `2 <= b <= 6`, `D` a nonempty subset of `{1..b-1}`: 57 designs, named `b:D`, so `2:1` is Thue-Morse and `6:135` counts the odd digits in base 6. Its sign is `x(n) = prod_i v(n_i)`, `v(d) = (-1)^[d in D]`.
- One `u64` word carries all 57 designs, Rudin-Shapiro (the parity of the `11` blocks) and a coin (the top bit of the `i`-th splitmix64 output). The word of `v` is read through digit tables of `3^13`, `4^10`, `5^9` and `6^8` entries, each entry the XOR of the design bits of its digits, folded over the chunks of `v`.
- A stream is read as every overlapping 8-block at every start `n < N`. The Walsh transform of the 256 counts gives every `T_S = sum_(n<N) prod_(i in S) y(n + i)`: the bias `z_B = T_{0}/sqrt(N)`, the lags `z_Ch = T_{0,h}/sqrt(N)` and the chi-square `X2 = sum_(S != 0) T_S^2/N`, split by the order `abs(S)`. For an iid stream `X2` has mean 255 and sd 38.7, as the census header prints.
- `phi(alpha) = b^-1 sum_(d<b) v(d) e(d alpha)`, `abs(phi)^2 = b^-2 (b + 2 sum_s a_s cos(2 pi s alpha))`, `a_s` the aperiodic autocorrelation of the tile `v`.

## VERBS

- `census K P`: the classes (copies, the folded class `n mod 2`, the twins `D` and `D xor odd` at odd bases, the mirrors `v*(d) = v(b-1) v(b-1-d)`), then the streams `y(n) = x(n^2)`, `n < 2^K`, and `x(p_i)` along the first `2^P` primes from a segmented sieve; the tables at the last checkpoint, the growth at every second power of 2, and the law table: the best periodic orbit of `abs(phi)` up to period `log_b 2^20`, its restriction to `q != 2 (mod 4)` and to squarefree `q`, the correlation dimension `D2`, the lit lags, the mirror pairs, Spearman ranks, the open class. `census 30 26` takes about 15 s on 8 threads.
- `eta`: for every design a bracket on the Fourier sup exponent `eta = -lim (1/L) log_b sup_alpha prod_(j<L) abs(phi(b^j alpha))`. The lower end is a max-plus certificate on the `b^r <= 2^17` cylinders of depth `r`: the sup of `abs(phi)^2` on a cylinder is bounded by its centre value plus `2 pi sum_s s abs(a_s)/(b^2 b^r)`, Howard's policy iteration picks a positive vector `w`, and `max_i G_i max_j w_j/w_i`, `j` over the `b` cylinders that `b I_i` covers, rounded up, bounds the growth rate. The upper end is the product along a periodic orbit `a/(b^k - 1)`, the best of every orbit up to period `log_b 2^20` and the optimal cycle of the certificate. Every float is an interval rounded outward after each operation, cosines by a Taylor series past `x^30` on `[0, pi/4]`, logarithms by the `atanh` series, and the ends are printed outward, the lower floored and the upper ceiled at 6 places. It also prints `D2 = 3 - log_b rho(T)`, `T` the integer `2 x 2` matrix that carries `(sum_h A_L(h)^2, sum_h A_L(h) A_L(h+1))` from level `L` to `L + 1`, `A_L` the aperiodic autocorrelation of the first `b^L` terms. About 6 s.
- `diagonal K P`: two models of the lags against the stream. Along the squares, `sum_(n<N) C_x(2hn + h^2)`, `C_x` the autocorrelation of the design from its `2 x 2` digit matrices `M_s = b^-1 [[al_s, al_(b-s)], [al_(s+1), al_(b-s-1)]]`, summed exactly by a digit walk over `n < N` that carries `2hn + h^2`; along the primes, `sum_i x(r) x((r + g) mod b) C_x(floor((r + g)/b))`, `r = p_i mod b`, `g = p_(i+h) - p_i`, which is exact when the digits above the last are averaged. `diagonal 24 22` takes about 4 s.

- `speed K`: the decay exponent of every strongly lit statistic along `n^2`: for each design with `eta < 1/4` and each of its bias and lags with `abs(z) >= 8` at `2^K`, the least-squares slope of `log2 abs(T_S/N)` over every `N = 2^k`, `k` in `18..K` and in `K-6..K`, against `-2 eta` and `-2 e_h`, with the mean and largest gap. `speed 30` takes about 13 s.

- `thue K`: Thue-Morse alone along `n^2`, the squares in 128 bits and updated by `+2n + 1`, the parity of `n^2` the popcount of its two halves XORed: the bias and the lags 3 and 6 at every `N = 2^k`, `k <= K`, and their decay slopes over `2^18..2^30`, `2^18..2^K` and `2^30..2^K`. It agrees with the tally to `2^22` (a test). `thue 40` takes about 8 minutes on 8 threads.

## RUN

```
cargo run --release --manifest-path research/random/lab/rs/squares/Cargo.toml -- census 30 26
cargo run --release --manifest-path research/random/lab/rs/squares/Cargo.toml -- eta
cargo run --release --manifest-path research/random/lab/rs/squares/Cargo.toml -- diagonal 24 22
cargo run --release --manifest-path research/random/lab/rs/squares/Cargo.toml -- speed 30
cargo run --release --manifest-path research/random/lab/rs/squares/Cargo.toml -- thue 40
cargo test --release --manifest-path research/random/lab/rs/squares/Cargo.toml
```
