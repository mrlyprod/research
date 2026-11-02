# readings

- The generator of the readings: one binary, five verbs, no dependencies, single-threaded.
- A reading is a design read along a polynomial or the consecutive primes: `y(n) = x(P(n))`, `x` one of Thue-Morse `t` (the parity of the binary ones), Rudin-Shapiro `r` (the parity of the adjacent pairs `11`) and the ternary parity `t3` (the parity of the ternary digits equal to `2`), as signs `+-1`; the primes are indexed from `p_0 = 2`.
- A battery pair is one of the 44 dilation pairs `(n, Bn + r)`, `2 <= B <= 9`, `r < B`, or the two kernel pairs `(bn, b^2 n)`, `b = 2, 3`. It is a near-head at the design's base `b` when `P(psi(n)) = b^L (P(phi(n)) + q) + s` with `0 <= s < b^L`, decided on the coefficients of the composed polynomials.
- The local constant of a near-head is the limit of the mean of `x(m + q) x(m)` along `m = P(phi(n))`: the classes of `n` modulo `b^J` (`2^(J+1)` for Rudin-Shapiro, whose rule reads the bit above the window) whose addition stays inside `J` digits give `kappa_J`, the others `eps_J`, and the limit lies in `[kappa_J - eps_J, kappa_J + eps_J]`; `eps_J` is printed rounded up.
- The battery follows the design ladder's protocol: Pearson correlations over `N` pairs, the shift test `max over 1 <= h <= 64`, the dilation test `max over the 44 pairs`, the kernel test `max over b = 2, 3`, the square `mean of y(n) y(2n) y(3n) y(6n)`, LZ78 as `c (log2 c + 1)/N` bits per bit for `c` phrases. Streams are bit-packed; the primes come from a segmented odd sieve.
- The head scan reads `y(phi^j n) y(phi^(j+1) n)` over `n < N` for `phi(n) = lambda n + mu`, `lambda in {2, 3, 4, 8, 9}`, `mu < lambda`, `j = 0, 1, 2`, and reports it when it is periodic in `n` with period at most 81: plain at period 1, twisted above.

## VERBS

- `identities`: the identities between readings at every `n < 2^24`; the kernel links `t((2^k n + 1)^2) = 1 + t(2^(k-1) n^2 + n)` and `t3((3^k n + 1)^2) = t3(3^k n^2 + 2n)`; the kernel witnesses `n = 2^(k-1) + 1` for every pair `3 <= k < k' <= 40` and `n = (3^k + 1)/2` for every pair `1 <= k < k' <= 24`. Runtime under 1 s.
- `battery`: three designs times eight readings (`n`, `n^2`, `n^2 + 1`, `n^2 + n`, `(2n + 1)^2`, `2n^2 + 1`, `n^3`, `p_n`) and six probes (`n^2 + 2`, `n^2 + 3`, `n^2 + 4`, `n^2 + 5`, `(n + 1)^2 + 1`, `n^3 + n`) at `N = 2^15`, with LZ78 at `2^14`, `2^16`, `2^18` and the head scan; the `splitmix64` control. Runtime about 13 s.
- `local`: every near-head battery pair of the 13 polynomial readings and the three designs, with `kappa_J` and `eps_J` at `J = 20` (`13` at base 3) beside the raw mean and the Pearson correlation at `N = 2^22` (`3^14`), a multiple of the class period, so the raw mean must lie in the enclosure; then the pair law `x(n^2 + c) x(n^2)`, `c = 1..8`, at `J = 26` (`16` at base 3), with the nearest rational `p/q`, `q <= 64`, and the bound `abs(kappa_J - p/q) + eps_J` on its distance from the limit. Runtime about 25 s.
- `falsify [e]`: the squares, the readings with no near-head among the battery pairs, the primes and the control at `N = 2^20, 2^22, ..., 2^e` (default `e = 24`), with `z = max|r| sqrt(N)` per test, the three largest dilation pairs, and the slopes of `log2 max|r|` against `log2 N`, overall and step by step. `falsify 26` takes about 130 s, of which 23 s is the sieve of 603979904 primes.
- `leak`: the sums `U_r(H)` of the Thue-Morse correlations `C_t(h)` over `h < H`, `h = r mod 3`, by the exact recursion `U(2H) = M U(H) + b` at `H = 2^k`, checked against direct sums for `k <= 24`, with `U_0(2^k) (2/3)^k` and the mean over the multiples of 3 to `k = 61`, and the closed form `(16/45)(3/2)^k + 8/15` at even `k`, `+ 7/15` at odd `k`, checked in exact integers on `6 . 2^k U(2^k)` for `k = 1..70`. Runtime under 1 s.

## RUN

```
cargo run --release --manifest-path research/random/lab/rs/readings/Cargo.toml -- identities
cargo run --release --manifest-path research/random/lab/rs/readings/Cargo.toml -- battery
cargo run --release --manifest-path research/random/lab/rs/readings/Cargo.toml -- local
cargo run --release --manifest-path research/random/lab/rs/readings/Cargo.toml -- falsify 26
cargo run --release --manifest-path research/random/lab/rs/readings/Cargo.toml -- leak
```
