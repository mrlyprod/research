# two-bases

- The generator of the moire limit: one binary, five verbs, no dependencies, single-threaded.
- `C_t` is the correlation of Thue-Morse, `C3` that of `t3`, the base-3 signed design of code `{2}`, and `C_z = C_t C3` that of the moire XOR `z = t . t3`.
- `C_t` is exact: the vector `w = 3 . 2^L (C_t(h), C_t(h + 1))` starts at `(3, -1)` and reads the `L` binary digits of `h` from the top, `w -> (2x, -(x + y))` on a `0` and `w -> (-(x + y), 2y)` on a `1`, so `C_t(h) = x/(3 . 2^L)`; past `2^64` it runs in a 192-bit two's-complement integer, which holds every value below `10^40` exactly.
- `C3` is exact: `C3(3q) = C3(q)`, `C3(3q + 1) = -C3(q + 1)/3`, `C3(3q + 2) = -C3(q)/3`, ending at `C3(0) = 1` or `C3(1) = 0`, so `C3(h)` is `0` or `+-3^(-m)`.
- `r_2(h)` is the number of runs of the binary word of `h`; `w_3(h)` the number of nonzero balanced-ternary digits.

## VERBS

- `lemma`: the form `Q(x, y) = 20x^2 + 20xy + 41y^2` and `rho^2 = 103/200` in exact rationals: the coefficients `a`, `d` and the discriminant `Delta = 4ad - B^2` of `rho^2 Q(u) - Q(A(t) u)` as polynomials in `t`, their concavity, their values at `t = -1/2` and `t = 1/4`, and their real roots on `[-1/2, 1/4]` counted by an exact Sturm sequence; the form norm of `A(t)` on a grid; `K`; the run pattern `(1, 1, 2)` and its spectral radius. Under 0.1 s.
- `decay`: every `h < 2^25` against `|C_t(h)| <= 1.0093 . 0.717636^(r_2(h))` and `|C3(h)| <= 3^(-w_3(h)/2)`; the maximum of `|C_t|` at each run count; the Cesaro means of `C_t^2` and `C_z^2` and the proved exponents `theta = log2(400/303)`, `log3(9/5)` and `1 - theta - log3(9/5)`; the family of run pattern `(1, 1, 2)` to 30 runs; `C_t(2^k - 1)` against its closed form. About 3 s.
- `windows`: the maximum of `|C_z|` on every `[2^k, 2^(k+1))`, `k <= 24`, with its argument and exact value; the lags `1 <= h < 2^25` with `abs C_z >= 10^(-3)`, their count, the last of them with its exact value, how many sit within `0.02` of the threshold, and the sup over `[2^22, 2^25)`. About 1.5 s.
- `changes`: every `1 <= h < 3^12` at `b = 2` and `b = 3` against `BC(h, b) <= 2 l_b(h) - 1`, `BC` the number of digit changes of the base-`b` word; `l_b` by the digit recursion `l(bq) = l(q)`, `l(bq + r) = min(r + l(q), b - r + l(q + 1))`, exact since each shift by `b` of the lowest coefficient toward `{r, r - b}` saves `b` and moves the rest by one unit, which costs at most `1`, against a breadth-first search over `[-3^13, 3^13]`; the equality cases; `BC(h, 2) = r_2(h) - 1`, `l_2 <= r_2 + 1` and `l_3 = w_3` asserted; `rho^4 < 1/3` and the constants of the rate. Under 0.3 s.
- `sparse R W D`: every `h < 10^D` with `r_2(h) <= R` and `w_3(h) <= W`. It walks every balanced-ternary word with at most `W` nonzero digits below the cut, keeps those whose low 64 bits have at most `R` digit changes (a necessary condition), then checks the run count, the cut and the weight exactly; balanced ternary is unique, so nothing is counted twice. It prints the list (its last 48 when longer than 20000), the counts per bit-length block, the largest element, its exact `C_z`, the maximum of `|C_z|` above each `2^k`, and what sits within `0.02` of the cut. `sparse 6 6 40` walks `1.45e10` words in about 10 s; `sparse 8 8 20` walks `2.08e10` in about 24 s.

## RUN

```
cargo run --release --manifest-path research/random/lab/rs/two-bases/Cargo.toml -- lemma
cargo run --release --manifest-path research/random/lab/rs/two-bases/Cargo.toml -- decay
cargo run --release --manifest-path research/random/lab/rs/two-bases/Cargo.toml -- windows
cargo run --release --manifest-path research/random/lab/rs/two-bases/Cargo.toml -- changes
cargo run --release --manifest-path research/random/lab/rs/two-bases/Cargo.toml -- sparse 6 6 40
cargo run --release --manifest-path research/random/lab/rs/two-bases/Cargo.toml -- sparse 8 8 20
```
